//! Inner distance calculation for paired-end RNA-seq reads.
//!
//! Reimplementation of RSeQC's `inner_distance.py`. Computes the inner distance
//! between read pairs, classifying each pair by gene model overlap and producing
//! histogram and R plot output.

use anyhow::{Context, Result};
use log::info;
use rust_htslib::bam::{self, Read as BamRead};
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Write};

// ============================================================================
// Data Structures
// ============================================================================

/// Result of inner distance analysis for a single BAM file.
#[derive(Debug)]
pub struct InnerDistanceResult {
    /// Per-pair detail records: (read_name, distance_or_NA, classification)
    pub pairs: Vec<PairRecord>,
    /// Histogram bin counts
    pub histogram: Vec<(i64, i64, u64)>,
    /// Total number of read pairs processed
    pub total_pairs: u64,
}

/// A single read pair's inner distance record.
#[derive(Debug)]
pub struct PairRecord {
    /// Read name
    pub name: String,
    /// Inner distance (None if different chromosomes)
    pub distance: Option<i64>,
    /// Classification string
    pub classification: String,
}

/// Exon bitset: tracks which genomic positions are exonic per chromosome.
/// Uses a simple sorted interval list with point-query via binary search.
#[derive(Debug, Default)]
struct ExonBitset {
    /// Sorted, non-overlapping intervals per chromosome
    intervals: HashMap<String, Vec<(u64, u64)>>,
}

impl ExonBitset {
    /// Add an exon interval for a chromosome.
    fn add(&mut self, chrom: &str, start: u64, end: u64) {
        self.intervals
            .entry(chrom.to_string())
            .or_default()
            .push((start, end));
    }

    /// Merge overlapping intervals after all have been added.
    fn merge(&mut self) {
        for intervals in self.intervals.values_mut() {
            intervals.sort();
            let mut merged: Vec<(u64, u64)> = Vec::new();
            for &(start, end) in intervals.iter() {
                if let Some(last) = merged.last_mut() {
                    if start <= last.1 {
                        last.1 = last.1.max(end);
                    } else {
                        merged.push((start, end));
                    }
                } else {
                    merged.push((start, end));
                }
            }
            *intervals = merged;
        }
    }

    /// Check if a chromosome exists in the bitset.
    fn has_chrom(&self, chrom: &str) -> bool {
        self.intervals.contains_key(chrom)
    }

    /// Count exonic bases in the range [start, end).
    fn count_exonic_bases(&self, chrom: &str, start: u64, end: u64) -> u64 {
        let intervals = match self.intervals.get(chrom) {
            Some(iv) => iv,
            None => return 0,
        };

        // Find the first interval that could overlap [start, end)
        let idx = match intervals.binary_search_by_key(&start, |&(s, _)| s) {
            Ok(i) => i,
            Err(i) => {
                if i > 0 {
                    i - 1
                } else {
                    0
                }
            }
        };

        let mut count = 0u64;
        for &(iv_start, iv_end) in &intervals[idx..] {
            if iv_start >= end {
                break;
            }
            let overlap_start = start.max(iv_start);
            let overlap_end = end.min(iv_end);
            if overlap_start < overlap_end {
                count += overlap_end - overlap_start;
            }
        }
        count
    }
}

/// Transcript interval tree: maps chromosome to sorted transcript ranges.
#[derive(Debug, Default)]
struct TranscriptTree {
    /// Per-chromosome sorted list of (start, end, name)
    transcripts: HashMap<String, Vec<(u64, u64, String)>>,
}

impl TranscriptTree {
    /// Add a transcript range. Note: RSeQC has a bug where the first transcript
    /// per chromosome is dropped (the dict is created but `add_interval` is only
    /// called in the `else` branch). We reproduce this bug for compatibility.
    fn add(&mut self, chrom: &str, start: u64, end: u64, name: &str, first_seen: bool) {
        if !first_seen {
            self.transcripts
                .entry(chrom.to_string())
                .or_default()
                .push((start, end, name.to_string()));
        } else {
            // Create entry but don't add the interval (RSeQC bug compatibility)
            self.transcripts.entry(chrom.to_string()).or_default();
        }
    }

    /// Sort intervals for efficient querying.
    fn sort(&mut self) {
        for transcripts in self.transcripts.values_mut() {
            transcripts.sort_by_key(|&(s, _, _)| s);
        }
    }

    /// Find transcript names that overlap a point.
    fn find_overlapping(&self, chrom: &str, point: u64) -> HashSet<String> {
        let mut result = HashSet::new();
        if let Some(transcripts) = self.transcripts.get(chrom) {
            for &(start, end, ref name) in transcripts {
                if start <= point && point < end {
                    result.insert(name.clone());
                }
            }
        }
        result
    }
}

// ============================================================================
// BED12 Parsing
// ============================================================================

/// Parse BED12 file and build exon bitsets and transcript tree.
fn parse_bed12(bed_path: &str) -> Result<(ExonBitset, TranscriptTree)> {
    let file = std::fs::File::open(bed_path)
        .with_context(|| format!("Failed to open BED file: {}", bed_path))?;
    let reader = BufReader::new(file);

    let mut exon_bitset = ExonBitset::default();
    let mut transcript_tree = TranscriptTree::default();
    let mut first_seen_chroms: HashSet<String> = HashSet::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            continue;
        }

        let chrom = fields[0].to_uppercase();
        let tx_start: u64 = fields[1].parse().context("Invalid txStart")?;
        let tx_end: u64 = fields[2].parse().context("Invalid txEnd")?;
        let name = fields[3];
        let block_count: usize = fields[9].parse().context("Invalid blockCount")?;

        // Skip single-exon transcripts (matches RSeQC: `if int(fields[9] == 1): continue`)
        if block_count == 1 {
            continue;
        }

        let block_sizes: Vec<u64> = fields[10]
            .trim_end_matches(',')
            .split(',')
            .map(|s| s.parse::<u64>())
            .collect::<std::result::Result<Vec<_>, _>>()
            .context("Invalid blockSizes")?;
        let block_starts: Vec<u64> = fields[11]
            .trim_end_matches(',')
            .split(',')
            .map(|s| s.parse::<u64>())
            .collect::<std::result::Result<Vec<_>, _>>()
            .context("Invalid blockStarts")?;

        // Add exons to bitset
        for i in 0..block_count.min(block_sizes.len()).min(block_starts.len()) {
            let exon_start = tx_start + block_starts[i];
            let exon_end = exon_start + block_sizes[i];
            exon_bitset.add(&chrom, exon_start, exon_end);
        }

        // Add transcript to tree (with RSeQC first-per-chrom bug)
        let is_first = first_seen_chroms.insert(chrom.clone());
        let tx_name = format!("{}:{}:{}-{}", name, chrom, tx_start, tx_end);
        transcript_tree.add(&chrom, tx_start, tx_end, &tx_name, is_first);
    }

    exon_bitset.merge();
    transcript_tree.sort();

    Ok((exon_bitset, transcript_tree))
}

// ============================================================================
// CIGAR Helpers
// ============================================================================

/// Extract intron blocks from CIGAR (N operations).
/// Matches RSeQC's `bam_cigar.fetch_intron()`: S operations do NOT advance position.
#[allow(dead_code)]
fn fetch_intron_blocks(record: &bam::Record) -> Vec<(u64, u64)> {
    let mut introns = Vec::new();
    let mut chrom_st = record.pos() as u64;

    for op in record.cigar().iter() {
        use rust_htslib::bam::record::Cigar::*;
        match op {
            Match(len) => chrom_st += *len as u64,
            Ins(_) => {}
            Del(len) => chrom_st += *len as u64,
            RefSkip(len) => {
                let intron_start = chrom_st;
                chrom_st += *len as u64;
                introns.push((intron_start, chrom_st));
            }
            SoftClip(_) => {} // Does NOT advance in fetch_intron
            _ => {}
        }
    }

    introns
}

/// Extract exon blocks from CIGAR (M operations).
/// Matches RSeQC's `bam_cigar.fetch_exon()`: S operations DO advance position.
fn fetch_exon_blocks(record: &bam::Record) -> Vec<(u64, u64)> {
    let mut exons = Vec::new();
    let mut chrom_st = record.pos() as u64;

    for op in record.cigar().iter() {
        use rust_htslib::bam::record::Cigar::*;
        match op {
            Match(len) => {
                let start = chrom_st;
                chrom_st += *len as u64;
                exons.push((start, chrom_st));
            }
            Ins(_) => {}
            Del(len) => chrom_st += *len as u64,
            RefSkip(len) => chrom_st += *len as u64,
            SoftClip(len) => chrom_st += *len as u64, // Advances in fetch_exon (RSeQC bug)
            _ => {}
        }
    }

    exons
}

// ============================================================================
// Main Analysis Function
// ============================================================================

/// Compute inner distances for paired-end reads in a BAM file.
///
/// # Arguments
/// * `bam_path` - Path to the BAM file
/// * `bed_path` - Path to the BED12 reference gene model
/// * `sample_size` - Maximum number of read pairs to process
/// * `mapq_cut` - Minimum MAPQ threshold
/// * `lower_bound` - Lower bound for histogram
/// * `upper_bound` - Upper bound for histogram
/// * `step` - Histogram bin width
///
/// # Returns
/// An `InnerDistanceResult` with per-pair details and histogram data.
#[allow(clippy::too_many_arguments)]
pub fn inner_distance(
    bam_path: &str,
    bed_path: &str,
    sample_size: u64,
    mapq_cut: u8,
    lower_bound: i64,
    upper_bound: i64,
    step: i64,
    reference: Option<&str>,
) -> Result<InnerDistanceResult> {
    info!("Processing inner distance for: {}", bam_path);

    // Parse reference gene model
    let (exon_bitset, transcript_tree) = parse_bed12(bed_path)?;

    // Open BAM/CRAM file
    let mut bam = if let Some(ref_path) = reference {
        let mut reader = bam::Reader::from_path(bam_path)
            .with_context(|| format!("Failed to open BAM file: {}", bam_path))?;
        reader.set_reference(ref_path)?;
        reader
    } else {
        bam::Reader::from_path(bam_path)
            .with_context(|| format!("Failed to open BAM file: {}", bam_path))?
    };
    let header = bam::Header::from_template(bam.header());

    // Build tid -> chromosome name mapping (uppercased)
    let tid_to_chrom: HashMap<i32, String> = {
        let header_view = bam.header();
        let mut map = HashMap::new();
        for (tid, name_bytes) in header_view.target_names().iter().enumerate() {
            let name = String::from_utf8_lossy(name_bytes).to_uppercase();
            map.insert(tid as i32, name);
        }
        map
    };

    let mut pairs: Vec<PairRecord> = Vec::new();
    let mut distances: Vec<i64> = Vec::new(); // for histogram
    let mut pair_num: u64 = 0;

    // Suppress the header unused warning
    let _ = header;

    for result in bam.records() {
        let record = result.context("Failed to read BAM record")?;

        // Flag filtering (matches RSeQC exactly)
        if record.flags() & 0x200 != 0 {
            continue; // QC fail
        }
        if record.flags() & 0x400 != 0 {
            continue; // Duplicate
        }
        if record.flags() & 0x100 != 0 {
            continue; // Secondary
        }
        if record.is_unmapped() {
            continue;
        }
        if !record.is_paired() {
            continue; // Single-end reads
        }
        if record.is_mate_unmapped() {
            continue;
        }
        if record.mapq() < mapq_cut {
            continue;
        }

        let read1_start = record.pos() as u64;
        let read2_start = record.mpos() as u64;

        // Skip if mate has lower position (already processed)
        if (record.mpos() as u64) < read1_start {
            continue;
        }

        // Skip if same position and this is read1
        if read2_start == read1_start && record.is_first_in_template() {
            continue;
        }

        pair_num += 1;
        if pair_num > sample_size {
            break;
        }

        let read_name = String::from_utf8_lossy(record.qname()).to_string();

        // Check for different chromosomes
        let read1_tid = record.tid();
        let read2_tid = record.mtid();
        if read1_tid != read2_tid {
            pairs.push(PairRecord {
                name: read_name,
                distance: None,
                classification: "sameChrom=No".to_string(),
            });
            continue;
        }

        let chrom = match tid_to_chrom.get(&read1_tid) {
            Some(c) => c.clone(),
            None => continue,
        };

        // Compute read1_end from CIGAR end position (accurate reference end).
        // RSeQC uses pos + qlen + intron_size, but cigar end_pos is more accurate
        // as it properly accounts for D/I operations.
        let read1_end = record.cigar().end_pos() as u64;

        // Step 1: Compute inner distance
        let inner_dist: i64 = if read2_start >= read1_end {
            // No overlap: positive or zero inner distance
            (read2_start - read1_end) as i64
        } else {
            // Overlap: negative inner distance
            // Count exonic positions of read1 in overlap region [read2_start, read1_end)
            let exon_blocks = fetch_exon_blocks(&record);
            let mut overlap_count: i64 = 0;
            for (ex_start, ex_end) in &exon_blocks {
                let ov_start = (*ex_start).max(read2_start);
                let ov_end = (*ex_end).min(read1_end);
                if ov_start < ov_end {
                    overlap_count += (ov_end - ov_start) as i64;
                }
            }
            -overlap_count
        };

        // Step 2: Check transcript membership (always, regardless of distance sign)
        // RSeQC uses find(read1_end-1, read1_end) for read1 and find(read2_start, read2_start+1)
        let read1_genes = transcript_tree.find_overlapping(&chrom, read1_end.saturating_sub(1));
        let read2_genes = transcript_tree.find_overlapping(&chrom, read2_start);
        let common_genes: HashSet<_> = read1_genes.intersection(&read2_genes).collect();

        let classification: String;

        if common_genes.is_empty() {
            // Different transcripts — use genomic distance (even if negative)
            classification = "sameTranscript=No,dist=genomic".to_string();
        } else if inner_dist > 0 {
            // Same transcript, positive distance: detailed exon analysis
            if !exon_bitset.has_chrom(&chrom) {
                classification = "unknownChromosome,dist=genomic".to_string();
            } else {
                let exonic_bases = exon_bitset.count_exonic_bases(&chrom, read1_end, read2_start);

                if exonic_bases as i64 == inner_dist {
                    classification = "sameTranscript=Yes,sameExon=Yes,dist=mRNA".to_string();
                } else if exonic_bases > 0 {
                    classification = "sameTranscript=Yes,sameExon=No,dist=mRNA".to_string();
                    // Use mRNA distance instead of genomic
                    let mrna_dist = exonic_bases as i64;
                    pairs.push(PairRecord {
                        name: read_name,
                        distance: Some(mrna_dist),
                        classification,
                    });
                    distances.push(mrna_dist);
                    continue;
                } else {
                    classification = "sameTranscript=Yes,nonExonic=Yes,dist=genomic".to_string();
                }
            }
        } else {
            // Same transcript, inner_distance <= 0: overlap
            classification = "readPairOverlap".to_string();
        }

        pairs.push(PairRecord {
            name: read_name,
            distance: Some(inner_dist),
            classification,
        });
        distances.push(inner_dist);
    }

    info!("Total read pairs processed: {}", pair_num);

    // Build histogram
    let histogram = build_histogram(&distances, lower_bound, upper_bound, step);

    Ok(InnerDistanceResult {
        pairs,
        histogram,
        total_pairs: pair_num,
    })
}

// ============================================================================
// Histogram
// ============================================================================

/// Build histogram bins for inner distances.
fn build_histogram(
    distances: &[i64],
    lower_bound: i64,
    upper_bound: i64,
    step: i64,
) -> Vec<(i64, i64, u64)> {
    let mut bins: Vec<(i64, i64, u64)> = Vec::new();
    let mut bin_start = lower_bound;
    while bin_start < upper_bound {
        let bin_end = bin_start + step;
        let count = distances
            .iter()
            .filter(|&&d| d > bin_start && d <= bin_end)
            .count() as u64;
        bins.push((bin_start, bin_end, count));
        bin_start = bin_end;
    }
    bins
}

// ============================================================================
// Output Writers
// ============================================================================

/// Write the per-pair detail file.
pub fn write_detail_file(result: &InnerDistanceResult, path: &str) -> Result<()> {
    let mut file = std::fs::File::create(path)
        .with_context(|| format!("Failed to create detail file: {}", path))?;

    for pair in &result.pairs {
        let dist_str = match pair.distance {
            Some(d) => d.to_string(),
            None => "NA".to_string(),
        };
        writeln!(file, "{}\t{}\t{}", pair.name, dist_str, pair.classification)?;
    }

    Ok(())
}

/// Write the frequency (histogram) file.
pub fn write_freq_file(result: &InnerDistanceResult, path: &str) -> Result<()> {
    let mut file = std::fs::File::create(path)
        .with_context(|| format!("Failed to create frequency file: {}", path))?;

    for &(bin_start, bin_end, count) in &result.histogram {
        writeln!(file, "{}\t{}\t{}", bin_start, bin_end, count)?;
    }

    Ok(())
}

/// Write the R script for plotting.
pub fn write_r_script(
    result: &InnerDistanceResult,
    prefix: &str,
    path: &str,
    step: i64,
) -> Result<()> {
    let mut file = std::fs::File::create(path)
        .with_context(|| format!("Failed to create R script: {}", path))?;

    // Build bin center and count vectors
    let bin_centers: Vec<String> = result
        .histogram
        .iter()
        .map(|&(start, _end, _count)| format!("{}", start as f64 + step as f64 / 2.0))
        .collect();
    let counts: Vec<String> = result
        .histogram
        .iter()
        .map(|&(_start, _end, count)| count.to_string())
        .collect();

    let num_bins = result.histogram.len();

    writeln!(file, "out_file = '{}'", prefix)?;
    writeln!(file, "pdf('{}.inner_distance_plot.pdf')", prefix)?;
    writeln!(
        file,
        "fragsize=rep(c({}),times=c({}))",
        bin_centers.join(","),
        counts.join(",")
    )?;
    writeln!(file, "frag_sd = sd(fragsize)")?;
    writeln!(file, "frag_mean = mean(fragsize)")?;
    writeln!(file, "frag_median = median(fragsize)")?;
    writeln!(
        file,
        "write(x=c(\"Name\",\"Mean\",\"Median\",\"sd\"), sep=\"\\t\", file=stdout(),ncolumns=4)"
    )?;
    writeln!(
        file,
        "write(c(out_file,frag_mean,frag_median,frag_sd),sep=\"\\t\", file=stdout(),ncolumns=4)"
    )?;
    writeln!(
        file,
        "hist(fragsize,probability=T,breaks={},xlab=\"mRNA insert size (bp)\",main=paste(c(\"Mean=\",frag_mean,\";\",\"SD=\",frag_sd),collapse=\"\"),border=\"blue\")",
        num_bins
    )?;
    writeln!(file, "lines(density(fragsize,bw={}),col='red')", 2 * step)?;
    writeln!(file, "dev.off()")?;

    Ok(())
}

/// Write the summary text file.
pub fn write_summary(result: &InnerDistanceResult, path: &str) -> Result<()> {
    let mut file = std::fs::File::create(path)
        .with_context(|| format!("Failed to create summary file: {}", path))?;

    // Count classifications
    let mut same_chrom_no = 0u64;
    let mut same_transcript_no = 0u64;
    let mut same_exon_yes = 0u64;
    let mut same_exon_no = 0u64;
    let mut non_exonic = 0u64;
    let mut unknown_chrom = 0u64;
    let mut overlap = 0u64;

    for pair in &result.pairs {
        match pair.classification.as_str() {
            "sameChrom=No" => same_chrom_no += 1,
            "sameTranscript=No,dist=genomic" => same_transcript_no += 1,
            "sameTranscript=Yes,sameExon=Yes,dist=mRNA" => same_exon_yes += 1,
            "sameTranscript=Yes,sameExon=No,dist=mRNA" => same_exon_no += 1,
            "sameTranscript=Yes,nonExonic=Yes,dist=genomic" => non_exonic += 1,
            "unknownChromosome,dist=genomic" => unknown_chrom += 1,
            "readPairOverlap" => overlap += 1,
            _ => {}
        }
    }

    writeln!(file, "Total read pairs: {}", result.total_pairs)?;
    writeln!(file, "Different chromosome: {}", same_chrom_no)?;
    writeln!(file, "Different transcript: {}", same_transcript_no)?;
    writeln!(file, "Same exon (mRNA distance): {}", same_exon_yes)?;
    writeln!(file, "Different exon (mRNA distance): {}", same_exon_no)?;
    writeln!(file, "Non-exonic (genomic distance): {}", non_exonic)?;
    writeln!(file, "Unknown chromosome: {}", unknown_chrom)?;
    writeln!(file, "Read pair overlap: {}", overlap)?;

    // Compute mean/median/sd from histogram
    if !result.histogram.is_empty() {
        let total: u64 = result.histogram.iter().map(|&(_, _, c)| c).sum();
        if total > 0 {
            let step = if result.histogram.len() > 1 {
                result.histogram[1].0 - result.histogram[0].0
            } else {
                5
            };
            let mean: f64 = result
                .histogram
                .iter()
                .map(|&(s, _, c)| (s as f64 + step as f64 / 2.0) * c as f64)
                .sum::<f64>()
                / total as f64;
            writeln!(file, "Mean inner distance: {:.2}", mean)?;
        }
    }

    Ok(())
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exon_bitset_merge() {
        let mut bitset = ExonBitset::default();
        bitset.add("CHR1", 100, 200);
        bitset.add("CHR1", 150, 300);
        bitset.add("CHR1", 400, 500);
        bitset.merge();

        let intervals = &bitset.intervals["CHR1"];
        assert_eq!(intervals.len(), 2);
        assert_eq!(intervals[0], (100, 300));
        assert_eq!(intervals[1], (400, 500));
    }

    #[test]
    fn test_exon_bitset_count() {
        let mut bitset = ExonBitset::default();
        bitset.add("CHR1", 100, 200);
        bitset.add("CHR1", 300, 400);
        bitset.merge();

        // Fully within first exon
        assert_eq!(bitset.count_exonic_bases("CHR1", 100, 200), 100);
        // Spanning intron
        assert_eq!(bitset.count_exonic_bases("CHR1", 150, 350), 100); // 50 + 50
                                                                      // Fully intronic
        assert_eq!(bitset.count_exonic_bases("CHR1", 200, 300), 0);
        // Unknown chrom
        assert_eq!(bitset.count_exonic_bases("CHR2", 100, 200), 0);
    }

    #[test]
    fn test_histogram() {
        let distances = vec![-5, -3, 0, 2, 7, 12];
        let hist = build_histogram(&distances, -10, 15, 5);
        assert_eq!(hist.len(), 5);
        assert_eq!(hist[0], (-10, -5, 1)); // -5 (d > -10 && d <= -5)
        assert_eq!(hist[1], (-5, 0, 2)); // -3, 0 (d > -5 && d <= 0)
        assert_eq!(hist[2], (0, 5, 1)); // 2 (d > 0 && d <= 5)
        assert_eq!(hist[3], (5, 10, 1)); // 7
        assert_eq!(hist[4], (10, 15, 1)); // 12
    }

    #[test]
    fn test_transcript_tree_first_per_chrom_bug() {
        let mut tree = TranscriptTree::default();
        tree.add("CHR1", 100, 200, "gene1", true); // First = dropped
        tree.add("CHR1", 150, 300, "gene2", false); // Second = added
        tree.sort();

        let genes = tree.find_overlapping("CHR1", 175);
        assert!(!genes.contains("gene1")); // First gene was dropped
        assert!(genes.contains("gene2"));
    }
}
