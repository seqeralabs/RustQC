//! Inner distance calculation for paired-end RNA-seq reads.
//!
//! Reimplementation of RSeQC's `inner_distance.py`. Computes the inner distance
//! between read pairs, classifying each pair by gene model overlap and producing
//! histogram and R plot output.

use crate::gtf::Gene;
use anyhow::{Context, Result};
use coitrees::{COITree, Interval, IntervalTree};
use indexmap::IndexMap;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, Write};

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
pub struct ExonBitset {
    /// Sorted, non-overlapping intervals per chromosome (0-based half-open)
    pub intervals: HashMap<String, Vec<(u64, u64)>>,
}

impl ExonBitset {
    /// Build an `ExonBitset` from parsed GTF gene annotations.
    ///
    /// Includes exons from ALL transcripts (including single-exon), converting
    /// from GTF 1-based inclusive coordinates to 0-based half-open coordinates.
    /// This matches upstream RSeQC's `getExon()` which does not filter by exon
    /// count. Chromosome names are uppercased.
    pub fn from_genes(genes: &IndexMap<String, Gene>) -> Self {
        let mut bitset = ExonBitset::default();
        for gene in genes.values() {
            for tx in &gene.transcripts {
                let chrom = tx.chrom.to_uppercase();
                for &(start, end) in &tx.exons {
                    // GTF: 1-based inclusive → BED-style 0-based half-open
                    bitset.add(&chrom, start - 1, end);
                }
            }
        }
        bitset.merge();
        bitset
    }

    /// Build an `ExonBitset` from a BED12 file.
    ///
    /// Parses the BED12 file and extracts exon blocks from ALL transcripts
    /// (including single-exon). This matches upstream RSeQC's `getExon()`
    /// which does not filter by block count.
    pub fn from_bed(bed_path: &str) -> Result<Self> {
        let reader = crate::io::open_reader(bed_path)
            .with_context(|| format!("Failed to open BED file: {}", bed_path))?;

        let mut bitset = ExonBitset::default();

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
            let block_count: usize = fields[9].parse().context("Invalid blockCount")?;

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

            for i in 0..block_count.min(block_sizes.len()).min(block_starts.len()) {
                let exon_start = tx_start + block_starts[i];
                let exon_end = exon_start + block_sizes[i];
                bitset.add(&chrom, exon_start, exon_end);
            }
        }

        bitset.merge();
        Ok(bitset)
    }

    /// Add an exon interval for a chromosome.
    pub fn add(&mut self, chrom: &str, start: u64, end: u64) {
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
    pub fn has_chrom(&self, chrom: &str) -> bool {
        self.intervals.contains_key(chrom)
    }

    /// Count exonic bases in the range [start, end).
    pub fn count_exonic_bases(&self, chrom: &str, start: u64, end: u64) -> u64 {
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

/// Transcript interval tree: maps chromosome to a cache-oblivious interval tree
/// for fast overlap queries. Replaces the previous linear-scan implementation.
#[derive(Default)]
pub struct TranscriptTree {
    /// Per-chromosome COITree with transcript name index as metadata.
    trees: HashMap<String, COITree<u32, u32>>,
    /// All transcript names, indexed by the metadata stored in the tree.
    names: Vec<String>,
    /// Temporary builder: per-chromosome interval list (only used during construction).
    pending: HashMap<String, Vec<(u64, u64, String)>>,
}

impl TranscriptTree {
    /// Build a `TranscriptTree` from parsed GTF gene annotations.
    ///
    /// Includes ALL transcripts (including single-exon) to match upstream
    /// RSeQC's `getTranscriptRanges()`. Converts from GTF 1-based inclusive
    /// coordinates to 0-based half-open coordinates. Chromosome names are
    /// uppercased. Transcript names use the format `"name:CHROM:start-end"`.
    ///
    /// Replicates the upstream RSeQC bug where the first transcript per
    /// chromosome is dropped (the `if/else` in the transcript_ranges builder
    /// creates the `Intersecter` but skips `add_interval` for the first entry).
    pub fn from_genes(genes: &IndexMap<String, Gene>) -> Self {
        let mut tree = TranscriptTree::default();
        // Collect all transcripts, then drop the first per chromosome
        let mut per_chrom: HashMap<String, Vec<(u64, u64, String)>> = HashMap::new();
        for gene in genes.values() {
            for tx in &gene.transcripts {
                let chrom = tx.chrom.to_uppercase();
                // GTF: 1-based inclusive → BED-style 0-based half-open
                let start = tx.start - 1;
                let end = tx.end;
                let name = format!("{}:{}:{}-{}", tx.transcript_id, chrom, start, end);
                per_chrom.entry(chrom).or_default().push((start, end, name));
            }
        }
        // Replicate the upstream bug: skip the first transcript per chromosome
        for (chrom, transcripts) in &per_chrom {
            for (start, end, name) in transcripts.iter().skip(1) {
                tree.add(chrom, *start, *end, name);
            }
        }
        tree.build();
        tree
    }

    /// Build a `TranscriptTree` from a BED12 file.
    ///
    /// Includes ALL transcripts (including single-exon) to match upstream
    /// RSeQC's `getTranscriptRanges()`. Replicates the upstream bug where
    /// the first transcript per chromosome is dropped.
    pub fn from_bed(bed_path: &str) -> Result<Self> {
        let reader = crate::io::open_reader(bed_path)
            .with_context(|| format!("Failed to open BED file: {}", bed_path))?;

        // Collect all transcripts first, then drop the first per chromosome
        let mut per_chrom: HashMap<String, Vec<(u64, u64, String)>> = HashMap::new();

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

            let tx_name = format!("{}:{}:{}-{}", name, chrom, tx_start, tx_end);
            per_chrom
                .entry(chrom)
                .or_default()
                .push((tx_start, tx_end, tx_name));
        }

        let mut tree = TranscriptTree::default();
        // Replicate the upstream bug: skip the first transcript per chromosome
        for (chrom, transcripts) in &per_chrom {
            for (start, end, name) in transcripts.iter().skip(1) {
                tree.add(chrom, *start, *end, name);
            }
        }

        tree.build();
        Ok(tree)
    }

    /// Add a transcript range to the pending list (must call `build()` after).
    pub fn add(&mut self, chrom: &str, start: u64, end: u64, name: &str) {
        self.pending
            .entry(chrom.to_string())
            .or_default()
            .push((start, end, name.to_string()));
    }

    /// Build the COITrees from accumulated pending intervals.
    /// Must be called once after all `add()` calls and before any `find_overlapping()`.
    pub fn build(&mut self) {
        let pending = std::mem::take(&mut self.pending);
        for (chrom, intervals) in pending {
            let coitree_intervals: Vec<Interval<u32>> = intervals
                .into_iter()
                .map(|(start, end, name)| {
                    let idx = self.names.len() as u32;
                    self.names.push(name);
                    // COITree uses end-inclusive i32 coordinates: [start, end-1]
                    debug_assert!(
                        start <= i32::MAX as u64 && end <= i32::MAX as u64,
                        "Coordinate overflow: TranscriptTree requires coordinates < 2^31"
                    );
                    Interval::new(start as i32, (end.saturating_sub(1)) as i32, idx)
                })
                .collect();
            self.trees.insert(chrom, COITree::new(&coitree_intervals));
        }
    }

    /// Find transcript names that overlap a point. Uses O(log T + hits) interval
    /// tree query instead of the previous O(T) linear scan.
    pub fn find_overlapping(&self, chrom: &str, point: u64) -> HashSet<String> {
        let mut result = HashSet::new();
        if let Some(tree) = self.trees.get(chrom) {
            let p = point as i32;
            tree.query(p, p, |node| {
                fn as_usize(v: impl std::borrow::Borrow<u32>) -> usize {
                    *v.borrow() as usize
                }
                let idx = as_usize(node.metadata);
                if idx < self.names.len() {
                    result.insert(self.names[idx].clone());
                }
            });
        }
        result
    }
}

// (BED12 parsing moved to ExonBitset::from_bed and TranscriptTree::from_bed)

// ============================================================================
// Histogram
// ============================================================================

/// Build histogram bins for inner distances.
///
/// Uses sort + single linear scan (O(n log n + m)) instead of per-bin filtering
/// (O(n × m)) where n = distances.len() and m = number of bins.
pub fn build_histogram(
    distances: &[i64],
    lower_bound: i64,
    upper_bound: i64,
    step: i64,
) -> Vec<(i64, i64, u64)> {
    let mut sorted = distances.to_vec();
    sorted.sort_unstable();

    let mut bins: Vec<(i64, i64, u64)> = Vec::new();
    let mut cursor = 0;
    let mut bin_start = lower_bound;
    while bin_start < upper_bound {
        let bin_end = bin_start + step;
        // Skip values <= bin_start
        while cursor < sorted.len() && sorted[cursor] <= bin_start {
            cursor += 1;
        }
        // Count values in (bin_start, bin_end]
        let start_cursor = cursor;
        while cursor < sorted.len() && sorted[cursor] <= bin_end {
            cursor += 1;
        }
        bins.push((bin_start, bin_end, (cursor - start_cursor) as u64));
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
    fn test_transcript_tree_first_per_chrom_dropped() {
        // When using from_bed / from_genes, the first transcript per chromosome
        // is dropped to replicate the upstream RSeQC bug. But when using add()
        // directly, all transcripts are included.
        let mut tree = TranscriptTree::default();
        tree.add("CHR1", 100, 200, "gene1");
        tree.add("CHR1", 150, 300, "gene2");
        tree.build();

        let genes = tree.find_overlapping("CHR1", 175);
        assert!(genes.contains("gene1")); // Direct add() includes all
        assert!(genes.contains("gene2"));
    }
}
