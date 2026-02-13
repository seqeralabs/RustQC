//! Read distribution across genomic features.
//!
//! Reimplementation of RSeQC's `read_distribution.py`: classifies BAM read tags
//! into CDS exons, 5'/3' UTRs, introns, and intergenic regions using a BED12
//! gene model. Tags (CIGAR M-blocks) are classified by midpoint with priority
//! CDS > UTR > Intron > Intergenic.

use std::collections::{BTreeMap, HashMap};
use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use log::{debug, info};
use rust_htslib::bam::{self, Read as BamRead};

// ===================================================================
// Region types
// ===================================================================

/// A genomic interval: [start, end) on a chromosome (0-based, half-open).
#[derive(Debug, Clone, Copy)]
struct Interval {
    start: u64,
    end: u64,
}

/// Merged intervals for a single chromosome, sorted and non-overlapping.
#[derive(Debug, Clone, Default)]
struct ChromIntervals {
    intervals: Vec<Interval>,
}

impl ChromIntervals {
    /// Add a raw interval (may overlap existing ones).
    fn add(&mut self, start: u64, end: u64) {
        if start < end {
            self.intervals.push(Interval { start, end });
        }
    }

    /// Sort and merge overlapping intervals.
    fn merge(&mut self) {
        if self.intervals.is_empty() {
            return;
        }
        self.intervals.sort_by_key(|i| (i.start, i.end));
        let mut merged = Vec::with_capacity(self.intervals.len());
        let mut current = self.intervals[0];
        for iv in &self.intervals[1..] {
            if iv.start <= current.end {
                current.end = current.end.max(iv.end);
            } else {
                merged.push(current);
                current = *iv;
            }
        }
        merged.push(current);
        self.intervals = merged;
    }

    /// Subtract another set of intervals from this one.
    /// Both must be pre-merged and sorted.
    fn subtract(&mut self, other: &ChromIntervals) {
        if self.intervals.is_empty() || other.intervals.is_empty() {
            return;
        }
        let mut result = Vec::new();
        let mut j = 0;
        for iv in &self.intervals {
            let mut start = iv.start;
            let end = iv.end;
            while j < other.intervals.len() && other.intervals[j].end <= start {
                j += 1;
            }
            let mut k = j;
            while k < other.intervals.len() && other.intervals[k].start < end {
                let sub = &other.intervals[k];
                if sub.start > start {
                    result.push(Interval {
                        start,
                        end: sub.start.min(end),
                    });
                }
                start = sub.end;
                k += 1;
            }
            if start < end {
                result.push(Interval { start, end });
            }
        }
        self.intervals = result;
    }

    /// Total bases covered by these intervals.
    fn total_bases(&self) -> u64 {
        self.intervals.iter().map(|i| i.end - i.start).sum()
    }

    /// Check if a point falls strictly inside any interval (binary search).
    /// Uses strict containment: `start < point < end`, matching RSeQC's
    /// bx-python `Intersecter.find(mid, mid)` which checks
    /// `(self.start < end) and (self.end > start)` with start==end==mid,
    /// i.e. `self.start < mid and self.end > mid`.
    fn contains(&self, point: u64) -> bool {
        self.intervals
            .binary_search_by(|iv| {
                if point <= iv.start {
                    std::cmp::Ordering::Greater
                } else if point >= iv.end {
                    std::cmp::Ordering::Less
                } else {
                    std::cmp::Ordering::Equal
                }
            })
            .is_ok()
    }
}

/// Per-chromosome region sets for all feature types.
#[derive(Debug, Default)]
struct RegionSets {
    cds_exon: HashMap<String, ChromIntervals>,
    utr_5: HashMap<String, ChromIntervals>,
    utr_3: HashMap<String, ChromIntervals>,
    intron: HashMap<String, ChromIntervals>,
    tss_up_1kb: HashMap<String, ChromIntervals>,
    tss_up_5kb: HashMap<String, ChromIntervals>,
    tss_up_10kb: HashMap<String, ChromIntervals>,
    tes_down_1kb: HashMap<String, ChromIntervals>,
    tes_down_5kb: HashMap<String, ChromIntervals>,
    tes_down_10kb: HashMap<String, ChromIntervals>,
}

// ===================================================================
// BED12 parsing and region extraction
// ===================================================================

/// A parsed BED12 transcript record.
#[derive(Debug)]
struct Bed12Record {
    chrom: String,
    tx_start: u64,
    tx_end: u64,
    strand: char,
    cds_start: u64,
    cds_end: u64,
    exon_starts: Vec<u64>,
    exon_ends: Vec<u64>,
}

impl Bed12Record {
    /// Parse a BED12 line. Returns None if malformed.
    fn parse(line: &str) -> Option<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            return None;
        }
        let chrom = fields[0].to_uppercase();
        let tx_start: u64 = fields[1].parse().ok()?;
        let tx_end: u64 = fields[2].parse().ok()?;
        let strand = fields[5].chars().next().unwrap_or('+');
        let cds_start: u64 = fields[6].parse().ok()?;
        let cds_end: u64 = fields[7].parse().ok()?;
        let block_count: usize = fields[9].parse().ok()?;
        let block_sizes: Vec<u64> = fields[10]
            .trim_end_matches(',')
            .split(',')
            .filter_map(|s| s.parse().ok())
            .collect();
        let block_starts: Vec<u64> = fields[11]
            .trim_end_matches(',')
            .split(',')
            .filter_map(|s| s.parse().ok())
            .collect();

        if block_sizes.len() < block_count || block_starts.len() < block_count {
            return None;
        }

        let exon_starts: Vec<u64> = block_starts
            .iter()
            .take(block_count)
            .map(|s| tx_start + s)
            .collect();
        let exon_ends: Vec<u64> = exon_starts
            .iter()
            .zip(block_sizes.iter().take(block_count))
            .map(|(s, sz)| s + sz)
            .collect();

        Some(Bed12Record {
            chrom,
            tx_start,
            tx_end,
            strand,
            cds_start,
            cds_end,
            exon_starts,
            exon_ends,
        })
    }

    /// Get CDS exon intervals (intersection of exon blocks with [cds_start, cds_end)).
    fn get_cds_exons(&self) -> Vec<(u64, u64)> {
        let mut result = Vec::new();
        for (start, end) in self.exon_starts.iter().zip(self.exon_ends.iter()) {
            if *end <= self.cds_start || *start >= self.cds_end {
                continue;
            }
            let s = (*start).max(self.cds_start);
            let e = (*end).min(self.cds_end);
            if s < e {
                result.push((s, e));
            }
        }
        result
    }

    /// Get 5' UTR intervals.
    fn get_utr_5(&self) -> Vec<(u64, u64)> {
        let mut result = Vec::new();
        for (start, end) in self.exon_starts.iter().zip(self.exon_ends.iter()) {
            match self.strand {
                '+' => {
                    if *start < self.cds_start {
                        result.push((*start, (*end).min(self.cds_start)));
                    }
                }
                '-' => {
                    if *end > self.cds_end {
                        result.push(((*start).max(self.cds_end), *end));
                    }
                }
                _ => {}
            }
        }
        result
    }

    /// Get 3' UTR intervals.
    fn get_utr_3(&self) -> Vec<(u64, u64)> {
        let mut result = Vec::new();
        for (start, end) in self.exon_starts.iter().zip(self.exon_ends.iter()) {
            match self.strand {
                '+' => {
                    if *end > self.cds_end {
                        result.push(((*start).max(self.cds_end), *end));
                    }
                }
                '-' => {
                    if *start < self.cds_start {
                        result.push((*start, (*end).min(self.cds_start)));
                    }
                }
                _ => {}
            }
        }
        result
    }

    /// Get intron intervals (gaps between exon blocks).
    fn get_introns(&self) -> Vec<(u64, u64)> {
        let mut result = Vec::new();
        for i in 0..self.exon_starts.len().saturating_sub(1) {
            let start = self.exon_ends[i];
            let end = self.exon_starts[i + 1];
            if start < end {
                result.push((start, end));
            }
        }
        result
    }

    /// Get TSS upstream region.
    fn get_tss_upstream(&self, size: u64) -> (u64, u64) {
        match self.strand {
            '-' => (self.tx_end, self.tx_end + size),
            _ => (self.tx_start.saturating_sub(size), self.tx_start),
        }
    }

    /// Get TES downstream region.
    fn get_tes_downstream(&self, size: u64) -> (u64, u64) {
        match self.strand {
            '-' => (self.tx_start.saturating_sub(size), self.tx_start),
            _ => (self.tx_end, self.tx_end + size),
        }
    }
}

// ===================================================================
// Build region sets from BED12 file
// ===================================================================

/// Build all region sets from a BED12 gene model file.
fn build_regions(bed_path: &str) -> Result<RegionSets> {
    let content = std::fs::read_to_string(bed_path)
        .with_context(|| format!("Failed to read BED file: {}", bed_path))?;

    let mut regions = RegionSets::default();

    for line in content.lines() {
        if line.starts_with('#') || line.starts_with("track") || line.is_empty() {
            continue;
        }
        let rec = match Bed12Record::parse(line) {
            Some(r) => r,
            None => continue,
        };

        // CDS exons
        for (s, e) in rec.get_cds_exons() {
            regions
                .cds_exon
                .entry(rec.chrom.clone())
                .or_default()
                .add(s, e);
        }

        // 5' UTR
        for (s, e) in rec.get_utr_5() {
            regions
                .utr_5
                .entry(rec.chrom.clone())
                .or_default()
                .add(s, e);
        }

        // 3' UTR
        for (s, e) in rec.get_utr_3() {
            regions
                .utr_3
                .entry(rec.chrom.clone())
                .or_default()
                .add(s, e);
        }

        // Introns
        for (s, e) in rec.get_introns() {
            regions
                .intron
                .entry(rec.chrom.clone())
                .or_default()
                .add(s, e);
        }

        // TSS upstream regions
        for (size, map) in [
            (1000u64, &mut regions.tss_up_1kb),
            (5000, &mut regions.tss_up_5kb),
            (10000, &mut regions.tss_up_10kb),
        ] {
            let (s, e) = rec.get_tss_upstream(size);
            map.entry(rec.chrom.clone()).or_default().add(s, e);
        }

        // TES downstream regions
        for (size, map) in [
            (1000u64, &mut regions.tes_down_1kb),
            (5000, &mut regions.tes_down_5kb),
            (10000, &mut regions.tes_down_10kb),
        ] {
            let (s, e) = rec.get_tes_downstream(size);
            map.entry(rec.chrom.clone()).or_default().add(s, e);
        }
    }

    // Merge all regions
    for map in [
        &mut regions.cds_exon,
        &mut regions.utr_5,
        &mut regions.utr_3,
        &mut regions.intron,
        &mut regions.tss_up_1kb,
        &mut regions.tss_up_5kb,
        &mut regions.tss_up_10kb,
        &mut regions.tes_down_1kb,
        &mut regions.tes_down_5kb,
        &mut regions.tes_down_10kb,
    ] {
        for intervals in map.values_mut() {
            intervals.merge();
        }
    }

    // Priority-based subtraction to make regions mutually exclusive
    // utr_5 -= cds_exon
    subtract_regions(&mut regions.utr_5, &regions.cds_exon);
    // utr_3 -= cds_exon
    subtract_regions(&mut regions.utr_3, &regions.cds_exon);
    // intron -= cds_exon, utr_5, utr_3
    subtract_regions(&mut regions.intron, &regions.cds_exon);
    subtract_regions(&mut regions.intron, &regions.utr_5);
    subtract_regions(&mut regions.intron, &regions.utr_3);
    // intergenic -= cds_exon, utr_5, utr_3, intron
    for intergenic in [
        &mut regions.tss_up_1kb,
        &mut regions.tss_up_5kb,
        &mut regions.tss_up_10kb,
        &mut regions.tes_down_1kb,
        &mut regions.tes_down_5kb,
        &mut regions.tes_down_10kb,
    ] {
        subtract_regions(intergenic, &regions.cds_exon);
        subtract_regions(intergenic, &regions.utr_5);
        subtract_regions(intergenic, &regions.utr_3);
        subtract_regions(intergenic, &regions.intron);
    }

    Ok(regions)
}

/// Subtract `sub` regions from `target` regions for all chromosomes.
fn subtract_regions(
    target: &mut HashMap<String, ChromIntervals>,
    sub: &HashMap<String, ChromIntervals>,
) {
    for (chrom, intervals) in target.iter_mut() {
        if let Some(sub_intervals) = sub.get(chrom) {
            intervals.subtract(sub_intervals);
        }
    }
}

/// Sum total bases across all chromosomes for a region set.
fn sum_bases(map: &HashMap<String, ChromIntervals>) -> u64 {
    map.values().map(|ci| ci.total_bases()).sum()
}

/// Check if a point on a chromosome falls within a region set.
fn point_in_region(map: &HashMap<String, ChromIntervals>, chrom: &str, point: u64) -> bool {
    map.get(chrom).map(|ci| ci.contains(point)).unwrap_or(false)
}

// ===================================================================
// BAM processing
// ===================================================================

/// Extract aligned blocks (M-operations) from CIGAR, returning (start, end) pairs.
/// Matches RSeQC's `fetch_exon()` behavior exactly:
/// - Only M (op 0) creates exon blocks
/// - D (op 2) and N (op 3) advance position without creating blocks
/// - S (op 4) advances position (bug in RSeQC, but replicated for compatibility)
/// - I (op 1), H (op 5), P (op 6) don't advance position
/// - = (op 7) and X (op 8) fall through to else branch in RSeQC's Python code,
///   which means they neither create blocks nor advance position (a known bug)
fn fetch_exon_blocks(record: &bam::Record) -> Vec<(u64, u64)> {
    let mut blocks = Vec::new();
    let mut pos = record.pos() as u64;
    for op in record.cigar().iter() {
        match op {
            rust_htslib::bam::record::Cigar::Match(len) => {
                let len = *len as u64;
                blocks.push((pos, pos + len));
                pos += len;
            }
            rust_htslib::bam::record::Cigar::Del(len)
            | rust_htslib::bam::record::Cigar::RefSkip(len) => {
                pos += *len as u64;
            }
            rust_htslib::bam::record::Cigar::SoftClip(len) => {
                // RSeQC's fetch_exon() advances position for soft clips.
                // This is arguably a bug, but we replicate the behavior.
                pos += *len as u64;
            }
            // RSeQC's fetch_exon() ignores = (SeqMatch), X (Diff), I (Ins),
            // H (HardClip), P (Pad) — none advance position or create blocks.
            _ => {}
        }
    }
    blocks
}

/// Result of read distribution analysis.
#[derive(Debug)]
pub struct ReadDistributionResult {
    /// Total reads processed (after filtering).
    pub total_reads: u64,
    /// Total tags (CIGAR M-blocks) across all reads.
    pub total_tags: u64,
    /// Region statistics: (name, total_bases, tag_count).
    pub regions: Vec<(String, u64, u64)>,
    /// Number of tags not assigned to any region.
    pub unassigned_tags: u64,
}

/// Run read distribution analysis on a BAM file.
///
/// # Arguments
/// * `bam_path` - Path to the BAM/CRAM file
/// * `bed_path` - Path to the BED12 reference gene model
/// * `reference` - Optional path to reference FASTA (required for CRAM)
pub fn read_distribution(
    bam_path: &str,
    bed_path: &str,
    reference: Option<&str>,
) -> Result<ReadDistributionResult> {
    info!("Building gene model from BED file: {}", bed_path);
    let regions = build_regions(bed_path)?;

    // Compute total bases for each region
    let cds_bases = sum_bases(&regions.cds_exon);
    let utr5_bases = sum_bases(&regions.utr_5);
    let utr3_bases = sum_bases(&regions.utr_3);
    let intron_bases = sum_bases(&regions.intron);
    let tss_1k_bases = sum_bases(&regions.tss_up_1kb);
    let tss_5k_bases = sum_bases(&regions.tss_up_5kb);
    let tss_10k_bases = sum_bases(&regions.tss_up_10kb);
    let tes_1k_bases = sum_bases(&regions.tes_down_1kb);
    let tes_5k_bases = sum_bases(&regions.tes_down_5kb);
    let tes_10k_bases = sum_bases(&regions.tes_down_10kb);

    debug!(
        "Region bases — CDS: {}, 5'UTR: {}, 3'UTR: {}, Intron: {}",
        cds_bases, utr5_bases, utr3_bases, intron_bases
    );

    // Iterate BAM and classify tags
    let mut bam = if let Some(ref_path) = reference {
        let mut r = bam::Reader::from_path(bam_path)
            .with_context(|| format!("Failed to open BAM file: {}", bam_path))?;
        r.set_reference(ref_path)
            .with_context(|| format!("Failed to set reference: {}", ref_path))?;
        r
    } else {
        bam::Reader::from_path(bam_path)
            .with_context(|| format!("Failed to open BAM file: {}", bam_path))?
    };

    // Build tid -> uppercase chrom name mapping
    let header = bam.header().clone();
    let tid_to_chrom: BTreeMap<u32, String> = (0..header.target_count())
        .filter_map(|tid| {
            let name = std::str::from_utf8(header.tid2name(tid)).ok()?;
            Some((tid, name.to_uppercase()))
        })
        .collect();

    let mut total_reads: u64 = 0;
    let mut total_tags: u64 = 0;
    let mut cds_tags: u64 = 0;
    let mut utr5_tags: u64 = 0;
    let mut utr3_tags: u64 = 0;
    let mut intron_tags: u64 = 0;
    let mut tss_1k_tags: u64 = 0;
    let mut tss_5k_tags: u64 = 0;
    let mut tss_10k_tags: u64 = 0;
    let mut tes_1k_tags: u64 = 0;
    let mut tes_5k_tags: u64 = 0;
    let mut tes_10k_tags: u64 = 0;
    let mut unassigned: u64 = 0;

    let mut record = bam::Record::new();
    while let Some(result) = bam.read(&mut record) {
        result.with_context(|| "Error reading BAM record")?;

        let flags = record.flags();
        // Skip: qc_fail, duplicate, secondary, unmapped
        if flags & 0x200 != 0 {
            continue;
        }
        if flags & 0x400 != 0 {
            continue;
        }
        if flags & 0x100 != 0 {
            continue;
        }
        if record.is_unmapped() {
            continue;
        }

        let tid = record.tid();
        if tid < 0 {
            continue;
        }
        let chrom = match tid_to_chrom.get(&(tid as u32)) {
            Some(c) => c.as_str(),
            None => continue,
        };

        let blocks = fetch_exon_blocks(&record);
        total_reads += 1;
        total_tags += blocks.len() as u64;

        for (block_start, block_end) in &blocks {
            let mid = block_start + (block_end - block_start) / 2;

            // Priority classification cascade
            if point_in_region(&regions.cds_exon, chrom, mid) {
                cds_tags += 1;
            } else {
                let in_utr5 = point_in_region(&regions.utr_5, chrom, mid);
                let in_utr3 = point_in_region(&regions.utr_3, chrom, mid);

                if in_utr5 && !in_utr3 {
                    utr5_tags += 1;
                } else if in_utr3 && !in_utr5 {
                    utr3_tags += 1;
                } else if in_utr5 && in_utr3 {
                    // Both UTRs — unassigned
                    unassigned += 1;
                } else if point_in_region(&regions.intron, chrom, mid) {
                    intron_tags += 1;
                } else {
                    // Intergenic classification
                    let in_tss_10k = point_in_region(&regions.tss_up_10kb, chrom, mid);
                    let in_tes_10k = point_in_region(&regions.tes_down_10kb, chrom, mid);

                    if in_tss_10k && in_tes_10k {
                        // Ambiguous intergenic — unassigned
                        unassigned += 1;
                    } else if in_tss_10k {
                        // Cumulative TSS counting
                        tss_10k_tags += 1;
                        if point_in_region(&regions.tss_up_5kb, chrom, mid) {
                            tss_5k_tags += 1;
                            if point_in_region(&regions.tss_up_1kb, chrom, mid) {
                                tss_1k_tags += 1;
                            }
                        }
                    } else if in_tes_10k {
                        // Cumulative TES counting
                        tes_10k_tags += 1;
                        if point_in_region(&regions.tes_down_5kb, chrom, mid) {
                            tes_5k_tags += 1;
                            if point_in_region(&regions.tes_down_1kb, chrom, mid) {
                                tes_1k_tags += 1;
                            }
                        }
                    } else {
                        unassigned += 1;
                    }
                }
            }
        }
    }

    let result_regions = vec![
        ("CDS_Exons".to_string(), cds_bases, cds_tags),
        ("5'UTR_Exons".to_string(), utr5_bases, utr5_tags),
        ("3'UTR_Exons".to_string(), utr3_bases, utr3_tags),
        ("Introns".to_string(), intron_bases, intron_tags),
        ("TSS_up_1kb".to_string(), tss_1k_bases, tss_1k_tags),
        ("TSS_up_5kb".to_string(), tss_5k_bases, tss_5k_tags),
        ("TSS_up_10kb".to_string(), tss_10k_bases, tss_10k_tags),
        ("TES_down_1kb".to_string(), tes_1k_bases, tes_1k_tags),
        ("TES_down_5kb".to_string(), tes_5k_bases, tes_5k_tags),
        ("TES_down_10kb".to_string(), tes_10k_bases, tes_10k_tags),
    ];

    Ok(ReadDistributionResult {
        total_reads,
        total_tags,
        regions: result_regions,
        unassigned_tags: unassigned,
    })
}

/// Write read distribution results to a file in RSeQC-compatible format.
pub fn write_read_distribution(result: &ReadDistributionResult, output_path: &Path) -> Result<()> {
    let mut writer = std::io::BufWriter::new(
        std::fs::File::create(output_path)
            .with_context(|| format!("Failed to create output file: {}", output_path.display()))?,
    );

    let assigned = result.total_tags - result.unassigned_tags;

    // Header stats
    writeln!(writer, "{:<30}{}", "Total Reads", result.total_reads)?;
    writeln!(writer, "{:<30}{}", "Total Tags", result.total_tags)?;
    writeln!(writer, "{:<30}{}", "Total Assigned Tags", assigned)?;

    // Separator
    writeln!(writer, "{}", "=".repeat(69))?;

    // Table header
    writeln!(
        writer,
        "{:<20}{:<20}{:<20}{:<20}",
        "Group", "Total_bases", "Tag_count", "Tags/Kb"
    )?;

    // Region rows
    for (name, bases, tags) in &result.regions {
        let tags_per_kb = *tags as f64 * 1000.0 / (*bases as f64 + 1.0);
        writeln!(
            writer,
            "{:<20}{:<20}{:<20}{:<18.2}",
            name, bases, tags, tags_per_kb
        )?;
    }

    // Footer separator
    writeln!(writer, "{}", "=".repeat(69))?;

    writer.flush()?;
    Ok(())
}

// ===================================================================
// Tests
// ===================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chrom_intervals_merge() {
        let mut ci = ChromIntervals::default();
        ci.add(10, 20);
        ci.add(15, 25);
        ci.add(30, 40);
        ci.merge();
        assert_eq!(ci.intervals.len(), 2);
        assert_eq!(ci.intervals[0].start, 10);
        assert_eq!(ci.intervals[0].end, 25);
        assert_eq!(ci.intervals[1].start, 30);
        assert_eq!(ci.intervals[1].end, 40);
    }

    #[test]
    fn test_chrom_intervals_subtract() {
        let mut a = ChromIntervals::default();
        a.add(10, 50);
        a.merge();

        let mut b = ChromIntervals::default();
        b.add(20, 30);
        b.merge();

        a.subtract(&b);
        assert_eq!(a.intervals.len(), 2);
        assert_eq!(a.intervals[0].start, 10);
        assert_eq!(a.intervals[0].end, 20);
        assert_eq!(a.intervals[1].start, 30);
        assert_eq!(a.intervals[1].end, 50);
    }

    #[test]
    fn test_chrom_intervals_contains() {
        let mut ci = ChromIntervals::default();
        ci.add(10, 20);
        ci.add(30, 40);
        ci.merge();
        // RSeQC uses strict containment: point > start && point < end
        // (matching bx-python IntervalTree find(mid, mid) semantics)
        assert!(!ci.contains(10)); // boundary: not contained
        assert!(ci.contains(11)); // inside [10, 20)
        assert!(ci.contains(19)); // inside [10, 20)
        assert!(!ci.contains(20)); // boundary: not contained
        assert!(ci.contains(35)); // inside [30, 40)
        assert!(!ci.contains(5)); // outside
        assert!(!ci.contains(25)); // outside (gap between intervals)
        assert!(!ci.contains(30)); // boundary: not contained
        assert!(ci.contains(31)); // inside [30, 40)
        assert!(!ci.contains(40)); // boundary: not contained
    }

    #[test]
    fn test_bed12_cds_exons() {
        // Simple: one exon, fully CDS
        let line = "chr1\t100\t200\tgene1\t0\t+\t100\t200\t0\t1\t100,\t0,";
        let rec = Bed12Record::parse(line).unwrap();
        let cds = rec.get_cds_exons();
        assert_eq!(cds, vec![(100, 200)]);
    }

    #[test]
    fn test_bed12_utr() {
        // Two exons: first is 5'UTR, second is CDS on + strand
        let line = "chr1\t100\t300\tgene1\t0\t+\t200\t300\t0\t2\t50,50,\t0,150,";
        let rec = Bed12Record::parse(line).unwrap();
        let utr5 = rec.get_utr_5();
        assert_eq!(utr5, vec![(100, 150)]); // first exon: 100-150, all before CDS at 200
        let utr3 = rec.get_utr_3();
        assert!(utr3.is_empty()); // no exon portions after CDS end (300)
    }

    #[test]
    fn test_total_bases() {
        let mut ci = ChromIntervals::default();
        ci.add(10, 20);
        ci.add(30, 50);
        ci.merge();
        assert_eq!(ci.total_bases(), 30); // 10 + 20
    }
}
