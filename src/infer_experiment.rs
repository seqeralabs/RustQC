//! Infer library strandedness from RNA-Seq alignments.
//!
//! Reimplementation of RSeQC's `infer_experiment.py`. Samples reads that overlap
//! gene models (BED12) and determines the fraction consistent with each strand
//! protocol.

use anyhow::{Context, Result};
use log::{debug, info, warn};
use rust_htslib::bam::{self, Read as BamRead};
use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

// ============================================================================
// BED12 gene model
// ============================================================================

/// A genomic interval with strand information, parsed from a BED12 file.
#[derive(Debug, Clone)]
struct TranscriptInterval {
    /// 0-based start position (BED start)
    start: u64,
    /// 0-based exclusive end position (BED end)
    end: u64,
    /// Strand: `b'+'` or `b'-'`
    strand: u8,
}

/// Per-chromosome collection of transcript intervals, sorted by start position.
#[derive(Debug, Default)]
pub struct GeneModel {
    /// Map from chromosome name to sorted list of transcript intervals.
    intervals: HashMap<String, Vec<TranscriptInterval>>,
}

impl GeneModel {
    /// Load a BED12 (or BED6+) file into the gene model.
    ///
    /// Extracts columns: chrom (0), start (1), end (2), strand (5).
    pub fn from_bed<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let file = fs::File::open(path)
            .with_context(|| format!("Failed to open BED file: {}", path.display()))?;
        let reader = BufReader::new(file);

        let mut model = GeneModel::default();
        let mut count: u64 = 0;

        for line in reader.lines() {
            let line = line.context("Failed to read BED line")?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') || line.starts_with("track") {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 6 {
                continue; // Need at least 6 columns for strand
            }

            let chrom = fields[0].to_string();
            let start: u64 = fields[1]
                .parse()
                .with_context(|| format!("Invalid BED start: {}", fields[1]))?;
            let end: u64 = fields[2]
                .parse()
                .with_context(|| format!("Invalid BED end: {}", fields[2]))?;
            let strand = match fields[5] {
                "+" => b'+',
                "-" => b'-',
                _ => continue, // Skip unknown strand
            };

            model
                .intervals
                .entry(chrom)
                .or_default()
                .push(TranscriptInterval { start, end, strand });
            count += 1;
        }

        // Sort intervals by start position for binary search
        for intervals in model.intervals.values_mut() {
            intervals.sort_by_key(|iv| iv.start);
        }

        info!("Loaded {} transcripts from BED file", count);
        Ok(model)
    }

    /// Find all strands of transcripts overlapping the query interval [qstart, qend).
    ///
    /// Returns a set of strand bytes found among overlapping transcripts.
    fn find_strands(&self, chrom: &str, qstart: u64, qend: u64) -> Vec<u8> {
        let mut strands = Vec::new();
        if let Some(intervals) = self.intervals.get(chrom) {
            // Binary search to find the first interval that could overlap
            // An interval overlaps if interval.start < qend && interval.end > qstart
            let idx = intervals.partition_point(|iv| iv.start < qend);
            // Check all intervals from 0..idx that could overlap
            // (their start < qend, but we also need end > qstart)
            // We scan backwards from idx since intervals are sorted by start
            // Actually, we need to scan from the beginning because intervals with
            // small start could still overlap if their end > qstart
            for iv in intervals.iter().take(idx) {
                if iv.end > qstart && !strands.contains(&iv.strand) {
                    strands.push(iv.strand);
                }
            }
        }
        strands
    }
}

// ============================================================================
// Strandedness inference
// ============================================================================

/// Result of strandedness inference.
#[derive(Debug)]
pub struct InferExperimentResult {
    /// Number of usable reads sampled.
    pub total_sampled: u64,
    /// Library type: "PairEnd", "SingleEnd", or "Mixture".
    pub library_type: String,
    /// Fraction of reads that failed to determine strandedness.
    pub frac_failed: f64,
    /// Fraction explained by protocol 1 (PE: 1++,1--,2+-,2-+; SE: ++,--)
    pub frac_protocol1: f64,
    /// Fraction explained by protocol 2 (PE: 1+-,1-+,2++,2--; SE: +-,-+)
    pub frac_protocol2: f64,
}

/// Infer library strandedness from a BAM file using a BED12 gene model.
///
/// # Arguments
///
/// * `bam_path` - Path to the BAM/SAM/CRAM file.
/// * `bed_path` - Path to the BED12 gene model file.
/// * `mapq_cut` - Minimum MAPQ score to consider a read (default 30).
/// * `sample_size` - Maximum number of usable reads to sample (default 200000).
/// * `reference` - Optional path to reference FASTA (for CRAM).
///
/// # Returns
///
/// An `InferExperimentResult` with the inferred strandedness fractions.
pub fn infer_experiment<P: AsRef<Path>>(
    bam_path: P,
    model: &GeneModel,
    mapq_cut: u8,
    sample_size: u64,
    reference: Option<&str>,
) -> Result<InferExperimentResult> {
    let bam_path = bam_path.as_ref();

    // Open BAM file
    let mut bam = bam::Reader::from_path(bam_path)
        .with_context(|| format!("Failed to open BAM file: {}", bam_path.display()))?;
    if let Some(ref_path) = reference {
        bam.set_reference(ref_path)
            .with_context(|| format!("Failed to set reference: {}", ref_path))?;
    }

    let header_view = bam.header();
    let tid_to_name: Vec<String> = (0..header_view.target_count())
        .map(|tid| {
            let raw = header_view.tid2name(tid);
            String::from_utf8_lossy(raw).to_string()
        })
        .collect();

    // Paired-end strandness counts: keys like "1++", "1--", "2+-", "2-+", etc.
    let mut p_strandness: HashMap<String, u64> = HashMap::new();
    // Single-end strandness counts: keys like "++", "--", "+-", "-+", etc.
    let mut s_strandness: HashMap<String, u64> = HashMap::new();

    let mut count: u64 = 0;
    let mut record = bam::Record::new();

    while count < sample_size {
        match bam.read(&mut record) {
            Some(Ok(())) => {}
            Some(Err(e)) => {
                warn!("Error reading BAM record: {}", e);
                continue;
            }
            None => break, // EOF
        }

        // Skip QC-failed, duplicate, secondary, unmapped
        let flags = record.flags();
        if record.is_quality_check_failed() {
            continue;
        }
        if record.is_duplicate() {
            continue;
        }
        if record.is_secondary() {
            continue;
        }
        if record.is_unmapped() {
            continue;
        }

        // Skip supplementary alignments
        if flags & 0x800 != 0 {
            continue;
        }

        // MAPQ filter
        if record.mapq() < mapq_cut {
            continue;
        }

        // Get chromosome name
        let tid = record.tid();
        if tid < 0 {
            continue;
        }
        let chrom = &tid_to_name[tid as usize];

        // Get read mapping strand
        let map_strand = if record.is_reverse() { '-' } else { '+' };

        // Get read span on the reference.
        // Note: RSeQC uses `pos + qlen` where pysam's qlen is the query alignment
        // length (M+I+=+X, excluding soft clips). This slightly overestimates the
        // reference span for reads with insertions (I consumes query, not reference)
        // and underestimates for reads with deletions (D consumes reference, not
        // query). We replicate this behavior to match RSeQC output.
        let read_start = record.pos() as u64;
        let qalen: u64 = record
            .cigar()
            .iter()
            .filter_map(|op| {
                use rust_htslib::bam::record::Cigar::*;
                match op {
                    Match(len) | Ins(len) | Equal(len) | Diff(len) => Some(*len as u64),
                    _ => None,
                }
            })
            .sum();
        let read_end = read_start + qalen;

        // Find overlapping gene strands
        let strands = model.find_strands(chrom, read_start, read_end);
        if strands.is_empty() {
            continue; // No overlapping genes — skip
        }

        // Build strand string from gene model
        let strand_str: String = strands
            .iter()
            .map(|&s| s as char)
            .collect::<Vec<char>>()
            .iter()
            .map(|c| c.to_string())
            .collect::<Vec<String>>()
            .join(":");

        if record.is_paired() {
            let read_id = if record.is_first_in_template() {
                "1"
            } else {
                "2"
            };
            let key = format!("{}{}{}", read_id, map_strand, strand_str);
            *p_strandness.entry(key).or_insert(0) += 1;
        } else {
            let key = format!("{}{}", map_strand, strand_str);
            *s_strandness.entry(key).or_insert(0) += 1;
        }

        count += 1;
    }

    if count < 1000 {
        warn!(
            "Too few usable reads ({}). Results may not be reliable.",
            count
        );
    }

    info!("Total {} usable reads were sampled", count);

    // Determine library type and calculate fractions
    let (library_type, frac_protocol1, frac_protocol2, frac_failed) =
        if !p_strandness.is_empty() && s_strandness.is_empty() {
            // Paired-end
            let total: u64 = p_strandness.values().sum();
            let total_f = total as f64;

            let spec1 = (*p_strandness.get("1++").unwrap_or(&0)
                + *p_strandness.get("1--").unwrap_or(&0)
                + *p_strandness.get("2+-").unwrap_or(&0)
                + *p_strandness.get("2-+").unwrap_or(&0)) as f64
                / total_f;

            let spec2 = (*p_strandness.get("1+-").unwrap_or(&0)
                + *p_strandness.get("1-+").unwrap_or(&0)
                + *p_strandness.get("2++").unwrap_or(&0)
                + *p_strandness.get("2--").unwrap_or(&0)) as f64
                / total_f;

            let other = (1.0 - spec1 - spec2).max(0.0);

            ("PairEnd".to_string(), spec1, spec2, other)
        } else if p_strandness.is_empty() && !s_strandness.is_empty() {
            // Single-end
            let total: u64 = s_strandness.values().sum();
            let total_f = total as f64;

            let spec1 = (*s_strandness.get("++").unwrap_or(&0)
                + *s_strandness.get("--").unwrap_or(&0)) as f64
                / total_f;

            let spec2 = (*s_strandness.get("+-").unwrap_or(&0)
                + *s_strandness.get("-+").unwrap_or(&0)) as f64
                / total_f;

            let other = (1.0 - spec1 - spec2).max(0.0);

            ("SingleEnd".to_string(), spec1, spec2, other)
        } else if !p_strandness.is_empty() && !s_strandness.is_empty() {
            // Mixture
            warn!("Data contains a mixture of paired-end and single-end reads.");
            ("Mixture".to_string(), 0.0, 0.0, 0.0)
        } else {
            // No usable reads
            warn!("No usable reads found overlapping gene models.");
            ("Unknown".to_string(), 0.0, 0.0, 0.0)
        };

    // Debug: print all strandness keys
    if !p_strandness.is_empty() {
        debug!("Paired-end strandness counts:");
        for (key, val) in &p_strandness {
            debug!("  {} = {}", key, val);
        }
    }
    if !s_strandness.is_empty() {
        debug!("Single-end strandness counts:");
        for (key, val) in &s_strandness {
            debug!("  {} = {}", key, val);
        }
    }

    Ok(InferExperimentResult {
        total_sampled: count,
        library_type,
        frac_failed,
        frac_protocol1,
        frac_protocol2,
    })
}

// ============================================================================
// Output
// ============================================================================

/// Write infer_experiment results to a file in RSeQC-compatible format.
///
/// # Arguments
///
/// * `result` - The inference result.
/// * `output_path` - Path to write the output file.
pub fn write_infer_experiment<P: AsRef<Path>>(
    result: &InferExperimentResult,
    output_path: P,
) -> Result<()> {
    let output_path = output_path.as_ref();
    let mut writer = fs::File::create(output_path)
        .with_context(|| format!("Failed to create output file: {}", output_path.display()))?;

    match result.library_type.as_str() {
        "PairEnd" => {
            writeln!(writer)?;
            writeln!(writer, "This is PairEnd Data")?;
            writeln!(
                writer,
                "Fraction of reads failed to determine: {:.4}",
                result.frac_failed
            )?;
            writeln!(
                writer,
                "Fraction of reads explained by \"1++,1--,2+-,2-+\": {:.4}",
                result.frac_protocol1
            )?;
            writeln!(
                writer,
                "Fraction of reads explained by \"1+-,1-+,2++,2--\": {:.4}",
                result.frac_protocol2
            )?;
        }
        "SingleEnd" => {
            writeln!(writer)?;
            writeln!(writer, "This is SingleEnd Data")?;
            writeln!(
                writer,
                "Fraction of reads failed to determine: {:.4}",
                result.frac_failed
            )?;
            writeln!(
                writer,
                "Fraction of reads explained by \"++,--\": {:.4}",
                result.frac_protocol1
            )?;
            writeln!(
                writer,
                "Fraction of reads explained by \"+-,-+\": {:.4}",
                result.frac_protocol2
            )?;
        }
        _ => {
            writeln!(writer)?;
            writeln!(writer, "Unknown Data type")?;
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
    fn test_bed_loading() {
        // Test that BED loading works with the small test file if available
        let bed_path = "../RustQC/benchmark/small/chr6.bed";
        if Path::new(bed_path).exists() {
            let model = GeneModel::from_bed(bed_path).unwrap();
            assert!(!model.intervals.is_empty(), "Should have loaded intervals");
            // chr6.bed should have 8443 transcripts
            let total: usize = model.intervals.values().map(|v| v.len()).sum();
            assert!(total > 0, "Should have transcripts");
        }
    }

    #[test]
    fn test_find_strands_no_overlap() {
        let model = GeneModel::default();
        let strands = model.find_strands("chr1", 100, 200);
        assert!(strands.is_empty());
    }

    #[test]
    fn test_find_strands_single_overlap() {
        let mut model = GeneModel::default();
        model.intervals.insert(
            "chr1".to_string(),
            vec![TranscriptInterval {
                start: 100,
                end: 500,
                strand: b'+',
            }],
        );
        let strands = model.find_strands("chr1", 200, 300);
        assert_eq!(strands, vec![b'+']);
    }

    #[test]
    fn test_find_strands_both_strands() {
        let mut model = GeneModel::default();
        model.intervals.insert(
            "chr1".to_string(),
            vec![
                TranscriptInterval {
                    start: 100,
                    end: 500,
                    strand: b'+',
                },
                TranscriptInterval {
                    start: 200,
                    end: 600,
                    strand: b'-',
                },
            ],
        );
        let strands = model.find_strands("chr1", 250, 350);
        assert_eq!(strands.len(), 2);
        assert!(strands.contains(&b'+'));
        assert!(strands.contains(&b'-'));
    }
}
