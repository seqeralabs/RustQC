//! Read duplication rate analysis (reimplementation of RSeQC's read_duplication.py).
//!
//! Computes position-based and sequence-based duplication histograms from
//! a BAM/SAM/CRAM file in a single pass.

use anyhow::{Context, Result};
use log::{debug, info};
use rust_htslib::bam::{self, Read as BamRead};
use std::collections::{BTreeMap, HashMap};
use std::io::Write;
use std::path::Path;
use std::time::Instant;

// ===================================================================
// Result structures
// ===================================================================

/// Histogram of duplication rates.
///
/// Maps occurrence count → number of distinct sequences/positions
/// observed at that count.
pub type DupHistogram = BTreeMap<u64, u64>;

/// Results from read duplication analysis.
#[derive(Debug)]
pub struct ReadDuplicationResult {
    /// Position-based duplication histogram.
    pub pos_histogram: DupHistogram,
    /// Sequence-based duplication histogram.
    pub seq_histogram: DupHistogram,
}

// ===================================================================
// Core algorithm
// ===================================================================

/// Build the position key for a read from its CIGAR alignment.
///
/// Constructs `{chrom}:{start}:{exon1_start}-{exon1_end}:{exon2_start}-{exon2_end}:...`
/// matching RSeQC's `fetch_exon` + position key logic.
///
/// # Arguments
/// * `chrom` - Reference sequence name
/// * `pos` - 0-based leftmost mapping position
/// * `cigar` - CIGAR view from the BAM record
fn build_position_key(chrom: &str, pos: i64, cigar: &bam::record::CigarStringView) -> String {
    use rust_htslib::bam::record::Cigar;

    let mut key = format!("{}:{}:", chrom, pos);
    let mut ref_pos = pos;

    for op in cigar.iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                // M/=/X: alignment match — creates an exon block
                let end = ref_pos + *len as i64;
                key.push_str(&format!("{}-{}:", ref_pos, end));
                ref_pos = end;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                // D/N: deletion or intron skip — advance reference position
                ref_pos += *len as i64;
            }
            Cigar::SoftClip(len) => {
                // S: RSeQC's fetch_exon advances chrom_st for soft clips
                // (arguably a bug, but we replicate for compatibility)
                ref_pos += *len as i64;
            }
            Cigar::Ins(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {
                // I/H/P: do not consume reference
            }
        }
    }

    key
}

/// Perform read duplication analysis on a BAM/SAM/CRAM file.
///
/// Single-pass iteration counting both position-based and sequence-based
/// duplicates. Filters: skip unmapped, QC-fail, and low MAPQ reads.
/// Does NOT skip duplicates or secondary alignments (intentional — we are
/// measuring duplication).
///
/// # Arguments
/// * `bam_path` - Path to the input alignment file
/// * `mapq_cut` - Minimum MAPQ threshold (reads below this are skipped)
/// * `reference` - Optional path to reference FASTA (for CRAM)
pub fn read_duplication(
    bam_path: &str,
    mapq_cut: u8,
    reference: Option<&str>,
) -> Result<ReadDuplicationResult> {
    let start = Instant::now();
    info!("Running read duplication analysis on {}", bam_path);

    // Open BAM file
    let mut bam = if let Some(ref_path) = reference {
        let mut reader = bam::Reader::from_path(bam_path)
            .with_context(|| format!("Failed to open BAM file: {}", bam_path))?;
        reader.set_reference(ref_path)?;
        reader
    } else {
        bam::Reader::from_path(bam_path)
            .with_context(|| format!("Failed to open BAM file: {}", bam_path))?
    };

    // Get reference names for position keys
    let header = bam.header().clone();
    let target_names: Vec<String> = (0..header.target_count())
        .map(|tid| String::from_utf8_lossy(header.tid2name(tid)).to_string())
        .collect();

    // Dictionaries for counting
    let mut seq_dup: HashMap<Vec<u8>, u64> = HashMap::new();
    let mut pos_dup: HashMap<String, u64> = HashMap::new();

    let mut total_processed = 0u64;

    for result in bam.records() {
        let record = result.context("Failed to read BAM record")?;

        // Filter: skip unmapped
        if record.is_unmapped() {
            continue;
        }

        // Filter: skip QC-fail
        if record.is_quality_check_failed() {
            continue;
        }

        // Filter: skip low MAPQ
        if record.mapq() < mapq_cut {
            continue;
        }

        total_processed += 1;

        // --- Sequence-based duplication ---
        // Use the full SEQ field, uppercased
        let seq = record.seq().as_bytes();
        let seq_upper: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
        *seq_dup.entry(seq_upper).or_insert(0) += 1;

        // --- Position-based duplication ---
        let tid = record.tid();
        if tid >= 0 && (tid as usize) < target_names.len() {
            let chrom = &target_names[tid as usize];
            let pos = record.pos();
            let cigar = record.cigar();
            let key = build_position_key(chrom, pos, &cigar);
            *pos_dup.entry(key).or_insert(0) += 1;
        }
    }

    debug!(
        "Read duplication: processed {} reads, {} distinct sequences, {} distinct positions",
        total_processed,
        seq_dup.len(),
        pos_dup.len()
    );

    // Aggregate into histograms
    let mut pos_histogram: DupHistogram = BTreeMap::new();
    for count in pos_dup.values() {
        *pos_histogram.entry(*count).or_insert(0) += 1;
    }

    let mut seq_histogram: DupHistogram = BTreeMap::new();
    for count in seq_dup.values() {
        *seq_histogram.entry(*count).or_insert(0) += 1;
    }

    let elapsed = start.elapsed();
    info!(
        "Read duplication analysis complete in {:.1}s ({} reads processed)",
        elapsed.as_secs_f64(),
        total_processed
    );

    Ok(ReadDuplicationResult {
        pos_histogram,
        seq_histogram,
    })
}

// ===================================================================
// Output
// ===================================================================

/// Write a duplication histogram to a TSV file in RSeQC format.
///
/// Format: `Occurrence\tUniqReadNumber\n`
///
/// # Arguments
/// * `histogram` - The duplication histogram to write
/// * `path` - Output file path
fn write_histogram(histogram: &DupHistogram, path: &Path) -> Result<()> {
    let mut file = std::fs::File::create(path)
        .with_context(|| format!("Failed to create output file: {}", path.display()))?;

    writeln!(file, "Occurrence\tUniqReadNumber")?;
    for (occurrence, uniq_count) in histogram {
        writeln!(file, "{}\t{}", occurrence, uniq_count)?;
    }

    Ok(())
}

/// Write read duplication results to output files.
///
/// Creates two files:
/// - `{stem}.pos.DupRate.xls` — position-based duplication histogram
/// - `{stem}.seq.DupRate.xls` — sequence-based duplication histogram
///
/// # Arguments
/// * `result` - The duplication analysis results
/// * `outdir` - Output directory
/// * `stem` - File name stem (typically BAM file name without extension)
pub fn write_read_duplication(
    result: &ReadDuplicationResult,
    outdir: &Path,
    stem: &str,
) -> Result<()> {
    let pos_path = outdir.join(format!("{}.pos.DupRate.xls", stem));
    let seq_path = outdir.join(format!("{}.seq.DupRate.xls", stem));

    write_histogram(&result.pos_histogram, &pos_path).with_context(|| {
        format!(
            "Failed to write position duplication file: {}",
            pos_path.display()
        )
    })?;
    info!(
        "Wrote position-based duplication rates to {}",
        pos_path.display()
    );

    write_histogram(&result.seq_histogram, &seq_path).with_context(|| {
        format!(
            "Failed to write sequence duplication file: {}",
            seq_path.display()
        )
    })?;
    info!(
        "Wrote sequence-based duplication rates to {}",
        seq_path.display()
    );

    Ok(())
}

// ===================================================================
// Tests
// ===================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_position_key_simple() {
        // A simple case: 50M read at chr1:1000
        use rust_htslib::bam::record::Cigar;
        use rust_htslib::bam::record::{CigarString, CigarStringView};

        let cigar_ops = vec![Cigar::Match(50)];
        let cigar_string = CigarString(cigar_ops);
        let cigar_view = CigarStringView::new(cigar_string, 1000);

        let key = build_position_key("chr1", 1000, &cigar_view);
        assert_eq!(key, "chr1:1000:1000-1050:");
    }

    #[test]
    fn test_build_position_key_spliced() {
        // Spliced read: 10M500N20M at chr1:1000
        use rust_htslib::bam::record::Cigar;
        use rust_htslib::bam::record::{CigarString, CigarStringView};

        let cigar_ops = vec![Cigar::Match(10), Cigar::RefSkip(500), Cigar::Match(20)];
        let cigar_string = CigarString(cigar_ops);
        let cigar_view = CigarStringView::new(cigar_string, 1000);

        let key = build_position_key("chr1", 1000, &cigar_view);
        assert_eq!(key, "chr1:1000:1000-1010:1510-1530:");
    }

    #[test]
    fn test_read_duplication_small() {
        // Integration test with the small test BAM
        let bam_path = "tests/data/test.bam";
        if !Path::new(bam_path).exists() {
            return; // Skip if test data not available
        }

        let result = read_duplication(bam_path, 30, None).unwrap();

        // Basic sanity: both histograms should be non-empty
        assert!(
            !result.pos_histogram.is_empty(),
            "Position histogram should not be empty"
        );
        assert!(
            !result.seq_histogram.is_empty(),
            "Sequence histogram should not be empty"
        );

        // The highest occurrence counts should be reasonable
        let max_pos_occ = *result.pos_histogram.keys().last().unwrap();
        assert!(max_pos_occ >= 1, "Should have at least occurrence=1");

        let max_seq_occ = *result.seq_histogram.keys().last().unwrap();
        assert!(max_seq_occ >= 1, "Should have at least occurrence=1");
    }
}
