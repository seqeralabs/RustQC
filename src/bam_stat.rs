//! Basic BAM alignment statistics.
//!
//! Reimplements RSeQC's `bam_stat.py`: a single-pass BAM reader that collects
//! fundamental alignment metrics (total records, QC-failed, duplicates, mapping
//! quality distribution, splice reads, proper pairs, etc.).

use anyhow::{Context, Result};
use log::{debug, info};
use rust_htslib::bam::{self, Read as BamRead};
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;

// ============================================================================
// BAM flag constants (matching SAM spec)
// ============================================================================

const BAM_FPAIRED: u16 = 0x1;
const BAM_FPROPER_PAIR: u16 = 0x2;
const BAM_FUNMAP: u16 = 0x4;
const BAM_FREVERSE: u16 = 0x10;
const BAM_FREAD1: u16 = 0x40;
const BAM_FREAD2: u16 = 0x80;
const BAM_FSECONDARY: u16 = 0x100;
const BAM_FQCFAIL: u16 = 0x200;
const BAM_FDUP: u16 = 0x400;

// ============================================================================
// Result struct
// ============================================================================

/// Holds all statistics collected during the BAM pass.
#[derive(Debug, Default)]
pub struct BamStatResult {
    /// Total number of records in the BAM file.
    pub total_records: u64,
    /// Number of records that failed QC (flag 0x200).
    pub qc_failed: u64,
    /// Number of PCR/optical duplicate records (flag 0x400).
    pub duplicates: u64,
    /// Number of non-primary (secondary) alignment records (flag 0x100).
    /// Note: supplementary alignments (0x800) are NOT counted here,
    /// matching RSeQC's behavior.
    pub non_primary: u64,
    /// Number of unmapped reads (flag 0x4).
    pub unmapped: u64,
    /// Mapped reads with MAPQ < cutoff (non-unique/multimapped).
    pub non_unique: u64,
    /// Mapped reads with MAPQ >= cutoff (uniquely mapped).
    pub unique: u64,
    /// Among unique reads: number that are read1 in a pair.
    pub read_1: u64,
    /// Among unique reads: number that are read2 in a pair.
    pub read_2: u64,
    /// Among unique reads: number mapping to forward strand.
    pub forward: u64,
    /// Among unique reads: number mapping to reverse strand.
    pub reverse: u64,
    /// Among unique reads: number with at least one splice junction (CIGAR N).
    pub splice: u64,
    /// Among unique reads: number without any splice junctions.
    pub non_splice: u64,
    /// Among unique reads: number in proper pairs (flag 0x2).
    pub proper_pairs: u64,
    /// Among proper-paired unique reads: number where mates map to different chromosomes.
    pub proper_pair_diff_chrom: u64,
    /// MAPQ distribution for all primary, non-QC-fail, non-dup, mapped reads.
    pub mapq_distribution: BTreeMap<u8, u64>,
}

// ============================================================================
// Core analysis function
// ============================================================================

/// Compute alignment statistics from a BAM/SAM/CRAM file in a single pass.
///
/// The counting logic mirrors RSeQC's `bam_stat.py` exactly:
/// 1. Every record increments `total_records`.
/// 2. QC-failed records (0x200) are counted and skipped.
/// 3. Duplicate records (0x400) are counted and skipped.
/// 4. Secondary alignments (0x100) are counted and skipped.
///    Supplementary alignments (0x800) are NOT counted as non-primary.
/// 5. Unmapped reads (0x4) are counted and skipped.
/// 6. Remaining mapped primary reads are split by MAPQ threshold.
/// 7. Reads with MAPQ >= cutoff are further classified by strand,
///    read number, splice status, and pairing.
///
/// # Arguments
/// * `bam_path` - Path to the BAM/SAM/CRAM file
/// * `mapq_cut` - MAPQ cutoff for unique/non-unique classification
/// * `reference` - Optional path to reference FASTA (for CRAM)
///
/// # Returns
/// A `BamStatResult` with all collected statistics.
pub fn bam_stat(bam_path: &str, mapq_cut: u8, reference: Option<&str>) -> Result<BamStatResult> {
    info!("Computing alignment statistics for {}", bam_path);

    let mut reader = if let Some(ref_path) = reference {
        let mut r = bam::Reader::from_path(bam_path)
            .with_context(|| format!("Failed to open BAM file: {}", bam_path))?;
        r.set_reference(ref_path)
            .with_context(|| format!("Failed to set reference: {}", ref_path))?;
        r
    } else {
        bam::Reader::from_path(bam_path)
            .with_context(|| format!("Failed to open BAM file: {}", bam_path))?
    };

    let mut result = BamStatResult::default();
    let mut record = bam::Record::new();

    while let Some(r) = reader.read(&mut record) {
        r.with_context(|| format!("Error reading record from {}", bam_path))?;

        let flags = record.flags();
        result.total_records += 1;

        // 1. QC-failed
        if flags & BAM_FQCFAIL != 0 {
            result.qc_failed += 1;
            continue;
        }

        // 2. Duplicate
        if flags & BAM_FDUP != 0 {
            result.duplicates += 1;
            continue;
        }

        // 3. Secondary (non-primary) — NOT supplementary
        if flags & BAM_FSECONDARY != 0 {
            result.non_primary += 1;
            continue;
        }

        // 4. Unmapped
        if flags & BAM_FUNMAP != 0 {
            result.unmapped += 1;
            continue;
        }

        // At this point: mapped, primary, non-QC-fail, non-duplicate
        let mapq = record.mapq();

        // Record MAPQ distribution
        *result.mapq_distribution.entry(mapq).or_insert(0) += 1;

        // 5. MAPQ classification
        if mapq < mapq_cut {
            result.non_unique += 1;
            continue;
        }

        // Uniquely mapped (MAPQ >= cutoff)
        result.unique += 1;

        // Read1 / Read2
        if flags & BAM_FREAD1 != 0 {
            result.read_1 += 1;
        }
        if flags & BAM_FREAD2 != 0 {
            result.read_2 += 1;
        }

        // Strand
        if flags & BAM_FREVERSE != 0 {
            result.reverse += 1;
        } else {
            result.forward += 1;
        }

        // Splice detection: look for 'N' (reference skip) in CIGAR
        let has_splice = record
            .cigar()
            .iter()
            .any(|op| matches!(op, rust_htslib::bam::record::Cigar::RefSkip(_)));
        if has_splice {
            result.splice += 1;
        } else {
            result.non_splice += 1;
        }

        // Proper pair analysis
        if flags & BAM_FPAIRED != 0 && flags & BAM_FPROPER_PAIR != 0 {
            result.proper_pairs += 1;
            // Check if mate maps to different chromosome
            if record.tid() != record.mtid() {
                result.proper_pair_diff_chrom += 1;
            }
        }
    }

    debug!(
        "bam_stat complete: {} total records, {} unique",
        result.total_records, result.unique
    );

    Ok(result)
}

// ============================================================================
// Output formatting
// ============================================================================

/// Write bam_stat results in the same format as RSeQC's bam_stat.py.
///
/// # Arguments
/// * `result` - The computed statistics
/// * `mapq_cut` - The MAPQ cutoff that was used
/// * `output_path` - Path to write the output file
pub fn write_bam_stat(result: &BamStatResult, _mapq_cut: u8, output_path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(output_path)
        .with_context(|| format!("Failed to create output file: {}", output_path.display()))?;

    // Header — matches RSeQC's exact format (50-char separator)
    writeln!(out, "\n#==================================================")?;
    writeln!(out, "#All numbers are READ count")?;
    writeln!(out, "#==================================================\n")?;

    // Main statistics — format matches RSeQC's %-40s%d
    // Note: RSeQC inserts blank lines after "Total records" and before "mapq >= ..."
    // Note: "Non primary hits" has no trailing colon in RSeQC
    writeln!(out, "{:<40}{}\n", "Total records:", result.total_records)?;
    writeln!(out, "{:<40}{}", "QC failed:", result.qc_failed)?;
    writeln!(out, "{:<40}{}", "Optical/PCR duplicate:", result.duplicates)?;
    writeln!(out, "{:<40}{}", "Non primary hits", result.non_primary)?;
    writeln!(out, "{:<40}{}", "Unmapped reads:", result.unmapped)?;
    writeln!(
        out,
        "{:<40}{}\n",
        "mapq < mapq_cut (non-unique):", result.non_unique
    )?;
    writeln!(out, "{:<40}{}", "mapq >= mapq_cut (unique):", result.unique)?;
    writeln!(out, "{:<40}{}", "Read-1:", result.read_1)?;
    writeln!(out, "{:<40}{}", "Read-2:", result.read_2)?;
    writeln!(out, "{:<40}{}", "Reads map to '+':", result.forward)?;
    writeln!(out, "{:<40}{}", "Reads map to '-':", result.reverse)?;
    writeln!(out, "{:<40}{}", "Non-splice reads:", result.non_splice)?;
    writeln!(out, "{:<40}{}", "Splice reads:", result.splice)?;
    writeln!(
        out,
        "{:<40}{}",
        "Reads mapped in proper pairs:", result.proper_pairs
    )?;
    writeln!(
        out,
        "Proper-paired reads map to different chrom:{}",
        result.proper_pair_diff_chrom
    )?;

    info!("Wrote bam_stat output to {}", output_path.display());
    Ok(())
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bam_stat_small() {
        // Uses the tiny test BAM from test data
        let result = bam_stat("tests/data/test.bam", 30, None);
        assert!(result.is_ok(), "bam_stat should succeed on test.bam");
        let r = result.unwrap();
        // Basic sanity checks
        assert!(r.total_records > 0, "Should have some records");
        assert_eq!(
            r.total_records,
            r.qc_failed + r.duplicates + r.non_primary + r.unmapped + r.non_unique + r.unique,
            "Counts should sum to total"
        );
    }
}
