//! Basic BAM alignment statistics.
//!
//! Reimplements RSeQC's `bam_stat.py`. Statistics are collected via the
//! `BamStatAccum` accumulator (in `accumulators.rs`) during the shared BAM pass,
//! then converted to `BamStatResult` for output.

use anyhow::{Context, Result};
use log::info;
use std::collections::{BTreeMap, HashMap};
use std::io::Write;
use std::path::Path;

// ============================================================================
// Result struct
// ============================================================================

/// Holds all statistics collected during the BAM pass.
///
/// Contains both the original RSeQC bam_stat fields and the additional
/// fields needed for samtools-compatible flagstat, idxstats, and stats output.
#[derive(Debug, Default)]
pub struct BamStatResult {
    // --- RSeQC bam_stat fields ---
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
    #[allow(dead_code)] // Kept for API completeness; not yet used by write_bam_stat
    pub mapq_distribution: BTreeMap<u8, u64>,

    // --- samtools flagstat fields ---
    /// Secondary alignments (0x100).
    pub secondary: u64,
    /// Supplementary alignments (0x800).
    pub supplementary: u64,
    /// All mapped records (not 0x4).
    pub mapped: u64,
    /// Paired reads (0x1).
    pub paired_flagstat: u64,
    /// Read1 in pair (0x40).
    pub read1_flagstat: u64,
    /// Read2 in pair (0x80).
    pub read2_flagstat: u64,
    /// First fragments for samtools stats (read1 or unpaired).
    pub first_fragments: u64,
    /// Last fragments for samtools stats (read2).
    pub last_fragments: u64,
    /// Properly paired reads (0x1 + 0x2).
    pub properly_paired: u64,
    /// Both mates mapped.
    pub both_mapped: u64,
    /// Singletons (this mapped, mate unmapped).
    pub singletons: u64,
    /// Paired, both mapped, different reference.
    pub mate_diff_chr: u64,
    /// Paired, both mapped, different reference, MAPQ >= 5.
    pub mate_diff_chr_mapq5: u64,

    // --- samtools idxstats fields ---
    /// Per-reference (tid) mapped and unmapped counts.
    pub chrom_counts: HashMap<i32, (u64, u64)>,
    /// Unmapped reads with no reference (tid < 0).
    pub unplaced_unmapped: u64,

    // --- samtools stats SN fields ---
    /// Sum of query sequence lengths for all primary reads.
    pub total_len: u64,
    /// Sum of first fragment (read1 or unpaired) sequence lengths.
    pub total_first_fragment_len: u64,
    /// Sum of last fragment (read2) sequence lengths.
    pub total_last_fragment_len: u64,
    /// Sum of query lengths for mapped primary reads.
    pub bases_mapped: u64,
    /// Sum of M/=/X CIGAR operations for mapped primary reads.
    pub bases_mapped_cigar: u64,
    /// Sum of query lengths for duplicate-flagged primary reads.
    pub bases_duplicated: u64,
    /// Maximum query sequence length (among primary reads).
    pub max_len: u64,
    /// Maximum first-fragment sequence length.
    pub max_first_fragment_len: u64,
    /// Maximum last-fragment sequence length.
    pub max_last_fragment_len: u64,
    /// Sum of average per-read base qualities.
    pub quality_sum: f64,
    /// Number of reads contributing to quality_sum.
    pub quality_count: u64,
    /// Sum of NM tag values across mapped primary reads.
    pub mismatches: u64,
    /// Sum of absolute TLEN for properly paired primary reads.
    pub insert_size_sum: f64,
    /// Sum of squared TLEN for insert size standard deviation.
    pub insert_size_sq_sum: f64,
    /// Number of reads contributing to insert size stats.
    pub insert_size_count: u64,
    /// Inward-oriented pairs (FR).
    pub inward_pairs: u64,
    /// Outward-oriented pairs (RF).
    pub outward_pairs: u64,
    /// Other orientation pairs (FF, RR).
    pub other_orientation: u64,
    /// Primary reads that are paired.
    #[allow(dead_code)]
    pub primary_paired: u64,
    /// Total primary reads (non-secondary, non-supplementary).
    pub primary_count: u64,
    /// Primary mapped reads count.
    pub primary_mapped: u64,
    /// Primary duplicate reads.
    pub primary_duplicates: u64,
    /// Number of mapped reads with mapping quality 0 (all reads, not just primary).
    pub reads_mq0: u64,
    /// Number of primary, non-QC-fail, mapped, paired reads whose mate is also mapped.
    pub reads_mapped_and_paired: u64,
}

// ============================================================================
// Output formatting
// ============================================================================

/// Write bam_stat results in the same format as RSeQC's bam_stat.py.
///
/// # Arguments
/// * `result` - The computed statistics
/// * `output_path` - Path to write the output file
pub fn write_bam_stat(result: &BamStatResult, output_path: &Path) -> Result<()> {
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
    use crate::rna::rseqc::accumulators::BamStatAccum;
    use rust_htslib::bam::{self, Read as BamRead};

    #[test]
    fn test_bam_stat_small() {
        // Uses the tiny test BAM from test data via the accumulator path
        let mut reader =
            bam::Reader::from_path("tests/data/test.bam").expect("Failed to open test.bam");
        let mut accum = BamStatAccum::default();
        let mut record = bam::Record::new();

        while let Some(res) = reader.read(&mut record) {
            res.expect("Error reading BAM record");
            accum.process_read(&record, 30);
        }

        let r = accum.into_result();
        // Basic sanity checks
        assert!(r.total_records > 0, "Should have some records");
        assert_eq!(
            r.total_records,
            r.qc_failed + r.duplicates + r.non_primary + r.unmapped + r.non_unique + r.unique,
            "Counts should sum to total"
        );
    }
}
