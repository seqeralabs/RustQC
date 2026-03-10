//! Basic BAM alignment statistics.
//!
//! Reimplements RSeQC's `bam_stat.py`. Statistics are collected via the
//! `BamStatAccum` accumulator (in `accumulators.rs`) during the shared BAM pass,
//! then converted to `BamStatResult` for output.

use anyhow::{Context, Result};
use log::debug;
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

// ============================================================================
// Result struct
// ============================================================================

/// Holds all statistics collected during the BAM pass.
///
/// Contains both the original RSeQC bam_stat fields and the additional
/// fields needed for samtools-compatible flagstat, idxstats, and stats output.
#[derive(Debug)]
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
    /// Sum of M/I/=/X CIGAR operations for mapped primary reads.
    pub bases_mapped_cigar: u64,
    /// Sum of query lengths for duplicate-flagged primary reads.
    pub bases_duplicated: u64,
    /// Maximum query sequence length (among primary reads).
    pub max_len: u64,
    /// Maximum first-fragment sequence length.
    pub max_first_fragment_len: u64,
    /// Maximum last-fragment sequence length.
    pub max_last_fragment_len: u64,
    /// Sum of individual base qualities across all primary non-QC-fail reads.
    pub quality_sum: f64,
    /// Total number of bases contributing to quality_sum.
    pub quality_count: u64,
    /// Sum of NM tag values across mapped primary reads.
    pub mismatches: u64,
    /// Insert size with orientation breakdown: abs(TLEN) → [total, inward, outward, other].
    /// Both mates contribute; divide by 2 at output to match samtools stats.
    pub is_hist: HashMap<u64, [u64; 4]>,
    /// Inward-oriented pairs (FR).
    pub inward_pairs: u64,
    /// Outward-oriented pairs (RF).
    pub outward_pairs: u64,
    /// Other orientation pairs (FF, RR).
    pub other_orientation: u64,
    /// Total primary reads (non-secondary, non-supplementary).
    pub primary_count: u64,
    /// Primary mapped reads count.
    pub primary_mapped: u64,
    /// Primary duplicate reads.
    pub primary_duplicates: u64,
    /// Number of primary mapped reads with mapping quality 0.
    pub reads_mq0: u64,
    /// Number of primary, non-QC-fail, mapped, paired reads whose mate is also mapped.
    pub reads_mapped_and_paired: u64,

    // --- samtools stats histogram/distribution fields ---
    /// Read length histogram (all primary reads): length → count.
    pub rl_hist: HashMap<u64, u64>,
    /// First fragment read length histogram: length → count.
    pub frl_hist: HashMap<u64, u64>,
    /// Last fragment read length histogram: length → count.
    pub lrl_hist: HashMap<u64, u64>,
    /// MAPQ histogram (primary, mapped, !qcfail, !dup): quality 0-255.
    pub mapq_hist: [u64; 256],
    /// Per-cycle quality for first fragments: outer = cycle, inner = quality bucket counts.
    pub ffq: Vec<[u64; 64]>,
    /// Per-cycle quality for last fragments: outer = cycle, inner = quality bucket counts.
    pub lfq: Vec<[u64; 64]>,
    /// GC content step-function for first fragments, 200 bins (samtools ngc=200).
    pub gcf: [u64; 200],
    /// GC content step-function for last fragments, 200 bins (samtools ngc=200).
    pub gcl: [u64; 200],
    /// Per-cycle base composition for first fragments: [A, C, G, T, N, Other] per cycle.
    #[allow(dead_code)]
    pub fbc: Vec<[u64; 6]>,
    /// Per-cycle base composition for last fragments: [A, C, G, T, N, Other] per cycle.
    #[allow(dead_code)]
    pub lbc: Vec<[u64; 6]>,
    /// Per-cycle base composition (read-oriented) for first fragments.
    pub fbc_ro: Vec<[u64; 6]>,
    /// Per-cycle base composition (read-oriented) for last fragments.
    pub lbc_ro: Vec<[u64; 6]>,
    /// Per-cycle base composition (reverse-complemented for reverse-strand reads,
    /// combined first+last fragments). Used for GCT output. [A, C, G, T] only.
    pub gcc_rc: Vec<[u64; 4]>,
    /// Total base counters for first fragments: [A, C, G, T, N].
    pub ftc: [u64; 5],
    /// Total base counters for last fragments: [A, C, G, T, N].
    pub ltc: [u64; 5],
    /// Indel distribution by size: length → [insertions, deletions].
    pub id_hist: HashMap<u64, [u64; 2]>,
    /// Indels per cycle: cycle → [ins_fwd, ins_rev, del_fwd, del_rev].
    pub ic: Vec<[u64; 4]>,
    /// CRC32 checksum sums: [names, sequences, qualities].
    pub chk: [u32; 3],
    /// Coverage distribution: depth → number of reference positions at that depth.
    pub cov_hist: HashMap<u32, u64>,
    /// GC-depth bins for the GCD section of samtools stats.
    pub gcd_bins: Vec<GcDepthBin>,
}

/// A single GC-depth bin, representing a `gcd_bin_size`-bp window of the genome.
#[derive(Debug, Clone)]
pub struct GcDepthBin {
    /// Accumulated GC fraction sum (sum of per-read `gc_count / seq_len`).
    /// Divided by `depth` during output to get the mean GC fraction.
    pub gc: f32,
    /// Number of reads whose start position falls in this bin.
    pub depth: u32,
}

impl Default for BamStatResult {
    fn default() -> Self {
        Self {
            total_records: 0,
            qc_failed: 0,
            duplicates: 0,
            non_primary: 0,
            unmapped: 0,
            non_unique: 0,
            unique: 0,
            read_1: 0,
            read_2: 0,
            forward: 0,
            reverse: 0,
            splice: 0,
            non_splice: 0,
            proper_pairs: 0,
            proper_pair_diff_chrom: 0,
            secondary: 0,
            supplementary: 0,
            mapped: 0,
            paired_flagstat: 0,
            read1_flagstat: 0,
            read2_flagstat: 0,
            first_fragments: 0,
            last_fragments: 0,
            properly_paired: 0,
            both_mapped: 0,
            singletons: 0,
            mate_diff_chr: 0,
            mate_diff_chr_mapq5: 0,
            chrom_counts: HashMap::new(),
            unplaced_unmapped: 0,
            total_len: 0,
            total_first_fragment_len: 0,
            total_last_fragment_len: 0,
            bases_mapped: 0,
            bases_mapped_cigar: 0,
            bases_duplicated: 0,
            max_len: 0,
            max_first_fragment_len: 0,
            max_last_fragment_len: 0,
            quality_sum: 0.0,
            quality_count: 0,
            mismatches: 0,
            is_hist: HashMap::new(),
            inward_pairs: 0,
            outward_pairs: 0,
            other_orientation: 0,
            primary_count: 0,
            primary_mapped: 0,
            primary_duplicates: 0,
            reads_mq0: 0,
            reads_mapped_and_paired: 0,
            rl_hist: HashMap::new(),
            frl_hist: HashMap::new(),
            lrl_hist: HashMap::new(),
            mapq_hist: [0u64; 256],
            ffq: Vec::new(),
            lfq: Vec::new(),
            gcf: [0u64; 200],
            gcl: [0u64; 200],
            fbc: Vec::new(),
            lbc: Vec::new(),
            fbc_ro: Vec::new(),
            lbc_ro: Vec::new(),
            gcc_rc: Vec::new(),
            ftc: [0u64; 5],
            ltc: [0u64; 5],
            id_hist: HashMap::new(),
            ic: Vec::new(),
            chk: [0u32; 3],
            cov_hist: HashMap::new(),
            gcd_bins: Vec::new(),
        }
    }
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

    debug!("Wrote bam_stat output to {}", output_path.display());
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
