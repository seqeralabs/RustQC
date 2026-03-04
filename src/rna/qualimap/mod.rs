//! Qualimap RNA-Seq QC module.
//!
//! Reimplements Qualimap's RNA-Seq quality control analysis with enclosure-based
//! gene assignment, per-transcript coverage tracking, and Qualimap-compatible output.

pub mod accumulator;
pub mod coverage;
pub mod index;
pub mod output;
pub mod plots;
pub mod report;

pub use accumulator::QualimapAccum;
pub use index::QualimapIndex;

use std::collections::HashMap;

// ============================================================================
// Result types
// ============================================================================

/// Final Qualimap result, produced from a merged `QualimapAccum`.
///
/// Holds all counters and coverage data needed to produce Qualimap-compatible
/// output files (coverage profiles, bias metrics, `rnaseq_qc_results.txt`).
#[derive(Debug)]
pub struct QualimapResult {
    // --- Read counters ---
    /// Total primary alignments (excluding unmapped, secondary, supplementary, QC-fail).
    pub primary_alignments: u64,
    /// Total secondary alignments seen.
    pub secondary_alignments: u64,
    /// Reads skipped because they were not aligned.
    pub not_aligned: u64,
    /// Reads skipped because NH > 1 (multi-mappers).
    pub alignment_not_unique: u64,

    // --- Gene assignment counters ---
    /// Reads assigned to exactly one gene (exonic).
    pub exonic_reads: u64,
    /// Reads where >1 gene encloses all blocks (ambiguous).
    pub ambiguous_reads: u64,
    /// Reads with no gene enclosing all blocks.
    pub no_feature: u64,
    /// Reads classified as intronic (M-block overlaps intron tree).
    pub intronic_reads: u64,
    /// Reads classified as intergenic (no exon or intron overlap).
    pub intergenic_reads: u64,
    /// Reads overlapping exons from multiple genes.
    pub overlapping_exon_reads: u64,

    // --- Fragment counters (PE) ---
    /// Total reads counted (PE: each mate counts as 1 read).
    pub read_count: u64,
    /// Total fragments counted (PE: 1 per paired fragment).
    #[allow(dead_code)] // populated but not yet consumed in output
    pub fragment_count: u64,
    /// Left-of-pair reads (first-of-pair in paired mode).
    pub left_proper_in_pair: u64,
    /// Right-of-pair reads (second-of-pair in paired mode).
    pub right_proper_in_pair: u64,
    /// Both mates with proper-pair flag (for numberOfMappedPairs = count / 2).
    pub both_proper_in_pair: u64,

    // --- Junction counters ---
    /// Reads containing at least one splice junction (N-op in CIGAR).
    pub reads_at_junctions: u64,
    /// Junction motif counts: canonical_motif_string -> count.
    pub junction_motifs: HashMap<String, u64>,

    // --- Per-transcript coverage ---
    /// Per-transcript coverage arrays (transcript_key -> per-base coverage).
    #[allow(dead_code)] // populated but not yet consumed in output
    pub transcript_coverage: HashMap<String, Vec<i32>>,

    /// Raw per-transcript coverage keyed by flat index (for bias computation).
    /// Retains the internal indexing used during accumulation so that bias
    /// computation can look up transcript metadata (strand, exonic length)
    /// efficiently via the `QualimapIndex`.
    pub transcript_coverage_raw: coverage::TranscriptCoverage,

    // --- Merged gene coverage ---
    /// Per-gene coverage using merged exon models (gene_idx -> per-base coverage).
    /// Qualimap's Picard Gene.Transcript merges all exons from all isoforms into
    /// a single non-redundant exon set per gene. This produces more accurate
    /// bias and profile values than per-transcript tracking.
    #[allow(dead_code)] // populated but not yet consumed in output
    pub merged_gene_coverage: coverage::MergedGeneCoverage,

    /// Forward strand estimation count for SSP.
    pub ssp_fwd: u64,
    /// Reverse strand estimation count for SSP.
    pub ssp_rev: u64,
}
