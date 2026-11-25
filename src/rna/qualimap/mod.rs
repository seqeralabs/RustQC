//! Qualimap RNA-Seq QC module.
//!
//! Reimplements Qualimap's RNA-Seq quality control analysis with enclosure-based
//! gene assignment, per-transcript coverage tracking, and Qualimap-compatible output.

pub mod accumulator;
pub mod coverage;
mod format;
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
    /// Mapped records not flagged secondary (supplementary and QC-fail included).
    pub primary_alignments: u64,
    /// Total secondary alignments seen.
    pub secondary_alignments: u64,
    /// Reads skipped because they were not aligned.
    pub not_aligned: u64,
    /// Reads skipped because NH > 1 (multi-mappers).
    pub alignment_not_unique: u64,
    /// Mapped reads on a sequence with no annotated exon ("Missing chromosome
    /// in annotation" in the HTML report).
    pub seq_not_found: u64,

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
    /// Mapped records on an annotated sequence (Qualimap `readCount`).
    pub read_count: u64,
    /// Left-of-pair reads (first-of-pair in paired mode).
    pub left_proper_in_pair: u64,
    /// Right-of-pair reads (second-of-pair in paired mode).
    pub right_proper_in_pair: u64,
    /// Both mates with proper-pair flag (for numberOfMappedPairs = count / 2).
    pub both_proper_in_pair: u64,

    // --- Junction counters ---
    /// Reads containing at least one splice junction (N-op in CIGAR).
    pub reads_at_junctions: u64,
    /// Junction events classified as "known" (both endpoints match known exon boundaries).
    pub known_junction_events: u64,
    /// Junction events classified as "partly known" (one endpoint matches).
    pub partly_known_junction_events: u64,
    /// Junction motif counts: canonical_motif_string -> count.
    pub junction_motifs: HashMap<String, u64>,

    // --- Per-transcript coverage ---
    /// Raw per-transcript coverage keyed by flat index (for bias computation).
    /// Retains the internal indexing used during accumulation so that bias
    /// computation can look up transcript metadata (strand, exonic length)
    /// efficiently via the `QualimapIndex`.
    pub transcript_coverage_raw: coverage::TranscriptCoverage,

    /// Forward strand estimation count for SSP.
    pub ssp_fwd: u64,
    /// Reverse strand estimation count for SSP.
    pub ssp_rev: u64,
}

// ============================================================================
// Junction motif ordering
// ============================================================================

/// Number of junction motifs listed in the report. Qualimap's output loop runs
/// while `count <= 10`, so it prints at most 11 motifs.
pub(crate) const TOP_JUNCTION_MOTIFS: usize = 11;

/// Order junction motifs for output: most frequent first, ties broken by motif
/// name (alphabetical).
///
/// Qualimap 2.3 keeps the motif counts in a `java.util.HashMap`, sorts them by
/// percentage with a stable sort and prints the list backwards
/// (`TranscriptDataHandler.computeSortedJunctionsMap` and the junction loops in
/// `RNASeqQCAnalysis`). Tied motifs therefore come out in reverse HashMap
/// iteration order. That order depends on the table size and, for keys that
/// share a bucket, on the order the motifs were first seen in the
/// (name-sorted, for paired-end) input. RustQC reads the BAM in a different
/// order, so it cannot reproduce Qualimap's tie order. It uses the motif name
/// instead, which keeps the output identical across runs and thread counts.
pub(crate) fn sorted_junction_motifs(motifs: &HashMap<String, u64>) -> Vec<(&str, u64)> {
    let mut sorted: Vec<(&str, u64)> = motifs.iter().map(|(m, &c)| (m.as_str(), c)).collect();
    sorted.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(b.0)));
    sorted
}

/// The junction motif rows Qualimap prints: up to [`TOP_JUNCTION_MOTIFS`]
/// motifs with their share of `reads_at_junctions`, computed as
/// `(count * 100.) / numReadsWithJunction` like
/// `TranscriptDataHandler.computeSortedJunctionsMap`. Empty when no read has a
/// junction.
pub(crate) fn junction_motif_rows(
    motifs: &HashMap<String, u64>,
    reads_at_junctions: u64,
) -> Vec<(&str, f64)> {
    if reads_at_junctions == 0 {
        return Vec::new();
    }
    sorted_junction_motifs(motifs)
        .into_iter()
        .take(TOP_JUNCTION_MOTIFS)
        .map(|(motif, count)| (motif, (count as f64 * 100.0) / reads_at_junctions as f64))
        .collect()
}

/// Qualimap's strand-specificity estimation (fwd, rev).
///
/// `ComputeCountsTask.getFwdStrandEstimation` divides in single precision,
/// `fwd / (float)(fwd + rev)`, so the value is computed as an `f32` and only
/// then widened for formatting. Both are NaN when no read overlapped a feature.
pub(crate) fn strand_estimation(fwd: u64, rev: u64) -> (f64, f64) {
    let total = (fwd + rev) as f32;
    (f64::from(fwd as f32 / total), f64::from(rev as f32 / total))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn motif_map(entries: &[(&str, u64)]) -> HashMap<String, u64> {
        entries.iter().map(|&(m, c)| (m.to_string(), c)).collect()
    }

    #[test]
    fn junction_motifs_sorted_by_count_descending() {
        let motifs = motif_map(&[("GTAG", 5), ("CTAC", 90), ("GCAG", 20)]);
        assert_eq!(
            sorted_junction_motifs(&motifs),
            vec![("CTAC", 90), ("GCAG", 20), ("GTAG", 5)]
        );
    }

    #[test]
    fn junction_motif_ties_broken_by_name() {
        // Four motifs tied at the same count, as seen in CI where their order
        // used to change from run to run.
        let motifs = motif_map(&[
            ("TGGT", 1),
            ("AGGT", 1),
            ("ACAA", 1),
            ("TTTT", 1),
            ("GTAG", 6),
        ]);
        assert_eq!(
            sorted_junction_motifs(&motifs),
            vec![
                ("GTAG", 6),
                ("ACAA", 1),
                ("AGGT", 1),
                ("TGGT", 1),
                ("TTTT", 1)
            ]
        );
    }

    #[test]
    fn junction_motif_order_independent_of_insertion_order() {
        let entries: Vec<(String, u64)> = (0..64u64)
            .map(|i| {
                let bases = b"ACGT";
                let motif: String = (0..4)
                    .map(|p| bases[((i >> (2 * p)) & 3) as usize] as char)
                    .collect();
                (motif, i % 3)
            })
            .collect();
        let forward: HashMap<String, u64> = entries.iter().cloned().collect();
        let mut reversed = HashMap::with_capacity(1);
        for (m, c) in entries.iter().rev() {
            reversed.insert(m.clone(), *c);
        }
        // Many ties (64 motifs, 3 distinct counts), so the top-11 cut falls
        // inside a tie group and must not depend on insertion order.
        let expected = sorted_junction_motifs(&forward);
        assert_eq!(expected, sorted_junction_motifs(&reversed));
        assert!(expected
            .windows(2)
            .all(|w| w[0].1 > w[1].1 || (w[0].1 == w[1].1 && w[0].0 < w[1].0)));
    }
}
