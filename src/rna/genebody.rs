//! Gene body coverage profiling and Qualimap-compatible output.
//!
//! Computes coverage profiles along transcript bodies, mapping read positions
//! to 100 percentile bins (5'→3'). Produces Qualimap rnaseq-compatible output
//! for MultiQC parsing.

use anyhow::{Context, Result};
use indexmap::IndexMap;
use log::{debug, info};
use std::io::Write;
use std::path::Path;

use crate::gtf::Gene;

/// Number of percentile bins for coverage profiling (0–99).
const NUM_BINS: usize = 100;

// ============================================================
// Gene exon offset map — precomputed for fast position mapping
// ============================================================

/// Per-gene exon structure for mapping genomic positions to transcript percentiles.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct GeneExonMap {
    /// Sorted exon intervals in 0-based half-open coordinates `(start, end)`.
    pub exons: Vec<(u64, u64)>,
    /// Cumulative start offsets for each exon (transcript-relative).
    /// `cumulative_starts[i]` = sum of lengths of exons 0..i.
    pub cumulative_starts: Vec<u64>,
    /// Total exonic length (sum of all exon lengths).
    pub total_exonic_length: u64,
    /// Strand ('+' or '-').
    pub strand: char,
}

/// Maps gene indices to their exon position data for fast percentile bin lookup.
#[derive(Debug)]
pub struct TranscriptPositionMap {
    /// Per-gene exon maps, indexed by gene position in the IndexMap.
    genes: Vec<Option<GeneExonMap>>,
}

impl TranscriptPositionMap {
    /// Build the position map from GTF gene annotations.
    ///
    /// Uses the longest transcript per gene for exon structure.
    /// Genes with no exons or zero exonic length are set to `None`.
    pub fn from_genes(genes: &IndexMap<String, Gene>) -> Self {
        let mut maps = Vec::with_capacity(genes.len());

        for gene in genes.values() {
            // Use the longest transcript (by total exonic length)
            let best_tx = gene.transcripts.iter().max_by_key(|tx| {
                tx.exons
                    .iter()
                    .map(|(s, e)| e.saturating_sub(*s))
                    .sum::<u64>()
            });

            let exon_map = best_tx.and_then(|tx| {
                // Convert GTF 1-based inclusive to 0-based half-open
                let mut exons: Vec<(u64, u64)> = tx
                    .exons
                    .iter()
                    .map(|(s, e)| (s.saturating_sub(1), *e))
                    .collect();

                // Sort by start position
                exons.sort_unstable_by_key(|(s, _)| *s);

                let total_exonic_length: u64 =
                    exons.iter().map(|(s, e)| e.saturating_sub(*s)).sum();

                if total_exonic_length == 0 || exons.is_empty() {
                    return None;
                }

                let mut cumulative_starts = Vec::with_capacity(exons.len());
                let mut cum = 0u64;
                for (s, e) in &exons {
                    cumulative_starts.push(cum);
                    cum += e - s;
                }

                Some(GeneExonMap {
                    exons,
                    cumulative_starts,
                    total_exonic_length,
                    strand: gene.strand,
                })
            });

            maps.push(exon_map);
        }

        Self { genes: maps }
    }

    /// Get the exon map for a gene by index.
    #[allow(dead_code)]
    pub fn get(&self, gene_idx: usize) -> Option<&GeneExonMap> {
        self.genes.get(gene_idx).and_then(|m| m.as_ref())
    }

    /// Number of genes in the map.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.genes.len()
    }
}

// ============================================================
// Coverage accumulator — accumulated during BAM pass
// ============================================================

/// Accumulator for gene body coverage profiling.
///
/// Tracks coverage bins and genomic origin counts during BAM processing.
#[derive(Debug, Clone)]
pub struct GenebodyCoverageAccum {
    /// Coverage profile bins (100 percentile bins, 5'→3').
    pub coverage_bins: [f64; NUM_BINS],
    /// Number of reads/fragments assigned to genes.
    pub aligned_to_genes: u64,
    /// Number of ambiguous reads (overlapping multiple genes).
    pub ambiguous: u64,
    /// Number of reads with no gene feature.
    pub no_feature: u64,
    /// Exonic base count (reads overlapping exons).
    pub exonic: u64,
    /// Intronic base count (reads in introns of assigned genes).
    pub intronic: u64,
    /// Intergenic base count (reads not overlapping any gene).
    pub intergenic: u64,
    /// Reads overlapping multiple exons of the same gene.
    pub overlapping_exon: u64,
}

impl Default for GenebodyCoverageAccum {
    fn default() -> Self {
        Self {
            coverage_bins: [0.0; NUM_BINS],
            aligned_to_genes: 0,
            ambiguous: 0,
            no_feature: 0,
            exonic: 0,
            intronic: 0,
            intergenic: 0,
            overlapping_exon: 0,
        }
    }
}

impl GenebodyCoverageAccum {
    /// Create a new empty accumulator.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record a read's gene body coverage by mapping aligned blocks to percentile bins.
    ///
    /// # Arguments
    /// * `gene_idx` - Index of the assigned gene
    /// * `aligned_blocks` - Read's aligned genomic blocks (0-based half-open)
    /// * `position_map` - Precomputed exon position map
    pub fn record_coverage(
        &mut self,
        gene_idx: usize,
        aligned_blocks: &[(u64, u64)],
        position_map: &TranscriptPositionMap,
    ) {
        let gene_map = match position_map.get(gene_idx) {
            Some(m) => m,
            None => return,
        };

        self.aligned_to_genes += 1;

        // Track exon overlap count for overlapping_exon metric
        let mut exons_hit = 0u32;

        // Map each aligned block to transcript-relative positions
        for &(block_start, block_end) in aligned_blocks {
            let mut exonic_bases: u64 = 0;

            for (i, &(exon_start, exon_end)) in gene_map.exons.iter().enumerate() {
                // Find overlap between aligned block and this exon
                let overlap_start = block_start.max(exon_start);
                let overlap_end = block_end.min(exon_end);

                if overlap_start >= overlap_end {
                    continue;
                }

                exonic_bases += overlap_end - overlap_start;
                exons_hit += 1;

                // Map overlapping bases to transcript-relative positions
                let tx_start = gene_map.cumulative_starts[i] + (overlap_start - exon_start);
                let tx_end = gene_map.cumulative_starts[i] + (overlap_end - exon_start);

                // Convert to percentile bins using range-based approach
                // instead of iterating per-base (O(bins_spanned) vs O(overlap_len))
                let total_len = gene_map.total_exonic_length as f64;
                let (range_start, range_end) = if gene_map.strand == '-' {
                    // Reverse strand: flip to 5'→3'
                    let rs = gene_map.total_exonic_length - tx_end;
                    let re = gene_map.total_exonic_length - tx_start;
                    (rs, re)
                } else {
                    (tx_start, tx_end)
                };

                // Compute start/end bins (clamped to valid range)
                let start_bin = ((range_start as f64 / total_len) * NUM_BINS as f64) as usize;
                let start_bin = start_bin.min(NUM_BINS - 1);
                let end_bin = (((range_end - 1) as f64 / total_len) * NUM_BINS as f64) as usize;
                let end_bin = end_bin.min(NUM_BINS - 1);

                if start_bin == end_bin {
                    // All bases fall in the same bin
                    self.coverage_bins[start_bin] += (range_end - range_start) as f64;
                } else {
                    // Bases span multiple bins -- compute exact counts per bin
                    let bin_width = total_len / NUM_BINS as f64;
                    // First bin: bases from range_start to the end of start_bin
                    let first_bin_end = ((start_bin + 1) as f64 * bin_width).ceil() as u64;
                    self.coverage_bins[start_bin] +=
                        (first_bin_end.min(range_end) - range_start) as f64;
                    // Middle bins: each gets a full bin_width of bases
                    for bin in (start_bin + 1)..end_bin {
                        let bin_start = (bin as f64 * bin_width).ceil() as u64;
                        let bin_end = ((bin + 1) as f64 * bin_width).ceil() as u64;
                        self.coverage_bins[bin] +=
                            (bin_end.min(range_end) - bin_start.max(range_start)) as f64;
                    }
                    // Last bin: bases from the start of end_bin to range_end
                    let last_bin_start = (end_bin as f64 * bin_width).ceil() as u64;
                    self.coverage_bins[end_bin] +=
                        (range_end - last_bin_start.max(range_start)) as f64;
                }
            }

            // Count exonic bases from actual overlap, remainder is intronic.
            // Cap exonic_bases at block_len in case overlapping exons
            // (from merged transcript models) double-count the same bases.
            let block_len = block_end - block_start;
            let exonic_capped = exonic_bases.min(block_len);
            self.exonic += exonic_capped;
            self.intronic += block_len - exonic_capped;
        }

        if exons_hit > 1 {
            self.overlapping_exon += 1;
        }
    }

    /// Record a read that was not assigned to any gene (intergenic).
    pub fn record_no_feature(&mut self, aligned_blocks: &[(u64, u64)]) {
        self.no_feature += 1;
        for &(start, end) in aligned_blocks {
            self.intergenic += end - start;
        }
    }

    /// Record an ambiguous read (overlapping multiple genes).
    pub fn record_ambiguous(&mut self) {
        self.ambiguous += 1;
    }

    /// Merge another accumulator into this one.
    pub fn merge(&mut self, other: &GenebodyCoverageAccum) {
        for i in 0..NUM_BINS {
            self.coverage_bins[i] += other.coverage_bins[i];
        }
        self.aligned_to_genes += other.aligned_to_genes;
        self.ambiguous += other.ambiguous;
        self.no_feature += other.no_feature;
        self.exonic += other.exonic;
        self.intronic += other.intronic;
        self.intergenic += other.intergenic;
        self.overlapping_exon += other.overlapping_exon;
    }
}

// ============================================================
// Results and output
// ============================================================

/// Gene body coverage results.
#[derive(Debug)]
pub struct GenebodyCoverageResult {
    /// Coverage profile bins (100 percentile bins).
    pub coverage_bins: [f64; NUM_BINS],
    /// Reads aligned to genes.
    pub aligned_to_genes: u64,
    /// Ambiguous reads.
    pub ambiguous: u64,
    /// Reads with no feature.
    pub no_feature: u64,
    /// Exonic bases.
    pub exonic: u64,
    /// Intronic bases.
    pub intronic: u64,
    /// Intergenic bases.
    pub intergenic: u64,
    /// Reads overlapping multiple exons.
    pub overlapping_exon: u64,
    /// 5' bias (mean of first 20 bins / overall mean).
    pub bias_5prime: f64,
    /// 3' bias (mean of last 20 bins / overall mean).
    pub bias_3prime: f64,
    /// 5'-3' bias (median of first 20 bins / median of last 20 bins).
    pub bias_5_to_3: f64,
}

impl GenebodyCoverageAccum {
    /// Convert accumulator to final results with bias calculations.
    pub fn into_result(self) -> GenebodyCoverageResult {
        let overall_mean: f64 = self.coverage_bins.iter().sum::<f64>() / NUM_BINS as f64;

        // 5' bias: mean of bins 0-19 / overall mean
        let first20_mean: f64 = self.coverage_bins[..20].iter().sum::<f64>() / 20.0;
        let bias_5prime = if overall_mean > 0.0 {
            first20_mean / overall_mean
        } else {
            0.0
        };

        // 3' bias: mean of bins 80-99 / overall mean
        let last20_mean: f64 = self.coverage_bins[80..].iter().sum::<f64>() / 20.0;
        let bias_3prime = if overall_mean > 0.0 {
            last20_mean / overall_mean
        } else {
            0.0
        };

        // 5'-3' bias: median of first 20 bins / median of last 20 bins
        let median_first20 = crate::io::median(&self.coverage_bins[..20]);
        let median_last20 = crate::io::median(&self.coverage_bins[80..]);
        let bias_5_to_3 = if median_last20 > 0.0 {
            median_first20 / median_last20
        } else {
            0.0
        };

        GenebodyCoverageResult {
            coverage_bins: self.coverage_bins,
            aligned_to_genes: self.aligned_to_genes,
            ambiguous: self.ambiguous,
            no_feature: self.no_feature,
            exonic: self.exonic,
            intronic: self.intronic,
            intergenic: self.intergenic,
            overlapping_exon: self.overlapping_exon,
            bias_5prime,
            bias_3prime,
            bias_5_to_3,
        }
    }
}

/// Write the coverage profile along genes in Qualimap format.
///
/// Format: TSV with columns `position` (0.0–99.0) and `coverage` (cumulative depth).
pub fn write_coverage_profile(result: &GenebodyCoverageResult, output_path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(output_path)
        .with_context(|| format!("Failed to create {}", output_path.display()))?;

    for (i, &val) in result.coverage_bins.iter().enumerate() {
        writeln!(f, "{:.1}\t{}", i as f64, val)?;
    }

    debug!("Wrote coverage profile to {}", output_path.display());
    Ok(())
}

/// Format a number with comma-separated thousands.
fn format_with_commas(n: u64) -> String {
    let s = n.to_string();
    let mut result = String::with_capacity(s.len() + s.len() / 3);
    for (i, c) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            result.push(',');
        }
        result.push(c);
    }
    result.chars().rev().collect()
}

/// Write Qualimap rnaseq-compatible results file.
///
/// Format matches `rnaseq_qc_results.txt` expected by MultiQC.
pub fn write_qualimap_results(
    result: &GenebodyCoverageResult,
    bam_path: &str,
    total_alignments: u64,
    secondary: u64,
    output_path: &Path,
) -> Result<()> {
    let mut f = std::fs::File::create(output_path)
        .with_context(|| format!("Failed to create {}", output_path.display()))?;

    let reads_aligned = total_alignments - secondary;
    let total_bases = result.exonic + result.intronic + result.intergenic;

    let pct = |n: u64, d: u64| -> String {
        if d == 0 {
            "0".to_string()
        } else {
            format!("{:.2}", n as f64 / d as f64 * 100.0)
        }
    };

    writeln!(f, ">>>>>>> Input")?;
    writeln!(f)?;
    writeln!(f, "     bam file = {}", bam_path)?;
    writeln!(f)?;

    writeln!(f, ">>>>>>> Reads alignment")?;
    writeln!(f)?;
    writeln!(
        f,
        "     reads aligned = {}",
        format_with_commas(reads_aligned)
    )?;
    writeln!(
        f,
        "     total alignments = {}",
        format_with_commas(total_alignments)
    )?;
    writeln!(
        f,
        "     secondary alignments = {}",
        format_with_commas(secondary)
    )?;
    writeln!(
        f,
        "     non-unique alignments = {}",
        format_with_commas(secondary)
    )?;
    writeln!(
        f,
        "     aligned to genes = {}",
        format_with_commas(result.aligned_to_genes)
    )?;
    writeln!(
        f,
        "     ambiguous alignments = {}",
        format_with_commas(result.ambiguous)
    )?;
    writeln!(
        f,
        "     no feature assigned = {}",
        format_with_commas(result.no_feature)
    )?;
    writeln!(f)?;

    writeln!(f, ">>>>>>> Reads genomic origin")?;
    writeln!(f)?;
    writeln!(
        f,
        "     exonic = {} ({}%)",
        format_with_commas(result.exonic),
        pct(result.exonic, total_bases)
    )?;
    writeln!(
        f,
        "     intronic = {} ({}%)",
        format_with_commas(result.intronic),
        pct(result.intronic, total_bases)
    )?;
    writeln!(
        f,
        "     intergenic = {} ({}%)",
        format_with_commas(result.intergenic),
        pct(result.intergenic, total_bases)
    )?;
    writeln!(
        f,
        "     overlapping exon = {} ({}%)",
        format_with_commas(result.overlapping_exon),
        pct(result.overlapping_exon, result.aligned_to_genes)
    )?;
    writeln!(f)?;

    writeln!(f, ">>>>>>> Transcript coverage profile")?;
    writeln!(f)?;
    writeln!(f, "     5' bias = {:.2}", result.bias_5prime)?;
    writeln!(f, "     3' bias = {:.2}", result.bias_3prime)?;
    writeln!(f, "     5'-3' bias = {:.2}", result.bias_5_to_3)?;

    info!("Wrote Qualimap results to {}", output_path.display());
    Ok(())
}

// ============================================================
// Unit tests
// ============================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_with_commas() {
        assert_eq!(format_with_commas(0), "0");
        assert_eq!(format_with_commas(999), "999");
        assert_eq!(format_with_commas(1000), "1,000");
        assert_eq!(format_with_commas(1234567), "1,234,567");
    }

    #[test]
    fn test_median() {
        assert!((crate::io::median(&[1.0, 2.0, 3.0]) - 2.0).abs() < 1e-10);
        assert!((crate::io::median(&[1.0, 2.0, 3.0, 4.0]) - 2.5).abs() < 1e-10);
        assert!((crate::io::median(&[]) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_gene_body_coverage_merge() {
        let mut a = GenebodyCoverageAccum::new();
        a.coverage_bins[0] = 10.0;
        a.coverage_bins[50] = 20.0;
        a.aligned_to_genes = 5;
        a.exonic = 100;

        let mut b = GenebodyCoverageAccum::new();
        b.coverage_bins[0] = 5.0;
        b.coverage_bins[50] = 15.0;
        b.aligned_to_genes = 3;
        b.intronic = 50;

        a.merge(&b);
        assert!((a.coverage_bins[0] - 15.0).abs() < 1e-10);
        assert!((a.coverage_bins[50] - 35.0).abs() < 1e-10);
        assert_eq!(a.aligned_to_genes, 8);
        assert_eq!(a.exonic, 100);
        assert_eq!(a.intronic, 50);
    }

    #[test]
    fn test_bias_calculation() {
        let mut accum = GenebodyCoverageAccum::new();
        // Uniform coverage
        for bin in &mut accum.coverage_bins {
            *bin = 100.0;
        }
        let result = accum.into_result();
        assert!((result.bias_5prime - 1.0).abs() < 1e-10);
        assert!((result.bias_3prime - 1.0).abs() < 1e-10);
        assert!((result.bias_5_to_3 - 1.0).abs() < 1e-10);
    }
}
