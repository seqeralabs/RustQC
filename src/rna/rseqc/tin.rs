//! TIN (Transcript Integrity Number) analysis.
//!
//! Measures transcript integrity via Shannon entropy of read coverage
//! uniformity across sampled exonic positions. Reimplementation of
//! RSeQC's `tin.py` tool.

use rust_htslib::bam;
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use indexmap::IndexMap;
use log::debug;

use crate::gtf::Gene;
use crate::io as rustqc_io;

// ===================================================================
// Data structures
// ===================================================================

/// A single transcript's sampling metadata.
#[derive(Debug, Clone)]
pub struct TranscriptSampling {
    /// Gene identifier (gene_id from GTF or BED name).
    pub gene_id: String,
    /// Chromosome name (original case, for output).
    pub chrom: String,
    /// Chromosome name uppercased (for lookup matching).
    pub chrom_upper: String,
    /// Transcript start position (0-based).
    pub tx_start: u64,
    /// Transcript end position (0-based, exclusive).
    pub tx_end: u64,
    /// Sorted exon blocks as (start, end) in 0-based half-open coords.
    #[allow(dead_code)]
    pub exon_regions: Vec<(u64, u64)>,
    /// Total exonic bases across all exon blocks.
    #[allow(dead_code)]
    pub exon_length: u64,
    /// Genomic positions sampled within exonic regions (sorted).
    pub sampled_positions: Vec<u64>,
}

/// Index for fast lookup of sampled positions during BAM processing.
///
/// Maps chromosome → sorted list of (genomic_pos, transcript_idx, slot_idx)
/// for binary search during per-read dispatch.
#[derive(Debug, Default)]
pub struct TinIndex {
    /// All transcripts with their sampling metadata.
    pub transcripts: Vec<TranscriptSampling>,
    /// Per-chromosome sorted position entries for binary search.
    /// Each entry: (genomic_position, transcript_index, slot_index).
    pub chrom_positions: HashMap<String, Vec<(u64, u32, u32)>>,
    /// Per-chromosome transcript spans for overlap pre-filtering.
    /// Each entry: (tx_start, tx_end, transcript_index), sorted by tx_start.
    pub chrom_spans: HashMap<String, Vec<(u64, u64, u32)>>,
}

/// Result of TIN computation for a single transcript.
#[derive(Debug)]
pub struct TinResult {
    /// Gene identifier.
    pub gene_id: String,
    /// Chromosome name.
    pub chrom: String,
    /// Transcript start (0-based).
    pub tx_start: u64,
    /// Transcript end (0-based, exclusive).
    pub tx_end: u64,
    /// Computed TIN score (0-100), or NaN if below coverage threshold.
    pub tin: f64,
}

/// Collection of all TIN results for a BAM file.
#[derive(Debug)]
pub struct TinResults {
    /// Per-transcript TIN scores.
    pub transcripts: Vec<TinResult>,
}

impl TinResults {
    /// Number of transcripts with computed TIN scores.
    pub fn len(&self) -> usize {
        self.transcripts.len()
    }

    /// Returns true if there are no results.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.transcripts.is_empty()
    }
}

// ===================================================================
// TIN index construction
// ===================================================================

impl TinIndex {
    /// Build TIN index from GTF gene annotations.
    ///
    /// Uses the longest transcript per gene (by exonic bases) as the
    /// representative transcript.
    pub fn from_genes(genes: &IndexMap<String, Gene>, sample_size: usize) -> Self {
        let mut index = TinIndex::default();

        for gene in genes.values() {
            // Find the longest transcript by total exonic bases
            let best_tx = gene.transcripts.iter().max_by_key(|tx| {
                tx.exons
                    .iter()
                    .map(|(s, e)| e.saturating_sub(*s))
                    .sum::<u64>()
            });

            if let Some(tx) = best_tx {
                let exon_regions: Vec<(u64, u64)> = tx
                    .exons
                    .iter()
                    .map(|&(s, e)| (s.saturating_sub(1), e)) // GTF 1-based to 0-based half-open
                    .collect();
                let exon_length: u64 = exon_regions.iter().map(|(s, e)| e - s).sum();

                if exon_length == 0 {
                    continue;
                }

                let sampled = sample_exonic_positions(&exon_regions, sample_size);
                if sampled.is_empty() {
                    continue;
                }

                let chrom = gene.chrom.clone();
                let chrom_upper = chrom.to_uppercase();
                let tx_start = exon_regions.first().map(|r| r.0).unwrap_or(0);
                let tx_end = exon_regions.last().map(|r| r.1).unwrap_or(0);

                index.add_transcript(TranscriptSampling {
                    gene_id: gene.gene_id.clone(),
                    chrom,
                    chrom_upper,
                    tx_start,
                    tx_end,
                    exon_regions,
                    exon_length,
                    sampled_positions: sampled,
                });
            }
        }

        index.build();
        index
    }

    /// Build TIN index from BED12 file.
    pub fn from_bed(bed_path: &str, sample_size: usize) -> Result<Self> {
        let content = rustqc_io::read_to_string(bed_path)
            .with_context(|| format!("Failed to read BED file: {bed_path}"))?;

        let mut index = TinIndex::default();

        for line in content.lines() {
            if line.starts_with('#') || line.starts_with("track") || line.is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 12 {
                continue;
            }

            let chrom = fields[0].to_string();
            let chrom_upper = chrom.to_uppercase();
            let tx_start: u64 = fields[1].parse().unwrap_or(0);
            let tx_end: u64 = fields[2].parse().unwrap_or(0);
            let name = fields[3].to_string();
            let block_count: usize = fields[9].parse().unwrap_or(0);

            if block_count == 0 {
                continue;
            }

            let block_sizes: Vec<u64> = fields[10]
                .split(',')
                .filter(|s| !s.is_empty())
                .filter_map(|s| s.parse().ok())
                .collect();
            let block_starts: Vec<u64> = fields[11]
                .split(',')
                .filter(|s| !s.is_empty())
                .filter_map(|s| s.parse().ok())
                .collect();

            if block_sizes.len() != block_count || block_starts.len() != block_count {
                continue;
            }

            let exon_regions: Vec<(u64, u64)> = block_starts
                .iter()
                .zip(block_sizes.iter())
                .map(|(&start, &size)| (tx_start + start, tx_start + start + size))
                .collect();

            let exon_length: u64 = exon_regions.iter().map(|(s, e)| e - s).sum();
            if exon_length == 0 {
                continue;
            }

            let sampled = sample_exonic_positions(&exon_regions, sample_size);
            if sampled.is_empty() {
                continue;
            }

            index.add_transcript(TranscriptSampling {
                gene_id: name,
                chrom,
                chrom_upper,
                tx_start,
                tx_end,
                exon_regions,
                exon_length,
                sampled_positions: sampled,
            });
        }

        index.build();
        Ok(index)
    }

    /// Add a transcript to the index (before calling `build()`).
    fn add_transcript(&mut self, tx: TranscriptSampling) {
        let tx_idx = self.transcripts.len() as u32;
        let chrom_upper = tx.chrom_upper.clone();
        let tx_start = tx.tx_start;
        let tx_end = tx.tx_end;

        // Add sampled positions to per-chromosome position list
        for (slot_idx, &pos) in tx.sampled_positions.iter().enumerate() {
            self.chrom_positions
                .entry(chrom_upper.clone())
                .or_default()
                .push((pos, tx_idx, slot_idx as u32));
        }

        // Add transcript span
        self.chrom_spans
            .entry(chrom_upper)
            .or_default()
            .push((tx_start, tx_end, tx_idx));

        self.transcripts.push(tx);
    }

    /// Sort per-chromosome data for binary search.
    fn build(&mut self) {
        for positions in self.chrom_positions.values_mut() {
            positions.sort_unstable();
        }
        for spans in self.chrom_spans.values_mut() {
            spans.sort_unstable();
        }
        debug!(
            "TIN index: {} transcripts, {} chromosomes",
            self.transcripts.len(),
            self.chrom_positions.len()
        );
    }
}

/// Sample `n` equally-spaced positions within the exonic regions.
///
/// Returns genomic coordinates of sampled positions (0-based).
fn sample_exonic_positions(exon_regions: &[(u64, u64)], n: usize) -> Vec<u64> {
    let total_bases: u64 = exon_regions.iter().map(|(s, e)| e - s).sum();
    if total_bases == 0 {
        return Vec::new();
    }

    let n = n as u64;
    let num_slots = n.min(total_bases);
    if num_slots == 0 {
        return Vec::new();
    }

    let mut positions = Vec::with_capacity(num_slots as usize);

    for i in 0..num_slots {
        // Map slot i to an offset within the concatenated exonic region
        let offset = if num_slots == 1 {
            total_bases / 2
        } else {
            (i * (total_bases - 1)) / (num_slots - 1)
        };

        // Convert offset to genomic coordinate
        let mut remaining = offset;
        for &(start, end) in exon_regions {
            let block_len = end - start;
            if remaining < block_len {
                positions.push(start + remaining);
                break;
            }
            remaining -= block_len;
        }
    }

    positions
}

// ===================================================================
// TIN accumulator
// ===================================================================

/// Per-read accumulator for TIN computation.
///
/// Tracks coverage at sampled positions and read start counts
/// per transcript during the BAM pass.
#[derive(Debug)]
pub struct TinAccum {
    /// Per-transcript per-slot coverage counts.
    /// Indexed as `coverage[tx_idx][slot_idx]`.
    pub coverage: Vec<Vec<u32>>,
    /// Per-transcript read start counts.
    pub read_starts: Vec<u32>,
    /// Number of sampled slots per transcript.
    #[allow(dead_code)]
    pub n_samples: Vec<u32>,
    /// Minimum MAPQ threshold.
    pub mapq_cut: u8,
    /// Minimum coverage threshold.
    pub min_cov: u32,
}

impl TinAccum {
    /// Create a new TIN accumulator for the given index.
    pub fn new(index: &TinIndex, mapq_cut: u8, min_cov: u32) -> Self {
        let n_transcripts = index.transcripts.len();
        let mut coverage = Vec::with_capacity(n_transcripts);
        let mut n_samples = Vec::with_capacity(n_transcripts);

        for tx in &index.transcripts {
            let n = tx.sampled_positions.len();
            coverage.push(vec![0u32; n]);
            n_samples.push(n as u32);
        }

        TinAccum {
            coverage,
            read_starts: vec![0u32; n_transcripts],
            n_samples,
            mapq_cut,
            min_cov,
        }
    }

    /// Process a single BAM record for TIN accumulation.
    ///
    /// Updates coverage at sampled positions that overlap the read's
    /// aligned blocks, and tracks read start positions.
    pub fn process_read(&mut self, record: &bam::Record, chrom_upper: &str, index: &TinIndex) {
        let flags = record.flags();

        // Skip unmapped, QC-fail, secondary, supplementary, duplicate
        if flags & 0x4 != 0
            || flags & 0x200 != 0
            || flags & 0x100 != 0
            || flags & 0x800 != 0
            || flags & 0x400 != 0
        {
            return;
        }

        // MAPQ filter
        if record.mapq() < self.mapq_cut {
            return;
        }

        // Get aligned blocks from CIGAR
        let blocks = get_aligned_blocks(record);
        if blocks.is_empty() {
            return;
        }

        let read_start = blocks[0].0;

        // Get chromosome positions for this chromosome
        let chrom_positions = match index.chrom_positions.get(chrom_upper) {
            Some(p) => p,
            None => return,
        };

        // Get transcript spans for this chromosome
        let chrom_spans = match index.chrom_spans.get(chrom_upper) {
            Some(s) => s,
            None => return,
        };

        // Find overlapping transcripts via read start
        let read_end = blocks.last().map(|b| b.1).unwrap_or(read_start);

        // Binary search for first transcript that could overlap
        let span_start = chrom_spans.partition_point(|s| s.1 <= read_start);

        for &(tx_start, tx_end, tx_idx) in &chrom_spans[span_start..] {
            if tx_start >= read_end {
                break;
            }
            // Read start falls within this transcript's exonic region
            if read_start >= tx_start && read_start < tx_end {
                self.read_starts[tx_idx as usize] += 1;
            }
        }

        // For each aligned block, find sampled positions that fall within it
        for &(block_start, block_end) in &blocks {
            // Binary search for first position >= block_start
            let pos_start = chrom_positions.partition_point(|p| p.0 < block_start);

            for &(pos, tx_idx, slot_idx) in &chrom_positions[pos_start..] {
                if pos >= block_end {
                    break;
                }
                // This sampled position is covered by this aligned block
                self.coverage[tx_idx as usize][slot_idx as usize] += 1;
            }
        }
    }

    /// Merge another accumulator into this one.
    pub fn merge(&mut self, other: TinAccum) {
        for (i, other_cov) in other.coverage.into_iter().enumerate() {
            for (j, count) in other_cov.into_iter().enumerate() {
                self.coverage[i][j] += count;
            }
            self.read_starts[i] += other.read_starts[i];
        }
    }

    /// Compute TIN scores for all transcripts.
    pub fn into_result(self, index: &TinIndex) -> TinResults {
        let mut transcripts = Vec::with_capacity(index.transcripts.len());

        for (i, tx) in index.transcripts.iter().enumerate() {
            let read_starts = self.read_starts[i];
            let coverage = &self.coverage[i];

            // Require minimum coverage threshold
            if read_starts < self.min_cov {
                transcripts.push(TinResult {
                    gene_id: tx.gene_id.clone(),
                    chrom: tx.chrom.clone(),
                    tx_start: tx.tx_start,
                    tx_end: tx.tx_end,
                    tin: 0.0,
                });
                continue;
            }

            let tin = compute_tin(coverage);

            transcripts.push(TinResult {
                gene_id: tx.gene_id.clone(),
                chrom: tx.chrom.clone(),
                tx_start: tx.tx_start,
                tx_end: tx.tx_end,
                tin,
            });
        }

        TinResults { transcripts }
    }
}

/// Compute TIN score from coverage vector using Shannon entropy.
///
/// TIN = 100 * exp(H) / n_nonzero
/// where H = -sum(P_i * ln(P_i)) for non-zero positions.
fn compute_tin(coverage: &[u32]) -> f64 {
    let nonzero: Vec<f64> = coverage
        .iter()
        .filter(|&&c| c > 0)
        .map(|&c| c as f64)
        .collect();

    let n = nonzero.len();
    if n <= 1 {
        return 0.0;
    }

    let total: f64 = nonzero.iter().sum();
    if total == 0.0 {
        return 0.0;
    }

    // Shannon entropy
    let entropy: f64 = nonzero
        .iter()
        .map(|&c| {
            let p = c / total;
            -p * p.ln()
        })
        .sum();

    // TIN = 100 * exp(H) / n_nonzero
    100.0 * entropy.exp() / n as f64
}

/// Get aligned blocks from a BAM record's CIGAR string.
///
/// Returns sorted (start, end) pairs of aligned blocks (0-based half-open).
fn get_aligned_blocks(record: &bam::Record) -> Vec<(u64, u64)> {
    use rust_htslib::bam::record::Cigar;

    let mut blocks = Vec::new();
    let mut pos = record.pos() as u64;

    for op in record.cigar().iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                blocks.push((pos, pos + *len as u64));
                pos += *len as u64;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                pos += *len as u64;
            }
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    blocks
}

// ===================================================================
// Output writers
// ===================================================================

/// Write per-transcript TIN scores in RSeQC format.
///
/// Format: TSV with header geneID\tchrom\ttx_start\ttx_end\tTIN
pub fn write_tin(results: &TinResults, output_path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(output_path)
        .with_context(|| format!("Failed to create TIN output: {}", output_path.display()))?;

    writeln!(f, "geneID\tchrom\ttx_start\ttx_end\tTIN")?;

    for r in &results.transcripts {
        if r.tin == 0.0 {
            continue; // Skip transcripts below coverage threshold
        }
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{:.2}",
            r.gene_id, r.chrom, r.tx_start, r.tx_end, r.tin
        )?;
    }

    Ok(())
}

/// Write TIN summary statistics.
///
/// Format: Bam_file\tTIN(mean)\tTIN(median)\tTIN(stdev)
pub fn write_tin_summary(results: &TinResults, bam_name: &str, output_path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(output_path)
        .with_context(|| format!("Failed to create TIN summary: {}", output_path.display()))?;

    let scores: Vec<f64> = results
        .transcripts
        .iter()
        .filter(|r| r.tin > 0.0)
        .map(|r| r.tin)
        .collect();

    let (mean, median, stdev) = if scores.is_empty() {
        (0.0, 0.0, 0.0)
    } else {
        let n = scores.len() as f64;
        let mean = scores.iter().sum::<f64>() / n;

        let median = crate::io::median(&scores);

        let variance = scores.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;
        let stdev = variance.sqrt();

        (mean, median, stdev)
    };

    writeln!(f, "Bam_file\tTIN(mean)\tTIN(median)\tTIN(stdev)")?;
    writeln!(f, "{}\t{:.2}\t{:.2}\t{:.2}", bam_name, mean, median, stdev)?;

    Ok(())
}

// ===================================================================
// Unit tests
// ===================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sample_exonic_positions() {
        // Single exon of 10 bases
        let exons = vec![(100, 110)];
        let positions = sample_exonic_positions(&exons, 5);
        assert_eq!(positions.len(), 5);
        // Positions should span the exonic region
        assert!(positions[0] >= 100 && positions[0] < 110);
        assert!(positions[4] >= 100 && positions[4] < 110);
        // First and last should be near boundaries
        assert_eq!(positions[0], 100);
        assert_eq!(positions[4], 109);
    }

    #[test]
    fn test_sample_multi_exon() {
        // Two exons
        let exons = vec![(100, 105), (200, 205)];
        let positions = sample_exonic_positions(&exons, 10);
        assert_eq!(positions.len(), 10);
        // First position in first exon, last in second exon
        assert_eq!(positions[0], 100);
        assert_eq!(positions[9], 204);
    }

    #[test]
    fn test_sample_tiny_exon() {
        // Exon smaller than sample size
        let exons = vec![(100, 103)];
        let positions = sample_exonic_positions(&exons, 100);
        assert_eq!(positions.len(), 3); // Clamped to exon length
    }

    #[test]
    fn test_compute_tin_uniform() {
        // Uniform coverage → max TIN
        let cov = vec![10, 10, 10, 10, 10, 10, 10, 10, 10, 10];
        let tin = compute_tin(&cov);
        assert!(
            (tin - 100.0).abs() < 0.01,
            "Uniform coverage TIN should be ~100: {tin}"
        );
    }

    #[test]
    fn test_compute_tin_degraded() {
        // 5' degradation: high at start, low at end
        let cov = vec![100, 80, 60, 40, 20, 10, 5, 2, 1, 1];
        let tin = compute_tin(&cov);
        assert!(
            tin > 0.0 && tin < 100.0,
            "Degraded TIN should be between 0 and 100: {tin}"
        );
    }

    #[test]
    fn test_compute_tin_all_zero() {
        let cov = vec![0, 0, 0, 0, 0];
        let tin = compute_tin(&cov);
        assert_eq!(tin, 0.0);
    }

    #[test]
    fn test_compute_tin_single_nonzero() {
        // Only one position covered → TIN = 0
        let cov = vec![0, 0, 100, 0, 0];
        let tin = compute_tin(&cov);
        assert_eq!(tin, 0.0, "Single position coverage should give TIN=0");
    }
}
