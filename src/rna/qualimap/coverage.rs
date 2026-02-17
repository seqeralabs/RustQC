//! Per-transcript and per-gene coverage tracking for Qualimap RNA-Seq QC.
//!
//! Maintains per-base coverage arrays for both individual transcripts and
//! merged gene models. Qualimap uses Picard's Gene.Transcript model which
//! merges all exons from all isoforms into a single non-redundant exon set
//! per gene. Coverage for bias/profile is tracked on this merged model.

use std::collections::HashMap;

use super::index::{MergedGeneModel, TranscriptInfo};

// ============================================================
// Per-transcript coverage accumulator
// ============================================================

/// Per-transcript per-base coverage tracker.
///
/// Stores a sparse map of transcript flat-index → per-base depth array.
/// Coverage is added for ALL transcripts of a gene when a read is assigned.
/// The depth array uses `i32` to allow efficient increment without overflow
/// concerns for typical RNA-Seq depths.
#[derive(Debug, Clone, Default)]
pub struct TranscriptCoverage {
    /// Per-base coverage arrays, keyed by flat transcript index in `QualimapIndex.transcripts`.
    /// Only transcripts that have been touched (at least one read assigned to their gene)
    /// get an entry here. The `Vec<i32>` has length = `transcript.length`.
    coverage: HashMap<u32, Vec<i32>>,
}

impl TranscriptCoverage {
    /// Create a new empty coverage tracker.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add coverage for all transcripts of a gene.
    ///
    /// For each transcript of the gene, maps the read's M-blocks (genomic coordinates,
    /// 0-based half-open) to transcript-relative positions and increments the depth array.
    ///
    /// # Arguments
    /// * `gene_idx` - Index of the assigned gene.
    /// * `m_blocks` - Read's M-only aligned blocks in 0-based half-open genomic coords.
    /// * `transcripts` - Slice of `TranscriptInfo` for this gene (from `QualimapIndex::gene_transcripts`).
    /// * `tx_flat_offset` - Flat index offset: the index of the first transcript of this gene
    ///   in the global transcripts vec.
    pub fn add_coverage(
        &mut self,
        m_blocks: &[(i32, i32)],
        transcripts: &[TranscriptInfo],
        tx_flat_offset: u32,
    ) {
        for (local_idx, tx_info) in transcripts.iter().enumerate() {
            let flat_idx = tx_flat_offset + local_idx as u32;
            let tx_len = tx_info.exonic_length as usize;
            if tx_len == 0 {
                continue;
            }

            // Get or create the depth array for this transcript.
            // Length = exonic_length (sum of exon lengths), matching Picard's
            // Gene.Transcript which stores coverage in exonic coordinates.
            let depth = self
                .coverage
                .entry(flat_idx)
                .or_insert_with(|| vec![0i32; tx_len]);

            // Map each M-block to exonic coordinates.
            // For each genomic position in an M-block that falls within an exon,
            // compute the exonic offset (position within concatenated exons)
            // and increment the depth at that position.
            // This matches Picard's Gene.Transcript.addCoverageCounts()
            // which calls getTranscriptCoordinate() per genomic position.
            for &(block_start, block_end) in m_blocks {
                // For each exon, compute overlap with this M-block.
                // exonic_offset tracks cumulative exonic bases from prior exons.
                let mut exonic_offset = 0usize;
                for &(exon_start, exon_end) in &tx_info.exons {
                    let exon_len = (exon_end - exon_start) as usize;

                    // Compute overlap between M-block [block_start, block_end)
                    // and exon [exon_start, exon_end) (both 0-based half-open)
                    let overlap_start = block_start.max(exon_start);
                    let overlap_end = block_end.min(exon_end);

                    if overlap_start < overlap_end {
                        // This M-block overlaps this exon.
                        // Map overlap to exonic coordinates:
                        // rel_start = exonic_offset + (overlap_start - exon_start)
                        let rel_start = exonic_offset + (overlap_start - exon_start) as usize;
                        let rel_end =
                            (exonic_offset + (overlap_end - exon_start) as usize).min(tx_len);
                        if rel_start < tx_len && rel_start < rel_end {
                            for d in &mut depth[rel_start..rel_end] {
                                *d += 1;
                            }
                        }
                    }

                    exonic_offset += exon_len;
                }
            }
        }
    }

    /// Merge another `TranscriptCoverage` into this one.
    ///
    /// For transcripts present in both, adds depth arrays element-wise.
    /// For transcripts only in `other`, moves them into `self`.
    pub fn merge(&mut self, other: TranscriptCoverage) {
        for (tx_idx, other_depth) in other.coverage {
            let entry = self.coverage.entry(tx_idx);
            match entry {
                std::collections::hash_map::Entry::Occupied(mut e) => {
                    let self_depth = e.get_mut();
                    for (i, &val) in other_depth.iter().enumerate() {
                        if i < self_depth.len() {
                            self_depth[i] += val;
                        }
                    }
                }
                std::collections::hash_map::Entry::Vacant(e) => {
                    e.insert(other_depth);
                }
            }
        }
    }

    /// Get the depth array for a transcript (by flat index).
    #[allow(dead_code)]
    pub fn get(&self, tx_flat_idx: u32) -> Option<&[i32]> {
        self.coverage.get(&tx_flat_idx).map(|v| v.as_slice())
    }

    /// Number of transcripts with coverage data.
    #[allow(dead_code)]
    pub fn num_transcripts_with_coverage(&self) -> usize {
        self.coverage.len()
    }

    /// Iterate over all transcripts with coverage data.
    pub fn iter(&self) -> impl Iterator<Item = (u32, &[i32])> {
        self.coverage.iter().map(|(&idx, v)| (idx, v.as_slice()))
    }

    /// Compute mean coverage for a transcript.
    #[allow(dead_code)]
    pub fn mean_coverage(&self, tx_flat_idx: u32) -> f64 {
        match self.coverage.get(&tx_flat_idx) {
            Some(depth) if !depth.is_empty() => {
                let sum: i64 = depth.iter().map(|&d| d as i64).sum();
                sum as f64 / depth.len() as f64
            }
            _ => 0.0,
        }
    }
}

// ============================================================
// Merged gene coverage accumulator
// ============================================================

/// Per-gene per-base coverage tracker using merged gene models.
///
/// Qualimap (via Picard's Gene.Transcript) merges all exons from all transcripts
/// of a gene into a single non-redundant exon set. Coverage is tracked on this
/// merged model, which is used for bias calculation and coverage profiles.
/// This produces more accurate bias values than per-transcript tracking because
/// coverage from reads mapping to different isoform exons is unified.
#[derive(Debug, Clone, Default)]
pub struct MergedGeneCoverage {
    /// Per-base coverage arrays, keyed by gene index.
    /// The `Vec<i32>` has length = `MergedGeneModel.exonic_length`.
    coverage: HashMap<u32, Vec<i32>>,
}

impl MergedGeneCoverage {
    /// Create a new empty merged gene coverage tracker.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add coverage for a gene using its merged model.
    ///
    /// Maps read M-blocks (genomic coordinates, 0-based half-open) to positions
    /// in the merged exonic coordinate space and increments the depth array.
    ///
    /// # Arguments
    /// * `gene_idx` - Index of the assigned gene.
    /// * `m_blocks` - Read's M-only aligned blocks in 0-based half-open genomic coords.
    /// * `model` - The merged gene model with non-redundant exon intervals.
    pub fn add_coverage(
        &mut self,
        gene_idx: u32,
        m_blocks: &[(i32, i32)],
        model: &MergedGeneModel,
    ) {
        let total_len = model.exonic_length as usize;
        if total_len == 0 {
            return;
        }

        let depth = self
            .coverage
            .entry(gene_idx)
            .or_insert_with(|| vec![0i32; total_len]);

        // Map each M-block to merged exonic coordinates.
        for &(block_start, block_end) in m_blocks {
            let mut exonic_offset = 0usize;
            for &(exon_start, exon_end) in &model.exons {
                let exon_len = (exon_end - exon_start) as usize;

                // Compute overlap between M-block and this merged exon
                let overlap_start = block_start.max(exon_start);
                let overlap_end = block_end.min(exon_end);

                if overlap_start < overlap_end {
                    let rel_start = exonic_offset + (overlap_start - exon_start) as usize;
                    let rel_end =
                        (exonic_offset + (overlap_end - exon_start) as usize).min(total_len);
                    if rel_start < total_len && rel_start < rel_end {
                        for d in &mut depth[rel_start..rel_end] {
                            *d += 1;
                        }
                    }
                }

                exonic_offset += exon_len;
            }
        }
    }

    /// Merge another `MergedGeneCoverage` into this one.
    ///
    /// For genes present in both, adds depth arrays element-wise.
    /// For genes only in `other`, moves them into `self`.
    #[allow(dead_code)]
    pub fn merge(&mut self, other: MergedGeneCoverage) {
        for (gene_idx, other_depth) in other.coverage {
            let entry = self.coverage.entry(gene_idx);
            match entry {
                std::collections::hash_map::Entry::Occupied(mut e) => {
                    let self_depth = e.get_mut();
                    for (i, &val) in other_depth.iter().enumerate() {
                        if i < self_depth.len() {
                            self_depth[i] += val;
                        }
                    }
                }
                std::collections::hash_map::Entry::Vacant(e) => {
                    e.insert(other_depth);
                }
            }
        }
    }

    /// Get the depth array for a gene (by gene index).
    #[allow(dead_code)]
    pub fn get(&self, gene_idx: u32) -> Option<&[i32]> {
        self.coverage.get(&gene_idx).map(|v| v.as_slice())
    }

    /// Iterate over all genes with coverage data.
    #[allow(dead_code)]
    pub fn iter(&self) -> impl Iterator<Item = (u32, &[i32])> {
        self.coverage.iter().map(|(&idx, v)| (idx, v.as_slice()))
    }

    /// Number of genes that have any coverage recorded.
    #[allow(dead_code)]
    pub fn num_genes_with_coverage(&self) -> usize {
        self.coverage.len()
    }
}

// ============================================================
// Unit tests
// ============================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn make_tx_info(start: i32, end: i32) -> TranscriptInfo {
        TranscriptInfo {
            gene_idx: 0,
            transcript_idx: 0,
            gene_id: "g1".to_string(),
            transcript_id: "t1".to_string(),
            chrom: "chr1".to_string(),
            start,
            end,
            strand: '+',
            length: (end - start) as u32,
            exons: vec![(start, end)],
            exonic_length: (end - start) as u32,
        }
    }

    #[test]
    fn test_add_coverage_single_block() {
        let mut cov = TranscriptCoverage::new();
        let tx = make_tx_info(100, 200);
        // M-block: 120..150 (0-based half-open)
        cov.add_coverage(&[(120, 150)], &[tx], 0);

        let depth = cov.get(0).unwrap();
        assert_eq!(depth.len(), 100);
        // Positions 20..50 (relative) should have depth 1
        assert_eq!(depth[19], 0);
        assert_eq!(depth[20], 1);
        assert_eq!(depth[49], 1);
        assert_eq!(depth[50], 0);
    }

    #[test]
    fn test_add_coverage_clipping() {
        let mut cov = TranscriptCoverage::new();
        let tx = make_tx_info(100, 200);
        // M-block extends beyond transcript: 180..220
        cov.add_coverage(&[(180, 220)], &[tx], 0);

        let depth = cov.get(0).unwrap();
        // Only 80..100 (relative) should be covered
        assert_eq!(depth[79], 0);
        assert_eq!(depth[80], 1);
        assert_eq!(depth[99], 1);
    }

    #[test]
    fn test_merge() {
        let mut cov1 = TranscriptCoverage::new();
        let tx = make_tx_info(100, 200);
        cov1.add_coverage(&[(120, 150)], &[tx.clone()], 0);

        let mut cov2 = TranscriptCoverage::new();
        cov2.add_coverage(&[(130, 160)], &[tx], 0);

        cov1.merge(cov2);
        let depth = cov1.get(0).unwrap();
        // 120..130 relative 20..30: depth 1
        assert_eq!(depth[25], 1);
        // 130..150 relative 30..50: depth 2
        assert_eq!(depth[35], 2);
        // 150..160 relative 50..60: depth 1
        assert_eq!(depth[55], 1);
    }

    #[test]
    fn test_mean_coverage() {
        let mut cov = TranscriptCoverage::new();
        let tx = make_tx_info(0, 10);
        // Cover all 10 bases with depth 1
        cov.add_coverage(&[(0, 10)], &[tx], 0);
        let mean = cov.mean_coverage(0);
        assert!((mean - 1.0).abs() < 1e-10);

        // Non-existent transcript
        assert!((cov.mean_coverage(99) - 0.0).abs() < 1e-10);
    }

    // --------------------------------------------------------
    // MergedGeneCoverage tests
    // --------------------------------------------------------

    #[test]
    fn test_merged_gene_single_exon() {
        let model = MergedGeneModel {
            gene_idx: 0,
            strand: '+',
            exons: vec![(100, 200)],
            exonic_length: 100,
        };
        let mut cov = MergedGeneCoverage::new();
        cov.add_coverage(0, &[(120, 150)], &model);

        let depth = cov.get(0).unwrap();
        assert_eq!(depth.len(), 100);
        assert_eq!(depth[19], 0);
        assert_eq!(depth[20], 1);
        assert_eq!(depth[49], 1);
        assert_eq!(depth[50], 0);
    }

    #[test]
    fn test_merged_gene_two_exons() {
        // Gene with two merged exons: (100,200) and (300,400)
        // Total exonic length = 200, offset: first exon 0..100, second 100..200
        let model = MergedGeneModel {
            gene_idx: 0,
            strand: '+',
            exons: vec![(100, 200), (300, 400)],
            exonic_length: 200,
        };
        let mut cov = MergedGeneCoverage::new();
        // M-block in second exon: 310..330 (genomic) → offset 100 + 10 = 110..130 (exonic)
        cov.add_coverage(0, &[(310, 330)], &model);

        let depth = cov.get(0).unwrap();
        assert_eq!(depth.len(), 200);
        assert_eq!(depth[109], 0);
        assert_eq!(depth[110], 1);
        assert_eq!(depth[129], 1);
        assert_eq!(depth[130], 0);
    }

    #[test]
    fn test_merged_gene_overlapping_transcripts() {
        // Simulates what happens when two transcripts have overlapping exons:
        // tx1: exons (100,200) and (300,400)
        // tx2: exons (150,250) and (350,400)
        // Merged: (100,250) and (300,400) → total 250
        let model = MergedGeneModel::from_transcripts(
            0,
            '+',
            &[(100, 200), (300, 400), (150, 250), (350, 400)],
        );
        assert_eq!(model.exons, vec![(100, 250), (300, 400)]);
        assert_eq!(model.exonic_length, 250);

        let mut cov = MergedGeneCoverage::new();
        // M-block spanning the merged first exon: 180..220 → offset 80..120
        cov.add_coverage(0, &[(180, 220)], &model);

        let depth = cov.get(0).unwrap();
        assert_eq!(depth.len(), 250);
        assert_eq!(depth[79], 0);
        assert_eq!(depth[80], 1);
        assert_eq!(depth[119], 1);
        assert_eq!(depth[120], 0);
    }

    #[test]
    fn test_merged_gene_merge() {
        let model = MergedGeneModel {
            gene_idx: 0,
            strand: '+',
            exons: vec![(100, 200)],
            exonic_length: 100,
        };

        let mut cov1 = MergedGeneCoverage::new();
        cov1.add_coverage(0, &[(120, 150)], &model);

        let mut cov2 = MergedGeneCoverage::new();
        cov2.add_coverage(0, &[(130, 160)], &model);

        cov1.merge(cov2);
        let depth = cov1.get(0).unwrap();
        // 120..130 → offset 20..30: depth 1
        assert_eq!(depth[25], 1);
        // 130..150 → offset 30..50: depth 2
        assert_eq!(depth[35], 2);
        // 150..160 → offset 50..60: depth 1
        assert_eq!(depth[55], 1);
    }
}
