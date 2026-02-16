//! Per-transcript coverage tracking for Qualimap RNA-Seq QC.
//!
//! Maintains per-base coverage arrays for all transcripts of assigned genes.
//! Qualimap tracks coverage across ALL transcripts (not just the longest),
//! using genomic coordinates mapped to transcript-relative positions.

use std::collections::HashMap;

use super::index::TranscriptInfo;

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
            let tx_len = tx_info.length as usize;
            if tx_len == 0 {
                continue;
            }

            // Get or create the depth array for this transcript
            let depth = self
                .coverage
                .entry(flat_idx)
                .or_insert_with(|| vec![0i32; tx_len]);

            // Map each M-block to transcript-relative coordinates
            let tx_start = tx_info.start;
            let tx_end = tx_info.end;
            for &(block_start, block_end) in m_blocks {
                // Clip to transcript span
                let clipped_start = block_start.max(tx_start);
                let clipped_end = block_end.min(tx_end);
                if clipped_start >= clipped_end {
                    continue;
                }

                // Convert to transcript-relative position
                let rel_start = (clipped_start - tx_start) as usize;
                let rel_end = (clipped_end - tx_start) as usize;
                let rel_end = rel_end.min(tx_len);

                for d in &mut depth[rel_start..rel_end] {
                    *d += 1;
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
}
