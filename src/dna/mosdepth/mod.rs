//! mosdepth-compatible depth of coverage results.
//!
//! [`ContigDepth::from_depths`] turns one contig's per-base depth vector into
//! everything the six mosdepth outputs need, in a single pass over the vector,
//! so the depth vector can be dropped as soon as the contig is done.
//!
//! # Output formats
//!
//! These were derived from mosdepth 0.3.14 output committed under
//! `tests/expected/dna/`, not from documentation, and every rule below was
//! checked against every row of those fixtures.
//!
//! `{prefix}.mosdepth.summary.txt` carries the header
//! `chrom length bases mean min max`, one row per contig, then one
//! `{contig}_region` row per contig when windows were requested, then `total`
//! and `total_region`. `mean` is `bases / length` to two decimals.
//!
//! `{prefix}.mosdepth.global.dist.txt` and `.region.dist.txt` carry
//! `chrom depth proportion` rows in descending depth order, where `proportion`
//! is the fraction at depth **at or above** `depth`, formatted to two
//! decimals, ending at depth 0 with `1.00`. Which depths get a row is the
//! non-obvious part:
//!
//! - depths 0 through [`DIST_DENSE_MAX`] always get a row, even when no base
//!   sits at that exact depth;
//! - above that, only depths that actually occur;
//! - the maximum observed depth never gets a row.
//!
//! The global distribution is over bases and their exact depth; the region
//! distribution is over windows and their **rounded** mean depth.

use std::collections::BTreeMap;

pub mod output;

/// Highest depth that always gets a distribution row, matching the size of
/// mosdepth's internal fixed depth array.
pub const DIST_DENSE_MAX: u32 = 300;

/// A run of consecutive bases sharing one depth, as written to `per-base.bed.gz`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DepthRun {
    /// Zero-based, inclusive start.
    pub start: u64,
    /// Zero-based, exclusive end.
    pub end: u64,
    /// Depth shared by every base in the run.
    pub depth: u32,
}

/// A fixed-width window and its mean depth, as written to `regions.bed.gz`.
#[derive(Debug, Clone, PartialEq)]
pub struct WindowDepth {
    /// Zero-based, inclusive start.
    pub start: u64,
    /// Zero-based, exclusive end.
    pub end: u64,
    /// Mean depth over the window.
    pub mean: f64,
}

/// One window's per-threshold counts, as written to `thresholds.bed.gz`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ThresholdRow {
    /// Zero-based, inclusive start.
    pub start: u64,
    /// Zero-based, exclusive end.
    pub end: u64,
    /// Bases at or above each requested threshold, in the requested order.
    pub counts: Vec<u64>,
}

/// Everything the mosdepth outputs need about one contig.
#[derive(Debug, Clone)]
pub struct ContigDepth {
    /// Contig name as it appears in the alignment header.
    pub name: String,
    /// Contig length in bases.
    pub length: u64,
    /// Sum of per-base depth over the contig.
    pub total_bases: u64,
    /// Lowest per-base depth seen.
    pub min: u32,
    /// Highest per-base depth seen.
    pub max: u32,
    /// Base count per exact depth.
    pub histogram: BTreeMap<u32, u64>,
    /// Collapsed runs of equal depth.
    pub runs: Vec<DepthRun>,
    /// Per-window mean depth; empty when no window size was requested.
    pub windows: Vec<WindowDepth>,
    /// Per-window threshold counts; empty when no thresholds were requested.
    pub thresholds: Vec<ThresholdRow>,
}

impl ContigDepth {
    /// Summarise one contig's per-base depths in a single pass.
    pub fn from_depths(
        name: &str,
        depths: &[u32],
        window_size: Option<u32>,
        thresholds: &[u32],
    ) -> Self {
        let length = depths.len() as u64;
        let mut histogram: BTreeMap<u32, u64> = BTreeMap::new();
        let mut runs: Vec<DepthRun> = Vec::new();
        let mut total_bases = 0u64;

        for (i, &depth) in depths.iter().enumerate() {
            total_bases += u64::from(depth);
            *histogram.entry(depth).or_insert(0) += 1;
            match runs.last_mut() {
                Some(run) if run.depth == depth => run.end = i as u64 + 1,
                _ => runs.push(DepthRun {
                    start: i as u64,
                    end: i as u64 + 1,
                    depth,
                }),
            }
        }

        let min = depths.iter().copied().min().unwrap_or(0);
        let max = depths.iter().copied().max().unwrap_or(0);

        let (windows, threshold_rows) = match window_size {
            Some(size) if size > 0 => Self::windowed(depths, u64::from(size), thresholds),
            _ => (Vec::new(), Vec::new()),
        };

        Self {
            name: name.to_string(),
            length,
            total_bases,
            min,
            max,
            histogram,
            runs,
            windows,
            thresholds: threshold_rows,
        }
    }

    /// Split the contig into fixed-width windows, computing each window's mean
    /// depth and its per-threshold base counts.
    fn windowed(
        depths: &[u32],
        size: u64,
        thresholds: &[u32],
    ) -> (Vec<WindowDepth>, Vec<ThresholdRow>) {
        let mut windows = Vec::new();
        let mut rows = Vec::new();
        for (index, chunk) in depths.chunks(size as usize).enumerate() {
            let start = index as u64 * size;
            let end = start + chunk.len() as u64;
            let sum: u64 = chunk.iter().map(|d| u64::from(*d)).sum();
            windows.push(WindowDepth {
                start,
                end,
                mean: sum as f64 / chunk.len() as f64,
            });
            if !thresholds.is_empty() {
                let counts = thresholds
                    .iter()
                    .map(|t| chunk.iter().filter(|d| *d >= t).count() as u64)
                    .collect();
                rows.push(ThresholdRow { start, end, counts });
            }
        }
        (windows, rows)
    }

    /// Mean depth over the contig.
    pub fn mean(&self) -> f64 {
        if self.length == 0 {
            0.0
        } else {
            self.total_bases as f64 / self.length as f64
        }
    }

    /// Histogram of window mean depths, rounded to the nearest integer, which
    /// is what the region distribution is built from.
    pub fn region_histogram(&self) -> BTreeMap<u32, u64> {
        let mut hist = BTreeMap::new();
        for window in &self.windows {
            let key = window.mean.round().max(0.0) as u32;
            *hist.entry(key).or_insert(0) += 1;
        }
        hist
    }
}

/// The mosdepth result for one alignment file.
#[derive(Debug, Clone)]
pub struct MosdepthResult {
    /// Per-contig results, in alignment-header order.
    pub contigs: Vec<ContigDepth>,
    /// Window size, when per-window output was requested.
    pub window_size: Option<u32>,
    /// Requested coverage thresholds, in the order they are reported.
    pub thresholds: Vec<u32>,
}

impl MosdepthResult {
    /// Total length across all contigs.
    pub fn total_length(&self) -> u64 {
        self.contigs.iter().map(|c| c.length).sum()
    }

    /// Total covered bases across all contigs.
    pub fn total_bases(&self) -> u64 {
        self.contigs.iter().map(|c| c.total_bases).sum()
    }

    /// Mean depth across all contigs.
    pub fn mean(&self) -> f64 {
        let length = self.total_length();
        if length == 0 {
            0.0
        } else {
            self.total_bases() as f64 / length as f64
        }
    }

    /// Lowest depth across all contigs.
    pub fn min(&self) -> u32 {
        self.contigs.iter().map(|c| c.min).min().unwrap_or(0)
    }

    /// Highest depth across all contigs.
    pub fn max(&self) -> u32 {
        self.contigs.iter().map(|c| c.max).max().unwrap_or(0)
    }
}

/// Merge histograms element-wise.
pub fn merge_histograms<'a>(
    parts: impl IntoIterator<Item = &'a BTreeMap<u32, u64>>,
) -> BTreeMap<u32, u64> {
    let mut merged = BTreeMap::new();
    for part in parts {
        for (depth, count) in part {
            *merged.entry(*depth).or_insert(0) += count;
        }
    }
    merged
}

/// The depths that get a distribution row, in descending order.
///
/// The rule was derived from the committed fixtures and holds for both
/// distribution files: every depth from 0 up to `min(DIST_DENSE_MAX, max)`
/// gets a row whether or not anything sits at it, and above
/// [`DIST_DENSE_MAX`] only depths that actually occur and lie strictly below
/// the maximum do.
///
/// The consequence worth stating plainly: the maximum observed depth gets a
/// row when it falls inside the dense range and no row when it does not. On
/// the project fixture the global distribution tops out at 866 with a maximum
/// of 867, while the region distribution does emit its maximum of 204.
pub fn dist_rows(histogram: &BTreeMap<u32, u64>) -> Vec<u32> {
    let observed_max = histogram
        .iter()
        .filter(|(_, count)| **count > 0)
        .map(|(depth, _)| *depth)
        .max()
        .unwrap_or(0);

    let mut depths: Vec<u32> = histogram
        .iter()
        .filter(|(depth, count)| **count > 0 && **depth > DIST_DENSE_MAX && **depth < observed_max)
        .map(|(depth, _)| *depth)
        .collect();
    depths.extend(0..=DIST_DENSE_MAX.min(observed_max));
    depths.sort_unstable_by(|a, b| b.cmp(a));
    depths.dedup();
    depths
}

/// Cumulative proportion at or above each depth in `rows`, given `histogram`
/// and a total to divide by.
pub fn dist_proportions(histogram: &BTreeMap<u32, u64>, rows: &[u32], total: u64) -> Vec<f64> {
    if total == 0 {
        return vec![0.0; rows.len()];
    }
    rows.iter()
        .map(|threshold| {
            let at_or_above: u64 = histogram
                .iter()
                .filter(|(depth, _)| *depth >= threshold)
                .map(|(_, count)| count)
                .sum();
            at_or_above as f64 / total as f64
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn runs_collapse_equal_neighbours() {
        let c = ContigDepth::from_depths("chr1", &[0, 0, 0, 2, 2, 1], None, &[]);
        assert_eq!(
            c.runs,
            vec![
                DepthRun {
                    start: 0,
                    end: 3,
                    depth: 0
                },
                DepthRun {
                    start: 3,
                    end: 5,
                    depth: 2
                },
                DepthRun {
                    start: 5,
                    end: 6,
                    depth: 1
                },
            ]
        );
    }

    #[test]
    fn summary_figures_are_computed_over_the_whole_contig() {
        let c = ContigDepth::from_depths("chr1", &[0, 0, 3, 5], None, &[]);
        assert_eq!(c.length, 4);
        assert_eq!(c.total_bases, 8);
        assert_eq!(c.min, 0);
        assert_eq!(c.max, 5);
        assert!((c.mean() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn windows_cover_the_tail_even_when_shorter_than_the_window() {
        let c = ContigDepth::from_depths("chr1", &[4, 4, 4, 4, 10], Some(4), &[]);
        assert_eq!(c.windows.len(), 2);
        assert_eq!(
            c.windows[0],
            WindowDepth {
                start: 0,
                end: 4,
                mean: 4.0
            }
        );
        assert_eq!(
            c.windows[1],
            WindowDepth {
                start: 4,
                end: 5,
                mean: 10.0
            }
        );
    }

    #[test]
    fn threshold_counts_are_at_or_above_each_threshold() {
        let c = ContigDepth::from_depths("chr1", &[0, 1, 5, 10], Some(4), &[1, 5, 20]);
        assert_eq!(c.thresholds.len(), 1);
        assert_eq!(c.thresholds[0].counts, vec![3, 2, 0]);
    }

    #[test]
    fn dist_rows_emit_a_maximum_that_falls_inside_the_dense_range() {
        let mut hist = BTreeMap::new();
        hist.insert(0u32, 10u64);
        hist.insert(204, 1); // the maximum, but below DIST_DENSE_MAX
        let rows = dist_rows(&hist);
        assert_eq!(
            rows.first(),
            Some(&204),
            "a maximum inside the dense range is emitted"
        );
        assert_eq!(rows.len(), 205, "0 through 204 inclusive");
    }

    #[test]
    fn dist_rows_skip_a_maximum_above_the_dense_range() {
        let mut hist = BTreeMap::new();
        hist.insert(0u32, 10u64);
        hist.insert(5, 2);
        hist.insert(400, 1);
        hist.insert(500, 1); // the maximum, never emitted
        let rows = dist_rows(&hist);
        assert!(!rows.contains(&500), "the maximum depth gets no row");
        assert!(
            rows.contains(&400),
            "an observed depth above the dense range does"
        );
        assert!(
            rows.contains(&7),
            "an unobserved depth inside the dense range does"
        );
        assert!(
            !rows.contains(&350),
            "an unobserved depth above the dense range does not"
        );
        assert_eq!(rows.first(), Some(&400), "descending order");
        assert_eq!(rows.last(), Some(&0), "down to zero");
    }

    #[test]
    fn dist_proportions_are_cumulative_from_the_top() {
        let mut hist = BTreeMap::new();
        hist.insert(0u32, 2u64);
        hist.insert(1, 1);
        hist.insert(3, 1);
        let rows = vec![3u32, 2, 1, 0];
        let props = dist_proportions(&hist, &rows, 4);
        assert!((props[0] - 0.25).abs() < 1e-12);
        assert!((props[1] - 0.25).abs() < 1e-12);
        assert!((props[2] - 0.50).abs() < 1e-12);
        assert!((props[3] - 1.00).abs() < 1e-12);
    }

    #[test]
    fn region_histogram_rounds_window_means() {
        let c = ContigDepth::from_depths("chr1", &[1, 2, 2, 3], Some(2), &[]);
        // Windows: mean 1.5 rounds to 2, mean 2.5 rounds to 3 (away from zero).
        let hist = c.region_histogram();
        assert_eq!(hist.get(&2), Some(&1));
        assert_eq!(hist.get(&3), Some(&1));
    }
}
