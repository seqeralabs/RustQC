//! BED interval parsing and merging for targeted mode.
//!
//! Picard consumes `.interval_list` files, RustQC accepts BED. The two differ
//! in a way that is easy to get wrong: BED is zero-based half-open, an
//! interval list is one-based inclusive, so `chr22 1 15000` in BED is
//! `chr22 2 15000` in an interval list. Everything here works in BED's
//! convention internally and converts only at the edges.

use std::collections::HashMap;
use std::path::Path;

use anyhow::{bail, Context, Result};

/// A half-open interval `[start, end)` on one contig, zero-based.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Interval {
    /// Zero-based inclusive start.
    pub start: u64,
    /// Zero-based exclusive end.
    pub end: u64,
}

impl Interval {
    /// Number of bases covered.
    pub fn len(&self) -> u64 {
        self.end.saturating_sub(self.start)
    }

    /// Whether the interval covers no bases.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Whether `position` falls inside.
    pub fn contains(&self, position: u64) -> bool {
        position >= self.start && position < self.end
    }
}

/// Merged, sorted intervals grouped by contig.
#[derive(Debug, Clone, Default)]
pub struct IntervalSet {
    /// Non-overlapping intervals per contig, ascending.
    by_contig: HashMap<String, Vec<Interval>>,
    /// A name for the set, used as `BAIT_SET` in the metrics.
    name: String,
}

impl IntervalSet {
    /// Read a BED file, merging any overlapping or touching intervals.
    ///
    /// Merging matters: overlapping targets would otherwise inflate the
    /// territory and double-count on-target bases.
    pub fn from_bed(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("Failed to read BED file: {}", path.display()))?;

        let name = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("targets")
            .to_string();

        let mut raw: HashMap<String, Vec<Interval>> = HashMap::new();
        for (number, line) in text.lines().enumerate() {
            let line = line.trim();
            if line.is_empty()
                || line.starts_with('#')
                || line.starts_with("track")
                || line.starts_with("browser")
            {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 3 {
                bail!(
                    "{}: line {} has {} fields, a BED interval needs at least 3",
                    path.display(),
                    number + 1,
                    fields.len()
                );
            }
            let start: u64 = fields[1].parse().with_context(|| {
                format!("{}: line {} has a bad start", path.display(), number + 1)
            })?;
            let end: u64 = fields[2].parse().with_context(|| {
                format!("{}: line {} has a bad end", path.display(), number + 1)
            })?;
            if end <= start {
                bail!(
                    "{}: line {} ends at or before it starts",
                    path.display(),
                    number + 1
                );
            }
            raw.entry(fields[0].to_string())
                .or_default()
                .push(Interval { start, end });
        }

        let by_contig = raw
            .into_iter()
            .map(|(contig, intervals)| (contig, merge(intervals)))
            .collect();

        Ok(Self { by_contig, name })
    }

    /// Build directly from intervals, for tests and for deriving one set from
    /// another.
    pub fn from_intervals(name: &str, by_contig: HashMap<String, Vec<Interval>>) -> Self {
        Self {
            by_contig: by_contig
                .into_iter()
                .map(|(contig, intervals)| (contig, merge(intervals)))
                .collect(),
            name: name.to_string(),
        }
    }

    /// The set's name, reported as `BAIT_SET`.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Total bases covered across every contig.
    pub fn territory(&self) -> u64 {
        self.by_contig
            .values()
            .flat_map(|intervals| intervals.iter())
            .map(|interval| interval.len())
            .sum()
    }

    /// Intervals on one contig, ascending, or an empty slice.
    pub fn on(&self, contig: &str) -> &[Interval] {
        self.by_contig
            .get(contig)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Total number of intervals.
    pub fn len(&self) -> usize {
        self.by_contig.values().map(|v| v.len()).sum()
    }

    /// Whether the set holds no intervals.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// A per-base membership mask for one contig, for fast position lookup in
    /// the inner loop.
    pub fn mask(&self, contig: &str, length: u64) -> Vec<bool> {
        let mut mask = vec![false; length as usize];
        for interval in self.on(contig) {
            let start = interval.start.min(length) as usize;
            let end = interval.end.min(length) as usize;
            mask[start..end].fill(true);
        }
        mask
    }
}

/// Sort and merge overlapping or adjacent intervals.
fn merge(mut intervals: Vec<Interval>) -> Vec<Interval> {
    intervals.sort();
    let mut merged: Vec<Interval> = Vec::with_capacity(intervals.len());
    for interval in intervals {
        match merged.last_mut() {
            Some(last) if interval.start <= last.end => {
                last.end = last.end.max(interval.end);
            }
            _ => merged.push(interval),
        }
    }
    merged
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write_bed(name: &str, contents: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join("rustqc-interval-tests");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(name);
        std::fs::write(&path, contents).unwrap();
        path
    }

    #[test]
    fn overlapping_intervals_are_merged() {
        let path = write_bed(
            "overlap.bed",
            "chr1\t100\t200\nchr1\t150\t300\nchr1\t400\t500\n",
        );
        let set = IntervalSet::from_bed(&path).unwrap();
        assert_eq!(
            set.on("chr1"),
            &[
                Interval {
                    start: 100,
                    end: 300
                },
                Interval {
                    start: 400,
                    end: 500
                },
            ]
        );
        assert_eq!(set.territory(), 300, "merged, not 100 + 150 + 100");
    }

    #[test]
    fn touching_intervals_are_merged_too() {
        let path = write_bed("touch.bed", "chr1\t100\t200\nchr1\t200\t300\n");
        let set = IntervalSet::from_bed(&path).unwrap();
        assert_eq!(
            set.on("chr1"),
            &[Interval {
                start: 100,
                end: 300
            }]
        );
    }

    #[test]
    fn comments_and_track_lines_are_ignored() {
        let path = write_bed(
            "comments.bed",
            "# a comment\ntrack name=x\nchr1\t10\t20\n\nbrowser position chr1\n",
        );
        let set = IntervalSet::from_bed(&path).unwrap();
        assert_eq!(set.len(), 1);
    }

    #[test]
    fn a_backwards_interval_is_an_error_rather_than_silently_empty() {
        let path = write_bed("backwards.bed", "chr1\t200\t100\n");
        let result = IntervalSet::from_bed(&path);
        assert!(result.is_err(), "an end before the start must be rejected");
    }

    #[test]
    fn the_mask_marks_exactly_the_covered_bases() {
        let path = write_bed("mask.bed", "chr1\t2\t5\n");
        let set = IntervalSet::from_bed(&path).unwrap();
        assert_eq!(
            set.mask("chr1", 8),
            vec![false, false, true, true, true, false, false, false]
        );
    }

    #[test]
    fn intervals_beyond_the_contig_end_do_not_overflow_the_mask() {
        let path = write_bed("beyond.bed", "chr1\t2\t100\n");
        let set = IntervalSet::from_bed(&path).unwrap();
        assert_eq!(set.mask("chr1", 4), vec![false, false, true, true]);
    }
}
