//! Picard `CollectInsertSizeMetrics` reimplementation.
//!
//! # Upstream semantics
//!
//! Every rule below was measured against Picard 3.4.0 output on
//! `tests/data/dna/test.dna.bam`, not recalled from documentation.
//!
//! A record contributes when it is paired, is neither secondary,
//! supplementary, duplicate-flagged nor unmapped, has a mapped mate, and
//! carries a positive `TLEN`. Taking only the positive `TLEN` of the two is
//! what counts each pair once. Proper-pair is deliberately **not** required:
//! requiring it drops one pair and shortens the maximum from 300 to 239 on the
//! project fixture.
//!
//! Pairs are grouped by orientation (`FR`, `RF`, `TANDEM`), each group
//! reported on its own row with its own histogram, exactly as Picard does.
//!
//! `MEAN_INSERT_SIZE` and `STANDARD_DEVIATION` are computed over the
//! histogram trimmed to `DEVIATIONS` median absolute deviations either side of
//! the median, and the standard deviation uses the `n - 1` denominator.
//! `MIN_INSERT_SIZE` and `MAX_INSERT_SIZE` are over the untrimmed set.
//!
//! `WIDTH_OF_XX_PERCENT` is the width of the smallest window centred on the
//! median that covers at least `XX` percent of pairs: grow `i` from zero until
//! the bins from `median - i` to `median + i` cover the target, then report
//! `2i + 1`.

use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use rust_htslib::bam;

use crate::common::bam_flags::*;

/// Percentiles Picard reports a width for, in output order.
pub const WIDTH_PERCENTILES: [u32; 11] = [10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99];

/// Picard's `DEVIATIONS` default: how many median absolute deviations either
/// side of the median survive trimming before the mean and standard deviation
/// are computed.
pub const DEFAULT_DEVIATIONS: f64 = 10.0;

/// Relative orientation of the two mates of a pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum PairOrientation {
    /// Forward-reverse, the usual Illumina paired-end arrangement.
    Fr,
    /// Reverse-forward, seen in mate-pair and some capture libraries.
    Rf,
    /// Both mates on the same strand.
    Tandem,
}

impl PairOrientation {
    /// The label Picard writes in the `PAIR_ORIENTATION` column.
    pub fn label(&self) -> &'static str {
        match self {
            PairOrientation::Fr => "FR",
            PairOrientation::Rf => "RF",
            PairOrientation::Tandem => "TANDEM",
        }
    }

    /// The prefix Picard uses for this orientation's histogram column.
    fn histogram_column(&self) -> &'static str {
        match self {
            PairOrientation::Fr => "fr",
            PairOrientation::Rf => "rf",
            PairOrientation::Tandem => "tandem",
        }
    }
}

/// Accumulates insert sizes, one histogram per orientation.
#[derive(Debug, Default)]
pub struct InsertSizeAccum {
    histograms: BTreeMap<PairOrientation, BTreeMap<u64, u64>>,
}

impl InsertSizeAccum {
    /// A new, empty accumulator.
    pub fn new() -> Self {
        Self::default()
    }

    /// Offer one record. Records that do not represent a countable pair are
    /// ignored.
    pub fn process_read(&mut self, record: &bam::Record) {
        let flags = record.flags();
        if flags & BAM_FPAIRED == 0 {
            return;
        }
        let excluded = BAM_FUNMAP | BAM_FMUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FDUP;
        if flags & excluded != 0 {
            return;
        }
        // Only the mate carrying the positive TLEN counts, so each pair is
        // counted once.
        let insert_size = record.insert_size();
        if insert_size <= 0 {
            return;
        }

        let orientation = orientation_of(flags);
        *self
            .histograms
            .entry(orientation)
            .or_default()
            .entry(insert_size as u64)
            .or_insert(0) += 1;
    }

    /// Fold another accumulator into this one.
    pub fn merge(&mut self, other: InsertSizeAccum) {
        for (orientation, histogram) in other.histograms {
            let target = self.histograms.entry(orientation).or_default();
            for (size, count) in histogram {
                *target.entry(size).or_insert(0) += count;
            }
        }
    }

    /// Summarise each orientation, in Picard's output order: most pairs first.
    pub fn into_result(self, deviations: f64) -> InsertSizeResult {
        let mut rows: Vec<InsertSizeRow> = self
            .histograms
            .into_iter()
            .map(|(orientation, histogram)| InsertSizeRow::new(orientation, histogram, deviations))
            .collect();
        rows.sort_by_key(|row| std::cmp::Reverse(row.read_pairs));
        InsertSizeResult { rows }
    }
}

/// Which orientation a record's flags describe.
fn orientation_of(flags: u16) -> PairOrientation {
    let read_reverse = flags & BAM_FREVERSE != 0;
    let mate_reverse = flags & BAM_FMREVERSE != 0;
    if read_reverse == mate_reverse {
        PairOrientation::Tandem
    } else if read_reverse {
        // This record carries the positive TLEN, so it is the leftmost mate.
        // Leftmost on the reverse strand means reverse-forward.
        PairOrientation::Rf
    } else {
        PairOrientation::Fr
    }
}

/// One orientation's metrics and histogram.
#[derive(Debug, Clone)]
pub struct InsertSizeRow {
    /// The orientation this row describes.
    pub orientation: PairOrientation,
    /// Insert size histogram, size to pair count.
    pub histogram: BTreeMap<u64, u64>,
    /// Number of pairs counted.
    pub read_pairs: u64,
    /// Median insert size.
    pub median: u64,
    /// Most frequent insert size; ties go to the smaller size.
    pub mode: u64,
    /// Median absolute deviation from the median.
    pub median_absolute_deviation: u64,
    /// Smallest insert size seen, before trimming.
    pub min: u64,
    /// Largest insert size seen, before trimming.
    pub max: u64,
    /// Mean over the trimmed histogram.
    pub mean: f64,
    /// Standard deviation over the trimmed histogram, `n - 1` denominator.
    pub standard_deviation: f64,
    /// Width of the smallest median-centred window covering each percentile,
    /// in the order of [`WIDTH_PERCENTILES`].
    pub widths: Vec<u64>,
}

impl InsertSizeRow {
    fn new(orientation: PairOrientation, histogram: BTreeMap<u64, u64>, deviations: f64) -> Self {
        let read_pairs: u64 = histogram.values().sum();
        let median = quantile(&histogram, read_pairs / 2);
        let mode = histogram
            .iter()
            .max_by_key(|(size, count)| (**count, std::cmp::Reverse(**size)))
            .map(|(size, _)| *size)
            .unwrap_or(0);

        // Median absolute deviation, itself a median over |size - median|.
        let mut deviation_histogram: BTreeMap<u64, u64> = BTreeMap::new();
        for (size, count) in &histogram {
            let deviation = size.abs_diff(median);
            *deviation_histogram.entry(deviation).or_insert(0) += count;
        }
        let median_absolute_deviation = quantile(&deviation_histogram, read_pairs / 2);

        let min = histogram.keys().copied().min().unwrap_or(0);
        let max = histogram.keys().copied().max().unwrap_or(0);

        // Trim to `deviations` MADs either side before the mean and SD.
        let span = deviations * median_absolute_deviation as f64;
        let low = (median as f64 - span).max(0.0);
        let high = median as f64 + span;
        let trimmed: Vec<(u64, u64)> = histogram
            .iter()
            .filter(|(size, _)| **size as f64 >= low && **size as f64 <= high)
            .map(|(size, count)| (*size, *count))
            .collect();

        let n: u64 = trimmed.iter().map(|(_, count)| count).sum();
        let mean = if n == 0 {
            0.0
        } else {
            trimmed
                .iter()
                .map(|(size, count)| *size as f64 * *count as f64)
                .sum::<f64>()
                / n as f64
        };
        let standard_deviation = if n < 2 {
            0.0
        } else {
            let variance = trimmed
                .iter()
                .map(|(size, count)| {
                    let diff = *size as f64 - mean;
                    diff * diff * *count as f64
                })
                .sum::<f64>()
                / (n - 1) as f64;
            variance.sqrt()
        };

        let widths = WIDTH_PERCENTILES
            .iter()
            .map(|pct| width_of_percent(&histogram, median, read_pairs, *pct))
            .collect();

        Self {
            orientation,
            histogram,
            read_pairs,
            median,
            mode,
            median_absolute_deviation,
            min,
            max,
            mean,
            standard_deviation,
            widths,
        }
    }
}

/// The value at `rank` when the histogram is expanded into a sorted list.
fn quantile(histogram: &BTreeMap<u64, u64>, rank: u64) -> u64 {
    let mut seen = 0u64;
    for (value, count) in histogram {
        seen += count;
        if seen > rank {
            return *value;
        }
    }
    histogram.keys().next_back().copied().unwrap_or(0)
}

/// Width of the smallest window centred on `median` covering `pct` percent of
/// `total` pairs.
fn width_of_percent(histogram: &BTreeMap<u64, u64>, median: u64, total: u64, pct: u32) -> u64 {
    if total == 0 {
        return 0;
    }
    let target = total as f64 * pct as f64 / 100.0;
    let mut covered = *histogram.get(&median).unwrap_or(&0) as f64;
    let mut i = 0u64;
    while covered < target {
        i += 1;
        covered += *histogram.get(&(median.saturating_sub(i))).unwrap_or(&0) as f64;
        covered += *histogram.get(&(median + i)).unwrap_or(&0) as f64;
        // Once the window spans the whole histogram there is nothing left to add.
        if median + i > *histogram.keys().next_back().unwrap_or(&0) && median < i {
            break;
        }
    }
    2 * i + 1
}

/// All orientations' metrics for one alignment file.
#[derive(Debug, Clone)]
pub struct InsertSizeResult {
    /// One row per orientation seen, most pairs first.
    pub rows: Vec<InsertSizeRow>,
}

/// Format a float the way Picard's metrics writer does: up to six decimals,
/// trailing zeros removed, and a bare integer when there is no fraction.
fn fmt_picard(value: f64) -> String {
    if value == value.trunc() && value.abs() < 1e15 {
        return format!("{}", value as i64);
    }
    let text = format!("{value:.6}");
    let trimmed = text.trim_end_matches('0').trim_end_matches('.');
    trimmed.to_string()
}

/// Write a Picard-compatible `insert_size_metrics.txt`.
///
/// The `## htsjdk...StringHeader` preamble Picard writes is omitted: it holds
/// only the command line and a start timestamp, both of which are noise in a
/// reproducible pipeline.
pub fn write_insert_size_metrics(result: &InsertSizeResult, path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .map(std::io::BufWriter::new)
        .with_context(|| format!("Failed to create insert size metrics: {}", path.display()))?;

    writeln!(out, "## METRICS CLASS\tpicard.analysis.InsertSizeMetrics")?;
    write!(
        out,
        "MEDIAN_INSERT_SIZE\tMODE_INSERT_SIZE\tMEDIAN_ABSOLUTE_DEVIATION\tMIN_INSERT_SIZE\t\
         MAX_INSERT_SIZE\tMEAN_INSERT_SIZE\tSTANDARD_DEVIATION\tREAD_PAIRS\tPAIR_ORIENTATION"
    )?;
    for pct in WIDTH_PERCENTILES {
        write!(out, "\tWIDTH_OF_{pct}_PERCENT")?;
    }
    writeln!(out, "\tSAMPLE\tLIBRARY\tREAD_GROUP")?;

    for row in &result.rows {
        write!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.median,
            row.mode,
            row.median_absolute_deviation,
            row.min,
            row.max,
            fmt_picard(row.mean),
            fmt_picard(row.standard_deviation),
            row.read_pairs,
            row.orientation.label(),
        )?;
        for width in &row.widths {
            write!(out, "\t{width}")?;
        }
        // Trailing SAMPLE, LIBRARY and READ_GROUP columns are empty at the
        // ALL_READS accumulation level, which is Picard's default.
        writeln!(out, "\t\t\t")?;
    }

    writeln!(out)?;
    writeln!(out, "## HISTOGRAM\tjava.lang.Integer")?;
    write!(out, "insert_size")?;
    for row in &result.rows {
        write!(
            out,
            "\tAll_Reads.{}_count",
            row.orientation.histogram_column()
        )?;
    }
    writeln!(out)?;

    // One row per insert size seen in any orientation, ascending.
    let mut sizes: Vec<u64> = result
        .rows
        .iter()
        .flat_map(|row| row.histogram.keys().copied())
        .collect();
    sizes.sort_unstable();
    sizes.dedup();
    for size in sizes {
        write!(out, "{size}")?;
        for row in &result.rows {
            write!(out, "\t{}", row.histogram.get(&size).copied().unwrap_or(0))?;
        }
        writeln!(out)?;
    }
    writeln!(out)?;

    out.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hist(pairs: &[(u64, u64)]) -> BTreeMap<u64, u64> {
        pairs.iter().copied().collect()
    }

    #[test]
    fn orientation_follows_the_strand_flags() {
        // The record carrying the positive TLEN is the leftmost mate, so its
        // own strand decides between FR and RF.
        assert_eq!(
            orientation_of(BAM_FPAIRED | BAM_FMREVERSE),
            PairOrientation::Fr,
            "leftmost forward, mate reverse"
        );
        assert_eq!(
            orientation_of(BAM_FPAIRED | BAM_FREVERSE),
            PairOrientation::Rf,
            "leftmost reverse, mate forward"
        );
        // Same strand either way round is tandem, including neither reversed.
        assert_eq!(
            orientation_of(BAM_FPAIRED),
            PairOrientation::Tandem,
            "both forward"
        );
        assert_eq!(
            orientation_of(BAM_FPAIRED | BAM_FREVERSE | BAM_FMREVERSE),
            PairOrientation::Tandem,
            "both reverse"
        );
    }

    #[test]
    fn mode_breaks_ties_towards_the_smaller_size() {
        let row = InsertSizeRow::new(
            PairOrientation::Fr,
            hist(&[(100, 5), (200, 5)]),
            DEFAULT_DEVIATIONS,
        );
        assert_eq!(row.mode, 100);
    }

    #[test]
    fn standard_deviation_uses_the_sample_denominator() {
        // Values 1, 2, 3: mean 2, sample variance 1, so SD is exactly 1.
        let row = InsertSizeRow::new(
            PairOrientation::Fr,
            hist(&[(1, 1), (2, 1), (3, 1)]),
            DEFAULT_DEVIATIONS,
        );
        assert!((row.mean - 2.0).abs() < 1e-12);
        assert!(
            (row.standard_deviation - 1.0).abs() < 1e-12,
            "got {}",
            row.standard_deviation
        );
    }

    #[test]
    fn trimming_excludes_outliers_beyond_the_deviation_span() {
        // Median 10, MAD 0, so a span of zero keeps only the median bin.
        let row = InsertSizeRow::new(
            PairOrientation::Fr,
            hist(&[(10, 9), (1000, 1)]),
            DEFAULT_DEVIATIONS,
        );
        assert_eq!(row.max, 1000, "the untrimmed maximum is still reported");
        assert!(
            (row.mean - 10.0).abs() < 1e-12,
            "the outlier must not reach the mean, got {}",
            row.mean
        );
    }

    #[test]
    fn width_grows_symmetrically_around_the_median() {
        // Ten pairs at the median, five either side one apart.
        let h = hist(&[(9, 5), (10, 10), (11, 5)]);
        assert_eq!(width_of_percent(&h, 10, 20, 50), 1, "the median bin alone");
        assert_eq!(width_of_percent(&h, 10, 20, 90), 3, "one bin either side");
    }

    #[test]
    fn picard_float_formatting_drops_trailing_zeros() {
        assert_eq!(fmt_picard(124.442269), "124.442269");
        assert_eq!(fmt_picard(3.5), "3.5");
        assert_eq!(fmt_picard(40001.0), "40001");
        assert_eq!(fmt_picard(0.0), "0");
    }

    #[test]
    fn rows_are_ordered_by_pair_count() {
        let mut accum = InsertSizeAccum::new();
        accum
            .histograms
            .insert(PairOrientation::Rf, hist(&[(100, 1)]));
        accum
            .histograms
            .insert(PairOrientation::Fr, hist(&[(100, 50)]));
        let result = accum.into_result(DEFAULT_DEVIATIONS);
        assert_eq!(result.rows[0].orientation, PairOrientation::Fr);
        assert_eq!(result.rows[1].orientation, PairOrientation::Rf);
    }
}
