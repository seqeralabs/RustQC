//! Read-level metrics, reproducing what seqkit and FastQC report.
//!
//! # Upstream semantics
//!
//! Derived by reproducing seqkit 2.13.0 and FastQC 0.12.1 on the project
//! fixture. Two of FastQC's conventions are worth stating:
//!
//! **Base positions are binned.** The first nine positions are reported
//! individually, then in groups of five, with the last group truncated to the
//! longest read. A 151 base read gives `1` through `9`, then `10-14` up to
//! `150-151`.
//!
//! **Quantiles come from per-position histograms**, not from sorting. Quality
//! scores are small integers, so a histogram per position gives exact
//! quantiles at bounded memory whatever the file size.

use std::collections::BTreeMap;

/// Highest Phred score tracked. Modern instruments top out well below this.
pub const MAX_PHRED: usize = 64;

/// Per-position accumulators.
#[derive(Debug, Clone)]
struct PositionStats {
    /// Count of each quality score seen at this position.
    quality: [u64; MAX_PHRED],
    /// Bases seen at this position, by base.
    bases: [u64; 5],
}

impl Default for PositionStats {
    fn default() -> Self {
        Self {
            quality: [0; MAX_PHRED],
            bases: [0; 5],
        }
    }
}

/// The error probability a Phred score encodes.
fn phred_to_probability(score: u8) -> f64 {
    10f64.powf(-f64::from(score) / 10.0)
}

/// Index of a base in the fixed `A, C, G, T, N` ordering.
fn base_index(base: u8) -> usize {
    match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

/// Everything the `reads` subcommand reports about one FASTQ.
#[derive(Debug, Clone, Default)]
pub struct ReadMetrics {
    /// Records seen.
    pub reads: u64,
    /// Bases across every record.
    pub bases: u64,
    /// `G` and `C` bases.
    pub gc_bases: u64,
    /// `N` bases.
    pub n_bases: u64,
    /// Bases at Phred 20 or better.
    pub q20_bases: u64,
    /// Bases at Phred 30 or better.
    pub q30_bases: u64,
    /// Summed quality across every base, for the arithmetic mean.
    quality_sum: u64,
    /// Summed error probability across every base, for the mean seqkit
    /// reports.
    error_probability_sum: f64,
    /// Read length histogram.
    pub length_histogram: BTreeMap<u64, u64>,
    /// Mean read quality histogram, keyed by the truncated mean.
    pub read_quality_histogram: BTreeMap<u8, u64>,
    /// Read GC percentage histogram, keyed by the rounded percentage.
    pub gc_histogram: BTreeMap<u8, u64>,
    /// Per-position accumulators, grown as longer reads arrive.
    positions: Vec<PositionStats>,
}

impl ReadMetrics {
    /// Fold one record in.
    pub fn observe(&mut self, sequence: &[u8], quality: &[u8]) {
        self.reads += 1;
        self.bases += sequence.len() as u64;
        *self
            .length_histogram
            .entry(sequence.len() as u64)
            .or_insert(0) += 1;

        if self.positions.len() < sequence.len() {
            self.positions
                .resize(sequence.len(), PositionStats::default());
        }

        let mut gc = 0u64;
        let mut read_quality = 0u64;
        for (index, (base, score)) in sequence.iter().zip(quality).enumerate() {
            let position = &mut self.positions[index];
            position.bases[base_index(*base)] += 1;
            position.quality[(*score as usize).min(MAX_PHRED - 1)] += 1;

            match base {
                b'G' | b'C' => {
                    gc += 1;
                    self.gc_bases += 1;
                }
                b'N' => self.n_bases += 1,
                _ => {}
            }
            if *score >= 20 {
                self.q20_bases += 1;
            }
            if *score >= 30 {
                self.q30_bases += 1;
            }
            read_quality += u64::from(*score);
            self.quality_sum += u64::from(*score);
            self.error_probability_sum += phred_to_probability(*score);
        }

        if !sequence.is_empty() {
            // FastQC bins a read by its mean quality, truncated.
            let mean = (read_quality / sequence.len() as u64) as u8;
            *self.read_quality_histogram.entry(mean).or_insert(0) += 1;

            // GC percentage per read, rounded to the nearest whole percent.
            let percent = (gc as f64 / sequence.len() as f64 * 100.0).round() as u8;
            *self.gc_histogram.entry(percent.min(100)).or_insert(0) += 1;
        }
    }

    /// Shortest read.
    pub fn min_length(&self) -> u64 {
        self.length_histogram.keys().next().copied().unwrap_or(0)
    }

    /// Longest read.
    pub fn max_length(&self) -> u64 {
        self.length_histogram
            .keys()
            .next_back()
            .copied()
            .unwrap_or(0)
    }

    /// Mean read length.
    pub fn mean_length(&self) -> f64 {
        if self.reads == 0 {
            0.0
        } else {
            self.bases as f64 / self.reads as f64
        }
    }

    /// Arithmetic mean of the Phred scores.
    ///
    /// Not what seqkit reports; see [`ReadMetrics::average_quality`].
    pub fn mean_phred(&self) -> f64 {
        if self.bases == 0 {
            0.0
        } else {
            self.quality_sum as f64 / self.bases as f64
        }
    }

    /// Average quality as seqkit reports it: the mean error probability,
    /// converted back to a Phred score.
    ///
    /// This is not the mean of the Phred scores, and the difference is large.
    /// On the project fixture the arithmetic mean is 34.00 while this is
    /// 25.60, because averaging in probability space is dominated by the worst
    /// bases. Reporting the arithmetic mean would flatter the run.
    pub fn average_quality(&self) -> f64 {
        if self.bases == 0 {
            return 0.0;
        }
        let mean_probability = self.error_probability_sum / self.bases as f64;
        if mean_probability <= 0.0 {
            return 0.0;
        }
        -10.0 * mean_probability.log10()
    }

    /// Percentage of bases at Phred 20 or better.
    pub fn q20_percent(&self) -> f64 {
        self.percent(self.q20_bases)
    }

    /// Percentage of bases at Phred 30 or better.
    pub fn q30_percent(&self) -> f64 {
        self.percent(self.q30_bases)
    }

    /// Percentage of bases that are `G` or `C`.
    pub fn gc_percent(&self) -> f64 {
        self.percent(self.gc_bases)
    }

    fn percent(&self, part: u64) -> f64 {
        if self.bases == 0 {
            0.0
        } else {
            part as f64 / self.bases as f64 * 100.0
        }
    }

    /// The position bins FastQC reports, as inclusive one-based ranges.
    ///
    /// The first nine positions stand alone, then groups of five, with the
    /// last truncated to the longest read.
    pub fn position_bins(&self) -> Vec<(u64, u64)> {
        let max = self.positions.len() as u64;
        let mut bins = Vec::new();
        let mut start = 1u64;
        while start <= max {
            if start <= 9 {
                bins.push((start, start));
                start += 1;
            } else {
                let end = (start + 4).min(max);
                bins.push((start, end));
                start = end + 1;
            }
        }
        bins
    }

    /// Quality statistics for one position bin.
    ///
    /// Every figure, the mean and all five quantiles, is the unweighted mean
    /// of that figure across the bin's positions. Pooling the positions into
    /// one histogram instead would give integer quantiles and weight the
    /// earlier positions more heavily once reads start running out; FastQC
    /// averages, which is why its tenth percentile for a bin can read 35.2
    /// when every position's own tenth percentile is a whole number.
    pub fn position_quality(&self, bin: (u64, u64)) -> QualitySummary {
        let mut summaries = Vec::new();
        for position in bin.0..=bin.1 {
            let index = position as usize - 1;
            if let Some(stats) = self.positions.get(index) {
                if stats.quality.iter().sum::<u64>() == 0 {
                    continue;
                }
                summaries.push(QualitySummary::from_histogram(&stats.quality));
            }
        }
        if summaries.is_empty() {
            return QualitySummary {
                mean: 0.0,
                median: 0.0,
                lower_quartile: 0.0,
                upper_quartile: 0.0,
                tenth: 0.0,
                ninetieth: 0.0,
            };
        }
        let n = summaries.len() as f64;
        let average =
            |f: fn(&QualitySummary) -> f64| -> f64 { summaries.iter().map(f).sum::<f64>() / n };
        QualitySummary {
            mean: average(|s| s.mean),
            median: average(|s| s.median),
            lower_quartile: average(|s| s.lower_quartile),
            upper_quartile: average(|s| s.upper_quartile),
            tenth: average(|s| s.tenth),
            ninetieth: average(|s| s.ninetieth),
        }
    }

    /// Percentage of `N` bases in one position bin.
    pub fn position_n_content(&self, bin: (u64, u64)) -> f64 {
        let mut n = 0u64;
        let mut total = 0u64;
        for position in bin.0..=bin.1 {
            if let Some(stats) = self.positions.get(position as usize - 1) {
                n += stats.bases[4];
                total += stats.bases.iter().sum::<u64>();
            }
        }
        if total == 0 {
            0.0
        } else {
            n as f64 / total as f64 * 100.0
        }
    }

    /// Base composition percentages for one position bin, as `A, C, G, T`.
    ///
    /// `N` bases are excluded from the denominator, so the four percentages
    /// sum to 100 even where the instrument called nothing. Including them
    /// would shift every percentage by the N rate, which on the project
    /// fixture is enough to move the fourth decimal.
    pub fn position_composition(&self, bin: (u64, u64)) -> [f64; 4] {
        let mut counts = [0u64; 5];
        for position in bin.0..=bin.1 {
            if let Some(stats) = self.positions.get(position as usize - 1) {
                for (index, count) in stats.bases.iter().enumerate() {
                    counts[index] += count;
                }
            }
        }
        let total: u64 = counts[..4].iter().sum();
        if total == 0 {
            return [0.0; 4];
        }
        [
            counts[0] as f64 / total as f64 * 100.0,
            counts[1] as f64 / total as f64 * 100.0,
            counts[2] as f64 / total as f64 * 100.0,
            counts[3] as f64 / total as f64 * 100.0,
        ]
    }
}

/// The six figures FastQC reports per position.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QualitySummary {
    /// Arithmetic mean.
    pub mean: f64,
    /// Median.
    pub median: f64,
    /// Lower quartile.
    pub lower_quartile: f64,
    /// Upper quartile.
    pub upper_quartile: f64,
    /// Tenth percentile.
    pub tenth: f64,
    /// Ninetieth percentile.
    pub ninetieth: f64,
}

impl QualitySummary {
    /// Summarise a quality histogram.
    fn from_histogram(histogram: &[u64; MAX_PHRED]) -> Self {
        let total: u64 = histogram.iter().sum();
        if total == 0 {
            return Self {
                mean: 0.0,
                median: 0.0,
                lower_quartile: 0.0,
                upper_quartile: 0.0,
                tenth: 0.0,
                ninetieth: 0.0,
            };
        }
        let sum: u64 = histogram
            .iter()
            .enumerate()
            .map(|(score, count)| score as u64 * count)
            .sum();
        let percentile = |fraction: f64| -> f64 {
            let target = (total as f64 * fraction).ceil() as u64;
            let mut seen = 0u64;
            for (score, count) in histogram.iter().enumerate() {
                seen += count;
                if seen >= target.max(1) {
                    return score as f64;
                }
            }
            0.0
        };
        Self {
            mean: sum as f64 / total as f64,
            median: percentile(0.5),
            lower_quartile: percentile(0.25),
            upper_quartile: percentile(0.75),
            tenth: percentile(0.10),
            ninetieth: percentile(0.90),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn position_bins_are_singles_then_fives_truncated_at_the_end() {
        let mut m = ReadMetrics::default();
        m.observe(&vec![b'A'; 12], &vec![30u8; 12]);
        assert_eq!(
            m.position_bins(),
            vec![
                (1, 1),
                (2, 2),
                (3, 3),
                (4, 4),
                (5, 5),
                (6, 6),
                (7, 7),
                (8, 8),
                (9, 9),
                (10, 12),
            ],
            "the last bin stops at the longest read, not at 14"
        );
    }

    #[test]
    fn a_short_read_gives_only_single_position_bins() {
        let mut m = ReadMetrics::default();
        m.observe(b"ACG", &[30, 30, 30]);
        assert_eq!(m.position_bins(), vec![(1, 1), (2, 2), (3, 3)]);
    }

    #[test]
    fn composition_excludes_n_from_its_denominator() {
        let mut m = ReadMetrics::default();
        m.observe(b"AN", &[30, 30]);
        let [a, c, g, t] = m.position_composition((1, 1));
        assert!(
            (a - 100.0).abs() < 1e-12,
            "the one called base is all of it"
        );
        assert_eq!([c, g, t], [0.0, 0.0, 0.0]);
        // Position 2 is nothing but N, so there is no composition to report.
        assert_eq!(m.position_composition((2, 2)), [0.0; 4]);
    }

    #[test]
    fn gc_and_n_are_counted_separately_and_not_as_each_other() {
        let mut m = ReadMetrics::default();
        m.observe(b"ACGTN", &[30, 30, 30, 30, 30]);
        assert_eq!(m.gc_bases, 2, "C and G");
        assert_eq!(m.n_bases, 1);
        assert_eq!(m.bases, 5);
        assert!((m.gc_percent() - 40.0).abs() < 1e-12);
    }

    #[test]
    fn quality_thresholds_are_at_or_above() {
        let mut m = ReadMetrics::default();
        m.observe(b"AAAA", &[19, 20, 29, 30]);
        assert_eq!(m.q20_bases, 3, "20 itself counts");
        assert_eq!(m.q30_bases, 1, "30 itself counts");
    }

    #[test]
    fn reads_of_different_lengths_grow_the_position_table() {
        let mut m = ReadMetrics::default();
        m.observe(b"AC", &[30, 30]);
        m.observe(b"ACGTA", &[30, 30, 30, 30, 30]);
        assert_eq!(m.max_length(), 5);
        assert_eq!(m.min_length(), 2);
        assert_eq!(m.position_bins().len(), 5);
        // Position 5 only saw one read.
        assert!((m.position_n_content((5, 5)) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn a_bin_averages_its_positions_quantiles_rather_than_pooling_them() {
        let mut m = ReadMetrics::default();
        // Positions 10 and 11 have medians of 10 and 40; the bin's median is
        // their average, 25, which no single base carries.
        m.observe(&[b'A'; 11], &[0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 40]);
        let q = m.position_quality((10, 14));
        assert!((q.median - 25.0).abs() < 1e-12, "got {}", q.median);
        assert!((q.mean - 25.0).abs() < 1e-12);
    }

    #[test]
    fn per_position_quantiles_come_out_of_the_histogram() {
        let mut m = ReadMetrics::default();
        // Four reads, one base each, qualities 10, 20, 30, 40.
        for score in [10u8, 20, 30, 40] {
            m.observe(b"A", &[score]);
        }
        let q = m.position_quality((1, 1));
        assert!((q.mean - 25.0).abs() < 1e-12);
        assert_eq!(q.median, 20.0, "the value at the halfway rank");
        assert_eq!(q.lower_quartile, 10.0);
        assert_eq!(q.upper_quartile, 30.0);
    }

    #[test]
    fn average_quality_is_computed_in_probability_space() {
        let mut m = ReadMetrics::default();
        // One base at Phred 40 and one at Phred 10: the arithmetic mean is 25,
        // but the mean error probability is (0.0001 + 0.1) / 2, which is Phred
        // 13.0 to one decimal.
        m.observe(b"AA", &[40, 10]);
        assert!((m.mean_phred() - 25.0).abs() < 1e-12);
        let average = m.average_quality();
        assert!(
            (average - 13.0).abs() < 0.05,
            "expected about 13, got {average}"
        );
        assert!(
            average < m.mean_phred(),
            "probability-space averaging is dominated by the worse base"
        );
    }

    #[test]
    fn an_empty_file_reports_zeroes_rather_than_dividing_by_zero() {
        let m = ReadMetrics::default();
        assert_eq!(m.reads, 0);
        assert_eq!(m.mean_length(), 0.0);
        assert_eq!(m.mean_phred(), 0.0);
        assert_eq!(m.average_quality(), 0.0);
        assert_eq!(m.gc_percent(), 0.0);
        assert!(m.position_bins().is_empty());
    }
}
