//! Length and composition statistics, reproducing `seqkit stats -a -T`.
//!
//! # Upstream semantics
//!
//! Derived by reproducing seqkit 2.13.0's output on the project fixtures.
//! Two of its conventions are worth stating, because neither is what a
//! statistics library gives you by default:
//!
//! **Quartiles are Tukey's halves.** `Q1` is the median of the lower half of
//! the sorted lengths, `Q3` the median of the upper half, with the overall
//! median excluded when the count is odd. On a ten-sequence file that gives
//! `Q1 = 157` where linear interpolation gives 165 and inverse-ECDF gives 157
//! but then disagrees on `Q2`.
//!
//! **Halves are rounded half-to-even.** A median of 235.5 is reported as 236,
//! but 376.5 and 516.5 are reported as 376 and 516. Ordinary rounding gets the
//! first right and the other two wrong.
//!
//! `N50` is the largest length `L` such that sequences of at least `L` cover
//! half the total residues; `N50_num` is how many sequences that takes.

use std::collections::BTreeMap;

use super::Record;

/// The statistics `seqkit stats -a` reports, plus composition.
#[derive(Debug, Clone)]
pub struct SequenceStats {
    /// Number of sequences.
    pub count: u64,
    /// Total residues.
    pub total: u64,
    /// Shortest sequence.
    pub min: u64,
    /// Longest sequence.
    pub max: u64,
    /// Mean length.
    pub mean: f64,
    /// Lower quartile, Tukey's lower half median.
    pub q1: u64,
    /// Median.
    pub q2: u64,
    /// Upper quartile, Tukey's upper half median.
    pub q3: u64,
    /// Gap characters, `-` and `.`, across all sequences.
    pub gaps: u64,
    /// N50 length.
    pub n50: u64,
    /// Number of distinct lengths needed to reach N50, which is what seqkit
    /// reports and is not the number of sequences.
    pub n50_num: u64,
    /// Residue counts, keyed by the residue character.
    pub composition: BTreeMap<u8, u64>,
}

impl SequenceStats {
    /// Summarise a set of records.
    pub fn from_records(records: &[Record]) -> Self {
        let mut lengths: Vec<u64> = records.iter().map(|r| r.len() as u64).collect();
        lengths.sort_unstable();

        let count = lengths.len() as u64;
        let total: u64 = lengths.iter().sum();
        let min = lengths.first().copied().unwrap_or(0);
        let max = lengths.last().copied().unwrap_or(0);
        let mean = if count == 0 {
            0.0
        } else {
            total as f64 / count as f64
        };

        let (q1, q2, q3) = tukey_quartiles(&lengths);
        let (n50, n50_num) = n50(&lengths, total);

        let mut composition: BTreeMap<u8, u64> = BTreeMap::new();
        let mut gaps = 0u64;
        for record in records {
            for residue in &record.residues {
                if *residue == b'-' || *residue == b'.' {
                    gaps += 1;
                }
                *composition.entry(*residue).or_insert(0) += 1;
            }
        }

        Self {
            count,
            total,
            min,
            max,
            mean,
            q1,
            q2,
            q3,
            gaps,
            n50,
            n50_num,
            composition,
        }
    }

    /// Fraction of residues that are the given one.
    pub fn fraction(&self, residue: u8) -> f64 {
        if self.total == 0 {
            0.0
        } else {
            self.composition.get(&residue).copied().unwrap_or(0) as f64 / self.total as f64
        }
    }
}

/// Median of an ascending slice, rounded half to even when the count is even.
fn median(sorted: &[u64]) -> u64 {
    match sorted.len() {
        0 => 0,
        n if n % 2 == 1 => sorted[n / 2],
        n => round_half_to_even((sorted[n / 2 - 1] + sorted[n / 2]) as f64 / 2.0),
    }
}

/// Round to the nearest integer, ties going to the even one.
///
/// seqkit reports 235.5 as 236 and 376.5 as 376, which only half-to-even
/// explains.
fn round_half_to_even(value: f64) -> u64 {
    let floor = value.floor();
    let fraction = value - floor;
    let rounded = if (fraction - 0.5).abs() < f64::EPSILON {
        if (floor as i64) % 2 == 0 {
            floor
        } else {
            floor + 1.0
        }
    } else {
        value.round()
    };
    rounded.max(0.0) as u64
}

/// Tukey's quartiles: the medians of the two halves, the overall median
/// excluded from both when the count is odd.
fn tukey_quartiles(sorted: &[u64]) -> (u64, u64, u64) {
    let n = sorted.len();
    if n == 0 {
        return (0, 0, 0);
    }
    let half = n / 2;
    let lower = &sorted[..half];
    let upper = &sorted[n - half..];
    (median(lower), median(sorted), median(upper))
}

/// N50 and the number of *distinct lengths* needed to reach it.
///
/// The second figure is not the number of sequences, which is the obvious
/// reading and the wrong one. seqkit walks the distinct lengths from longest
/// down and counts how many it consumed: three sequences of lengths 10, 10 and
/// 3 give `N50_num` of 1, not 2, because the two tens are one length. A file
/// whose lengths are all distinct hides the difference entirely, which is why
/// the project's protein fixtures did not catch it.
fn n50(sorted_ascending: &[u64], total: u64) -> (u64, u64) {
    let mut cumulative = 0u64;
    let mut distinct = 0u64;
    let mut previous: Option<u64> = None;
    for length in sorted_ascending.iter().rev() {
        if previous != Some(*length) {
            distinct += 1;
            previous = Some(*length);
        }
        cumulative += length;
        if cumulative * 2 >= total {
            return (*length, distinct);
        }
    }
    (0, 0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn records(lengths: &[usize]) -> Vec<Record> {
        lengths
            .iter()
            .enumerate()
            .map(|(i, n)| Record {
                id: format!("s{i}"),
                residues: vec![b'A'; *n],
            })
            .collect()
    }

    #[test]
    fn quartiles_are_tukeys_halves_not_interpolated() {
        // The project's yeast fixture: linear interpolation gives 165 for Q1.
        let stats = SequenceStats::from_records(&records(&[
            81, 140, 157, 189, 198, 273, 381, 526, 584, 755,
        ]));
        assert_eq!(stats.q1, 157);
        assert_eq!(stats.q2, 236);
        assert_eq!(stats.q3, 526);
    }

    #[test]
    fn halves_round_to_even() {
        assert_eq!(round_half_to_even(235.5), 236, "236 is even");
        assert_eq!(round_half_to_even(376.5), 376, "376 is even");
        assert_eq!(round_half_to_even(516.5), 516, "516 is even");
        assert_eq!(round_half_to_even(2.4), 2, "not a tie, ordinary rounding");
        assert_eq!(round_half_to_even(2.6), 3);
    }

    #[test]
    fn an_odd_count_excludes_the_median_from_both_halves() {
        // Five values: halves are [1, 2] and [4, 5], so Q1 is 2 and Q3 is 4
        // after half-to-even rounding of 1.5 and 4.5.
        let stats = SequenceStats::from_records(&records(&[1, 2, 3, 4, 5]));
        assert_eq!(stats.q2, 3, "the median itself");
        assert_eq!(stats.q1, 2, "median of [1, 2] is 1.5, rounded to 2");
        assert_eq!(stats.q3, 4, "median of [4, 5] is 4.5, rounded to 4");
    }

    #[test]
    fn n50_num_counts_distinct_lengths_not_sequences() {
        // Two sequences of 10 and one of 3: two sequences are needed to cover
        // half of 23, but they share a length, so seqkit reports 1. Every
        // length in the project fixtures is distinct, so only a case like this
        // separates the two definitions.
        let stats = SequenceStats::from_records(&records(&[3, 10, 10]));
        assert_eq!(stats.n50, 10);
        assert_eq!(stats.n50_num, 1, "not 2");
    }

    #[test]
    fn n50_is_the_length_covering_half_the_residues() {
        // 755 + 584 + 526 = 1865, more than half of 3284.
        let stats = SequenceStats::from_records(&records(&[
            81, 140, 157, 189, 198, 273, 381, 526, 584, 755,
        ]));
        assert_eq!(stats.n50, 526);
        assert_eq!(stats.n50_num, 3);
    }

    #[test]
    fn gaps_are_counted_but_still_appear_in_the_composition() {
        let records = vec![Record {
            id: "a".into(),
            residues: b"MK-T.A".to_vec(),
        }];
        let stats = SequenceStats::from_records(&records);
        assert_eq!(stats.gaps, 2);
        assert_eq!(
            stats.total, 6,
            "gaps count towards the length, as seqkit does"
        );
        assert_eq!(stats.composition.get(&b'M'), Some(&1));
    }

    #[test]
    fn an_empty_input_does_not_divide_by_zero() {
        let stats = SequenceStats::from_records(&[]);
        assert_eq!(stats.count, 0);
        assert_eq!(stats.mean, 0.0);
        assert_eq!(stats.n50, 0);
        assert_eq!(stats.fraction(b'A'), 0.0);
    }
}
