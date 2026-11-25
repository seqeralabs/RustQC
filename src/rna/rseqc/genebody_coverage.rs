//! RSeQC `geneBody_coverage.py` reimplementation.
//!
//! Measures how evenly reads cover the gene body from the 5' to the 3' end,
//! which is the standard way RNA degradation and 3' bias show up.
//!
//! # Upstream semantics
//!
//! Taken from RSeQC 5.0.4's `geneBody_coverage.py`, whose behaviour is not
//! obvious from its output:
//!
//! Each transcript's mRNA bases are laid out end to end, exon by exon, in
//! one-based genome coordinates. Transcripts shorter than
//! [`DEFAULT_MIN_MRNA_LENGTH`] are skipped entirely. The remainder are reduced
//! to exactly 100 positions by [`percentile_positions`], which is a linear
//! interpolation between neighbouring bases, rounded. Those are *positions*,
//! not windows: coverage is read at 100 single bases, and everything between
//! them is never looked at.
//!
//! A base counts a read when the read is not deleted at that position and is
//! not QC-failed, secondary, unmapped or duplicate-flagged.
//!
//! A transcript on the minus strand has its 100 values reversed, so index 0 is
//! always the 5' end. The reversal happens per transcript, before summing, so
//! a gene set mixing strands still aggregates 5' to 3'.

use std::collections::BTreeMap;

/// Transcripts with fewer mRNA bases than this are skipped.
pub const DEFAULT_MIN_MRNA_LENGTH: usize = 100;

/// Number of points each transcript is reduced to.
pub const POINTS: usize = 100;

/// One transcript reduced to the positions its coverage is read at.
#[derive(Debug, Clone)]
pub struct TranscriptPoints {
    /// Contig the transcript sits on.
    pub chrom: String,
    /// `+` or `-`.
    pub reverse: bool,
    /// One-based genome positions, in the order the exons were concatenated.
    /// Ascending for a normal transcript, but not for one whose exons
    /// overlap, because RSeQC never sorts them.
    pub positions: Vec<i64>,
}

/// Round to the nearest integer, ties going to the even one.
///
/// RSeQC interpolates with Python's `round`, which rounds half to even.
/// Rounding away from zero instead moves an interpolated position by one base
/// wherever the interpolation lands exactly halfway, which on the project
/// fixture changes the 23rd percentile's coverage by a whole read.
fn round_half_to_even(value: f64) -> i64 {
    let floor = value.floor();
    if (value - floor - 0.5).abs() < 1e-9 {
        if (floor as i64) % 2 == 0 {
            floor as i64
        } else {
            floor as i64 + 1
        }
    } else {
        value.round() as i64
    }
}

/// Reduce a sorted list of values to 100 percentile points.
///
/// This is RSeQC's `mystat.percentile_list`. A list shorter than 100 is
/// returned unchanged rather than padded, which is why transcripts below the
/// minimum length are dropped beforehand: they would otherwise contribute
/// fewer than 100 values and shift every later index.
pub fn percentile_positions(sorted: &[i64]) -> Vec<i64> {
    if sorted.is_empty() {
        return Vec::new();
    }
    if sorted.len() < POINTS {
        return sorted.to_vec();
    }
    let last = (sorted.len() - 1) as f64;
    (1..=POINTS)
        .map(|i| {
            let k = last * i as f64 / 100.0;
            let floor = k.floor();
            let ceil = k.ceil();
            if (floor - ceil).abs() < f64::EPSILON {
                sorted[k as usize]
            } else {
                let low = sorted[floor as usize] as f64 * (ceil - k);
                let high = sorted[ceil as usize] as f64 * (k - floor);
                round_half_to_even(low + high)
            }
        })
        .collect()
}

/// Build the percentile points for one transcript from its exons.
///
/// Exons are half-open zero-based, as BED and the GTF parser give them; the
/// positions returned are one-based, as RSeQC works in, and are in the order
/// the exons were given rather than sorted.
pub fn transcript_points(
    chrom: &str,
    reverse: bool,
    exons: &[(i64, i64)],
    min_length: usize,
) -> Option<TranscriptPoints> {
    // Exons are concatenated in the order given and deliberately **not**
    // sorted. RSeQC does not sort either, and its percentile interpolation
    // assumes a sorted list, so a transcript with overlapping exons feeds it a
    // list that dips backwards. The interpolated positions then repeat, and
    // the curve for that transcript is computed over fewer than 100 distinct
    // bases. Sorting here would be the more sensible thing and would disagree
    // with the tool this reproduces; see `GeneBodyCoverage::add_transcript`
    // for what the repeats then do.
    let mut bases: Vec<i64> = Vec::new();
    for (start, end) in exons {
        bases.extend((*start + 1)..=*end);
    }
    if bases.len() < min_length {
        return None;
    }
    Some(TranscriptPoints {
        chrom: chrom.to_string(),
        reverse,
        positions: percentile_positions(&bases),
    })
}

/// Aggregated coverage across every transcript, 5' to 3'.
#[derive(Debug, Clone)]
pub struct GeneBodyCoverage {
    /// One total per percentile point.
    pub totals: Vec<u64>,
    /// Transcripts that contributed.
    pub transcripts: u64,
}

impl Default for GeneBodyCoverage {
    fn default() -> Self {
        Self {
            totals: vec![0; POINTS],
            transcripts: 0,
        }
    }
}

impl GeneBodyCoverage {
    /// Fold one transcript's per-position coverage in.
    ///
    /// Only the **distinct** positions contribute, which matters on short
    /// transcripts where the interpolation repeats a base. RSeQC keeps its
    /// coverage in a dictionary keyed by position, so repeats collapse and the
    /// list it aggregates is shorter than 100; it then adds those values at
    /// indices 0 upwards, which quietly compresses that transcript's curve
    /// towards the 5' end.
    ///
    /// That is a defect in RSeQC rather than a design decision, and it is
    /// reproduced here on purpose: the point of this module is to agree with
    /// the tool people already have. On the project fixture one gene of eight
    /// has 100 points across 88 distinct positions, and correcting the
    /// behaviour moves the aggregate curve. Anyone deciding to diverge should
    /// do so knowingly, which is why this is spelled out rather than hidden.
    pub fn add_transcript(&mut self, points: &TranscriptPoints, coverage: &BTreeMap<i64, u64>) {
        let mut distinct: Vec<i64> = points.positions.clone();
        distinct.sort_unstable();
        distinct.dedup();

        let mut values: Vec<u64> = distinct
            .iter()
            .map(|p| coverage.get(p).copied().unwrap_or(0))
            .collect();
        if points.reverse {
            values.reverse();
        }
        for (index, value) in values.iter().enumerate() {
            if index < self.totals.len() {
                self.totals[index] += value;
            }
        }
        self.transcripts += 1;
    }

    /// Pearson's moment coefficient of skewness over the aggregated curve,
    /// which RSeQC reports alongside it.
    pub fn skewness(&self) -> f64 {
        let n = self.totals.len();
        if n < 2 {
            return 0.0;
        }
        let values: Vec<f64> = self.totals.iter().map(|v| *v as f64).collect();
        let mean = values.iter().sum::<f64>() / n as f64;
        let variance = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n as f64;
        let sd = variance.sqrt();
        if sd == 0.0 {
            return 0.0;
        }
        values
            .iter()
            .map(|v| ((v - mean) / sd).powi(3))
            .sum::<f64>()
            / n as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_list_shorter_than_a_hundred_is_returned_unchanged() {
        let values = vec![1i64, 2, 3];
        assert_eq!(percentile_positions(&values), values);
    }

    #[test]
    fn a_hundred_values_map_to_themselves_offset_by_one() {
        // With exactly 100 values the interpolation lands on each in turn,
        // starting at the second, because the points run from 1 to 100 rather
        // than 0 to 99.
        let values: Vec<i64> = (1..=100).collect();
        let points = percentile_positions(&values);
        assert_eq!(points.len(), 100);
        assert_eq!(points[0], 2);
        assert_eq!(points[99], 100);
    }

    #[test]
    fn interpolation_rounds_between_neighbours() {
        let values: Vec<i64> = (0..200).map(|i| i * 10).collect();
        let points = percentile_positions(&values);
        assert_eq!(points.len(), 100);
        // Ascending and inside the input range.
        assert!(points.windows(2).all(|w| w[0] <= w[1]));
        assert!(*points.first().unwrap() >= 0);
        assert!(*points.last().unwrap() <= 1990);
    }

    #[test]
    fn exons_are_laid_out_end_to_end_in_one_based_coordinates() {
        // Two exons of 60 bases each: 120 mRNA bases, so the transcript is
        // kept and reduced to 100 points.
        let points = transcript_points("chr1", false, &[(0, 60), (100, 160)], 100).unwrap();
        assert_eq!(points.positions.len(), 100);
        assert_eq!(*points.positions.first().unwrap(), 2, "one-based");
        // Non-overlapping exons concatenate in ascending order, so these do
        // come out sorted; overlapping ones do not.
        assert!(points.positions.windows(2).all(|w| w[0] <= w[1]));
    }

    #[test]
    fn overlapping_exons_leave_the_positions_unsorted() {
        // Two exons covering 1..=600 and 401..=900. Concatenated in order the
        // base list dips backwards at the junction, and because RSeQC never
        // sorts it the interpolated positions dip too. That is the mechanism
        // behind the repeated positions seen on real annotations; whether a
        // given transcript actually repeats one depends on where the
        // interpolation lands.
        let points = transcript_points("chr1", false, &[(0, 600), (400, 900)], 100).unwrap();
        assert_eq!(points.positions.len(), 100);
        assert!(
            points.positions.windows(2).any(|w| w[0] > w[1]),
            "the concatenation is not sorted, so the points should dip"
        );
    }

    #[test]
    fn interpolation_can_land_inside_an_intron() {
        // Two exons either side of a 40 base intron. The percentile points are
        // interpolated between neighbouring *mRNA* bases, and the pair
        // straddling the junction interpolates between genome positions 60 and
        // 101, so a point can fall in the intron and coverage is then read at
        // a base the transcript does not contain.
        //
        // This is RSeQC's behaviour, not a defect introduced here, and
        // reproducing it is the point. It is asserted so that anyone
        // "correcting" it later sees the intent.
        let points = transcript_points("chr1", false, &[(0, 60), (100, 160)], 100).unwrap();
        let intronic: Vec<i64> = points
            .positions
            .iter()
            .copied()
            .filter(|p| (61..=100).contains(p))
            .collect();
        assert!(
            !intronic.is_empty(),
            "expected at least one interpolated point in the intron"
        );
    }

    #[test]
    fn a_transcript_below_the_minimum_length_is_dropped() {
        assert!(transcript_points("chr1", false, &[(0, 50)], 100).is_none());
        assert!(transcript_points("chr1", false, &[(0, 100)], 100).is_some());
    }

    #[test]
    fn a_minus_strand_transcript_is_reversed_before_aggregating() {
        let points = TranscriptPoints {
            chrom: "chr1".into(),
            reverse: true,
            positions: vec![10, 20, 30],
        };
        let coverage: BTreeMap<i64, u64> = [(10, 1), (20, 2), (30, 3)].into_iter().collect();
        let mut aggregate = GeneBodyCoverage {
            totals: vec![0; 3],
            transcripts: 0,
        };
        aggregate.add_transcript(&points, &coverage);
        assert_eq!(
            aggregate.totals,
            vec![3, 2, 1],
            "index 0 is the 5' end, which for a minus-strand gene is the last position"
        );
    }

    #[test]
    fn a_plus_strand_transcript_keeps_its_order() {
        let points = TranscriptPoints {
            chrom: "chr1".into(),
            reverse: false,
            positions: vec![10, 20, 30],
        };
        let coverage: BTreeMap<i64, u64> = [(10, 1), (20, 2), (30, 3)].into_iter().collect();
        let mut aggregate = GeneBodyCoverage {
            totals: vec![0; 3],
            transcripts: 0,
        };
        aggregate.add_transcript(&points, &coverage);
        assert_eq!(aggregate.totals, vec![1, 2, 3]);
    }

    #[test]
    fn an_uncovered_position_contributes_zero_rather_than_being_skipped() {
        let points = TranscriptPoints {
            chrom: "chr1".into(),
            reverse: false,
            positions: vec![10, 20, 30],
        };
        let coverage: BTreeMap<i64, u64> = [(10, 5)].into_iter().collect();
        let mut aggregate = GeneBodyCoverage {
            totals: vec![0; 3],
            transcripts: 0,
        };
        aggregate.add_transcript(&points, &coverage);
        assert_eq!(aggregate.totals, vec![5, 0, 0]);
    }

    #[test]
    fn interpolation_ties_round_to_the_even_position() {
        assert_eq!(round_half_to_even(10.5), 10, "10 is even");
        assert_eq!(round_half_to_even(11.5), 12, "12 is even");
        assert_eq!(round_half_to_even(10.4), 10);
        assert_eq!(round_half_to_even(10.6), 11);
    }

    #[test]
    fn repeated_positions_collapse_and_shorten_the_curve() {
        // Three points across two distinct positions: RSeQC aggregates two
        // values at indices 0 and 1, leaving index 2 untouched.
        let points = TranscriptPoints {
            chrom: "chr1".into(),
            reverse: false,
            positions: vec![10, 10, 20],
        };
        let coverage: BTreeMap<i64, u64> = [(10, 5), (20, 7)].into_iter().collect();
        let mut aggregate = GeneBodyCoverage {
            totals: vec![0; 3],
            transcripts: 0,
        };
        aggregate.add_transcript(&points, &coverage);
        assert_eq!(
            aggregate.totals,
            vec![5, 7, 0],
            "the repeat collapses rather than contributing twice"
        );
    }

    #[test]
    fn a_flat_curve_has_no_skew() {
        let aggregate = GeneBodyCoverage {
            totals: vec![10; 100],
            transcripts: 1,
        };
        assert_eq!(aggregate.skewness(), 0.0);
    }
}
