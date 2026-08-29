//! RSeQC `read_GC.py` reimplementation.
//!
//! Reports the distribution of GC content across reads, which is how library
//! preparation bias and contamination show up before any annotation is
//! involved.
//!
//! # Upstream semantics
//!
//! Taken from RSeQC 5.0.4's `qcmodule.SAM.readGC`. A read is skipped when it
//! is unmapped, QC-failed, or below the mapping quality cutoff, which defaults
//! to 30. Secondary alignments are not excluded.
//!
//! GC percentage is `(C + G)` over the **full sequence length**, so an `N` pads
//! the denominator without ever counting as GC and pulls the percentage down.
//! Dividing by the called bases instead is the more defensible measurement and
//! the wrong one to report under this name. The result is formatted to two
//! decimals, and reads sharing a formatted value share a bin.

use std::collections::BTreeMap;

use rust_htslib::bam;

use crate::common::bam_flags::*;

/// The histogram key for a percentage: the value scaled by 100 and rounded
/// half to even.
///
/// RSeQC formats with `%4.2f`, which rounds ties to the even digit. A read
/// with 9 of 32 bases G or C is exactly 28.125 percent and lands in the 28.12
/// bin, not 28.13. Rounding away from zero puts it in the wrong bin and, on
/// the project fixture, moves fourteen reads.
fn bin_key(percent: f64) -> u32 {
    let scaled = percent * 100.0;
    let floor = scaled.floor();
    let value = if (scaled - floor - 0.5).abs() < 1e-9 {
        if (floor as i64) % 2 == 0 {
            floor
        } else {
            floor + 1.0
        }
    } else {
        scaled.round()
    };
    value.max(0.0) as u32
}

/// Accumulates the read GC distribution.
///
/// Percentages are keyed by their value scaled by 100 and rounded, so the map
/// is exact rather than relying on float keys.
#[derive(Debug, Default, Clone)]
pub struct ReadGcAccum {
    /// Read counts keyed by GC percentage times 100.
    counts: BTreeMap<u32, u64>,
    /// Mapped reads seen.
    pub reads: u64,
    /// Reads whose sequence held no called base at all.
    pub reads_without_bases: u64,
}

impl ReadGcAccum {
    /// A new, empty accumulator.
    pub fn new() -> Self {
        Self::default()
    }

    /// Offer one record.
    pub fn process_read(&mut self, record: &bam::Record, mapq_cut: u8) {
        let flags = record.flags();
        if flags & BAM_FUNMAP != 0 || flags & BAM_FQCFAIL != 0 || record.mapq() < mapq_cut {
            return;
        }
        self.reads += 1;

        let sequence = record.seq().as_bytes();
        if sequence.is_empty() {
            self.reads_without_bases += 1;
            return;
        }
        let gc = sequence
            .iter()
            .filter(|b| matches!(b.to_ascii_uppercase(), b'G' | b'C'))
            .count() as f64;

        let percent = gc / sequence.len() as f64 * 100.0;
        *self.counts.entry(bin_key(percent)).or_insert(0) += 1;
    }

    /// Fold another accumulator in.
    pub fn merge(&mut self, other: ReadGcAccum) {
        for (key, count) in other.counts {
            *self.counts.entry(key).or_insert(0) += count;
        }
        self.reads += other.reads;
        self.reads_without_bases += other.reads_without_bases;
    }

    /// The distribution as `(percentage, count)` pairs, ascending.
    pub fn distribution(&self) -> Vec<(f64, u64)> {
        self.counts
            .iter()
            .map(|(key, count)| (f64::from(*key) / 100.0, *count))
            .collect()
    }

    /// Reads that landed in the histogram.
    pub fn counted(&self) -> u64 {
        self.counts.values().sum()
    }

    /// Mean GC percentage across the counted reads.
    pub fn mean(&self) -> f64 {
        let total = self.counted();
        if total == 0 {
            return 0.0;
        }
        self.counts
            .iter()
            .map(|(key, count)| f64::from(*key) / 100.0 * *count as f64)
            .sum::<f64>()
            / total as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString, Record};

    /// RSeQC's default mapping quality cutoff.
    const CUTOFF: u8 = 30;

    fn record(sequence: &[u8], flags: u16, mapq: u8) -> Record {
        let mut r = Record::new();
        let cigar = CigarString(vec![Cigar::Match(sequence.len() as u32)]);
        let quality = vec![30u8; sequence.len()];
        r.set(b"q", Some(&cigar), sequence, &quality);
        r.set_flags(flags);
        r.set_mapq(mapq);
        r
    }

    #[test]
    fn n_pads_the_denominator_rather_than_being_excluded() {
        let mut accum = ReadGcAccum::new();
        // Two of five bases are G or C, and the N counts towards the length,
        // so this is 40 percent rather than the 50 that dropping it would give.
        accum.process_read(&record(b"ACGTN", 0, 60), CUTOFF);
        assert_eq!(accum.distribution(), vec![(40.0, 1)]);
    }

    #[test]
    fn unmapped_and_qc_failed_reads_are_skipped() {
        let mut accum = ReadGcAccum::new();
        accum.process_read(&record(b"GGGG", BAM_FUNMAP, 60), CUTOFF);
        accum.process_read(&record(b"GGGG", BAM_FQCFAIL, 60), CUTOFF);
        assert_eq!(accum.reads, 0);
        assert!(accum.distribution().is_empty());
    }

    #[test]
    fn reads_below_the_mapping_quality_cutoff_are_skipped() {
        let mut accum = ReadGcAccum::new();
        accum.process_read(&record(b"GGCC", 0, 29), CUTOFF);
        assert_eq!(accum.reads, 0, "29 is below the default cutoff of 30");
        accum.process_read(&record(b"GGCC", 0, 30), CUTOFF);
        assert_eq!(accum.reads, 1, "30 itself passes");
    }

    #[test]
    fn secondary_alignments_still_count() {
        let mut accum = ReadGcAccum::new();
        accum.process_read(&record(b"GGCC", BAM_FSECONDARY, 60), CUTOFF);
        assert_eq!(
            accum.reads, 1,
            "RSeQC excludes neither secondaries nor duplicates"
        );
        assert_eq!(accum.distribution(), vec![(100.0, 1)]);
    }

    #[test]
    fn ties_round_to_the_even_digit_as_printf_does() {
        // 9 of 32 bases is exactly 28.125 percent, which "%4.2f" writes as
        // 28.12; rounding away from zero would give 28.13.
        let mut sequence = vec![b'G'; 9];
        sequence.extend(std::iter::repeat_n(b'A', 23));
        let mut accum = ReadGcAccum::new();
        accum.process_read(&record(&sequence, 0, 60), CUTOFF);
        assert_eq!(accum.distribution(), vec![(28.12, 1)]);
    }

    #[test]
    fn percentages_are_kept_to_two_decimals() {
        let mut accum = ReadGcAccum::new();
        // One G in three bases is 33.333..., reported as 33.33.
        accum.process_read(&record(b"AGT", 0, 60), CUTOFF);
        assert_eq!(accum.distribution(), vec![(33.33, 1)]);
    }

    #[test]
    fn merging_adds_the_distributions() {
        let mut a = ReadGcAccum::new();
        a.process_read(&record(b"GGCC", 0, 60), CUTOFF);
        let mut b = ReadGcAccum::new();
        b.process_read(&record(b"GGCC", 0, 60), CUTOFF);
        b.process_read(&record(b"AATT", 0, 60), CUTOFF);
        a.merge(b);
        assert_eq!(a.distribution(), vec![(0.0, 1), (100.0, 2)]);
        assert_eq!(a.reads, 3);
    }

    #[test]
    fn the_mean_weights_by_read_count() {
        let mut accum = ReadGcAccum::new();
        accum.process_read(&record(b"GGCC", 0, 60), CUTOFF);
        accum.process_read(&record(b"AATT", 0, 60), CUTOFF);
        accum.process_read(&record(b"AATT", 0, 60), CUTOFF);
        assert!((accum.mean() - 100.0 / 3.0).abs() < 1e-9);
    }
}
