//! Per-run and per-level spectrum metrics.

use std::collections::BTreeMap;

/// Metrics for one MS level.
#[derive(Debug, Clone, Default)]
pub struct LevelMetrics {
    /// Spectra acquired at this level.
    pub spectra: u64,
    /// Peaks across those spectra.
    pub peaks: u64,
    /// Summed total ion current.
    pub total_ion_current: f64,
    /// Fewest peaks in any one spectrum.
    pub min_peaks: u64,
    /// Most peaks in any one spectrum.
    pub max_peaks: u64,
}

impl LevelMetrics {
    /// Mean peaks per spectrum.
    pub fn mean_peaks(&self) -> f64 {
        if self.spectra == 0 {
            0.0
        } else {
            self.peaks as f64 / self.spectra as f64
        }
    }
}

/// Everything the `spectra` mode reports about one run.
#[derive(Debug, Clone, Default)]
pub struct SpectraMetrics {
    /// Per-level metrics, keyed by MS level.
    pub levels: BTreeMap<u8, LevelMetrics>,
    /// Earliest retention time seen, in minutes.
    pub rt_min: f64,
    /// Latest retention time seen, in minutes.
    pub rt_max: f64,
    /// Charge state distribution of the precursors selected for fragmentation.
    pub precursor_charges: BTreeMap<i32, u64>,
    /// Precursors whose charge the instrument did not assign.
    pub precursors_without_charge: u64,
    /// Lowest precursor m/z selected.
    pub precursor_mz_min: f64,
    /// Highest precursor m/z selected.
    pub precursor_mz_max: f64,
    /// Precursors selected in total.
    pub precursors: u64,
    /// Whether anything has been observed yet, so the ranges know to
    /// initialise rather than compare against a default.
    seen: bool,
}

impl SpectraMetrics {
    /// Fold one spectrum in.
    pub fn observe(
        &mut self,
        level: u8,
        peaks: u64,
        total_ion_current: f64,
        start_time: f64,
        precursor: Option<(f64, Option<i32>)>,
    ) {
        let entry = self.levels.entry(level).or_default();
        entry.spectra += 1;
        entry.peaks += peaks;
        entry.total_ion_current += total_ion_current;
        if entry.spectra == 1 {
            entry.min_peaks = peaks;
            entry.max_peaks = peaks;
        } else {
            entry.min_peaks = entry.min_peaks.min(peaks);
            entry.max_peaks = entry.max_peaks.max(peaks);
        }

        if !self.seen {
            self.rt_min = start_time;
            self.rt_max = start_time;
            self.seen = true;
        } else {
            self.rt_min = self.rt_min.min(start_time);
            self.rt_max = self.rt_max.max(start_time);
        }

        if let Some((mz, charge)) = precursor {
            if self.precursors == 0 {
                self.precursor_mz_min = mz;
                self.precursor_mz_max = mz;
            } else {
                self.precursor_mz_min = self.precursor_mz_min.min(mz);
                self.precursor_mz_max = self.precursor_mz_max.max(mz);
            }
            self.precursors += 1;
            match charge {
                Some(z) => *self.precursor_charges.entry(z).or_insert(0) += 1,
                None => self.precursors_without_charge += 1,
            }
        }
    }

    /// Called once every spectrum has been offered.
    pub fn finish(&mut self) {
        if !self.seen {
            self.rt_min = 0.0;
            self.rt_max = 0.0;
        }
    }

    /// Spectra across every level.
    pub fn total_spectra(&self) -> u64 {
        self.levels.values().map(|l| l.spectra).sum()
    }

    /// Peaks across every level.
    pub fn total_peaks(&self) -> u64 {
        self.levels.values().map(|l| l.peaks).sum()
    }

    /// Retention time span, in minutes.
    pub fn rt_span(&self) -> f64 {
        self.rt_max - self.rt_min
    }

    /// Ratio of fragmentation spectra to survey spectra, the usual measure of
    /// how hard the instrument was working.
    ///
    /// `None` when there are no MS1 spectra to divide by, which is what an
    /// MS2-only file gives.
    pub fn ms2_per_ms1(&self) -> Option<f64> {
        let ms1 = self.levels.get(&1)?.spectra;
        if ms1 == 0 {
            return None;
        }
        let ms2 = self.levels.get(&2).map(|l| l.spectra).unwrap_or(0);
        Some(ms2 as f64 / ms1 as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn levels_accumulate_separately() {
        let mut m = SpectraMetrics::default();
        m.observe(1, 100, 10.0, 0.5, None);
        m.observe(2, 20, 2.0, 0.6, Some((500.0, Some(2))));
        m.observe(1, 200, 20.0, 0.7, None);
        m.finish();

        assert_eq!(m.levels[&1].spectra, 2);
        assert_eq!(m.levels[&1].peaks, 300);
        assert_eq!(m.levels[&2].spectra, 1);
        assert_eq!(m.total_spectra(), 3);
        assert_eq!(m.total_peaks(), 320);
    }

    #[test]
    fn peak_bounds_start_from_the_first_spectrum_not_from_zero() {
        let mut m = SpectraMetrics::default();
        m.observe(1, 50, 1.0, 0.0, None);
        m.observe(1, 90, 1.0, 0.1, None);
        m.finish();
        assert_eq!(m.levels[&1].min_peaks, 50, "not 0");
        assert_eq!(m.levels[&1].max_peaks, 90);
    }

    #[test]
    fn the_retention_time_range_spans_what_was_seen() {
        let mut m = SpectraMetrics::default();
        m.observe(1, 1, 1.0, 5.0, None);
        m.observe(1, 1, 1.0, 2.0, None);
        m.observe(1, 1, 1.0, 9.0, None);
        m.finish();
        assert_eq!(m.rt_min, 2.0, "not 0.0 from the default");
        assert_eq!(m.rt_max, 9.0);
        assert_eq!(m.rt_span(), 7.0);
    }

    #[test]
    fn precursors_without_a_charge_are_counted_apart() {
        let mut m = SpectraMetrics::default();
        m.observe(2, 1, 1.0, 0.0, Some((400.0, Some(2))));
        m.observe(2, 1, 1.0, 0.1, Some((600.0, None)));
        m.observe(2, 1, 1.0, 0.2, Some((500.0, Some(2))));
        m.finish();
        assert_eq!(m.precursors, 3);
        assert_eq!(m.precursor_charges[&2], 2);
        assert_eq!(m.precursors_without_charge, 1);
        assert_eq!(m.precursor_mz_min, 400.0);
        assert_eq!(m.precursor_mz_max, 600.0);
    }

    #[test]
    fn the_fragmentation_ratio_needs_survey_spectra_to_divide_by() {
        let mut m = SpectraMetrics::default();
        m.observe(2, 1, 1.0, 0.0, None);
        m.finish();
        assert_eq!(m.ms2_per_ms1(), None, "an MS2-only file has no ratio");

        let mut m = SpectraMetrics::default();
        m.observe(1, 1, 1.0, 0.0, None);
        m.observe(2, 1, 1.0, 0.1, None);
        m.observe(2, 1, 1.0, 0.2, None);
        m.finish();
        assert_eq!(m.ms2_per_ms1(), Some(2.0));
    }

    #[test]
    fn an_empty_run_reports_zeroes_rather_than_a_default_range() {
        let mut m = SpectraMetrics::default();
        m.finish();
        assert_eq!(m.total_spectra(), 0);
        assert_eq!(m.rt_span(), 0.0);
        assert_eq!(m.ms2_per_ms1(), None);
    }
}
