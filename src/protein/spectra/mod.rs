//! Mass spectrometry run quality control from mzML.
//!
//! Answers the question a QC report on a raw run should answer: how many
//! spectra were acquired at each level, how much signal they carry, over what
//! retention time window, and what the precursors selected for fragmentation
//! looked like.
//!
//! Reading is delegated to [`mzdata`], which is why this module sits behind
//! the `proteomics` cargo feature. Building without it drops the `spectra`
//! mode rather than failing at run time.

pub mod metrics;
pub mod output;

use std::path::Path;

use anyhow::{Context, Result};
use mzdata::prelude::*;
use mzdata::MZReader;

use metrics::SpectraMetrics;

/// Read an mzML file and summarise it.
pub fn analyse(path: &Path) -> Result<SpectraMetrics> {
    let reader = MZReader::open_path(path)
        .with_context(|| format!("Failed to open spectra file: {}", path.display()))?;

    let mut metrics = SpectraMetrics::default();
    for spectrum in reader {
        let level = spectrum.ms_level();
        let peaks = spectrum.peaks().len() as u64;

        // Total ion current as the sum of the peak intensities actually
        // present, rather than whatever the instrument wrote in the header.
        // The two usually agree; where they do not, the header can describe a
        // profile spectrum that has since been centroided, and the peaks are
        // the honest answer for a file as it stands.
        let tic = f64::from(spectrum.peaks().tic());

        let start_time = spectrum.start_time();

        let precursor = spectrum
            .precursor()
            .and_then(|p| p.ion())
            .map(|ion| (ion.mz, ion.charge));

        metrics.observe(level, peaks, tic, start_time, precursor);
    }

    metrics.finish();
    Ok(metrics)
}
