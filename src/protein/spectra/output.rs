//! Writer for the mass spectrometry run report.

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};

use super::metrics::SpectraMetrics;

/// Write the run report.
pub fn write_report(file: &str, metrics: &SpectraMetrics, path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .map(std::io::BufWriter::new)
        .with_context(|| format!("Failed to create spectra report: {}", path.display()))?;

    writeln!(out, "# RustQC protein spectra report")?;
    writeln!(out, "# file\t{file}")?;
    writeln!(out)?;

    writeln!(out, "## Run")?;
    writeln!(out, "metric\tvalue")?;
    writeln!(out, "spectra\t{}", metrics.total_spectra())?;
    writeln!(out, "peaks\t{}", metrics.total_peaks())?;
    writeln!(out, "rt_min\t{:.6}", metrics.rt_min)?;
    writeln!(out, "rt_max\t{:.6}", metrics.rt_max)?;
    writeln!(out, "rt_span\t{:.6}", metrics.rt_span())?;
    match metrics.ms2_per_ms1() {
        Some(ratio) => writeln!(out, "ms2_per_ms1\t{ratio:.6}")?,
        // Written as NA rather than 0, which would read as "no fragmentation"
        // when the truth is "no survey scans to divide by".
        None => writeln!(out, "ms2_per_ms1\tNA")?,
    }
    writeln!(out)?;

    writeln!(out, "## Levels")?;
    writeln!(
        out,
        "ms_level\tspectra\tpeaks\tmean_peaks\tmin_peaks\tmax_peaks\ttotal_ion_current"
    )?;
    for (level, m) in &metrics.levels {
        writeln!(
            out,
            "{level}\t{}\t{}\t{:.4}\t{}\t{}\t{:.4}",
            m.spectra,
            m.peaks,
            m.mean_peaks(),
            m.min_peaks,
            m.max_peaks,
            m.total_ion_current,
        )?;
    }
    writeln!(out)?;

    writeln!(out, "## Precursors")?;
    writeln!(out, "metric\tvalue")?;
    writeln!(out, "selected\t{}", metrics.precursors)?;
    writeln!(out, "without_charge\t{}", metrics.precursors_without_charge)?;
    if metrics.precursors > 0 {
        writeln!(out, "mz_min\t{:.4}", metrics.precursor_mz_min)?;
        writeln!(out, "mz_max\t{:.4}", metrics.precursor_mz_max)?;
    }
    writeln!(out)?;

    if !metrics.precursor_charges.is_empty() {
        writeln!(out, "## Charge states")?;
        writeln!(out, "charge\tcount\tfraction")?;
        let assigned: u64 = metrics.precursor_charges.values().sum();
        for (charge, count) in &metrics.precursor_charges {
            writeln!(
                out,
                "{charge}\t{count}\t{:.6}",
                *count as f64 / assigned as f64
            )?;
        }
        writeln!(out)?;
    }

    out.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scratch(name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join("rustqc-spectra-output-tests");
        std::fs::create_dir_all(&dir).unwrap();
        dir.join(name)
    }

    #[test]
    fn a_run_without_survey_scans_reports_na_rather_than_zero() {
        let mut metrics = SpectraMetrics::default();
        metrics.observe(2, 10, 1.0, 0.0, None);
        metrics.finish();
        let path = scratch("ms2only.txt");
        write_report("f.mzML", &metrics, &path).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(
            text.contains("ms2_per_ms1\tNA"),
            "0 would read as no fragmentation, which is the opposite of the truth"
        );
    }

    #[test]
    fn the_charge_section_is_omitted_when_nothing_was_assigned() {
        let mut metrics = SpectraMetrics::default();
        metrics.observe(2, 10, 1.0, 0.0, Some((500.0, None)));
        metrics.finish();
        let path = scratch("nocharge.txt");
        write_report("f.mzML", &metrics, &path).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(text.contains("without_charge\t1"));
        assert!(!text.contains("## Charge states"));
    }
}
