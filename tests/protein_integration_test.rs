//! Parity tests for the protein pipeline against the committed reference
//! outputs.
//!
//! The fixtures under `tests/expected/protein/` are the output of the upstream
//! tools themselves, at the versions pinned in `VERSIONS.txt`. A failure here
//! is a defect in RustQC, not a reason to regenerate the fixture.

use std::path::{Path, PathBuf};

use rustqc::protein::sequence::{self, defects, output, stats::SequenceStats};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/expected/protein")
        .join(name)
}

fn input(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/protein")
        .join(name)
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join("rustqc-protein-parity");
    std::fs::create_dir_all(&dir).unwrap();
    dir.join(name)
}

/// The fixtures record which seqkit produced them; a mismatch is a version
/// skew rather than a parity failure, so it is reported as such.
#[test]
fn fixture_tool_versions_are_the_pinned_ones() {
    let versions = std::fs::read_to_string(fixture("VERSIONS.txt")).unwrap();
    assert!(
        versions.contains("seqkit\t2.13.0"),
        "unexpected seqkit fixture version: {versions}"
    );
}

/// The whole `seqkit stats -a -T` table, both fixtures, byte for byte.
#[test]
fn sequence_stats_match_seqkit() {
    for name in ["yeast_UPS_mini.fasta", "protein_mini_with_cazymes.faa"] {
        let records = sequence::read_fasta(&input(name)).unwrap();
        let stats = SequenceStats::from_records(&records);
        let path = scratch(&format!("{name}.tsv"));
        output::write_seqkit_stats(&[(name.to_string(), stats)], &path).unwrap();

        let stem = name.rsplit_once('.').unwrap().0;
        let got = std::fs::read_to_string(&path).unwrap();
        let want = std::fs::read_to_string(fixture(&format!("{stem}.seqkit.tsv"))).unwrap();

        let got_lines: Vec<&str> = got.lines().collect();
        let want_lines: Vec<&str> = want.lines().collect();
        assert_eq!(got_lines.len(), want_lines.len(), "{name}: line count");
        for (i, (a, b)) in got_lines.iter().zip(&want_lines).enumerate() {
            assert_eq!(a, b, "{name}: line {} differs", i + 1);
        }
    }
}

/// The two conventions that had to be recovered from seqkit's behaviour, since
/// neither is what a statistics library gives by default.
#[test]
fn quartiles_use_tukeys_halves_with_bankers_rounding() {
    let records = sequence::read_fasta(&input("yeast_UPS_mini.fasta")).unwrap();
    let stats = SequenceStats::from_records(&records);
    // Lengths are 81 140 157 189 198 273 381 526 584 755.
    assert_eq!(stats.q1, 157, "linear interpolation would give 165");
    assert_eq!(
        stats.q2, 236,
        "the median is 235.5, rounded to the even 236"
    );
    assert_eq!(stats.q3, 526);

    let records = sequence::read_fasta(&input("protein_mini_with_cazymes.faa")).unwrap();
    let stats = SequenceStats::from_records(&records);
    assert_eq!(stats.q2, 376, "376.5 rounds down, to the even 376");
    assert_eq!(stats.q3, 516, "516.5 likewise");
}

/// Composition is what seqkit does not report, so it is checked against the
/// input directly: every residue counted once, summing to the total length.
#[test]
fn composition_accounts_for_every_residue() {
    for name in ["yeast_UPS_mini.fasta", "protein_mini_with_cazymes.faa"] {
        let records = sequence::read_fasta(&input(name)).unwrap();
        let stats = SequenceStats::from_records(&records);
        let counted: u64 = stats.composition.values().sum();
        assert_eq!(
            counted, stats.total,
            "{name}: composition must account for every residue"
        );
        let fractions: f64 = stats.composition.keys().map(|r| stats.fraction(*r)).sum();
        assert!(
            (fractions - 1.0).abs() < 1e-12,
            "{name}: fractions sum to {fractions}"
        );
    }
}

/// Both fixtures are real reference proteomes, so they should be clean.
#[test]
fn the_reference_proteomes_carry_no_defects() {
    for name in ["yeast_UPS_mini.fasta", "protein_mini_with_cazymes.faa"] {
        let records = sequence::read_fasta(&input(name)).unwrap();
        let found = defects::inspect(&records, false);
        assert!(
            found.is_clean(),
            "{name}: expected a clean proteome, found {} defects: {found:?}",
            found.count()
        );
    }
}

/// The report is RustQC's own format, so it is checked for structure and for
/// carrying the figures the stats table does not.
#[test]
fn the_report_carries_composition_and_defects() {
    let name = "yeast_UPS_mini.fasta";
    let records = sequence::read_fasta(&input(name)).unwrap();
    let stats = SequenceStats::from_records(&records);
    let found = defects::inspect(&records, false);
    let path = scratch("report.txt");
    output::write_report(name, &stats, &found, &path).unwrap();

    let text = std::fs::read_to_string(&path).unwrap();
    assert!(text.contains("## Composition"));
    assert!(text.contains("## Defects"));
    assert!(text.contains("internal_stop\t0"));
    // Methionine starts every protein, so it cannot be absent.
    let methionine = text
        .lines()
        .find(|l| l.starts_with("M\t"))
        .expect("a row for methionine");
    let count: u64 = methionine.split('\t').nth(1).unwrap().parse().unwrap();
    assert!(count > 0, "methionine should appear in a real proteome");
}

// ===================================================================
// Mass spectrometry
// ===================================================================

/// The reference figures pyteomics produced, keyed by metric name.
#[cfg(feature = "proteomics")]
fn pyteomics_reference() -> std::collections::HashMap<String, String> {
    std::fs::read_to_string(fixture("small.pyteomics.tsv"))
        .unwrap()
        .lines()
        .skip(1)
        .filter_map(|l| l.split_once('\t'))
        .map(|(k, v)| (k.to_string(), v.to_string()))
        .collect()
}

/// Every counter mzdata and pyteomics can both report, on the same file.
///
/// Two independent readers agreeing on 305213 peaks across 48 spectra is a
/// stronger statement than either agreeing with itself.
#[cfg(feature = "proteomics")]
#[test]
fn spectra_metrics_match_pyteomics() {
    use rustqc::protein::spectra;

    let metrics = spectra::analyse(&input("small.mzML")).unwrap();
    let want = pyteomics_reference();
    let integer = |key: &str| -> u64 { want[key].parse().unwrap() };
    let float = |key: &str| -> f64 { want[key].parse().unwrap() };

    assert_eq!(metrics.total_spectra(), integer("spectra"));
    assert_eq!(metrics.total_peaks(), integer("peaks"));
    assert_eq!(metrics.precursors, integer("precursors"));
    assert_eq!(
        metrics.precursors_without_charge,
        integer("precursors_without_charge"),
        "this file annotates no charge states, and that must be visible"
    );

    assert!(
        (metrics.rt_min - float("rt_min")).abs() < 1e-6,
        "rt_min was {}",
        metrics.rt_min
    );
    assert!(
        (metrics.rt_max - float("rt_max")).abs() < 1e-6,
        "rt_max was {}",
        metrics.rt_max
    );
    assert!((metrics.precursor_mz_min - float("precursor_mz_min")).abs() < 1e-4);
    assert!((metrics.precursor_mz_max - float("precursor_mz_max")).abs() < 1e-4);

    for level in [1u8, 2] {
        let m = metrics
            .levels
            .get(&level)
            .unwrap_or_else(|| panic!("no MS{level} spectra"));
        assert_eq!(m.spectra, integer(&format!("ms{level}_spectra")));
        assert_eq!(m.peaks, integer(&format!("ms{level}_peaks")));
        assert_eq!(m.min_peaks, integer(&format!("ms{level}_min_peaks")));
        assert_eq!(m.max_peaks, integer(&format!("ms{level}_max_peaks")));

        // Intensities are 32-bit in the file, so the two readers accumulate
        // them at different precision; a relative tolerance is the honest
        // comparison rather than an exact one.
        let want_tic = float(&format!("ms{level}_total_ion_current"));
        let relative = (m.total_ion_current - want_tic).abs() / want_tic;
        assert!(
            relative < 1e-6,
            "MS{level} total ion current: got {}, want {want_tic}, relative {relative}",
            m.total_ion_current
        );
    }
}

/// The derived figures the report leans on.
#[cfg(feature = "proteomics")]
#[test]
fn spectra_derived_figures_are_consistent() {
    use rustqc::protein::spectra;

    let metrics = spectra::analyse(&input("small.mzML")).unwrap();
    assert_eq!(
        metrics.ms2_per_ms1(),
        Some(34.0 / 14.0),
        "34 fragmentation scans for 14 survey scans"
    );
    assert!(metrics.rt_span() > 0.0, "the run spans some time");
    assert_eq!(
        metrics.total_peaks(),
        metrics.levels.values().map(|l| l.peaks).sum::<u64>(),
        "the total must be the sum of the levels"
    );
}

/// The report is RustQC's own format, so it is checked for structure.
#[cfg(feature = "proteomics")]
#[test]
fn the_spectra_report_carries_every_section() {
    use rustqc::protein::spectra;

    let metrics = spectra::analyse(&input("small.mzML")).unwrap();
    let path = scratch("spectra_report.txt");
    spectra::output::write_report("small.mzML", &metrics, &path).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();

    assert!(text.contains("## Run"));
    assert!(text.contains("## Levels"));
    assert!(text.contains("## Precursors"));
    assert!(text.contains("spectra\t48"));
    assert!(text.contains("without_charge\t34"));
    assert!(
        !text.contains("## Charge states"),
        "this file annotates no charges, so the section is omitted"
    );
}
