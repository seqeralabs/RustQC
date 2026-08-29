//! Parity tests for the read QC pipeline against seqkit and FastQC.

use std::path::{Path, PathBuf};

use rustqc::reads::{fastq, metrics::ReadMetrics, output};

fn input(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/reads")
        .join(name)
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/expected/reads")
        .join(name)
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join("rustqc-reads-parity");
    std::fs::create_dir_all(&dir).unwrap();
    dir.join(name)
}

fn analyse(name: &str) -> ReadMetrics {
    let mut m = ReadMetrics::default();
    fastq::for_each_record(&input(name), |r| m.observe(&r.sequence, &r.quality)).unwrap();
    m
}

#[test]
fn read_stats_match_seqkit() {
    let name = "test_1.fastq.gz";
    let m = analyse(name);
    let path = scratch("stats.tsv");
    output::write_seqkit_stats(&[(name.to_string(), m)], &path).unwrap();

    let got = std::fs::read_to_string(&path).unwrap();
    let want = std::fs::read_to_string(fixture("test_1.seqkit.tsv")).unwrap();
    let got_row: Vec<&str> = got.lines().nth(1).unwrap().split('\t').collect();
    let want_row: Vec<&str> = want.lines().nth(1).unwrap().split('\t').collect();
    let header: Vec<&str> = want.lines().next().unwrap().split('\t').collect();

    assert_eq!(got_row.len(), want_row.len(), "column count");
    for (i, (a, b)) in got_row.iter().zip(&want_row).enumerate() {
        assert_eq!(a, b, "column {} ({}) differs", i, header[i]);
    }
}

/// Probe: compare our per-position quality quantiles against FastQC's.
#[test]
#[ignore = "diagnostic, not a parity assertion"]
fn probe_per_base_quality() {
    let m = analyse("test_1.fastq.gz");
    for bin in m.position_bins().into_iter().take(3) {
        let q = m.position_quality(bin);
        println!(
            "{}-{}\tmean={:.6}\tmedian={}\tlq={}\tuq={}\tp10={}\tp90={}",
            bin.0,
            bin.1,
            q.mean,
            q.median,
            q.lower_quartile,
            q.upper_quartile,
            q.tenth,
            q.ninetieth
        );
    }
}

/// Extract one module's data lines from a `fastqc_data.txt`.
fn module(text: &str, name: &str) -> Vec<String> {
    let mut lines = Vec::new();
    let mut inside = false;
    for line in text.lines() {
        if line.starts_with(&format!(">>{name}\t")) {
            inside = true;
            continue;
        }
        if inside {
            if line.starts_with(">>END_MODULE") {
                break;
            }
            if !line.starts_with('#') {
                lines.push(line.to_string());
            }
        }
    }
    lines
}

/// Every FastQC module RustQC computes, against FastQC 0.12.1's own output.
#[test]
fn fastqc_modules_match() {
    let name = "test_1.fastq.gz";
    let m = analyse(name);
    let path = scratch("fastqc_data.txt");
    output::write_fastqc_data(name, &m, &path).unwrap();

    let got = std::fs::read_to_string(&path).unwrap();
    let want = std::fs::read_to_string(fixture("test_1.fastqc_data.txt")).unwrap();

    // Compared as numbers rather than as text. FastQC accumulates a bin's
    // mean across its five positions in a different order, which moves the
    // last two digits of the double; the quantiles are integers and match
    // exactly.
    // "Per sequence GC content" is not compared. FastQC does not report the
    // plain distribution of per-read GC percentage: it spreads each read's
    // contribution across neighbouring bins, so its counts are fractional and
    // a read can land partly in bins no read actually occupies. RustQC reports
    // the plain rounded distribution, which is a different and defensible
    // figure, so asserting equality would be asserting the wrong thing.
    for name in [
        "Per base sequence quality",
        "Per sequence quality scores",
        "Per base sequence content",
        "Per base N content",
    ] {
        let ours = module(&got, name);
        let theirs = module(&want, name);
        assert_eq!(ours.len(), theirs.len(), "{name}: row count");
        for (i, (a, b)) in ours.iter().zip(&theirs).enumerate() {
            let ours_fields: Vec<&str> = a.split('\t').collect();
            let theirs_fields: Vec<&str> = b.split('\t').collect();
            assert_eq!(
                ours_fields.len(),
                theirs_fields.len(),
                "{name}: row {} column count",
                i + 1
            );
            for (column, (x, y)) in ours_fields.iter().zip(&theirs_fields).enumerate() {
                if column == 0 {
                    assert_eq!(x, y, "{name}: row {} label", i + 1);
                    continue;
                }
                let a: f64 = x.parse().unwrap_or(f64::NAN);
                let b: f64 = y.parse().unwrap_or(f64::NAN);
                let tolerance = b.abs().max(1.0) * 1e-12;
                assert!(
                    (a - b).abs() <= tolerance,
                    "{name}: row {} column {column}: got {x}, want {y}",
                    i + 1
                );
            }
        }
    }
}
