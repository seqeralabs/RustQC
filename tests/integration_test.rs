//! Integration tests comparing RustQC output against R dupRadar reference output.
//!
//! These tests run RustQC on the same test BAM/GTF data used to generate
//! the R reference outputs in tests/expected/.

use std::fs;
use std::path::Path;
use std::path::PathBuf;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use rust_htslib::bam;
use rust_htslib::bam::header::HeaderRecord;

/// Helper: get the path to the rustqc binary.
///
/// Uses the `CARGO_BIN_EXE_rustqc` env var when available (set by cargo test),
/// otherwise falls back to building and using `target/debug/rustqc`.
fn rustqc_binary() -> String {
    // cargo test sets this for integration tests automatically
    if let Ok(path) = std::env::var("CARGO_BIN_EXE_rustqc") {
        return path;
    }

    // Fallback: build and use debug binary
    let status = Command::new("cargo")
        .args(["build"])
        .status()
        .expect("Failed to run cargo build");
    assert!(status.success(), "cargo build failed");

    let path = Path::new("target/debug/rustqc");
    assert!(path.exists(), "Binary not found at {:?}", path);
    path.to_str().unwrap().to_string()
}

/// Helper: run rustqc rna on test data and return the output directory
fn run_rustqc(outdir: &str) -> std::process::Output {
    run_rustqc_with_input(outdir, "tests/data/test.bam")
}

/// Helper: run rustqc rna with a specified input file and return the output.
/// When `skip_dup_check` is true (the default for most tests), passes
/// `--skip-dup-check` so the test doesn't depend on duplicate flags.
fn run_rustqc_impl(outdir: &str, input: &str, skip_dup_check: bool) -> std::process::Output {
    let binary = rustqc_binary();
    let mut args = vec![
        "rna",
        input,
        "--gtf",
        "tests/data/test.gtf",
        "--outdir",
        outdir,
    ];
    if skip_dup_check {
        args.push("--skip-dup-check");
    }
    Command::new(&binary)
        .args(&args)
        .output()
        .expect("Failed to execute rustqc")
}

/// Helper: run rustqc rna with custom arguments and environment.
fn run_rustqc_custom(
    outdir: &Path,
    input: &Path,
    gtf: &Path,
    threads: usize,
    skip_dup_check: bool,
    extra_env: &[(&str, &str)],
) -> std::process::Output {
    let binary = rustqc_binary();
    let mut args = vec![
        "rna".to_string(),
        input.display().to_string(),
        "--gtf".to_string(),
        gtf.display().to_string(),
        "--outdir".to_string(),
        outdir.display().to_string(),
        "--threads".to_string(),
        threads.to_string(),
    ];
    if skip_dup_check {
        args.push("--skip-dup-check".to_string());
    }

    let mut cmd = Command::new(&binary);
    cmd.args(&args);
    for (key, value) in extra_env {
        cmd.env(key, value);
    }
    cmd.output().expect("Failed to execute rustqc")
}

fn unique_test_dir(label: &str) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir =
        std::env::temp_dir().join(format!("rustqc-{}-{}-{}", label, std::process::id(), nonce));
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn write_test_gtf(path: &Path, chroms: &[&str]) {
    let mut gtf = String::new();
    for (i, chrom) in chroms.iter().enumerate() {
        let gene_id = format!("GENE_{}", i + 1);
        gtf.push_str(&format!(
            "{chrom}\ttest\tgene\t1\t100\t.\t+\t.\tgene_id \"{gene_id}\"; gene_name \"{gene_id}\";\n"
        ));
        gtf.push_str(&format!(
            "{chrom}\ttest\texon\t1\t100\t.\t+\t.\tgene_id \"{gene_id}\"; gene_name \"{gene_id}\"; exon_number \"1\";\n"
        ));
    }
    fs::write(path, gtf).unwrap();
}

fn write_bam_fixture(path: &Path, chroms: &[(&str, u64)], sam_lines: &[&str]) {
    let mut header = bam::Header::new();
    header.push_record(
        HeaderRecord::new(b"HD")
            .push_tag(b"VN", "1.6")
            .push_tag(b"SO", "coordinate"),
    );
    for (chrom, len) in chroms {
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", *chrom)
                .push_tag(b"LN", *len as i64),
        );
    }

    let header_view = bam::HeaderView::from_header(&header);
    let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam).unwrap();
    for line in sam_lines {
        let record = bam::Record::from_sam(&header_view, line.as_bytes()).unwrap();
        writer.write(&record).unwrap();
    }
    drop(writer);

    let bai_path = PathBuf::from(format!("{}.bai", path.display()));
    bam::index::build(path, Some(&bai_path), bam::index::Type::Bai, 1).unwrap();
}

fn build_late_duplicate_fixture(root: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let bam_path = root.join("late-dup.bam");
    let gtf_path = root.join("late-dup.gtf");
    let outdir = root.join("out");
    fs::create_dir_all(&outdir).unwrap();

    write_test_gtf(&gtf_path, &["chrLate"]);
    write_bam_fixture(
        &bam_path,
        &[("chrLate", 1000)],
        &[
            "late1\t0\tchrLate\t1\t60\t4M\t*\t0\t0\tACGT\tIIII",
            "late2\t0\tchrLate\t10\t60\t4M\t*\t0\t0\tTGCA\tIIII",
            "late3\t1024\tchrLate\t20\t60\t4M\t*\t0\t0\tGATC\tIIII",
        ],
    );

    (bam_path, gtf_path, outdir)
}

fn build_parallel_duplicate_fixture(root: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let bam_path = root.join("parallel-dup.bam");
    let gtf_path = root.join("parallel-dup.gtf");
    let outdir = root.join("out");
    fs::create_dir_all(&outdir).unwrap();

    write_test_gtf(&gtf_path, &["chrDup", "chrNoDup"]);
    write_bam_fixture(
        &bam_path,
        &[("chrDup", 1000), ("chrNoDup", 1000)],
        &[
            "dup1\t1024\tchrDup\t1\t60\t4M\t*\t0\t0\tACGT\tIIII",
            "nodup1\t0\tchrNoDup\t1\t60\t4M\t*\t0\t0\tTGCA\tIIII",
            "nodup2\t0\tchrNoDup\t10\t60\t4M\t*\t0\t0\tGATC\tIIII",
        ],
    );

    (bam_path, gtf_path, outdir)
}

/// Convenience: run with a custom input and --skip-dup-check.
fn run_rustqc_with_input(outdir: &str, input: &str) -> std::process::Output {
    run_rustqc_impl(outdir, input, true)
}

/// Helper: run rustqc rna with --flat-output flag (all files in outdir, no subdirectories)
fn run_rustqc_flat(outdir: &str) -> std::process::Output {
    let binary = rustqc_binary();
    Command::new(&binary)
        .args([
            "rna",
            "tests/data/test.bam",
            "--gtf",
            "tests/data/test.gtf",
            "--outdir",
            outdir,
            "--skip-dup-check",
            "--flat-output",
        ])
        .output()
        .expect("Failed to execute rustqc")
}

/// Parse a dupMatrix TSV file into a Vec of (gene_id, HashMap<col_name, value_string>)
fn parse_dup_matrix(path: &str) -> Vec<(String, Vec<String>)> {
    let content = fs::read_to_string(path).expect("Failed to read dup matrix");
    let mut lines = content.lines();
    let _header = lines.next().expect("Missing header");
    let mut rows = Vec::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        let gene_id = fields[0].clone();
        rows.push((gene_id, fields));
    }
    rows
}

/// Parse the intercept/slope file (format: "label\tvalue\n")
fn parse_intercept_slope(path: &str) -> (f64, f64) {
    let content = fs::read_to_string(path).expect("Failed to read intercept/slope");
    let mut intercept = 0.0;
    let mut slope = 0.0;
    for line in content.lines() {
        let line = line.trim();
        if let Some(val) = line.strip_prefix("intercept\t") {
            // Legacy tab-separated format
            intercept = val.trim().parse().unwrap();
        } else if let Some(val) = line.strip_prefix("slope\t") {
            slope = val.trim().parse().unwrap();
        } else if line.contains("dupRadar Int") {
            // Upstream R dupRadar format: "sample - dupRadar Int (...): <value>"
            if let Some(val) = line.rsplit(": ").next() {
                intercept = val.trim().parse().unwrap();
            }
        } else if line.contains("dupRadar Sl") {
            if let Some(val) = line.rsplit(": ").next() {
                slope = val.trim().parse().unwrap();
            }
        }
    }
    (intercept, slope)
}

/// Helper: run rustqc rna on test data with multiple BAM files
fn run_rustqc_multi(outdir: &str) -> std::process::Output {
    let binary = rustqc_binary();

    // Ensure output directory exists before copying files into it
    fs::create_dir_all(outdir).expect("Failed to create output directory");

    // Create a copy of test.bam with a different stem so outputs don't collide
    let copy_bam = format!("{}/test_copy.bam", outdir);
    let copy_bai = format!("{}/test_copy.bam.bai", outdir);
    fs::copy("tests/data/test.bam", &copy_bam).expect("Failed to copy BAM");
    fs::copy("tests/data/test.bam.bai", &copy_bai).expect("Failed to copy BAI");

    Command::new(&binary)
        .args([
            "rna",
            "tests/data/test.bam",
            &copy_bam,
            "--gtf",
            "tests/data/test.gtf",
            "--outdir",
            outdir,
            "--skip-dup-check",
        ])
        .output()
        .expect("Failed to execute rustqc")
}

/// Compare two floating-point values with a relative tolerance
fn approx_eq(a: f64, b: f64, rel_tol: f64) -> bool {
    if a.is_nan() && b.is_nan() {
        return true;
    }
    if a == b {
        return true;
    }
    let max_abs = a.abs().max(b.abs());
    if max_abs == 0.0 {
        return (a - b).abs() < 1e-15;
    }
    (a - b).abs() / max_abs < rel_tol
}

#[test]
fn test_dup_matrix_exact_match() {
    let outdir = "tests/output_integration";
    let _ = fs::remove_dir_all(outdir);
    fs::create_dir_all(outdir).unwrap();

    let output = run_rustqc(outdir);
    assert!(
        output.status.success(),
        "rustqc failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Compare dup matrix files
    let r_matrix = parse_dup_matrix("tests/expected/dupMatrix.txt");
    let rust_matrix = parse_dup_matrix(&format!("{}/dupradar/test_dupMatrix.txt", outdir));

    assert_eq!(
        r_matrix.len(),
        rust_matrix.len(),
        "Different number of genes: R={}, Rust={}",
        r_matrix.len(),
        rust_matrix.len()
    );

    for (r_row, rust_row) in r_matrix.iter().zip(rust_matrix.iter()) {
        assert_eq!(
            r_row.0, rust_row.0,
            "Gene order mismatch: R={}, Rust={}",
            r_row.0, rust_row.0
        );

        assert_eq!(
            r_row.1.len(),
            rust_row.1.len(),
            "Different number of columns for gene {}",
            r_row.0
        );

        for (col_idx, (r_val, rust_val)) in r_row.1.iter().zip(rust_row.1.iter()).enumerate() {
            // Gene ID and gene length are exact
            if col_idx <= 1 {
                assert_eq!(
                    r_val, rust_val,
                    "Column {} mismatch for gene {}: R='{}', Rust='{}'",
                    col_idx, r_row.0, r_val, rust_val
                );
                continue;
            }

            // For numeric columns, compare as floats with tolerance
            if r_val == "NA" || rust_val == "NA" {
                assert_eq!(
                    r_val, rust_val,
                    "NA mismatch at col {} for gene {}: R='{}', Rust='{}'",
                    col_idx, r_row.0, r_val, rust_val
                );
                continue;
            }

            let r_f: f64 = r_val.parse().unwrap_or_else(|_| {
                panic!(
                    "Failed to parse R value '{}' at col {} for gene {}",
                    r_val, col_idx, r_row.0
                )
            });
            let rust_f: f64 = rust_val.parse().unwrap_or_else(|_| {
                panic!(
                    "Failed to parse Rust value '{}' at col {} for gene {}",
                    rust_val, col_idx, r_row.0
                )
            });

            assert!(
                approx_eq(r_f, rust_f, 1e-6),
                "Value mismatch at col {} for gene {}: R={}, Rust={} (diff={})",
                col_idx,
                r_row.0,
                r_f,
                rust_f,
                (r_f - rust_f).abs()
            );
        }
    }

    // Cleanup
    let _ = fs::remove_dir_all(outdir);
}

#[test]
fn test_intercept_slope_match() {
    let outdir = "tests/output_integration_fit";
    let _ = fs::remove_dir_all(outdir);
    fs::create_dir_all(outdir).unwrap();

    let output = run_rustqc(outdir);
    assert!(
        output.status.success(),
        "rustqc failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let (r_intercept, r_slope) = parse_intercept_slope("tests/expected/intercept_slope.txt");
    let (rust_intercept, rust_slope) =
        parse_intercept_slope(&format!("{}/dupradar/test_intercept_slope.txt", outdir));

    // Allow slightly more tolerance for the IRLS implementation differences
    assert!(
        approx_eq(r_intercept, rust_intercept, 1e-6),
        "Intercept mismatch: R={}, Rust={} (diff={})",
        r_intercept,
        rust_intercept,
        (r_intercept - rust_intercept).abs()
    );
    assert!(
        approx_eq(r_slope, rust_slope, 1e-6),
        "Slope mismatch: R={}, Rust={} (diff={})",
        r_slope,
        rust_slope,
        (r_slope - rust_slope).abs()
    );

    // Cleanup
    let _ = fs::remove_dir_all(outdir);
}

#[test]
fn test_all_output_files_generated() {
    let outdir = "tests/output_integration_files";
    let _ = fs::remove_dir_all(outdir);
    fs::create_dir_all(outdir).unwrap();

    let output = run_rustqc(outdir);
    assert!(
        output.status.success(),
        "rustqc failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Check all expected output files exist
    let expected_files = vec![
        "dupradar/test_dupMatrix.txt",
        "dupradar/test_intercept_slope.txt",
        "dupradar/test_duprateExpDens.png",
        "dupradar/test_duprateExpDens.svg",
        "dupradar/test_duprateExpBoxplot.png",
        "dupradar/test_duprateExpBoxplot.svg",
        "dupradar/test_expressionHist.png",
        "dupradar/test_expressionHist.svg",
        "dupradar/test_dup_intercept_mqc.txt",
        "dupradar/test_duprateExpDensCurve_mqc.txt",
        "featurecounts/test.featureCounts.tsv",
        "featurecounts/test.featureCounts.tsv.summary",
        "rseqc/bam_stat/test.bam_stat.txt",
        "rseqc/infer_experiment/test.infer_experiment.txt",
        "rseqc/read_duplication/test.pos.DupRate.xls",
        "rseqc/read_duplication/test.seq.DupRate.xls",
        "rseqc/read_distribution/test.read_distribution.txt",
        "rseqc/junction_annotation/test.junction_annotation.log",
        "rseqc/junction_saturation/test.junctionSaturation_summary.txt",
        "rseqc/inner_distance/test.inner_distance.txt",
    ];

    // Files that may be empty with a small test dataset (e.g. no proper pairs
    // for inner distance calculation)
    let allow_empty: Vec<&str> = vec!["rseqc/inner_distance/test.inner_distance.txt"];

    for file in &expected_files {
        let path = format!("{}/{}", outdir, file);
        assert!(Path::new(&path).exists(), "Missing output file: {}", file);

        // Check file is non-empty (unless it is expected to be empty for the small test data)
        if !allow_empty.contains(file) {
            let metadata = fs::metadata(&path).unwrap();
            assert!(metadata.len() > 0, "Output file is empty: {}", file);
        }
    }

    // Cleanup
    let _ = fs::remove_dir_all(outdir);
}

#[test]
fn test_dup_matrix_header() {
    let outdir = "tests/output_integration_header";
    let _ = fs::remove_dir_all(outdir);
    fs::create_dir_all(outdir).unwrap();

    let output = run_rustqc(outdir);
    assert!(output.status.success());

    let content = fs::read_to_string(format!("{}/dupradar/test_dupMatrix.txt", outdir)).unwrap();
    let header = content.lines().next().unwrap();
    let expected_cols = vec![
        "ID",
        "geneLength",
        "allCountsMulti",
        "filteredCountsMulti",
        "dupRateMulti",
        "dupsPerIdMulti",
        "RPKMulti",
        "PKMMulti",
        "allCounts",
        "filteredCounts",
        "dupRate",
        "dupsPerId",
        "RPK",
        "RPKM",
    ];

    let actual_cols: Vec<&str> = header.split('\t').collect();
    assert_eq!(
        actual_cols.len(),
        expected_cols.len(),
        "Wrong number of columns: expected {}, got {}",
        expected_cols.len(),
        actual_cols.len()
    );

    for (i, (expected, actual)) in expected_cols.iter().zip(actual_cols.iter()).enumerate() {
        assert_eq!(
            expected, actual,
            "Column {} name mismatch: expected '{}', got '{}'",
            i, expected, actual
        );
    }

    // Cleanup
    let _ = fs::remove_dir_all(outdir);
}

#[test]
fn test_mqc_intercept_format() {
    let outdir = "tests/output_integration_mqc";
    let _ = fs::remove_dir_all(outdir);
    fs::create_dir_all(outdir).unwrap();

    let output = run_rustqc(outdir);
    assert!(output.status.success());

    // Check MultiQC intercept file format
    let content =
        fs::read_to_string(format!("{}/dupradar/test_dup_intercept_mqc.txt", outdir)).unwrap();
    // Skip YAML comment lines (starting with #)
    let data_lines: Vec<&str> = content.lines().filter(|l| !l.starts_with('#')).collect();
    assert!(data_lines.len() >= 2, "MultiQC intercept file too short");

    // First data line should be a header with "Sample" and "dupRadar_intercept"
    assert!(
        data_lines[0].contains("Sample"),
        "Missing 'Sample' in MQC header"
    );
    assert!(
        data_lines[0].contains("dupRadar_intercept"),
        "Missing 'dupRadar_intercept' in MQC header"
    );

    // Second data line should have sample name and numeric value
    let parts: Vec<&str> = data_lines[1].split('\t').collect();
    assert_eq!(parts.len(), 2, "MQC data line should have 2 columns");
    let _intercept: f64 = parts[1]
        .trim()
        .parse()
        .expect("MQC intercept should be numeric");

    // Cleanup
    let _ = fs::remove_dir_all(outdir);
}

#[test]
fn test_mqc_curve_format() {
    let outdir = "tests/output_integration_curve";
    let _ = fs::remove_dir_all(outdir);
    fs::create_dir_all(outdir).unwrap();

    let output = run_rustqc(outdir);
    assert!(output.status.success());

    // Check MultiQC curve file format
    let content = fs::read_to_string(format!(
        "{}/dupradar/test_duprateExpDensCurve_mqc.txt",
        outdir
    ))
    .unwrap();
    // Skip YAML comment lines (starting with #)
    let data_lines: Vec<&str> = content.lines().filter(|l| !l.starts_with('#')).collect();
    // Header + at least some data points (101 evenly spaced + header = 102)
    assert!(
        data_lines.len() >= 3,
        "MultiQC curve file too short: {} lines",
        data_lines.len()
    );

    // Should have a header line with 2 columns
    let header_parts: Vec<&str> = data_lines[0].split('\t').collect();
    assert_eq!(header_parts.len(), 2, "Curve header should have 2 columns");

    // Data lines should have numeric values
    for line in &data_lines[1..] {
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        assert_eq!(parts.len(), 2, "Curve data should have 2 columns");
        let _x: f64 = parts[0].parse().expect("X should be numeric");
        let _y: f64 = parts[1].parse().expect("Y should be numeric");
    }

    // Cleanup
    let _ = fs::remove_dir_all(outdir);
}

#[test]
fn test_gene_order_preserved() {
    let outdir = "tests/output_integration_order";
    let _ = fs::remove_dir_all(outdir);
    fs::create_dir_all(outdir).unwrap();

    let output = run_rustqc(outdir);
    assert!(output.status.success());

    // Both R and Rust should output genes in the same order (GTF order)
    let r_matrix = parse_dup_matrix("tests/expected/dupMatrix.txt");
    let rust_matrix = parse_dup_matrix(&format!("{}/dupradar/test_dupMatrix.txt", outdir));

    let r_genes: Vec<&str> = r_matrix.iter().map(|(id, _)| id.as_str()).collect();
    let rust_genes: Vec<&str> = rust_matrix.iter().map(|(id, _)| id.as_str()).collect();

    assert_eq!(
        r_genes, rust_genes,
        "Gene order differs:\nR: {:?}\nRust: {:?}",
        r_genes, rust_genes
    );

    // Cleanup
    let _ = fs::remove_dir_all(outdir);
}

#[test]
fn test_count_values_exact() {
    let outdir = "tests/output_integration_counts";
    let _ = fs::remove_dir_all(outdir);
    fs::create_dir_all(outdir).unwrap();

    let output = run_rustqc(outdir);
    assert!(output.status.success());

    let r_matrix = parse_dup_matrix("tests/expected/dupMatrix.txt");
    let rust_matrix = parse_dup_matrix(&format!("{}/dupradar/test_dupMatrix.txt", outdir));

    // Integer count columns should match EXACTLY (columns 2-3, 5, 8-9, 11)
    // col indices: 2=allCountsMulti, 3=filteredCountsMulti, 5=dupsPerIdMulti,
    //              8=allCounts, 9=filteredCounts, 11=dupsPerId
    let count_cols = [2, 3, 5, 8, 9, 11];

    for (r_row, rust_row) in r_matrix.iter().zip(rust_matrix.iter()) {
        for &col in &count_cols {
            assert_eq!(
                r_row.1[col], rust_row.1[col],
                "Count mismatch at col {} for gene {}: R='{}', Rust='{}'",
                col, r_row.0, r_row.1[col], rust_row.1[col]
            );
        }
    }

    // Cleanup
    let _ = fs::remove_dir_all(outdir);
}

// ===================================================================
// Multi-BAM tests
// ===================================================================

/// Multiple BAM files should each produce their own output files.
#[test]
fn test_multiple_bam_files() {
    let outdir = "tests/output_multi_bam";
    let output = run_rustqc_multi(outdir);

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "rustqc failed with multiple BAMs:\n{}",
        stderr
    );

    // Both BAM stems should have produced output files
    let expected_files = [
        "dupradar/test_dupMatrix.txt",
        "dupradar/test_duprateExpDens.png",
        "dupradar/test_duprateExpBoxplot.png",
        "dupradar/test_expressionHist.png",
        "dupradar/test_copy_dupMatrix.txt",
        "dupradar/test_copy_duprateExpDens.png",
        "dupradar/test_copy_duprateExpBoxplot.png",
        "dupradar/test_copy_expressionHist.png",
    ];

    for filename in &expected_files {
        let path = format!("{}/{}", outdir, filename);
        assert!(
            std::path::Path::new(&path).exists(),
            "Expected output file missing: {}",
            path
        );
    }

    // Both dupMatrix files should match the R reference (same input data)
    let r_matrix = parse_dup_matrix("tests/expected/dupMatrix.txt");
    for stem in &["test", "test_copy"] {
        let rust_matrix = parse_dup_matrix(&format!("{}/dupradar/{}_dupMatrix.txt", outdir, stem));
        assert_eq!(
            r_matrix.len(),
            rust_matrix.len(),
            "Row count mismatch for {} matrix",
            stem
        );
    }

    // Cleanup
    let _ = fs::remove_dir_all(outdir);
}

/// Passing the same BAM file twice should fail with a duplicate-stem error.
#[test]
fn test_duplicate_bam_stems_rejected() {
    let bin = rustqc_binary();
    let output = Command::new(&bin)
        .args([
            "rna",
            "tests/data/test.bam",
            "tests/data/test.bam",
            "--gtf",
            "tests/data/test.gtf",
            "--outdir",
            "tests/output_dup_stems",
        ])
        .output()
        .expect("Failed to execute rustqc");

    assert!(
        !output.status.success(),
        "rustqc should fail with duplicate BAM stems"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("duplicate") || stderr.contains("Duplicate"),
        "Error message should mention duplicate stems, got:\n{}",
        stderr
    );

    // Cleanup (may not exist if it failed early)
    let _ = fs::remove_dir_all("tests/output_dup_stems");
}

#[test]
fn test_sam_input_matches_bam() {
    // SAM input should produce identical results to BAM input.
    // SAM files have no index, so this exercises the single-threaded code path.
    let outdir_bam = "tests/output_integration_sam_bam";
    let outdir_sam = "tests/output_integration_sam_sam";
    let _ = fs::remove_dir_all(outdir_bam);
    let _ = fs::remove_dir_all(outdir_sam);
    fs::create_dir_all(outdir_bam).unwrap();
    fs::create_dir_all(outdir_sam).unwrap();

    let output_bam = run_rustqc_with_input(outdir_bam, "tests/data/test.bam");
    assert!(
        output_bam.status.success(),
        "rustqc with BAM input failed: {}",
        String::from_utf8_lossy(&output_bam.stderr)
    );

    let output_sam = run_rustqc_with_input(outdir_sam, "tests/data/test.sam");
    assert!(
        output_sam.status.success(),
        "rustqc with SAM input failed: {}",
        String::from_utf8_lossy(&output_sam.stderr)
    );

    // Compare dup matrices — SAM and BAM should produce identical output
    let bam_matrix = parse_dup_matrix(&format!("{}/dupradar/test_dupMatrix.txt", outdir_bam));
    let sam_matrix = parse_dup_matrix(&format!("{}/dupradar/test_dupMatrix.txt", outdir_sam));

    assert_eq!(
        bam_matrix.len(),
        sam_matrix.len(),
        "Different number of genes: BAM={}, SAM={}",
        bam_matrix.len(),
        sam_matrix.len()
    );

    for (bam_row, sam_row) in bam_matrix.iter().zip(sam_matrix.iter()) {
        assert_eq!(
            bam_row.0, sam_row.0,
            "Gene order mismatch: BAM={}, SAM={}",
            bam_row.0, sam_row.0
        );

        for (col_idx, (bam_val, sam_val)) in bam_row.1.iter().zip(sam_row.1.iter()).enumerate() {
            assert_eq!(
                bam_val, sam_val,
                "Value mismatch at col {} for gene {}: BAM='{}', SAM='{}'",
                col_idx, bam_row.0, bam_val, sam_val
            );
        }
    }

    // Cleanup
    let _ = fs::remove_dir_all(outdir_bam);
    let _ = fs::remove_dir_all(outdir_sam);
}

// ===================================================================
// Flat output mode tests
// ===================================================================

/// With --flat-output, all files should be written directly in the output directory
/// with no subdirectories.
#[test]
fn test_flat_output_no_subdirectories() {
    let outdir = "tests/output_integration_flat";
    let _ = fs::remove_dir_all(outdir);
    fs::create_dir_all(outdir).unwrap();

    let output = run_rustqc_flat(outdir);
    assert!(
        output.status.success(),
        "rustqc --flat-output failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // These files should exist directly in outdir (no subdirectories)
    let expected_flat_files = vec![
        "test_dupMatrix.txt",
        "test_intercept_slope.txt",
        "test_duprateExpDens.png",
        "test.featureCounts.tsv",
        "test.featureCounts.tsv.summary",
        "test.bam_stat.txt",
        "test.infer_experiment.txt",
        "test.pos.DupRate.xls",
        "test.seq.DupRate.xls",
        "test.read_distribution.txt",
        "test.junction_annotation.log",
        "test.junctionSaturation_summary.txt",
        "test.inner_distance.txt",
    ];

    for file in &expected_flat_files {
        let path = format!("{}/{}", outdir, file);
        assert!(
            Path::new(&path).exists(),
            "Missing flat output file: {}",
            file
        );
    }

    // Subdirectories should NOT exist in flat mode
    let nested_dirs = vec!["dupradar", "featurecounts", "rseqc", "samtools"];
    for dir in &nested_dirs {
        let path = format!("{}/{}", outdir, dir);
        assert!(
            !Path::new(&path).is_dir(),
            "Subdirectory '{}' should not exist in flat output mode",
            dir
        );
    }

    // Cleanup
    let _ = fs::remove_dir_all(outdir);
}

// ===================================================================
// Duplicate-flag detection tests (https://github.com/seqeralabs/RustQC/issues/71)
//
// PR #84 replaced the fragile @PG header check with runtime inspection of
// the actual 0x400 duplicate flag. These tests exercise that logic.
// ===================================================================

/// Convenience: run WITHOUT --skip-dup-check (the default user experience).
fn run_rustqc_dup_check(outdir: &str, input: &str) -> std::process::Output {
    run_rustqc_impl(outdir, input, false)
}

/// A BAM with duplicates marked (0x400) should succeed without --skip-dup-check.
/// This is the normal happy path — tools like Picard, samblaster, or Parabricks
/// set the duplicate flag, and RustQC should accept the file.
#[test]
fn test_dup_check_passes_with_marked_duplicates() {
    let outdir = "tests/output_dup_check_pass";
    let _ = fs::remove_dir_all(outdir);
    fs::create_dir_all(outdir).unwrap();

    // test.bam has 137/488 reads with 0x400 set
    let output = run_rustqc_dup_check(outdir, "tests/data/test.bam");
    assert!(
        output.status.success(),
        "rustqc should succeed when BAM has duplicate-flagged reads:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Cleanup
    let _ = fs::remove_dir_all(outdir);
}

/// A BAM with NO duplicates flagged should fail with a clear error message.
/// test_nodup.bam is test.bam with all 0x400 flags cleared.
#[test]
fn test_dup_check_fails_without_marked_duplicates() {
    let outdir = "tests/output_dup_check_fail";
    let _ = fs::remove_dir_all(outdir);
    fs::create_dir_all(outdir).unwrap();

    let output = run_rustqc_dup_check(outdir, "tests/data/test_nodup.bam");
    assert!(
        !output.status.success(),
        "rustqc should fail when BAM has no duplicate-flagged reads"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("No duplicate-flagged reads found"),
        "Error should mention missing duplicate flags, got:\n{}",
        stderr
    );
    assert!(
        stderr.contains("--skip-dup-check"),
        "Error should suggest --skip-dup-check, got:\n{}",
        stderr
    );

    // Cleanup
    let _ = fs::remove_dir_all(outdir);
}

/// --skip-dup-check should bypass the duplicate flag requirement,
/// allowing a no-duplicate BAM to be processed successfully.
#[test]
fn test_skip_dup_check_bypasses_failure() {
    let outdir = "tests/output_dup_check_skip";
    let _ = fs::remove_dir_all(outdir);
    fs::create_dir_all(outdir).unwrap();

    let output = run_rustqc_with_input(outdir, "tests/data/test_nodup.bam");
    assert!(
        output.status.success(),
        "rustqc --skip-dup-check should succeed even without duplicate flags:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Cleanup
    let _ = fs::remove_dir_all(outdir);
}

/// SAM input (single-threaded path) should also detect missing duplicates.
/// This exercises the non-parallel code path in count_reads().
#[test]
fn test_dup_check_fails_for_sam_without_duplicates() {
    let outdir = "tests/output_dup_check_sam";
    let _ = fs::remove_dir_all(outdir);
    fs::create_dir_all(outdir).unwrap();

    let output = run_rustqc_dup_check(outdir, "tests/data/test_nodup.sam");
    assert!(
        !output.status.success(),
        "rustqc should fail for SAM input without duplicate flags"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("No duplicate-flagged reads found"),
        "Error should mention missing duplicate flags for SAM input, got:\n{}",
        stderr
    );

    // Cleanup
    let _ = fs::remove_dir_all(outdir);
}

/// A valid BAM should still pass if duplicate-flagged reads appear later than
/// the early-check threshold. The full file contains duplicates, so the
/// validation must consider the complete input rather than a prefix.
#[test]
fn test_dup_check_allows_late_duplicate_flags() {
    let root = unique_test_dir("late-dup");
    let (bam_path, gtf_path, outdir) = build_late_duplicate_fixture(&root);

    let output = run_rustqc_custom(
        &outdir,
        &bam_path,
        &gtf_path,
        1,
        false,
        &[("RUSTQC_DUP_CHECK_THRESHOLD", "2")],
    );

    assert!(
        output.status.success(),
        "rustqc should succeed when duplicate-flagged reads appear later in the file:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let _ = fs::remove_dir_all(root);
}

/// Parallel duplicate validation must use the global duplicate count, not a
/// per-worker count. One chromosome batch here has duplicate flags and the
/// other does not; the BAM should still be accepted.
#[test]
fn test_dup_check_parallel_uses_global_duplicate_state() {
    let root = unique_test_dir("parallel-dup");
    let (bam_path, gtf_path, outdir) = build_parallel_duplicate_fixture(&root);

    let output = run_rustqc_custom(
        &outdir,
        &bam_path,
        &gtf_path,
        2,
        false,
        &[("RUSTQC_DUP_CHECK_THRESHOLD", "2")],
    );

    assert!(
        output.status.success(),
        "rustqc should succeed in parallel mode when duplicates exist globally:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let _ = fs::remove_dir_all(root);
}

// ============================================================================
// bigWig coverage track tests (bedtools genomecov + UCSC bedGraphToBigWig)
// ============================================================================

/// Parse a 4-column bedGraph file into tuples.
fn read_bedgraph(path: &Path) -> Vec<(String, u32, u32, u32)> {
    fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read bedGraph {}: {e}", path.display()))
        .lines()
        .filter(|l| !l.is_empty())
        .map(|line| {
            let parts: Vec<&str> = line.split('\t').collect();
            (
                parts[0].to_string(),
                parts[1].parse().unwrap(),
                parts[2].parse().unwrap(),
                parts[3].parse::<f64>().unwrap() as u32,
            )
        })
        .collect()
}

fn which_bedtools() -> Option<String> {
    for candidate in ["/tmp/bedtools2/bin/bedtools", "bedtools"] {
        if Command::new(candidate)
            .arg("--version")
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .status()
            .map(|s| s.success())
            .unwrap_or(false)
        {
            return Some(candidate.to_string());
        }
    }
    None
}

fn which_ucsc_tool(name: &str) -> Option<String> {
    [format!("/tmp/{name}"), name.to_string()]
        .into_iter()
        .find(|candidate| Path::new(candidate).exists())
}

fn bigwig_to_bedgraph(bw: &Path, out: &Path) -> bool {
    let bgtobg = match which_ucsc_tool("bigWigToBedGraph") {
        Some(p) => p,
        None => return false,
    };
    Command::new(&bgtobg)
        .args([bw.to_str().unwrap(), out.to_str().unwrap()])
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

fn write_test_chrom_sizes(path: &Path) {
    use rust_htslib::bam::Read;
    let r = rust_htslib::bam::Reader::from_path("tests/data/test.bam").unwrap();
    let h = r.header();
    let mut sizes = String::new();
    for tid in 0..h.target_count() {
        let name = String::from_utf8_lossy(h.tid2name(tid as u32));
        let len = h.target_len(tid as u32).unwrap_or(0);
        sizes.push_str(&format!("{name}\t{len}\n"));
    }
    fs::write(path, sizes).unwrap();
}

fn sort_bedgraph_file(path: &Path) {
    let sorted = fs::read_to_string(path).unwrap();
    let mut lines: Vec<String> = sorted.lines().map(String::from).collect();
    lines.sort_by(|a, b| {
        let pa: Vec<&str> = a.split('\t').collect();
        let pb: Vec<&str> = b.split('\t').collect();
        pa[0].cmp(pb[0]).then_with(|| {
            pa[1]
                .parse::<u64>()
                .unwrap()
                .cmp(&pb[1].parse::<u64>().unwrap())
        })
    });
    fs::write(path, lines.join("\n") + "\n").unwrap();
}

#[test]
fn test_bigwig_combined_matches_bedtools_genomecov() {
    if which_bedtools().is_none() || which_ucsc_tool("bigWigToBedGraph").is_none() {
        eprintln!("Skipping bigWig cross-check: upstream tools not found");
        return;
    }

    let root = unique_test_dir("bigwig");
    fs::create_dir_all(&root).unwrap();
    let chrom_sizes_path = root.join("chrom.sizes");
    write_test_chrom_sizes(&chrom_sizes_path);

    let bedtools = which_bedtools().unwrap();
    let bg = root.join("combined.bedGraph");
    let clip = root.join("combined.clip.bedGraph");
    let ref_bw = root.join("ref.bigWig");

    let status = Command::new(&bedtools)
        .args(["genomecov", "-ibam", "tests/data/test.bam", "-split", "-bg"])
        .stdout(std::process::Stdio::from(
            std::fs::File::create(&bg).expect("create bedGraph"),
        ))
        .status()
        .expect("bedtools");
    assert!(status.success());
    sort_bedgraph_file(&bg);

    Command::new(which_ucsc_tool("bedClip").unwrap())
        .args([
            bg.to_str().unwrap(),
            chrom_sizes_path.to_str().unwrap(),
            clip.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    Command::new(which_ucsc_tool("bedGraphToBigWig").unwrap())
        .args([
            clip.to_str().unwrap(),
            chrom_sizes_path.to_str().unwrap(),
            ref_bw.to_str().unwrap(),
        ])
        .status()
        .unwrap();

    let outdir = root.join("rustqc");
    let output = run_rustqc_impl(outdir.to_str().unwrap(), "tests/data/test.bam", true);
    assert!(
        output.status.success(),
        "rustqc failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let rust_bw = outdir.join("bigwig/test.bigWig");
    assert!(rust_bw.exists(), "RustQC bigWig missing");

    let ref_bg = root.join("ref.bg");
    let rust_bg = root.join("rust.bg");
    assert!(bigwig_to_bedgraph(&ref_bw, &ref_bg));
    assert!(bigwig_to_bedgraph(&rust_bw, &rust_bg));
    assert_eq!(read_bedgraph(&ref_bg), read_bedgraph(&rust_bg));

    let _ = fs::remove_dir_all(root);
}

#[test]
fn test_bigwig_stranded_forward_matches_bedtools() {
    if which_bedtools().is_none() || which_ucsc_tool("bigWigToBedGraph").is_none() {
        eprintln!("Skipping stranded bigWig cross-check: upstream tools not found");
        return;
    }

    let root = unique_test_dir("bigwig-stranded");
    fs::create_dir_all(&root).unwrap();
    let chrom_sizes_path = root.join("chrom.sizes");
    write_test_chrom_sizes(&chrom_sizes_path);

    let bedtools = which_bedtools().unwrap();
    let bg = root.join("forward.bedGraph");
    let clip = root.join("forward.clip.bedGraph");
    let ref_bw = root.join("ref.forward.bigWig");

    let status = Command::new(&bedtools)
        .args([
            "genomecov",
            "-ibam",
            "tests/data/test.bam",
            "-split",
            "-du",
            "-strand",
            "+",
            "-bg",
        ])
        .stdout(std::process::Stdio::from(
            std::fs::File::create(&bg).expect("create bedGraph"),
        ))
        .status()
        .expect("bedtools");
    assert!(status.success());
    sort_bedgraph_file(&bg);

    Command::new(which_ucsc_tool("bedClip").unwrap())
        .args([
            bg.to_str().unwrap(),
            chrom_sizes_path.to_str().unwrap(),
            clip.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    Command::new(which_ucsc_tool("bedGraphToBigWig").unwrap())
        .args([
            clip.to_str().unwrap(),
            chrom_sizes_path.to_str().unwrap(),
            ref_bw.to_str().unwrap(),
        ])
        .status()
        .unwrap();

    let outdir = root.join("rustqc");
    let binary = rustqc_binary();
    let output = Command::new(&binary)
        .args([
            "rna",
            "tests/data/test.bam",
            "--gtf",
            "tests/data/test.gtf",
            "--outdir",
            outdir.to_str().unwrap(),
            "--stranded",
            "forward",
            "--skip-dup-check",
        ])
        .output()
        .expect("rustqc");
    assert!(
        output.status.success(),
        "rustqc failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let rust_bw = outdir.join("bigwig/test.forward.bigWig");
    assert!(rust_bw.exists());

    let ref_bg = root.join("ref.bg");
    let rust_bg = root.join("rust.bg");
    assert!(bigwig_to_bedgraph(&ref_bw, &ref_bg));
    assert!(bigwig_to_bedgraph(&rust_bw, &rust_bg));
    assert_eq!(read_bedgraph(&ref_bg), read_bedgraph(&rust_bg));

    let _ = fs::remove_dir_all(root);
}

/// Build a synthetic BAM exercising every bigWig strand/CIGAR/flag edge case,
/// and confirm the combined and per-strand tracks match `bedtools genomecov`
/// exactly. The reads deliberately cover:
///
/// - read-1 forward/reverse (never strand-flipped under `-du`)
/// - read-2 forward/reverse with mate **mapped** (flipped under `-du`)
/// - read-2 forward/reverse with mate **unmapped** (NOT flipped — the subtle
///   case: bedtools only flips when the mate is mapped)
/// - read-1 with an unmapped mate, and single-end reads (no flip)
/// - spliced (`N`) and deletion (`D`) CIGARs (`-split` block boundaries)
/// - a read overhanging the chromosome end (clamping)
/// - secondary, duplicate, QC-fail and supplementary alignments (all counted)
/// - a fully unmapped read (ignored)
/// - a second chromosome (exercises the parallel per-chromosome merge)
///
/// Run at one and four threads to cover the sequential and parallel paths.
#[test]
fn test_bigwig_strand_edge_cases_match_bedtools() {
    if which_bedtools().is_none()
        || which_ucsc_tool("bigWigToBedGraph").is_none()
        || which_ucsc_tool("bedClip").is_none()
        || which_ucsc_tool("bedGraphToBigWig").is_none()
    {
        eprintln!("Skipping bigWig strand edge-case cross-check: upstream tools not found");
        return;
    }

    let root = unique_test_dir("bigwig-strand-edge");
    let bam = root.join("strand.bam");
    let gtf = root.join("strand.gtf");
    let chrom_sizes = root.join("chrom.sizes");

    // SEQ/QUAL are 10 bases for every record (query length is 10 for 10M,
    // 5M5N5M and 5M2D5M alike), so a single constant pair suffices.
    write_bam_fixture(
        &bam,
        &[("chrTest", 300), ("chrTest2", 200)],
        &[
            // chrTest — coordinate-sorted
            "r1_fwd_mm\t65\tchrTest\t11\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            "r1_rev_mm\t81\tchrTest\t31\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            "r2_fwd_mm\t129\tchrTest\t51\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            "r2_rev_mm\t145\tchrTest\t71\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            "r2_fwd_mu\t137\tchrTest\t91\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            "r2_rev_mu\t153\tchrTest\t111\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            "r1_fwd_mu\t73\tchrTest\t131\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            "se_fwd\t0\tchrTest\t151\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            "secondary_fwd\t256\tchrTest\t151\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            "se_rev\t16\tchrTest\t171\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            "spliced_fwd\t65\tchrTest\t191\t60\t5M5N5M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            "del_fwd\t0\tchrTest\t211\t60\t5M2D5M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            "dup_fwd\t1024\tchrTest\t231\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            "qcfail_fwd\t512\tchrTest\t241\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            "supp_fwd\t2048\tchrTest\t251\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            "clipend_fwd\t0\tchrTest\t295\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            // chrTest2
            "t2_fwd\t0\tchrTest2\t11\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            "t2_rev\t16\tchrTest2\t51\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            // unmapped read (sorted last) — must be ignored by both tools
            "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
        ],
    );
    write_test_gtf(&gtf, &["chrTest", "chrTest2"]);
    fs::write(&chrom_sizes, "chrTest\t300\nchrTest2\t200\n").unwrap();

    let bedtools = which_bedtools().unwrap();
    let binary = rustqc_binary();

    // (track filename suffix, bedtools genomecov strand args)
    let tracks: [(&str, &[&str]); 3] = [
        ("strand.bigWig", &["-split"]),
        ("strand.forward.bigWig", &["-split", "-du", "-strand", "+"]),
        ("strand.reverse.bigWig", &["-split", "-du", "-strand", "-"]),
    ];

    for threads in ["1", "4"] {
        let outdir = root.join(format!("out_t{threads}"));
        let output = Command::new(&binary)
            .args([
                "rna",
                bam.to_str().unwrap(),
                "--gtf",
                gtf.to_str().unwrap(),
                "--outdir",
                outdir.to_str().unwrap(),
                "--stranded",
                "forward",
                "--threads",
                threads,
                "--skip-dup-check",
            ])
            .output()
            .expect("rustqc");
        assert!(
            output.status.success(),
            "rustqc failed (threads={threads}):\n{}",
            String::from_utf8_lossy(&output.stderr)
        );

        for (suffix, strand_args) in &tracks {
            let label = format!("threads={threads} track={suffix}");

            // Reference: bedtools genomecov -> bedClip -> bedGraphToBigWig.
            let raw = root.join(format!("ref.{threads}.{suffix}.bg"));
            let mut args: Vec<&str> = vec!["genomecov", "-ibam", bam.to_str().unwrap()];
            args.extend_from_slice(strand_args);
            args.push("-bg");
            let status = Command::new(&bedtools)
                .args(&args)
                .stdout(std::process::Stdio::from(
                    std::fs::File::create(&raw).expect("create bedGraph"),
                ))
                .status()
                .expect("bedtools");
            assert!(status.success(), "bedtools failed for {label}");
            sort_bedgraph_file(&raw);

            let clip = root.join(format!("ref.{threads}.{suffix}.clip.bg"));
            Command::new(which_ucsc_tool("bedClip").unwrap())
                .args([
                    raw.to_str().unwrap(),
                    chrom_sizes.to_str().unwrap(),
                    clip.to_str().unwrap(),
                ])
                .status()
                .unwrap();
            let ref_bw = root.join(format!("ref.{threads}.{suffix}.bw"));
            Command::new(which_ucsc_tool("bedGraphToBigWig").unwrap())
                .args([
                    clip.to_str().unwrap(),
                    chrom_sizes.to_str().unwrap(),
                    ref_bw.to_str().unwrap(),
                ])
                .status()
                .unwrap();

            // bedtools emits no per-strand track when there is zero coverage;
            // RustQC likewise skips the file. Treat "no reference intervals"
            // and "no RustQC file" as consistent.
            let ref_bg = root.join(format!("ref.{threads}.{suffix}.dec.bg"));
            let has_ref =
                bigwig_to_bedgraph(&ref_bw, &ref_bg) && !read_bedgraph(&ref_bg).is_empty();

            if !has_ref {
                continue;
            }
            let rust_bw = outdir.join(format!("bigwig/{suffix}"));
            assert!(rust_bw.exists(), "missing RustQC bigWig for {label}");
            let rust_bg = root.join(format!("rust.{threads}.{suffix}.dec.bg"));
            assert!(
                bigwig_to_bedgraph(&rust_bw, &rust_bg),
                "decode failed {label}"
            );
            assert_eq!(
                read_bedgraph(&ref_bg),
                read_bedgraph(&rust_bg),
                "bedGraph mismatch for {label}"
            );
        }
    }

    let _ = fs::remove_dir_all(root);
}
