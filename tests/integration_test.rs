//! Integration tests comparing RustQC output against R dupRadar reference output.
//!
//! These tests run RustQC on the same test BAM/GTF data used to generate
//! the R reference outputs in tests/expected/.

use std::fs;
use std::path::Path;
use std::path::PathBuf;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use rustqc::rna::bam;
use rustqc::rna::bam::header::{HeaderRecord, HeaderView};

/// Helper: get the path to the rustqc binary.
///
/// Cargo sets `CARGO_BIN_EXE_rustqc` at compile time for integration tests, so
/// the path is always available. (At run time the variable is only set by newer
/// Cargo versions, which made these tests fail on the MSRV toolchain.)
fn rustqc_binary() -> String {
    env!("CARGO_BIN_EXE_rustqc").to_string()
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
                .push_tag(b"LN", len.to_string()),
        );
    }

    let header_view = HeaderView::from_header(&header);
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

    // Same layout as the nf-core dupradar module's file: a YAML header in `#`
    // lines, then space-separated "Sample dupRadar_intercept" and value lines.
    let content =
        fs::read_to_string(format!("{}/dupradar/test_dup_intercept_mqc.txt", outdir)).unwrap();
    assert!(content.starts_with("#id: DupInt\n#plot_type: 'generalstats'\n"));
    let data_lines: Vec<&str> = content.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(data_lines.len(), 2, "MultiQC intercept file: {content}");
    assert_eq!(data_lines[0], "Sample dupRadar_intercept");
    let parts: Vec<&str> = data_lines[1].split(' ').collect();
    assert_eq!(parts.len(), 2, "MQC data line should have 2 columns");
    assert_eq!(parts[0], "test");
    let intercept: f64 = parts[1].parse().expect("MQC intercept should be numeric");

    // The same intercept as in the R fixture.
    let (r_intercept, _) = parse_intercept_slope("tests/expected/intercept_slope.txt");
    assert!(
        approx_eq(r_intercept, intercept, 1e-6),
        "MQC intercept {intercept} vs R {r_intercept}"
    );

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

    // Same layout as the nf-core dupradar module's file: a YAML header in `#`
    // lines, then "RPK duplication-percent" rows (space-separated, no header).
    let content = fs::read_to_string(format!(
        "{}/dupradar/test_duprateExpDensCurve_mqc.txt",
        outdir
    ))
    .unwrap();
    assert!(content.starts_with("#id: dupradar\n#plot_type: 'linegraph'\n"));
    let data_lines: Vec<&str> = content.lines().filter(|l| !l.starts_with('#')).collect();
    // One row per 10 expressed genes (RPK > 0), starting with the least expressed.
    let matrix = fs::read_to_string(format!("{}/dupradar/test_dupMatrix.txt", outdir)).unwrap();
    let rpk_col = matrix
        .lines()
        .next()
        .unwrap()
        .split('\t')
        .position(|h| h == "RPK")
        .unwrap();
    let expressed = matrix
        .lines()
        .skip(1)
        .filter(|l| l.split('\t').nth(rpk_col).unwrap() != "0")
        .count();
    assert!(expressed > 0);
    assert_eq!(data_lines.len(), expressed.div_ceil(10), "{content}");
    let mut prev_x = 0.0;
    for line in &data_lines {
        let parts: Vec<&str> = line.split(' ').collect();
        assert_eq!(parts.len(), 2, "Curve data should have 2 columns: {line}");
        let x: f64 = parts[0].parse().expect("X should be numeric");
        let y: f64 = parts[1].parse().expect("Y should be numeric");
        // Rows are genes sorted by expression; y is a percentage.
        assert!(x > 0.0 && x >= prev_x, "{line}");
        assert!((0.0..=100.0).contains(&y), "{line}");
        prev_x = x;
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

// ===================================================================
// split_bam — parity with RSeQC 5.0.4 split_bam.py
// ===================================================================
//
// tests/expected/split_bam/ holds the real split_bam.py output for the
// tests/data/split_bam_test.{bam,bed} fixture; regenerate both with
// tests/create_split_bam_test_data.sh.

const SPLIT_BAM_FIXTURE: &str = "tests/data/split_bam_test.bam";
const SPLIT_BAM_BED: &str = "tests/data/split_bam_test.bed";
const SPLIT_BAM_EXPECTED: &str = "tests/expected/split_bam";

/// Run rustqc on the split_bam fixture with the given extra arguments.
fn run_split_bam(outdir: &Path, extra: &[&str]) -> std::process::Output {
    let mut args = vec![
        "rna",
        SPLIT_BAM_FIXTURE,
        "--gtf",
        "tests/data/test.gtf",
        "--skip-dup-check",
        "--outdir",
        outdir.to_str().unwrap(),
    ];
    args.extend_from_slice(extra);
    Command::new(rustqc_binary())
        .args(&args)
        .output()
        .expect("Failed to execute rustqc")
}

/// Read a BAM's header text and all records.
fn read_bam(path: &Path) -> (noodles::sam::Header, Vec<bam::Record>) {
    use rustqc::rna::bam::Read;
    let mut reader = bam::Reader::from_path(path)
        .unwrap_or_else(|e| panic!("cannot open {}: {e}", path.display()));
    let header = reader.header().noodles_header().clone();
    let mut records = Vec::new();
    let mut record = bam::Record::new();
    while let Some(result) = reader.read(&mut record) {
        result.unwrap();
        records.push(record.clone());
    }
    (header, records)
}

/// Assert the split_bam summary is byte-identical to RSeQC's stdout.
fn assert_split_bam_stats_match(split_dir: &Path) {
    let got = fs::read_to_string(split_dir.join("split_bam_test.split_bam.txt"))
        .expect("split_bam summary not written");
    let expected =
        fs::read_to_string(Path::new(SPLIT_BAM_EXPECTED).join("split_bam_test.split_bam.txt"))
            .unwrap();
    assert_eq!(got, expected, "split_bam summary differs from RSeQC");
}

/// Counts and the three BAMs must equal RSeQC's, record for record and in
/// order, on both the sequential (-t 1) and parallel indexed (-t 2) paths.
#[test]
fn test_split_bam_matches_rseqc() {
    for threads in ["1", "2"] {
        let root = unique_test_dir(&format!("split-bam-t{threads}"));
        let output = run_split_bam(
            &root,
            &[
                "--split-bam-bed",
                SPLIT_BAM_BED,
                "--split-bam-write-bams",
                "--threads",
                threads,
            ],
        );
        assert!(
            output.status.success(),
            "rustqc failed (-t {threads}):\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let split_dir = root.join("rseqc").join("split_bam");
        assert_split_bam_stats_match(&split_dir);

        for class in ["in", "ex", "junk"] {
            let name = format!("split_bam_test.{class}.bam");
            let (got_header, got) = read_bam(&split_dir.join(&name));
            let (exp_header, expected) = read_bam(&Path::new(SPLIT_BAM_EXPECTED).join(&name));
            assert_eq!(
                got_header, exp_header,
                "{name} header differs (-t {threads})"
            );
            assert_eq!(
                got.len(),
                expected.len(),
                "{name} record count differs (-t {threads})"
            );
            for (i, (g, e)) in got.iter().zip(expected.iter()).enumerate() {
                assert!(
                    g == e,
                    "{name} record {i} differs (-t {threads}): got {}, expected {}",
                    String::from_utf8_lossy(g.qname()),
                    String::from_utf8_lossy(e.qname())
                );
            }
        }
        let _ = fs::remove_dir_all(root);
    }
}

/// Config-file activation, no BAMs by default, and BED names matched against
/// alignment names even when a GTF chromosome mapping renames a contig
/// (chr2 → gtf_chr2 for GTF lookups; the BED's chr2 blocks must still apply).
#[test]
fn test_split_bam_config_and_counts_only() {
    let root = unique_test_dir("split-bam-config");
    let config_path = root.join("rustqc.yml");
    fs::write(
        &config_path,
        format!(
            "rna:\n  chromosome_mapping:\n    gtf_chr2: chr2\n  split_bam:\n    bed: {}\n",
            fs::canonicalize(SPLIT_BAM_BED).unwrap().display()
        ),
    )
    .unwrap();
    let outdir = root.join("out");
    let output = run_split_bam(&outdir, &["-c", config_path.to_str().unwrap()]);
    assert!(
        output.status.success(),
        "rustqc failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let split_dir = outdir.join("rseqc").join("split_bam");
    assert_split_bam_stats_match(&split_dir);
    for class in ["in", "ex", "junk"] {
        assert!(
            !split_dir
                .join(format!("split_bam_test.{class}.bam"))
                .exists(),
            "{class}.bam must not be written unless requested"
        );
    }
    let _ = fs::remove_dir_all(root);
}

/// split_bam stays off when no BED file is configured.
#[test]
fn test_split_bam_off_without_bed() {
    let root = unique_test_dir("split-bam-off");
    let output = run_split_bam(&root, &[]);
    assert!(
        output.status.success(),
        "rustqc failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        !root.join("rseqc").join("split_bam").exists(),
        "split_bam output written without a BED file"
    );
    let _ = fs::remove_dir_all(root);
}

// ===================================================================
// samtools flagstat / stats parity for QC-failed and dual-flagged records
//
// Cases adapted from seqeralabs/RustQC#175 (Nico DeVeaux) and #131
// (Benjamin Demaille). Expected text is samtools 1.22.1 output.
// ===================================================================

/// The two plain mapped reads every flagstat case sits between, so that some
/// read lands in a gene and `rustqc rna` produces its outputs.
const FLAGSTAT_BACKGROUND: [&str; 2] = [
    "r0\t0\tchr1\t1001\t60\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\t??????????????????????????????????????????????????",
    "r1\t0\tchr1\t1201\t60\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\t??????????????????????????????????????????????????",
];

/// Run `rustqc rna` on the background reads plus one record with `flag` (and
/// mate field `rnext`) and return its flagstat and stats outputs.
fn samtools_for_one_record(label: &str, flag: u16, rnext: &str) -> (String, String) {
    let root = unique_test_dir(label);
    let bam_path = root.join("s.bam");
    let gtf_path = root.join("a.gtf");
    let outdir = root.join("out");
    fs::create_dir_all(&outdir).unwrap();
    fs::write(
        &gtf_path,
        "chr1\tt\tgene\t900\t3000\t.\t+\t.\tgene_id \"g1\"; gene_name \"g1\";\n\
         chr1\tt\texon\t900\t3000\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\"; exon_number \"1\";\n",
    )
    .unwrap();
    let pnext = if rnext == "*" { 0 } else { 1001 };
    let case = format!(
        "r2\t{flag}\tchr1\t1001\t60\t50M\t{rnext}\t{pnext}\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\t??????????????????????????????????????????????????"
    );
    write_bam_fixture(
        &bam_path,
        &[("chr1", 100_000), ("chr2", 100_000)],
        &[FLAGSTAT_BACKGROUND[0], &case, FLAGSTAT_BACKGROUND[1]],
    );

    let output = run_rustqc_custom(&outdir, &bam_path, &gtf_path, 1, true, &[]);
    assert!(
        output.status.success(),
        "rustqc failed on flag {flag:#06x}:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let flagstat = fs::read_to_string(outdir.join("samtools/s.flagstat")).unwrap();
    let stats = fs::read_to_string(outdir.join("samtools/s.stats")).unwrap();
    let _ = fs::remove_dir_all(root);
    (flagstat, stats)
}

fn flagstat_for_one_record(label: &str, flag: u16, rnext: &str) -> String {
    samtools_for_one_record(label, flag, rnext).0
}

/// A QC-failed read (0x200) belongs in samtools' second column.
#[test]
fn test_flagstat_qc_failed_column() {
    let expected = "\
2 + 1 in total (QC-passed reads + QC-failed reads)
2 + 1 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
2 + 1 mapped (100.00% : 100.00%)
2 + 1 primary mapped (100.00% : 100.00%)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
";
    assert_eq!(flagstat_for_one_record("fs-qcfail", 0x200, "*"), expected);
}

/// "properly paired" needs the read itself mapped (samtools bam_stat.c),
/// so PAIRED|PROPER_PAIR|UNMAP|MUNMAP|READ1 (0x4f) does not count.
#[test]
fn test_flagstat_properly_paired_requires_mapped_read() {
    let expected = "\
3 + 0 in total (QC-passed reads + QC-failed reads)
3 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
2 + 0 mapped (66.67% : N/A)
2 + 0 primary mapped (66.67% : N/A)
1 + 0 paired in sequencing
1 + 0 read1
0 + 0 read2
0 + 0 properly paired (0.00% : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
";
    assert_eq!(
        flagstat_for_one_record("fs-proper-unmapped", 0x4f, "="),
        expected
    );
}

/// A QC-failed duplicate (0x600) is a duplicate in the QC-failed column, so
/// "primary duplicates" never exceeds "duplicates".
#[test]
fn test_flagstat_qc_failed_duplicate_counts_as_duplicate() {
    let expected = "\
2 + 1 in total (QC-passed reads + QC-failed reads)
2 + 1 primary
0 + 0 secondary
0 + 0 supplementary
0 + 1 duplicates
0 + 1 primary duplicates
2 + 1 mapped (100.00% : 100.00%)
2 + 1 primary mapped (100.00% : 100.00%)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
";
    assert_eq!(
        flagstat_for_one_record("fs-qcfail-dup", 0x600, "*"),
        expected
    );
}

/// SECONDARY takes priority over SUPPLEMENTARY, in flagstat and in the
/// `stats` summary numbers alike (samtools 1.22.1 `flagstat_loop` and
/// `stats.c`): a record with both flags is secondary only.
#[test]
fn test_secondary_supplementary_overlap_counts_secondary_only() {
    let (flagstat, stats) = samtools_for_one_record("fs-sec-supp", 0x900, "*");
    let expected = "\
3 + 0 in total (QC-passed reads + QC-failed reads)
2 + 0 primary
1 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
3 + 0 mapped (100.00% : N/A)
2 + 0 primary mapped (100.00% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
";
    assert_eq!(flagstat, expected);
    assert!(
        stats.contains("SN\tnon-primary alignments:\t1\n"),
        "{stats}"
    );
    assert!(
        stats.contains("SN\tsupplementary alignments:\t0\n"),
        "{stats}"
    );
}

/// A supplementary-only record is not a "non-primary alignment" in `stats`.
#[test]
fn test_stats_supplementary_is_not_non_primary() {
    let (_, stats) = samtools_for_one_record("stats-supp-only", 0x800, "*");
    assert!(
        stats.contains("SN\tnon-primary alignments:\t0\n"),
        "{stats}"
    );
    assert!(
        stats.contains("SN\tsupplementary alignments:\t1\n"),
        "{stats}"
    );
}

/// A dual-flagged QC-failed duplicate is secondary in the failed column,
/// stays a duplicate, and contributes to neither supplementary nor primary.
#[test]
fn test_flagstat_secondary_supplementary_qcfail_duplicate() {
    let (flagstat, stats) = samtools_for_one_record("fs-sec-supp-qcfail-dup", 0xf00, "*");
    let expected = "\
2 + 1 in total (QC-passed reads + QC-failed reads)
2 + 0 primary
0 + 1 secondary
0 + 0 supplementary
0 + 1 duplicates
0 + 0 primary duplicates
2 + 1 mapped (100.00% : 100.00%)
2 + 0 primary mapped (100.00% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
";
    assert_eq!(flagstat, expected);
    assert!(
        stats.contains("SN\tnon-primary alignments:\t1\n"),
        "{stats}"
    );
    assert!(
        stats.contains("SN\tsupplementary alignments:\t0\n"),
        "{stats}"
    );
}

/// SAM lines for every 12-bit flag in three placements: mate on the same
/// contig, mate on the other contig (MAPQ 60), and mate on the other contig
/// with MAPQ 3. `tests/expected/samtools/all_flags.*` is samtools 1.22.1's
/// output for exactly these records (see the README there).
fn all_flags_sam_lines() -> Vec<String> {
    let mut chr1 = Vec::new();
    let mut chr2 = Vec::new();
    for flag in 0..4096u16 {
        chr1.push(format!(
            "f{flag}a\t{flag}\tchr1\t10\t60\t4M\t=\t10\t0\tACGT\tIIII"
        ));
        chr1.push(format!(
            "f{flag}b\t{flag}\tchr1\t10\t60\t4M\tchr2\t10\t0\tACGT\tIIII"
        ));
        chr2.push(format!(
            "f{flag}c\t{flag}\tchr2\t10\t3\t4M\tchr1\t10\t0\tACGT\tIIII"
        ));
    }
    chr1.extend(chr2);
    chr1
}

/// Every flag, against the unmodified samtools 1.22.1 golden output, on both
/// the serial and per-contig parallel paths: all of `flagstat` and the `SN`
/// section of `stats`.
#[test]
fn test_samtools_flagstat_and_stats_match_on_every_flag() {
    let expected_flagstat =
        fs::read_to_string("tests/expected/samtools/all_flags.flagstat").unwrap();
    let expected_sn = fs::read_to_string("tests/expected/samtools/all_flags.stats.sn").unwrap();

    let lines = all_flags_sam_lines();
    let lines: Vec<&str> = lines.iter().map(String::as_str).collect();
    for threads in [1, 2] {
        let root = unique_test_dir(&format!("fs-all-flags-t{threads}"));
        let bam_path = root.join("all_flags.bam");
        let gtf_path = root.join("all_flags.gtf");
        let outdir = root.join("out");
        fs::create_dir_all(&outdir).unwrap();
        write_test_gtf(&gtf_path, &["chr1", "chr2"]);
        write_bam_fixture(&bam_path, &[("chr1", 1000), ("chr2", 1000)], &lines);

        let output = run_rustqc_custom(&outdir, &bam_path, &gtf_path, threads, true, &[]);
        assert!(
            output.status.success(),
            "rustqc failed with --threads {threads}:\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let flagstat = fs::read_to_string(outdir.join("samtools/all_flags.flagstat")).unwrap();
        assert_eq!(
            flagstat, expected_flagstat,
            "flagstat with --threads {threads}"
        );
        let stats = fs::read_to_string(outdir.join("samtools/all_flags.stats")).unwrap();
        let sn: String = stats
            .lines()
            .filter(|l| l.starts_with("SN\t"))
            .map(|l| format!("{l}\n"))
            .collect();
        assert_eq!(sn, expected_sn, "stats SN section with --threads {threads}");
        let _ = fs::remove_dir_all(root);
    }
}

// ===================================================================
// Unmapped sweep on the parallel path
// ===================================================================

/// A placed unmapped mate (0x4 with a reference) in the last reference's final
/// 16 kb, plus two unplaced reads. The parallel path's unmapped sweep must not
/// count the placed mate a second time. Expected text is `samtools flagstat`
/// 1.22.1 output for this input.
#[test]
fn test_parallel_unmapped_sweep_counts_placed_unmapped_mate_once() {
    let expected = "\
6 + 0 in total (QC-passed reads + QC-failed reads)
6 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
3 + 0 mapped (50.00% : N/A)
3 + 0 primary mapped (50.00% : N/A)
2 + 0 paired in sequencing
1 + 0 read1
1 + 0 read2
0 + 0 properly paired (0.00% : N/A)
0 + 0 with itself and mate mapped
1 + 0 singletons (50.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
";
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
    let qual = "??????????????????????????????????????????????????";
    let lines: Vec<String> = [
        "r0\t0\tchr1\t1001\t60\t50M\t*\t0\t0",
        "r1\t0\tchr1\t1201\t60\t50M\t*\t0\t0",
        "r2\t73\tchr2\t99001\t60\t50M\t=\t99001\t0",
        "r3\t133\tchr2\t99001\t0\t*\t=\t99001\t0",
        "u1\t4\t*\t0\t0\t*\t*\t0\t0",
        "u2\t4\t*\t0\t0\t*\t*\t0\t0",
    ]
    .iter()
    .map(|fields| format!("{fields}\t{seq}\t{qual}"))
    .collect();
    let lines: Vec<&str> = lines.iter().map(String::as_str).collect();

    for threads in [1, 4] {
        let root = unique_test_dir(&format!("unmapped-sweep-t{threads}"));
        let bam_path = root.join("s.bam");
        let gtf_path = root.join("a.gtf");
        let outdir = root.join("out");
        fs::create_dir_all(&outdir).unwrap();
        fs::write(
            &gtf_path,
            "chr1\tt\tgene\t900\t3000\t.\t+\t.\tgene_id \"g1\"; gene_name \"g1\";\n\
             chr1\tt\texon\t900\t3000\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\"; exon_number \"1\";\n",
        )
        .unwrap();
        write_bam_fixture(&bam_path, &[("chr1", 100_000), ("chr2", 100_000)], &lines);

        let output = run_rustqc_custom(&outdir, &bam_path, &gtf_path, threads, true, &[]);
        assert!(
            output.status.success(),
            "rustqc failed with --threads {threads}:\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let flagstat = fs::read_to_string(outdir.join("samtools/s.flagstat")).unwrap();
        assert_eq!(flagstat, expected, "flagstat with --threads {threads}");
        let _ = fs::remove_dir_all(root);
    }
}

/// The same records as a samtools-written CRAM holding all of them in one
/// multi-reference slice: each per-reference fetch on the parallel path must
/// return only its own reference's records. Expected text is `samtools
/// flagstat` 1.22.1 output for this file.
#[test]
fn test_parallel_cram_multi_reference_slice_counts_each_record_once() {
    let expected = "\
6 + 0 in total (QC-passed reads + QC-failed reads)
6 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
3 + 0 mapped (50.00% : N/A)
3 + 0 primary mapped (50.00% : N/A)
2 + 0 paired in sequencing
1 + 0 read1
1 + 0 read2
0 + 0 properly paired (0.00% : N/A)
0 + 0 with itself and mate mapped
1 + 0 singletons (50.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
";
    let data = Path::new("tests/data/unmapped_sweep");
    for threads in [1, 4] {
        let root = unique_test_dir(&format!("cram-multi-ref-t{threads}"));
        let gtf_path = root.join("a.gtf");
        let outdir = root.join("out");
        fs::create_dir_all(&outdir).unwrap();
        fs::write(
            &gtf_path,
            "chr1\tt\tgene\t1\t1000\t.\t+\t.\tgene_id \"g1\"; gene_name \"g1\";\n\
             chr1\tt\texon\t1\t1000\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\"; exon_number \"1\";\n",
        )
        .unwrap();

        let output = Command::new(rustqc_binary())
            .args(["rna", "--skip-dup-check", "--threads", &threads.to_string()])
            .arg("--reference")
            .arg(data.join("sweep.fa"))
            .arg("--gtf")
            .arg(&gtf_path)
            .arg("--outdir")
            .arg(&outdir)
            .arg(data.join("multi_ref_slice.cram"))
            .output()
            .expect("Failed to execute rustqc");
        assert!(
            output.status.success(),
            "rustqc failed with --threads {threads}:\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let flagstat =
            fs::read_to_string(outdir.join("samtools/multi_ref_slice.flagstat")).unwrap();
        assert_eq!(flagstat, expected, "flagstat with --threads {threads}");
        let _ = fs::remove_dir_all(root);
    }
}

/// Multiple CRAI slice entries for one reference/container must not duplicate reads.
#[test]
fn test_parallel_cram_same_reference_multi_slice_counts_each_record_once() {
    let expected = "\
6 + 0 in total (QC-passed reads + QC-failed reads)
6 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
4 + 0 mapped (66.67% : N/A)
4 + 0 primary mapped (66.67% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
";
    let data = Path::new("tests/data/unmapped_sweep");
    for threads in [1, 4] {
        let root = unique_test_dir(&format!("cram-same-reference-multi-slice-t{threads}"));
        let gtf_path = root.join("a.gtf");
        let outdir = root.join("out");
        fs::write(&gtf_path,
            "chr1\tt\tgene\t1\t1000\t.\t+\t.\tgene_id \"g1\"; gene_name \"g1\";\n\
             chr1\tt\texon\t1\t1000\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\"; exon_number \"1\";\n"
        ).unwrap();
        let output = Command::new(rustqc_binary())
            .args(["rna", "--skip-dup-check", "--threads", &threads.to_string()])
            .arg("--reference")
            .arg(data.join("sweep.fa"))
            .arg("--gtf")
            .arg(&gtf_path)
            .arg("--outdir")
            .arg(&outdir)
            .arg(data.join("same_reference_multi_slice.cram"))
            .output()
            .expect("Failed to execute rustqc");
        assert!(
            output.status.success(),
            "rustqc failed with --threads {threads}:\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let flagstat =
            fs::read_to_string(outdir.join("samtools/same_reference_multi_slice.flagstat"))
                .unwrap();
        assert_eq!(flagstat, expected, "flagstat with --threads {threads}");
        fs::remove_dir_all(root).unwrap();
    }
}

/// A CRAM with no unplaced records must succeed on the parallel sweep path too.
#[test]
fn test_parallel_cram_with_no_unplaced_reads() {
    let expected = "\
4 + 0 in total (QC-passed reads + QC-failed reads)
4 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
3 + 0 mapped (75.00% : N/A)
3 + 0 primary mapped (75.00% : N/A)
2 + 0 paired in sequencing
1 + 0 read1
1 + 0 read2
0 + 0 properly paired (0.00% : N/A)
0 + 0 with itself and mate mapped
1 + 0 singletons (50.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
";
    let data = Path::new("tests/data/unmapped_sweep");
    for threads in [1, 4] {
        let root = unique_test_dir(&format!("cram-no-unplaced-t{threads}"));
        let gtf_path = root.join("a.gtf");
        let outdir = root.join("out");
        fs::write(&gtf_path,
            "chr1\tt\tgene\t1\t1000\t.\t+\t.\tgene_id \"g1\"; gene_name \"g1\";\n\
             chr1\tt\texon\t1\t1000\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\"; exon_number \"1\";\n"
        ).unwrap();
        let output = Command::new(rustqc_binary())
            .args(["rna", "--skip-dup-check", "--threads", &threads.to_string()])
            .arg("--reference")
            .arg(data.join("sweep.fa"))
            .arg("--gtf")
            .arg(&gtf_path)
            .arg("--outdir")
            .arg(&outdir)
            .arg(data.join("placed_only.cram"))
            .output()
            .expect("Failed to execute rustqc");
        assert!(
            output.status.success(),
            "rustqc failed with --threads {threads}:\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let flagstat = fs::read_to_string(outdir.join("samtools/placed_only.flagstat")).unwrap();
        assert_eq!(flagstat, expected, "flagstat with --threads {threads}");
        fs::remove_dir_all(root).unwrap();
    }
}
