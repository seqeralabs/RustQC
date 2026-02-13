//! Integration tests comparing RustQC output against R dupRadar reference output.
//!
//! These tests run RustQC on the same test BAM/GTF data used to generate
//! the R reference outputs in tests/expected/.

use std::fs;
use std::path::Path;
use std::process::Command;

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

/// Helper: run rustqc rna with a specified input file and return the output
fn run_rustqc_with_input(outdir: &str, input: &str) -> std::process::Output {
    let binary = rustqc_binary();
    Command::new(&binary)
        .args([
            "rna",
            input,
            "--gtf",
            "tests/data/test.gtf",
            "--outdir",
            outdir,
            "--skip-dup-check",
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
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() == 2 {
            match parts[0] {
                "intercept" => intercept = parts[1].parse().unwrap(),
                "slope" => slope = parts[1].parse().unwrap(),
                _ => {}
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
    let rust_matrix = parse_dup_matrix(&format!("{}/test_dupMatrix.txt", outdir));

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
        parse_intercept_slope(&format!("{}/test_intercept_slope.txt", outdir));

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
        "test_dupMatrix.txt",
        "test_intercept_slope.txt",
        "test_duprateExpDens.png",
        "test_duprateExpDens.svg",
        "test_duprateExpBoxplot.png",
        "test_duprateExpBoxplot.svg",
        "test_expressionHist.png",
        "test_expressionHist.svg",
        "test_dup_intercept_mqc.txt",
        "test_duprateExpDensCurve_mqc.txt",
    ];

    for file in &expected_files {
        let path = format!("{}/{}", outdir, file);
        assert!(Path::new(&path).exists(), "Missing output file: {}", file);

        // Check file is non-empty
        let metadata = fs::metadata(&path).unwrap();
        assert!(metadata.len() > 0, "Output file is empty: {}", file);
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

    let content = fs::read_to_string(format!("{}/test_dupMatrix.txt", outdir)).unwrap();
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
    let content = fs::read_to_string(format!("{}/test_dup_intercept_mqc.txt", outdir)).unwrap();
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
    let content =
        fs::read_to_string(format!("{}/test_duprateExpDensCurve_mqc.txt", outdir)).unwrap();
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
    let rust_matrix = parse_dup_matrix(&format!("{}/test_dupMatrix.txt", outdir));

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
    let rust_matrix = parse_dup_matrix(&format!("{}/test_dupMatrix.txt", outdir));

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
        "test_dupMatrix.txt",
        "test_duprateExpDens.png",
        "test_duprateExpBoxplot.png",
        "test_expressionHist.png",
        "test_copy_dupMatrix.txt",
        "test_copy_duprateExpDens.png",
        "test_copy_duprateExpBoxplot.png",
        "test_copy_expressionHist.png",
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
        let rust_matrix = parse_dup_matrix(&format!("{}/{}_dupMatrix.txt", outdir, stem));
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
    let bin = env!("CARGO_BIN_EXE_rustqc");
    let output = Command::new(bin)
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
    let bam_matrix = parse_dup_matrix(&format!("{}/test_dupMatrix.txt", outdir_bam));
    let sam_matrix = parse_dup_matrix(&format!("{}/test_dupMatrix.txt", outdir_sam));

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
