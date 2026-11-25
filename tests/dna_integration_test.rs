//! Parity tests for the DNA pipeline against the committed reference outputs.
//!
//! The fixtures under `tests/expected/dna/` are the output of the upstream
//! tools themselves, at the versions pinned in `VERSIONS.txt`. A failure here
//! is a defect in RustQC, not a reason to regenerate the fixture.
//!
//! Compressed outputs are compared on their decompressed bytes. Two bgzf
//! writers at the same compression level need not emit identical compressed
//! bytes, so comparing the `.gz` files directly would test the compressor
//! rather than this code.

use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};

use rust_htslib::bam::Read as BamRead;
use rust_htslib::{bam, bgzf};

use rustqc::dna::depth::{DepthAccum, MOSDEPTH_DEFAULT_EXCLUDE};
use rustqc::dna::insert_size::{self, InsertSizeAccum};
use rustqc::dna::mosdepth::{output, ContigDepth, MosdepthResult};
use rustqc::dna::wgs_metrics::{self, WgsAccum, WgsMetricsResult};

/// Window size and thresholds the fixtures were generated with.
const WINDOW_SIZE: u32 = 500;
const THRESHOLDS: [u32; 7] = [1, 5, 10, 15, 20, 30, 50];

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/expected/dna")
        .join(name)
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join("rustqc-dna-parity");
    std::fs::create_dir_all(&dir).unwrap();
    dir.join(name)
}

/// Run the depth engine over the committed test BAM and summarise it exactly
/// as the fixtures were generated.
fn compute() -> MosdepthResult {
    let bam_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/dna/test.dna.bam");
    let reader = bam::Reader::from_path(&bam_path).unwrap();
    let header = reader.header().to_owned();

    let mut contigs = Vec::new();
    for tid in 0..header.target_count() {
        let name = String::from_utf8(header.tid2name(tid).to_vec()).unwrap();
        let length = header.target_len(tid).unwrap();
        let mut accum = DepthAccum::new(length, 0, MOSDEPTH_DEFAULT_EXCLUDE);

        let mut record = bam::Record::new();
        let mut per_contig = bam::Reader::from_path(&bam_path).unwrap();
        while let Some(result) = per_contig.read(&mut record) {
            result.unwrap();
            if record.tid() == tid as i32 {
                accum.process_read(&record);
            }
        }
        let depths = accum.into_depths();
        contigs.push(ContigDepth::from_depths(
            &name,
            &depths,
            Some(WINDOW_SIZE),
            &THRESHOLDS,
        ));
    }

    MosdepthResult {
        contigs,
        window_size: Some(WINDOW_SIZE),
        thresholds: THRESHOLDS.to_vec(),
    }
}

fn read_bgzf(path: &Path) -> String {
    let mut reader = bgzf::Reader::from_path(path).unwrap();
    let mut buf = Vec::new();
    reader.read_to_end(&mut buf).unwrap();
    String::from_utf8(buf).unwrap()
}

/// Compare line by line so a failure names the offending row.
fn assert_same_lines(actual: &str, expected: &str, what: &str) {
    let a: Vec<&str> = actual.lines().collect();
    let e: Vec<&str> = expected.lines().collect();
    for (i, (got, want)) in a.iter().zip(e.iter()).enumerate() {
        assert_eq!(got, want, "{what}: line {} differs", i + 1);
    }
    assert_eq!(a.len(), e.len(), "{what}: line count differs");
}

#[test]
fn fixture_tool_versions_are_the_pinned_ones() {
    let versions = std::fs::read_to_string(fixture("VERSIONS.txt")).unwrap();
    assert!(
        versions.contains("mosdepth\t0.3.14"),
        "unexpected mosdepth fixture version: {versions}"
    );
    assert!(
        versions.contains("samtools\t1.24"),
        "unexpected samtools fixture version: {versions}"
    );
}

#[test]
fn summary_matches_mosdepth() {
    let path = scratch("test.mosdepth.summary.txt");
    output::write_summary(&compute(), &path).unwrap();
    assert_same_lines(
        &std::fs::read_to_string(&path).unwrap(),
        &std::fs::read_to_string(fixture("test.mosdepth.summary.txt")).unwrap(),
        "summary",
    );
}

#[test]
fn global_dist_matches_mosdepth() {
    let path = scratch("test.mosdepth.global.dist.txt");
    output::write_global_dist(&compute(), &path).unwrap();
    assert_same_lines(
        &std::fs::read_to_string(&path).unwrap(),
        &std::fs::read_to_string(fixture("test.mosdepth.global.dist.txt")).unwrap(),
        "global dist",
    );
}

#[test]
fn region_dist_matches_mosdepth() {
    let path = scratch("test.mosdepth.region.dist.txt");
    output::write_region_dist(&compute(), &path).unwrap();
    assert_same_lines(
        &std::fs::read_to_string(&path).unwrap(),
        &std::fs::read_to_string(fixture("test.mosdepth.region.dist.txt")).unwrap(),
        "region dist",
    );
}

#[test]
fn per_base_matches_mosdepth() {
    let path = scratch("test.per-base.bed.gz");
    output::write_per_base(&compute(), &path).unwrap();
    assert_same_lines(
        &read_bgzf(&path),
        &read_bgzf(&fixture("test.per-base.bed.gz")),
        "per-base",
    );
}

#[test]
fn regions_match_mosdepth() {
    let path = scratch("test.regions.bed.gz");
    output::write_regions(&compute(), &path).unwrap();
    assert_same_lines(
        &read_bgzf(&path),
        &read_bgzf(&fixture("test.regions.bed.gz")),
        "regions",
    );
}

#[test]
fn thresholds_match_mosdepth() {
    let path = scratch("test.thresholds.bed.gz");
    output::write_thresholds(&compute(), &path).unwrap();
    assert_same_lines(
        &read_bgzf(&path),
        &read_bgzf(&fixture("test.thresholds.bed.gz")),
        "thresholds",
    );
}

/// The depth histogram is the input to both distribution files, so pinning it
/// separately makes a distribution failure easy to attribute.
#[test]
fn depth_histogram_matches_the_per_base_fixture() {
    let result = compute();
    let mut expected: BTreeMap<u32, u64> = BTreeMap::new();
    for line in read_bgzf(&fixture("test.per-base.bed.gz")).lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        let start: u64 = fields[1].parse().unwrap();
        let end: u64 = fields[2].parse().unwrap();
        let depth: u32 = fields[3].parse().unwrap();
        *expected.entry(depth).or_insert(0) += end - start;
    }
    assert_eq!(result.contigs[0].histogram, expected);
}

// ===================================================================
// End-to-end parity: the binary, not just the library
// ===================================================================

/// Run `rustqc dna` once into a scratch directory shared by every end-to-end
/// test, with the same window size and thresholds the fixtures were made with.
fn run_binary() -> &'static Path {
    static OUTDIR: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
    OUTDIR.get_or_init(|| {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"));
        let outdir = std::env::temp_dir().join("rustqc-dna-e2e");
        let _ = std::fs::remove_dir_all(&outdir);
        std::fs::create_dir_all(&outdir).unwrap();

        let status = std::process::Command::new(env!("CARGO_BIN_EXE_rustqc"))
            .arg("dna")
            .arg(root.join("tests/data/dna/test.dna.bam"))
            .arg("--outdir")
            .arg(&outdir)
            .arg("--window-size")
            .arg(WINDOW_SIZE.to_string())
            .arg("--reference")
            .arg(root.join("tests/data/dna/genome.fasta"))
            .arg("--quiet")
            .status()
            .expect("failed to run the rustqc binary");
        assert!(status.success(), "rustqc dna exited with {status}");
        outdir
    })
}

/// The sample name is the BAM file stem, dots included.
const SAMPLE: &str = "test.dna";

fn produced(subdir: &str, name: &str) -> PathBuf {
    run_binary().join(subdir).join(name)
}

#[test]
fn binary_writes_every_mosdepth_output_byte_for_byte() {
    for suffix in [
        "mosdepth.summary.txt",
        "mosdepth.global.dist.txt",
        "mosdepth.region.dist.txt",
    ] {
        let got = std::fs::read_to_string(produced("mosdepth", &format!("{SAMPLE}.{suffix}")))
            .unwrap_or_else(|e| panic!("reading {suffix}: {e}"));
        let want = std::fs::read_to_string(fixture(&format!("test.{suffix}"))).unwrap();
        assert_same_lines(&got, &want, suffix);
    }
    for suffix in ["per-base.bed.gz", "regions.bed.gz", "thresholds.bed.gz"] {
        let got = read_bgzf(&produced("mosdepth", &format!("{SAMPLE}.{suffix}")));
        let want = read_bgzf(&fixture(&format!("test.{suffix}")));
        assert_same_lines(&got, &want, suffix);
    }
}

#[test]
fn binary_writes_flagstat_and_idxstats_byte_for_byte() {
    for suffix in ["flagstat", "idxstats"] {
        let got = std::fs::read_to_string(produced("samtools", &format!("{SAMPLE}.{suffix}.txt")))
            .unwrap();
        let want = std::fs::read_to_string(fixture(&format!("test.{suffix}.txt"))).unwrap();
        assert_eq!(got, want, "{suffix} must match samtools exactly");
    }
}

/// `samtools stats` output is compared on its data lines only. RustQC writes
/// its own `#` header, naming itself rather than reproducing samtools' command
/// line and version banner, which is deliberate and shared with the `rna`
/// pipeline. Everything below the header must match exactly.
#[test]
fn binary_writes_samtools_stats_data_lines_byte_for_byte() {
    let got =
        std::fs::read_to_string(produced("samtools", &format!("{SAMPLE}.stats.txt"))).unwrap();
    let want = std::fs::read_to_string(fixture("test.stats.txt")).unwrap();
    let strip = |s: &str| {
        s.lines()
            .filter(|l| !l.starts_with('#'))
            .collect::<Vec<_>>()
            .join("\n")
    };
    assert_same_lines(&strip(&got), &strip(&want), "samtools stats data lines");
}

#[test]
fn the_stats_header_does_not_claim_the_wrong_subcommand() {
    let got =
        std::fs::read_to_string(produced("samtools", &format!("{SAMPLE}.stats.txt"))).unwrap();
    assert!(
        !got.contains("rustqc rna"),
        "the dna pipeline must not label its output as rna output"
    );
}

#[test]
fn binary_refuses_input_without_duplicate_marks() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let outdir = std::env::temp_dir().join("rustqc-dna-nodup");
    let _ = std::fs::remove_dir_all(&outdir);
    let output = std::process::Command::new(env!("CARGO_BIN_EXE_rustqc"))
        .arg("dna")
        .arg(root.join("tests/data/test_nodup.bam"))
        .arg("--outdir")
        .arg(&outdir)
        .arg("--json-summary")
        .arg("-")
        .output()
        .expect("failed to run the rustqc binary");
    let combined = format!(
        "{}{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        combined.contains("duplicate-flagged") || combined.contains("failed"),
        "expected a duplicate-marking complaint, got: {combined}"
    );
}

/// The JSON summary is the machine-readable face of a run, so its DNA block is
/// pinned against the same figures the mosdepth fixtures carry.
#[test]
fn json_summary_carries_the_dna_block() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let outdir = std::env::temp_dir().join("rustqc-dna-json");
    let _ = std::fs::remove_dir_all(&outdir);
    std::fs::create_dir_all(&outdir).unwrap();
    let json_path = outdir.join("summary.json");

    let status = std::process::Command::new(env!("CARGO_BIN_EXE_rustqc"))
        .arg("dna")
        .arg(root.join("tests/data/dna/test.dna.bam"))
        .arg("--outdir")
        .arg(&outdir)
        .arg("--window-size")
        .arg(WINDOW_SIZE.to_string())
        .arg("--json-summary")
        .arg(&json_path)
        .arg("--quiet")
        .status()
        .expect("failed to run the rustqc binary");
    assert!(status.success());

    let text = std::fs::read_to_string(&json_path).unwrap();
    // Checked as text rather than parsed: the point is that these exact
    // figures reach the summary, and pulling in a JSON parser for one test
    // would not make the assertion any stronger.
    for needle in [
        "\"genome_length\": 40001",
        "\"covered_bases\": 247878",
        "\"max_coverage\": 867",
        "\"total_reads\": 5644",
        "\"duplicates\": 1656",
    ] {
        assert!(
            text.contains(needle),
            "summary is missing {needle}:\n{text}"
        );
    }
    assert!(
        !text.contains("\"dupradar\""),
        "a dna run must not emit the rna summary blocks"
    );
}

/// The citations file names the tools this pipeline actually replicated.
#[test]
fn citations_name_the_dna_tools_only() {
    let citations = std::fs::read_to_string(run_binary().join("CITATIONS.md")).unwrap();
    assert!(citations.contains("mosdepth"), "mosdepth must be cited");
    assert!(citations.contains("Samtools"), "samtools must be cited");
    assert!(
        !citations.contains("dupRadar") && !citations.contains("RSeQC"),
        "a dna run must not cite the rna-only tools"
    );
}

/// The `.csi` companion indexes are not compared byte for byte: an index is
/// binary metadata over the compressed blocks, and two writers answering the
/// same queries need not produce the same bytes. What matters is that a region
/// query returns the same rows through our index as through mosdepth's.
///
/// The query goes through the `tabix` binary rather than rust-htslib's tabix
/// reader, which ends a fetched region by yielding a `TabixTruncatedRecord`
/// instead of stopping, and does so at different points for the two files. The
/// test is skipped where `tabix` is not installed, the same way the fixtures
/// themselves depend on the upstream tools being present.
#[test]
fn csi_indexes_answer_region_queries_like_mosdepths() {
    if std::process::Command::new("tabix")
        .arg("--version")
        .output()
        .is_err()
    {
        eprintln!("skipping: tabix is not installed");
        return;
    }

    let query = |path: &Path| -> String {
        let out = std::process::Command::new("tabix")
            .arg(path)
            .arg("chr22:2000-2500")
            .output()
            .unwrap_or_else(|e| panic!("querying {}: {e}", path.display()));
        assert!(
            out.status.success(),
            "tabix failed on {}: {}",
            path.display(),
            String::from_utf8_lossy(&out.stderr)
        );
        String::from_utf8(out.stdout).unwrap()
    };

    for suffix in ["per-base.bed.gz", "regions.bed.gz", "thresholds.bed.gz"] {
        let ours = produced("mosdepth", &format!("{SAMPLE}.{suffix}"));
        let index = ours.with_file_name(format!("{SAMPLE}.{suffix}.csi"));
        assert!(
            index.exists(),
            "{suffix} must have a .csi companion at {}",
            index.display()
        );

        let mine = query(&ours);
        let theirs = query(&fixture(&format!("test.{suffix}")));
        assert!(!mine.is_empty(), "{suffix}: the query returned nothing");
        assert_eq!(mine, theirs, "{suffix}: region query results differ");
    }
}

// ===================================================================
// Picard CollectInsertSizeMetrics
// ===================================================================

/// Drive the insert size accumulator over the whole test BAM.
fn insert_size_result() -> insert_size::InsertSizeResult {
    let bam_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/dna/test.dna.bam");
    let mut reader = bam::Reader::from_path(&bam_path).unwrap();
    let mut accum = InsertSizeAccum::new();
    let mut record = bam::Record::new();
    while let Some(result) = reader.read(&mut record) {
        result.unwrap();
        accum.process_read(&record);
    }
    accum.into_result(insert_size::DEFAULT_DEVIATIONS)
}

#[test]
fn insert_size_metrics_match_picard() {
    let path = scratch("test.insert_size_metrics.txt");
    insert_size::write_insert_size_metrics(&insert_size_result(), &path).unwrap();
    assert_same_lines(
        &std::fs::read_to_string(&path).unwrap(),
        &std::fs::read_to_string(fixture("test.insert_size_metrics.txt")).unwrap(),
        "insert size metrics",
    );
}

/// The headline figures, pinned separately so a failure in the metrics row is
/// easy to tell apart from a failure in the histogram below it.
#[test]
fn insert_size_headline_figures_match_picard() {
    let result = insert_size_result();
    let fr = result
        .rows
        .iter()
        .find(|r| r.orientation == insert_size::PairOrientation::Fr)
        .expect("the fixture library is FR");
    assert_eq!(fr.read_pairs, 1992, "read pairs");
    assert_eq!(fr.median, 122, "median insert size");
    assert_eq!(fr.mode, 96, "mode");
    assert_eq!(fr.median_absolute_deviation, 23, "MAD");
    assert_eq!(fr.min, 32, "minimum");
    assert_eq!(fr.max, 300, "maximum");
    assert!((fr.mean - 124.442269).abs() < 1e-6, "mean was {}", fr.mean);
    assert!(
        (fr.standard_deviation - 32.720214).abs() < 1e-6,
        "standard deviation was {}",
        fr.standard_deviation
    );
    assert_eq!(
        fr.widths,
        vec![9, 19, 27, 37, 47, 57, 69, 83, 103, 127, 181],
        "the eleven percentile widths"
    );
}

// ===================================================================
// Picard CollectWgsMetrics
// ===================================================================

fn wgs_result() -> WgsMetricsResult {
    let bam_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/dna/test.dna.bam");
    let mut reader = bam::Reader::from_path(&bam_path).unwrap();
    let header = reader.header().to_owned();
    let length = header.target_len(0).unwrap();

    let mut accum = WgsAccum::new(
        length,
        wgs_metrics::DEFAULT_MIN_MAPPING_QUALITY,
        wgs_metrics::DEFAULT_MIN_BASE_QUALITY,
    );
    let mut record = bam::Record::new();
    while let Some(result) = reader.read(&mut record) {
        result.unwrap();
        accum.process_read(&record);
    }
    let (counters, depths) = accum.into_parts();
    // The fixture reference carries no N bases, so the territory is its length.
    WgsMetricsResult::new(&depths, counters, length, wgs_metrics::DEFAULT_COVERAGE_CAP)
}

/// The exclusion breakdown is the heart of this tool: it is what separates
/// Picard's coverage from a plain depth count, and each fraction is a
/// different rule. They are pinned individually so a failure names the rule
/// that broke.
#[test]
fn wgs_exclusion_fractions_match_picard() {
    let result = wgs_result();
    assert_eq!(
        result.counters.total_aligned_bases, 670_989,
        "the denominator is every reference-aligned base of every primary mapped record"
    );
    let [dupe, mapq, unpaired, baseq, overlap, capped, total] = result.exclusion_fractions();
    let close = |got: f64, want: f64, what: &str| {
        assert!((got - want).abs() < 1e-6, "{what}: got {got}, want {want}");
    };
    close(dupe, 0.299737, "PCT_EXC_DUPE");
    close(mapq, 0.0, "PCT_EXC_MAPQ");
    close(unpaired, 0.0, "PCT_EXC_UNPAIRED");
    close(baseq, 0.007352, "PCT_EXC_BASEQ");
    close(overlap, 0.324694, "PCT_EXC_OVERLAP");
    close(capped, 0.157699, "PCT_EXC_CAPPED");
    close(total, 0.789481, "PCT_EXC_TOTAL");
}

#[test]
fn wgs_headline_figures_match_picard() {
    let result = wgs_result();
    assert_eq!(result.genome_territory, 40_001);
    assert_eq!(result.median_coverage, 0);
    assert_eq!(result.mad_coverage, 0);
    assert!(
        (result.mean_coverage - 3.531312).abs() < 1e-6,
        "mean was {}",
        result.mean_coverage
    );
    assert!(
        (result.sd_coverage - 27.339314).abs() < 1e-6,
        "standard deviation was {}",
        result.sd_coverage
    );
    let f = result.coverage_fractions();
    assert!((f[0] - 0.029124).abs() < 1e-6, "PCT_1X was {}", f[0]);
    assert!((f[13] - 0.01505).abs() < 1e-6, "PCT_100X was {}", f[13]);
}

/// The whole file, except the five columns RustQC does not compute.
///
/// `FOLD_80/90/95_BASE_PENALTY` are `?` in the fixture too, because Picard
/// could not compute them on this data. `HET_SNP_SENSITIVITY` and `HET_SNP_Q`
/// come from a Monte Carlo simulation that is out of scope, so RustQC writes
/// `?` where Picard writes a sampled value. Those two positions are the only
/// permitted difference.
#[test]
fn wgs_metrics_file_matches_picard_except_the_simulated_columns() {
    let path = scratch("test.wgs_metrics.txt");
    wgs_metrics::write_wgs_metrics(&wgs_result(), &path).unwrap();
    let got = std::fs::read_to_string(&path).unwrap();
    let want = std::fs::read_to_string(fixture("test.wgs_metrics.txt")).unwrap();

    let got_lines: Vec<&str> = got.lines().collect();
    let want_lines: Vec<&str> = want.lines().collect();
    assert_eq!(
        got_lines.len(),
        want_lines.len(),
        "line count differs: {} versus {}",
        got_lines.len(),
        want_lines.len()
    );

    for (i, (a, b)) in got_lines.iter().zip(want_lines.iter()).enumerate() {
        if i == 2 {
            // The metrics row: compare every column but the last two.
            let ours: Vec<&str> = a.split('\t').collect();
            let theirs: Vec<&str> = b.split('\t').collect();
            assert_eq!(ours.len(), theirs.len(), "column count differs");
            let simulated = ours.len() - 2;
            for (col, (x, y)) in ours.iter().zip(theirs.iter()).enumerate() {
                if col >= simulated {
                    continue;
                }
                assert_eq!(x, y, "column {col} of the metrics row differs");
            }
            assert_eq!(
                &ours[simulated..],
                &["?", "?"],
                "the simulated columns must be written as ?"
            );
        } else {
            assert_eq!(a, b, "line {} differs", i + 1);
        }
    }
}

#[test]
fn binary_writes_insert_size_metrics_byte_for_byte() {
    let got = std::fs::read_to_string(produced(
        "picard/insert_size",
        &format!("{SAMPLE}.insert_size_metrics.txt"),
    ))
    .unwrap();
    let want = std::fs::read_to_string(fixture("test.insert_size_metrics.txt")).unwrap();
    assert_same_lines(&got, &want, "insert size metrics from the binary");
}

/// As with the library-level check, the two Monte Carlo columns are the only
/// permitted difference.
#[test]
fn binary_writes_wgs_metrics_bar_the_simulated_columns() {
    let got = std::fs::read_to_string(produced(
        "picard/wgs_metrics",
        &format!("{SAMPLE}.wgs_metrics.txt"),
    ))
    .unwrap();
    let want = std::fs::read_to_string(fixture("test.wgs_metrics.txt")).unwrap();

    let got_lines: Vec<&str> = got.lines().collect();
    let want_lines: Vec<&str> = want.lines().collect();
    assert_eq!(got_lines.len(), want_lines.len(), "line count differs");
    for (i, (a, b)) in got_lines.iter().zip(want_lines.iter()).enumerate() {
        if i == 2 {
            let ours: Vec<&str> = a.split('\t').collect();
            let theirs: Vec<&str> = b.split('\t').collect();
            let simulated = ours.len() - 2;
            assert_eq!(&ours[..simulated], &theirs[..simulated], "metrics row");
        } else {
            assert_eq!(a, b, "line {} differs", i + 1);
        }
    }
}

/// Without a reference there is no way to size the genome territory, so the
/// analysis is skipped rather than reported against a wrong denominator.
#[test]
fn wgs_metrics_are_skipped_without_a_reference() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let outdir = std::env::temp_dir().join("rustqc-dna-noref");
    let _ = std::fs::remove_dir_all(&outdir);
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_rustqc"))
        .arg("dna")
        .arg(root.join("tests/data/dna/test.dna.bam"))
        .arg("--outdir")
        .arg(&outdir)
        .arg("--quiet")
        .status()
        .unwrap();
    assert!(status.success());
    assert!(
        !outdir.join("picard/wgs_metrics").exists(),
        "no reference means no WGS metrics"
    );
    assert!(
        outdir.join("picard/insert_size").exists(),
        "insert size needs no reference and must still be written"
    );
}
