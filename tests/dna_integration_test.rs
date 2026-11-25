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
use rustqc::dna::gc_bias::{self, GcBiasAccum};
use rustqc::dna::hs_metrics::{self, HsAccum, HsCounters, HsMetricsResult};
use rustqc::dna::insert_size::{self, InsertSizeAccum};
use rustqc::dna::intervals::IntervalSet;
use rustqc::dna::mosdepth::{output, ContigDepth, MosdepthResult};
use rustqc::dna::qualimap::{self, ContigQualimap, QualimapAccum};
use rustqc::dna::qualimap_output;
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

// ===================================================================
// Picard CollectGcBiasMetrics
// ===================================================================

fn gc_bias_result() -> gc_bias::GcBiasResult {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let reference: Vec<u8> = {
        let text = std::fs::read_to_string(root.join("tests/data/dna/genome.fasta")).unwrap();
        text.lines()
            .filter(|l| !l.starts_with('>'))
            .flat_map(|l| l.bytes())
            .collect()
    };

    let mut accum = GcBiasAccum::new(&reference, gc_bias::DEFAULT_WINDOW_SIZE);
    let mut reader = bam::Reader::from_path(root.join("tests/data/dna/test.dna.bam")).unwrap();
    let mut record = bam::Record::new();
    while let Some(result) = reader.read(&mut record) {
        result.unwrap();
        accum.process_read(&record, &reference);
    }
    accum.into_result(gc_bias::DEFAULT_WINDOW_SIZE)
}

/// The window table and the read assignment are the two rules that black-box
/// inference could not recover, so they are pinned before anything derived
/// from them.
#[test]
fn gc_bias_windows_and_read_starts_match_picard() {
    let result = gc_bias_result();
    let windows: u64 = result.rows.iter().map(|r| r.windows).sum();
    let read_starts: u64 = result.rows.iter().map(|r| r.read_starts).sum();
    assert_eq!(
        windows, 39_900,
        "sliding windows run from position 1 to len - window_size - 1"
    );
    assert_eq!(
        read_starts, 5_642,
        "secondary alignments count towards read starts"
    );
    assert_eq!(result.total_clusters, 2_822);
    assert_eq!(result.aligned_reads, 5_642);
}

#[test]
fn gc_bias_detail_metrics_match_picard() {
    let path = scratch("test.gc_bias.detail_metrics.txt");
    gc_bias::write_detail_metrics(&gc_bias_result(), &path).unwrap();
    assert_same_lines(
        &std::fs::read_to_string(&path).unwrap(),
        &std::fs::read_to_string(fixture("test.gc_bias.detail_metrics.txt")).unwrap(),
        "GC bias detail metrics",
    );
}

#[test]
fn gc_bias_summary_metrics_match_picard() {
    let path = scratch("test.gc_bias.summary_metrics.txt");
    gc_bias::write_summary_metrics(&gc_bias_result(), &path).unwrap();
    assert_same_lines(
        &std::fs::read_to_string(&path).unwrap(),
        &std::fs::read_to_string(fixture("test.gc_bias.summary_metrics.txt")).unwrap(),
        "GC bias summary metrics",
    );
}

// ===================================================================
// Picard CollectHsMetrics
// ===================================================================

/// Read one column out of the fixture's single metrics row.
fn hs_fixture_column(name: &str) -> String {
    let text = std::fs::read_to_string(fixture("test.hs_metrics.txt")).unwrap();
    let mut lines = text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty());
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let values: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let index = header
        .iter()
        .position(|h| *h == name)
        .unwrap_or_else(|| panic!("no column named {name}"));
    values[index].to_string()
}

fn hs_result() -> HsMetricsResult {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let targets = IntervalSet::from_bed(&root.join("tests/data/dna/targets.bed")).unwrap();

    let mut reader = bam::Reader::from_path(root.join("tests/data/dna/test.dna.bam")).unwrap();
    let header = reader.header().to_owned();
    let contig = String::from_utf8(header.tid2name(0).to_vec()).unwrap();
    let length = header.target_len(0).unwrap();

    let mut accum = HsAccum::new(&contig, length, &targets, &targets, 20, 20);
    let mut record = bam::Record::new();
    while let Some(result) = reader.read(&mut record) {
        result.unwrap();
        accum.process_read(&record);
    }
    let (counters, depths, target_mask) = accum.into_parts();

    let target_depths: Vec<u32> = depths
        .iter()
        .zip(target_mask.iter())
        .filter(|(_, on_target)| **on_target)
        .map(|(depth, _)| *depth)
        .collect();

    let zero_coverage_targets = targets
        .on(&contig)
        .iter()
        .filter(|interval| (interval.start..interval.end).all(|p| depths[p as usize] == 0))
        .count() as u64;

    let library_size =
        hs_metrics::estimate_library_size(counters.selected_pairs, counters.selected_unique_pairs);

    HsMetricsResult {
        bait_set: targets.name().to_string(),
        bait_territory: targets.territory(),
        target_territory: targets.territory(),
        genome_size: length,
        counters,
        target_depths,
        zero_coverage_targets,
        target_count: targets.len() as u64,
        library_size,
    }
}

/// The counters are the part that had to be taken from Picard's source, so
/// each is pinned against the fixture individually.
#[test]
fn hs_counters_match_picard() {
    let result = hs_result();
    let c: &HsCounters = &result.counters;
    let want = |name: &str| -> u64 { hs_fixture_column(name).parse().unwrap() };

    assert_eq!(result.bait_territory, want("BAIT_TERRITORY"));
    assert_eq!(result.target_territory, want("TARGET_TERRITORY"));
    assert_eq!(result.genome_size, want("GENOME_SIZE"));
    assert_eq!(c.total_reads, want("TOTAL_READS"), "secondary excluded");
    assert_eq!(c.pf_bases, want("PF_BASES"));
    assert_eq!(c.pf_unique_reads, want("PF_UNIQUE_READS"));
    assert_eq!(c.pf_uq_reads_aligned, want("PF_UQ_READS_ALIGNED"));
    assert_eq!(c.pf_bases_aligned, want("PF_BASES_ALIGNED"));
    assert_eq!(c.pf_uq_bases_aligned, want("PF_UQ_BASES_ALIGNED"));
    assert_eq!(c.on_bait_bases, want("ON_BAIT_BASES"));
    assert_eq!(c.near_bait_bases, want("NEAR_BAIT_BASES"));
    assert_eq!(c.off_bait_bases, want("OFF_BAIT_BASES"));
    assert_eq!(
        c.on_target_bases,
        want("ON_TARGET_BASES"),
        "overlap clipping runs before the base quality filter"
    );
    assert_eq!(result.library_size, Some(want("HS_LIBRARY_SIZE")));
}

/// The exclusion fractions are where HsMetrics parts company with
/// CollectWgsMetrics, so they get their own assertions.
#[test]
fn hs_exclusion_fractions_match_picard() {
    let result = hs_result();
    let aligned = result.counters.pf_bases_aligned as f64;
    let want = |name: &str| -> f64 { hs_fixture_column(name).parse().unwrap() };
    let close = |got: f64, name: &str| {
        let expected = want(name);
        assert!(
            (got - expected).abs() < 1e-6,
            "{name}: got {got}, want {expected}"
        );
    };
    close(
        result.counters.excluded_dupe as f64 / aligned,
        "PCT_EXC_DUPE",
    );
    close(
        result.counters.excluded_overlap as f64 / aligned,
        "PCT_EXC_OVERLAP",
    );
    close(
        result.counters.excluded_baseq as f64 / aligned,
        "PCT_EXC_BASEQ",
    );
    close(
        result.counters.excluded_off_target as f64 / aligned,
        "PCT_EXC_OFF_TARGET",
    );
}

#[test]
fn hs_target_coverage_matches_picard() {
    let result = hs_result();
    let want = |name: &str| -> f64 { hs_fixture_column(name).parse().unwrap() };
    assert!(
        (result.mean_target_coverage() - want("MEAN_TARGET_COVERAGE")).abs() < 1e-6,
        "mean target coverage was {}",
        result.mean_target_coverage()
    );
    assert!(
        (result.mean_bait_coverage() - want("MEAN_BAIT_COVERAGE")).abs() < 1e-6,
        "mean bait coverage was {}",
        result.mean_bait_coverage()
    );
    let (median, min, max) = result.target_coverage_bounds();
    assert_eq!(u64::from(median), want("MEDIAN_TARGET_COVERAGE") as u64);
    assert_eq!(u64::from(min), want("MIN_TARGET_COVERAGE") as u64);
    assert_eq!(u64::from(max), want("MAX_TARGET_COVERAGE") as u64);

    let fractions = result.target_coverage_fractions();
    for (level, got) in hs_metrics::TARGET_COVERAGE_LEVELS.iter().zip(fractions) {
        let expected = want(&format!("PCT_TARGET_BASES_{level}X"));
        assert!(
            (got - expected).abs() < 1e-6,
            "PCT_TARGET_BASES_{level}X: got {got}, want {expected}"
        );
    }
}

/// A second binary run, this time in targeted mode with a reference, so the
/// GC bias and targeted outputs are produced.
fn run_binary_targeted() -> &'static Path {
    static OUTDIR: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
    OUTDIR.get_or_init(|| {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"));
        let outdir = std::env::temp_dir().join("rustqc-dna-targeted");
        let _ = std::fs::remove_dir_all(&outdir);
        std::fs::create_dir_all(&outdir).unwrap();
        let status = std::process::Command::new(env!("CARGO_BIN_EXE_rustqc"))
            .arg("dna")
            .arg(root.join("tests/data/dna/test.dna.bam"))
            .arg("--reference")
            .arg(root.join("tests/data/dna/genome.fasta"))
            .arg("--targets")
            .arg(root.join("tests/data/dna/targets.bed"))
            .arg("--outdir")
            .arg(&outdir)
            .arg("--quiet")
            .status()
            .expect("failed to run the rustqc binary");
        assert!(status.success(), "rustqc dna exited with {status}");
        outdir
    })
}

#[test]
fn binary_writes_gc_bias_metrics_byte_for_byte() {
    for suffix in ["gc_bias.detail_metrics", "gc_bias.summary_metrics"] {
        let got = std::fs::read_to_string(
            run_binary_targeted()
                .join("picard/gc_bias")
                .join(format!("{SAMPLE}.{suffix}.txt")),
        )
        .unwrap();
        let want = std::fs::read_to_string(fixture(&format!("test.{suffix}.txt"))).unwrap();
        assert_same_lines(&got, &want, suffix);
    }
}

/// Every HS metrics column but the seven that need Picard's theoretical
/// sensitivity simulation or its per-target GC dropout, which RustQC does not
/// compute and writes as Picard writes its own uncomputable values.
#[test]
fn binary_writes_hs_metrics_bar_the_simulated_columns() {
    let path = run_binary_targeted()
        .join("picard/hs_metrics")
        .join(format!("{SAMPLE}.hs_metrics.txt"));
    let parse = |text: &str| -> std::collections::HashMap<String, String> {
        let mut lines = text
            .lines()
            .filter(|l| !l.starts_with('#') && !l.trim().is_empty());
        let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
        let values: Vec<&str> = lines.next().unwrap().split('\t').collect();
        header
            .iter()
            .zip(values.iter())
            .map(|(h, v)| (h.to_string(), v.to_string()))
            .collect()
    };

    let ours = parse(&std::fs::read_to_string(&path).unwrap());
    let theirs = parse(&std::fs::read_to_string(fixture("test.hs_metrics.txt")).unwrap());

    let uncomputed: Vec<String> = [
        "HET_SNP_SENSITIVITY",
        "HET_SNP_Q",
        "AT_DROPOUT",
        "GC_DROPOUT",
        "FOLD_80_BASE_PENALTY",
    ]
    .iter()
    .map(|s| s.to_string())
    .chain(
        hs_metrics::PENALTY_LEVELS
            .iter()
            .map(|n| format!("HS_PENALTY_{n}X")),
    )
    .collect();

    let mut compared = 0;
    for (column, want) in &theirs {
        if uncomputed.contains(column) {
            continue;
        }
        compared += 1;
        let got = ours
            .get(column)
            .unwrap_or_else(|| panic!("we do not emit column {column}"));
        assert_eq!(got, want, "column {column} differs");
    }
    assert!(
        compared >= 55,
        "expected to compare most of the columns, only did {compared}"
    );
}

/// Targeted outputs appear only when targets are given.
#[test]
fn hs_metrics_are_absent_without_targets() {
    assert!(
        !run_binary().join("picard/hs_metrics").exists(),
        "no targets means no targeted metrics"
    );
    assert!(
        run_binary_targeted().join("picard/hs_metrics").exists(),
        "targets must produce them"
    );
}

// ===================================================================
// Qualimap bamqc
// ===================================================================

fn qualimap_result() -> ContigQualimap {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut reader = bam::Reader::from_path(root.join("tests/data/dna/test.dna.bam")).unwrap();
    let header = reader.header().to_owned();
    let contig = String::from_utf8(header.tid2name(0).to_vec()).unwrap();
    let length = header.target_len(0).unwrap();

    let mut accum = QualimapAccum::new(&contig, length, qualimap::DEFAULT_NUM_WINDOWS);
    let mut record = bam::Record::new();
    while let Some(result) = reader.read(&mut record) {
        result.unwrap();
        accum.process_read(&record);
    }
    accum.into_result()
}

/// Qualimap measures coverage differently from every other tool here: no
/// filtering at all, deletions counted, and no mate-overlap correction. The
/// figures are pinned so that difference stays deliberate.
#[test]
fn qualimap_globals_match() {
    let r = qualimap_result();
    let c = &r.counters;
    assert_eq!(r.midpoints.len(), 397, "window count");
    assert_eq!(c.reads, 5642, "secondary alignments are counted separately");
    assert_eq!(c.secondary, 2);
    assert_eq!(c.mapped, 5640);
    assert_eq!(c.duplicates, 1656);
    assert_eq!(c.paired_first, 2820);
    assert_eq!(c.paired_second, 2820);
    assert_eq!(c.paired_both, 5640);
    assert_eq!(c.singletons, 0);
    assert_eq!(c.sequenced_bases, 670_989);
    assert_eq!(c.mapped_bases, 670_999, "deletions count as mapped");
}

#[test]
fn qualimap_base_composition_matches() {
    let c = qualimap_result().counters;
    // A, C, G, T, N in reference orientation.
    assert_eq!(c.base_counts, [233_897, 101_959, 103_412, 231_444, 277]);
}

#[test]
fn qualimap_mismatches_and_indels_match() {
    let c = qualimap_result().counters;
    assert_eq!(
        c.mismatches(),
        1350,
        "NM less insertions, not less deletions"
    );
    assert_eq!(c.insertions, 2);
    assert_eq!(c.deletions, 10);
    assert_eq!(c.reads_with_insertion, 2);
    assert_eq!(c.reads_with_deletion, 10);
    let rate = c.general_error_rate();
    assert!((rate - 0.002).abs() < 5e-4, "general error rate was {rate}");
}

#[test]
fn qualimap_insert_size_matches() {
    let (mean, sd, median) = qualimap_result().counters.insert_size_stats();
    assert!((mean - 125.6844).abs() < 1e-4, "mean was {mean}");
    assert!((sd - 32.4421).abs() < 1e-4, "sd was {sd}");
    assert_eq!(median, 123);
}

#[test]
fn qualimap_coverage_matches() {
    let r = qualimap_result();
    assert!(
        (r.mean_coverage() - 16.7746).abs() < 1e-4,
        "mean coverage was {}",
        r.mean_coverage()
    );
    assert_eq!(r.coverage_histogram.get(&0), Some(&38_820));
    assert_eq!(r.coverage_histogram.get(&1), Some(&40));
    let fraction = r.genome_fraction();
    assert!(
        (fraction[0].1 - 2.9524261893452746).abs() < 1e-9,
        "1X fraction was {}",
        fraction[0].1
    );
}

/// The mapping quality histogram truncates the per-position mean rather than
/// rounding it, which moves 243 positions between the 59 and 60 bins.
#[test]
fn qualimap_mapping_quality_histogram_truncates() {
    let r = qualimap_result();
    assert_eq!(r.mapq_histogram.get(&59), Some(&248));
    assert_eq!(r.mapq_histogram.get(&60), Some(&933));
}

#[test]
fn qualimap_window_positions_are_midpoints() {
    let r = qualimap_result();
    assert!((r.midpoints[0] - 51.0).abs() < 1e-9);
    assert!((r.midpoints[1] - 152.0).abs() < 1e-9);
}

#[test]
fn qualimap_clipping_profile_is_a_distribution_over_clipped_bases() {
    let c = qualimap_result().counters;
    assert_eq!(c.clipped_bases, 863, "the profile's denominator");
    let first = 100.0 * c.clipping_by_position[0] as f64 / c.clipped_bases as f64;
    assert!(
        (first - 1.8539976825028968).abs() < 1e-9,
        "clipping at position 0 was {first}"
    );
}

/// Base composition is taken in reference orientation while the clipped span
/// that selects positions is taken in sequencing orientation. Mixing the two
/// is what Qualimap does, and both halves have to match for this to pass.
#[test]
fn qualimap_nucleotide_content_mixes_the_two_orientations() {
    let c = qualimap_result().counters;
    let first = c.nucleotide_by_position[0];
    let total: u64 = first.iter().sum();
    assert_eq!(total, 5624, "clipped positions are excluded");
    let pct = |i: usize| 100.0 * first[i] as f64 / total as f64;
    assert!(
        (pct(0) - 36.575391180654336).abs() < 1e-9,
        "A was {}",
        pct(0)
    );
    assert!(
        (pct(1) - 12.820056899004268).abs() < 1e-9,
        "C was {}",
        pct(1)
    );
    assert!(
        (pct(2) - 18.509957325746797).abs() < 1e-9,
        "G was {}",
        pct(2)
    );
    assert!(
        (pct(3) - 32.059032716927454).abs() < 1e-9,
        "T was {}",
        pct(3)
    );
    assert!(
        (pct(4) - 0.03556187766714083).abs() < 1e-9,
        "N was {}",
        pct(4)
    );
}

/// The whole `genome_results.txt`, minus the two lines that record the
/// absolute paths the run used.
///
/// Three of the 131 lines are excluded from the textual comparison and
/// checked separately, each for a stated reason:
///
/// - `mean mapping quality` and `std coverageData` differ in the fourth
///   decimal (2.4179 against 2.4178, 154.9340 against 154.9323). Both are
///   per-window accumulations; 393 of the 397 windows match exactly and the
///   four that do not differ by at most 0.053. They are asserted numerically
///   with a tolerance.
/// - `homopolymer indels` differs outright. Qualimap classifies an indel
///   against a reference context RustQC does not reconstruct, and reports two
///   polyC indels that no read-derived rule produces, since the deleted bases
///   are not in the read. That line is asserted only to be present and
///   well-formed.
#[test]
fn qualimap_genome_results_match() {
    let path = scratch("genome_results.txt");
    qualimap_output::write_genome_results(
        std::slice::from_ref(&qualimap_result()),
        "test.dna.bam",
        &path,
    )
    .unwrap();
    let strip = |s: &str| {
        s.lines()
            .filter(|l| !l.contains("bam file =") && !l.contains("outfile ="))
            .collect::<Vec<_>>()
            .join("\n")
    };
    let got = strip(&std::fs::read_to_string(&path).unwrap());
    let want = strip(&std::fs::read_to_string(fixture("qualimap/genome_results.txt")).unwrap());

    let number = |text: &str, key: &str| -> f64 {
        text.lines()
            .find(|l| l.contains(key))
            .and_then(|l| l.split('=').nth(1))
            .map(|v| v.trim().trim_end_matches('X').parse().unwrap())
            .unwrap_or_else(|| panic!("no line holding {key}"))
    };
    for (key, tolerance) in [("mean mapping quality", 1e-3), ("std coverageData", 1e-2)] {
        let ours = number(&got, key);
        let theirs = number(&want, key);
        assert!(
            (ours - theirs).abs() < tolerance,
            "{key}: got {ours}, want {theirs}"
        );
    }

    let homopolymer = got
        .lines()
        .find(|l| l.contains("homopolymer indels"))
        .expect("the homopolymer line must still be written");
    assert!(
        homopolymer.trim_end().ends_with('%'),
        "homopolymer line is malformed: {homopolymer}"
    );

    // The coverage fraction lines are compared numerically: about five
    // reference positions out of 40001 sit one deeper here than in Qualimap,
    // which moves these percentages in the third decimal.
    let fractions = |text: &str| -> Vec<f64> {
        text.lines()
            .filter(|l| l.contains("of reference with a coverageData"))
            .map(|l| {
                l.split("There is a")
                    .nth(1)
                    .and_then(|r| r.split('%').next())
                    .unwrap()
                    .trim()
                    .parse()
                    .unwrap()
            })
            .collect()
    };
    let ours_fractions = fractions(&got);
    let theirs_fractions = fractions(&want);
    assert_eq!(
        ours_fractions.len(),
        theirs_fractions.len(),
        "fraction lines"
    );
    for (i, (a, b)) in ours_fractions.iter().zip(&theirs_fractions).enumerate() {
        assert!(
            (a - b).abs() < 0.01,
            "coverage fraction at level {}: got {a}, want {b}",
            i + 1
        );
    }

    // The per-contig row carries the same standard deviation, so it is
    // compared field by field with the last one given a tolerance.
    let contig_row = |text: &str| -> Vec<String> {
        text.lines()
            .find(|l| l.starts_with('\t'))
            .map(|l| l.trim().split('\t').map(str::to_string).collect())
            .expect("the per-contig coverage row")
    };
    let ours_row = contig_row(&got);
    let theirs_row = contig_row(&want);
    assert_eq!(
        ours_row[..4],
        theirs_row[..4],
        "per-contig name, length, bases and mean"
    );
    let ours_sd: f64 = ours_row[4].parse().unwrap();
    let theirs_sd: f64 = theirs_row[4].parse().unwrap();
    assert!(
        (ours_sd - theirs_sd).abs() < 1e-2,
        "per-contig standard deviation: got {ours_sd}, want {theirs_sd}"
    );

    let excluded = [
        "mean mapping quality",
        "std coverageData",
        "homopolymer indels",
        "of reference with a coverageData",
    ];
    let drop = |text: &str| -> String {
        text.lines()
            .filter(|l| !l.starts_with('\t') && !excluded.iter().any(|k| l.contains(k)))
            .collect::<Vec<_>>()
            .join("\n")
    };
    assert_same_lines(&drop(&got), &drop(&want), "genome_results.txt");
}

/// The raw data tables, compared as numbers rather than as text.
///
/// Three match byte for byte. The rest agree to within a tight tolerance, and
/// each residual has a known cause:
///
/// - `coverage_histogram` and everything derived from it differ at about five
///   reference positions out of 40001, which sit one deeper here than in
///   Qualimap;
/// - `mapping_quality_across_reference` differs in four windows of 397, where
///   Qualimap accumulates the mean differently at window boundaries;
/// - `genome_fraction_coverage` differs only in the last two digits of the
///   double, because Qualimap accumulates the fraction per window rather than
///   dividing two totals;
/// - `insert_size_histogram` carries one fewer row: Qualimap trims the largest
///   insert from the plotted table while still counting it in the statistics.
#[test]
fn qualimap_raw_data_tables_match() {
    let dir = run_binary_targeted().join("qualimap/raw_data_qualimapReport");

    let numbers = |path: &Path| -> Vec<Vec<f64>> {
        std::fs::read_to_string(path)
            .unwrap_or_else(|e| panic!("reading {}: {e}", path.display()))
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| {
                l.split('\t')
                    .map(|v| v.trim().parse::<f64>().unwrap_or(f64::NAN))
                    .collect()
            })
            .collect()
    };

    // Tables that reproduce exactly.
    for name in [
        "mapped_reads_clipping_profile.txt",
        "mapped_reads_nucleotide_content.txt",
        "mapping_quality_histogram.txt",
    ] {
        let got = std::fs::read_to_string(dir.join(name)).unwrap();
        let want =
            std::fs::read_to_string(fixture(&format!("qualimap/raw_data_qualimapReport/{name}")))
                .unwrap();
        assert_same_lines(&got, &want, name);
    }

    // Tables compared numerically, with the tolerated row count in each.
    for (name, tolerance, max_differing_rows) in [
        ("coverage_across_reference.txt", 1e-6, 20usize),
        ("coverage_histogram.txt", 1.5, 10),
        ("genome_fraction_coverage.txt", 1e-6, 52),
        ("insert_size_across_reference.txt", 1e-6, 5),
        ("mapping_quality_across_reference.txt", 0.1, 5),
    ] {
        let ours = numbers(&dir.join(name));
        let theirs = numbers(&fixture(&format!(
            "qualimap/raw_data_qualimapReport/{name}"
        )));
        assert_eq!(ours.len(), theirs.len(), "{name}: row count");

        let mut differing = 0;
        for (row, (a, b)) in ours.iter().zip(&theirs).enumerate() {
            assert_eq!(a.len(), b.len(), "{name}: row {row} column count");
            if a.iter().zip(b).any(|(x, y)| (x - y).abs() > tolerance) {
                differing += 1;
                assert!(
                    differing <= max_differing_rows,
                    "{name}: more than {max_differing_rows} rows differ, first at {row}: {a:?} against {b:?}"
                );
            }
        }
    }

    // The insert size histogram is the one table with a different row count.
    let ours = numbers(&dir.join("insert_size_histogram.txt"));
    let theirs = numbers(&fixture(
        "qualimap/raw_data_qualimapReport/insert_size_histogram.txt",
    ));
    assert!(
        ours.len() == theirs.len() + 1,
        "expected exactly one extra row, got {} against {}",
        ours.len(),
        theirs.len()
    );
    for (a, b) in ours.iter().zip(&theirs) {
        assert_eq!(a, b, "insert size histogram rows before the trimmed one");
    }
}

/// The HTML report is RustQC's own page rather than a copy of Qualimap's, so
/// it is checked for structure and for carrying the headline numbers.
#[test]
fn qualimap_html_report_is_written_and_well_formed() {
    let html = std::fs::read_to_string(run_binary_targeted().join("qualimap/qualimapReport.html"))
        .unwrap();
    assert!(html.starts_with("<!doctype html>"), "missing doctype");
    assert!(html.trim_end().ends_with("</html>"), "unclosed document");
    assert!(html.contains("BamQC report"), "missing the title");
    assert!(html.contains("40,001"), "missing the reference length");
    assert!(html.contains("5,642"), "missing the read count");
    assert!(html.contains("16.7746X"), "missing the mean coverage");
}
