//! CLI regression for enabled preseq analyses that cannot estimate complexity.
//! A skipped analysis must leave an explicit status, without a numeric curve.

use rust_htslib::bam;
use rust_htslib::bam::header::HeaderRecord;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

struct Fixture {
    root: PathBuf,
    bam: PathBuf,
    gtf: PathBuf,
}

impl Fixture {
    fn new(multiplicities: &[usize]) -> Self {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "rustqc-preseq-output-{}-{nonce}",
            std::process::id()
        ));
        fs::create_dir_all(&root).unwrap();
        let bam_path = root.join("reads.bam");
        let gtf = root.join("genes.gtf");
        fs::write(
            &gtf,
            "chr1\ttest\tgene\t1\t100000\t.\t+\t.\tgene_id \"g\"; gene_name \"g\";\n\
             chr1\ttest\texon\t1\t100000\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\"; gene_name \"g\";\n",
        )
        .unwrap();
        let mut header = bam::Header::new();
        header.push_record(
            HeaderRecord::new(b"HD")
                .push_tag(b"VN", "1.6")
                .push_tag(b"SO", "coordinate"),
        );
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 100000),
        );
        let view = bam::HeaderView::from_header(&header);
        let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
        for (molecule, &copies) in multiplicities.iter().enumerate() {
            for copy in 0..copies {
                let sam = format!(
                    "m{molecule}_{copy}\t0\tchr1\t{}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII",
                    molecule * 100 + 1
                );
                writer
                    .write(&bam::Record::from_sam(&view, sam.as_bytes()).unwrap())
                    .unwrap();
            }
        }
        drop(writer);
        bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();
        Self {
            root,
            bam: bam_path,
            gtf,
        }
    }

    fn run(&self, enabled: bool, flat: bool) -> Output {
        let config = self.root.join("config.yaml");
        let mut yaml = String::from(
            "rna:\n  dupradar:\n    dup_matrix: false\n    intercept_slope: false\n    density_scatter_plot: false\n    boxplot: false\n    expression_histogram: false\n    multiqc_intercept: false\n    multiqc_curve: false\n  featurecounts:\n    counts_file: false\n    summary_file: false\n    biotype_summary_file: false\n    biotype_counts: false\n    biotype_counts_mqc: false\n    biotype_rrna_mqc: false\n",
        );
        for tool in [
            "bam_stat",
            "infer_experiment",
            "read_duplication",
            "read_distribution",
            "junction_annotation",
            "junction_saturation",
            "inner_distance",
            "tin",
            "flagstat",
            "idxstats",
            "samtools_stats",
            "qualimap",
        ] {
            yaml.push_str(&format!("  {tool}:\n    enabled: false\n"));
        }
        yaml.push_str(&format!(
            "  preseq:\n    enabled: {enabled}\n    n_bootstraps: 0\n    max_extrap: 1000\n    step_size: 100\n    max_terms: 4\n"
        ));
        fs::write(&config, yaml).unwrap();
        let mut command = Command::new(env!("CARGO_BIN_EXE_rustqc"));
        command
            .arg("rna")
            .arg(&self.bam)
            .arg("--gtf")
            .arg(&self.gtf)
            .arg("--config")
            .arg(config)
            .arg("--outdir")
            .arg(self.root.join("out"))
            .args([
                "--sample-name",
                "sample",
                "--skip-dup-check",
                "--threads",
                "1",
            ])
            .env("XDG_CONFIG_HOME", self.root.join("config-home"));
        for variable in [
            "RUSTQC_SKIP_PRESEQ",
            "RUSTQC_PRESEQ_SEED",
            "RUSTQC_PRESEQ_MAX_EXTRAP",
            "RUSTQC_PRESEQ_STEP_SIZE",
            "RUSTQC_PRESEQ_N_BOOTSTRAPS",
            "RUSTQC_PRESEQ_SEG_LEN",
            "RUSTQC_PAIRED",
            "RUSTQC_FLAT_OUTPUT",
        ] {
            command.env_remove(variable);
        }
        if flat {
            command.arg("--flat-output");
        }
        command.output().unwrap()
    }

    fn output_dir(&self, flat: bool) -> PathBuf {
        if flat {
            self.root.join("out")
        } else {
            self.root.join("out/preseq")
        }
    }
}

fn messages(output: &Output) -> String {
    format!(
        "{}\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    )
}

fn assert_skipped(fixture: &Fixture, reason: &str, flat: bool) {
    let output = fixture.run(true, flat);
    assert!(output.status.success(), "{}", messages(&output));
    assert!(messages(&output).contains(reason), "{}", messages(&output));
    let directory = fixture.output_dir(flat);
    assert!(!directory.join("sample.lc_extrap.txt").exists());
    let status_path = directory.join("sample.preseq_skipped.txt");
    assert!(
        status_path.is_file(),
        "enabled but unestimable preseq must write a skipped-status artifact: {}",
        status_path.display()
    );
    let status = fs::read_to_string(status_path).unwrap();
    let rows: Vec<_> = status.lines().collect();
    assert_eq!(rows.len(), 2, "{status}");
    assert_eq!(rows[0], "STATUS\tREASON");
    assert!(rows[1].starts_with("skipped\t"), "{status}");
    assert!(rows[1].contains(reason), "{status}");
    assert!(!status.contains("EXPECTED_DISTINCT"));
    assert!(!status.contains("distinct_reads"));
}

#[test]
fn test_preseq_too_few_molecules_writes_status_without_curve() {
    let fixture = Fixture::new(&[1, 1]);
    assert_skipped(&fixture, "Too few distinct molecules (2)", false);
    fs::remove_dir_all(fixture.root).unwrap();
}

#[test]
fn test_preseq_negative_good_toulmin_writes_status_without_curve() {
    let fixture = Fixture::new(&[2; 30]);
    assert_skipped(
        &fixture,
        "Good-Toulmin 2x extrapolation check failed (value=-30.00)",
        false,
    );
    fs::remove_dir_all(fixture.root).unwrap();
}

#[test]
fn test_preseq_skipped_status_respects_flat_output() {
    let fixture = Fixture::new(&[1, 1]);
    assert_skipped(&fixture, "Too few distinct molecules (2)", true);
    assert!(!fixture.root.join("out/preseq").exists());
    fs::remove_dir_all(fixture.root).unwrap();
}

#[test]
fn test_preseq_disabled_does_not_write_curve_or_status() {
    let fixture = Fixture::new(&[1, 1]);
    let output = fixture.run(false, false);
    assert!(output.status.success(), "{}", messages(&output));
    assert!(!fixture.output_dir(false).exists());
    fs::remove_dir_all(fixture.root).unwrap();
}

#[test]
fn test_preseq_skipped_status_write_error_is_not_swallowed() {
    let fixture = Fixture::new(&[1, 1]);
    fs::create_dir_all(fixture.output_dir(false).join("sample.preseq_skipped.txt")).unwrap();
    let output = fixture.run(true, false);
    assert!(!output.status.success(), "{}", messages(&output));
    assert!(
        messages(&output).contains("Failed to create preseq skipped-status output"),
        "{}",
        messages(&output)
    );
    assert!(!fixture
        .output_dir(false)
        .join("sample.lc_extrap.txt")
        .exists());
    fs::remove_dir_all(fixture.root).unwrap();
}

fn estimable_fixture() -> Fixture {
    let multiplicities: Vec<_> = [(1, 100), (2, 50), (3, 25), (4, 12), (5, 6), (6, 3), (7, 1)]
        .into_iter()
        .flat_map(|(copies, molecules)| std::iter::repeat_n(copies, molecules))
        .collect();
    Fixture::new(&multiplicities)
}

fn assert_numeric_curve(directory: &Path) {
    let table = fs::read_to_string(directory.join("sample.lc_extrap.txt")).unwrap();
    assert_eq!(
        table.lines().next(),
        Some("TOTAL_READS\tEXPECTED_DISTINCT\tLOWER_0.95CI\tUPPER_0.95CI")
    );
    assert!(table.lines().count() > 2, "{table}");
    let mut positive_prediction = false;
    for line in table.lines().skip(1) {
        let fields: Vec<_> = line.split('\t').collect();
        assert_eq!(fields.len(), 4, "{line}");
        for value in &fields[..2] {
            assert!(value.parse::<f64>().unwrap().is_finite(), "{line}");
        }
        positive_prediction |= fields[1].parse::<f64>().unwrap() > 0.0;
        for value in &fields[2..] {
            assert!(
                *value == "-" || value.parse::<f64>().unwrap().is_finite(),
                "{line}"
            );
        }
    }
    assert!(positive_prediction, "{table}");
}

#[test]
fn test_preseq_estimable_sample_keeps_numeric_output() {
    let fixture = estimable_fixture();
    let output = fixture.run(true, false);
    assert!(output.status.success(), "{}", messages(&output));
    let directory = fixture.output_dir(false);
    assert_numeric_curve(&directory);
    assert!(!directory.join("sample.preseq_skipped.txt").exists());
    fs::remove_dir_all(fixture.root).unwrap();
}

fn replace_input(target: &Fixture, source: &Fixture) {
    fs::copy(&source.bam, &target.bam).unwrap();
    bam::index::build(&target.bam, None, bam::index::Type::Bai, 1).unwrap();
}

fn write_neighbors(directory: &Path) {
    fs::write(directory.join("other.lc_extrap.txt"), "other curve\n").unwrap();
    fs::write(directory.join("other.preseq_skipped.txt"), "other status\n").unwrap();
}

fn assert_neighbors(directory: &Path) {
    assert_eq!(
        fs::read_to_string(directory.join("other.lc_extrap.txt")).unwrap(),
        "other curve\n"
    );
    assert_eq!(
        fs::read_to_string(directory.join("other.preseq_skipped.txt")).unwrap(),
        "other status\n"
    );
}

#[test]
fn test_preseq_success_then_skip_removes_only_own_curve() {
    for flat in [false, true] {
        let fixture = estimable_fixture();
        let skipped = Fixture::new(&[1, 1]);
        let output = fixture.run(true, flat);
        assert!(output.status.success(), "{}", messages(&output));
        let directory = fixture.output_dir(flat);
        assert_numeric_curve(&directory);
        write_neighbors(&directory);
        replace_input(&fixture, &skipped);
        assert_skipped(&fixture, "Too few distinct molecules (2)", flat);
        assert_neighbors(&directory);
        fs::remove_dir_all(fixture.root).unwrap();
        fs::remove_dir_all(skipped.root).unwrap();
    }
}

#[test]
fn test_preseq_skip_then_success_removes_only_own_status() {
    for flat in [false, true] {
        let fixture = Fixture::new(&[1, 1]);
        let estimable = estimable_fixture();
        assert_skipped(&fixture, "Too few distinct molecules (2)", flat);
        let directory = fixture.output_dir(flat);
        write_neighbors(&directory);
        replace_input(&fixture, &estimable);
        let output = fixture.run(true, flat);
        assert!(output.status.success(), "{}", messages(&output));
        assert_numeric_curve(&directory);
        assert!(!directory.join("sample.preseq_skipped.txt").exists());
        assert_neighbors(&directory);
        fs::remove_dir_all(fixture.root).unwrap();
        fs::remove_dir_all(estimable.root).unwrap();
    }
}

#[test]
fn test_preseq_stale_curve_removal_error_is_not_swallowed() {
    for flat in [false, true] {
        let fixture = Fixture::new(&[1, 1]);
        let directory = fixture.output_dir(flat);
        let blocked = directory.join("sample.lc_extrap.txt");
        fs::create_dir_all(&blocked).unwrap();
        fs::write(blocked.join("keep"), "do not remove recursively\n").unwrap();
        let output = fixture.run(true, flat);
        assert!(!output.status.success(), "{}", messages(&output));
        assert!(
            messages(&output).contains("Failed to remove stale preseq output"),
            "{}",
            messages(&output)
        );
        assert!(directory.join("sample.preseq_skipped.txt").is_file());
        assert_eq!(
            fs::read_to_string(blocked.join("keep")).unwrap(),
            "do not remove recursively\n"
        );
        fs::remove_dir_all(fixture.root).unwrap();
    }
}

#[test]
fn test_preseq_stale_status_removal_error_is_not_swallowed() {
    for flat in [false, true] {
        let fixture = estimable_fixture();
        let directory = fixture.output_dir(flat);
        let blocked = directory.join("sample.preseq_skipped.txt");
        fs::create_dir_all(&blocked).unwrap();
        fs::write(blocked.join("keep"), "do not remove recursively\n").unwrap();
        let output = fixture.run(true, flat);
        assert!(!output.status.success(), "{}", messages(&output));
        assert!(
            messages(&output).contains("Failed to remove stale preseq output"),
            "{}",
            messages(&output)
        );
        assert_numeric_curve(&directory);
        assert_eq!(
            fs::read_to_string(blocked.join("keep")).unwrap(),
            "do not remove recursively\n"
        );
        fs::remove_dir_all(fixture.root).unwrap();
    }
}

#[test]
fn test_preseq_numeric_output_write_error_is_not_swallowed() {
    let fixture = estimable_fixture();
    let directory = fixture.output_dir(false);
    fs::create_dir_all(directory.join("sample.lc_extrap.txt")).unwrap();
    let output = fixture.run(true, false);
    assert!(!output.status.success(), "{}", messages(&output));
    assert!(
        messages(&output).contains("Failed to create preseq output"),
        "{}",
        messages(&output)
    );
    assert!(!directory.join("sample.preseq_skipped.txt").exists());
    fs::remove_dir_all(fixture.root).unwrap();
}

#[test]
fn test_preseq_disabled_keeps_existing_outputs() {
    let fixture = Fixture::new(&[1, 1]);
    let directory = fixture.output_dir(false);
    fs::create_dir_all(&directory).unwrap();
    let curve = directory.join("sample.lc_extrap.txt");
    let status = directory.join("sample.preseq_skipped.txt");
    fs::write(&curve, "existing curve\n").unwrap();
    fs::write(&status, "existing status\n").unwrap();
    let output = fixture.run(false, false);
    assert!(output.status.success(), "{}", messages(&output));
    assert_eq!(fs::read_to_string(curve).unwrap(), "existing curve\n");
    assert_eq!(fs::read_to_string(status).unwrap(), "existing status\n");
    fs::remove_dir_all(fixture.root).unwrap();
}

#[test]
fn test_preseq_skipped_status_sanitizes_record_delimiters() {
    let fixture = Fixture::new(&[]);
    let path = fixture.root.join("skipped.txt");
    rustqc::rna::preseq::write_skipped_output(&path, "bad\treason\ncontinued\rline").unwrap();
    assert_eq!(
        fs::read_to_string(path).unwrap(),
        "STATUS\tREASON\nskipped\tbad reason continued line\n"
    );
    fs::remove_dir_all(fixture.root).unwrap();
}
