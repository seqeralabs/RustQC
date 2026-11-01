//! Native regression controls for explicit, file-order strandedness sampling.
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use rust_htslib::bam;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::HeaderView;

fn root(label: &str) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = std::env::temp_dir().join(format!(
        "rustqc-sampling-{label}-{}-{nonce}",
        std::process::id()
    ));
    fs::create_dir(&path).unwrap();
    path
}

fn write_alignment(
    path: &Path,
    chroms: &[(&str, u64)],
    rows: &[(String, u64)],
    format: bam::Format,
    reference: Option<&Path>,
) {
    let mut header = bam::Header::new();
    header.push_record(
        HeaderRecord::new(b"HD")
            .push_tag(b"VN", &"1.6")
            .push_tag(b"SO", &"coordinate"),
    );
    for (chrom, length) in chroms {
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", chrom)
                .push_tag(b"LN", length),
        );
    }
    let view = HeaderView::from_header(&header);
    let mut writer = bam::Writer::from_path(path, &header, format).unwrap();
    if let Some(reference) = reference {
        writer.set_reference(reference).unwrap();
    }
    for (fields, count) in rows {
        let line = format!("{fields}\t{}\t{}", "A".repeat(50), "I".repeat(50));
        let record = bam::Record::from_sam(&view, line.as_bytes()).unwrap();
        for _ in 0..*count {
            writer.write(&record).unwrap();
        }
    }
    drop(writer);
    bam::index::build(path, None, bam::index::Type::Bai, 1).unwrap();
}

fn write_bam(path: &Path, chroms: &[(&str, u64)], rows: &[(String, u64)]) {
    write_alignment(path, chroms, rows, bam::Format::Bam, None);
}

fn run(
    root: &Path,
    threads: usize,
    cli_cap: Option<u64>,
    yaml_cap: Option<u64>,
    additions: &str,
    expected: [f64; 3],
) -> String {
    let folder = root.join(format!("t{threads}-cli{cli_cap:?}-yaml{yaml_cap:?}"));
    fs::create_dir(&folder).unwrap();
    let mut config = include_str!("data/infer_exp_sampling.yaml").to_string();
    if let Some(cap) = yaml_cap {
        config = config.replace(
            "  infer_experiment:\n    enabled: true\n",
            &format!("  infer_experiment:\n    enabled: true\n    sample_size: {cap}\n"),
        );
    }
    config.push_str(additions);
    let config_path = folder.join("config.yaml");
    fs::write(&config_path, config).unwrap();
    let mut command = Command::new(env!("CARGO_BIN_EXE_rustqc"));
    command.env_remove("RUSTQC_INFER_EXPERIMENT_SAMPLE_SIZE");
    command
        .args(["rna", "-q", "--skip-dup-check", "--config"])
        .arg(&config_path)
        .arg("--gtf")
        .arg(root.join("a.gtf"))
        .arg("--outdir")
        .arg(&folder)
        .args(["--threads", &threads.to_string()])
        .arg(root.join("s.bam"));
    if let Some(cap) = cli_cap {
        command.args(["--infer-experiment-sample-size", &cap.to_string()]);
    }
    let output = command.output().unwrap();
    assert!(
        output.status.success(),
        "sampling run failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text =
        fs::read_to_string(folder.join("rseqc/infer_experiment/s.infer_experiment.txt")).unwrap();
    let observed: Vec<f64> = text
        .lines()
        .filter(|line| line.starts_with("Fraction of reads"))
        .map(|line| line.rsplit_once(": ").unwrap().1.parse().unwrap())
        .collect();
    assert_eq!(
        observed, expected,
        "sampling t{threads} CLI={cli_cap:?} YAML={yaml_cap:?}:\n{text}"
    );
    fs::read_to_string(folder.join("samtools/s.flagstat")).unwrap()
}

fn write_gtf(path: &Path, rows: &[(&str, u64, u64, &str, &str)]) {
    let mut annotation = String::new();
    for (chrom, start, end, strand, id) in rows {
        annotation.push_str(&format!("{chrom}\tt\texon\t{start}\t{end}\t.\t{strand}\t.\tgene_id \"{id}\"; transcript_id \"{id}\";\n"));
    }
    fs::write(path, annotation).unwrap();
}

#[test]
fn test_explicit_caps_file_order_default_and_yaml() {
    let root = root("large");
    write_gtf(
        &root.join("a.gtf"),
        &[("chr1", 1, 100, "+", "one"), ("chr2", 1, 100, "+", "two")],
    );
    let rows: Vec<_> = [
        ("sense1", "chr1", 0, 156_000),
        ("anti1", "chr1", 16, 44_000),
        ("sense2", "chr2", 0, 100_000),
    ]
    .into_iter()
    .map(|(id, chrom, flags, count)| {
        (
            format!("{id}\t{flags}\t{chrom}\t11\t60\t50M\t*\t0\t0"),
            count,
        )
    })
    .collect();
    write_bam(
        &root.join("s.bam"),
        &[("chr1", 100_000), ("chr2", 200_000)],
        &rows,
    );
    let mut flagstat = None;
    for threads in [1, 2, 4, 8] {
        for (cli, yaml, expected) in [
            (None, None, [0.0, 0.8533, 0.1467]),
            (Some(u64::MAX), None, [0.0, 0.8533, 0.1467]),
            (Some(1), None, [0.0, 1.0, 0.0]),
            (Some(200_000), None, [0.0, 0.78, 0.22]),
            (Some(400_000), None, [0.0, 0.8533, 0.1467]),
            (None, Some(200_000), [0.0, 0.78, 0.22]),
            (Some(1), Some(200_000), [0.0, 1.0, 0.0]),
        ] {
            let observed = run(&root, threads, cli, yaml, "", expected);
            assert!(
                observed.starts_with("300000 + 0 in total"),
                "sampling must not truncate other analyses"
            );
            if let Some(ref baseline) = flagstat {
                assert_eq!(&observed, baseline);
            } else {
                flagstat = Some(observed);
            }
        }
    }
    let _ = fs::remove_dir_all(root);
}

#[test]
fn test_explicit_caps_filtering_and_ambiguous_overlap() {
    let root = root("eligible");
    write_gtf(
        &root.join("a.gtf"),
        &[
            ("chr1", 1001, 1400, "+", "plus"),
            ("chr1", 1101, 1200, "-", "minus"),
            ("chr2", 1001, 2000, "-", "other"),
        ],
    );
    let rows: Vec<_> = [
        ("no_overlap", 0, "chr1", 11, 60),
        ("qcfail", 512, "chr1", 1001, 60),
        ("duplicate", 1024, "chr1", 1002, 60),
        ("secondary", 256, "chr1", 1003, 60),
        ("low_mapq", 0, "chr1", 1004, 29),
        ("unmapped", 4, "chr1", 1005, 60),
        ("ambiguous", 0, "chr1", 1110, 60),
        ("sense", 0, "chr1", 1250, 60),
        ("supplementary", 2048, "chr1", 1251, 60),
        ("antisense", 0, "chr2", 1001, 60),
    ]
    .into_iter()
    .map(|(id, flags, chrom, pos, mapq)| {
        (
            format!("{id}\t{flags}\t{chrom}\t{pos}\t{mapq}\t50M\t*\t0\t0"),
            1,
        )
    })
    .collect();
    write_bam(
        &root.join("s.bam"),
        &[("chr1", 100_000), ("chr2", 200_000)],
        &rows,
    );
    let baseline = run(&root, 1, None, None, "", [0.25, 0.5, 0.25]);
    for threads in [1, 2, 4, 8] {
        for (cap, expected) in [
            (1, [1.0, 0.0, 0.0]),
            (2, [0.5, 0.5, 0.0]),
            (3, [0.3333, 0.6667, 0.0]),
            (10, [0.25, 0.5, 0.25]),
        ] {
            let flagstat = run(&root, threads, Some(cap), None, "", expected);
            assert_eq!(
                flagstat, baseline,
                "sampling must preserve uncapped flagstat"
            );
        }
    }
    let _ = fs::remove_dir_all(root);
}

#[test]
fn test_explicit_caps_resolve_chromosome_names() {
    for additions in [
        "  chromosome_prefix: chr\n",
        "  chromosome_mapping:\n    chr1: \"1\"\n  chromosome_prefix: ignored\n",
    ] {
        let root = root("chromosome-names");
        write_gtf(&root.join("a.gtf"), &[("chr1", 1, 100, "+", "one")]);
        write_bam(
            &root.join("s.bam"),
            &[("1", 100_000)],
            &[
                ("sense\t0\t1\t11\t60\t50M\t*\t0\t0".into(), 1),
                ("anti\t16\t1\t12\t60\t50M\t*\t0\t0".into(), 1),
            ],
        );
        for threads in [1, 2, 4, 8] {
            assert!(
                run(&root, threads, Some(1), None, additions, [0.0, 1.0, 0.0])
                    .starts_with("2 + 0 in total")
            );
        }
        let _ = fs::remove_dir_all(root);
    }
}

#[test]
fn test_explicit_caps_keep_cram_reference_configuration() {
    let root = root("cram-reference");
    write_gtf(
        &root.join("a.gtf"),
        &[("chr1", 1, 1000, "+", "one"), ("chr2", 1, 1000, "+", "two")],
    );
    let reference = root.join("reference.fa");
    fs::write(
        &reference,
        format!(">chr1\n{}\n>chr2\n{}\n", "A".repeat(1000), "A".repeat(1000)),
    )
    .unwrap();
    rust_htslib::faidx::build(&reference).unwrap();
    let rows: Vec<_> = [
        "first\t0\tchr1\t101\t60\t50M\t*\t0\t0",
        "anti1\t16\tchr1\t201\t60\t50M\t*\t0\t0",
        "anti2\t16\tchr2\t101\t60\t50M\t*\t0\t0",
        "anti3\t16\tchr2\t151\t60\t50M\t*\t0\t0",
        "anti4\t16\tchr2\t201\t60\t50M\t*\t0\t0",
        "anti5\t16\tchr2\t251\t60\t50M\t*\t0\t0",
    ]
    .into_iter()
    .map(|line| (line.to_owned(), 1))
    .collect();
    let cram = root.join("s.cram");
    write_alignment(
        &cram,
        &[("chr1", 1000), ("chr2", 1000)],
        &rows,
        bam::Format::Cram,
        Some(&reference),
    );
    let config = root.join("config.yaml");
    fs::write(&config, include_str!("data/infer_exp_sampling.yaml")).unwrap();
    for threads in [1, 2, 4, 8] {
        let outdir = root.join(format!("t{threads}"));
        let output = Command::new(env!("CARGO_BIN_EXE_rustqc"))
            .env_remove("RUSTQC_INFER_EXPERIMENT_SAMPLE_SIZE")
            .args(["rna", "-q", "--skip-dup-check", "--config"])
            .arg(&config)
            .arg("--gtf")
            .arg(root.join("a.gtf"))
            .arg("--outdir")
            .arg(&outdir)
            .args([
                "--threads",
                &threads.to_string(),
                "--infer-experiment-sample-size",
                "1",
                "--reference",
            ])
            .arg(&reference)
            .arg(&cram)
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "CRAM sampling failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let text = fs::read_to_string(outdir.join("rseqc/infer_experiment/s.infer_experiment.txt"))
            .unwrap();
        assert!(text.contains("This is SingleEnd Data"));
        assert!(
            text.contains("Fraction of reads explained by \"++,--\": 1.0000"),
            "{text}"
        );
        assert!(
            text.contains("Fraction of reads explained by \"+-,-+\": 0.0000"),
            "{text}"
        );
        assert!(fs::read_to_string(outdir.join("samtools/s.flagstat"))
            .unwrap()
            .starts_with("6 + 0 in total"));
    }
    let _ = fs::remove_dir_all(root);
}

#[test]
fn test_explicit_caps_sample_paired_records_in_file_order() {
    let root = root("paired-prefix");
    write_gtf(
        &root.join("a.gtf"),
        &[("chr1", 1, 1000, "+", "one"), ("chr2", 1, 1000, "+", "two")],
    );
    let rows: Vec<_> = [
        "a\t99\tchr1\t101\t60\t50M\t=\t151\t100",
        "a\t147\tchr1\t151\t60\t50M\t=\t101\t-100",
        "b\t163\tchr1\t201\t60\t50M\t=\t251\t100",
        "b\t83\tchr1\t251\t60\t50M\t=\t201\t-100",
        "c\t163\tchr2\t101\t60\t50M\t=\t151\t100",
        "c\t83\tchr2\t151\t60\t50M\t=\t101\t-100",
    ]
    .into_iter()
    .map(|line| (line.to_string(), 1))
    .collect();
    write_bam(
        &root.join("s.bam"),
        &[("chr1", 100_000), ("chr2", 200_000)],
        &rows,
    );
    let config = root.join("config.yaml");
    fs::write(&config, include_str!("data/infer_exp_sampling.yaml")).unwrap();
    for threads in [1, 4] {
        for (cap, expected) in [
            (1, [0.0, 1.0, 0.0]),
            (3, [0.0, 0.6667, 0.3333]),
            (10, [0.0, 0.3333, 0.6667]),
        ] {
            let outdir = root.join(format!("t{threads}-cap{cap}"));
            let output = Command::new(env!("CARGO_BIN_EXE_rustqc"))
                .env_remove("RUSTQC_INFER_EXPERIMENT_SAMPLE_SIZE")
                .args(["rna", "-q", "--skip-dup-check", "--paired", "--config"])
                .arg(&config)
                .arg("--gtf")
                .arg(root.join("a.gtf"))
                .arg("--outdir")
                .arg(&outdir)
                .args([
                    "--threads",
                    &threads.to_string(),
                    "--infer-experiment-sample-size",
                    &cap.to_string(),
                ])
                .arg(root.join("s.bam"))
                .output()
                .unwrap();
            assert!(
                output.status.success(),
                "paired sampling failed: {}",
                String::from_utf8_lossy(&output.stderr)
            );
            let text =
                fs::read_to_string(outdir.join("rseqc/infer_experiment/s.infer_experiment.txt"))
                    .unwrap();
            assert!(text.contains("This is PairEnd Data"));
            let actual: Vec<f64> = text
                .lines()
                .filter(|line| line.starts_with("Fraction of reads"))
                .map(|line| line.rsplit_once(": ").unwrap().1.parse().unwrap())
                .collect();
            assert_eq!(actual, expected, "t{threads} cap{cap}:\n{text}");
            assert!(fs::read_to_string(outdir.join("samtools/s.flagstat"))
                .unwrap()
                .starts_with("6 + 0 in total"));
        }
    }
    let _ = fs::remove_dir_all(root);
}

#[test]
fn test_explicit_zero_cap_retains_unknown_frame() {
    let root = root("zero-frame");
    write_gtf(&root.join("a.gtf"), &[("chr1", 1, 1000, "+", "one")]);
    write_bam(
        &root.join("s.bam"),
        &[("chr1", 100_000)],
        &[("sense\t0\tchr1\t11\t60\t50M\t*\t0\t0".into(), 1)],
    );
    for threads in [1, 4] {
        for cli in [false, true] {
            let folder = root.join(format!("t{threads}-cli{cli}"));
            fs::create_dir(&folder).unwrap();
            let mut config = include_str!("data/infer_exp_sampling.yaml").to_string();
            if !cli {
                config = config.replace(
                    "  infer_experiment:\n    enabled: true\n",
                    "  infer_experiment:\n    enabled: true\n    sample_size: 0\n",
                );
            }
            let config_path = folder.join("config.yaml");
            fs::write(&config_path, config).unwrap();
            let mut command = Command::new(env!("CARGO_BIN_EXE_rustqc"));
            command.env_remove("RUSTQC_INFER_EXPERIMENT_SAMPLE_SIZE");
            command
                .args(["rna", "-q", "--skip-dup-check", "--config"])
                .arg(&config_path)
                .arg("--gtf")
                .arg(root.join("a.gtf"))
                .arg("--outdir")
                .arg(&folder)
                .args(["--threads", &threads.to_string()])
                .arg(root.join("s.bam"));
            if cli {
                command.args(["--infer-experiment-sample-size", "0"]);
            }
            let output = command.output().unwrap();
            assert!(
                output.status.success(),
                "zero-cap run failed: {}",
                String::from_utf8_lossy(&output.stderr)
            );
            assert_eq!(
                fs::read(folder.join("rseqc/infer_experiment/s.infer_experiment.txt")).unwrap(),
                b"Unknown Data type\n".to_vec()
            );
            assert!(fs::read_to_string(folder.join("samtools/s.flagstat"))
                .unwrap()
                .starts_with("1 + 0 in total"));
        }
    }
    let _ = fs::remove_dir_all(root);
}

#[test]
fn test_explicit_caps_reset_for_each_input_in_either_order() {
    let root = root("multiple-inputs");
    write_gtf(&root.join("a.gtf"), &[("chr1", 1, 1000, "+", "one")]);
    let sense = root.join("sense_first.bam");
    let anti = root.join("anti_first.bam");
    for (path, first, rest) in [(&sense, 0, 16), (&anti, 16, 0)] {
        let rows: Vec<_> = [first, rest, rest, rest]
            .into_iter()
            .enumerate()
            .map(|(index, flags)| {
                (
                    format!(
                        "r{index}\t{flags}\tchr1\t{}\t60\t50M\t*\t0\t0",
                        101 + index * 100
                    ),
                    1,
                )
            })
            .collect();
        write_bam(path, &[("chr1", 100_000), ("chr2", 200_000)], &rows);
    }
    let config = root.join("config.yaml");
    fs::write(&config, include_str!("data/infer_exp_sampling.yaml")).unwrap();
    for threads in [1, 4] {
        for reverse in [false, true] {
            let outdir = root.join(format!("t{threads}-reversed{reverse}"));
            let (first, second) = if reverse {
                (&anti, &sense)
            } else {
                (&sense, &anti)
            };
            let output = Command::new(env!("CARGO_BIN_EXE_rustqc"))
                .env_remove("RUSTQC_INFER_EXPERIMENT_SAMPLE_SIZE")
                .args(["rna", "-q", "--skip-dup-check", "--config"])
                .arg(&config)
                .arg("--gtf")
                .arg(root.join("a.gtf"))
                .arg("--outdir")
                .arg(&outdir)
                .args([
                    "--threads",
                    &threads.to_string(),
                    "--infer-experiment-sample-size",
                    "1",
                ])
                .arg(first)
                .arg(second)
                .output()
                .unwrap();
            assert!(
                output.status.success(),
                "multiple-input run failed: {}",
                String::from_utf8_lossy(&output.stderr)
            );
            for (name, expected) in [
                ("sense_first", [0.0, 1.0, 0.0]),
                ("anti_first", [0.0, 0.0, 1.0]),
            ] {
                let text = fs::read_to_string(outdir.join(format!(
                    "rseqc/infer_experiment/{name}.infer_experiment.txt"
                )))
                .unwrap();
                let actual: Vec<f64> = text
                    .lines()
                    .filter(|line| line.starts_with("Fraction of reads"))
                    .map(|line| line.rsplit_once(": ").unwrap().1.parse().unwrap())
                    .collect();
                assert_eq!(
                    actual, expected,
                    "{name}, t{threads}, reversed{reverse}:\n{text}"
                );
                assert!(
                    fs::read_to_string(outdir.join(format!("samtools/{name}.flagstat")))
                        .unwrap()
                        .starts_with("4 + 0 in total")
                );
            }
        }
    }
    let _ = fs::remove_dir_all(root);
}
