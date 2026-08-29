//! RustQC - Fast quality control tools for sequencing data
//!
//! A collection of Rust-based QC tools for bioinformatics.
//! The `rna` subcommand runs all RNA-Seq QC analyses in a single pass:
//! dupRadar duplication rate analysis, featureCounts-compatible gene counting,
//! 8 RSeQC-equivalent tools (bam_stat, infer_experiment, read_duplication,
//! read_distribution, junction_annotation, junction_saturation, inner_distance, TIN),
//! preseq library complexity extrapolation, samtools-compatible outputs
//! (flagstat, idxstats, stats), and Qualimap gene body coverage profiling.
//! Individual tools can be disabled via the YAML config file.

mod citations;
mod cli;
mod ui;

use anyhow::{ensure, Context, Result};
use indexmap::IndexMap;
use log::debug;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use rustqc::io::{format_count, format_duration, format_pct};
use rustqc::{common, config, cpu, gtf, rna, summary};

use ui::{Ui, Verbosity};

use rust_htslib::bam;
use rust_htslib::bam::Read as BamRead;

use rna::rseqc::accumulators::{RseqcAccumulators, RseqcAnnotations, RseqcConfig};

/// Common BAM filename suffixes added by alignment and duplicate-marking tools.
///
/// Return the current UTC time formatted as ISO-8601 (e.g. `2026-03-07T12:34:56Z`).
///
/// Uses only `std::time` — no external chrono dependency.
fn format_utc_now() -> String {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();

    // Convert epoch seconds to a UTC date-time the hard way.
    let days = (secs / 86400) as i64;
    let time_of_day = secs % 86400;
    let hh = time_of_day / 3600;
    let mm = (time_of_day % 3600) / 60;
    let ss = time_of_day % 60;

    // Days since 1970-01-01 → (year, month, day) using the civil-from-days algorithm.
    // Ref: Howard Hinnant, chrono-Compatible Low-Level Date Algorithms
    // <http://howardhinnant.github.io/date_algorithms.html>
    let z = days + 719468;
    let era = (if z >= 0 { z } else { z - 146096 }) / 146097;
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };

    format!("{:04}-{:02}-{:02}T{:02}:{:02}:{:02}Z", y, m, d, hh, mm, ss)
}

fn main() -> Result<()> {
    // Guard against running a SIMD-optimized binary on incompatible hardware.
    // Must run before any auto-vectorized code to prevent SIGILL.
    cpu::check_cpu_compat()?;

    let cli = cli::parse_args();

    // Determine verbosity from CLI flags
    let (quiet, verbose) = match &cli.command {
        cli::Commands::Rna(args) => (args.quiet, args.verbose),
        cli::Commands::Dna(args) => (args.quiet, args.verbose),
        cli::Commands::Reads(args) => (args.quiet, args.verbose),
        cli::Commands::Protein(args) => match &args.mode {
            cli::ProteinMode::Sequence(args) => (args.quiet, args.verbose),
            #[cfg(feature = "proteomics")]
            cli::ProteinMode::Spectra(args) => (args.quiet, args.verbose),
        },
    };
    let verbosity = match (quiet, verbose) {
        (true, _) => Verbosity::Quiet,
        (_, true) => Verbosity::Verbose,
        _ => Verbosity::Normal,
    };

    // Initialize env_logger: only for debug/trace (user-facing output goes through Ui).
    // In verbose mode, lower the threshold so debug!() messages are visible too.
    let log_level = match verbosity {
        Verbosity::Quiet => "warn",
        Verbosity::Normal => "warn",
        Verbosity::Verbose => "debug",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level))
        .format_timestamp(None)
        .init();

    let ui = Ui::new(verbosity);

    match cli.command {
        cli::Commands::Rna(args) => run_rna(args, &ui),
        cli::Commands::Dna(args) => run_dna(args, &ui),
        cli::Commands::Reads(args) => run_reads(args, &ui),
        cli::Commands::Protein(args) => run_protein(args, &ui),
    }
}

/// Run `reads`: raw FASTQ quality control.
///
/// Records stream through the accumulator and are never all held at once, so
/// memory does not grow with the file.
fn run_reads(args: cli::ReadsArgs, ui: &Ui) -> Result<()> {
    use rustqc::reads::{fastq, metrics::ReadMetrics, output};

    let run_start = Instant::now();
    let timestamp_start = format_utc_now();

    let (_merged, config_paths) = config::load_merged_config(args.config.as_deref())?;
    let outdir = Path::new(&args.outdir);
    std::fs::create_dir_all(outdir)
        .with_context(|| format!("Failed to create output directory: {}", outdir.display()))?;

    ui.header(
        env!("CARGO_PKG_VERSION"),
        env!("GIT_SHORT_HASH"),
        env!("BUILD_TIMESTAMP"),
        Some(&rustqc::cpu::cpu_info_line()),
    );
    for (path, source) in &config_paths {
        ui.config("Config", &format!("{} ({source})", path.display()));
    }
    ui.config("Output dir", &args.outdir);

    let dir = if args.flat_output {
        outdir.to_path_buf()
    } else {
        outdir.join("reads")
    };
    std::fs::create_dir_all(&dir)?;

    let mut table_rows = Vec::new();
    let mut inputs = Vec::new();

    for path in &args.input {
        let file_start = Instant::now();
        let name = Path::new(path)
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or(path.as_str())
            .to_string();

        let mut metrics = ReadMetrics::default();
        match fastq::for_each_record(Path::new(path), |r| {
            metrics.observe(&r.sequence, &r.quality)
        }) {
            Ok(()) => {
                let sample_name = args.sample_name.clone().unwrap_or_else(|| {
                    // Strip both extensions of a .fastq.gz, not just the last.
                    let stem = Path::new(path)
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("sample");
                    Path::new(stem)
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or(stem)
                        .to_string()
                });

                let data = dir.join(format!("{sample_name}_fastqc_data.txt"));
                output::write_fastqc_data(&name, &metrics, &data)?;
                ui.output_item("reads", &data.display().to_string());

                ui.detail(&format!(
                    "{name}: {} reads, {} bases, mean quality {:.2}",
                    metrics.reads,
                    metrics.bases,
                    metrics.average_quality()
                ));

                inputs.push(summary::InputSummary {
                    bam_file: path.clone(),
                    status: "success".to_string(),
                    error: None,
                    runtime_seconds: file_start.elapsed().as_secs_f64(),
                    counting: None,
                    dupradar: None,
                    dna: None,
                    outputs: vec![summary::OutputFile {
                        tool: "reads".to_string(),
                        path: data.display().to_string(),
                    }],
                });
                table_rows.push((name, metrics));
            }
            Err(e) => {
                ui.bam_result_err(&name, &format!("{e:#}"));
                inputs.push(summary::InputSummary {
                    bam_file: path.clone(),
                    status: "failed".to_string(),
                    error: Some(format!("{e:#}")),
                    runtime_seconds: file_start.elapsed().as_secs_f64(),
                    counting: None,
                    dupradar: None,
                    dna: None,
                    outputs: Vec::new(),
                });
            }
        }
    }

    if !table_rows.is_empty() {
        let stats_path = dir.join("read_stats.tsv");
        output::write_seqkit_stats(&table_rows, &stats_path)?;
        ui.output_item("reads", &stats_path.display().to_string());
    }

    if let Some(ref json_path) = args.json_summary {
        let summary = summary::RunSummary {
            version: env!("CARGO_PKG_VERSION").to_string(),
            commit: env!("GIT_SHORT_HASH").to_string(),
            binary_target: cpu::binary_target().to_string(),
            cpu_features: cpu::detected_features()
                .iter()
                .map(|s| s.to_string())
                .collect(),
            timestamp_start,
            timestamp_end: format_utc_now(),
            runtime_seconds: run_start.elapsed().as_secs_f64(),
            inputs,
        };
        let json = serde_json::to_string_pretty(&summary)?;
        if json_path == "-" {
            println!("{json}");
        } else {
            let path = if json_path.is_empty() {
                outdir.join("rustqc_summary.json")
            } else {
                PathBuf::from(json_path)
            };
            std::fs::write(&path, json)
                .with_context(|| format!("Failed to write JSON summary: {}", path.display()))?;
        }
    }

    ui.finish("Read QC", run_start.elapsed());
    Ok(())
}

/// Run the protein QC pipeline, dispatching on the chosen mode.
fn run_protein(args: cli::ProteinArgs, ui: &Ui) -> Result<()> {
    match args.mode {
        cli::ProteinMode::Sequence(args) => run_protein_sequence(args, ui),
        #[cfg(feature = "proteomics")]
        cli::ProteinMode::Spectra(args) => run_protein_spectra(args, ui),
    }
}

/// Run `protein spectra`: mass spectrometry run QC from mzML.
#[cfg(feature = "proteomics")]
fn run_protein_spectra(args: cli::ProteinSpectraArgs, ui: &Ui) -> Result<()> {
    use rustqc::protein::spectra;

    let run_start = Instant::now();
    let timestamp_start = format_utc_now();

    let (merged, config_paths) = config::load_merged_config(args.config.as_deref())?;
    let config = merged.protein;
    let flat_output = args.flat_output || config.flat_output;

    let outdir = Path::new(&args.outdir);
    std::fs::create_dir_all(outdir)
        .with_context(|| format!("Failed to create output directory: {}", outdir.display()))?;

    ui.header(
        env!("CARGO_PKG_VERSION"),
        env!("GIT_SHORT_HASH"),
        env!("BUILD_TIMESTAMP"),
        Some(&rustqc::cpu::cpu_info_line()),
    );
    for (path, source) in &config_paths {
        ui.config("Config", &format!("{} ({source})", path.display()));
    }
    ui.config("Output dir", &args.outdir);

    let dir = if flat_output {
        outdir.to_path_buf()
    } else {
        outdir.join("spectra")
    };
    std::fs::create_dir_all(&dir)?;

    let mut inputs = Vec::new();
    for path in &args.input {
        let file_start = Instant::now();
        let name = Path::new(path)
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or(path.as_str())
            .to_string();

        match spectra::analyse(Path::new(path)) {
            Ok(metrics) => {
                let sample_name = args.sample_name.clone().unwrap_or_else(|| {
                    Path::new(path)
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("sample")
                        .to_string()
                });
                let report = dir.join(format!("{sample_name}.spectra_report.txt"));
                spectra::output::write_report(&name, &metrics, &report)?;
                ui.output_item("protein spectra", &report.display().to_string());
                ui.detail(&format!(
                    "{name}: {} spectra, {} peaks",
                    metrics.total_spectra(),
                    metrics.total_peaks()
                ));
                if metrics.precursors > 0 && metrics.precursors_without_charge == metrics.precursors
                {
                    ui.warn(&format!(
                        "{name}: no precursor charge states are annotated, so charge metrics are unavailable"
                    ));
                }

                inputs.push(summary::InputSummary {
                    bam_file: path.clone(),
                    status: "success".to_string(),
                    error: None,
                    runtime_seconds: file_start.elapsed().as_secs_f64(),
                    counting: None,
                    dupradar: None,
                    dna: None,
                    outputs: vec![summary::OutputFile {
                        tool: "protein spectra".to_string(),
                        path: report.display().to_string(),
                    }],
                });
            }
            Err(e) => {
                ui.bam_result_err(&name, &format!("{e:#}"));
                inputs.push(summary::InputSummary {
                    bam_file: path.clone(),
                    status: "failed".to_string(),
                    error: Some(format!("{e:#}")),
                    runtime_seconds: file_start.elapsed().as_secs_f64(),
                    counting: None,
                    dupradar: None,
                    dna: None,
                    outputs: Vec::new(),
                });
            }
        }
    }

    if let Some(ref json_path) = args.json_summary {
        let summary = summary::RunSummary {
            version: env!("CARGO_PKG_VERSION").to_string(),
            commit: env!("GIT_SHORT_HASH").to_string(),
            binary_target: cpu::binary_target().to_string(),
            cpu_features: cpu::detected_features()
                .iter()
                .map(|s| s.to_string())
                .collect(),
            timestamp_start,
            timestamp_end: format_utc_now(),
            runtime_seconds: run_start.elapsed().as_secs_f64(),
            inputs,
        };
        let json = serde_json::to_string_pretty(&summary)?;
        if json_path == "-" {
            println!("{json}");
        } else {
            let path = if json_path.is_empty() {
                outdir.join("rustqc_summary.json")
            } else {
                PathBuf::from(json_path)
            };
            std::fs::write(&path, json)
                .with_context(|| format!("Failed to write JSON summary: {}", path.display()))?;
        }
    }

    ui.finish("Protein spectra QC", run_start.elapsed());
    Ok(())
}

/// Run `protein sequence`: FASTA statistics, composition and defects.
///
/// Each input file is summarised on its own, and the seqkit-compatible table
/// carries one row per file, which is how seqkit reports several files too.
fn run_protein_sequence(args: cli::ProteinSequenceArgs, ui: &Ui) -> Result<()> {
    use rustqc::protein::sequence::output;

    let run_start = Instant::now();
    let timestamp_start = format_utc_now();

    let (merged, config_paths) = config::load_merged_config(args.config.as_deref())?;
    let config = merged.protein;
    let flat_output = args.flat_output || config.flat_output;
    let expect_stop = args.expect_stop || config.sequence.expect_stop;
    let min_length = if args.min_length > 0 {
        args.min_length
    } else {
        config.sequence.min_length
    };

    let outdir = Path::new(&args.outdir);
    std::fs::create_dir_all(outdir)
        .with_context(|| format!("Failed to create output directory: {}", outdir.display()))?;

    ui.header(
        env!("CARGO_PKG_VERSION"),
        env!("GIT_SHORT_HASH"),
        env!("BUILD_TIMESTAMP"),
        Some(&rustqc::cpu::cpu_info_line()),
    );
    for (path, source) in &config_paths {
        ui.config("Config", &format!("{} ({source})", path.display()));
    }
    ui.config("Output dir", &args.outdir);
    if min_length > 0 {
        ui.config("Min length", &min_length.to_string());
    }

    let dir = if flat_output {
        outdir.to_path_buf()
    } else {
        outdir.join("sequence")
    };
    std::fs::create_dir_all(&dir)?;

    let mut table_rows = Vec::new();
    let mut inputs = Vec::new();

    for path in &args.input {
        let file_start = Instant::now();
        let name = Path::new(path)
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or(path.as_str())
            .to_string();

        match summarise_fasta(path, min_length, expect_stop) {
            Ok((stats, found, kept, skipped)) => {
                if skipped > 0 {
                    ui.detail(&format!(
                        "{name}: skipped {skipped} sequences shorter than {min_length}"
                    ));
                }
                let sample_name = args.sample_name.clone().unwrap_or_else(|| {
                    Path::new(path)
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("sample")
                        .to_string()
                });

                let report = dir.join(format!("{sample_name}.sequence_report.txt"));
                output::write_report(&name, &stats, &found, &report)?;
                ui.output_item("protein sequence", &report.display().to_string());

                if !found.is_clean() {
                    ui.warn(&format!(
                        "{name}: {} defective sequences, see the report",
                        found.count()
                    ));
                }

                inputs.push(summary::InputSummary {
                    bam_file: path.clone(),
                    status: "success".to_string(),
                    error: None,
                    runtime_seconds: file_start.elapsed().as_secs_f64(),
                    counting: None,
                    dupradar: None,
                    dna: None,
                    outputs: vec![summary::OutputFile {
                        tool: "protein sequence".to_string(),
                        path: report.display().to_string(),
                    }],
                });
                ui.detail(&format!("{name}: {kept} sequences"));
                table_rows.push((name, stats));
            }
            Err(e) => {
                ui.bam_result_err(&name, &format!("{e:#}"));
                inputs.push(summary::InputSummary {
                    bam_file: path.clone(),
                    status: "failed".to_string(),
                    error: Some(format!("{e:#}")),
                    runtime_seconds: file_start.elapsed().as_secs_f64(),
                    counting: None,
                    dupradar: None,
                    dna: None,
                    outputs: Vec::new(),
                });
            }
        }
    }

    if !table_rows.is_empty() {
        let stats_path = dir.join("sequence_stats.tsv");
        output::write_seqkit_stats(&table_rows, &stats_path)?;
        ui.output_item("protein sequence", &stats_path.display().to_string());
    }

    if let Some(ref json_path) = args.json_summary {
        let summary = summary::RunSummary {
            version: env!("CARGO_PKG_VERSION").to_string(),
            commit: env!("GIT_SHORT_HASH").to_string(),
            binary_target: cpu::binary_target().to_string(),
            cpu_features: cpu::detected_features()
                .iter()
                .map(|s| s.to_string())
                .collect(),
            timestamp_start,
            timestamp_end: format_utc_now(),
            runtime_seconds: run_start.elapsed().as_secs_f64(),
            inputs,
        };
        let json = serde_json::to_string_pretty(&summary)?;
        if json_path == "-" {
            println!("{json}");
        } else {
            let path = if json_path.is_empty() {
                outdir.join("rustqc_summary.json")
            } else {
                PathBuf::from(json_path)
            };
            std::fs::write(&path, json)
                .with_context(|| format!("Failed to write JSON summary: {}", path.display()))?;
        }
    }

    ui.finish("Protein sequence QC", run_start.elapsed());
    Ok(())
}

/// Read one FASTA and summarise it, returning the statistics, the defects, and
/// how many sequences were kept and skipped by the length filter.
fn summarise_fasta(
    path: &str,
    min_length: usize,
    expect_stop: bool,
) -> Result<(
    rustqc::protein::sequence::stats::SequenceStats,
    rustqc::protein::sequence::defects::Defects,
    usize,
    usize,
)> {
    use rustqc::protein::sequence::{self, defects, stats::SequenceStats};

    let all = sequence::read_fasta(Path::new(path))?;
    let total = all.len();
    let kept: Vec<_> = all.into_iter().filter(|r| r.len() >= min_length).collect();
    let skipped = total - kept.len();
    let stats = SequenceStats::from_records(&kept);
    let found = defects::inspect(&kept, expect_stop);
    Ok((stats, found, kept.len(), skipped))
}

/// Run the DNA QC pipeline: depth of coverage, samtools-compatible outputs
/// and library complexity estimation in a single pass over each input.
///
/// Contigs are processed in parallel, one worker per contig, each holding its
/// own depth array. Input files are processed one after another so that the
/// per-contig parallelism gets the whole thread budget.
fn run_dna(args: cli::DnaArgs, ui: &Ui) -> Result<()> {
    let run_start = Instant::now();
    let timestamp_start = format_utc_now();

    let (merged, config_paths) = config::load_merged_config(args.config.as_deref())?;
    let mut config = merged.dna;

    // CLI flags override the configuration file.
    if !args.depth_thresholds.is_empty() {
        config.mosdepth.thresholds = args.depth_thresholds.clone();
    }
    if let Some(window) = args.window_size {
        config.mosdepth.window_size = Some(window);
    }
    if args.skip_per_base {
        config.mosdepth.skip_per_base = true;
    }
    if args.skip_preseq {
        config.preseq.enabled = false;
    }
    if let Some(seed) = args.preseq_seed {
        config.preseq.seed = seed;
    }
    if let Some(val) = args.preseq_max_extrap {
        config.preseq.max_extrap = val;
    }
    if let Some(val) = args.preseq_step_size {
        config.preseq.step_size = val;
    }
    if let Some(val) = args.preseq_n_bootstraps {
        config.preseq.n_bootstraps = val;
    }
    if let Some(val) = args.preseq_seg_len {
        config.preseq.max_segment_length = val;
    }

    let flat_output = args.flat_output || config.flat_output;
    let outdir = Path::new(&args.outdir);
    std::fs::create_dir_all(outdir)
        .with_context(|| format!("Failed to create output directory: {}", outdir.display()))?;

    ui.header(
        env!("CARGO_PKG_VERSION"),
        env!("GIT_SHORT_HASH"),
        env!("BUILD_TIMESTAMP"),
        Some(&rustqc::cpu::cpu_info_line()),
    );
    for (path, source) in &config_paths {
        ui.config("Config", &format!("{} ({source})", path.display()));
    }
    ui.config("Output dir", &args.outdir);
    ui.config("Threads", &args.threads.to_string());
    if let Some(ref targets) = args.targets {
        ui.config("Targets", targets);
    }
    if let Some(ref baits) = args.baits {
        ui.config("Baits", baits);
    }

    let mut inputs = Vec::new();
    for bam_path in &args.input {
        let bam_start = Instant::now();
        let name = Path::new(bam_path)
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or(bam_path.as_str())
            .to_string();

        match process_single_dna_bam(bam_path, &args, &config, outdir, flat_output, ui) {
            Ok(mut summary) => {
                summary.runtime_seconds = bam_start.elapsed().as_secs_f64();
                ui.bam_result_ok(&name, bam_start.elapsed());
                inputs.push(summary);
            }
            Err(e) => {
                ui.bam_result_err(&name, &format!("{e:#}"));
                inputs.push(summary::InputSummary {
                    bam_file: bam_path.clone(),
                    status: "failed".to_string(),
                    error: Some(format!("{e:#}")),
                    runtime_seconds: bam_start.elapsed().as_secs_f64(),
                    counting: None,
                    dupradar: None,
                    dna: None,
                    outputs: Vec::new(),
                });
            }
        }
    }

    if let Some(ref json_path) = args.json_summary {
        let summary = summary::RunSummary {
            version: env!("CARGO_PKG_VERSION").to_string(),
            commit: env!("GIT_SHORT_HASH").to_string(),
            binary_target: cpu::binary_target().to_string(),
            cpu_features: cpu::detected_features()
                .iter()
                .map(|s| s.to_string())
                .collect(),
            timestamp_start,
            timestamp_end: format_utc_now(),
            runtime_seconds: run_start.elapsed().as_secs_f64(),
            inputs,
        };
        let json = serde_json::to_string_pretty(&summary)?;
        if json_path == "-" {
            println!("{json}");
        } else {
            let path = if json_path.is_empty() {
                outdir.join("rustqc_summary.json")
            } else {
                PathBuf::from(json_path)
            };
            std::fs::write(&path, json)
                .with_context(|| format!("Failed to write JSON summary: {}", path.display()))?;
        }
    }

    let citations_path = outdir.join("CITATIONS.md");
    citations::write_dna_citations(
        &citations_path,
        &config,
        env!("CARGO_PKG_VERSION"),
        env!("GIT_SHORT_HASH"),
    )?;
    ui.output_item("citations", &citations_path.display().to_string());

    ui.finish("DNA QC", run_start.elapsed());
    Ok(())
}

/// Process one alignment file through the DNA pipeline.
fn process_single_dna_bam(
    bam_path: &str,
    args: &cli::DnaArgs,
    config: &config::DnaConfig,
    outdir: &Path,
    flat_output: bool,
    ui: &Ui,
) -> Result<summary::InputSummary> {
    use rustqc::common::bam_stat_accum::BamStatAccum;
    use rustqc::common::preseq::PreseqAccum;
    use rustqc::dna::depth::{DepthAccum, MOSDEPTH_DEFAULT_EXCLUDE};
    use rustqc::dna::gc_bias::{self, GcBiasAccum};
    use rustqc::dna::hs_metrics::{self, HsAccum, HsCounters, HsMetricsResult};
    use rustqc::dna::insert_size::{self, InsertSizeAccum};
    use rustqc::dna::intervals::IntervalSet;
    use rustqc::dna::mosdepth::{output as mos_out, ContigDepth, MosdepthResult};
    use rustqc::dna::qualimap::{ContigQualimap, QualimapAccum};
    use rustqc::dna::qualimap_output;
    use rustqc::dna::wgs_metrics::{self, WgsAccum, WgsCounters, WgsMetricsResult};

    let sample_name = args
        .sample_name
        .clone()
        .or_else(|| config.sample_name.clone())
        .unwrap_or_else(|| {
            Path::new(bam_path)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("sample")
                .to_string()
        });

    let is_cram = bam_path.ends_with(".cram");
    ensure!(
        !is_cram || args.reference.is_some(),
        "CRAM input requires --reference"
    );

    // Read the header once to learn the contigs.
    let header = {
        let reader = bam::IndexedReader::from_path(bam_path)
            .with_context(|| format!("Failed to open alignment file: {bam_path}"))?;
        reader.header().to_owned()
    };
    let mut contigs: Vec<(u32, String, u64)> = (0..header.target_count())
        .map(|tid| {
            let name = String::from_utf8_lossy(header.tid2name(tid)).to_string();
            let len = header.target_len(tid).unwrap_or(0);
            (tid, name, len)
        })
        .collect();
    // Longest first, so the biggest depth arrays are allocated while the pool
    // is emptiest.
    contigs.sort_by_key(|contig| std::cmp::Reverse(contig.2));

    let largest = contigs.first().map(|c| c.2).unwrap_or(0);
    let workers = depth_worker_budget(args.threads, args.max_depth_workers, largest);
    ui.config("Depth workers", &workers.to_string());

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(workers)
        .build()
        .context("Failed to build rayon thread pool")?;

    let thresholds = config.mosdepth.thresholds.clone();
    let window_size = config.mosdepth.window_size;
    let preseq_enabled = config.preseq.enabled;
    let seg_len = config.preseq.max_segment_length;
    let mapq_cut = args.mapq_cut;
    // CollectWgsMetrics needs the reference to count non-N bases, so without
    // one it is skipped rather than reported against a wrong territory.
    let wgs_enabled = config.wgs_metrics.enabled && args.reference.is_some();
    if config.wgs_metrics.enabled && args.reference.is_none() {
        ui.warn("CollectWgsMetrics needs --reference to size the genome territory, skipping");
    }
    let insert_size_enabled = config.insert_size.enabled;
    // GC bias bins reference windows, so it needs the reference just as the
    // WGS metrics do.
    let gc_bias_enabled = config.gc_bias.enabled && args.reference.is_some();
    if config.gc_bias.enabled && args.reference.is_none() {
        ui.warn("CollectGcBiasMetrics needs --reference to bin the genome, skipping");
    }
    let gc_window = config.gc_bias.window_size;

    // Targeted mode is switched on by --targets alone; --baits defaults to it.
    let targets = match args.targets.as_deref() {
        Some(path) => Some(IntervalSet::from_bed(Path::new(path))?),
        None => None,
    };
    let baits = match args.baits.as_deref() {
        Some(path) => Some(IntervalSet::from_bed(Path::new(path))?),
        None => targets.clone(),
    };
    let hs_enabled = config.hs_metrics.enabled && targets.is_some();
    let hs_min_mapq = config.hs_metrics.min_mapping_quality;
    let hs_min_baseq = config.hs_metrics.min_base_quality;
    let qualimap_enabled = config.qualimap.enabled;
    let qualimap_windows = config.qualimap.num_windows;
    let wgs_min_mapq = config.wgs_metrics.min_mapping_quality;
    let wgs_min_baseq = config.wgs_metrics.min_base_quality;
    let coverage_cap = config.wgs_metrics.coverage_cap;

    /// What one contig worker hands back: its depth summary, the read-level
    /// counters, and the optional per-tool accumulators.
    type ContigOutput = (
        ContigDepth,
        BamStatAccum,
        Option<PreseqAccum>,
        Option<(WgsCounters, Vec<u32>)>,
        Option<InsertSizeAccum>,
        Option<GcBiasAccum>,
        Option<(HsCounters, Vec<u32>, Vec<bool>, String)>,
        Option<ContigQualimap>,
    );

    let results: Vec<Result<ContigOutput>> = pool.install(|| {
        contigs
            .par_iter()
            .map(|(tid, name, len)| -> Result<ContigOutput> {
                let mut reader = bam::IndexedReader::from_path(bam_path)
                    .with_context(|| format!("Failed to open alignment file: {bam_path}"))?;
                if let Some(reference) = args.reference.as_deref() {
                    reader
                        .set_reference(reference)
                        .with_context(|| format!("Failed to set reference: {reference}"))?;
                }
                reader
                    .fetch(*tid)
                    .with_context(|| format!("Failed to fetch contig {name}"))?;

                let mut depth = DepthAccum::new(*len, mapq_cut, MOSDEPTH_DEFAULT_EXCLUDE);
                let mut bam_stat = BamStatAccum::default();
                let mut preseq = preseq_enabled.then(|| PreseqAccum::new(seg_len));
                // Picard filters differently from mosdepth, so its coverage
                // needs its own accumulator rather than a correction applied
                // to a shared one.
                let mut wgs = wgs_enabled.then(|| WgsAccum::new(*len, wgs_min_mapq, wgs_min_baseq));
                let mut insert_sizes = insert_size_enabled.then(InsertSizeAccum::new);
                let mut qualimap =
                    qualimap_enabled.then(|| QualimapAccum::new(name, *len, qualimap_windows));

                // GC bias and the targeted metrics both need per-contig
                // context, fetched once here rather than per record.
                let reference_bases: Option<Vec<u8>> = if gc_bias_enabled {
                    let reader = rust_htslib::faidx::Reader::from_path(
                        args.reference.as_deref().unwrap_or_default(),
                    )
                    .with_context(|| "Failed to open the reference FASTA index")?;
                    let length = reader.fetch_seq_len(name) as usize;
                    Some(
                        reader
                            .fetch_seq(name, 0, length.saturating_sub(1))
                            .map(|s| s.to_vec())
                            .with_context(|| format!("Failed to read reference for {name}"))?,
                    )
                } else {
                    None
                };
                let mut gc = reference_bases
                    .as_ref()
                    .map(|bases| GcBiasAccum::new(bases, gc_window));
                let mut hs = hs_enabled.then(|| {
                    HsAccum::new(
                        name,
                        *len,
                        baits.as_ref().unwrap_or_else(|| targets.as_ref().unwrap()),
                        targets.as_ref().unwrap(),
                        hs_min_mapq,
                        hs_min_baseq,
                    )
                });

                let mut record = bam::Record::new();
                while let Some(result) = reader.read(&mut record) {
                    result.context("Failed to read record")?;
                    depth.process_read(&record);
                    bam_stat.process_read(&record, mapq_cut);
                    if let Some(accum) = preseq.as_mut() {
                        accum.process_read(&record);
                    }
                    if let Some(accum) = wgs.as_mut() {
                        accum.process_read(&record);
                    }
                    if let Some(accum) = insert_sizes.as_mut() {
                        accum.process_read(&record);
                    }
                    if let (Some(accum), Some(bases)) = (gc.as_mut(), reference_bases.as_ref()) {
                        accum.process_read(&record, bases);
                    }
                    if let Some(accum) = hs.as_mut() {
                        accum.process_read(&record);
                    }
                    if let Some(accum) = qualimap.as_mut() {
                        accum.process_read(&record);
                    }
                }

                let depths = depth.into_depths();
                let contig = ContigDepth::from_depths(name, &depths, window_size, &thresholds);
                Ok((
                    contig,
                    bam_stat,
                    preseq,
                    wgs.map(|accum| accum.into_parts()),
                    insert_sizes,
                    gc,
                    hs.map(|accum| {
                        let (counters, depths, mask) = accum.into_parts();
                        (counters, depths, mask, name.clone())
                    }),
                    qualimap.map(|accum| accum.into_result()),
                ))
            })
            .collect()
    });

    let mut per_contig = Vec::new();
    let mut bam_stat_total = BamStatAccum::default();
    let mut preseq_total: Option<PreseqAccum> = None;
    let mut wgs_counters = WgsCounters::default();
    let mut wgs_depths: Vec<u32> = Vec::new();
    let mut saw_wgs = false;
    let mut insert_size_total: Option<InsertSizeAccum> = None;
    let mut gc_total: Option<GcBiasAccum> = None;
    let mut hs_counters = HsCounters::default();
    let mut hs_target_depths: Vec<u32> = Vec::new();
    let mut hs_target_count = 0u64;
    let mut hs_zero_targets = 0u64;
    let mut qualimap_contigs: Vec<ContigQualimap> = Vec::new();
    for result in results {
        let (contig, bam_stat, preseq, wgs, insert_sizes, gc, hs, qualimap) = result?;
        per_contig.push(contig);
        bam_stat_total.merge(bam_stat);
        match (preseq_total.as_mut(), preseq) {
            (Some(total), Some(part)) => total.merge(part),
            (None, part) => preseq_total = part,
            _ => {}
        }
        if let Some((counters, depths)) = wgs {
            saw_wgs = true;
            wgs_counters.merge(&counters);
            // Depths concatenate rather than merge: each worker owns a
            // distinct contig and the metrics span all of them.
            wgs_depths.extend(depths);
        }
        match (insert_size_total.as_mut(), insert_sizes) {
            (Some(total), Some(part)) => total.merge(part),
            (None, part) => insert_size_total = part,
            _ => {}
        }
        match (gc_total.as_mut(), gc) {
            (Some(total), Some(part)) => total.merge(&part),
            (None, part) => gc_total = part,
            _ => {}
        }
        if let Some(part) = qualimap {
            qualimap_contigs.push(part);
        }
        if let Some((counters, depths, mask, contig_name)) = hs {
            hs_counters.merge(&counters);
            for (depth, on_target) in depths.iter().zip(mask.iter()) {
                if *on_target {
                    hs_target_depths.push(*depth);
                }
            }
            if let Some(set) = targets.as_ref() {
                for interval in set.on(&contig_name) {
                    hs_target_count += 1;
                    if (interval.start..interval.end)
                        .all(|p| depths.get(p as usize).copied().unwrap_or(0) == 0)
                    {
                        hs_zero_targets += 1;
                    }
                }
            }
        }
    }

    // Unmapped records carry no contig, so they need their own pass; flagstat
    // and idxstats both report them.
    {
        let mut reader = bam::IndexedReader::from_path(bam_path)
            .with_context(|| format!("Failed to open alignment file: {bam_path}"))?;
        if let Some(reference) = args.reference.as_deref() {
            reader.set_reference(reference).ok();
        }
        if reader.fetch(bam::FetchDefinition::Unmapped).is_ok() {
            // Unmapped records reach no contig worker, yet they still count
            // towards several metrics: flagstat and idxstats report them, HS
            // metrics count them in TOTAL_READS and PF_BASES, and GC bias
            // counts them as clusters. Both accumulators short-circuit on an
            // unmapped record, so an empty contig is enough context here.
            let mut hs_unmapped = hs_enabled.then(|| {
                HsAccum::new(
                    "",
                    0,
                    baits.as_ref().unwrap_or_else(|| targets.as_ref().unwrap()),
                    targets.as_ref().unwrap(),
                    hs_min_mapq,
                    hs_min_baseq,
                )
            });
            let mut gc_unmapped = gc_bias_enabled.then(|| GcBiasAccum::new(&[], gc_window));
            let mut qualimap_unmapped =
                qualimap_enabled.then(|| QualimapAccum::new("", 0, qualimap_windows));

            let mut record = bam::Record::new();
            while let Some(result) = reader.read(&mut record) {
                result.context("Failed to read unmapped record")?;
                bam_stat_total.process_read(&record, mapq_cut);
                if let Some(accum) = hs_unmapped.as_mut() {
                    accum.process_read(&record);
                }
                if let Some(accum) = gc_unmapped.as_mut() {
                    accum.process_read(&record, &[]);
                }
                if let Some(accum) = qualimap_unmapped.as_mut() {
                    accum.process_read(&record);
                }
            }

            if let Some(accum) = hs_unmapped {
                let (counters, _, _) = accum.into_parts();
                hs_counters.merge(&counters);
            }
            if let (Some(total), Some(part)) = (gc_total.as_mut(), gc_unmapped) {
                total.merge(&part);
            }
            if let (Some(first), Some(part)) = (qualimap_contigs.first_mut(), qualimap_unmapped) {
                // The unmapped pass only moves read counters, so folding it
                // into the first contig keeps the totals right without
                // inventing a contig for reads that have none.
                first.counters.merge(&part.into_result().counters);
            }
        }
    }

    // Workers ran longest-contig-first; outputs go out in header order.
    let order: Vec<String> = (0..header.target_count())
        .map(|tid| String::from_utf8_lossy(header.tid2name(tid)).to_string())
        .collect();
    per_contig.sort_by_key(|contig| {
        order
            .iter()
            .position(|name| name == &contig.name)
            .unwrap_or(usize::MAX)
    });

    let result = MosdepthResult {
        contigs: per_contig,
        window_size,
        thresholds: thresholds.clone(),
    };

    let bam_stat_result = bam_stat_total.into_result();
    ensure!(
        args.skip_dup_check || bam_stat_result.duplicates > 0,
        "No duplicate-flagged reads found in {bam_path}. RustQC expects \
         duplicate-marked (not removed) input. Pass --skip-dup-check to override."
    );

    let dir = |name: &str| -> PathBuf {
        if flat_output {
            outdir.to_path_buf()
        } else {
            outdir.join(name)
        }
    };
    let mut written: Vec<summary::OutputFile> = Vec::new();
    let mut record_output = |tool: &str, path: PathBuf| {
        ui.output_item(tool, &path.display().to_string());
        written.push(summary::OutputFile {
            tool: tool.to_string(),
            path: path.display().to_string(),
        });
    };

    if config.mosdepth.enabled {
        let mos_dir = dir("mosdepth");
        std::fs::create_dir_all(&mos_dir)?;
        // Built with format! rather than with_extension: a sample name that
        // contains a dot (test.dna, say) would otherwise lose its last segment.
        let prefix = |suffix: &str| mos_dir.join(format!("{sample_name}.{suffix}"));

        let path = prefix("mosdepth.summary.txt");
        mos_out::write_summary(&result, &path)?;
        record_output("mosdepth", path);

        let path = prefix("mosdepth.global.dist.txt");
        mos_out::write_global_dist(&result, &path)?;
        record_output("mosdepth", path);

        if !config.mosdepth.skip_per_base {
            let path = prefix("per-base.bed.gz");
            mos_out::write_per_base(&result, &path)?;
            record_output("mosdepth", path);
        }

        if window_size.is_some() {
            let path = prefix("mosdepth.region.dist.txt");
            mos_out::write_region_dist(&result, &path)?;
            record_output("mosdepth", path);

            let path = prefix("regions.bed.gz");
            mos_out::write_regions(&result, &path)?;
            record_output("mosdepth", path);

            if !thresholds.is_empty() {
                let path = prefix("thresholds.bed.gz");
                mos_out::write_thresholds(&result, &path)?;
                record_output("mosdepth", path);
            }
        }
    }

    if config.samtools.enabled {
        let sam_dir = dir("samtools");
        std::fs::create_dir_all(&sam_dir)?;

        let path = sam_dir.join(format!("{sample_name}.stats.txt"));
        common::samtools::stats::write_stats(&bam_stat_result, &path)?;
        record_output("samtools stats", path);

        let path = sam_dir.join(format!("{sample_name}.flagstat.txt"));
        common::samtools::flagstat::write_flagstat(&bam_stat_result, &path)?;
        record_output("samtools flagstat", path);

        let refs: Vec<(String, u64)> = (0..header.target_count())
            .map(|tid| {
                (
                    String::from_utf8_lossy(header.tid2name(tid)).to_string(),
                    header.target_len(tid).unwrap_or(0),
                )
            })
            .collect();
        let path = sam_dir.join(format!("{sample_name}.idxstats.txt"));
        common::samtools::idxstats::write_idxstats(&bam_stat_result, &refs, &path)?;
        record_output("samtools idxstats", path);
    }

    if saw_wgs {
        let territory = match args.reference.as_deref() {
            Some(reference) => genome_territory(reference)?,
            // Unreachable: saw_wgs implies a reference was given.
            None => wgs_depths.len() as u64,
        };
        let result = WgsMetricsResult::new(&wgs_depths, wgs_counters, territory, coverage_cap);
        let dir_path = dir("picard").join("wgs_metrics");
        std::fs::create_dir_all(&dir_path)?;
        let path = dir_path.join(format!("{sample_name}.wgs_metrics.txt"));
        wgs_metrics::write_wgs_metrics(&result, &path)?;
        record_output("picard CollectWgsMetrics", path);
    }

    if let Some(accum) = insert_size_total {
        let result = accum.into_result(config.insert_size.deviations);
        if result.rows.is_empty() {
            ui.warn("no paired records with a usable insert size, skipping insert size metrics");
        } else {
            let dir_path = dir("picard").join("insert_size");
            std::fs::create_dir_all(&dir_path)?;
            let path = dir_path.join(format!("{sample_name}.insert_size_metrics.txt"));
            insert_size::write_insert_size_metrics(&result, &path)?;
            record_output("picard CollectInsertSizeMetrics", path);
        }
    }

    if let Some(accum) = gc_total {
        let result = accum.into_result(gc_window);
        let dir_path = dir("picard").join("gc_bias");
        std::fs::create_dir_all(&dir_path)?;
        let detail = dir_path.join(format!("{sample_name}.gc_bias.detail_metrics.txt"));
        gc_bias::write_detail_metrics(&result, &detail)?;
        record_output("picard CollectGcBiasMetrics", detail);
        let summary = dir_path.join(format!("{sample_name}.gc_bias.summary_metrics.txt"));
        gc_bias::write_summary_metrics(&result, &summary)?;
        record_output("picard CollectGcBiasMetrics", summary);
    }

    if hs_enabled {
        let target_set = targets.as_ref().expect("hs_enabled implies --targets");
        let bait_set = baits.as_ref().unwrap_or(target_set);
        let library_size = hs_metrics::estimate_library_size(
            hs_counters.selected_pairs,
            hs_counters.selected_unique_pairs,
        );
        let result = HsMetricsResult {
            bait_set: bait_set.name().to_string(),
            bait_territory: bait_set.territory(),
            target_territory: target_set.territory(),
            genome_size: contigs.iter().map(|(_, _, len)| len).sum(),
            counters: hs_counters,
            target_depths: hs_target_depths,
            zero_coverage_targets: hs_zero_targets,
            target_count: hs_target_count,
            library_size,
        };
        let dir_path = dir("picard").join("hs_metrics");
        std::fs::create_dir_all(&dir_path)?;
        let path = dir_path.join(format!("{sample_name}.hs_metrics.txt"));
        hs_metrics::write_hs_metrics(&result, &path)?;
        record_output("picard CollectHsMetrics", path);
    }

    if !qualimap_contigs.is_empty() {
        // Restore header order, since the workers ran longest contig first.
        qualimap_contigs.sort_by_key(|c| {
            order
                .iter()
                .position(|name| *name == c.name)
                .unwrap_or(usize::MAX)
        });
        let dir_path = dir("qualimap");
        std::fs::create_dir_all(&dir_path)?;
        let results = dir_path.join("genome_results.txt");
        qualimap_output::write_genome_results(&qualimap_contigs, bam_path, &results)?;
        record_output("qualimap", results);
        qualimap_output::write_raw_data(
            &qualimap_contigs,
            &dir_path.join("raw_data_qualimapReport"),
        )?;
        record_output("qualimap", dir_path.join("raw_data_qualimapReport"));
        let report = dir_path.join("qualimapReport.html");
        qualimap_output::write_html_report(&qualimap_contigs, &sample_name, &report)?;
        record_output("qualimap", report);
    }

    if let Some(mut accum) = preseq_total {
        let preseq_dir = dir("preseq");
        std::fs::create_dir_all(&preseq_dir)?;
        accum.finalize();
        let total_reads = accum.total_fragments;
        let n_distinct = accum.n_distinct();
        let histogram = accum.into_histogram();
        match common::preseq::estimate_complexity(
            &histogram,
            total_reads,
            n_distinct,
            &config.preseq,
        ) {
            Ok(preseq_result) => {
                let path = preseq_dir.join(format!("{sample_name}.lc_extrap.txt"));
                common::preseq::write_output(
                    &preseq_result,
                    &path,
                    config.preseq.confidence_level,
                )?;
                record_output("preseq", path);
            }
            Err(e) => ui.warn(&format!("preseq: {e:#}")),
        }
    }

    Ok(summary::InputSummary {
        bam_file: bam_path.to_string(),
        status: "success".to_string(),
        error: None,
        runtime_seconds: 0.0,
        counting: None,
        dupradar: None,
        dna: Some(dna_summary(&result, &bam_stat_result, &thresholds)),
        outputs: written,
    })
}

/// Build the JSON summary block for a `dna` run.
fn dna_summary(
    result: &rustqc::dna::mosdepth::MosdepthResult,
    bam_stat: &rustqc::common::bam_stat::BamStatResult,
    thresholds: &[u32],
) -> summary::DnaSummary {
    let genome_length = result.total_length();
    let histogram =
        rustqc::dna::mosdepth::merge_histograms(result.contigs.iter().map(|c| &c.histogram));

    let coverage_thresholds = thresholds
        .iter()
        .map(|threshold| {
            let at_or_above: u64 = histogram
                .iter()
                .filter(|(depth, _)| *depth >= threshold)
                .map(|(_, count)| count)
                .sum();
            summary::CoverageThreshold {
                threshold: *threshold,
                pct_bases: if genome_length == 0 {
                    0.0
                } else {
                    at_or_above as f64 * 100.0 / genome_length as f64
                },
            }
        })
        .collect();

    // Median: walk the depth histogram until half the reference is behind us.
    let mut seen = 0u64;
    let mut median = 0u32;
    for (depth, count) in &histogram {
        seen += count;
        if seen * 2 >= genome_length {
            median = *depth;
            break;
        }
    }

    let duplicate_pct = if bam_stat.total_records == 0 {
        0.0
    } else {
        bam_stat.duplicates as f64 * 100.0 / bam_stat.total_records as f64
    };

    summary::DnaSummary {
        genome_length,
        covered_bases: result.total_bases(),
        mean_coverage: result.mean(),
        median_coverage: median,
        max_coverage: result.max(),
        coverage_thresholds,
        total_reads: bam_stat.total_records,
        duplicates: bam_stat.duplicates,
        duplicate_pct,
    }
}

/// Count the reference's non-N bases, which is Picard's `GENOME_TERRITORY`.
///
/// The whole reference is read once. Picard does the same, and the figure
/// cannot be taken from the alignment header, which records contig lengths
/// including their N runs.
fn genome_territory(reference: &str) -> Result<u64> {
    use std::io::BufRead;

    let reader = rustqc::io::open_reader(reference)
        .with_context(|| format!("Failed to open reference FASTA: {reference}"))?;
    let mut territory = 0u64;
    for line in reader.lines() {
        let line = line.with_context(|| format!("Failed to read reference FASTA: {reference}"))?;
        if line.starts_with('>') {
            continue;
        }
        territory += line.bytes().filter(|b| !matches!(b, b'N' | b'n')).count() as u64;
    }
    Ok(territory)
}

/// How many contig depth arrays may be live at once.
///
/// Each worker holds four bytes per base of its contig, so the largest contig
/// sets the per-worker cost: about 1 GB for GRCh38 chr1. There is no portable
/// way to ask the operating system how much memory is free, so the budget is a
/// fixed 4 GB unless the user overrides it with `--max-depth-workers`.
fn depth_worker_budget(threads: usize, override_value: Option<usize>, largest: u64) -> usize {
    const BUDGET_BYTES: u64 = 4 * 1024 * 1024 * 1024;
    if let Some(value) = override_value {
        return value.max(1);
    }
    let per_worker = largest.saturating_mul(4).max(1);
    let affordable = (BUDGET_BYTES / per_worker).max(1) as usize;
    threads.min(affordable).max(1)
}

/// Reconstruct the command line for the featureCounts-compatible header comment.
fn reconstruct_command_line(args: &cli::RnaArgs) -> String {
    let mut parts = vec![format!(
        "rustqc rna {}",
        args.input
            .iter()
            .map(|s| shell_escape(s))
            .collect::<Vec<_>>()
            .join(" "),
    )];
    parts.push(format!("--gtf {}", shell_escape(&args.gtf)));
    if let Some(s) = args.stranded {
        parts.push(format!("-s {}", s));
    }
    if args.paired {
        parts.push("-p".to_string());
    }
    if args.threads != 1 {
        parts.push(format!("-t {}", args.threads));
    }
    if args.outdir != "." {
        parts.push(format!("-o {}", shell_escape(&args.outdir)));
    }
    if let Some(ref config_path) = args.config {
        parts.push(format!("-c {}", shell_escape(config_path)));
    }
    if let Some(ref biotype) = args.biotype_attribute {
        parts.push(format!("--biotype-attribute {}", shell_escape(biotype)));
    }
    if let Some(ref reference) = args.reference {
        parts.push(format!("--reference {}", shell_escape(reference)));
    }
    if args.flat_output {
        parts.push("--flat-output".to_string());
    }
    if args.skip_dup_check {
        parts.push("--skip-dup-check".to_string());
    }
    parts.join(" ")
}

/// Shell-escaping: wrap in single quotes if the string contains shell metacharacters.
///
/// Single quotes prevent all shell interpretation. Any embedded single quotes
/// are escaped using the `'\''` pattern (end quote, escaped quote, restart quote).
fn shell_escape(s: &str) -> String {
    if s.contains(|c: char| c.is_whitespace() || "\"'\\$`!#&|;(){}[]<>?*~".contains(c)) {
        format!("'{}'", s.replace('\'', "'\\''"))
    } else {
        s.to_string()
    }
}

/// Run the full RNA-Seq QC pipeline: dupRadar + featureCounts + RSeQC analyses.
///
/// When multiple BAM files are provided, they are processed in parallel using
/// rayon. The GTF annotation is parsed once and shared across all BAM files.
/// Available threads are distributed across the parallel BAM processing jobs.
fn run_rna(args: cli::RnaArgs, ui: &Ui) -> Result<()> {
    // Validate thread count
    ensure!(
        args.threads >= 1,
        "--threads must be at least 1 (got {})",
        args.threads
    );

    // Load and merge configuration files from all sources:
    // XDG system → XDG user → RUSTQC_CONFIG env → explicit -c flag.
    // Each layer overrides only the leaf fields it sets.
    let (cfg, config_sources) = config::load_merged_config(args.config.as_deref())?;
    let mut config = cfg.rna;

    if !config_sources.is_empty() {
        for (path, source) in &config_sources {
            ui.detail(&format!("Loaded config: {} ({})", path.display(), source));
        }
        if config.has_chromosome_mapping() {
            ui.detail(&format!(
                "Chromosome name mapping: {} entries",
                config.chromosome_mapping.len()
            ));
        }
    }

    // Apply CLI overrides to skip flags
    if args.skip_tin {
        config.tin.enabled = false;
    }
    if args.skip_read_duplication {
        config.read_duplication.enabled = false;
    }
    if args.skip_preseq {
        config.preseq.enabled = false;
    }
    if let Some(val) = args.preseq_max_extrap {
        config.preseq.max_extrap = val;
    }
    if let Some(val) = args.preseq_step_size {
        config.preseq.step_size = val;
    }
    if let Some(val) = args.preseq_n_bootstraps {
        config.preseq.n_bootstraps = val;
    }
    if let Some(val) = args.preseq_seg_len {
        config.preseq.max_segment_length = val;
    }

    // Apply CLI overrides for RSeQC tool parameters
    if let Some(val) = args.infer_experiment_sample_size {
        config.infer_experiment.sample_size = Some(val);
    }
    if let Some(val) = args.min_intron {
        config.junction_annotation.min_intron = Some(val);
    }
    if let Some(val) = args.junction_saturation_min_coverage {
        config.junction_saturation.min_coverage = Some(val);
    }
    if let Some(val) = args.junction_saturation_percentile_floor {
        config.junction_saturation.percentile_floor = Some(val);
    }
    if let Some(val) = args.junction_saturation_percentile_ceiling {
        config.junction_saturation.percentile_ceiling = Some(val);
    }
    if let Some(val) = args.junction_saturation_percentile_step {
        config.junction_saturation.percentile_step = Some(val);
    }
    if let Some(val) = args.inner_distance_sample_size {
        config.inner_distance.sample_size = Some(val);
    }
    if let Some(val) = args.inner_distance_lower_bound {
        config.inner_distance.lower_bound = Some(val);
    }
    if let Some(val) = args.inner_distance_upper_bound {
        config.inner_distance.upper_bound = Some(val);
    }
    if let Some(val) = args.inner_distance_step {
        config.inner_distance.step = Some(val);
    }

    // Apply per-tool seed overrides from CLI flags
    if let Some(seed) = args.preseq_seed {
        config.preseq.seed = seed;
    }
    if let Some(seed) = args.tin_seed {
        config.tin.seed = Some(seed);
    }
    if let Some(seed) = args.junction_saturation_seed {
        config.junction_saturation.seed = Some(seed);
    }

    // Warn early if CRAM input is likely but no reference is provided
    if args.input.iter().any(|f| f.ends_with(".cram"))
        && args.reference.is_none()
        && std::env::var("REF_PATH").is_err()
        && std::env::var("REF_CACHE").is_err()
    {
        anyhow::bail!(
            "CRAM input requires a reference FASTA. \
             Pass --reference <FASTA> or set the REF_PATH environment variable."
        );
    }

    // Validate all input alignment files before expensive GTF parsing
    for bam_path in &args.input {
        let mut reader = rust_htslib::bam::Reader::from_path(bam_path)
            .with_context(|| format!("Cannot open alignment file '{}'", bam_path))?;
        if let Some(ref reference) = args.reference {
            reader
                .set_reference(reference)
                .with_context(|| format!("Cannot set reference for CRAM file '{}'", bam_path))?;
        }
        let _header = reader.header().clone();
        // Reader dropped — just validating the file is openable
    }

    let start = Instant::now();
    let start_time = format_utc_now();
    let n_bams = args.input.len();

    // --sample-name (or config sample_name) only makes sense for a single BAM
    let effective_sample_name = args
        .sample_name
        .as_deref()
        .or(config.sample_name.as_deref());
    if n_bams > 1 && effective_sample_name.is_some() {
        let source = if args.sample_name.is_some() {
            "--sample-name flag"
        } else {
            "config file sample_name"
        };
        anyhow::bail!(
            "{source} cannot be used with multiple input files \
             (would produce identical output filenames)"
        );
    }

    let cpu_target = cpu::binary_target();
    let cpu_features = cpu::detected_features();
    let cpu_info = cpu::cpu_info_line();
    ui.header(
        env!("CARGO_PKG_VERSION"),
        env!("GIT_SHORT_HASH"),
        env!("BUILD_TIMESTAMP"),
        Some(&cpu_info),
    );

    // Resolve effective stranded/paired early for display (config loaded above)
    let effective_stranded = args
        .stranded
        .or(config.stranded)
        .unwrap_or(rustqc::Strandedness::Unstranded);
    let effective_paired = args.paired || config.paired.unwrap_or(false);

    if n_bams == 1 {
        ui.config("Input", &args.input[0]);
    } else {
        ui.config("Input", &format!("{} BAM files", n_bams));
        for f in &args.input {
            ui.detail(&format!("  {f}"));
        }
    }
    ui.config("Annotation", &args.gtf);
    if let Some(ref reference) = args.reference {
        ui.config("Reference", reference);
    }
    ui.config("Stranded", &effective_stranded.to_string());
    ui.config("Paired", &effective_paired.to_string());
    ui.config("CPU Threads", &args.threads.to_string());
    let outdir_display = if std::path::Path::new(&args.outdir).is_relative() {
        format!("./{}", args.outdir)
    } else {
        args.outdir.clone()
    };
    ui.config("Output dir", &outdir_display);

    // Determine biotype attribute name (CLI overrides config, with auto-detection fallback)
    let configured_biotype = args
        .biotype_attribute
        .clone()
        .unwrap_or_else(|| config.featurecounts.biotype_attribute.clone());

    // Determine which extra GTF attributes we need, and parse GTF if provided
    let mut extra_attributes: Vec<String> = Vec::new();
    let mut biotype_attribute = configured_biotype.clone();
    let need_biotype = config.any_biotype_output();

    // Detect biotype attributes in GTF
    let gtf_path = &args.gtf;
    if need_biotype {
        // Check if the configured biotype attribute exists in the GTF.
        // If not found and the user didn't explicitly set it, try common alternatives:
        // Ensembl GTFs use "gene_biotype", GENCODE GTFs use "gene_type".
        let user_explicit = args.biotype_attribute.is_some();
        if gtf::attribute_exists_in_gtf(gtf_path, &biotype_attribute, 1000) {
            extra_attributes.push(biotype_attribute.clone());
            ui.detail(&format!("Biotype attribute: {}", biotype_attribute));
        } else if !user_explicit {
            // Auto-detect: try known alternatives
            let alternatives = if biotype_attribute == "gene_biotype" {
                vec!["gene_type"]
            } else if biotype_attribute == "gene_type" {
                vec!["gene_biotype"]
            } else {
                vec!["gene_biotype", "gene_type"]
            };
            let mut found = false;
            for alt in &alternatives {
                if gtf::attribute_exists_in_gtf(gtf_path, alt, 1000) {
                    ui.detail(&format!(
                        "Biotype attribute '{}' not found, using '{}'",
                        biotype_attribute, alt
                    ));
                    biotype_attribute = alt.to_string();
                    extra_attributes.push(biotype_attribute.clone());
                    found = true;
                    break;
                }
            }
            if !found {
                let tried: Vec<_> = std::iter::once(configured_biotype.as_str())
                    .chain(alternatives.iter().copied())
                    .collect();
                let names = tried
                    .iter()
                    .map(|a| format!("'{a}'"))
                    .collect::<Vec<_>>()
                    .join(" and ");
                ui.warn(&format!(
                    "Biotype attributes {} not found in GTF, skipping biotype outputs \
                         (use --biotype-attribute to specify)",
                    names
                ));
            }
        } else {
            ui.warn(&format!(
                "Biotype attribute '{}' not found in GTF, skipping biotype outputs",
                biotype_attribute
            ));
        }
    }

    // Step 1: Parse GTF annotation (shared across all BAM files)
    ui.blank();
    ui.step("Parsing GTF annotation...");
    let gtf_start = Instant::now();
    let genes = gtf::parse_gtf(gtf_path, &extra_attributes)?;
    ui.detail(&format!(
        "Parsed {} genes in {}",
        format_count(genes.len() as u64),
        format_duration(gtf_start.elapsed()),
    ));

    // Build RSeQC data structures from GTF annotation.
    // These are built once and shared across all BAM files.
    // Each tool's data is only built when enabled in the config.
    let gene_model = if config.infer_experiment.enabled {
        ui.detail("Building gene model for infer_experiment...");
        Some(rna::rseqc::infer_experiment::GeneModel::from_genes(&genes))
    } else {
        None
    };

    let ref_junctions = if config.junction_annotation.enabled {
        ui.detail("Building reference junctions...");
        Some(rna::rseqc::common::build_reference_junctions_from_genes(
            &genes,
        ))
    } else {
        None
    };

    let known_junctions = if config.junction_saturation.enabled {
        ui.detail("Building known junction set...");
        Some(rna::rseqc::common::build_known_junctions_from_genes(&genes))
    } else {
        None
    };

    let rd_regions = if config.read_distribution.enabled {
        ui.detail("Building genomic region sets...");
        Some(rna::rseqc::read_distribution::build_regions_from_genes(
            &genes,
        ))
    } else {
        None
    };

    let exon_bitset = if config.inner_distance.enabled {
        ui.detail("Building exon bitset...");
        Some(rna::rseqc::inner_distance::ExonBitset::from_genes(&genes))
    } else {
        None
    };

    let transcript_tree = if config.inner_distance.enabled {
        ui.detail("Building transcript tree...");
        Some(rna::rseqc::inner_distance::TranscriptTree::from_genes(
            &genes,
        ))
    } else {
        None
    };

    let tin_sample_size = config.tin.sample_size.unwrap_or(100) as usize;
    let tin_index = if config.tin.enabled {
        ui.detail("Building TIN index...");
        Some(rna::rseqc::tin::TinIndex::from_genes(
            &genes,
            tin_sample_size,
        ))
    } else {
        None
    };

    let chrom_mapping = config.alignment_to_gtf_mapping();
    let chrom_prefix = config.chromosome_prefix().map(|s| s.to_owned());

    // Reconstruct command line for featureCounts-compatible header
    let command_line = reconstruct_command_line(&args);

    // Validate that input file stems are unique (otherwise outputs would collide).
    if n_bams > 1 {
        let mut seen_stems = HashSet::new();
        for bam_path in &args.input {
            let stem = Path::new(bam_path)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or(bam_path);
            anyhow::ensure!(
                seen_stems.insert(stem.to_owned()),
                "Duplicate BAM file stem '{}': multiple BAM files with the same \
                 filename would produce conflicting output files. Rename or \
                 reorganise input files so each has a unique filename.",
                stem
            );
        }
    }

    // Determine thread allocation for parallel BAM processing.
    // When processing multiple BAMs, we run BAMs in parallel and divide threads
    // among them. Each BAM's count_reads() creates its own rayon pool internally.
    // Note: the outer pool threads are mostly blocked waiting on inner pools, so
    // actual CPU-active threads stay close to `args.threads`. However, the total
    // OS thread count may briefly exceed `--threads` due to the outer pool threads
    // and temporary plot-generation threads (3 per BAM via std::thread::scope).
    // Integer division may leave up to `n_parallel - 1` threads unused.
    let n_parallel = n_bams.min(args.threads).max(1);
    let threads_per_bam = (args.threads / n_parallel).max(1);
    if n_bams > 1 {
        ui.detail(&format!(
            "Processing {} BAM files ({} in parallel, {} threads each)",
            n_bams, n_parallel, threads_per_bam
        ));
    }

    // Create output directory
    let outdir = Path::new(&args.outdir);
    std::fs::create_dir_all(outdir)
        .with_context(|| format!("Failed to create output directory: {}", outdir.display()))?;

    // Determine if biotype attribute was found in the GTF
    let biotype_in_gtf = extra_attributes.contains(&biotype_attribute);

    // Effective flat_output: true if enabled by either CLI flag or config file
    let flat_output = args.flat_output || config.flat_output;

    // Build the shared parameters struct for process_single_bam
    let shared = SharedParams {
        ui,
        stranded: effective_stranded,
        paired: effective_paired,
        chrom_mapping: &chrom_mapping,
        chrom_prefix: chrom_prefix.as_deref(),
        outdir,
        flat_output,
        reference: args.reference.as_deref(),
        skip_dup_check: args.skip_dup_check,
        config: &config,
        biotype_attribute: &biotype_attribute,
        biotype_in_gtf,
        command_line: &command_line,
        gene_model: gene_model.as_ref(),
        ref_junctions: ref_junctions.as_ref(),
        known_junctions: known_junctions.as_ref(),
        rd_regions: rd_regions.as_ref(),
        exon_bitset: exon_bitset.as_ref(),
        transcript_tree: transcript_tree.as_ref(),
        mapq_cut: args.mapq_cut,
        infer_experiment_sample_size: config.infer_experiment.sample_size.unwrap_or(200_000),
        min_intron: config.junction_annotation.min_intron.unwrap_or(50),
        junction_saturation_min_coverage: config.junction_saturation.min_coverage.unwrap_or(1),
        junction_saturation_percentile_floor: config
            .junction_saturation
            .percentile_floor
            .unwrap_or(5),
        junction_saturation_percentile_ceiling: config
            .junction_saturation
            .percentile_ceiling
            .unwrap_or(100),
        junction_saturation_percentile_step: config
            .junction_saturation
            .percentile_step
            .unwrap_or(5),
        junction_saturation_seed: config.junction_saturation.seed.unwrap_or(42),
        inner_distance_sample_size: config.inner_distance.sample_size.unwrap_or(1_000_000),
        inner_distance_lower_bound: config.inner_distance.lower_bound.unwrap_or(-250),
        inner_distance_upper_bound: config.inner_distance.upper_bound.unwrap_or(250),
        inner_distance_step: config.inner_distance.step.unwrap_or(5),
        tin_index: tin_index.as_ref(),
        tin_sample_size,
        tin_min_coverage: config.tin.min_coverage.unwrap_or(10),
        gtf_path: &args.gtf,
        sample_name_override: effective_sample_name,
    };

    // Step 2: Process all alignment files (in parallel when multiple)
    let bam_results: Vec<(String, Result<BamResult>)> = if n_bams == 1 {
        // Single file: use all threads directly, no outer rayon pool needed
        vec![(
            args.input[0].clone(),
            process_single_bam(&args.input[0], &genes, args.threads, &shared),
        )]
    } else {
        // Multiple files: process in parallel with a dedicated rayon pool
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(n_parallel)
            .build()
            .context("Failed to create rayon thread pool for parallel BAM processing")?;

        pool.install(|| {
            args.input
                .par_iter()
                .map(|bam_path| {
                    (
                        bam_path.clone(),
                        process_single_bam(bam_path, &genes, threads_per_bam, &shared),
                    )
                })
                .collect()
        })
    };

    // Collect results for summary
    let mut n_err = 0;
    let mut input_summaries: Vec<summary::InputSummary> = Vec::new();

    for (bam_path, result) in &bam_results {
        match result {
            Ok(bam_result) => {
                input_summaries.push(bam_result.to_input_summary(bam_path));
            }
            Err(e) => {
                n_err += 1;
                ui.error(&format!("Failed to process {}: {:#}", bam_path, e));
                input_summaries.push(summary::InputSummary {
                    bam_file: bam_path.clone(),
                    status: "failed".to_string(),
                    error: Some(format!("{:#}", e)),
                    runtime_seconds: 0.0,
                    counting: None,
                    dupradar: None,
                    dna: None,
                    outputs: vec![],
                });
            }
        }
    }

    // Multi-BAM summary
    if n_bams > 1 {
        ui.blank();
        ui.section(&format!(
            "Processed {} files in {}",
            n_bams,
            format_duration(start.elapsed()),
        ));
        for (bam_path, result) in &bam_results {
            let name = Path::new(bam_path)
                .file_name()
                .and_then(|s| s.to_str())
                .unwrap_or(bam_path);
            match result {
                Ok(bam_result) => {
                    ui.bam_result_ok(name, bam_result.duration);
                }
                Err(e) => {
                    let msg = format!("{:#}", e);
                    // Truncate long error messages for the summary line
                    let short = console::truncate_str(&msg, 60, "…");
                    ui.bam_result_err(name, &short);
                }
            }
        }
    }

    let elapsed = start.elapsed();
    let end_time = format_utc_now();

    // Write JSON summary if requested
    if let Some(ref json_path) = args.json_summary {
        let run_summary = summary::RunSummary {
            version: env!("CARGO_PKG_VERSION").to_string(),
            commit: env!("GIT_SHORT_HASH").to_string(),
            binary_target: cpu_target.to_string(),
            cpu_features: cpu_features.iter().map(|s| s.to_string()).collect(),
            timestamp_start: start_time.clone(),
            timestamp_end: end_time.clone(),
            runtime_seconds: elapsed.as_secs_f64(),
            inputs: input_summaries,
        };
        let json = serde_json::to_string_pretty(&run_summary)
            .context("Failed to serialize JSON summary")?;
        if json_path == "-" {
            println!("{json}");
            ui.detail("JSON summary written to stdout");
        } else {
            let path = if json_path.is_empty() {
                outdir.join("rustqc_summary.json")
            } else {
                Path::new(json_path).to_path_buf()
            };
            std::fs::write(&path, &json)
                .with_context(|| format!("Failed to write JSON summary to {}", path.display()))?;
            ui.output_item("JSON summary", &path.display().to_string());
        }
    }

    // Write citations file
    let citations_path = outdir.join("CITATIONS.md");
    citations::write_citations(
        &citations_path,
        &config,
        env!("CARGO_PKG_VERSION"),
        env!("GIT_SHORT_HASH"),
    )?;
    ui.output_item("Citations", &citations_path.display().to_string());

    // Check for strandedness mismatch between user-specified and inferred values
    for (bam_path, result) in &bam_results {
        if let Ok(bam_result) = result {
            if let Some(ref ie_result) = bam_result.infer_experiment {
                if let Some((inferred, suggestion)) =
                    rna::rseqc::infer_experiment::check_strandedness_mismatch(
                        ie_result,
                        effective_stranded,
                    )
                {
                    let bam_name = Path::new(bam_path)
                        .file_name()
                        .and_then(|s| s.to_str())
                        .unwrap_or(bam_path);
                    let line1 = "Strandedness mismatch detected!".to_string();
                    let line2 = format!(
                        "{} - you specified '--stranded {}' but infer_experiment suggests '{}'",
                        bam_name, effective_stranded, inferred,
                    );
                    let line3 = format!(
                        "(forward fraction: {:.4}, reverse fraction: {:.4})",
                        ie_result.frac_protocol1, ie_result.frac_protocol2,
                    );
                    let line4 = format!("Consider re-running with '--stranded {}'", suggestion,);
                    ui.blank();
                    ui.warn_box(&[&line1, &line2, &line3, &line4]);
                }
            }
        }
    }

    ui.finish("RustQC run finished", elapsed);

    if n_err > 0 {
        anyhow::bail!("{} file(s) failed to process", n_err);
    }

    Ok(())
}

// ============================================================================
// Shared parameters
// ============================================================================

/// Parameters shared across all BAM files in a run.
///
/// Bundles the read-only configuration, annotation data, and tool parameters
/// that are computed once in `run_rna()` and passed to each `process_single_bam()`.
/// This avoids a long parameter list on the processing function.
struct SharedParams<'a> {
    /// Terminal UI handle.
    ui: &'a Ui,
    /// Library strandedness.
    stranded: rustqc::Strandedness,
    /// Whether the library is paired-end.
    paired: bool,
    /// Alignment-to-GTF chromosome name mapping.
    chrom_mapping: &'a HashMap<String, String>,
    /// Optional chromosome name prefix.
    chrom_prefix: Option<&'a str>,
    /// Output directory for results.
    outdir: &'a Path,
    /// When true, write all files directly to outdir (no subdirectories).
    flat_output: bool,
    /// Optional override for the sample name used in output filenames.
    sample_name_override: Option<&'a str>,
    /// Optional reference FASTA for CRAM files.
    reference: Option<&'a str>,
    /// Whether to skip duplicate-marking validation.
    skip_dup_check: bool,
    /// Configuration for conditional outputs.
    config: &'a config::RnaConfig,
    /// GTF attribute name for biotype counting.
    biotype_attribute: &'a str,
    /// Whether the biotype attribute was found in the GTF.
    biotype_in_gtf: bool,
    /// Reconstructed command line for featureCounts header.
    command_line: &'a str,
    /// Pre-built gene model for infer_experiment (from GTF).
    gene_model: Option<&'a rna::rseqc::infer_experiment::GeneModel>,
    /// Pre-built reference junctions for junction_annotation (from GTF).
    ref_junctions: Option<&'a rna::rseqc::common::ReferenceJunctions>,
    /// Pre-built known junction set for junction_saturation (from GTF).
    known_junctions: Option<&'a rna::rseqc::common::KnownJunctionSet>,
    /// Pre-built genomic region sets for read_distribution (from GTF).
    rd_regions: Option<&'a rna::rseqc::read_distribution::RegionSets>,
    /// Pre-built exon bitset for inner_distance (from GTF).
    exon_bitset: Option<&'a rna::rseqc::inner_distance::ExonBitset>,
    /// Pre-built transcript tree for inner_distance (from GTF).
    transcript_tree: Option<&'a rna::rseqc::inner_distance::TranscriptTree>,
    /// MAPQ cutoff for read quality filtering.
    mapq_cut: u8,
    /// Maximum reads to sample for strandedness inference.
    infer_experiment_sample_size: u64,
    /// Minimum intron size for junction filtering.
    min_intron: u64,
    /// Minimum coverage for junction saturation.
    junction_saturation_min_coverage: u64,
    /// Sampling start percentage for junction saturation.
    junction_saturation_percentile_floor: u64,
    /// Sampling end percentage for junction saturation.
    junction_saturation_percentile_ceiling: u64,
    /// Sampling step percentage for junction saturation.
    junction_saturation_percentile_step: u64,
    /// Random seed for junction saturation shuffle.
    junction_saturation_seed: u64,
    /// Maximum read pairs to sample for inner distance.
    inner_distance_sample_size: u64,
    /// Lower bound of inner distance histogram.
    inner_distance_lower_bound: i64,
    /// Upper bound of inner distance histogram.
    inner_distance_upper_bound: i64,
    /// Bin width for inner distance histogram.
    inner_distance_step: i64,
    /// Pre-built TIN index for transcript integrity analysis (from GTF).
    tin_index: Option<&'a rna::rseqc::tin::TinIndex>,
    /// Number of equally-spaced positions to sample per transcript for TIN.
    tin_sample_size: usize,
    /// Minimum read-start count per transcript to compute TIN.
    tin_min_coverage: u32,
    /// Path to GTF file (for Qualimap report output).
    gtf_path: &'a str,
}

// ============================================================================
// Per-BAM result
// ============================================================================

/// Collected results from processing a single BAM file, used for the summary
/// box display and JSON output.
#[derive(Debug, Default)]
struct BamResult {
    /// Wall-clock processing time.
    duration: std::time::Duration,
    /// Total reads in the BAM.
    total_reads: u64,
    /// Total fragments (pairs or single reads).
    total_fragments: u64,
    /// Mapped reads.
    mapped_reads: u64,
    /// Fragments assigned to a gene.
    assigned: u64,
    /// Fragments with no overlapping gene.
    no_features: u64,
    /// Fragments overlapping multiple genes.
    ambiguous: u64,
    /// Duplicate-flagged reads.
    duplicates: u64,
    /// Multimapping reads.
    multimappers: u64,
    /// dupRadar fit intercept (if available).
    dupradar_intercept: Option<f64>,
    /// dupRadar fit slope (if available).
    dupradar_slope: Option<f64>,
    /// dupRadar gene stats.
    dupradar_total_genes: u64,
    /// Genes with reads.
    dupradar_genes_with_reads: u64,
    /// Genes with duplicates.
    dupradar_genes_with_dups: u64,
    /// Output files written (tool, path).
    outputs: Vec<(String, String)>,
    /// infer_experiment result (if the tool was enabled).
    infer_experiment: Option<rna::rseqc::infer_experiment::InferExperimentResult>,
}

impl BamResult {
    /// Convert to a JSON-serializable InputSummary.
    fn to_input_summary(&self, bam_path: &str) -> summary::InputSummary {
        let counting = Some(summary::CountingSummary {
            total_reads: self.total_reads,
            mapped_reads: self.mapped_reads,
            fragments: self.total_fragments,
            assigned: self.assigned,
            no_features: self.no_features,
            ambiguous: self.ambiguous,
            duplicates: self.duplicates,
            multimappers: self.multimappers,
            assigned_pct: if self.total_fragments > 0 {
                self.assigned as f64 / self.total_fragments as f64 * 100.0
            } else {
                0.0
            },
            duplicate_pct: if self.mapped_reads > 0 {
                self.duplicates as f64 / self.mapped_reads as f64 * 100.0
            } else {
                0.0
            },
        });

        let dupradar = if self.dupradar_total_genes > 0 {
            Some(summary::DupradarSummary {
                total_genes: self.dupradar_total_genes,
                genes_with_reads: self.dupradar_genes_with_reads,
                genes_with_duplication: self.dupradar_genes_with_dups,
                intercept: self.dupradar_intercept,
                slope: self.dupradar_slope,
            })
        } else {
            None
        };

        summary::InputSummary {
            bam_file: bam_path.to_string(),
            status: "success".to_string(),
            error: None,
            runtime_seconds: self.duration.as_secs_f64(),
            counting,
            dupradar,
            dna: None,
            outputs: self
                .outputs
                .iter()
                .map(|(tool, path)| summary::OutputFile {
                    tool: tool.clone(),
                    path: path.clone(),
                })
                .collect(),
        }
    }
}

// ============================================================================
// Per-BAM processing
// ============================================================================

/// Process a single alignment file through the full analysis pipeline.
///
/// This runs the complete analysis for one file: counting, featureCounts output,
/// biotype counting, duplication matrix, model fitting, plotting, MultiQC output,
/// and all enabled RSeQC analyses.
///
/// # Arguments
///
/// * `bam_path` - Path to the duplicate-marked alignment file (SAM/BAM/CRAM)
/// * `genes` - Parsed GTF gene annotations
/// * `threads` - Number of threads for this file's read counting
/// * `params` - Shared parameters (config, annotations, tool settings)
fn process_single_bam(
    bam_path: &str,
    genes: &IndexMap<String, gtf::Gene>,
    threads: usize,
    params: &SharedParams,
) -> Result<BamResult> {
    let ui = params.ui;
    let bam_stem = Path::new(bam_path)
        .file_stem()
        .context("Input path has no filename")?
        .to_str()
        .context("Input filename is not valid UTF-8")?;
    let sample_name = params
        .sample_name_override
        .map(|s| s.to_string())
        .unwrap_or_else(|| bam_stem.to_string());

    let bam_start = Instant::now();
    let config = params.config;
    let outdir = params.outdir;

    // Track output files for the summary
    let mut written_outputs: Vec<(String, String)> = Vec::new();

    // === Build RSeQC config and annotations ===
    let rseqc_config = RseqcConfig {
        mapq_cut: params.mapq_cut,
        infer_experiment_sample_size: params.infer_experiment_sample_size,
        min_intron: params.min_intron,
        junction_saturation_min_coverage: params.junction_saturation_min_coverage as u32,
        junction_saturation_sample_start: params.junction_saturation_percentile_floor as u32,
        junction_saturation_sample_end: params.junction_saturation_percentile_ceiling as u32,
        junction_saturation_sample_step: params.junction_saturation_percentile_step as u32,
        inner_distance_sample_size: params.inner_distance_sample_size,
        inner_distance_lower_bound: params.inner_distance_lower_bound,
        inner_distance_upper_bound: params.inner_distance_upper_bound,
        inner_distance_step: params.inner_distance_step,
        bam_stat_enabled: config.bam_stat.enabled
            || config.flagstat.enabled
            || config.idxstats.enabled
            || config.samtools_stats.enabled,
        infer_experiment_enabled: config.infer_experiment.enabled && params.gene_model.is_some(),
        read_duplication_enabled: config.read_duplication.enabled,
        read_distribution_enabled: config.read_distribution.enabled && params.rd_regions.is_some(),
        junction_annotation_enabled: config.junction_annotation.enabled
            && params.ref_junctions.is_some(),
        junction_saturation_enabled: config.junction_saturation.enabled
            && params.known_junctions.is_some(),
        inner_distance_enabled: config.inner_distance.enabled
            && params.exon_bitset.is_some()
            && params.transcript_tree.is_some(),
        tin_enabled: config.tin.enabled && params.tin_index.is_some(),
        tin_sample_size: params.tin_sample_size,
        tin_min_coverage: params.tin_min_coverage,
        tin_seed: config.tin.seed,
        junction_saturation_seed: config.junction_saturation.seed.unwrap_or(42),
        preseq_enabled: config.preseq.enabled,
        preseq_max_segment_length: config.preseq.max_segment_length,
    };

    let rseqc_annotations = RseqcAnnotations {
        gene_model: params.gene_model,
        ref_junctions: params.ref_junctions,
        rd_regions: params.rd_regions,
        exon_bitset: params.exon_bitset,
        transcript_tree: params.transcript_tree,
        tin_index: params.tin_index,
    };

    let any_rseqc_enabled = rseqc_config.bam_stat_enabled
        || rseqc_config.infer_experiment_enabled
        || rseqc_config.read_duplication_enabled
        || rseqc_config.read_distribution_enabled
        || rseqc_config.junction_annotation_enabled
        || rseqc_config.junction_saturation_enabled
        || rseqc_config.inner_distance_enabled
        || rseqc_config.preseq_enabled
        || rseqc_config.tin_enabled;

    // === Build Qualimap exon index (if enabled) ===
    let qualimap_index = if params.config.qualimap.enabled {
        Some(rna::qualimap::QualimapIndex::from_genes(genes))
    } else {
        None
    };

    // === dupRadar counting ===
    ui.blank();
    ui.section(&format!("Processing {}", bam_path));
    let pb = ui.progress_bar();
    let count_start = Instant::now();
    let mut count_result = rna::dupradar::counting::count_reads(
        bam_path,
        genes,
        params.stranded,
        params.paired,
        threads,
        params.chrom_mapping,
        params.chrom_prefix,
        params.reference,
        params.skip_dup_check,
        params.biotype_attribute,
        if any_rseqc_enabled {
            Some(&rseqc_config)
        } else {
            None
        },
        if any_rseqc_enabled {
            Some(&rseqc_annotations)
        } else {
            None
        },
        qualimap_index.as_ref(),
        Some(&pb),
    )?;
    let count_duration = count_start.elapsed();
    ui.finish_progress(&pb, count_result.stat_total_reads, count_duration);

    // Extract RSeQC accumulators from count_result.
    let rseqc_accums = count_result.rseqc.take();

    // Summary stats for the box
    let total_mapped = count_result.stat_total_mapped;
    let total_dup = count_result.stat_total_dup;

    // Summary box
    ui.summary_box(
        &format!("{} — Counting Summary", sample_name),
        &[
            (
                "Total reads:",
                format_count(count_result.stat_total_reads),
                format!(
                    "({} fragments)",
                    format_count(count_result.stat_total_fragments)
                ),
            ),
            (
                "Assigned:",
                format_count(count_result.stat_assigned),
                format_pct(
                    count_result.stat_assigned,
                    count_result.stat_total_fragments,
                ),
            ),
            (
                "No features:",
                format_count(count_result.stat_no_features),
                format_pct(
                    count_result.stat_no_features,
                    count_result.stat_total_fragments,
                ),
            ),
            (
                "Ambiguous:",
                format_count(count_result.stat_ambiguous),
                format_pct(
                    count_result.stat_ambiguous,
                    count_result.stat_total_fragments,
                ),
            ),
            (
                "Duplicates:",
                format_count(total_dup),
                format_pct(total_dup, total_mapped),
            ),
            (
                "Multimappers:",
                format_count(count_result.fc_multimapping),
                format_pct(count_result.fc_multimapping, count_result.stat_total_reads),
            ),
        ],
    );

    // Output directories: nested subdirectories by default, flat if requested
    let fc_dir = if params.flat_output {
        outdir.to_path_buf()
    } else {
        outdir.join("featurecounts")
    };
    let dr_dir = if params.flat_output {
        outdir.to_path_buf()
    } else {
        outdir.join("dupradar")
    };

    // === featureCounts outputs ===
    ui.section("Writing outputs...");
    if config.any_featurecounts_output() {
        std::fs::create_dir_all(&fc_dir).with_context(|| {
            format!(
                "Failed to create featurecounts output directory: {}",
                fc_dir.display()
            )
        })?;

        if config.featurecounts.counts_file {
            let counts_path = fc_dir.join(format!("{}.featureCounts.tsv", sample_name));
            rna::featurecounts::output::write_counts_file(
                &counts_path,
                genes,
                &count_result,
                bam_path,
                params.command_line,
            )?;
            let p = counts_path.display().to_string();
            ui.output_item("featureCounts", &p);
            written_outputs.push(("featureCounts".into(), p));
        }

        if config.featurecounts.summary_file {
            let summary_path = fc_dir.join(format!("{}.featureCounts.tsv.summary", sample_name));
            rna::featurecounts::output::write_summary_file(&summary_path, &count_result, bam_path)?;
            let p = summary_path.display().to_string();
            ui.output_detail(&format!("Summary: {p}"));
            written_outputs.push(("featureCounts summary".into(), p));
        }

        // Biotype outputs (only if attribute was found in GTF)
        if params.biotype_in_gtf && config.any_biotype_output() {
            let biotype_counts =
                rna::featurecounts::output::aggregate_biotype_counts(&count_result);
            ui.detail(&format!(
                "Biotype counting: {} biotypes found",
                biotype_counts.len()
            ));

            if config.featurecounts.biotype_summary_file {
                let biotype_summary_path =
                    fc_dir.join(format!("{}.featureCounts.biotype.tsv.summary", sample_name));
                rna::featurecounts::output::write_biotype_summary_file(
                    &biotype_summary_path,
                    &count_result,
                    bam_path,
                )?;
                let p = biotype_summary_path.display().to_string();
                ui.output_detail(&format!("Biotype summary: {p}"));
                written_outputs.push(("biotype summary".into(), p));
            }

            if config.featurecounts.biotype_counts {
                let biotype_path = fc_dir.join(format!("{}.biotype_counts.tsv", sample_name));
                rna::featurecounts::output::write_biotype_counts(&biotype_path, &biotype_counts)?;
                let p = biotype_path.display().to_string();
                ui.output_detail(&format!("Biotype counts: {p}"));
                written_outputs.push(("biotype counts".into(), p));
            }

            if config.featurecounts.biotype_counts_mqc {
                let mqc_biotype_path =
                    fc_dir.join(format!("{}.biotype_counts_mqc.tsv", sample_name));
                rna::featurecounts::output::write_biotype_counts_mqc(
                    &mqc_biotype_path,
                    &biotype_counts,
                )?;
                let p = mqc_biotype_path.display().to_string();
                ui.output_detail(&format!("Biotype MultiQC: {p}"));
                written_outputs.push(("biotype MultiQC".into(), p));
            }

            if config.featurecounts.biotype_rrna_mqc {
                let mqc_rrna_path =
                    fc_dir.join(format!("{}.biotype_counts_rrna_mqc.tsv", sample_name));
                rna::featurecounts::output::write_biotype_rrna_mqc(
                    &mqc_rrna_path,
                    &biotype_counts,
                    count_result.fc_biotype_assigned,
                    &sample_name,
                )?;
                let p = mqc_rrna_path.display().to_string();
                ui.output_detail(&format!("rRNA MultiQC: {p}"));
                written_outputs.push(("rRNA MultiQC".into(), p));
            }
        }
    }

    // === dupRadar outputs ===
    let mut dupradar_intercept: Option<f64> = None;
    let mut dupradar_slope: Option<f64> = None;
    let mut dupradar_total_genes: u64 = 0;
    let mut dupradar_genes_with_reads: u64 = 0;
    let mut dupradar_genes_with_dups: u64 = 0;

    if config.any_dupradar_output() {
        std::fs::create_dir_all(&dr_dir).with_context(|| {
            format!(
                "Failed to create dupradar output directory: {}",
                dr_dir.display()
            )
        })?;
        let dup_matrix = rna::dupradar::dupmatrix::DupMatrix::build(genes, &count_result);

        let stats = dup_matrix.get_stats();
        dupradar_total_genes = stats.n_regions as u64;
        dupradar_genes_with_reads = stats.n_regions_covered as u64;
        dupradar_genes_with_dups = stats.n_regions_duplication as u64;

        // Write duplication matrix
        if config.dupradar.dup_matrix {
            let matrix_path = dr_dir.join(format!("{}_dupMatrix.txt", sample_name));
            dup_matrix.write_tsv(&matrix_path)?;
            let p = matrix_path.display().to_string();
            written_outputs.push(("dupRadar matrix".into(), p));
        }

        // Fit logistic regression model (needed for intercept/slope, density plot, and MultiQC)
        let need_fit = config.dupradar.intercept_slope
            || config.dupradar.density_scatter_plot
            || config.dupradar.multiqc_intercept
            || config.dupradar.multiqc_curve;

        let fit_ok = if need_fit {
            let rpk_values: Vec<f64> = dup_matrix.rows.iter().map(|r| r.rpk).collect();
            let dup_rate_values: Vec<f64> = dup_matrix.rows.iter().map(|r| r.dup_rate).collect();

            let fit_result = rna::dupradar::fitting::duprate_exp_fit(&rpk_values, &dup_rate_values);
            match &fit_result {
                Ok(fit) => {
                    dupradar_intercept = Some(fit.intercept);
                    dupradar_slope = Some(fit.slope);
                    ui.output_detail(&format!(
                        "Model fit: intercept={:.6}, slope={:.6}",
                        fit.intercept, fit.slope
                    ));
                    if config.dupradar.intercept_slope {
                        let fit_path = dr_dir.join(format!("{}_intercept_slope.txt", sample_name));
                        rna::dupradar::plots::write_intercept_slope(fit, &sample_name, &fit_path)?;
                        let p = fit_path.display().to_string();
                        ui.output_detail(&format!("Fit results: {p}"));
                        written_outputs.push(("dupRadar fit".into(), p));
                    }
                    Some(fit.clone())
                }
                Err(e) => {
                    ui.warn(&format!("Could not fit dupRadar model: {}", e));
                    None
                }
            }
        } else {
            None
        };

        // Generate plots (in parallel — all plots read shared immutable data)
        let any_plot = config.dupradar.density_scatter_plot
            || config.dupradar.boxplot
            || config.dupradar.expression_histogram;

        if any_plot {
            let rpkm_threshold = 0.5;
            let rpkm_threshold_rpk = fit_ok.as_ref().and_then(|_fit| {
                let rpk_values: Vec<f64> = dup_matrix.rows.iter().map(|r| r.rpk).collect();
                let rpkm_values: Vec<f64> = dup_matrix.rows.iter().map(|r| r.rpkm).collect();
                rna::dupradar::fitting::compute_rpkm_threshold_rpk(
                    &rpk_values,
                    &rpkm_values,
                    rpkm_threshold,
                )
            });

            let density_path = dr_dir.join(format!("{}_duprateExpDens.png", sample_name));
            let boxplot_path = dr_dir.join(format!("{}_duprateExpBoxplot.png", sample_name));
            let histogram_path = dr_dir.join(format!("{}_expressionHist.png", sample_name));

            let sample_name_str = sample_name.as_str();
            std::thread::scope(|s| -> Result<()> {
                // Density scatter plot (only if fit succeeded and enabled)
                let density_handle = if config.dupradar.density_scatter_plot {
                    fit_ok.as_ref().map(|fit| {
                        let dm_ref = &dup_matrix;
                        let thresh = rpkm_threshold_rpk;
                        let path = &density_path;
                        s.spawn(move || {
                            rna::dupradar::plots::density_scatter_plot(
                                dm_ref,
                                fit,
                                thresh,
                                rpkm_threshold,
                                sample_name_str,
                                path,
                            )
                        })
                    })
                } else {
                    None
                };

                // Boxplot
                let boxplot_handle = if config.dupradar.boxplot {
                    let dm_ref = &dup_matrix;
                    let path = &boxplot_path;
                    Some(s.spawn(move || {
                        rna::dupradar::plots::duprate_boxplot(dm_ref, sample_name_str, path)
                    }))
                } else {
                    None
                };

                // Histogram
                let histogram_handle = if config.dupradar.expression_histogram {
                    let dm_ref = &dup_matrix;
                    let path = &histogram_path;
                    Some(s.spawn(move || {
                        rna::dupradar::plots::expression_histogram(dm_ref, sample_name_str, path)
                    }))
                } else {
                    None
                };

                // Collect results
                if let Some(handle) = density_handle {
                    handle
                        .join()
                        .map_err(|_| anyhow::anyhow!("density scatter plot thread panicked"))??;
                }
                if let Some(handle) = boxplot_handle {
                    handle
                        .join()
                        .map_err(|_| anyhow::anyhow!("boxplot thread panicked"))??;
                }
                if let Some(handle) = histogram_handle {
                    handle
                        .join()
                        .map_err(|_| anyhow::anyhow!("histogram thread panicked"))??;
                }

                Ok(())
            })?;

            let plots_dir = dr_dir.display().to_string();
            written_outputs.push(("dupRadar plots".into(), plots_dir));
        }

        // Write MultiQC-compatible output files
        if let Some(ref fit) = fit_ok {
            if config.dupradar.multiqc_intercept {
                let mqc_intercept_path =
                    dr_dir.join(format!("{}_dup_intercept_mqc.txt", sample_name));
                rna::dupradar::plots::write_mqc_intercept(fit, &sample_name, &mqc_intercept_path)?;
                written_outputs.push((
                    "dupRadar MultiQC".into(),
                    mqc_intercept_path.display().to_string(),
                ));
            }

            if config.dupradar.multiqc_curve {
                let mqc_curve_path =
                    dr_dir.join(format!("{}_duprateExpDensCurve_mqc.txt", sample_name));
                rna::dupradar::plots::write_mqc_curve(fit, &dup_matrix, &mqc_curve_path)?;
                written_outputs.push((
                    "dupRadar MultiQC curve".into(),
                    mqc_curve_path.display().to_string(),
                ));
            }
        }

        // Consolidated dupRadar output line
        ui.output_item("dupRadar", &format!("{}/*", dr_dir.display()));
        ui.output_detail(&format!(
            "{} genes, {} with reads",
            format_count(stats.n_regions as u64),
            format_count(stats.n_regions_covered as u64),
        ));
        if let (Some(intercept), Some(slope)) = (dupradar_intercept, dupradar_slope) {
            ui.output_detail(&format!(
                "Model fit: intercept={:.6}, slope={:.6}",
                intercept, slope,
            ));
        }
    }
    // === Qualimap RNA-Seq QC output ===
    if let (Some(ref qm_result), Some(ref qm_index)) = (&count_result.qualimap, &qualimap_index) {
        let qm_dir = if params.flat_output {
            outdir.to_path_buf()
        } else {
            outdir.join("qualimap")
        };

        rna::qualimap::output::write_qualimap_results(
            qm_result,
            qm_index,
            bam_path,
            params.gtf_path,
            params.stranded,
            &qm_dir,
            &sample_name,
        )?;
        let p = qm_dir.display().to_string();
        ui.output_item("Qualimap", &format!("{p}/*"));
        written_outputs.push(("Qualimap".into(), p));
    }

    // === RSeQC analyses (post-processing of single-pass accumulators) ===
    let rseqc_accums = rseqc_accums.unwrap_or_else(|| {
        ui.detail("No RSeQC tools enabled, skipping");
        RseqcAccumulators::empty()
    });
    // Extract BAM header info (reference names + lengths) for samtools-compatible outputs
    let bam_header_refs = {
        let reader = rust_htslib::bam::Reader::from_path(bam_path)
            .with_context(|| format!("Failed to open BAM for header: {}", bam_path))?;
        let header = reader.header();
        (0..header.target_count())
            .map(|tid| {
                let name = String::from_utf8_lossy(header.tid2name(tid)).to_string();
                let len = header.target_len(tid).unwrap_or(0);
                (name, len)
            })
            .collect::<Vec<(String, u64)>>()
    };
    let rseqc_outputs = write_rseqc_outputs(
        bam_path,
        &sample_name,
        params,
        rseqc_accums,
        &bam_header_refs,
    )?;
    written_outputs.extend(rseqc_outputs.written);

    let bam_duration = bam_start.elapsed();
    ui.finish(bam_stem, bam_duration);

    Ok(BamResult {
        duration: bam_duration,
        total_reads: count_result.stat_total_reads,
        total_fragments: count_result.stat_total_fragments,
        mapped_reads: total_mapped,
        assigned: count_result.stat_assigned,
        no_features: count_result.stat_no_features,
        ambiguous: count_result.stat_ambiguous,
        duplicates: total_dup,
        multimappers: count_result.fc_multimapping,
        dupradar_intercept,
        dupradar_slope,
        dupradar_total_genes,
        dupradar_genes_with_reads,
        dupradar_genes_with_dups,
        outputs: written_outputs,
        infer_experiment: rseqc_outputs.infer_experiment,
    })
}

// ============================================================================
// RSeQC output writing (post-processing of single-pass accumulators)
// ============================================================================

/// Results returned from `write_rseqc_outputs`, bundling output file paths
/// with the optional infer_experiment result for downstream strandedness checks.
struct RseqcOutputs {
    /// Output files written (tool, path).
    written: Vec<(String, String)>,
    /// infer_experiment result, if the tool was enabled and produced data.
    infer_experiment: Option<rna::rseqc::infer_experiment::InferExperimentResult>,
}

/// Write all RSeQC outputs from the single-pass accumulators.
///
/// Converts accumulated data to tool-specific result types and writes all
/// output files, plots, and summaries. Returns the written output paths and
/// the infer_experiment result (if available) for strandedness mismatch checking.
fn write_rseqc_outputs(
    bam_path: &str,
    sample_name: &str,
    params: &SharedParams,
    accums: RseqcAccumulators,
    bam_header_refs: &[(String, u64)],
) -> Result<RseqcOutputs> {
    let ui = params.ui;
    let outdir = params.outdir;
    let mut written: Vec<(String, String)> = Vec::new();
    let mut infer_experiment_result: Option<rna::rseqc::infer_experiment::InferExperimentResult> =
        None;

    // Build tool-specific output directories: nested subdirectories by default, flat if requested
    let flat = params.flat_output;
    let rseqc_bam_stat_dir = if flat {
        outdir.to_path_buf()
    } else {
        outdir.join("rseqc").join("bam_stat")
    };
    let rseqc_read_dup_dir = if flat {
        outdir.to_path_buf()
    } else {
        outdir.join("rseqc").join("read_duplication")
    };
    let rseqc_infer_exp_dir = if flat {
        outdir.to_path_buf()
    } else {
        outdir.join("rseqc").join("infer_experiment")
    };
    let rseqc_read_dist_dir = if flat {
        outdir.to_path_buf()
    } else {
        outdir.join("rseqc").join("read_distribution")
    };
    let rseqc_junc_annot_dir = if flat {
        outdir.to_path_buf()
    } else {
        outdir.join("rseqc").join("junction_annotation")
    };
    let rseqc_junc_sat_dir = if flat {
        outdir.to_path_buf()
    } else {
        outdir.join("rseqc").join("junction_saturation")
    };
    let rseqc_inner_dist_dir = if flat {
        outdir.to_path_buf()
    } else {
        outdir.join("rseqc").join("inner_distance")
    };

    // --- samtools-compatible outputs (flagstat, idxstats, stats) ---
    let samtools_dir = if flat {
        outdir.to_path_buf()
    } else {
        outdir.join("samtools")
    };

    // Compute bam_stat result once — used by both samtools and RSeQC outputs.
    let bam_stat_result = accums.bam_stat.map(|accum| accum.into_result());

    if let Some(ref result) = bam_stat_result {
        let has_samtools = params.config.flagstat.enabled
            || params.config.idxstats.enabled
            || params.config.samtools_stats.enabled;
        if has_samtools {
            ui.output_group("samtools");
        }

        // --- flagstat ---
        if params.config.flagstat.enabled {
            std::fs::create_dir_all(&samtools_dir)?;
            let flagstat_path = samtools_dir.join(format!("{}.flagstat", sample_name));
            common::samtools::flagstat::write_flagstat(result, &flagstat_path)?;
            let p = flagstat_path.display().to_string();
            ui.output_item("flagstat", &p);
            written.push(("flagstat".into(), p));
        }

        // --- idxstats ---
        if params.config.idxstats.enabled {
            std::fs::create_dir_all(&samtools_dir)?;
            let idxstats_path = samtools_dir.join(format!("{}.idxstats", sample_name));
            common::samtools::idxstats::write_idxstats(result, bam_header_refs, &idxstats_path)?;
            let p = idxstats_path.display().to_string();
            ui.output_item("idxstats", &p);
            written.push(("idxstats".into(), p));
        }

        // --- samtools stats SN ---
        if params.config.samtools_stats.enabled {
            std::fs::create_dir_all(&samtools_dir)?;
            let stats_path = samtools_dir.join(format!("{}.stats", sample_name));
            common::samtools::stats::write_stats(result, &stats_path)?;
            let p = stats_path.display().to_string();
            ui.output_item("stats", &p);
            written.push(("samtools stats".into(), p));
        }
    }

    // --- RSeQC outputs ---
    {
        let has_rseqc = (bam_stat_result.is_some() && params.config.bam_stat.enabled)
            || accums.read_dup.is_some()
            || accums.infer_exp.is_some()
            || accums.read_dist.is_some()
            || accums.junc_annot.is_some()
            || accums.junc_sat.is_some()
            || accums.inner_dist.is_some()
            || accums.tin.is_some();
        if has_rseqc {
            ui.output_group("RSeQC");
        }
    }

    // --- bam_stat ---
    if let Some(ref result) = bam_stat_result {
        if params.config.bam_stat.enabled {
            std::fs::create_dir_all(&rseqc_bam_stat_dir)?;
            let output_path = rseqc_bam_stat_dir.join(format!("{}.bam_stat.txt", sample_name));
            rna::rseqc::bam_stat::write_bam_stat(result, &output_path)?;
            let p = output_path.display().to_string();
            ui.output_item("bam_stat", &p);
            written.push(("bam_stat".into(), p));
        }
    }

    // --- read_duplication ---
    if let Some(accum) = accums.read_dup {
        std::fs::create_dir_all(&rseqc_read_dup_dir)?;
        let result = accum.into_result();
        rna::rseqc::read_duplication::write_read_duplication(
            &result,
            &rseqc_read_dup_dir,
            sample_name,
        )?;
        let plot_path = rseqc_read_dup_dir.join(format!("{}.DupRate_plot.png", sample_name));
        rna::rseqc::plots::read_duplication_plot(&result, sample_name, &plot_path)?;
        let p = rseqc_read_dup_dir.display().to_string();
        ui.output_item("read_duplication", &format!("{p}/{sample_name}.*"));
        written.push(("read_duplication".into(), p));
    }

    // --- infer_experiment ---
    if let Some(accum) = accums.infer_exp {
        std::fs::create_dir_all(&rseqc_infer_exp_dir)?;
        let result = accum.into_result();
        let output_path = rseqc_infer_exp_dir.join(format!("{}.infer_experiment.txt", sample_name));
        rna::rseqc::infer_experiment::write_infer_experiment(&result, &output_path)?;
        let p = output_path.display().to_string();
        ui.output_item("infer_experiment", &p);
        ui.output_detail(&format!(
            "{} usable reads sampled",
            format_count(result.total_sampled),
        ));
        written.push(("infer_experiment".into(), p));
        infer_experiment_result = Some(result);
    }

    // --- read_distribution ---
    if let Some(accum) = accums.read_dist {
        std::fs::create_dir_all(&rseqc_read_dist_dir)?;
        let rd_regions = params
            .rd_regions
            .context("rd_regions must be Some when read_distribution accumulator exists")?;
        let result = accum.into_result(rd_regions);
        let output_path =
            rseqc_read_dist_dir.join(format!("{}.read_distribution.txt", sample_name));
        rna::rseqc::read_distribution::write_read_distribution(&result, &output_path)?;
        let p = output_path.display().to_string();
        ui.output_item("read_distribution", &p);
        ui.output_detail(&format!(
            "{} reads, {} tags, {} assigned",
            format_count(result.total_reads),
            format_count(result.total_tags),
            format_count(result.total_tags - result.unassigned_tags),
        ));
        written.push(("read_distribution".into(), p));
    }

    // --- junction_annotation ---
    if let Some(accum) = accums.junc_annot {
        std::fs::create_dir_all(&rseqc_junc_annot_dir)?;
        let prefix = rseqc_junc_annot_dir
            .join(sample_name)
            .to_string_lossy()
            .to_string();
        let results = accum.into_result(bam_header_refs);

        let xls_path = rseqc_junc_annot_dir.join(format!("{}.junction.xls", sample_name));
        rna::rseqc::junction_annotation::write_junction_xls(&results, &xls_path)?;

        let bed_out_path = rseqc_junc_annot_dir.join(format!("{}.junction.bed", sample_name));
        rna::rseqc::junction_annotation::write_junction_bed(&results, &bed_out_path)?;

        let interact_path =
            rseqc_junc_annot_dir.join(format!("{}.junction.Interact.bed", sample_name));
        rna::rseqc::junction_annotation::write_junction_interact_bed(
            &results,
            bam_path,
            &interact_path,
        )?;

        let r_path = rseqc_junc_annot_dir.join(format!("{}.junction_plot.r", sample_name));
        rna::rseqc::junction_annotation::write_junction_plot_r(&results, &prefix, &r_path)?;

        rna::rseqc::plots::junction_annotation_plot(&results, &prefix, sample_name)?;

        let summary_path =
            rseqc_junc_annot_dir.join(format!("{}.junction_annotation.log", sample_name));
        rna::rseqc::junction_annotation::write_summary(&results, &summary_path, params.gtf_path)?;

        // Only print the detailed junction summary in verbose mode
        if ui.is_verbose() {
            rna::rseqc::junction_annotation::print_summary(&results);
        }

        let p = rseqc_junc_annot_dir.display().to_string();
        ui.output_item("junction_annotation", &format!("{p}/{sample_name}.*"));
        written.push(("junction_annotation".into(), p));
    }

    // --- junction_saturation ---
    if let Some(accum) = accums.junc_sat {
        std::fs::create_dir_all(&rseqc_junc_sat_dir)?;
        let prefix = rseqc_junc_sat_dir
            .join(sample_name)
            .to_string_lossy()
            .to_string();
        let known_junctions = params
            .known_junctions
            .context("known_junctions must be Some when junction_saturation accumulator exists")?;
        let results = accum.into_result(
            known_junctions,
            params.junction_saturation_percentile_floor as u32,
            params.junction_saturation_percentile_ceiling as u32,
            params.junction_saturation_percentile_step as u32,
            params.junction_saturation_min_coverage as u32,
            params.junction_saturation_seed,
        );

        rna::rseqc::junction_saturation::write_r_script(&results, &prefix)?;

        let plot_path =
            rseqc_junc_sat_dir.join(format!("{}.junctionSaturation_plot.png", sample_name));
        rna::rseqc::plots::junction_saturation_plot(&results, sample_name, &plot_path)?;

        let summary_path = format!("{prefix}.junctionSaturation_summary.txt");
        rna::rseqc::junction_saturation::write_summary(&results, &summary_path)?;

        let p = rseqc_junc_sat_dir.display().to_string();
        ui.output_item("junction_saturation", &format!("{p}/{sample_name}.*"));
        written.push(("junction_saturation".into(), p));
    }

    // --- inner_distance ---
    if let Some(accum) = accums.inner_dist {
        std::fs::create_dir_all(&rseqc_inner_dist_dir)?;
        let prefix = rseqc_inner_dist_dir
            .join(sample_name)
            .to_string_lossy()
            .to_string();
        let results = accum.into_result(
            params.inner_distance_lower_bound,
            params.inner_distance_upper_bound,
            params.inner_distance_step,
        )?;

        let detail_path = format!("{prefix}.inner_distance.txt");
        rna::rseqc::inner_distance::write_detail_file(&results, &detail_path)?;

        let freq_path = format!("{prefix}.inner_distance_freq.txt");
        rna::rseqc::inner_distance::write_freq_file(&results, &freq_path)?;

        let r_path = format!("{prefix}.inner_distance_plot.r");
        rna::rseqc::inner_distance::write_r_script(
            &results,
            &prefix,
            &r_path,
            params.inner_distance_step,
        )?;

        let plot_path =
            rseqc_inner_dist_dir.join(format!("{}.inner_distance_plot.png", sample_name));
        rna::rseqc::plots::inner_distance_plot(
            &results,
            params.inner_distance_step,
            params.inner_distance_lower_bound,
            params.inner_distance_upper_bound,
            sample_name,
            &plot_path,
        )?;

        let summary_path = format!("{prefix}.inner_distance_summary.txt");
        rna::rseqc::inner_distance::write_summary(&results, &summary_path)?;

        let mean_path = format!("{prefix}.inner_distance_mean.txt");
        rna::rseqc::inner_distance::write_mean_file(&results, sample_name, &mean_path)?;

        let p = rseqc_inner_dist_dir.display().to_string();
        ui.output_item("inner_distance", &format!("{p}/{sample_name}.*"));
        ui.output_detail(&format!(
            "{} read pairs processed",
            format_count(results.total_pairs),
        ));
        written.push(("inner_distance".into(), p));
    }

    // --- TIN ---
    if let Some(accum) = accums.tin {
        let rseqc_tin_dir = if flat {
            outdir.to_path_buf()
        } else {
            outdir.join("rseqc").join("tin")
        };
        std::fs::create_dir_all(&rseqc_tin_dir)?;
        let prefix = rseqc_tin_dir
            .join(sample_name)
            .to_string_lossy()
            .to_string();
        let tin_index = params
            .tin_index
            .as_ref()
            .context("TIN index must exist when TIN accumulator is present")?;
        let results = accum.into_result(tin_index);
        rna::rseqc::tin::write_tin(&results, Path::new(&format!("{prefix}.tin.xls")))?;
        rna::rseqc::tin::write_tin_summary(
            &results,
            bam_path,
            Path::new(&format!("{prefix}.summary.txt")),
        )?;
        let p = format!("{prefix}.tin.xls");
        ui.output_item("TIN", &p);
        ui.output_detail(&format!(
            "{} transcripts",
            format_count(results.len() as u64)
        ));
        written.push(("TIN".into(), p));
    }

    // --- preseq (library complexity) ---
    if let Some(mut accum) = accums.preseq {
        let preseq_dir = if flat {
            outdir.to_path_buf()
        } else {
            outdir.join("preseq")
        };
        std::fs::create_dir_all(&preseq_dir)?;
        let output_path = preseq_dir.join(format!("{}.lc_extrap.txt", sample_name));
        accum.finalize();
        let total_reads = accum.total_fragments;
        let n_distinct = accum.n_distinct();
        let histogram = accum.into_histogram();
        debug!(
            "preseq: {} histogram bins, {} total reads, {} distinct",
            histogram.len(),
            total_reads,
            n_distinct,
        );
        match rna::preseq::estimate_complexity(
            &histogram,
            total_reads,
            n_distinct,
            &params.config.preseq,
        ) {
            Ok(result) => {
                rna::preseq::write_output(
                    &result,
                    &output_path,
                    params.config.preseq.confidence_level,
                )?;
                let p = output_path.display().to_string();
                ui.output_item("preseq", &p);
                ui.output_detail(&format!("{} extrapolation points", result.curve.len(),));
                written.push(("preseq".into(), p));
            }
            Err(e) => {
                ui.warn(&format!("preseq: skipped — {}", e));
            }
        }
    }

    Ok(RseqcOutputs {
        written,
        infer_experiment: infer_experiment_result,
    })
}
