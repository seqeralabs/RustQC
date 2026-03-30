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

mod cli;
mod config;
mod gtf;
mod io;
mod rna;
mod summary;
mod ui;

use anyhow::{ensure, Context, Result};
use indexmap::IndexMap;
use log::debug;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use ui::{format_count, format_duration, format_pct, Ui, Verbosity};

use rust_htslib::bam::Read as BamRead;

use rna::rseqc::accumulators::{RseqcAccumulators, RseqcAnnotations, RseqcConfig};

/// Common BAM filename suffixes added by alignment and duplicate-marking tools.
///
/// These are stripped (case-insensitively) from the BAM file stem to derive a
/// clean sample name for output content (e.g., intercept/slope tables). Longer
/// / more specific suffixes come first so they are matched greedily.
const BAM_SUFFIXES_TO_STRIP: &[&str] = &[
    ".Aligned.sortedByCoord.out.markDups",
    ".Aligned.sortedByCoord.out",
    ".Aligned.toTranscriptome.out",
    ".markdup.sorted",
    ".sorted.markdup",
    ".markDups.sorted",
    ".sorted.markDups",
    ".mrkdup.sorted",
    ".sorted.mrkdup",
    ".dedup.sorted",
    ".sorted.dedup",
    ".markdup",
    ".markDups",
    ".mrkdup",
    ".dedup",
    ".sorted",
];

/// Derive a clean sample name from a BAM file stem.
///
/// Strips common alignment / duplicate-marking suffixes (case-insensitive) so
/// that e.g. `GM12878_REP1.markdup.sorted` becomes `GM12878_REP1`. If no known
/// suffix matches, the stem is returned unchanged.
fn clean_sample_name(bam_stem: &str) -> String {
    for suffix in BAM_SUFFIXES_TO_STRIP {
        if bam_stem.len() >= suffix.len() {
            let split = bam_stem.len() - suffix.len();
            if bam_stem[split..].eq_ignore_ascii_case(suffix) {
                return bam_stem[..split].to_string();
            }
        }
    }
    bam_stem.to_string()
}

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
    let cli = cli::parse_args();

    // Determine verbosity from CLI flags
    let verbosity = match &cli.command {
        cli::Commands::Rna(args) if args.quiet => Verbosity::Quiet,
        cli::Commands::Rna(args) if args.verbose => Verbosity::Verbose,
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
    }
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

    // Load configuration file if provided, then extract the `rna` section
    let mut config = if let Some(ref config_path) = args.config {
        let cfg = config::Config::from_file(Path::new(config_path))?;
        ui.detail(&format!("Loaded config from: {}", config_path));
        if cfg.rna.has_chromosome_mapping() {
            ui.detail(&format!(
                "Chromosome name mapping: {} entries",
                cfg.rna.chromosome_mapping.len()
            ));
        }
        cfg.rna
    } else {
        config::RnaConfig::default()
    };

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
    if let Some(val) = args.preseq_seed {
        config.preseq.seed = val;
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

    ui.header(
        env!("CARGO_PKG_VERSION"),
        env!("GIT_SHORT_HASH"),
        env!("BUILD_TIMESTAMP"),
    );

    // Resolve effective stranded/paired early for display (config loaded above)
    let effective_stranded = args
        .stranded
        .or(config.stranded)
        .unwrap_or(cli::Strandedness::Unstranded);
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
        inner_distance_sample_size: config.inner_distance.sample_size.unwrap_or(1_000_000),
        inner_distance_lower_bound: config.inner_distance.lower_bound.unwrap_or(-250),
        inner_distance_upper_bound: config.inner_distance.upper_bound.unwrap_or(250),
        inner_distance_step: config.inner_distance.step.unwrap_or(5),
        tin_index: tin_index.as_ref(),
        tin_sample_size,
        tin_min_coverage: config.tin.min_coverage.unwrap_or(10),
        gtf_path: &args.gtf,
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
    stranded: cli::Strandedness,
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
    let sample_name = clean_sample_name(bam_stem);

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
        &format!("{} — Counting Summary", bam_stem),
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
            let counts_path = fc_dir.join(format!("{}.featureCounts.tsv", bam_stem));
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
            let summary_path = fc_dir.join(format!("{}.featureCounts.tsv.summary", bam_stem));
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

            if config.featurecounts.biotype_counts {
                let biotype_path = fc_dir.join(format!("{}.biotype_counts.tsv", bam_stem));
                rna::featurecounts::output::write_biotype_counts(&biotype_path, &biotype_counts)?;
                let p = biotype_path.display().to_string();
                ui.output_detail(&format!("Biotype counts: {p}"));
                written_outputs.push(("biotype counts".into(), p));
            }

            if config.featurecounts.biotype_counts_mqc {
                let mqc_biotype_path = fc_dir.join(format!("{}.biotype_counts_mqc.tsv", bam_stem));
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
                    fc_dir.join(format!("{}.biotype_counts_rrna_mqc.tsv", bam_stem));
                rna::featurecounts::output::write_biotype_rrna_mqc(
                    &mqc_rrna_path,
                    &biotype_counts,
                    count_result.fc_biotype_assigned,
                    bam_stem,
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
            let matrix_path = dr_dir.join(format!("{}_dupMatrix.txt", bam_stem));
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
                        let fit_path = dr_dir.join(format!("{}_intercept_slope.txt", bam_stem));
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

            let density_path = dr_dir.join(format!("{}_duprateExpDens.png", bam_stem));
            let boxplot_path = dr_dir.join(format!("{}_duprateExpBoxplot.png", bam_stem));
            let histogram_path = dr_dir.join(format!("{}_expressionHist.png", bam_stem));

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
                                bam_stem,
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
                        rna::dupradar::plots::duprate_boxplot(dm_ref, bam_stem, path)
                    }))
                } else {
                    None
                };

                // Histogram
                let histogram_handle = if config.dupradar.expression_histogram {
                    let dm_ref = &dup_matrix;
                    let path = &histogram_path;
                    Some(s.spawn(move || {
                        rna::dupradar::plots::expression_histogram(dm_ref, bam_stem, path)
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
                let mqc_intercept_path = dr_dir.join(format!("{}_dup_intercept_mqc.txt", bam_stem));
                rna::dupradar::plots::write_mqc_intercept(fit, &sample_name, &mqc_intercept_path)?;
                written_outputs.push((
                    "dupRadar MultiQC".into(),
                    mqc_intercept_path.display().to_string(),
                ));
            }

            if config.dupradar.multiqc_curve {
                let mqc_curve_path =
                    dr_dir.join(format!("{}_duprateExpDensCurve_mqc.txt", bam_stem));
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
            bam_stem,
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
    let rseqc_outputs =
        write_rseqc_outputs(bam_path, bam_stem, params, rseqc_accums, &bam_header_refs)?;
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
    bam_stem: &str,
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
            let flagstat_path = samtools_dir.join(format!("{}.flagstat", bam_stem));
            rna::rseqc::flagstat::write_flagstat(result, &flagstat_path)?;
            let p = flagstat_path.display().to_string();
            ui.output_item("flagstat", &p);
            written.push(("flagstat".into(), p));
        }

        // --- idxstats ---
        if params.config.idxstats.enabled {
            std::fs::create_dir_all(&samtools_dir)?;
            let idxstats_path = samtools_dir.join(format!("{}.idxstats", bam_stem));
            rna::rseqc::idxstats::write_idxstats(result, bam_header_refs, &idxstats_path)?;
            let p = idxstats_path.display().to_string();
            ui.output_item("idxstats", &p);
            written.push(("idxstats".into(), p));
        }

        // --- samtools stats SN ---
        if params.config.samtools_stats.enabled {
            std::fs::create_dir_all(&samtools_dir)?;
            let stats_path = samtools_dir.join(format!("{}.stats", bam_stem));
            rna::rseqc::stats::write_stats(result, &stats_path)?;
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
            let output_path = rseqc_bam_stat_dir.join(format!("{}.bam_stat.txt", bam_stem));
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
            bam_stem,
        )?;
        let plot_path = rseqc_read_dup_dir.join(format!("{}.DupRate_plot.png", bam_stem));
        let sample_name = clean_sample_name(bam_stem);
        rna::rseqc::plots::read_duplication_plot(&result, &sample_name, &plot_path)?;
        let p = rseqc_read_dup_dir.display().to_string();
        ui.output_item("read_duplication", &format!("{p}/{bam_stem}.*"));
        written.push(("read_duplication".into(), p));
    }

    // --- infer_experiment ---
    if let Some(accum) = accums.infer_exp {
        std::fs::create_dir_all(&rseqc_infer_exp_dir)?;
        let result = accum.into_result();
        let output_path = rseqc_infer_exp_dir.join(format!("{}.infer_experiment.txt", bam_stem));
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
        let output_path = rseqc_read_dist_dir.join(format!("{}.read_distribution.txt", bam_stem));
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
            .join(bam_stem)
            .to_string_lossy()
            .to_string();
        let results = accum.into_result(bam_header_refs);

        let xls_path = rseqc_junc_annot_dir.join(format!("{}.junction.xls", bam_stem));
        rna::rseqc::junction_annotation::write_junction_xls(&results, &xls_path)?;

        let bed_out_path = rseqc_junc_annot_dir.join(format!("{}.junction.bed", bam_stem));
        rna::rseqc::junction_annotation::write_junction_bed(&results, &bed_out_path)?;

        let interact_path =
            rseqc_junc_annot_dir.join(format!("{}.junction.Interact.bed", bam_stem));
        rna::rseqc::junction_annotation::write_junction_interact_bed(
            &results,
            bam_path,
            &interact_path,
        )?;

        let r_path = rseqc_junc_annot_dir.join(format!("{}.junction_plot.r", bam_stem));
        rna::rseqc::junction_annotation::write_junction_plot_r(&results, &prefix, &r_path)?;

        let sample_name = clean_sample_name(bam_stem);
        rna::rseqc::plots::junction_annotation_plot(&results, &prefix, &sample_name)?;

        let summary_path =
            rseqc_junc_annot_dir.join(format!("{}.junction_annotation.txt", bam_stem));
        rna::rseqc::junction_annotation::write_summary(&results, &summary_path)?;

        // Only print the detailed junction summary in verbose mode
        if ui.is_verbose() {
            rna::rseqc::junction_annotation::print_summary(&results);
        }

        let p = rseqc_junc_annot_dir.display().to_string();
        ui.output_item("junction_annotation", &format!("{p}/{bam_stem}.*"));
        written.push(("junction_annotation".into(), p));
    }

    // --- junction_saturation ---
    if let Some(accum) = accums.junc_sat {
        std::fs::create_dir_all(&rseqc_junc_sat_dir)?;
        let prefix = rseqc_junc_sat_dir
            .join(bam_stem)
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
        );

        rna::rseqc::junction_saturation::write_r_script(&results, &prefix)?;

        let plot_path =
            rseqc_junc_sat_dir.join(format!("{}.junctionSaturation_plot.png", bam_stem));
        let sample_name = clean_sample_name(bam_stem);
        rna::rseqc::plots::junction_saturation_plot(&results, &sample_name, &plot_path)?;

        let summary_path = format!("{prefix}.junctionSaturation_summary.txt");
        rna::rseqc::junction_saturation::write_summary(&results, &summary_path)?;

        let p = rseqc_junc_sat_dir.display().to_string();
        ui.output_item("junction_saturation", &format!("{p}/{bam_stem}.*"));
        written.push(("junction_saturation".into(), p));
    }

    // --- inner_distance ---
    if let Some(accum) = accums.inner_dist {
        std::fs::create_dir_all(&rseqc_inner_dist_dir)?;
        let prefix = rseqc_inner_dist_dir
            .join(bam_stem)
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

        let plot_path = rseqc_inner_dist_dir.join(format!("{}.inner_distance_plot.png", bam_stem));
        let sample_name = clean_sample_name(bam_stem);
        rna::rseqc::plots::inner_distance_plot(
            &results,
            params.inner_distance_step,
            params.inner_distance_lower_bound,
            params.inner_distance_upper_bound,
            &sample_name,
            &plot_path,
        )?;

        let summary_path = format!("{prefix}.inner_distance_summary.txt");
        rna::rseqc::inner_distance::write_summary(&results, &summary_path)?;

        let mean_path = format!("{prefix}.inner_distance_mean.txt");
        let clean_name = clean_sample_name(bam_stem);
        rna::rseqc::inner_distance::write_mean_file(&results, &clean_name, &mean_path)?;

        let p = rseqc_inner_dist_dir.display().to_string();
        ui.output_item("inner_distance", &format!("{p}/{bam_stem}.*"));
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
        let prefix = rseqc_tin_dir.join(bam_stem).to_string_lossy().to_string();
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
        let output_path = preseq_dir.join(format!("{}.lc_extrap.txt", bam_stem));
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_clean_sample_name_markdup_sorted() {
        assert_eq!(
            clean_sample_name("GM12878_REP1.markdup.sorted"),
            "GM12878_REP1"
        );
    }

    #[test]
    fn test_clean_sample_name_sorted_markdup() {
        assert_eq!(clean_sample_name("SAMPLE.sorted.markdup"), "SAMPLE");
    }

    #[test]
    fn test_clean_sample_name_star_suffix() {
        assert_eq!(
            clean_sample_name("MySample.Aligned.sortedByCoord.out.markDups"),
            "MySample"
        );
    }

    #[test]
    fn test_clean_sample_name_star_no_markdups() {
        assert_eq!(
            clean_sample_name("MySample.Aligned.sortedByCoord.out"),
            "MySample"
        );
    }

    #[test]
    fn test_clean_sample_name_case_insensitive() {
        assert_eq!(clean_sample_name("SAMPLE.MARKDUP.SORTED"), "SAMPLE");
    }

    #[test]
    fn test_clean_sample_name_no_suffix() {
        assert_eq!(clean_sample_name("CleanName"), "CleanName");
    }

    #[test]
    fn test_clean_sample_name_only_sorted() {
        assert_eq!(clean_sample_name("data.sorted"), "data");
    }

    #[test]
    fn test_clean_sample_name_dedup() {
        assert_eq!(clean_sample_name("run1.dedup.sorted"), "run1");
    }

    #[test]
    fn test_clean_sample_name_mrkdup() {
        assert_eq!(clean_sample_name("abc.mrkdup"), "abc");
    }

    #[test]
    fn test_clean_sample_name_empty_after_strip() {
        assert_eq!(clean_sample_name(".sorted"), "");
    }
}
