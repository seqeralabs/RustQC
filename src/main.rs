//! RustQC - Fast quality control tools for sequencing data
//!
//! A collection of Rust-based QC tools for bioinformatics.
//! The `rna` subcommand runs all RNA-Seq QC analyses in a single pass:
//! dupRadar duplication rate analysis, featureCounts-compatible output,
//! and RSeQC-equivalent metrics (bam_stat, infer_experiment, read_duplication,
//! read_distribution, junction_annotation, junction_saturation, inner_distance).
//! Individual tools can be disabled via the YAML config file.

mod cli;
mod config;
mod gtf;
mod io;
mod rna;

use anyhow::{ensure, Context, Result};
use indexmap::IndexMap;
use log::{info, warn};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::time::Instant;

use rust_htslib::bam::Read as BamRead;

use rna::rseqc::accumulators::{RseqcAccumulators, RseqcAnnotations, RseqcConfig};

fn main() -> Result<()> {
    // Initialize logging
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
        .format_timestamp(None)
        .init();

    let cli = cli::parse_args();

    match cli.command {
        cli::Commands::Rna(args) => run_rna(args),
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
    if let Some(ref gtf) = args.gtf {
        parts.push(format!("--gtf {}", shell_escape(gtf)));
    }
    if let Some(ref bed) = args.bed {
        parts.push(format!("--bed {}", shell_escape(bed)));
    }
    parts.push(format!("-s {}", args.stranded));
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

/// Simple shell-escaping: wrap in quotes if the string contains spaces.
fn shell_escape(s: &str) -> String {
    if s.contains(' ') {
        format!("\"{}\"", s)
    } else {
        s.to_string()
    }
}

/// Run the full RNA-Seq QC pipeline: dupRadar + featureCounts + RSeQC analyses.
///
/// When multiple BAM files are provided, they are processed in parallel using
/// rayon. The GTF annotation is parsed once and shared across all BAM files.
/// Available threads are distributed across the parallel BAM processing jobs.
fn run_rna(args: cli::RnaArgs) -> Result<()> {
    // Validate thread count
    ensure!(
        args.threads >= 1,
        "--threads must be at least 1 (got {})",
        args.threads
    );

    // Load configuration file if provided
    let mut config = if let Some(ref config_path) = args.config {
        let cfg = config::Config::from_file(Path::new(config_path))?;
        info!("Loaded config from: {}", config_path);
        if cfg.has_chromosome_mapping() {
            info!(
                "Chromosome name mapping: {} entries",
                cfg.chromosome_mapping.len()
            );
        }
        cfg
    } else {
        config::Config::default()
    };

    // Apply CLI overrides to preseq config
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

    let start = Instant::now();
    let n_bams = args.input.len();
    info!("RustQC rna v{}", env!("CARGO_PKG_VERSION"));
    info!(
        "Input file{}: {}",
        if n_bams > 1 { "s" } else { "" },
        args.input.join(", ")
    );
    if let Some(ref gtf) = args.gtf {
        info!("GTF file: {}", gtf);
    }
    if let Some(ref bed) = args.bed {
        info!("BED file: {}", bed);
        warn!(
            "BED-only mode: dupRadar and featureCounts analyses will be skipped. \
             Use --gtf instead to run all analyses."
        );
    }
    if let Some(ref reference) = args.reference {
        info!("Reference FASTA: {}", reference);
    }
    info!(
        "Stranded: {}",
        match args.stranded {
            0 => "unstranded",
            1 => "forward",
            2 => "reverse",
            _ => "unknown",
        }
    );
    info!("Paired: {}", args.paired);
    info!("Threads: {}", args.threads);

    // Determine biotype attribute name (CLI overrides config, with auto-detection fallback)
    let configured_biotype = args
        .biotype_attribute
        .clone()
        .unwrap_or_else(|| config.featurecounts.biotype_attribute.clone());

    // Determine which extra GTF attributes we need, and parse GTF if provided
    let mut extra_attributes: Vec<String> = Vec::new();
    let mut biotype_attribute = configured_biotype.clone();
    let need_biotype = config.any_biotype_output();

    // GTF mode: parse GTF and detect biotype attributes
    let genes: Option<IndexMap<String, gtf::Gene>> = if let Some(ref gtf_path) = args.gtf {
        if need_biotype {
            // Check if the configured biotype attribute exists in the GTF.
            // If not found and the user didn't explicitly set it, try common alternatives:
            // Ensembl GTFs use "gene_biotype", GENCODE GTFs use "gene_type".
            let user_explicit = args.biotype_attribute.is_some();
            if gtf::attribute_exists_in_gtf(gtf_path, &biotype_attribute, 1000) {
                extra_attributes.push(biotype_attribute.clone());
                info!("Biotype attribute: {}", biotype_attribute);
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
                for alt in alternatives {
                    if gtf::attribute_exists_in_gtf(gtf_path, alt, 1000) {
                        info!(
                            "Biotype attribute '{}' not found in GTF, using '{}' instead",
                            biotype_attribute, alt
                        );
                        biotype_attribute = alt.to_string();
                        extra_attributes.push(biotype_attribute.clone());
                        found = true;
                        break;
                    }
                }
                if !found {
                    warn!(
                        "Biotype attribute '{}' not found in GTF, skipping biotype outputs",
                        configured_biotype
                    );
                }
            } else {
                warn!(
                    "Biotype attribute '{}' not found in GTF, skipping biotype outputs",
                    biotype_attribute
                );
            }
        }

        // Step 1: Parse GTF annotation (shared across all BAM files)
        info!("Parsing GTF annotation...");
        let gtf_start = Instant::now();
        let parsed = gtf::parse_gtf(gtf_path, &extra_attributes)?;
        info!(
            "Parsed {} genes in {:.2}s",
            parsed.len(),
            gtf_start.elapsed().as_secs_f64()
        );
        Some(parsed)
    } else {
        None
    };

    // Build RSeQC data structures from annotation (GTF or BED).
    // These are built once and shared across all BAM files.
    // Each tool's data is only built when enabled in the config.
    macro_rules! build_rseqc_data {
        ($enabled:expr, $gtf_msg:expr, $gtf_builder:expr, $bed_msg:expr, $bed_builder:expr) => {
            if $enabled {
                if let Some(ref genes) = genes {
                    info!($gtf_msg);
                    Some($gtf_builder(genes))
                } else if let Some(ref bed_path) = args.bed {
                    info!($bed_msg);
                    Some($bed_builder(bed_path)?)
                } else {
                    None
                }
            } else {
                None
            }
        };
    }

    let gene_model = build_rseqc_data!(
        config.infer_experiment.enabled,
        "Building gene model from GTF for infer_experiment...",
        rna::rseqc::infer_experiment::GeneModel::from_genes,
        "Loading gene model from BED file for infer_experiment...",
        rna::rseqc::infer_experiment::GeneModel::from_bed
    );

    let ref_junctions = build_rseqc_data!(
        config.junction_annotation.enabled,
        "Building reference junctions from GTF...",
        rna::rseqc::common::build_reference_junctions_from_genes,
        "Parsing reference junctions from BED...",
        rna::rseqc::common::parse_reference_junctions_from_bed
    );

    let known_junctions = build_rseqc_data!(
        config.junction_saturation.enabled,
        "Building known junction set from GTF...",
        rna::rseqc::common::build_known_junctions_from_genes,
        "Parsing known junction set from BED...",
        rna::rseqc::common::parse_known_junctions_from_bed
    );

    let rd_regions = build_rseqc_data!(
        config.read_distribution.enabled,
        "Building genomic region sets from GTF...",
        rna::rseqc::read_distribution::build_regions_from_genes,
        "Building genomic region sets from BED...",
        rna::rseqc::read_distribution::build_regions_from_bed
    );

    let exon_bitset = build_rseqc_data!(
        config.inner_distance.enabled,
        "Building exon bitset from GTF...",
        rna::rseqc::inner_distance::ExonBitset::from_genes,
        "Building exon bitset from BED...",
        rna::rseqc::inner_distance::ExonBitset::from_bed
    );

    let transcript_tree = build_rseqc_data!(
        config.inner_distance.enabled,
        "Building transcript tree from GTF...",
        rna::rseqc::inner_distance::TranscriptTree::from_genes,
        "Building transcript tree from BED...",
        rna::rseqc::inner_distance::TranscriptTree::from_bed
    );

    let tin_sample_size = config.tin.sample_size.unwrap_or(100) as usize;
    let tin_index = build_rseqc_data!(
        config.tin.enabled,
        "Building TIN index from GTF...",
        |genes| rna::rseqc::tin::TinIndex::from_genes(genes, tin_sample_size),
        "Building TIN index from BED...",
        |bed_path| rna::rseqc::tin::TinIndex::from_bed(bed_path, tin_sample_size)
    );

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
        info!(
            "Processing {} BAM files ({} in parallel, {} threads each)",
            n_bams, n_parallel, threads_per_bam
        );
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
        stranded: args.stranded,
        paired: args.paired,
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
        infer_experiment_sample_size: args.infer_experiment_sample_size,
        min_intron: args.min_intron,
        junction_saturation_min_coverage: args.junction_saturation_min_coverage,
        junction_saturation_percentile_floor: args.junction_saturation_percentile_floor,
        junction_saturation_percentile_ceiling: args.junction_saturation_percentile_ceiling,
        junction_saturation_percentile_step: args.junction_saturation_percentile_step,
        inner_distance_sample_size: args.inner_distance_sample_size,
        inner_distance_lower_bound: args.inner_distance_lower_bound,
        inner_distance_upper_bound: args.inner_distance_upper_bound,
        inner_distance_step: args.inner_distance_step,
        tin_index: tin_index.as_ref(),
        tin_sample_size: config.tin.sample_size.unwrap_or(100) as usize,
        tin_min_coverage: config.tin.min_coverage.unwrap_or(10),
        gtf_path: args.gtf.as_deref(),
    };

    // Step 2: Process all alignment files (in parallel when multiple)
    let bam_results: Vec<Result<()>> = if n_bams == 1 {
        // Single file: use all threads directly, no outer rayon pool needed
        vec![process_single_bam(
            &args.input[0],
            genes.as_ref(),
            args.threads,
            &shared,
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
                    process_single_bam(bam_path, genes.as_ref(), threads_per_bam, &shared)
                })
                .collect()
        })
    };

    // Report results
    let mut n_ok = 0;
    let mut n_err = 0;
    for (bam_path, result) in args.input.iter().zip(bam_results.into_iter()) {
        match result {
            Ok(()) => n_ok += 1,
            Err(e) => {
                n_err += 1;
                log::error!("Failed to process {}: {:#}", bam_path, e);
            }
        }
    }

    if n_bams > 1 {
        info!(
            "Processed {} file{}: {} succeeded, {} failed",
            n_bams, "s", n_ok, n_err,
        );
    }

    info!("Total runtime: {:.2}s", start.elapsed().as_secs_f64());

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
    /// Library strandedness (0/1/2).
    stranded: u8,
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
    config: &'a config::Config,
    /// GTF attribute name for biotype counting.
    biotype_attribute: &'a str,
    /// Whether the biotype attribute was found in the GTF.
    biotype_in_gtf: bool,
    /// Reconstructed command line for featureCounts header.
    command_line: &'a str,
    /// Pre-built gene model for infer_experiment (from GTF or BED).
    gene_model: Option<&'a rna::rseqc::infer_experiment::GeneModel>,
    /// Pre-built reference junctions for junction_annotation (from GTF or BED).
    ref_junctions: Option<&'a rna::rseqc::common::ReferenceJunctions>,
    /// Pre-built known junction set for junction_saturation (from GTF or BED).
    known_junctions: Option<&'a rna::rseqc::common::KnownJunctionSet>,
    /// Pre-built genomic region sets for read_distribution (from GTF or BED).
    rd_regions: Option<&'a rna::rseqc::read_distribution::RegionSets>,
    /// Pre-built exon bitset for inner_distance (from GTF or BED).
    exon_bitset: Option<&'a rna::rseqc::inner_distance::ExonBitset>,
    /// Pre-built transcript tree for inner_distance (from GTF or BED).
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
    /// Pre-built TIN index for transcript integrity analysis (from GTF or BED).
    tin_index: Option<&'a rna::rseqc::tin::TinIndex>,
    /// Number of equally-spaced positions to sample per transcript for TIN.
    tin_sample_size: usize,
    /// Minimum read-start count per transcript to compute TIN.
    tin_min_coverage: u32,
    /// Path to GTF file (for Qualimap report output).
    gtf_path: Option<&'a str>,
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
/// When `genes` is `None` (BED-only mode), the dupRadar and featureCounts
/// analyses are skipped entirely and only RSeQC tools are run.
///
/// # Arguments
///
/// * `bam_path` - Path to the duplicate-marked alignment file (SAM/BAM/CRAM)
/// * `genes` - Parsed GTF gene annotations (None in BED-only mode)
/// * `threads` - Number of threads for this file's read counting
/// * `params` - Shared parameters (config, annotations, tool settings)
fn process_single_bam(
    bam_path: &str,
    genes: Option<&IndexMap<String, gtf::Gene>>,
    threads: usize,
    params: &SharedParams,
) -> Result<()> {
    let bam_stem = Path::new(bam_path)
        .file_stem()
        .context("Input path has no filename")?
        .to_str()
        .context("Input filename is not valid UTF-8")?;

    let bam_start = Instant::now();
    let config = params.config;
    let outdir = params.outdir;

    // === Build RSeQC config and annotations ===
    let ref_chroms: HashSet<String> = params
        .known_junctions
        .map(|kj| {
            kj.junctions
                .iter()
                .map(|s| s.split(':').next().unwrap_or("").to_uppercase())
                .collect()
        })
        .unwrap_or_default();

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
    };

    let rseqc_annotations = RseqcAnnotations {
        gene_model: params.gene_model,
        ref_junctions: params.ref_junctions,
        known_junctions: params.known_junctions,
        rd_regions: params.rd_regions,
        exon_bitset: params.exon_bitset,
        transcript_tree: params.transcript_tree,
        ref_chroms: Some(&ref_chroms),
        tin_index: params.tin_index,
    };

    let any_rseqc_enabled = rseqc_config.bam_stat_enabled
        || rseqc_config.infer_experiment_enabled
        || rseqc_config.read_duplication_enabled
        || rseqc_config.read_distribution_enabled
        || rseqc_config.junction_annotation_enabled
        || rseqc_config.junction_saturation_enabled
        || rseqc_config.inner_distance_enabled
        || rseqc_config.preseq_enabled;

    // === Build gene body coverage position map (if enabled) ===
    let genebody_position_map = if params.config.genebody_coverage.enabled {
        genes.map(rna::genebody::TranscriptPositionMap::from_genes)
    } else {
        None
    };

    // === Build Qualimap exon index (if enabled, GTF-only) ===
    let qualimap_index = if params.config.qualimap.enabled {
        genes.map(rna::qualimap::QualimapIndex::from_genes)
    } else {
        None
    };

    // === dupRadar counting (requires GTF) ===
    let mut count_result = if let Some(genes) = genes {
        info!(
            "[{}] Counting reads (dupRadar + featureCounts + RSeQC)...",
            bam_stem
        );
        let count_start = Instant::now();
        let result = rna::dupradar::counting::count_reads(
            bam_path,
            genes,
            params.stranded,
            params.paired,
            threads,
            params.chrom_mapping,
            params.chrom_prefix,
            params.reference,
            params.skip_dup_check,
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
            genebody_position_map.as_ref(),
            qualimap_index.as_ref(),
        )?;
        info!(
            "[{}] Counting complete in {:.2}s",
            bam_stem,
            count_start.elapsed().as_secs_f64()
        );
        Some(result)
    } else {
        None
    };

    // Extract RSeQC accumulators from count_result (available in GTF mode).
    // In BED-only mode, we'll run a separate single-pass below.
    let rseqc_accums = count_result.as_mut().and_then(|cr| cr.rseqc.take());

    // === dupRadar + featureCounts outputs (GTF-only mode) ===
    // When genes are available (GTF mode), we have count_result and can produce
    // all counting-based outputs. In BED-only mode both are None and we skip.
    if let Some(count_result) = count_result {
        // genes is guaranteed to be Some when count_result is Some
        let genes = genes.expect("genes must be Some when count_result is Some");

        // Log counting summary stats
        info!(
            "[{}] Total reads processed: {} ({} fragments)",
            bam_stem, count_result.stat_total_reads, count_result.stat_total_fragments
        );
        info!("[{}] Assigned: {}", bam_stem, count_result.stat_assigned);
        info!(
            "[{}] No features: {}",
            bam_stem, count_result.stat_no_features
        );
        info!("[{}] Ambiguous: {}", bam_stem, count_result.stat_ambiguous);

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
        if config.any_featurecounts_output() {
            std::fs::create_dir_all(&fc_dir).with_context(|| {
                format!(
                    "Failed to create featurecounts output directory: {}",
                    fc_dir.display()
                )
            })?;
            info!(
                "[{}] Writing featureCounts-compatible output files...",
                bam_stem
            );

            if config.featurecounts.counts_file {
                let counts_path = fc_dir.join(format!("{}.featureCounts.tsv", bam_stem));
                rna::featurecounts::output::write_counts_file(
                    &counts_path,
                    genes,
                    &count_result,
                    bam_path,
                    params.command_line,
                )?;
                info!(
                    "[{}] Counts file written to {}",
                    bam_stem,
                    counts_path.display()
                );
            }

            if config.featurecounts.summary_file {
                let summary_path = fc_dir.join(format!("{}.featureCounts.tsv.summary", bam_stem));
                rna::featurecounts::output::write_summary_file(
                    &summary_path,
                    &count_result,
                    bam_path,
                )?;
                info!(
                    "[{}] Summary file written to {}",
                    bam_stem,
                    summary_path.display()
                );
            }

            // Biotype outputs (only if attribute was found in GTF)
            if params.biotype_in_gtf && config.any_biotype_output() {
                let biotype_counts = rna::featurecounts::output::aggregate_biotype_counts(
                    genes,
                    &count_result,
                    params.biotype_attribute,
                );
                info!(
                    "[{}] Biotype counting: {} biotypes found",
                    bam_stem,
                    biotype_counts.len()
                );

                if config.featurecounts.biotype_counts {
                    let biotype_path = fc_dir.join(format!("{}.biotype_counts.tsv", bam_stem));
                    rna::featurecounts::output::write_biotype_counts(
                        &biotype_path,
                        &biotype_counts,
                    )?;
                    info!(
                        "[{}] Biotype counts written to {}",
                        bam_stem,
                        biotype_path.display()
                    );
                }

                if config.featurecounts.biotype_counts_mqc {
                    let mqc_biotype_path =
                        fc_dir.join(format!("{}.biotype_counts_mqc.tsv", bam_stem));
                    rna::featurecounts::output::write_biotype_counts_mqc(
                        &mqc_biotype_path,
                        &biotype_counts,
                    )?;
                    info!(
                        "[{}] Biotype MultiQC file written to {}",
                        bam_stem,
                        mqc_biotype_path.display()
                    );
                }

                if config.featurecounts.biotype_rrna_mqc {
                    let mqc_rrna_path =
                        fc_dir.join(format!("{}.biotype_counts_rrna_mqc.tsv", bam_stem));
                    rna::featurecounts::output::write_biotype_rrna_mqc(
                        &mqc_rrna_path,
                        &biotype_counts,
                        count_result.fc_assigned,
                        bam_stem,
                    )?;
                    info!(
                        "[{}] rRNA MultiQC file written to {}",
                        bam_stem,
                        mqc_rrna_path.display()
                    );
                }
            }
        }

        // === dupRadar outputs ===
        if config.any_dupradar_output() {
            std::fs::create_dir_all(&dr_dir).with_context(|| {
                format!(
                    "Failed to create dupradar output directory: {}",
                    dr_dir.display()
                )
            })?;
            info!("[{}] Building duplication matrix...", bam_stem);
            let dup_matrix = rna::dupradar::dupmatrix::DupMatrix::build(genes, &count_result);

            let stats = dup_matrix.get_stats();
            info!(
                "[{}] Matrix built for {} genes ({} with reads)",
                bam_stem, stats.n_regions, stats.n_regions_covered
            );

            // Write duplication matrix
            if config.dupradar.dup_matrix {
                let matrix_path = dr_dir.join(format!("{}_dupMatrix.txt", bam_stem));
                dup_matrix.write_tsv(&matrix_path)?;
                info!(
                    "[{}] Duplication matrix written to {}",
                    bam_stem,
                    matrix_path.display()
                );
            }

            // Fit logistic regression model (needed for intercept/slope, density plot, and MultiQC)
            let need_fit = config.dupradar.intercept_slope
                || config.dupradar.density_scatter_plot
                || config.dupradar.multiqc_intercept
                || config.dupradar.multiqc_curve;

            let fit_ok = if need_fit {
                info!("[{}] Fitting duplication rate model...", bam_stem);
                let rpk_values: Vec<f64> = dup_matrix.rows.iter().map(|r| r.rpk).collect();
                let dup_rate_values: Vec<f64> =
                    dup_matrix.rows.iter().map(|r| r.dup_rate).collect();

                let fit_result =
                    rna::dupradar::fitting::duprate_exp_fit(&rpk_values, &dup_rate_values);
                match &fit_result {
                    Ok(fit) => {
                        info!(
                            "[{}] Model fit: intercept={:.6}, slope={:.6}",
                            bam_stem, fit.intercept, fit.slope
                        );
                        if config.dupradar.intercept_slope {
                            let fit_path = dr_dir.join(format!("{}_intercept_slope.txt", bam_stem));
                            rna::dupradar::plots::write_intercept_slope(fit, &fit_path)?;
                            info!(
                                "[{}] Fit results written to {}",
                                bam_stem,
                                fit_path.display()
                            );
                        }
                        Some(fit.clone())
                    }
                    Err(e) => {
                        warn!("[{}] Could not fit model: {}", bam_stem, e);
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
                info!("[{}] Generating plots...", bam_stem);
                let plot_start = Instant::now();

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
                        handle.join().map_err(|_| {
                            anyhow::anyhow!("density scatter plot thread panicked")
                        })??;
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

                info!(
                    "[{}] Plots generated in {:.2}s",
                    bam_stem,
                    plot_start.elapsed().as_secs_f64()
                );
            }

            // Write MultiQC-compatible output files
            if let Some(ref fit) = fit_ok {
                if config.dupradar.multiqc_intercept {
                    let mqc_intercept_path =
                        dr_dir.join(format!("{}_dup_intercept_mqc.txt", bam_stem));
                    rna::dupradar::plots::write_mqc_intercept(fit, bam_stem, &mqc_intercept_path)?;
                }

                if config.dupradar.multiqc_curve {
                    let mqc_curve_path =
                        dr_dir.join(format!("{}_duprateExpDensCurve_mqc.txt", bam_stem));
                    rna::dupradar::plots::write_mqc_curve(fit, &dup_matrix, &mqc_curve_path)?;
                }

                if config.dupradar.multiqc_intercept || config.dupradar.multiqc_curve {
                    info!("[{}] MultiQC output files written", bam_stem);
                }
            }

            // Print dupRadar summary statistics
            info!("[{}] --- dupRadar Summary ---", bam_stem);
            info!("[{}] Total genes: {}", bam_stem, stats.n_regions);
            info!(
                "[{}] Genes with reads: {} ({:.1}%)",
                bam_stem,
                stats.n_regions_covered,
                stats.f_regions_covered * 100.0
            );
            info!(
                "[{}] Genes with duplication: {} ({:.1}%)",
                bam_stem,
                stats.n_regions_duplication,
                stats.f_regions_duplication * 100.0
            );
        }
        // === Qualimap RNA-Seq QC output ===
        if let (Some(ref qm_result), Some(ref qm_index)) = (&count_result.qualimap, &qualimap_index)
        {
            let qm_dir = if params.flat_output {
                outdir.to_path_buf()
            } else {
                outdir.join("qualimap")
            };

            // Extract junction event counts from the RSeQC accumulator (if available)
            let junction_counts = rseqc_accums
                .as_ref()
                .and_then(|a| a.junc_annot.as_ref())
                .map(|ja| {
                    (
                        ja.known_events,
                        ja.partial_novel_events,
                        ja.complete_novel_events,
                    )
                });

            rna::qualimap::output::write_qualimap_results(
                qm_result,
                qm_index,
                bam_path,
                params.gtf_path.unwrap_or(""),
                params.stranded,
                &qm_dir,
                bam_stem,
                junction_counts,
            )?;
            info!(
                "[{}] Qualimap RNA-Seq QC results written to {:?}",
                bam_stem, qm_dir
            );
        }
    } // end if let Some(count_result) — GTF-only outputs

    // === RSeQC analyses (post-processing of single-pass accumulators) ===
    let rseqc_accums = if let Some(accums) = rseqc_accums {
        accums
    } else if genes.is_none() {
        // BED-only mode: run a standalone single-pass for RSeQC
        run_rseqc_single_pass(bam_path, bam_stem, params, threads)?
    } else {
        // No RSeQC tools enabled — skip
        info!("[{}] No RSeQC tools enabled, skipping", bam_stem);
        RseqcAccumulators::empty()
    };
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
    write_rseqc_outputs(bam_path, bam_stem, params, rseqc_accums, &bam_header_refs)?;

    info!(
        "[{}] Completed in {:.2}s",
        bam_stem,
        bam_start.elapsed().as_secs_f64()
    );

    Ok(())
}

// ============================================================================
// RSeQC output writing (post-processing of single-pass accumulators)
// ============================================================================

/// Write all RSeQC outputs from the single-pass accumulators.
///
/// Converts accumulated data to tool-specific result types and writes all
/// output files, plots, and summaries.
fn write_rseqc_outputs(
    bam_path: &str,
    bam_stem: &str,
    params: &SharedParams,
    accums: RseqcAccumulators,
    bam_header_refs: &[(String, u64)],
) -> Result<()> {
    let outdir = params.outdir;

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
    // These consume data from BamStatAccum, so they are written before
    // the main bam_stat output to share the same result.
    let samtools_dir = if flat {
        outdir.to_path_buf()
    } else {
        outdir.join("samtools")
    };

    // --- bam_stat ---
    if let Some(accum) = accums.bam_stat {
        info!("[{}] Writing bam_stat results...", bam_stem);
        std::fs::create_dir_all(&rseqc_bam_stat_dir)?;
        let result = accum.into_result();
        if params.config.bam_stat.enabled {
            let output_path = rseqc_bam_stat_dir.join(format!("{}.bam_stat.txt", bam_stem));
            rna::rseqc::bam_stat::write_bam_stat(&result, &output_path)?;
        }

        // --- flagstat ---
        if params.config.flagstat.enabled {
            info!("[{}] Writing flagstat results...", bam_stem);
            std::fs::create_dir_all(&samtools_dir)?;
            let flagstat_path = samtools_dir.join(format!("{}.flagstat", bam_stem));
            rna::rseqc::flagstat::write_flagstat(&result, &flagstat_path)?;
        }

        // --- idxstats ---
        if params.config.idxstats.enabled {
            info!("[{}] Writing idxstats results...", bam_stem);
            std::fs::create_dir_all(&samtools_dir)?;
            let idxstats_path = samtools_dir.join(format!("{}.idxstats", bam_stem));
            rna::rseqc::idxstats::write_idxstats(&result, bam_header_refs, &idxstats_path)?;
        }

        // --- samtools stats SN ---
        if params.config.samtools_stats.enabled {
            info!("[{}] Writing samtools stats results...", bam_stem);
            std::fs::create_dir_all(&samtools_dir)?;
            let stats_path = samtools_dir.join(format!("{}.stats", bam_stem));
            rna::rseqc::stats::write_stats(&result, &stats_path)?;
        }
    }

    // --- read_duplication ---
    if let Some(accum) = accums.read_dup {
        info!("[{}] Writing read_duplication results...", bam_stem);
        std::fs::create_dir_all(&rseqc_read_dup_dir)?;
        let result = accum.into_result();
        rna::rseqc::read_duplication::write_read_duplication(
            &result,
            &rseqc_read_dup_dir,
            bam_stem,
        )?;
        let plot_path = rseqc_read_dup_dir.join(format!("{}.DupRate_plot.png", bam_stem));
        rna::rseqc::plots::read_duplication_plot(&result, &plot_path)?;
    }

    // --- infer_experiment ---
    if let Some(accum) = accums.infer_exp {
        info!("[{}] Writing infer_experiment results...", bam_stem);
        std::fs::create_dir_all(&rseqc_infer_exp_dir)?;
        let result = accum.into_result();
        let output_path = rseqc_infer_exp_dir.join(format!("{}.infer_experiment.txt", bam_stem));
        rna::rseqc::infer_experiment::write_infer_experiment(&result, &output_path)?;
        info!(
            "[{}] Total {} usable reads were sampled",
            bam_stem, result.total_sampled
        );
    }

    // --- read_distribution ---
    if let Some(accum) = accums.read_dist {
        info!("[{}] Writing read_distribution results...", bam_stem);
        std::fs::create_dir_all(&rseqc_read_dist_dir)?;
        // Safety: rd_regions is Some when read_distribution_enabled is true,
        // which is required for this accumulator to exist.
        let result = accum.into_result(params.rd_regions.unwrap());
        let output_path = rseqc_read_dist_dir.join(format!("{}.read_distribution.txt", bam_stem));
        rna::rseqc::read_distribution::write_read_distribution(&result, &output_path)?;
        info!(
            "[{}] Total {} reads, {} tags ({} assigned)",
            bam_stem,
            result.total_reads,
            result.total_tags,
            result.total_tags - result.unassigned_tags
        );
    }

    // --- junction_annotation ---
    if let Some(accum) = accums.junc_annot {
        info!("[{}] Writing junction_annotation results...", bam_stem);
        std::fs::create_dir_all(&rseqc_junc_annot_dir)?;
        let prefix = rseqc_junc_annot_dir
            .join(bam_stem)
            .to_string_lossy()
            .to_string();
        let results = accum.into_result();

        let xls_path = rseqc_junc_annot_dir.join(format!("{}.junction.xls", bam_stem));
        rna::rseqc::junction_annotation::write_junction_xls(&results, &xls_path)?;

        let bed_out_path = rseqc_junc_annot_dir.join(format!("{}.junction.bed", bam_stem));
        rna::rseqc::junction_annotation::write_junction_bed(&results, &bed_out_path)?;

        let r_path = rseqc_junc_annot_dir.join(format!("{}.junction_plot.r", bam_stem));
        rna::rseqc::junction_annotation::write_junction_plot_r(&results, &prefix, &r_path)?;

        rna::rseqc::plots::junction_annotation_plot(&results, &prefix)?;

        let summary_path =
            rseqc_junc_annot_dir.join(format!("{}.junction_annotation.txt", bam_stem));
        rna::rseqc::junction_annotation::write_summary(&results, &summary_path)?;

        rna::rseqc::junction_annotation::print_summary(&results);
    }

    // --- junction_saturation ---
    if let Some(accum) = accums.junc_sat {
        info!("[{}] Writing junction_saturation results...", bam_stem);
        std::fs::create_dir_all(&rseqc_junc_sat_dir)?;
        let prefix = rseqc_junc_sat_dir
            .join(bam_stem)
            .to_string_lossy()
            .to_string();
        // Safety: known_junctions is Some when junction_saturation_enabled is true,
        // which is required for this accumulator to exist.
        let results = accum.into_result(
            params.known_junctions.unwrap(),
            params.junction_saturation_percentile_floor as u32,
            params.junction_saturation_percentile_ceiling as u32,
            params.junction_saturation_percentile_step as u32,
            params.junction_saturation_min_coverage as u32,
        );

        rna::rseqc::junction_saturation::write_r_script(&results, &prefix)?;

        let plot_path =
            rseqc_junc_sat_dir.join(format!("{}.junctionSaturation_plot.png", bam_stem));
        rna::rseqc::plots::junction_saturation_plot(&results, &plot_path)?;

        let summary_path = format!("{prefix}.junctionSaturation_summary.txt");
        rna::rseqc::junction_saturation::write_summary(&results, &summary_path)?;
    }

    // --- inner_distance ---
    if let Some(accum) = accums.inner_dist {
        info!("[{}] Writing inner_distance results...", bam_stem);
        std::fs::create_dir_all(&rseqc_inner_dist_dir)?;
        let prefix = rseqc_inner_dist_dir
            .join(bam_stem)
            .to_string_lossy()
            .to_string();
        let results = accum.into_result(
            params.inner_distance_lower_bound,
            params.inner_distance_upper_bound,
            params.inner_distance_step,
        );

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
        rna::rseqc::plots::inner_distance_plot(&results, params.inner_distance_step, &plot_path)?;

        let summary_path = format!("{prefix}.inner_distance_summary.txt");
        rna::rseqc::inner_distance::write_summary(&results, &summary_path)?;

        info!(
            "[{}] inner_distance: {} read pairs processed",
            bam_stem, results.total_pairs
        );
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
            .expect("TIN index must exist when TIN accumulator is present");
        let results = accum.into_result(tin_index);
        info!(
            "[{}] Writing TIN results ({} transcripts)...",
            bam_stem,
            results.len()
        );
        rna::rseqc::tin::write_tin(&results, Path::new(&format!("{prefix}.tin.xls")))?;
        rna::rseqc::tin::write_tin_summary(
            &results,
            bam_path,
            Path::new(&format!("{prefix}.summary.txt")),
        )?;
    }

    // --- preseq (library complexity) ---
    if let Some(accum) = accums.preseq {
        let preseq_dir = if flat {
            outdir.to_path_buf()
        } else {
            outdir.join("preseq")
        };
        std::fs::create_dir_all(&preseq_dir)?;
        let output_path = preseq_dir.join(format!("{}.lc_extrap.txt", bam_stem));
        let total_reads = accum.total_fragments;
        let n_distinct = accum.n_distinct();
        let histogram = accum.into_histogram();
        info!(
            "[{}] Running preseq complexity estimation ({} histogram bins, {} total reads, {} distinct)...",
            bam_stem,
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
                info!(
                    "[{}] preseq: {} extrapolation points written to {}",
                    bam_stem,
                    result.curve.len(),
                    output_path.display()
                );
            }
            Err(e) => {
                log::warn!("[{}] preseq: skipped — {}", bam_stem, e);
            }
        }
    }

    Ok(())
}

// ============================================================================
// BED-only mode: standalone single-pass for RSeQC
// ============================================================================

/// Run RSeQC analyses via a standalone single-pass BAM read (BED-only mode).
///
/// When no GTF is provided, dupRadar/featureCounts are skipped but RSeQC tools
/// still need to run. This function does a single sequential BAM pass using the
/// accumulator infrastructure, then returns the results for output writing.
fn run_rseqc_single_pass(
    bam_path: &str,
    bam_stem: &str,
    params: &SharedParams,
    threads: usize,
) -> Result<RseqcAccumulators> {
    use rust_htslib::bam::{self, Read as BamRead};

    info!("[{}] Running RSeQC single-pass (BED mode)...", bam_stem);
    let pass_start = Instant::now();

    // Build config and annotations
    let ref_chroms: Option<std::collections::HashSet<String>> = params.known_junctions.map(|kj| {
        kj.junctions
            .iter()
            .map(|j| j.split(':').next().unwrap_or("").to_uppercase())
            .collect()
    });

    let rseqc_config = RseqcConfig {
        mapq_cut: params.mapq_cut,
        min_intron: params.min_intron,
        infer_experiment_sample_size: params.infer_experiment_sample_size,
        junction_saturation_min_coverage: params.junction_saturation_min_coverage as u32,
        junction_saturation_sample_start: params.junction_saturation_percentile_floor as u32,
        junction_saturation_sample_end: params.junction_saturation_percentile_ceiling as u32,
        junction_saturation_sample_step: params.junction_saturation_percentile_step as u32,
        inner_distance_sample_size: params.inner_distance_sample_size,
        inner_distance_lower_bound: params.inner_distance_lower_bound,
        inner_distance_upper_bound: params.inner_distance_upper_bound,
        inner_distance_step: params.inner_distance_step,
        bam_stat_enabled: params.config.bam_stat.enabled
            || params.config.flagstat.enabled
            || params.config.idxstats.enabled
            || params.config.samtools_stats.enabled,
        infer_experiment_enabled: params.config.infer_experiment.enabled
            && params.gene_model.is_some(),
        read_duplication_enabled: params.config.read_duplication.enabled,
        read_distribution_enabled: params.config.read_distribution.enabled
            && params.rd_regions.is_some(),
        junction_annotation_enabled: params.config.junction_annotation.enabled
            && params.ref_junctions.is_some(),
        junction_saturation_enabled: params.config.junction_saturation.enabled
            && params.known_junctions.is_some(),
        inner_distance_enabled: params.config.inner_distance.enabled
            && params.exon_bitset.is_some()
            && params.transcript_tree.is_some(),
        tin_enabled: params.config.tin.enabled && params.tin_index.is_some(),
        tin_sample_size: params.tin_sample_size,
        tin_min_coverage: params.tin_min_coverage,
        preseq_enabled: params.config.preseq.enabled,
    };

    let annotations = RseqcAnnotations {
        gene_model: params.gene_model,
        ref_junctions: params.ref_junctions,
        known_junctions: params.known_junctions,
        rd_regions: params.rd_regions,
        exon_bitset: params.exon_bitset,
        transcript_tree: params.transcript_tree,
        ref_chroms: ref_chroms.as_ref(),
        tin_index: params.tin_index,
    };

    let mut accums = RseqcAccumulators::new(&rseqc_config, Some(&annotations));

    // Open BAM reader
    let mut bam = bam::Reader::from_path(bam_path)
        .with_context(|| format!("Failed to open BAM file: {bam_path}"))?;
    if let Some(ref_path) = params.reference {
        bam.set_reference(ref_path)
            .with_context(|| format!("Failed to set reference: {ref_path}"))?;
    }
    // Use all available threads for decompression in single-pass mode
    let _ = bam.set_threads(threads);

    let header = bam.header().clone();
    let tid_to_name: Vec<String> = (0..header.target_count())
        .map(|tid| String::from_utf8_lossy(header.tid2name(tid)).to_string())
        .collect();

    // Pre-compute resolved chromosome names (apply prefix/mapping like counting.rs)
    let tid_to_rseqc_chrom: Vec<String> = tid_to_name
        .iter()
        .map(|name| {
            if let Some(mapped) = params.chrom_mapping.get(name.as_str()) {
                mapped.clone()
            } else if let Some(prefix) = params.chrom_prefix {
                format!("{prefix}{name}")
            } else {
                name.clone()
            }
        })
        .collect();
    let tid_to_rseqc_upper: Vec<String> = tid_to_rseqc_chrom
        .iter()
        .map(|s| s.to_uppercase())
        .collect();

    let mut record = bam::Record::new();
    while let Some(read_result) = bam.read(&mut record) {
        read_result.context("Error reading BAM record")?;
        let tid = record.tid();
        if tid < 0 {
            // Process unmapped reads through accumulators (bam_stat wants them)
            accums.process_read(&record, "", "", &annotations, &rseqc_config);
            continue;
        }
        let idx = tid as usize;
        let rseqc_chrom = &tid_to_rseqc_chrom[idx];
        let rseqc_upper = &tid_to_rseqc_upper[idx];
        accums.process_read(
            &record,
            rseqc_chrom,
            rseqc_upper,
            &annotations,
            &rseqc_config,
        );
    }

    info!(
        "[{}] RSeQC single-pass completed in {:.2}s",
        bam_stem,
        pass_start.elapsed().as_secs_f64()
    );

    Ok(accums)
}
