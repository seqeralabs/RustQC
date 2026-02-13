//! RustQC - Fast quality control tools for sequencing data
//!
//! A collection of Rust-based QC tools for bioinformatics.
//! The `rna` subcommand is a reimplementation of the dupRadar Bioconductor
//! R package, analysing PCR duplicate rates as a function of gene expression
//! level in RNA-Seq datasets. It also produces featureCounts-compatible output
//! files and biotype count summaries in a single pass.

mod cli;
mod config;
mod gtf;
mod rna;

use anyhow::{ensure, Context, Result};
use indexmap::IndexMap;
use log::{info, warn};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::time::Instant;

fn main() -> Result<()> {
    // Initialize logging
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
        .format_timestamp(None)
        .init();

    let cli = cli::parse_args();

    match cli.command {
        cli::Commands::Rna(args) => run_rna(args),
        cli::Commands::BamStat(args) => run_bam_stat(args),
        cli::Commands::InferExperiment(args) => run_infer_experiment(args),
        cli::Commands::ReadDuplication(args) => run_read_duplication(args),
        cli::Commands::ReadDistribution(args) => run_read_distribution(args),
        cli::Commands::JunctionAnnotation(args) => run_junction_annotation(args),
        cli::Commands::JunctionSaturation(args) => run_junction_saturation(args),
        cli::Commands::InnerDistance(args) => run_inner_distance(args),
    }
}

/// Reconstruct the command line for the featureCounts-compatible header comment.
fn reconstruct_command_line(args: &cli::RnaArgs) -> String {
    let mut parts = vec![format!(
        "rustqc rna {} --gtf {}",
        args.input
            .iter()
            .map(|s| shell_escape(s))
            .collect::<Vec<_>>()
            .join(" "),
        shell_escape(&args.gtf)
    )];
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

/// Run the RNA-Seq duplication rate analysis (dupRadar equivalent)
/// and featureCounts-compatible output generation.
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
    let config = if let Some(ref config_path) = args.config {
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
    info!("GTF file: {}", args.gtf);
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

    // Determine which extra GTF attributes we need
    let mut extra_attributes: Vec<String> = Vec::new();
    let mut biotype_attribute = configured_biotype.clone();
    let need_biotype = config.any_biotype_output();
    if need_biotype {
        // Check if the configured biotype attribute exists in the GTF.
        // If not found and the user didn't explicitly set it, try common alternatives:
        // Ensembl GTFs use "gene_biotype", GENCODE GTFs use "gene_type".
        let user_explicit = args.biotype_attribute.is_some();
        if gtf::attribute_exists_in_gtf(&args.gtf, &biotype_attribute, 1000) {
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
                if gtf::attribute_exists_in_gtf(&args.gtf, alt, 1000) {
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
    let genes = gtf::parse_gtf(&args.gtf, &extra_attributes)?;
    info!(
        "Parsed {} genes in {:.2}s",
        genes.len(),
        gtf_start.elapsed().as_secs_f64()
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

    // Step 2: Process all alignment files (in parallel when multiple)
    let bam_results: Vec<Result<()>> = if n_bams == 1 {
        // Single file: use all threads directly, no outer rayon pool needed
        vec![process_single_bam(
            &args.input[0],
            &genes,
            args.stranded,
            args.paired,
            args.threads,
            &chrom_mapping,
            chrom_prefix.as_deref(),
            outdir,
            args.reference.as_deref(),
            args.skip_dup_check,
            &config,
            &biotype_attribute,
            biotype_in_gtf,
            &command_line,
        )]
    } else {
        // Multiple files: process in parallel with a dedicated rayon pool
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(n_parallel)
            .build()
            .context("Failed to create rayon thread pool for parallel BAM processing")?;

        let reference = args.reference.as_deref();
        let skip_dup_check = args.skip_dup_check;

        pool.install(|| {
            args.input
                .par_iter()
                .map(|bam_path| {
                    process_single_bam(
                        bam_path,
                        &genes,
                        args.stranded,
                        args.paired,
                        threads_per_bam,
                        &chrom_mapping,
                        chrom_prefix.as_deref(),
                        outdir,
                        reference,
                        skip_dup_check,
                        &config,
                        &biotype_attribute,
                        biotype_in_gtf,
                        &command_line,
                    )
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

/// Process a single alignment file through the full analysis pipeline.
///
/// This runs the complete analysis for one file: counting, featureCounts output,
/// biotype counting, duplication matrix, model fitting, plotting, and MultiQC output.
///
/// # Arguments
///
/// * `bam_path` - Path to the duplicate-marked alignment file (SAM/BAM/CRAM)
/// * `genes` - Parsed GTF gene annotations (shared, read-only)
/// * `stranded` - Library strandedness (0/1/2)
/// * `paired` - Whether the library is paired-end
/// * `threads` - Number of threads for this file's read counting
/// * `chrom_mapping` - Alignment-to-GTF chromosome name mapping
/// * `chrom_prefix` - Optional chromosome name prefix
/// * `outdir` - Output directory for results
/// * `reference` - Optional reference FASTA for CRAM files
/// * `skip_dup_check` - Whether to skip duplicate-marking validation
/// * `config` - Configuration for conditional outputs
/// * `biotype_attribute` - GTF attribute name for biotype counting
/// * `biotype_in_gtf` - Whether the biotype attribute was found in the GTF
/// * `command_line` - Reconstructed command line for featureCounts header
#[allow(clippy::too_many_arguments)]
fn process_single_bam(
    bam_path: &str,
    genes: &IndexMap<String, gtf::Gene>,
    stranded: u8,
    paired: bool,
    threads: usize,
    chrom_mapping: &HashMap<String, String>,
    chrom_prefix: Option<&str>,
    outdir: &Path,
    reference: Option<&str>,
    skip_dup_check: bool,
    config: &config::Config,
    biotype_attribute: &str,
    biotype_in_gtf: bool,
    command_line: &str,
) -> Result<()> {
    let bam_stem = Path::new(bam_path)
        .file_stem()
        .context("Input path has no filename")?
        .to_str()
        .context("Input filename is not valid UTF-8")?;

    let bam_start = Instant::now();
    info!("[{}] Counting reads across 4 modes...", bam_stem);
    let count_start = Instant::now();
    let count_result = rna::dupradar::counting::count_reads(
        bam_path,
        genes,
        stranded,
        paired,
        threads,
        chrom_mapping,
        chrom_prefix,
        reference,
        skip_dup_check,
        biotype_attribute,
    )?;
    info!(
        "[{}] Counting complete in {:.2}s",
        bam_stem,
        count_start.elapsed().as_secs_f64()
    );

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

    // === featureCounts outputs ===
    if config.any_featurecounts_output() {
        info!(
            "[{}] Writing featureCounts-compatible output files...",
            bam_stem
        );

        if config.featurecounts.counts_file {
            let counts_path = outdir.join(format!("{}.featureCounts.tsv", bam_stem));
            rna::featurecounts::output::write_counts_file(
                &counts_path,
                genes,
                &count_result,
                bam_path,
                command_line,
            )?;
            info!(
                "[{}] Counts file written to {}",
                bam_stem,
                counts_path.display()
            );
        }

        if config.featurecounts.summary_file {
            let summary_path = outdir.join(format!("{}.featureCounts.tsv.summary", bam_stem));
            rna::featurecounts::output::write_summary_file(&summary_path, &count_result, bam_path)?;
            info!(
                "[{}] Summary file written to {}",
                bam_stem,
                summary_path.display()
            );
        }

        // Biotype outputs (only if attribute was found in GTF)
        if biotype_in_gtf && config.any_biotype_output() {
            let biotype_counts = rna::featurecounts::output::aggregate_biotype_counts(
                genes,
                &count_result,
                biotype_attribute,
            );
            info!(
                "[{}] Biotype counting: {} biotypes found",
                bam_stem,
                biotype_counts.len()
            );

            if config.featurecounts.biotype_counts {
                let biotype_path = outdir.join(format!("{}.biotype_counts.tsv", bam_stem));
                rna::featurecounts::output::write_biotype_counts(&biotype_path, &biotype_counts)?;
                info!(
                    "[{}] Biotype counts written to {}",
                    bam_stem,
                    biotype_path.display()
                );
            }

            if config.featurecounts.biotype_counts_mqc {
                let mqc_biotype_path = outdir.join(format!("{}.biotype_counts_mqc.tsv", bam_stem));
                rna::featurecounts::output::write_biotype_counts_mqc(
                    &mqc_biotype_path,
                    &biotype_counts,
                    bam_stem,
                )?;
                info!(
                    "[{}] Biotype MultiQC file written to {}",
                    bam_stem,
                    mqc_biotype_path.display()
                );
            }

            if config.featurecounts.biotype_rrna_mqc {
                let mqc_rrna_path =
                    outdir.join(format!("{}.biotype_counts_rrna_mqc.tsv", bam_stem));
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
        info!("[{}] Building duplication matrix...", bam_stem);
        let dup_matrix = rna::dupradar::dupmatrix::DupMatrix::build(genes, &count_result);

        let stats = dup_matrix.get_stats();
        info!(
            "[{}] Matrix built for {} genes ({} with reads)",
            bam_stem, stats.n_regions, stats.n_regions_covered
        );

        // Write duplication matrix
        if config.dupradar.dup_matrix {
            let matrix_path = outdir.join(format!("{}_dupMatrix.txt", bam_stem));
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
            let dup_rate_values: Vec<f64> = dup_matrix.rows.iter().map(|r| r.dup_rate).collect();

            let fit_result = rna::dupradar::fitting::duprate_exp_fit(&rpk_values, &dup_rate_values);
            match &fit_result {
                Ok(fit) => {
                    info!(
                        "[{}] Model fit: intercept={:.6}, slope={:.6}",
                        bam_stem, fit.intercept, fit.slope
                    );
                    if config.dupradar.intercept_slope {
                        let fit_path = outdir.join(format!("{}_intercept_slope.txt", bam_stem));
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

            let density_path = outdir.join(format!("{}_duprateExpDens.png", bam_stem));
            let boxplot_path = outdir.join(format!("{}_duprateExpBoxplot.png", bam_stem));
            let histogram_path = outdir.join(format!("{}_expressionHist.png", bam_stem));

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

            info!(
                "[{}] Plots generated in {:.2}s",
                bam_stem,
                plot_start.elapsed().as_secs_f64()
            );
        }

        // Write MultiQC-compatible output files
        if let Some(ref fit) = fit_ok {
            if config.dupradar.multiqc_intercept {
                let mqc_intercept_path = outdir.join(format!("{}_dup_intercept_mqc.txt", bam_stem));
                rna::dupradar::plots::write_mqc_intercept(fit, bam_stem, &mqc_intercept_path)?;
            }

            if config.dupradar.multiqc_curve {
                let mqc_curve_path =
                    outdir.join(format!("{}_duprateExpDensCurve_mqc.txt", bam_stem));
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

    info!(
        "[{}] Completed in {:.2}s",
        bam_stem,
        bam_start.elapsed().as_secs_f64()
    );

    Ok(())
}

// ============================================================================
// bam-stat subcommand
// ============================================================================

/// Run the bam-stat analysis for one or more BAM files.
fn run_bam_stat(args: cli::BamStatArgs) -> Result<()> {
    let start = Instant::now();
    let outdir = Path::new(&args.outdir);
    std::fs::create_dir_all(outdir)
        .with_context(|| format!("Failed to create output directory: {}", args.outdir))?;

    for bam_path in &args.input {
        let bam_start = Instant::now();
        let bam_stem = Path::new(bam_path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("output");

        info!("[{}] Running bam_stat analysis...", bam_stem);

        let result =
            rna::rseqc::bam_stat::bam_stat(bam_path, args.mapq_cut, args.reference.as_deref())?;

        let output_path = outdir.join(format!("{}.bam_stat.txt", bam_stem));
        rna::rseqc::bam_stat::write_bam_stat(&result, args.mapq_cut, &output_path)?;

        info!(
            "[{}] bam_stat completed in {:.2}s",
            bam_stem,
            bam_start.elapsed().as_secs_f64()
        );
    }

    if args.input.len() > 1 {
        info!(
            "All BAM files processed in {:.2}s",
            start.elapsed().as_secs_f64()
        );
    }

    Ok(())
}

// ============================================================================

/// Run the infer-experiment analysis for one or more BAM files.
fn run_infer_experiment(args: cli::InferExperimentArgs) -> Result<()> {
    let start = Instant::now();
    let outdir = Path::new(&args.outdir);
    std::fs::create_dir_all(outdir)
        .with_context(|| format!("Failed to create output directory: {}", args.outdir))?;

    info!("Loading gene model from BED file: {}", args.bed);
    let gene_model = rna::rseqc::infer_experiment::GeneModel::from_bed(&args.bed)?;

    for bam_path in &args.input {
        let bam_start = Instant::now();
        let bam_stem = Path::new(bam_path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("output");

        info!("[{}] Running infer_experiment analysis...", bam_stem);

        let result = rna::rseqc::infer_experiment::infer_experiment(
            bam_path,
            &gene_model,
            args.mapq_cut,
            args.sample_size,
            args.reference.as_deref(),
        )?;

        let output_path = outdir.join(format!("{}.infer_experiment.txt", bam_stem));
        rna::rseqc::infer_experiment::write_infer_experiment(&result, &output_path)?;

        info!(
            "[{}] Total {} usable reads were sampled",
            bam_stem, result.total_sampled
        );
        info!(
            "[{}] infer_experiment completed in {:.2}s",
            bam_stem,
            bam_start.elapsed().as_secs_f64()
        );
    }

    if args.input.len() > 1 {
        info!(
            "All BAM files processed in {:.2}s",
            start.elapsed().as_secs_f64()
        );
    }

    Ok(())
}

/// Run read duplication analysis for one or more BAM files.
fn run_read_duplication(args: cli::ReadDuplicationArgs) -> Result<()> {
    let start = Instant::now();
    let outdir = Path::new(&args.outdir);
    std::fs::create_dir_all(outdir).context("Failed to create output directory")?;

    for bam_path in &args.input {
        let bam_start = Instant::now();
        let bam_stem = Path::new(bam_path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("output")
            .to_string();

        info!("[{}] Running read duplication analysis...", bam_stem);

        let result = rna::rseqc::read_duplication::read_duplication(
            bam_path,
            args.mapq_cut,
            args.reference.as_deref(),
        )?;

        rna::rseqc::read_duplication::write_read_duplication(&result, outdir, &bam_stem)?;

        info!(
            "[{}] read_duplication completed in {:.2}s",
            bam_stem,
            bam_start.elapsed().as_secs_f64()
        );
    }

    if args.input.len() > 1 {
        info!(
            "All BAM files processed in {:.2}s",
            start.elapsed().as_secs_f64()
        );
    }

    Ok(())
}

// ============================================================================
// read-distribution subcommand
// ============================================================================

/// Run read distribution analysis for one or more BAM files.
fn run_read_distribution(args: cli::ReadDistributionArgs) -> Result<()> {
    let start = Instant::now();
    let outdir = Path::new(&args.outdir);
    std::fs::create_dir_all(outdir)
        .with_context(|| format!("Failed to create output directory: {}", args.outdir))?;

    for bam_path in &args.input {
        let bam_start = Instant::now();
        let bam_stem = Path::new(bam_path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("output");

        info!("[{}] Running read_distribution analysis...", bam_stem);

        let result = rna::rseqc::read_distribution::read_distribution(
            bam_path,
            &args.bed,
            args.reference.as_deref(),
        )?;

        let output_path = outdir.join(format!("{}.read_distribution.txt", bam_stem));
        rna::rseqc::read_distribution::write_read_distribution(&result, &output_path)?;

        info!(
            "[{}] Total {} reads, {} tags ({} assigned)",
            bam_stem,
            result.total_reads,
            result.total_tags,
            result.total_tags - result.unassigned_tags
        );
        info!(
            "[{}] read_distribution completed in {:.2}s",
            bam_stem,
            bam_start.elapsed().as_secs_f64()
        );
    }

    if args.input.len() > 1 {
        info!(
            "All BAM files processed in {:.2}s",
            start.elapsed().as_secs_f64()
        );
    }

    Ok(())
}

// ============================================================================
// junction-annotation subcommand
// ============================================================================

/// Run junction annotation analysis for one or more BAM files.
fn run_junction_annotation(args: cli::JunctionAnnotationArgs) -> Result<()> {
    let start = Instant::now();
    let outdir = Path::new(&args.outdir);
    std::fs::create_dir_all(outdir)
        .with_context(|| format!("Failed to create output directory: {}", args.outdir))?;

    for bam_path in &args.input {
        let bam_start = Instant::now();
        let bam_stem = Path::new(bam_path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("output");

        info!("[{}] Running junction_annotation analysis...", bam_stem);

        let results = rna::rseqc::junction_annotation::junction_annotation(
            bam_path,
            &args.bed,
            args.min_intron,
            args.mapq_cut,
            args.reference.as_deref(),
        )?;

        // Write output files
        let xls_path = outdir.join(format!("{}.junction.xls", bam_stem));
        rna::rseqc::junction_annotation::write_junction_xls(&results, &xls_path)?;

        let bed_path = outdir.join(format!("{}.junction.bed", bam_stem));
        rna::rseqc::junction_annotation::write_junction_bed(&results, &bed_path)?;

        let r_path = outdir.join(format!("{}.junction_plot.r", bam_stem));
        let r_prefix = outdir.join(bam_stem).to_string_lossy().to_string();
        rna::rseqc::junction_annotation::write_junction_plot_r(&results, &r_prefix, &r_path)?;

        let summary_path = outdir.join(format!("{}.junction_annotation.txt", bam_stem));
        rna::rseqc::junction_annotation::write_summary(&results, &summary_path)?;

        // Print summary to stderr
        rna::rseqc::junction_annotation::print_summary(&results);

        info!(
            "[{}] junction_annotation completed in {:.2}s",
            bam_stem,
            bam_start.elapsed().as_secs_f64()
        );
    }

    if args.input.len() > 1 {
        info!(
            "All BAM files processed in {:.2}s",
            start.elapsed().as_secs_f64()
        );
    }

    Ok(())
}

// ============================================================================
// junction-saturation subcommand
// ============================================================================

/// Run the junction_saturation analysis on one or more BAM files.
fn run_junction_saturation(args: cli::JunctionSaturationArgs) -> Result<()> {
    let start = Instant::now();
    let outdir = Path::new(&args.outdir);
    std::fs::create_dir_all(outdir)
        .with_context(|| format!("Failed to create output directory: {}", args.outdir))?;

    for bam_path in &args.input {
        let bam_start = Instant::now();
        let bam_stem = Path::new(bam_path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("output");

        info!("[{}] Running junction_saturation analysis...", bam_stem);

        let results = rna::rseqc::junction_saturation::junction_saturation(
            bam_path,
            &args.bed,
            args.min_intron,
            args.mapq_cut,
            args.min_coverage as u32,
            args.percentile_floor as u32,
            args.percentile_ceiling as u32,
            args.percentile_step as u32,
            args.reference.as_deref(),
        )?;

        // Write R script and summary
        let prefix = outdir.join(bam_stem).to_string_lossy().to_string();
        rna::rseqc::junction_saturation::write_r_script(&results, &prefix)?;
        let summary_path = format!("{prefix}.junctionSaturation_summary.txt");
        rna::rseqc::junction_saturation::write_summary(&results, &summary_path)?;

        info!(
            "[{}] junction_saturation completed in {:.2}s",
            bam_stem,
            bam_start.elapsed().as_secs_f64()
        );
    }

    if args.input.len() > 1 {
        info!(
            "All BAM files processed in {:.2}s",
            start.elapsed().as_secs_f64()
        );
    }

    Ok(())
}

// ===================================================================
// Inner Distance
// ===================================================================

/// Run the inner_distance subcommand.
fn run_inner_distance(args: cli::InnerDistanceArgs) -> Result<()> {
    let start = Instant::now();
    let outdir = &args.outdir;
    std::fs::create_dir_all(outdir)
        .with_context(|| format!("Failed to create output directory: {}", outdir))?;

    for bam_path in &args.input {
        let bam_start = Instant::now();
        let bam_stem = Path::new(bam_path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("output");

        let results = rna::rseqc::inner_distance::inner_distance(
            bam_path,
            &args.bed,
            args.sample_size,
            args.mapq_cut,
            args.lower_bound,
            args.upper_bound,
            args.step,
            args.reference.as_deref(),
        )?;

        // Write output files
        let prefix = Path::new(outdir)
            .join(bam_stem)
            .to_string_lossy()
            .to_string();
        let detail_path = format!("{prefix}.inner_distance.txt");
        rna::rseqc::inner_distance::write_detail_file(&results, &detail_path)?;

        let freq_path = format!("{prefix}.inner_distance_freq.txt");
        rna::rseqc::inner_distance::write_freq_file(&results, &freq_path)?;

        let r_path = format!("{prefix}.inner_distance_plot.r");
        rna::rseqc::inner_distance::write_r_script(&results, &prefix, &r_path, args.step)?;

        let summary_path = format!("{prefix}.inner_distance_summary.txt");
        rna::rseqc::inner_distance::write_summary(&results, &summary_path)?;

        info!(
            "[{}] inner_distance completed in {:.2}s — {} read pairs processed",
            bam_stem,
            bam_start.elapsed().as_secs_f64(),
            results.total_pairs
        );
    }

    if args.input.len() > 1 {
        info!(
            "All BAM files processed in {:.2}s",
            start.elapsed().as_secs_f64()
        );
    }

    Ok(())
}
