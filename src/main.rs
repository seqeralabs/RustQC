//! RustQC - Fast quality control tools for sequencing data
//!
//! A collection of Rust-based QC tools for bioinformatics.
//! The `rna` subcommand is a reimplementation of the dupRadar Bioconductor
//! R package, analysing PCR duplicate rates as a function of gene expression
//! level in RNA-Seq datasets.

mod cli;
mod config;
mod counting;
mod dupmatrix;
mod fitting;
mod gtf;
mod plots;

use anyhow::{Context, Result};
use indexmap::IndexMap;
use log::info;
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
    }
}

/// Run the RNA-Seq duplication rate analysis (dupRadar equivalent).
///
/// When multiple BAM files are provided, they are processed in parallel using
/// rayon. The GTF annotation is parsed once and shared across all BAM files.
/// Available threads are distributed across the parallel BAM processing jobs.
fn run_rna(args: cli::RnaArgs) -> Result<()> {
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

    // Step 1: Parse GTF annotation (shared across all BAM files)
    info!("Parsing GTF annotation...");
    let gtf_start = Instant::now();
    let genes = gtf::parse_gtf(&args.gtf)?;
    info!(
        "Parsed {} genes in {:.2}s",
        genes.len(),
        gtf_start.elapsed().as_secs_f64()
    );

    let chrom_mapping = config.alignment_to_gtf_mapping();
    let chrom_prefix = config.chromosome_prefix().map(|s| s.to_owned());

    // Validate that BAM file stems are unique (otherwise outputs would collide).
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
    std::fs::create_dir_all(outdir)?;

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
                log::error!("Failed to process {}: {:?}", bam_path, e);
            }
        }
    }

    if n_bams > 1 {
        info!(
            "Processed {} file{}: {} succeeded, {} failed",
            n_bams,
            if n_bams > 1 { "s" } else { "" },
            n_ok,
            n_err,
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
/// This runs the complete dupRadar-equivalent analysis for one file:
/// counting, matrix construction, model fitting, plotting, and MultiQC output.
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
) -> Result<()> {
    let bam_stem = Path::new(bam_path)
        .file_stem()
        .context("Input path has no filename")?
        .to_str()
        .context("Input filename is not valid UTF-8")?;

    let bam_start = Instant::now();
    info!("[{}] Counting reads across 4 modes...", bam_stem);
    let count_start = Instant::now();
    let count_result = counting::count_reads(
        bam_path,
        genes,
        stranded,
        paired,
        threads,
        chrom_mapping,
        chrom_prefix,
        reference,
        skip_dup_check,
    )?;
    info!(
        "[{}] Counting complete in {:.2}s",
        bam_stem,
        count_start.elapsed().as_secs_f64()
    );

    // Build duplication matrix
    info!("[{}] Building duplication matrix...", bam_stem);
    let dup_matrix = dupmatrix::DupMatrix::build(genes, &count_result);

    let stats = dup_matrix.get_stats();
    info!(
        "[{}] Matrix built for {} genes ({} with reads)",
        bam_stem, stats.n_regions, stats.n_regions_covered
    );

    // Write duplication matrix
    let matrix_path = outdir.join(format!("{}_dupMatrix.txt", bam_stem));
    dup_matrix.write_tsv(&matrix_path)?;
    info!(
        "[{}] Duplication matrix written to {}",
        bam_stem,
        matrix_path.display()
    );

    // Fit logistic regression model
    info!("[{}] Fitting duplication rate model...", bam_stem);
    let rpk_values: Vec<f64> = dup_matrix.rows.iter().map(|r| r.rpk).collect();
    let dup_rate_values: Vec<f64> = dup_matrix.rows.iter().map(|r| r.dup_rate).collect();

    let fit_result = fitting::duprate_exp_fit(&rpk_values, &dup_rate_values);
    let fit_ok = match &fit_result {
        Ok(fit) => {
            info!(
                "[{}] Model fit: intercept={:.6}, slope={:.6}",
                bam_stem, fit.intercept, fit.slope
            );
            let fit_path = outdir.join(format!("{}_intercept_slope.txt", bam_stem));
            plots::write_intercept_slope(fit, &fit_path)?;
            info!(
                "[{}] Fit results written to {}",
                bam_stem,
                fit_path.display()
            );
            Some(fit.clone())
        }
        Err(e) => {
            info!("[{}] Warning: Could not fit model: {}", bam_stem, e);
            None
        }
    };

    // Generate plots (in parallel — all plots read shared immutable data)
    info!("[{}] Generating plots...", bam_stem);
    let plot_start = Instant::now();

    let rpkm_threshold_rpk = fit_ok.as_ref().and_then(|_fit| {
        let rpkm_values: Vec<f64> = dup_matrix.rows.iter().map(|r| r.rpkm).collect();
        fitting::compute_rpkm_threshold_rpk(&rpk_values, &rpkm_values, 0.5)
    });

    let density_path = outdir.join(format!("{}_duprateExpDens.png", bam_stem));
    let boxplot_path = outdir.join(format!("{}_duprateExpBoxplot.png", bam_stem));
    let histogram_path = outdir.join(format!("{}_expressionHist.png", bam_stem));

    std::thread::scope(|s| -> Result<()> {
        // Density scatter plot (only if fit succeeded)
        let density_handle = fit_ok.as_ref().map(|fit| {
            let dm_ref = &dup_matrix;
            let thresh = rpkm_threshold_rpk;
            let path = &density_path;
            s.spawn(move || plots::density_scatter_plot(dm_ref, fit, thresh, bam_stem, path))
        });

        // Boxplot
        let boxplot_handle = {
            let dm_ref = &dup_matrix;
            let path = &boxplot_path;
            s.spawn(move || plots::duprate_boxplot(dm_ref, bam_stem, path))
        };

        // Histogram
        let histogram_handle = {
            let dm_ref = &dup_matrix;
            let path = &histogram_path;
            s.spawn(move || plots::expression_histogram(dm_ref, bam_stem, path))
        };

        // Collect results
        if let Some(handle) = density_handle {
            handle
                .join()
                .expect("density scatter plot thread panicked")?;
        }
        boxplot_handle.join().expect("boxplot thread panicked")?;
        histogram_handle
            .join()
            .expect("histogram thread panicked")?;

        Ok(())
    })?;

    info!(
        "[{}] Plots generated in {:.2}s",
        bam_stem,
        plot_start.elapsed().as_secs_f64()
    );

    // Write MultiQC-compatible output files
    if let Some(ref fit) = fit_ok {
        let mqc_intercept_path = outdir.join(format!("{}_dup_intercept_mqc.txt", bam_stem));
        plots::write_mqc_intercept(fit, bam_stem, &mqc_intercept_path)?;

        let mqc_curve_path = outdir.join(format!("{}_duprateExpDensCurve_mqc.txt", bam_stem));
        plots::write_mqc_curve(fit, &dup_matrix, &mqc_curve_path)?;
        info!("[{}] MultiQC output files written", bam_stem);
    }

    // Print summary statistics
    info!("[{}] --- Summary ---", bam_stem);
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
    info!(
        "[{}] Completed in {:.2}s",
        bam_stem,
        bam_start.elapsed().as_secs_f64()
    );

    Ok(())
}
