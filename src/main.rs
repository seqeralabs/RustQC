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
use log::info;
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

    let start = Instant::now();
    info!("RustQC rna v{}", env!("CARGO_PKG_VERSION"));
    info!("BAM file: {}", args.bam);
    info!("GTF file: {}", args.gtf);
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

    // Step 1: Parse GTF annotation
    info!("Parsing GTF annotation...");
    let gtf_start = Instant::now();
    let genes = gtf::parse_gtf(&args.gtf)?;
    info!(
        "Parsed {} genes in {:.2}s",
        genes.len(),
        gtf_start.elapsed().as_secs_f64()
    );

    // Step 2: Count reads with featureCounts-compatible logic
    info!("Counting reads across 4 modes...");
    let count_start = Instant::now();
    let chrom_mapping = config.bam_to_gtf_mapping();
    let count_result = counting::count_reads(
        &args.bam,
        &genes,
        args.stranded,
        args.paired,
        args.threads,
        &chrom_mapping,
        config.chromosome_prefix(),
    )?;
    info!(
        "Counting complete in {:.2}s",
        count_start.elapsed().as_secs_f64()
    );

    // Step 3: Build duplication matrix
    info!("Building duplication matrix...");
    let dup_matrix = dupmatrix::DupMatrix::build(&genes, &count_result);

    let stats = dup_matrix.get_stats();
    info!(
        "Matrix built for {} genes ({} with reads)",
        stats.n_regions, stats.n_regions_covered
    );

    // Step 4: Write duplication matrix
    let outdir = Path::new(&args.outdir);
    std::fs::create_dir_all(outdir)?;
    let bam_stem = Path::new(&args.bam)
        .file_stem()
        .context("BAM path has no filename")?
        .to_str()
        .context("BAM filename is not valid UTF-8")?;

    let matrix_path = outdir.join(format!("{}_dupMatrix.txt", bam_stem));
    dup_matrix.write_tsv(&matrix_path)?;
    info!("Duplication matrix written to {}", matrix_path.display());

    // Step 5: Fit logistic regression model
    info!("Fitting duplication rate model...");
    let rpk_values: Vec<f64> = dup_matrix.rows.iter().map(|r| r.rpk).collect();
    let dup_rate_values: Vec<f64> = dup_matrix.rows.iter().map(|r| r.dup_rate).collect();

    let fit_result = fitting::duprate_exp_fit(&rpk_values, &dup_rate_values);
    let fit_ok = match &fit_result {
        Ok(fit) => {
            info!(
                "Model fit: intercept={:.6}, slope={:.6}",
                fit.intercept, fit.slope
            );
            let fit_path = outdir.join(format!("{}_intercept_slope.txt", bam_stem));
            plots::write_intercept_slope(fit, &fit_path)?;
            info!("Fit results written to {}", fit_path.display());
            Some(fit.clone())
        }
        Err(e) => {
            info!("Warning: Could not fit model: {}", e);
            None
        }
    };

    // Step 6: Generate plots (in parallel — all plots read shared immutable data)
    info!("Generating plots...");
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
        "Plots generated in {:.2}s",
        plot_start.elapsed().as_secs_f64()
    );

    // Step 7: Write MultiQC-compatible output files
    if let Some(ref fit) = fit_ok {
        let mqc_intercept_path = outdir.join(format!("{}_dup_intercept_mqc.txt", bam_stem));
        plots::write_mqc_intercept(fit, bam_stem, &mqc_intercept_path)?;

        let mqc_curve_path = outdir.join(format!("{}_duprateExpDensCurve_mqc.txt", bam_stem));
        plots::write_mqc_curve(fit, &dup_matrix, &mqc_curve_path)?;
        info!("MultiQC output files written");
    }

    // Step 8: Print summary statistics
    info!("--- Summary ---");
    info!("Total genes: {}", stats.n_regions);
    info!(
        "Genes with reads: {} ({:.1}%)",
        stats.n_regions_covered,
        stats.f_regions_covered * 100.0
    );
    info!(
        "Genes with duplication: {} ({:.1}%)",
        stats.n_regions_duplication,
        stats.f_regions_duplication * 100.0
    );

    info!("Total runtime: {:.2}s", start.elapsed().as_secs_f64());

    Ok(())
}
