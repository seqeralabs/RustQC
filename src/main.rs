//! RustQC - Fast quality control tools for sequencing data
//!
//! A collection of Rust-based QC tools for bioinformatics.
//! The `rna` subcommand is a reimplementation of the dupRadar Bioconductor
//! R package, analysing PCR duplicate rates as a function of gene expression
//! level in RNA-Seq datasets. It also produces featureCounts-compatible output
//! files and biotype count summaries in a single pass.

mod cli;
mod config;
mod counting;
mod dupmatrix;
mod featurecounts;
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

/// Reconstruct the command line for the featureCounts-compatible header comment.
fn reconstruct_command_line(args: &cli::RnaArgs) -> String {
    let mut parts = vec![format!(
        "rustqc rna {} {}",
        shell_escape(&args.bam),
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
                info!(
                    "Warning: biotype attribute '{}' not found in GTF, skipping biotype outputs",
                    configured_biotype
                );
            }
        } else {
            info!(
                "Warning: biotype attribute '{}' not found in GTF, skipping biotype outputs",
                biotype_attribute
            );
        }
    }

    // Step 1: Parse GTF annotation
    info!("Parsing GTF annotation...");
    let gtf_start = Instant::now();
    let genes = gtf::parse_gtf(&args.gtf, &extra_attributes)?;
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

    // Log counting summary stats
    info!(
        "Total reads processed: {} ({} fragments)",
        count_result.stat_total_reads, count_result.stat_total_fragments
    );
    info!("Assigned: {}", count_result.stat_assigned);
    info!("No features: {}", count_result.stat_no_features);
    info!("Ambiguous: {}", count_result.stat_ambiguous);

    // Step 3: Build duplication matrix (only if any dupradar outputs requested)
    let outdir = Path::new(&args.outdir);
    std::fs::create_dir_all(outdir)?;
    let bam_stem = Path::new(&args.bam)
        .file_stem()
        .context("BAM path has no filename")?
        .to_str()
        .context("BAM filename is not valid UTF-8")?;

    // === featureCounts outputs ===
    if config.any_featurecounts_output() {
        info!("Writing featureCounts-compatible output files...");
        let command_line = reconstruct_command_line(&args);

        if config.featurecounts.counts_file {
            let counts_path = outdir.join(format!("{}.featureCounts.tsv", bam_stem));
            featurecounts::write_counts_file(
                &counts_path,
                &genes,
                &count_result,
                &args.bam,
                &command_line,
            )?;
            info!("Counts file written to {}", counts_path.display());
        }

        if config.featurecounts.summary_file {
            let summary_path = outdir.join(format!("{}.featureCounts.tsv.summary", bam_stem));
            featurecounts::write_summary_file(&summary_path, &count_result, &args.bam)?;
            info!("Summary file written to {}", summary_path.display());
        }

        // Biotype outputs (only if attribute was found in GTF)
        let biotype_in_gtf = extra_attributes.contains(&biotype_attribute);
        if biotype_in_gtf && config.any_biotype_output() {
            let biotype_counts =
                featurecounts::aggregate_biotype_counts(&genes, &count_result, &biotype_attribute);
            info!("Biotype counting: {} biotypes found", biotype_counts.len());

            if config.featurecounts.biotype_counts {
                let biotype_path = outdir.join(format!("{}.biotype_counts.tsv", bam_stem));
                featurecounts::write_biotype_counts(&biotype_path, &biotype_counts)?;
                info!("Biotype counts written to {}", biotype_path.display());
            }

            if config.featurecounts.biotype_counts_mqc {
                let mqc_biotype_path = outdir.join(format!("{}.biotype_counts_mqc.tsv", bam_stem));
                featurecounts::write_biotype_counts_mqc(
                    &mqc_biotype_path,
                    &biotype_counts,
                    bam_stem,
                )?;
                info!(
                    "Biotype MultiQC file written to {}",
                    mqc_biotype_path.display()
                );
            }

            if config.featurecounts.biotype_rrna_mqc {
                let mqc_rrna_path =
                    outdir.join(format!("{}.biotype_counts_rrna_mqc.tsv", bam_stem));
                featurecounts::write_biotype_rrna_mqc(
                    &mqc_rrna_path,
                    &biotype_counts,
                    count_result.stat_assigned,
                    bam_stem,
                )?;
                info!("rRNA MultiQC file written to {}", mqc_rrna_path.display());
            }
        }
    }

    // === dupRadar outputs ===
    if config.any_dupradar_output() {
        // Build duplication matrix
        info!("Building duplication matrix...");
        let dup_matrix = dupmatrix::DupMatrix::build(&genes, &count_result);

        let stats = dup_matrix.get_stats();
        info!(
            "Matrix built for {} genes ({} with reads)",
            stats.n_regions, stats.n_regions_covered
        );

        // Write duplication matrix
        if config.dupradar.dup_matrix {
            let matrix_path = outdir.join(format!("{}_dupMatrix.txt", bam_stem));
            dup_matrix.write_tsv(&matrix_path)?;
            info!("Duplication matrix written to {}", matrix_path.display());
        }

        // Fit logistic regression model (needed for intercept/slope, density plot, and MultiQC)
        let need_fit = config.dupradar.intercept_slope
            || config.dupradar.density_scatter_plot
            || config.dupradar.multiqc_intercept
            || config.dupradar.multiqc_curve;

        let fit_ok = if need_fit {
            let rpk_values: Vec<f64> = dup_matrix.rows.iter().map(|r| r.rpk).collect();
            let dup_rate_values: Vec<f64> = dup_matrix.rows.iter().map(|r| r.dup_rate).collect();

            let fit_result = fitting::duprate_exp_fit(&rpk_values, &dup_rate_values);
            match &fit_result {
                Ok(fit) => {
                    info!(
                        "Model fit: intercept={:.6}, slope={:.6}",
                        fit.intercept, fit.slope
                    );
                    if config.dupradar.intercept_slope {
                        let fit_path = outdir.join(format!("{}_intercept_slope.txt", bam_stem));
                        plots::write_intercept_slope(fit, &fit_path)?;
                        info!("Fit results written to {}", fit_path.display());
                    }
                    Some(fit.clone())
                }
                Err(e) => {
                    info!("Warning: Could not fit model: {}", e);
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
            info!("Generating plots...");
            let plot_start = Instant::now();

            let rpkm_threshold_rpk = fit_ok.as_ref().and_then(|_fit| {
                let rpk_values: Vec<f64> = dup_matrix.rows.iter().map(|r| r.rpk).collect();
                let rpkm_values: Vec<f64> = dup_matrix.rows.iter().map(|r| r.rpkm).collect();
                fitting::compute_rpkm_threshold_rpk(&rpk_values, &rpkm_values, 0.5)
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
                        s.spawn(move || plots::density_scatter_plot(dm_ref, fit, thresh, path))
                    })
                } else {
                    None
                };

                // Boxplot
                let boxplot_handle = if config.dupradar.boxplot {
                    let dm_ref = &dup_matrix;
                    let path = &boxplot_path;
                    Some(s.spawn(move || plots::duprate_boxplot(dm_ref, path)))
                } else {
                    None
                };

                // Histogram
                let histogram_handle = if config.dupradar.expression_histogram {
                    let dm_ref = &dup_matrix;
                    let path = &histogram_path;
                    Some(s.spawn(move || plots::expression_histogram(dm_ref, path)))
                } else {
                    None
                };

                // Collect results
                if let Some(handle) = density_handle {
                    handle
                        .join()
                        .expect("density scatter plot thread panicked")?;
                }
                if let Some(handle) = boxplot_handle {
                    handle.join().expect("boxplot thread panicked")?;
                }
                if let Some(handle) = histogram_handle {
                    handle.join().expect("histogram thread panicked")?;
                }

                Ok(())
            })?;

            info!(
                "Plots generated in {:.2}s",
                plot_start.elapsed().as_secs_f64()
            );
        }

        // Write MultiQC-compatible output files
        if let Some(ref fit) = fit_ok {
            if config.dupradar.multiqc_intercept {
                let mqc_intercept_path = outdir.join(format!("{}_dup_intercept_mqc.txt", bam_stem));
                plots::write_mqc_intercept(fit, bam_stem, &mqc_intercept_path)?;
            }

            if config.dupradar.multiqc_curve {
                let mqc_curve_path =
                    outdir.join(format!("{}_duprateExpDensCurve_mqc.txt", bam_stem));
                plots::write_mqc_curve(fit, &dup_matrix, &mqc_curve_path)?;
            }

            if config.dupradar.multiqc_intercept || config.dupradar.multiqc_curve {
                info!("MultiQC output files written");
            }
        }

        // Print summary statistics
        info!("--- dupRadar Summary ---");
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
    }

    info!("Total runtime: {:.2}s", start.elapsed().as_secs_f64());

    Ok(())
}
