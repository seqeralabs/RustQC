//! Junction saturation analysis for RNA-seq data.
//!
//! Reimplements RSeQC's `junction_saturation.py`: subsamples splice junction
//! observations at increasing percentages of total reads and reports how many
//! known / novel / total unique junctions are detected at each level.

use anyhow::{Context, Result};
use log::info;
use std::io::Write;

// ============================= Saturation result =============================

/// Results from junction saturation analysis at each sampling percentage.
#[derive(Debug)]
pub struct SaturationResult {
    /// Sampling percentages (e.g., 5, 10, ..., 100).
    pub percentages: Vec<u32>,
    /// Number of known junctions detected at each percentage.
    pub known_counts: Vec<usize>,
    /// Number of novel junctions detected at each percentage.
    pub novel_counts: Vec<usize>,
    /// Number of all unique junctions detected at each percentage.
    pub all_counts: Vec<usize>,
}

// ============================= Output formatting =============================

/// Write the junction saturation R script.
///
/// Produces an R script matching RSeQC's `junctionSaturation_plot.r` format.
pub fn write_r_script(result: &SaturationResult, prefix: &str) -> Result<()> {
    let r_path = format!("{prefix}.junctionSaturation_plot.r");
    let mut f =
        std::fs::File::create(&r_path).with_context(|| format!("creating R script: {r_path}"))?;

    let pdf_path = format!("{prefix}.junctionSaturation_plot.pdf");

    // x values
    let x_str: Vec<String> = result.percentages.iter().map(|p| p.to_string()).collect();

    // y = known, z = all, w = novel
    let y_str: Vec<String> = result.known_counts.iter().map(|c| c.to_string()).collect();
    let z_str: Vec<String> = result.all_counts.iter().map(|c| c.to_string()).collect();
    let w_str: Vec<String> = result.novel_counts.iter().map(|c| c.to_string()).collect();

    // m = max of last values / 1000, n = min of first values / 1000
    // R's max()/min() computes the result — we emit known, all, novel in that order
    // to match RSeQC's output format.
    let known_last_k = *result.known_counts.last().unwrap_or(&0) as i64 / 1000;
    let all_last_k = *result.all_counts.last().unwrap_or(&0) as i64 / 1000;
    let novel_last_k = *result.novel_counts.last().unwrap_or(&0) as i64 / 1000;
    let known_first_k = *result.known_counts.first().unwrap_or(&0) as i64 / 1000;
    let all_first_k = *result.all_counts.first().unwrap_or(&0) as i64 / 1000;
    let novel_first_k = *result.novel_counts.first().unwrap_or(&0) as i64 / 1000;

    writeln!(f, "pdf('{pdf_path}')")?;
    writeln!(f, "x=c({})", x_str.join(","))?;
    writeln!(f, "y=c({})", y_str.join(","))?;
    writeln!(f, "z=c({})", z_str.join(","))?;
    writeln!(f, "w=c({})", w_str.join(","))?;
    writeln!(f, "m=max({known_last_k},{all_last_k},{novel_last_k})")?;
    writeln!(f, "n=min({known_first_k},{all_first_k},{novel_first_k})")?;
    writeln!(
        f,
        "plot(x,z/1000,xlab='percent of total reads',ylab='Number of splicing junctions (x1000)',type='o',col='blue',ylim=c(n,m))"
    )?;
    writeln!(f, "points(x,y/1000,type='o',col='red')")?;
    writeln!(f, "points(x,w/1000,type='o',col='green')")?;
    writeln!(
        f,
        "legend(5,{}, legend=c(\"All junctions\",\"known junctions\", \"novel junctions\"),col=c(\"blue\",\"red\",\"green\"),lwd=1,pch=1)",
        all_last_k
    )?;
    writeln!(f, "dev.off()")?;

    info!("Wrote R script: {r_path}");
    Ok(())
}

/// Write a summary text file with junction saturation statistics.
pub fn write_summary(result: &SaturationResult, path: &str) -> Result<()> {
    let mut f = std::fs::File::create(path).with_context(|| format!("creating summary: {path}"))?;

    writeln!(
        f,
        "Percent\tAll_Junctions\tKnown_Junctions\tNovel_Junctions"
    )?;
    for i in 0..result.percentages.len() {
        writeln!(
            f,
            "{}\t{}\t{}\t{}",
            result.percentages[i],
            result.all_counts[i],
            result.known_counts[i],
            result.novel_counts[i],
        )?;
    }

    info!("Wrote summary: {path}");
    Ok(())
}

// ============================= Tests =============================

// Note: CIGAR intron extraction tests and BED junction parsing tests have been
// moved to the common module (rseqc::common::tests). Tests here focus on
// saturation-specific logic.
