//! Qualimap `bamqc`-compatible output files.
//!
//! Writes `genome_results.txt` plus the `raw_data_qualimapReport/` tables that
//! MultiQC's Qualimap BamQC module parses.

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use log::debug;

use super::accumulator::{BamqcResult, COVERAGE_THRESHOLDS};

/// Write all bamqc output files under `outdir`.
///
/// Produces:
/// - `genome_results.txt`
/// - `raw_data_qualimapReport/coverage_histogram.txt`
/// - `raw_data_qualimapReport/genome_fraction_coverage.txt`
/// - `raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt`
/// - `raw_data_qualimapReport/mapped_reads_nucleotide_content.txt`
/// - `raw_data_qualimapReport/insert_size_histogram.txt`
///
/// # Arguments
/// * `result` - The accumulated metrics
/// * `sample_name` - Sample name recorded in the report
/// * `bam_path` - Input path recorded in the report
/// * `outdir` - Directory to write into
///
/// # Returns
/// The paths written, in the order listed above
pub fn write_all(
    result: &BamqcResult,
    sample_name: &str,
    bam_path: &str,
    outdir: &Path,
) -> Result<Vec<std::path::PathBuf>> {
    std::fs::create_dir_all(outdir)
        .with_context(|| format!("Failed to create output directory: {}", outdir.display()))?;
    let raw_dir = outdir.join("raw_data_qualimapReport");
    std::fs::create_dir_all(&raw_dir)
        .with_context(|| format!("Failed to create directory: {}", raw_dir.display()))?;

    let mut written = Vec::new();

    let genome_results = outdir.join("genome_results.txt");
    write_genome_results(result, sample_name, bam_path, &genome_results)?;
    written.push(genome_results);

    let cov_hist = raw_dir.join("coverage_histogram.txt");
    write_coverage_histogram(result, &cov_hist)?;
    written.push(cov_hist);

    let genome_fraction = raw_dir.join("genome_fraction_coverage.txt");
    write_genome_fraction(result, &genome_fraction)?;
    written.push(genome_fraction);

    let gc = raw_dir.join("mapped_reads_gc-content_distribution.txt");
    write_gc_distribution(result, &gc)?;
    written.push(gc);

    let nucleotide = raw_dir.join("mapped_reads_nucleotide_content.txt");
    write_nucleotide_content(result, &nucleotide)?;
    written.push(nucleotide);

    let insert = raw_dir.join("insert_size_histogram.txt");
    write_insert_size_histogram(result, &insert)?;
    written.push(insert);

    debug!("Wrote bamqc outputs to {}", outdir.display());
    Ok(written)
}

/// Format a count with Qualimap's thousands separators (e.g. `1,234,567`).
fn commas(value: u64) -> String {
    let digits = value.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, ch) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(ch);
    }
    out
}

/// Percentage of `part` in `whole`, guarding against a zero denominator.
fn pct(part: u64, whole: u64) -> f64 {
    if whole == 0 {
        0.0
    } else {
        100.0 * part as f64 / whole as f64
    }
}

/// Write the main `genome_results.txt` summary.
fn write_genome_results(
    result: &BamqcResult,
    sample_name: &str,
    bam_path: &str,
    path: &Path,
) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .with_context(|| format!("Failed to create {}", path.display()))?;

    writeln!(out, "BamQC report")?;
    writeln!(out, "-----------------------------------")?;
    writeln!(out)?;
    writeln!(out, ">>>>>>> Input")?;
    writeln!(out)?;
    writeln!(out, "     bam file = {bam_path}")?;
    writeln!(out, "     outfile = genome_results.txt")?;
    writeln!(out, "     sample = {sample_name}")?;
    writeln!(out)?;
    writeln!(out, ">>>>>>> Reference")?;
    writeln!(out)?;
    writeln!(
        out,
        "     number of bases = {} bp",
        commas(result.reference_bases)
    )?;
    writeln!(out, "     number of contigs = {}", result.num_contigs)?;
    writeln!(out)?;
    writeln!(out, ">>>>>>> Globals")?;
    writeln!(out)?;
    writeln!(out, "     number of reads = {}", commas(result.total_reads))?;
    writeln!(
        out,
        "     number of mapped reads = {} ({:.2}%)",
        commas(result.mapped_reads),
        pct(result.mapped_reads, result.total_reads)
    )?;
    writeln!(
        out,
        "     number of secondary alignments = {}",
        commas(result.secondary)
    )?;
    writeln!(
        out,
        "     number of supplementary alignments = {}",
        commas(result.supplementary)
    )?;
    writeln!(
        out,
        "     number of mapped paired reads (first in pair) = {}",
        commas(result.first_in_pair)
    )?;
    writeln!(
        out,
        "     number of mapped paired reads (second in pair) = {}",
        commas(result.second_in_pair)
    )?;
    writeln!(
        out,
        "     number of mapped paired reads (both in pair) = {}",
        commas(result.both_mates_mapped)
    )?;
    writeln!(
        out,
        "     number of mapped paired reads (singletons) = {}",
        commas(result.singletons)
    )?;
    writeln!(
        out,
        "     number of duplicated reads (flagged) = {}",
        commas(result.duplicates)
    )?;
    writeln!(
        out,
        "     number of QC-failed reads = {}",
        commas(result.qc_failed)
    )?;
    writeln!(out)?;
    writeln!(
        out,
        "     number of mapped bases = {} bp",
        commas(result.aligned_bases)
    )?;
    writeln!(
        out,
        "     number of sequenced bases = {} bp",
        commas(result.sequenced_bases)
    )?;
    writeln!(out)?;
    writeln!(out, "     read min length = {}", result.read_length_min)?;
    writeln!(out, "     read max length = {}", result.read_length_max)?;
    writeln!(
        out,
        "     read mean length = {:.2}",
        result.read_length_mean
    )?;
    writeln!(out)?;
    writeln!(out, ">>>>>>> Insert size")?;
    writeln!(out)?;
    writeln!(out, "     mean insert size = {:.4}", result.insert_mean)?;
    writeln!(out, "     std insert size = {:.4}", result.insert_std)?;
    writeln!(out, "     median insert size = {}", result.insert_median)?;
    writeln!(out)?;
    writeln!(out, ">>>>>>> Mapping quality")?;
    writeln!(out)?;
    writeln!(out, "     mean mapping quality = {:.4}", result.mapq_mean)?;
    writeln!(out)?;
    writeln!(out, ">>>>>>> ACTG content")?;
    writeln!(out)?;
    let acgt: u64 = result.base_counts[..4].iter().sum();
    let labels = ["A", "C", "G", "T", "N"];
    for (i, label) in labels.iter().enumerate() {
        if i < 4 {
            writeln!(
                out,
                "     number of {}'s = {} bp ({:.2}%)",
                label,
                commas(result.base_counts[i]),
                pct(result.base_counts[i], acgt)
            )?;
        } else {
            writeln!(
                out,
                "     number of {}'s = {} bp",
                label,
                commas(result.base_counts[i])
            )?;
        }
    }
    writeln!(out, "     GC percentage = {:.2}%", result.gc_percentage())?;
    writeln!(out)?;
    writeln!(out, ">>>>>>> Coverage")?;
    writeln!(out)?;
    writeln!(out, "     mean coverageData = {:.4}X", result.mean_coverage)?;
    writeln!(out, "     std coverageData = {:.4}X", result.std_coverage)?;
    writeln!(out)?;
    for threshold in COVERAGE_THRESHOLDS {
        writeln!(
            out,
            "     There is a {:.2}% of reference with a coverageData >= {}X",
            result.genome_fraction_at(threshold),
            threshold
        )?;
    }
    writeln!(out)?;
    writeln!(out, ">>>>>>> Coverage per contig")?;
    writeln!(out)?;
    for contig in &result.contigs {
        writeln!(
            out,
            "\t{}\t{}\t{}\t{:.10}\t{:.10}",
            contig.name,
            contig.length,
            contig.mapped_bases,
            contig.mean_coverage,
            contig.std_coverage
        )?;
    }

    Ok(())
}

/// Write the depth histogram (`coverage_histogram.txt`).
fn write_coverage_histogram(result: &BamqcResult, path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .with_context(|| format!("Failed to create {}", path.display()))?;
    writeln!(out, "#Coverage\tNumber of genomic locations")?;

    let mut depths: Vec<(&u32, &u64)> = result.coverage_hist.iter().collect();
    depths.sort_unstable();
    for (&depth, &count) in depths {
        writeln!(out, "{}.0\t{}.0", depth, count)?;
    }
    Ok(())
}

/// Write the cumulative genome fraction curve (`genome_fraction_coverage.txt`).
fn write_genome_fraction(result: &BamqcResult, path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .with_context(|| format!("Failed to create {}", path.display()))?;
    writeln!(out, "#Coverage (X)\tPercentage of reference (%)")?;

    let max_depth = result.coverage_hist.keys().copied().max().unwrap_or(0);
    // Cap the reported curve the way Qualimap does, so a handful of very deep
    // positions cannot produce a multi-million-row file.
    let limit = max_depth.min(1000);
    for depth in 1..=limit.max(1) {
        writeln!(out, "{}.0\t{:.10}", depth, result.genome_fraction_at(depth))?;
    }
    Ok(())
}

/// Write the read GC distribution (`mapped_reads_gc-content_distribution.txt`).
fn write_gc_distribution(result: &BamqcResult, path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .with_context(|| format!("Failed to create {}", path.display()))?;
    writeln!(out, "#GC Content (%)\tSample")?;

    let total: u64 = result.gc_hist.iter().sum();
    for (pct_bin, &count) in result.gc_hist.iter().enumerate() {
        let fraction = if total == 0 {
            0.0
        } else {
            count as f64 / total as f64
        };
        writeln!(out, "{}.0\t{:.10}", pct_bin, fraction)?;
    }
    Ok(())
}

/// Write the overall base composition (`mapped_reads_nucleotide_content.txt`).
fn write_nucleotide_content(result: &BamqcResult, path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .with_context(|| format!("Failed to create {}", path.display()))?;
    writeln!(out, "#Base\tCount\tFraction")?;

    let total: u64 = result.base_counts.iter().sum();
    for (label, &count) in ["A", "C", "G", "T", "N"].iter().zip(&result.base_counts) {
        let fraction = if total == 0 {
            0.0
        } else {
            count as f64 / total as f64
        };
        writeln!(out, "{label}\t{count}\t{fraction:.10}")?;
    }
    Ok(())
}

/// Write the insert size histogram (`insert_size_histogram.txt`).
fn write_insert_size_histogram(result: &BamqcResult, path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .with_context(|| format!("Failed to create {}", path.display()))?;
    writeln!(out, "#Insert size\tNumber of reads")?;
    for (&size, &count) in &result.insert_hist {
        writeln!(out, "{}.0\t{}.0", size, count)?;
    }
    Ok(())
}

// ===================================================================
// Tests
// ===================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_commas() {
        assert_eq!(commas(0), "0");
        assert_eq!(commas(999), "999");
        assert_eq!(commas(1000), "1,000");
        assert_eq!(commas(1234567), "1,234,567");
    }
}
