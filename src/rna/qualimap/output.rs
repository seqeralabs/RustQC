//! Qualimap RNA-Seq QC output generation.
//!
//! Produces `rnaseq_qc_results.txt` and coverage profile TSVs in exact
//! Qualimap-compatible format, parseable by MultiQC.

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use log::info;

use super::coverage::TranscriptCoverage;
use super::index::QualimapIndex;
use super::plots;
use super::QualimapResult;

// ============================= Constants =======================================

/// Number of percentile bins for coverage profiles.
const NUM_BINS: usize = 100;

/// Number of transcripts for the "high" expression tier.
const HIGH_TIER_COUNT: usize = 500;

/// Number of transcripts for the "low" expression tier.
const LOW_TIER_COUNT: usize = 500;

/// Number of top transcripts used for 5'-3' bias calculation.
const NUM_TRANSCRIPTS_FOR_BIAS: usize = 1000;

/// Number of bases at 5'/3' ends used for bias calculation.
const NUM_PRIME_BASES: usize = 100;

/// Minimum transcript length for bias calculation.
const MIN_TRANSCRIPT_LENGTH_FOR_BIAS: usize = 500;

// ============================= Number Formatting ===============================

/// Format an integer with comma separators (e.g., 1234567 → "1,234,567").
use crate::io::format_with_commas;

/// Format a percentage with Qualimap-style precision.
///
/// Qualimap uses Java's default `String.format("%.2f", val)` but with
/// trailing-zero trimming on the decimal part. Looking at real output:
/// - 68.3% (not 68.30%)
/// - 29.69%
/// - 2.01%
/// - 9.31%
///
/// The pattern is: 2 decimal places, then strip trailing zeros (but keep
/// at least one decimal digit).
fn format_percentage(val: f64) -> String {
    let s = format!("{:.2}", val);
    // Strip trailing zeros after decimal point, but keep at least one decimal
    if s.contains('.') {
        let trimmed = s.trim_end_matches('0');
        let trimmed = if trimmed.ends_with('.') {
            // e.g., "68." → keep at least "68.0"
            // Actually Qualimap shows "68.3" not "68.0" — but if we get exactly
            // an integer like 68.00, it would show "68.0" based on the trim logic.
            // However looking at real data, Qualimap never shows ".0" — it shows
            // integers as "68" etc. Let's check: the real data shows "68.3%", "29.69%"
            // etc. When the value is exactly integer, Qualimap shows e.g. "9.31%".
            // There's no case of exactly-integer percentages in the reference.
            // Be safe: keep "X.0" format.
            &s[..trimmed.len() + 1]
        } else {
            trimmed
        };
        trimmed.to_string()
    } else {
        s
    }
}

/// Format a bias value with Qualimap-style precision.
///
/// Qualimap uses 2 decimal places with trailing-zero trimming.
/// Examples: 0.71, 0.57, 1.3, 1.47, 0.31, 0.22
fn format_bias(val: f64) -> String {
    let s = format!("{:.2}", val);
    if s.contains('.') {
        let trimmed = s.trim_end_matches('0');
        if trimmed.ends_with('.') {
            format!("{}.0", trimmed.trim_end_matches('.'))
        } else {
            trimmed.to_string()
        }
    } else {
        s
    }
}

// ============================= Coverage Profile ================================

/// Per-transcript coverage data used for profile computation.
struct TranscriptCoverageEntry {
    /// Per-base depth array.
    coverage: Vec<i32>,
    /// Mean coverage across the array.
    mean_coverage: f64,
    /// Strand: '+' or '-'.
    strand: char,
    /// Gene index for best-per-gene selection.
    gene_idx: u32,
    /// Flat transcript index (for diagnostics).
    #[allow(dead_code)]
    flat_idx: u32,
}

/// Compute coverage profiles using Qualimap's GenericHistogram algorithm.
///
/// For each transcript with mean coverage > 0:
/// 1. Bin the per-base coverage into `NUM_BINS` percentile bins
/// 2. Normalize by transcript length
/// 3. Accumulate across all transcripts
///
/// Returns a 100-element array of f64 values.
fn compute_coverage_profile(entries: &[&TranscriptCoverageEntry]) -> [f64; NUM_BINS] {
    let mut hist = [0.0f64; NUM_BINS];

    for entry in entries {
        let data = &entry.coverage;
        let len = data.len();
        if len == 0 {
            continue;
        }

        // Bin coverage into NUM_BINS (100) bins using step = ceil(len / 100).
        // Matches Qualimap's GenericHistogram behavior where each bin accumulates
        // coverage from `step` consecutive positions, normalized by transcript length.
        let step = len.div_ceil(NUM_BINS).max(1);
        let norm = len as f64;
        let mut bin = 0;
        let mut bin_sum = 0.0f64;
        let mut count = 0;

        for (i, &val) in data.iter().enumerate() {
            bin_sum += val as f64;
            count += 1;
            if count == step || i == len - 1 {
                hist[bin] += bin_sum / norm;
                bin += 1;
                bin_sum = 0.0;
                count = 0;
            }
        }
    }

    hist
}

/// Compute 5'-3' bias values from per-transcript coverage.
///
/// For each qualifying transcript (length >= 500, mean >= 1.0):
/// - Compute mean coverage of first `NUM_PRIME_BASES` bases (5' end)
/// - Compute mean coverage of last `NUM_PRIME_BASES` bases (3' end)
/// - Compute mean coverage of entire transcript
/// - 5' bias = 5' mean / whole mean
/// - 3' bias = 3' mean / whole mean
/// - 5'-3' bias = 5' mean / 3' mean
///
/// For negative-strand transcripts, the coverage array is reversed before
/// computing 5'/3' regions.
///
/// Returns (5' bias, 3' bias, 5'-3' bias) as medians across all qualifying
/// transcripts.
fn compute_bias(entries: &[TranscriptCoverageEntry]) -> (f64, f64, f64) {
    use std::collections::HashMap;

    // Step 1: Pick best transcript per gene (highest mean coverage).
    // Qualimap selects one transcript per gene, then filters by length >= 500
    // and mean >= 1.0, then takes the top 1000 by mean coverage.
    let mut best_per_gene: HashMap<u32, &TranscriptCoverageEntry> = HashMap::new();
    for entry in entries {
        if entry.coverage.len() < MIN_TRANSCRIPT_LENGTH_FOR_BIAS || entry.mean_coverage < 1.0 {
            continue;
        }
        let gene = entry.gene_idx;
        let current_best = best_per_gene.get(&gene);
        if current_best.is_none() || entry.mean_coverage > current_best.unwrap().mean_coverage {
            best_per_gene.insert(gene, entry);
        }
    }

    // Collect qualifying transcripts and sort by mean coverage descending
    let mut qualifying: Vec<(f64, &TranscriptCoverageEntry)> = best_per_gene
        .values()
        .map(|e| (e.mean_coverage, *e))
        .collect();

    if qualifying.is_empty() {
        return (f64::NAN, f64::NAN, f64::NAN);
    }

    // Sort by mean coverage descending, take top N
    qualifying.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    qualifying.truncate(NUM_TRANSCRIPTS_FOR_BIAS);

    let mut five_prime_biases = Vec::with_capacity(qualifying.len());
    let mut three_prime_biases = Vec::with_capacity(qualifying.len());
    let mut five_three_biases = Vec::with_capacity(qualifying.len());

    for (_, entry) in &qualifying {
        let cov = &entry.coverage;
        let len = cov.len();
        if len < NUM_PRIME_BASES {
            continue;
        }

        // For negative strand, reverse the array conceptually
        let (five_prime_sum, three_prime_sum, whole_sum) = if entry.strand == '-' {
            // 5' end is at the END of the array for negative strand
            let five_sum: f64 = cov[len - NUM_PRIME_BASES..].iter().map(|&v| v as f64).sum();
            let three_sum: f64 = cov[..NUM_PRIME_BASES].iter().map(|&v| v as f64).sum();
            let whole_sum: f64 = cov.iter().map(|&v| v as f64).sum();
            (five_sum, three_sum, whole_sum)
        } else {
            let five_sum: f64 = cov[..NUM_PRIME_BASES].iter().map(|&v| v as f64).sum();
            let three_sum: f64 = cov[len - NUM_PRIME_BASES..].iter().map(|&v| v as f64).sum();
            let whole_sum: f64 = cov.iter().map(|&v| v as f64).sum();
            (five_sum, three_sum, whole_sum)
        };

        let five_mean = five_prime_sum / NUM_PRIME_BASES as f64;
        let three_mean = three_prime_sum / NUM_PRIME_BASES as f64;
        let whole_mean = whole_sum / len as f64;

        // Qualimap pushes bias values for all qualifying transcripts unconditionally
        // (lines 247-249 of TranscriptDataHandler.java). Since pickTranscripts()
        // already ensures mean >= 1.0, whole_mean is always positive. The 5'-3'
        // ratio may produce infinity or NaN when three_mean == 0; Java's
        // StatUtils.percentile handles these via Arrays.sort which places NaN
        // at the end. We filter NaN values since they don't contribute meaningful
        // information to the median.
        five_prime_biases.push(five_mean / whole_mean);
        three_prime_biases.push(three_mean / whole_mean);
        let ratio = five_mean / three_mean;
        if !ratio.is_nan() {
            five_three_biases.push(ratio);
        }
    }

    let five = median(&mut five_prime_biases);
    let three = median(&mut three_prime_biases);
    let five_three = median(&mut five_three_biases);

    (five, three, five_three)
}

/// Compute the median of a slice of f64 values.
fn median(values: &mut [f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = values.len();
    if n.is_multiple_of(2) {
        (values[n / 2 - 1] + values[n / 2]) / 2.0
    } else {
        values[n / 2]
    }
}

// ============================= Coverage Histogram ==============================

/// Number of bins for the mean coverage histogram (0 to 50, inclusive).
const COVERAGE_HIST_BINS: usize = 51;

/// Compute the mean transcript coverage histogram.
///
/// For each transcript, computes `mean_coverage = sum(per_base_depth) / length`
/// and bins it into 51 buckets (0-50, with bin 50 collecting all transcripts
/// with mean coverage >= 50). Matches Qualimap's `computeMeanTranscriptCoverageHist`.
fn compute_mean_coverage_histogram(
    entries: &[TranscriptCoverageEntry],
) -> [u64; COVERAGE_HIST_BINS] {
    let mut hist = [0u64; COVERAGE_HIST_BINS];

    for entry in entries {
        if entry.coverage.is_empty() {
            hist[0] += 1;
            continue;
        }
        let mean_level = entry.mean_coverage as usize;
        let bin = mean_level.min(COVERAGE_HIST_BINS - 1);
        hist[bin] += 1;
    }

    hist
}

// ============================= Output Writers ==================================

/// Write the coverage profile TSV file.
///
/// Format: tab-separated, header line starting with `#`, 100 rows (0.0–99.0).
fn write_coverage_profile(profile: &[f64; NUM_BINS], path: &Path) -> Result<()> {
    use std::fs::File;

    let mut f = File::create(path)
        .with_context(|| format!("Failed to create coverage profile: {}", path.display()))?;

    writeln!(f, "#Transcript position\tTranscript coverage profile")?;

    for (i, &val) in profile.iter().enumerate() {
        writeln!(f, "{:.1}\t{}", i as f64, val)?;
    }

    Ok(())
}

/// Write the complete Qualimap RNA-Seq QC results: text files, coverage profiles,
/// and all plots (coverage profile, histogram, reads genomic origin, junction analysis).
///
/// Produces:
/// - `rnaseq_qc_results.txt` in exact Qualimap format
/// - `raw_data_qualimapReport/` with 3 coverage profile TSVs
/// - `images_qualimapReport/` with 6 PNG + SVG chart images
///
/// # Arguments
/// * `junction_counts` - Optional (known, partly_known, novel) event counts for the
///   Junction Analysis pie chart. Pass `None` to skip that chart (e.g., when junction
///   annotation is disabled).
#[allow(clippy::too_many_arguments)]
pub fn write_qualimap_results(
    result: &QualimapResult,
    index: &QualimapIndex,
    bam_path: &str,
    gtf_path: &str,
    stranded: u8,
    output_dir: &Path,
    sample_name: &str,
    junction_counts: Option<(u64, u64, u64)>,
) -> Result<()> {
    // Auto-detect paired mode from the data
    let paired = result.left_proper_in_pair > 0 || result.right_proper_in_pair > 0;
    use std::fs;

    fs::create_dir_all(output_dir)
        .with_context(|| format!("Failed to create output dir: {}", output_dir.display()))?;

    // --- Build per-transcript coverage entries ---
    let entries = build_transcript_entries(&result.transcript_coverage_raw, index);

    // --- Compute coverage profiles ---
    // Qualimap uses ALL transcripts with coverage > 0, sorted by mean coverage.
    // The profile is a SUM of per-transcript length-normalized histograms.
    let active_entries: Vec<&TranscriptCoverageEntry> =
        entries.iter().filter(|e| e.mean_coverage > 0.0).collect();
    let total_profile = compute_coverage_profile(&active_entries);

    log::debug!(
        "QM_PROFILE: {} total entries, {} active (mean>0)",
        entries.len(),
        active_entries.len()
    );

    // Sort by mean coverage for high/low tiers
    let mut sorted_by_mean: Vec<&TranscriptCoverageEntry> = active_entries.clone();
    sorted_by_mean.sort_by(|a, b| {
        a.mean_coverage
            .partial_cmp(&b.mean_coverage)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Low: bottom N transcripts (ascending order)
    let low_entries: Vec<&TranscriptCoverageEntry> = sorted_by_mean
        .iter()
        .take(LOW_TIER_COUNT)
        .copied()
        .collect();
    let low_profile = compute_coverage_profile(&low_entries);

    // High: top N transcripts (descending order)
    let high_entries: Vec<&TranscriptCoverageEntry> = sorted_by_mean
        .iter()
        .rev()
        .take(HIGH_TIER_COUNT)
        .copied()
        .collect();
    let high_profile = compute_coverage_profile(&high_entries);

    // --- Compute bias ---
    let (five_bias, three_bias, five_three_bias) = compute_bias(&entries);

    // --- Compute mean coverage histogram ---
    let coverage_histogram = compute_mean_coverage_histogram(&entries);

    // --- Write coverage profile TSVs ---
    let raw_data_dir = output_dir.join("raw_data_qualimapReport");
    fs::create_dir_all(&raw_data_dir)?;

    write_coverage_profile(
        &total_profile,
        &raw_data_dir.join("coverage_profile_along_genes_(total).txt"),
    )?;
    write_coverage_profile(
        &high_profile,
        &raw_data_dir.join("coverage_profile_along_genes_(high).txt"),
    )?;
    write_coverage_profile(
        &low_profile,
        &raw_data_dir.join("coverage_profile_along_genes_(low).txt"),
    )?;

    info!("Wrote coverage profiles to {}", raw_data_dir.display());

    // --- Write rnaseq_qc_results.txt ---
    let results_path = output_dir.join("rnaseq_qc_results.txt");
    write_results_file(
        result,
        bam_path,
        gtf_path,
        paired,
        stranded,
        five_bias,
        three_bias,
        five_three_bias,
        &results_path,
    )?;

    info!("Wrote {}", results_path.display());

    // --- Generate plots ---
    let images_dir = output_dir.join("images_qualimapReport");
    fs::create_dir_all(&images_dir)?;

    // Coverage Profile Along Genes (Total / High / Low)
    plots::coverage_profile_plot(
        &total_profile,
        "Coverage Profile Along Genes (Total)",
        sample_name,
        &images_dir.join("Coverage Profile Along Genes (Total).png"),
    )?;

    plots::coverage_profile_plot(
        &high_profile,
        "Coverage Profile Along Genes (High)",
        sample_name,
        &images_dir.join("Coverage Profile Along Genes (High).png"),
    )?;

    plots::coverage_profile_plot(
        &low_profile,
        "Coverage Profile Along Genes (Low)",
        sample_name,
        &images_dir.join("Coverage Profile Along Genes (Low).png"),
    )?;

    // Coverage Histogram (0-50X)
    plots::coverage_histogram_plot(
        &coverage_histogram,
        sample_name,
        &images_dir.join("Transcript coverage histogram.png"),
    )?;

    // Reads Genomic Origin pie chart
    plots::reads_genomic_origin_plot(
        result.exonic_reads,
        result.intronic_reads,
        result.intergenic_reads,
        sample_name,
        &images_dir.join("Reads Genomic Origin.png"),
    )?;

    // Junction Analysis pie chart (optional — requires non-zero junction counts)
    if let Some((known, partly_known, novel)) = junction_counts {
        if known + partly_known + novel > 0 {
            plots::junction_analysis_plot(
                known,
                partly_known,
                novel,
                sample_name,
                &images_dir.join("Junction Analysis.png"),
            )?;
        }
    }

    info!("Wrote Qualimap plots to {}", images_dir.display());

    // --- HTML report ---
    let report_data = super::report::ReportData {
        sample_name,
        bam_path,
        gtf_path,
        paired,
        stranded,
        left_proper: result.left_proper_in_pair,
        right_proper: result.right_proper_in_pair,
        both_proper: result.both_proper_in_pair,
        read_count: result.read_count,
        primary_alignments: result.primary_alignments,
        secondary_alignments: result.secondary_alignments,
        alignment_not_unique: result.alignment_not_unique,
        exonic_reads: result.exonic_reads,
        ambiguous_reads: result.ambiguous_reads,
        no_feature: result.no_feature,
        not_aligned: result.not_aligned,
        intronic_reads: result.intronic_reads,
        intergenic_reads: result.intergenic_reads,
        overlapping_exon_reads: result.overlapping_exon_reads,
        ssp_fwd: result.ssp_fwd,
        ssp_rev: result.ssp_rev,
        reads_at_junctions: result.reads_at_junctions,
        junction_motifs: &result.junction_motifs,
        five_bias,
        three_bias,
        five_three_bias,
        junction_counts,
    };
    super::report::write_html_report(&report_data, output_dir)?;

    Ok(())
}

/// Write the rnaseq_qc_results.txt file in exact Qualimap format.
#[allow(clippy::too_many_arguments)]
fn write_results_file(
    result: &QualimapResult,
    bam_path: &str,
    gtf_path: &str,
    paired: bool,
    stranded: u8,
    five_bias: f64,
    three_bias: f64,
    five_three_bias: f64,
    path: &Path,
) -> Result<()> {
    use std::fs::File;

    let mut f =
        File::create(path).with_context(|| format!("Failed to create {}", path.display()))?;

    // Header
    writeln!(f, "RNA-Seq QC report")?;
    writeln!(f, "-----------------------------------")?;
    writeln!(f)?;

    // >>>>>>> Input
    writeln!(f, ">>>>>>> Input")?;
    writeln!(f)?;
    writeln!(f, "    bam file = {bam_path}")?;
    writeln!(f, "    gff file = {gtf_path}")?;
    writeln!(f, "    counting algorithm = uniquely-mapped-reads")?;
    let protocol = match stranded {
        1 => "strand-specific-forward",
        2 => "strand-specific-reverse",
        _ => "non-strand-specific",
    };
    writeln!(f, "    protocol = {protocol}")?;
    writeln!(f, "    5'-3' bias region size = {NUM_PRIME_BASES}")?;
    writeln!(
        f,
        "    5'-3' bias number of top transcripts = {NUM_TRANSCRIPTS_FOR_BIAS}"
    )?;
    writeln!(f)?;
    writeln!(f)?;

    // >>>>>>> Reads alignment
    writeln!(f, ">>>>>>> Reads alignment")?;
    writeln!(f)?;

    if paired {
        writeln!(
            f,
            "    reads aligned (left/right) = {} / {}",
            format_with_commas(result.left_proper_in_pair),
            format_with_commas(result.right_proper_in_pair)
        )?;
        writeln!(
            f,
            "    read pairs aligned  = {}",
            format_with_commas(result.both_proper_in_pair / 2)
        )?;
    } else {
        writeln!(
            f,
            "    reads aligned = {}",
            format_with_commas(result.read_count)
        )?;
    }
    writeln!(
        f,
        "    total alignments = {}",
        format_with_commas(result.primary_alignments + result.secondary_alignments)
    )?;
    writeln!(
        f,
        "    secondary alignments = {}",
        format_with_commas(result.secondary_alignments)
    )?;
    writeln!(
        f,
        "    non-unique alignments = {}",
        format_with_commas(result.alignment_not_unique)
    )?;
    writeln!(
        f,
        "    aligned to genes  = {}",
        format_with_commas(result.exonic_reads)
    )?;
    writeln!(
        f,
        "    ambiguous alignments = {}",
        format_with_commas(result.ambiguous_reads)
    )?;
    writeln!(
        f,
        "    no feature assigned = {}",
        format_with_commas(result.no_feature)
    )?;
    writeln!(
        f,
        "    not aligned = {}",
        format_with_commas(result.not_aligned)
    )?;

    // SSP estimation — only written when protocol was not pre-specified (stranded == 0)
    // When the user specifies --stranded, Qualimap does not output this line.
    if stranded == 0 {
        let ssp_total = result.ssp_fwd + result.ssp_rev;
        if ssp_total > 0 {
            let fwd = result.ssp_fwd as f64 / ssp_total as f64;
            let rev = result.ssp_rev as f64 / ssp_total as f64;
            writeln!(f, "    SSP estimation (fwd/rev) = {:.2} / {:.2}", fwd, rev)?;
        }
    }
    writeln!(f)?;
    writeln!(f)?;

    // >>>>>>> Reads genomic origin
    writeln!(f, ">>>>>>> Reads genomic origin")?;
    writeln!(f)?;

    let total_classified = result.exonic_reads + result.intronic_reads + result.intergenic_reads;
    if total_classified > 0 {
        let exonic_pct = result.exonic_reads as f64 / total_classified as f64 * 100.0;
        let intronic_pct = result.intronic_reads as f64 / total_classified as f64 * 100.0;
        let intergenic_pct = result.intergenic_reads as f64 / total_classified as f64 * 100.0;

        // overlapping exon percentage is relative to total classified
        let overlapping_pct =
            result.overlapping_exon_reads as f64 / total_classified as f64 * 100.0;

        writeln!(
            f,
            "    exonic =  {} ({}%)",
            format_with_commas(result.exonic_reads),
            format_percentage(exonic_pct)
        )?;
        writeln!(
            f,
            "    intronic = {} ({}%)",
            format_with_commas(result.intronic_reads),
            format_percentage(intronic_pct)
        )?;
        writeln!(
            f,
            "    intergenic = {} ({}%)",
            format_with_commas(result.intergenic_reads),
            format_percentage(intergenic_pct)
        )?;
        writeln!(
            f,
            "    overlapping exon = {} ({}%)",
            format_with_commas(result.overlapping_exon_reads),
            format_percentage(overlapping_pct)
        )?;
    }
    writeln!(f)?;
    writeln!(f)?;

    // >>>>>>> Transcript coverage profile
    writeln!(f, ">>>>>>> Transcript coverage profile")?;
    writeln!(f)?;
    writeln!(f, "    5' bias = {}", format_bias(five_bias))?;
    writeln!(f, "    3' bias = {}", format_bias(three_bias))?;
    writeln!(f, "    5'-3' bias = {}", format_bias(five_three_bias))?;
    writeln!(f)?;
    writeln!(f)?;

    // >>>>>>> Junction analysis
    writeln!(f, ">>>>>>> Junction analysis")?;
    writeln!(f)?;
    writeln!(
        f,
        "    reads at junctions = {}",
        format_with_commas(result.reads_at_junctions)
    )?;
    writeln!(f)?;

    // Sort junction motifs by percentage descending
    if result.reads_at_junctions > 0 {
        let total_junctions = result.reads_at_junctions as f64;
        let mut motif_pcts: Vec<(String, f64)> = result
            .junction_motifs
            .iter()
            .map(|(motif, &count)| (motif.clone(), count as f64 * 100.0 / total_junctions))
            .collect();
        motif_pcts.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

        // Qualimap shows top 11 motifs (count <= 10 in their loop)
        for (motif, pct) in motif_pcts.iter().take(11) {
            writeln!(f, "    {} : {:.2}%", motif, pct)?;
        }
    }
    writeln!(f)?;

    Ok(())
}

/// Build transcript coverage entries from the raw per-transcript coverage.
///
/// Uses the raw `TranscriptCoverage` (keyed by flat transcript index) to avoid
/// the overhead and potential collisions of string-keyed lookup.
fn build_transcript_entries(
    coverage: &TranscriptCoverage,
    index: &QualimapIndex,
) -> Vec<TranscriptCoverageEntry> {
    let mut entries = Vec::new();

    for (flat_idx, tx_info) in index.transcripts.iter().enumerate() {
        let flat_idx = flat_idx as u32;
        match coverage.get(flat_idx) {
            Some(depth) => {
                let sum: f64 = depth.iter().map(|&v| v as f64).sum();
                let mean = if depth.is_empty() {
                    0.0
                } else {
                    sum / depth.len() as f64
                };
                entries.push(TranscriptCoverageEntry {
                    coverage: depth.to_vec(),
                    mean_coverage: mean,
                    strand: tx_info.strand,
                    gene_idx: tx_info.gene_idx,
                    flat_idx,
                });
            }
            None => {
                // Transcript with no coverage — still include for total count
                entries.push(TranscriptCoverageEntry {
                    coverage: vec![0; tx_info.exonic_length as usize],
                    mean_coverage: 0.0,
                    strand: tx_info.strand,
                    gene_idx: tx_info.gene_idx,
                    flat_idx,
                });
            }
        }
    }

    entries
}

// ============================= Tests ===========================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_with_commas() {
        assert_eq!(format_with_commas(0), "0");
        assert_eq!(format_with_commas(999), "999");
        assert_eq!(format_with_commas(1000), "1,000");
        assert_eq!(format_with_commas(1234567), "1,234,567");
        assert_eq!(format_with_commas(87559392), "87,559,392");
    }

    #[test]
    fn test_format_percentage() {
        assert_eq!(format_percentage(68.30), "68.3");
        assert_eq!(format_percentage(29.69), "29.69");
        assert_eq!(format_percentage(2.01), "2.01");
        assert_eq!(format_percentage(9.31), "9.31");
        assert_eq!(format_percentage(100.00), "100.0");
    }

    #[test]
    fn test_format_bias() {
        assert_eq!(format_bias(0.71), "0.71");
        assert_eq!(format_bias(0.57), "0.57");
        assert_eq!(format_bias(1.30), "1.3");
        assert_eq!(format_bias(1.47), "1.47");
        assert_eq!(format_bias(0.22), "0.22");
    }

    #[test]
    fn test_median() {
        assert_eq!(median(&mut [1.0, 3.0, 2.0]), 2.0);
        assert_eq!(median(&mut [1.0, 2.0, 3.0, 4.0]), 2.5);
        assert!(median(&mut []).is_nan());
    }

    /// Helper: build a single-transcript slice for compute_coverage_profile.
    fn single_entry(coverage: Vec<i32>) -> TranscriptCoverageEntry {
        let mean = if coverage.is_empty() {
            0.0
        } else {
            coverage.iter().map(|&v| v as f64).sum::<f64>() / coverage.len() as f64
        };
        TranscriptCoverageEntry {
            coverage,
            mean_coverage: mean,
            strand: '+',
            gene_idx: 0,
            flat_idx: 0,
        }
    }

    // ---- compute_coverage_profile tests matching Qualimap's GenericHistogram ----

    /// Uniform array of length 100: step=1, each bin = 1.0/100 = 0.01.
    #[test]
    fn test_coverage_profile_binning_uniform_100() {
        let entry = single_entry(vec![1; 100]);
        let profile = compute_coverage_profile(&[&entry]);
        for (i, &val) in profile.iter().enumerate() {
            assert!(
                (val - 0.01).abs() < 1e-10,
                "bin {i}: expected 0.01, got {val}"
            );
        }
    }

    /// Uniform array of length 200: step=2, each bin = 2.0/200 = 0.01.
    #[test]
    fn test_coverage_profile_binning_uniform_200() {
        let entry = single_entry(vec![1; 200]);
        let profile = compute_coverage_profile(&[&entry]);
        for (i, &val) in profile.iter().enumerate() {
            assert!(
                (val - 0.01).abs() < 1e-10,
                "bin {i}: expected 0.01, got {val}"
            );
        }
    }

    /// Uniform array of length 150: step=2, 75 bins filled (each ≈ 0.01333), bins 75-99 = 0.0.
    #[test]
    fn test_coverage_profile_binning_uniform_150() {
        let entry = single_entry(vec![1; 150]);
        let profile = compute_coverage_profile(&[&entry]);
        let expected_filled = 2.0 / 150.0; // 0.01333...
        for (i, &val) in profile.iter().enumerate() {
            if i < 75 {
                assert!(
                    (val - expected_filled).abs() < 1e-10,
                    "bin {i}: expected {expected_filled}, got {val}"
                );
            } else {
                assert!(
                    val.abs() < 1e-10,
                    "bin {i}: expected 0.0 (unfilled), got {val}"
                );
            }
        }
    }

    /// Ramp array [0, 1, ..., 99], length 100: step=1, bin i = i/100.0.
    #[test]
    fn test_coverage_profile_binning_ramp_100() {
        let coverage: Vec<i32> = (0..100).collect();
        let entry = single_entry(coverage);
        let profile = compute_coverage_profile(&[&entry]);
        for (i, &val) in profile.iter().enumerate() {
            let expected = i as f64 / 100.0;
            assert!(
                (val - expected).abs() < 1e-10,
                "bin {i}: expected {expected}, got {val}"
            );
        }
    }

    /// Short uniform array of length 50: step=1 (clamped), 50 bins filled (each = 1.0/50 = 0.02),
    /// bins 50-99 = 0.0.
    #[test]
    fn test_coverage_profile_binning_short_50() {
        let entry = single_entry(vec![1; 50]);
        let profile = compute_coverage_profile(&[&entry]);
        let expected_filled = 1.0 / 50.0; // 0.02
        for (i, &val) in profile.iter().enumerate() {
            if i < 50 {
                assert!(
                    (val - expected_filled).abs() < 1e-10,
                    "bin {i}: expected {expected_filled}, got {val}"
                );
            } else {
                assert!(
                    val.abs() < 1e-10,
                    "bin {i}: expected 0.0 (unfilled), got {val}"
                );
            }
        }
    }

    #[test]
    fn test_coverage_profile_single_transcript() {
        // Simple transcript: 200 bases, uniform coverage of 10
        let entry = TranscriptCoverageEntry {
            coverage: vec![10; 200],
            mean_coverage: 10.0,
            strand: '+',
            gene_idx: 0,
            flat_idx: 0,
        };
        let entries: Vec<&TranscriptCoverageEntry> = vec![&entry];
        let profile = compute_coverage_profile(&entries);

        // Each bin should have roughly 10.0 (normalized by length)
        // With 200 bases and 100 bins: step=2, each bin sums 2 values of 10 = 20
        // Divided by norm (200) = 0.1 per bin? No...
        // Actually: bin_coverage = sum of values in bin = 20
        // hist[bin] += bin_coverage / norm = 20 / 200 = 0.1
        for &val in profile.iter() {
            assert!((val - 0.1).abs() < 1e-10, "Expected 0.1, got {val}");
        }
    }
}
