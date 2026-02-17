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
fn format_with_commas(n: u64) -> String {
    let s = n.to_string();
    let bytes = s.as_bytes();
    let len = bytes.len();
    if len <= 3 {
        return s;
    }
    let mut result = String::with_capacity(len + (len - 1) / 3);
    for (i, &b) in bytes.iter().enumerate() {
        if i > 0 && (len - i).is_multiple_of(3) {
            result.push(',');
        }
        result.push(b as char);
    }
    result
}

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

        // Qualimap's GenericHistogram: normalize=true means divide by data.length
        let norm = len as f64;
        let step = len.div_ceil(NUM_BINS); // ceil division

        let mut bin_coverage = 0.0f64;
        let mut bin_index = 0usize;
        let mut count = step; // next boundary

        for (i, &val) in data.iter().enumerate() {
            bin_coverage += val as f64;

            if i + 1 == count || i == len - 1 {
                if bin_index < NUM_BINS {
                    hist[bin_index] += bin_coverage / norm;
                }
                bin_coverage = 0.0;
                bin_index += 1;
                count += step;
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

        if whole_mean > 0.0 {
            five_prime_biases.push(five_mean / whole_mean);
        }
        if whole_mean > 0.0 {
            three_prime_biases.push(three_mean / whole_mean);
        }
        if three_mean > 0.0 {
            five_three_biases.push(five_mean / three_mean);
        }
    }

    // Diagnostic logging for bias analysis
    let zero_five = five_prime_biases.iter().filter(|&&v| v == 0.0).count();
    let zero_three = three_prime_biases.iter().filter(|&&v| v == 0.0).count();
    log::debug!(
        "QM_BIAS qualifying={} zero_5'={} zero_3'={} total_entries={}",
        qualifying.len(),
        zero_five,
        zero_three,
        entries.len()
    );
    // Print details of first 5 zero-5' transcripts
    let mut zero_five_count = 0;
    for (mean, entry) in &qualifying {
        let cov = &entry.coverage;
        let len = cov.len();
        if len < NUM_PRIME_BASES {
            continue;
        }
        let first_100: f64 = if entry.strand == '-' {
            cov[len - NUM_PRIME_BASES..].iter().map(|&v| v as f64).sum()
        } else {
            cov[..NUM_PRIME_BASES].iter().map(|&v| v as f64).sum()
        };
        if first_100 == 0.0 && zero_five_count < 5 {
            let first_nonzero = cov.iter().position(|&v| v > 0).unwrap_or(len);
            let last_nonzero = cov.iter().rposition(|&v| v > 0).unwrap_or(0);
            log::debug!(
                "QM_BIAS_ZERO_5 gene_idx={} strand={} len={} mean={:.2} first_nonzero={} last_nonzero={} first_5cov=[{},{},{},{},{}] last_5cov=[{},{},{},{},{}]",
                entry.gene_idx,
                entry.strand,
                len,
                mean,
                first_nonzero,
                last_nonzero,
                cov[0], cov[1], cov[2], cov[3], cov[4],
                cov[len-5], cov[len-4], cov[len-3], cov[len-2], cov[len-1],
            );
            zero_five_count += 1;
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

/// Write the complete Qualimap RNA-Seq QC results file.
///
/// Produces `rnaseq_qc_results.txt` in exact Qualimap format with all sections:
/// Input, Reads alignment, Reads genomic origin, Transcript coverage profile,
/// and Junction analysis.
#[allow(clippy::too_many_arguments)]
pub fn write_qualimap_results(
    result: &QualimapResult,
    index: &QualimapIndex,
    bam_path: &str,
    gtf_path: &str,
    output_dir: &Path,
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
        five_bias,
        three_bias,
        five_three_bias,
        &results_path,
    )?;

    info!("Wrote {}", results_path.display());

    Ok(())
}

/// Write the rnaseq_qc_results.txt file in exact Qualimap format.
#[allow(clippy::too_many_arguments)]
fn write_results_file(
    result: &QualimapResult,
    bam_path: &str,
    gtf_path: &str,
    paired: bool,
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
    writeln!(f, "    protocol = non-strand-specific")?;
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

    // SSP estimation — from strand comparison during enclosure check
    let ssp_total = result.ssp_fwd + result.ssp_rev;
    if ssp_total > 0 {
        let fwd = result.ssp_fwd as f64 / ssp_total as f64;
        let rev = result.ssp_rev as f64 / ssp_total as f64;
        writeln!(f, "    SSP estimation (fwd/rev) = {:.2} / {:.2}", fwd, rev)?;
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
        assert_eq!(median(&mut vec![1.0, 3.0, 2.0]), 2.0);
        assert_eq!(median(&mut vec![1.0, 2.0, 3.0, 4.0]), 2.5);
        assert!(median(&mut vec![]).is_nan());
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
