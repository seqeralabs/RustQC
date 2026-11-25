//! Qualimap RNA-Seq QC output generation.
//!
//! Produces `rnaseq_qc_results.txt` and coverage profile TSVs in exact
//! Qualimap-compatible format, parseable by MultiQC.

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use log::debug;

use super::coverage::TranscriptCoverage;
use super::format::{
    format_decimal, format_long, format_percentage, java_double_to_string, percentage_of,
};
use super::index::QualimapIndex;
use super::plots;
use super::QualimapResult;
use crate::Strandedness;

// ============================= Constants =======================================

/// Number of percentile bins for coverage profiles.
const NUM_BINS: usize = 100;

/// Number of transcripts for the "high" expression tier.
const HIGH_TIER_COUNT: usize = 500;

/// Number of transcripts for the "low" expression tier.
const LOW_TIER_COUNT: usize = 500;

/// Number of top transcripts used for 5'-3' bias calculation, as Qualimap
/// prints it. Its selection code actually keeps the top 1001 (see
/// [`compute_bias`]).
pub(super) const NUM_TRANSCRIPTS_FOR_BIAS: usize = 1000;

/// Number of bases at 5'/3' ends used for bias calculation.
pub(super) const NUM_PRIME_BASES: usize = 100;

/// Minimum transcript length for bias calculation.
const MIN_TRANSCRIPT_LENGTH_FOR_BIAS: usize = 500;

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
/// Follows Qualimap's `TranscriptDataHandler.calculateCoverageBias()` and
/// `pickTranscripts()` (themselves ported from Picard's `RnaMetricsCollector`):
///
/// 1. Per gene, pick the transcript with the highest mean coverage among
///    those at least 500 bp long with mean coverage >= 1. Transcripts are
///    visited in Java `HashMap` order (`order`) and only a strictly higher
///    mean replaces the current best, so on ties the first one visited wins.
/// 2. Keep every picked transcript whose mean is at least the 1001st highest
///    mean (`coverages[max(0, n - 1001)]` in the sorted means), so up to 1001
///    transcripts plus ties, not 1000.
/// 3. For each, with the coverage reversed for minus-strand genes:
///    5' bias = mean(first 100 bases) / mean, 3' bias = mean(last 100) / mean,
///    5'-3' bias = mean(first 100) / mean(last 100). No value is filtered:
///    a transcript with no coverage in either end gives 0/0 = NaN.
/// 4. Report `StatUtils.percentile(values, 50)` of each list (Commons Math
///    2.0, see [`commons_math_median`]).
fn compute_bias(entries: &[TranscriptCoverageEntry], order: &[usize]) -> (f64, f64, f64) {
    use std::collections::HashMap;

    // Step 1: best transcript per gene.
    let mut best_per_gene: HashMap<u32, usize> = HashMap::new();
    for &idx in order {
        let entry = &entries[idx];
        if entry.coverage.len() < MIN_TRANSCRIPT_LENGTH_FOR_BIAS || entry.mean_coverage < 1.0 {
            continue;
        }
        let better = match best_per_gene.get(&entry.gene_idx) {
            None => true,
            Some(&best) => entry.mean_coverage > entries[best].mean_coverage,
        };
        if better {
            best_per_gene.insert(entry.gene_idx, idx);
        }
    }
    if best_per_gene.is_empty() {
        return (f64::NAN, f64::NAN, f64::NAN);
    }

    // Step 2: keep transcripts at or above the 1001st best mean coverage.
    let mut means: Vec<f64> = best_per_gene
        .values()
        .map(|&i| entries[i].mean_coverage)
        .collect();
    means.sort_by(java_double_cmp);
    let min_mean = means[means.len().saturating_sub(NUM_TRANSCRIPTS_FOR_BIAS + 1)];

    let mut five_prime_biases = Vec::new();
    let mut three_prime_biases = Vec::new();
    let mut five_three_biases = Vec::new();

    for &idx in best_per_gene.values() {
        let entry = &entries[idx];
        if entry.mean_coverage < min_mean {
            continue;
        }
        let cov = &entry.coverage;
        let len = cov.len();
        // Picard's MathUtil.mean: sequential sum / count. For minus-strand
        // genes the array is reversed, so the 5' end is at the end.
        let mean_of = |range: std::ops::Range<usize>| -> f64 {
            let n = range.len() as f64;
            range.map(|i| cov[i] as f64).sum::<f64>() / n
        };
        let (five_mean, three_mean) = if entry.strand == '-' {
            (
                mean_of(len - NUM_PRIME_BASES..len),
                mean_of(0..NUM_PRIME_BASES),
            )
        } else {
            (
                mean_of(0..NUM_PRIME_BASES),
                mean_of(len - NUM_PRIME_BASES..len),
            )
        };
        let whole_mean = mean_of(0..len);

        five_prime_biases.push(five_mean / whole_mean);
        three_prime_biases.push(three_mean / whole_mean);
        five_three_biases.push(five_mean / three_mean);
    }

    debug!(
        "Qualimap 5'-3' bias: {} transcripts, {} with no coverage at either end (NaN)",
        five_three_biases.len(),
        five_three_biases.iter().filter(|v| v.is_nan()).count()
    );

    (
        commons_math_median(&mut five_prime_biases),
        commons_math_median(&mut three_prime_biases),
        commons_math_median(&mut five_three_biases),
    )
}

/// Order doubles as Java's `Arrays.sort(double[])` does: NaN sorts after
/// +Infinity (whatever its sign bit) and -0.0 before 0.0.
fn java_double_cmp(a: &f64, b: &f64) -> std::cmp::Ordering {
    match (a.is_nan(), b.is_nan()) {
        (false, false) => a.total_cmp(b),
        (a_nan, b_nan) => a_nan.cmp(&b_nan),
    }
}

/// Commons Math 2.0 `StatUtils.percentile(values, 50)`, as Qualimap uses it.
///
/// The values are sorted with `Arrays.sort` (NaN last, not removed), then
/// `pos = p * (n + 1) / 100` is interpolated between its neighbours as
/// `lower + (pos - floor(pos)) * (upper - lower)`.
fn commons_math_median(values: &mut [f64]) -> f64 {
    const P: f64 = 50.0;
    let length = values.len();
    match length {
        0 => return f64::NAN,
        1 => return values[0],
        _ => {}
    }
    let n = length as f64;
    let pos = P * (n + 1.0) / 100.0;
    let fpos = pos.floor();
    let int_pos = fpos as usize;
    let dif = pos - fpos;
    values.sort_by(java_double_cmp);
    if pos < 1.0 {
        return values[0];
    }
    if pos >= n {
        return values[length - 1];
    }
    let lower = values[int_pos - 1];
    let upper = values[int_pos];
    lower + dif * (upper - lower)
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
/// Qualimap writes both columns by string concatenation, i.e. with
/// `Double.toString` (`BamQCChart.exportData`).
fn write_coverage_profile(profile: &[f64; NUM_BINS], path: &Path) -> Result<()> {
    use std::fs::File;

    let mut f = File::create(path)
        .with_context(|| format!("Failed to create coverage profile: {}", path.display()))?;

    writeln!(f, "#Transcript position\tTranscript coverage profile")?;

    for (i, &val) in profile.iter().enumerate() {
        writeln!(
            f,
            "{}\t{}",
            java_double_to_string(i as f64),
            java_double_to_string(val)
        )?;
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
#[allow(clippy::too_many_arguments)]
pub fn write_qualimap_results(
    result: &QualimapResult,
    index: &QualimapIndex,
    bam_path: &str,
    gtf_path: &str,
    stranded: Strandedness,
    output_dir: &Path,
    sample_name: &str,
) -> Result<()> {
    // Auto-detect paired mode from the data
    let paired = result.left_proper_in_pair > 0 || result.right_proper_in_pair > 0;
    use std::fs;

    fs::create_dir_all(output_dir)
        .with_context(|| format!("Failed to create output dir: {}", output_dir.display()))?;

    // --- Build per-transcript coverage entries ---
    let entries = build_transcript_entries(&result.transcript_coverage_raw, index);

    // --- Compute coverage profiles ---
    // Qualimap uses a TreeMap<Double, int[]> keyed by mean coverage for profile
    // computation. Java's TreeMap deduplicates by key: if two transcripts have
    // the exact same mean coverage (as a double), only the last one inserted
    // survives. We replicate this by using a BTreeMap keyed by the f64 bit
    // representation (which preserves the same total ordering as Java's
    // Double.compareTo). This produces identical profile values.
    use std::collections::BTreeMap;
    let mut deduped_map: BTreeMap<u64, usize> = BTreeMap::new();
    // Qualimap iterates genes via HashMap<String,Gene>, then transcripts via
    // the Gene's internal HashMap<String,Transcript>. To match this iteration
    // order exactly (and thus match which transcript "wins" in the TreeMap when
    // two transcripts share the same mean coverage), we simulate Java's
    // HashMap iteration: iterate buckets 0..capacity-1, within each bucket
    // entries are in reverse insertion order (Java 7/pre-8 HashMap).
    // We iterate entries sorted by java_hashmap_order(gene_id, transcript_id)
    // and insert into the BTreeMap with last-wins (matching TreeMap.put).
    let java_ordered_indices = java_hashmap_order_indices(&entries, index);
    for &idx in &java_ordered_indices {
        let entry = &entries[idx];
        if entry.mean_coverage > 0.0 {
            let key = entry.mean_coverage.to_bits();
            deduped_map.insert(key, idx);
        }
    }
    // Collect deduplicated entries in ascending mean-coverage order (BTreeMap iteration order)
    let sorted_by_mean: Vec<&TranscriptCoverageEntry> =
        deduped_map.values().map(|&idx| &entries[idx]).collect();

    // Log dedup collision details for debugging
    {
        let active_count = entries.iter().filter(|e| e.mean_coverage > 0.0).count();
        let dedup_count = sorted_by_mean.len();
        let num_collisions = active_count - dedup_count;
        log::debug!(
            "QM_PROFILE: {} total entries, {} active (mean>0), {} after dedup (TreeMap emulation), {} collisions",
            entries.len(),
            active_count,
            dedup_count,
            num_collisions
        );
        if num_collisions > 0 {
            // Report the mean values that collided
            let mut mean_counts: std::collections::HashMap<u64, Vec<usize>> =
                std::collections::HashMap::new();
            for (idx, entry) in entries.iter().enumerate() {
                if entry.mean_coverage > 0.0 {
                    mean_counts
                        .entry(entry.mean_coverage.to_bits())
                        .or_default()
                        .push(idx);
                }
            }
            for (bits, indices) in &mean_counts {
                if indices.len() > 1 {
                    let mean_val = f64::from_bits(*bits);
                    log::debug!(
                        "  COLLISION: mean={:.15} ({} transcripts, indices: {:?})",
                        mean_val,
                        indices.len(),
                        indices
                    );
                }
            }
        }
    }

    // Total: all deduplicated entries (ascending order, like Qualimap)
    let total_profile = compute_coverage_profile(&sorted_by_mean);

    // Low: bottom N transcripts (ascending order — first N from BTreeMap)
    let low_entries: Vec<&TranscriptCoverageEntry> = sorted_by_mean
        .iter()
        .take(LOW_TIER_COUNT)
        .copied()
        .collect();
    let low_profile = compute_coverage_profile(&low_entries);

    // High: top N transcripts (descending order — last N from BTreeMap, reversed)
    let high_entries: Vec<&TranscriptCoverageEntry> = sorted_by_mean
        .iter()
        .rev()
        .take(HIGH_TIER_COUNT)
        .copied()
        .collect();
    let high_profile = compute_coverage_profile(&high_entries);

    // --- Compute bias ---
    // TranscriptDataHandler.calculateCoverageBias returns early when no
    // transcript received any coverage, leaving the three medians at 0.0.
    let (five_bias, three_bias, five_three_bias) = if result
        .transcript_coverage_raw
        .num_transcripts_with_coverage()
        == 0
    {
        (0.0, 0.0, 0.0)
    } else {
        compute_bias(&entries, &java_ordered_indices)
    };

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

    debug!("Wrote coverage profiles to {}", raw_data_dir.display());

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

    debug!("Wrote {}", results_path.display());

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

    // Junction Analysis pie chart — uses Qualimap-native junction classification.
    // Novel = total N-operations - known - partly_known (matching Qualimap Java).
    let junction_counts = if result.reads_at_junctions > 0 {
        let known = result.known_junction_events;
        let partly_known = result.partly_known_junction_events;
        let novel = result.reads_at_junctions - known - partly_known;
        plots::junction_analysis_plot(
            known,
            partly_known,
            novel,
            sample_name,
            &images_dir.join("Junction Analysis.png"),
        )?;
        Some((known, partly_known, novel))
    } else {
        None
    };

    debug!("Wrote Qualimap plots to {}", images_dir.display());

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
        seq_not_found: result.seq_not_found,
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
    stranded: Strandedness,
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
        Strandedness::Forward => "strand-specific-forward",
        Strandedness::Reverse => "strand-specific-reverse",
        Strandedness::Unstranded => "non-strand-specific",
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
            format_long(result.left_proper_in_pair),
            format_long(result.right_proper_in_pair)
        )?;
        writeln!(
            f,
            "    read pairs aligned  = {}",
            format_long(result.both_proper_in_pair / 2)
        )?;
    } else {
        writeln!(
            f,
            "    reads aligned  = {}",
            format_long(result.primary_alignments)
        )?;
    }
    writeln!(
        f,
        "    total alignments = {}",
        format_long(result.primary_alignments + result.secondary_alignments)
    )?;
    writeln!(
        f,
        "    secondary alignments = {}",
        format_long(result.secondary_alignments)
    )?;
    writeln!(
        f,
        "    non-unique alignments = {}",
        format_long(result.alignment_not_unique)
    )?;
    writeln!(
        f,
        "    aligned to genes  = {}",
        format_long(result.exonic_reads)
    )?;
    writeln!(
        f,
        "    ambiguous alignments = {}",
        format_long(result.ambiguous_reads)
    )?;
    writeln!(
        f,
        "    no feature assigned = {}",
        format_long(result.no_feature)
    )?;
    writeln!(f, "    not aligned = {}", format_long(result.not_aligned))?;

    // SSP estimation: only when the protocol is not strand-specific.
    if stranded == Strandedness::Unstranded {
        let (fwd, rev) = super::strand_estimation(result.ssp_fwd, result.ssp_rev);
        writeln!(
            f,
            "    SSP estimation (fwd/rev) = {} / {}",
            format_decimal(fwd),
            format_decimal(rev)
        )?;
    }
    writeln!(f)?;
    writeln!(f)?;

    // >>>>>>> Reads genomic origin
    writeln!(f, ">>>>>>> Reads genomic origin")?;
    writeln!(f)?;

    // Qualimap divides by exonic + "no feature" reads (RNASeqQCAnalysis.writeReport).
    let total_reads = result.exonic_reads + result.no_feature;
    for (label, count) in [
        ("exonic = ", result.exonic_reads),
        ("intronic =", result.intronic_reads),
        ("intergenic =", result.intergenic_reads),
        ("overlapping exon =", result.overlapping_exon_reads),
    ] {
        writeln!(
            f,
            "    {label} {} ({})",
            format_long(count),
            format_percentage(percentage_of(count, total_reads))
        )?;
    }
    writeln!(f)?;
    writeln!(f)?;

    // >>>>>>> Transcript coverage profile
    writeln!(f, ">>>>>>> Transcript coverage profile")?;
    writeln!(f)?;
    writeln!(f, "    5' bias = {}", format_decimal(five_bias))?;
    writeln!(f, "    3' bias = {}", format_decimal(three_bias))?;
    writeln!(f, "    5'-3' bias = {}", format_decimal(five_three_bias))?;
    writeln!(f)?;
    writeln!(f)?;

    // >>>>>>> Junction analysis
    writeln!(f, ">>>>>>> Junction analysis")?;
    writeln!(f)?;
    writeln!(
        f,
        "    reads at junctions = {}",
        format_long(result.reads_at_junctions)
    )?;
    writeln!(f)?;

    // Junction motifs by percentage descending, ties broken by motif name
    for (motif, pct) in
        super::junction_motif_rows(&result.junction_motifs, result.reads_at_junctions)
    {
        writeln!(f, "    {} : {}", motif, format_percentage(pct))?;
    }

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

// ============================= Java HashMap Order Simulation ===================

/// Compute Java's `String.hashCode()` for a UTF-8 string.
///
/// Java's algorithm: `s[0]*31^(n-1) + s[1]*31^(n-2) + ... + s[n-1]`
/// using wrapping 32-bit signed arithmetic.
fn java_string_hashcode(s: &str) -> i32 {
    let mut h: i32 = 0;
    for &b in s.as_bytes() {
        h = h.wrapping_mul(31).wrapping_add(b as i32);
    }
    h
}

/// Simulate Java HashMap bucket assignment for a given capacity.
///
/// Java 7 HashMap uses: `indexFor(hash, length) = hash & (length - 1)`
/// where `hash` is the result of `hash(key.hashCode())` and `length` is
/// the table capacity (always a power of 2).
///
/// The supplemental hash function for Java 8+ is:
/// ```java
/// static final int hash(Object key) {
///     int h;
///     return (key == null) ? 0 : (h = key.hashCode()) ^ (h >>> 16);
/// }
/// ```
fn java8_hash(h: i32) -> i32 {
    h ^ ((h as u32 >> 16) as i32)
}

/// Compute the bucket index for a Java 8+ HashMap with the given capacity.
fn java8_bucket(key_hashcode: i32, capacity: usize) -> usize {
    let h = java8_hash(key_hashcode);
    (h as u32 as usize) & (capacity - 1)
}

/// Determine Java HashMap capacity for `n` insertions.
///
/// Java HashMap starts at capacity 16 (default) with load factor 0.75.
/// It doubles when `size >= capacity * 0.75` (threshold).
fn java_hashmap_capacity(n: usize) -> usize {
    let mut capacity: usize = 16;
    let load_factor = 0.75f64;
    while n > (capacity as f64 * load_factor) as usize {
        capacity *= 2;
    }
    capacity
}

/// Return transcript entry indices in the order Qualimap's Java code would
/// iterate them (HashMap<String,Gene> iteration order → Gene's internal
/// HashMap<String,Transcript> iteration order).
///
/// This simulates Java 7 HashMap iteration: iterate buckets 0..capacity-1,
/// within each bucket follow the linked list in insertion order
/// (oldest entry first, matching Java 8+'s tail-insertion).
fn java_hashmap_order_indices(
    entries: &[TranscriptCoverageEntry],
    index: &QualimapIndex,
) -> Vec<usize> {
    // Step 1: Group entry indices by gene_id, preserving within-gene order
    let mut gene_transcripts: indexmap::IndexMap<String, Vec<usize>> = indexmap::IndexMap::new();
    for (idx, _entry) in entries.iter().enumerate() {
        let tx_info = &index.transcripts[idx];
        gene_transcripts
            .entry(tx_info.gene_id.clone())
            .or_default()
            .push(idx);
    }

    // Step 2: Compute Java HashMap iteration order for gene_ids
    let gene_ids: Vec<String> = gene_transcripts.keys().cloned().collect();
    let gene_capacity = java_hashmap_capacity(gene_ids.len());

    // Build bucket lists.
    // Java 8+ HashMap inserts at TAIL of the bucket's linked list (or tree),
    // so entries within a bucket are in insertion order.
    let mut gene_buckets: Vec<Vec<String>> = vec![Vec::new(); gene_capacity];
    for gene_id in &gene_ids {
        let hc = java_string_hashcode(gene_id);
        let bucket = java8_bucket(hc, gene_capacity);
        // Java 8+ inserts at TAIL
        gene_buckets[bucket].push(gene_id.clone());
    }

    // Step 3: For each gene (in Java HashMap bucket order), compute transcript order
    let mut result = Vec::with_capacity(entries.len());

    for bucket in &gene_buckets {
        for gene_id in bucket {
            let tx_indices = &gene_transcripts[gene_id];

            // Compute Java HashMap order for transcript_ids within this gene
            let tx_ids: Vec<(&str, usize)> = tx_indices
                .iter()
                .map(|&idx| (index.transcripts[idx].transcript_id.as_str(), idx))
                .collect();

            let tx_capacity = java_hashmap_capacity(tx_ids.len());
            let mut tx_buckets: Vec<Vec<usize>> = vec![Vec::new(); tx_capacity];

            for &(tx_id, idx) in &tx_ids {
                let hc = java_string_hashcode(tx_id);
                let bucket_idx = java8_bucket(hc, tx_capacity);
                // Java 8+ inserts at TAIL
                tx_buckets[bucket_idx].push(idx);
            }

            for tx_bucket in &tx_buckets {
                for &idx in tx_bucket {
                    result.push(idx);
                }
            }
        }
    }

    result
}

// ============================= Tests ===========================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_commons_math_median() {
        assert_eq!(commons_math_median(&mut [1.0, 3.0, 2.0]), 2.0);
        assert_eq!(commons_math_median(&mut [1.0, 2.0, 3.0, 4.0]), 2.5);
        assert_eq!(commons_math_median(&mut [7.0]), 7.0);
        assert!(commons_math_median(&mut []).is_nan());
        // NaN is kept and sorts last (whatever its sign), like Java's
        // Arrays.sort, which moves the median up; +Infinity sorts before NaN.
        assert_eq!(commons_math_median(&mut [f64::NAN, 1.0, 2.0, 3.0]), 2.5);
        assert_eq!(
            commons_math_median(&mut [-f64::NAN, f64::INFINITY, 1.0, 2.0]),
            f64::INFINITY
        );
        // `lower + 0 * (NaN - lower)` is NaN, as in Java.
        assert!(commons_math_median(&mut [f64::NAN, 1.0, 2.0]).is_nan());
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
