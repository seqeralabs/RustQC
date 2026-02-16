//! Library complexity estimation (preseq lc_extrap reimplementation).
//!
//! Estimates the expected number of distinct molecules as a function of
//! sequencing depth using the Good-Toulmin rational function extrapolation
//! method, matching the behavior of preseq v3.

use anyhow::{bail, Context, Result};
use log::{debug, info, warn};
use rand::prelude::*;
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

use crate::config::PreseqConfig;

// ============================================================================
// Constants
// ============================================================================

/// Lanczos approximation coefficients (g=7, n=9) for ln_gamma.
#[allow(clippy::excessive_precision)]
const LANCZOS_COEFFS: [f64; 9] = [
    0.999_999_999_999_809_93,
    676.520_368_121_885_1,
    -1_259.139_216_722_403,
    771.323_428_777_653_1,
    -176.615_029_162_140_6,
    12.507_343_278_686_905,
    -0.138_571_095_265_720_12,
    9.984_369_578_019_572e-6,
    1.505_632_735_149_311_6e-7,
];

/// Lanczos g parameter.
const LANCZOS_G: f64 = 7.0;

/// Tolerance for numerical comparisons.
const TOLERANCE: f64 = 1e-20;

/// Minimum number of observed distinct molecules for analysis.
const MIN_DISTINCT: u64 = 10;

// ============================================================================
// PreseqAccum — fragment counting accumulator
// ============================================================================

/// Accumulator that counts fragment occurrences for library complexity estimation.
///
/// Counts how many times each distinct fragment is observed in the BAM file.
/// Matching preseq behavior: only unmapped reads (tid == -1) are skipped.
/// Secondary, supplementary, and QC-fail reads are all counted.
///
/// Uses a hash table to count ALL occurrences of each fragment globally. This
/// correctly groups duplicates regardless of their order in the BAM file.
/// Multi-thread safe (accumulators can be merged).
///
/// For paired-end data, a fragment is identified by `(tid, pos, mtid, mpos)` —
/// no strand, no normalization. Both reads of a pair produce separate entries.
/// For single-end data, the key is `(tid, pos)`.
#[derive(Debug)]
pub struct PreseqAccum {
    /// Maps fragment key (hash) to observation count.
    fragment_counts: HashMap<u64, u32>,
    /// Total number of fragments counted.
    pub total_fragments: u64,
}

impl PreseqAccum {
    /// Create a new accumulator.
    pub fn new() -> Self {
        Self {
            fragment_counts: HashMap::new(),
            total_fragments: 0,
        }
    }

    /// Process a fragment for library complexity estimation.
    ///
    /// Callers are responsible for fragment-level counting:
    /// - **PE data**: Count read1 of primary proper pairs + unpaired primary reads.
    ///   Compute fragment boundaries: `frag_start = pos`, `frag_end = pos + tlen`.
    ///   This matches preseq v3.2.0's `merge_mates()` behavior which creates
    ///   fragments keyed by `(chrom, fragment_start, fragment_end)`.
    /// - **SE data**: Count primary mapped reads with `frag_start = pos`, `frag_end = 0`.
    ///
    /// Each fragment is hashed and stored in a hash table. This correctly groups
    /// all duplicates regardless of BAM order.
    ///
    /// Only unmapped reads (tid < 0) are skipped internally.
    ///
    /// # Arguments
    /// * `tid` - Target (chromosome) ID.
    /// * `frag_start` - Fragment start position (leftmost alignment position).
    /// * `frag_end` - Fragment end position (for PE: pos + tlen; for SE: 0).
    pub fn process_read(&mut self, tid: i32, frag_start: i64, frag_end: i64) {
        // Only skip unmapped reads (tid == -1), matching preseq behavior
        if tid < 0 {
            return;
        }

        let key = compute_hash(tid as u64, frag_start as u64, frag_end as u64, 0);
        *self.fragment_counts.entry(key).or_insert(0) += 1;
        self.total_fragments += 1;
    }

    /// Merge another accumulator into this one by combining hash tables.
    pub fn merge(&mut self, other: PreseqAccum) {
        for (key, count) in other.fragment_counts {
            *self.fragment_counts.entry(key).or_insert(0) += count;
        }
        self.total_fragments += other.total_fragments;
    }

    /// Convert fragment counts into a frequency-of-frequencies histogram.
    ///
    /// Returns a sorted `Vec<(u64, u64)>` where each entry is `(j, n_j)`:
    /// `n_j` distinct molecules were observed exactly `j` times.
    pub fn into_histogram(self) -> Vec<(u64, u64)> {
        let mut freq_of_freq: HashMap<u64, u64> = HashMap::new();
        for &count in self.fragment_counts.values() {
            *freq_of_freq.entry(count as u64).or_insert(0) += 1;
        }
        let mut hist: Vec<(u64, u64)> = freq_of_freq.into_iter().collect();
        hist.sort_by_key(|&(j, _)| j);
        hist
    }

    /// Get the number of distinct fragments observed.
    pub fn n_distinct(&self) -> u64 {
        self.fragment_counts.len() as u64
    }
}

/// Compute a hash for a fragment key.
///
/// Uses FNV-1a hashing for speed with small keys.
fn compute_hash(a: u64, b: u64, c: u64, d: u64) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325; // FNV offset basis
    let prime: u64 = 0x100000001b3; // FNV prime
    for val in [a, b, c, d] {
        for byte in val.to_le_bytes() {
            h ^= byte as u64;
            h = h.wrapping_mul(prime);
        }
    }
    h
}

// ============================================================================
// Mathematical helpers
// ============================================================================

/// Natural log of the gamma function via Lanczos approximation.
///
/// Uses the Lanczos approximation with g=7, n=9 coefficients.
/// Accurate to ~15 significant digits for positive arguments.
fn ln_gamma(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::INFINITY;
    }
    let x = x - 1.0;
    let mut sum = LANCZOS_COEFFS[0];
    for (i, &coeff) in LANCZOS_COEFFS.iter().enumerate().skip(1) {
        sum += coeff / (x + i as f64);
    }
    let t = x + LANCZOS_G + 0.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (t.ln() * (x + 0.5)) - t + sum.ln()
}

/// Log binomial coefficient: ln(C(n, k)).
///
/// Uses the identity ln(C(n,k)) = ln_gamma(n+1) - ln_gamma(k+1) - ln_gamma(n-k+1).
fn ln_binomial(n: u64, k: u64) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    if k == 0 || k == n {
        return 0.0;
    }
    ln_gamma(n as f64 + 1.0) - ln_gamma(k as f64 + 1.0) - ln_gamma((n - k) as f64 + 1.0)
}

// ============================================================================
// Interpolation — Heck (1975) formula
// ============================================================================

/// Estimate expected distinct molecules at a target sample size ≤ N (interpolation).
///
/// Uses the Heck (1975) formula:
/// E[D(n)] = S - Σ_{j=1}^{max_freq} n_j * exp(ln_binom(N-j, n) - ln_binom(N, n))
///
/// # Arguments
/// * `histogram` - Frequency-of-frequencies: `(j, n_j)` pairs.
/// * `total_reads` - Total number of reads observed (N).
/// * `n_distinct` - Total distinct molecules observed (S).
/// * `target` - Target sample size (must be ≤ total_reads).
///
/// # Returns
/// Expected number of distinct molecules at the target sample size.
fn interpolate(histogram: &[(u64, u64)], total_reads: u64, n_distinct: u64, target: f64) -> f64 {
    let n = total_reads;
    let t = target as u64;
    if t >= n {
        return n_distinct as f64;
    }
    let mut expected = n_distinct as f64;
    for &(j, n_j) in histogram {
        if j > n {
            continue;
        }
        let log_prob = ln_binomial(n - j, t) - ln_binomial(n, t);
        expected -= n_j as f64 * log_prob.exp();
    }
    expected
}

// ============================================================================
// Extrapolation — Good-Toulmin power series + continued fraction
// ============================================================================

/// Compute the power series coefficients from the histogram.
///
/// For the Good-Toulmin estimator, the coefficients are:
/// ps_coeffs[j] = (-1)^(j+1) * n_{j+1}
///
/// # Arguments
/// * `histogram` - Frequency-of-frequencies: `(j, n_j)` pairs.
/// * `max_terms` - Maximum number of terms to compute.
///
/// # Returns
/// Vector of power series coefficients.
fn power_series_coeffs(histogram: &[(u64, u64)], max_terms: usize) -> Vec<f64> {
    // Build a lookup from frequency j -> count n_j
    let max_j = histogram.iter().map(|&(j, _)| j).max().unwrap_or(0) as usize;
    let mut nj = vec![0u64; max_j + 2];
    for &(j, n_j) in histogram {
        if (j as usize) < nj.len() {
            nj[j as usize] = n_j;
        }
    }

    let n_terms = max_terms.min(max_j);
    let mut coeffs = Vec::with_capacity(n_terms);
    for j in 0..n_terms {
        let freq = j + 1; // n_{j+1}
        let n_val = if freq < nj.len() { nj[freq] } else { 0 };
        let sign = if j % 2 == 0 { 1.0 } else { -1.0 };
        coeffs.push(sign * n_val as f64);
    }
    coeffs
}

/// Quotient-Difference (QD) algorithm to convert power series to CF coefficients.
///
/// Faithfully follows the preseq C++ implementation (`quotdiff_algorithm`
/// in `continued_fraction.cpp`). Takes power series coefficients and computes
/// the corresponding continued fraction representation using 2D QD tables.
///
/// # Arguments
/// * `ps_coeffs` - Power series coefficients `[a_0, a_1, ...]`.
///
/// # Returns
/// Vector of continued fraction coefficients, or `None` if the algorithm fails.
fn qd_algorithm(ps_coeffs: &[f64]) -> Option<Vec<f64>> {
    let n = ps_coeffs.len();
    if n < 2 {
        return ps_coeffs.first().map(|&c| vec![c]);
    }

    // Matches preseq quotdiff_algorithm exactly:
    // depth = ps_coeffs.size(), tables are depth × (depth+1)
    let depth = n;
    let width = depth + 1;

    // 2D QD tables (row-major flat vectors)
    let mut q_table = vec![0.0f64; depth * width];
    let mut e_table = vec![0.0f64; depth * width];

    let idx = |row: usize, col: usize| -> usize { row * width + col };

    // Initialize first row of quotients: q[1][j] = ps[j+1]/ps[j]
    for j in 0..depth.saturating_sub(1) {
        if ps_coeffs[j].abs() < TOLERANCE {
            return None;
        }
        if j + 1 < n {
            q_table[idx(1, j)] = ps_coeffs[j + 1] / ps_coeffs[j];
        }
    }
    // e_table[0][j] = 0 (already initialized)

    // First row of e-values: e[1][j] = q[1][j+1] - q[1][j] + e[0][j+1]
    for j in 0..depth.saturating_sub(1) {
        e_table[idx(1, j)] = q_table[idx(1, j + 1)] - q_table[idx(1, j)] + e_table[idx(0, j + 1)];
    }

    // Build remaining QD rows (matching preseq: i = 2..depth-1)
    for i in 2..depth {
        // q[i][j] = q[i-1][j+1] * e[i-1][j+1] / e[i-1][j]
        for j in 0..depth {
            let denom = e_table[idx(i - 1, j)];
            if denom.abs() < TOLERANCE {
                break;
            }
            q_table[idx(i, j)] = q_table[idx(i - 1, j + 1)] * e_table[idx(i - 1, j + 1)] / denom;
        }
        // e[i][j] = q[i][j+1] - q[i][j] + e[i-1][j+1]
        for j in 0..depth {
            e_table[idx(i, j)] =
                q_table[idx(i, j + 1)] - q_table[idx(i, j)] + e_table[idx(i - 1, j + 1)];
        }
    }

    // Extract CF coefficients (matching preseq convention):
    // cf has exactly `depth` entries.
    // cf[0] = ps_coeffs[0]
    // cf[i odd]  = -q_table[(i+1)/2][0]
    // cf[i even] = -e_table[i/2][0]
    let mut cf_coeffs = Vec::with_capacity(depth);
    cf_coeffs.push(ps_coeffs[0]);
    for i in 1..depth {
        if i % 2 != 0 {
            cf_coeffs.push(-q_table[idx(i.div_ceil(2), 0)]);
        } else {
            cf_coeffs.push(-e_table[idx(i / 2, 0)]);
        }
    }

    Some(cf_coeffs)
}

/// Evaluate a continued fraction using Euler's forward recurrence with rescaling.
///
/// Evaluates the CF representation as used by preseq:
/// `cf[0] / (1 + cf[1]*t / (1 + cf[2]*t / (1 + ...)))`
///
/// Uses the forward recurrence (Wallis recursion) which is numerically stable
/// with periodic rescaling to avoid overflow/underflow.
///
/// # Arguments
/// * `cf_coeffs` - Continued fraction coefficients from the QD algorithm.
/// * `t` - The evaluation point.
/// * `degree` - Number of terms to use.
///
/// # Returns
/// The value of the continued fraction, or None if evaluation fails.
fn evaluate_cf(cf_coeffs: &[f64], t: f64, degree: usize) -> Option<f64> {
    let n = degree.min(cf_coeffs.len());
    if n == 0 {
        return None;
    }

    // Euler's forward recurrence (matches preseq evaluate_on_diagonal):
    // Euler's forward recurrence, matching preseq's evaluate_on_diagonal:
    //   prev_num = 0, curr_num = cf[0]
    //   prev_den = 1, curr_den = 1
    //   curr_num_i = curr_num_{i-1} + cf[i]*val * prev_num_{i-1}
    //   curr_den_i = curr_den_{i-1} + cf[i]*val * prev_den_{i-1}
    //   result = curr_num / curr_den
    let mut prev_num = 0.0_f64;
    let mut curr_num = cf_coeffs[0];
    let mut prev_den = 1.0_f64;
    let mut curr_den = 1.0_f64;

    let limit = n.min(cf_coeffs.len());
    for &cf_coeff in cf_coeffs.iter().take(limit).skip(1) {
        let cf_val = cf_coeff * t;
        let new_num = curr_num + cf_val * prev_num;
        let new_den = curr_den + cf_val * prev_den;

        // Shift
        prev_num = curr_num;
        curr_num = new_num;
        prev_den = curr_den;
        curr_den = new_den;

        // Rescale to avoid overflow/underflow (matching preseq's threshold logic)
        let rescale_factor = curr_num.abs() + curr_den.abs();
        if rescale_factor >= 1.0 / TOLERANCE {
            curr_num /= rescale_factor;
            curr_den /= rescale_factor;
            prev_num /= rescale_factor;
            prev_den /= rescale_factor;
        }
        if rescale_factor <= TOLERANCE {
            curr_num /= rescale_factor;
            curr_den /= rescale_factor;
            prev_num /= rescale_factor;
            prev_den /= rescale_factor;
        }

        if !curr_num.is_finite() || !curr_den.is_finite() {
            return None;
        }
    }

    if curr_den.abs() < TOLERANCE {
        return None;
    }

    let val = curr_num / curr_den;
    if val.is_finite() {
        Some(val)
    } else {
        None
    }
}

/// Select the optimal degree for the continued fraction approximation.
///
/// Matches preseq's `optimal_cont_frac_distinct` in `continued_fraction.cpp`.
/// Tries degrees starting at `7 + (max_terms % 2 == 0)`, stepping by 2, up to
/// `max_terms`. For each candidate degree, evaluates `t * CF(t)` on a fine grid
/// (step=0.05 up to search_max=100) and checks for stability: non-negative,
/// finite, monotonically increasing, and concave (negative second derivative).
///
/// Falls back to trying `max_terms` if it's in [3, 6].
///
/// # Arguments
/// * `cf_coeffs` - Continued fraction coefficients.
/// * `total_reads` - Total observed reads (N).
/// * `max_terms` - Maximum number of power series terms (determines search range).
/// * `max_extrap` - Maximum extrapolation point.
///
/// # Returns
/// The optimal degree, or None if no suitable degree is found.
fn select_degree(
    cf_coeffs: &[f64],
    total_reads: u64,
    max_terms: usize,
    max_extrap: f64,
) -> Option<usize> {
    let n = total_reads as f64;

    // preseq uses search_max = min(max_extrap/n, 100), step = 0.05
    let search_max = (max_extrap / n).min(100.0);
    let search_step = 0.05f64;
    let n_test = (search_max / search_step).ceil() as usize;

    let check_degree = |degree: usize| -> bool {
        if degree > cf_coeffs.len() || degree == 0 {
            return false;
        }

        // Evaluate t * CF(t) on the grid (preseq's extrapolate_distinct)
        let mut prev_val = 0.0f64;
        let mut prev_deriv = 0.0f64;

        for i in 1..=n_test {
            let t = search_step * i as f64;

            let cf_val = match evaluate_cf(cf_coeffs, t, degree) {
                Some(v) => v,
                None => return false,
            };

            let val = t * cf_val;

            // Check non-negative and finite
            if val < 0.0 || !val.is_finite() {
                return false;
            }

            // Check monotonically increasing
            if val < prev_val {
                return false;
            }

            // Check concavity (second derivative <= 0)
            if i >= 2 {
                let curr_deriv = val - prev_val;
                if curr_deriv > prev_deriv {
                    return false;
                }
                prev_deriv = curr_deriv;
            } else {
                prev_deriv = val - prev_val;
            }

            prev_val = val;
        }

        true
    };

    // For small max_terms (3-6), just check if the full CF works
    if (3..=6).contains(&max_terms) {
        if check_degree(max_terms) {
            debug!(
                "Selected CF degree {} for extrapolation (small max_terms)",
                max_terms
            );
            return Some(max_terms);
        }
        return None;
    }

    // preseq: for i = 7 + (max_terms % 2 == 0); i <= max_terms; i += 2
    let start = 7 + if max_terms.is_multiple_of(2) { 1 } else { 0 };
    let mut degree = start;
    while degree <= max_terms {
        if check_degree(degree) {
            debug!("Selected CF degree {} for extrapolation", degree);
            return Some(degree);
        }
        degree += 2;
    }

    None
}

/// Extrapolate the expected number of distinct molecules beyond the observed data.
///
/// Uses the Good-Toulmin continued fraction approximation. The formula is:
/// `E[D(N+t)] = S + fold * CF(fold)`
/// where `fold = (target - N) / N` and CF is the continued fraction evaluation.
///
/// This matches preseq's `extrapolate_distinct` which computes `t * CF(t)`.
///
/// # Arguments
/// * `cf_coeffs` - Continued fraction coefficients from the QD algorithm.
/// * `degree` - Number of CF terms to use.
/// * `total_reads` - Total observed reads (N).
/// * `n_distinct` - Total observed distinct molecules (S).
/// * `target` - Target total reads (must be > total_reads).
///
/// # Returns
/// Expected number of distinct molecules at the target sample size.
fn extrapolate(
    cf_coeffs: &[f64],
    degree: usize,
    total_reads: u64,
    n_distinct: u64,
    target: f64,
) -> Option<f64> {
    let n = total_reads as f64;
    let s = n_distinct as f64;

    if target <= n {
        return Some(s);
    }

    let fold = (target - n) / n;
    let cf_val = evaluate_cf(cf_coeffs, fold, degree)?;
    // preseq formula: initial_distinct + fold * CF(fold)
    let result = s + fold * cf_val;

    if result.is_finite() && result >= 0.0 {
        Some(result)
    } else {
        None
    }
}

// ============================================================================
// Extrapolation with defects model
// ============================================================================

/// Compute power series coefficients using the "defects" formulation.
///
/// Instead of using n_j directly, this adjusts for the expected number of
/// unseen species, providing better extrapolation when the library has
/// low complexity or high duplication.
///
/// # Arguments
/// * `histogram` - Frequency-of-frequencies.
/// * `total_reads` - Total observed reads.
/// * `max_terms` - Maximum terms to compute.
///
/// # Returns
/// Power series coefficients for the defects model.
fn power_series_coeffs_defects(
    histogram: &[(u64, u64)],
    total_reads: u64,
    max_terms: usize,
) -> Vec<f64> {
    let n = total_reads as f64;

    // Build lookup
    let max_j = histogram.iter().map(|&(j, _)| j).max().unwrap_or(0) as usize;
    let mut nj = vec![0u64; max_j + 2];
    for &(j, n_j) in histogram {
        if (j as usize) < nj.len() {
            nj[j as usize] = n_j;
        }
    }

    let n_terms = max_terms.min(max_j);
    let mut coeffs = Vec::with_capacity(n_terms);

    for j in 0..n_terms {
        let freq = j + 1;
        let n_val = if freq < nj.len() {
            nj[freq] as f64
        } else {
            0.0
        };
        let sign = if j % 2 == 0 { 1.0 } else { -1.0 };
        // Adjust coefficient by (j+1)/N factor for defects model
        let adjusted = sign * n_val * ((freq + 1) as f64) / n;
        coeffs.push(adjusted);
    }
    coeffs
}

// ============================================================================
// Bootstrap confidence intervals
// ============================================================================

/// Resample a histogram using categorical (multinomial) sampling to create a bootstrap
/// replicate, matching preseq's `resample_hist` implementation.
///
/// preseq builds a "distinct-counts histogram" — the histogram of multiplicities
/// of the original histogram — and uses `std::discrete_distribution` (categorical
/// sampling) to draw `n_distinct` samples. Each draw picks a frequency bin
/// weighted by how many distinct molecules have that frequency, then the
/// selected bin's count is incremented.
///
/// # Arguments
/// * `histogram` - Original frequency-of-frequencies.
/// * `rng` - Random number generator.
///
/// # Returns
/// A new histogram from the bootstrap sample, along with the total reads and
/// number of distinct molecules in the resample.
fn bootstrap_resample(
    histogram: &[(u64, u64)],
    rng: &mut ChaCha8Rng,
) -> (Vec<(u64, u64)>, u64, u64) {
    use rand::distr::weighted::WeightedIndex;

    // Build parallel vectors of histogram indices and their counts (n_j values).
    // These are the indices into the original histogram that have nonzero n_j.
    // In preseq this is "orig_hist_distinct_counts" (the indices) and
    // "distinct_orig_hist" (the positive n_j values).
    let mut hist_indices: Vec<u64> = Vec::new();
    let mut hist_weights: Vec<u64> = Vec::new();
    for &(j, n_j) in histogram {
        if n_j > 0 {
            hist_indices.push(j);
            hist_weights.push(n_j);
        }
    }

    // Total distinct molecules = sum of all n_j values
    let n_distinct: u64 = hist_weights.iter().sum();

    if n_distinct == 0 || hist_weights.is_empty() {
        return (Vec::new(), 0, 0);
    }

    // Categorical sampling: draw n_distinct samples weighted by hist_weights.
    // Each draw selects a frequency bin; we increment that bin's count in the
    // resampled histogram.
    let weights_f64: Vec<f64> = hist_weights.iter().map(|&w| w as f64).collect();
    let dist = match WeightedIndex::new(&weights_f64) {
        Ok(d) => d,
        Err(_) => return (Vec::new(), 0, 0),
    };

    // sample[k] = how many times frequency bin k was drawn
    let mut sample = vec![0u64; hist_indices.len()];
    for _ in 0..n_distinct {
        let idx = dist.sample(rng);
        sample[idx] += 1;
    }

    // Reconstruct the histogram from the resampled data.
    // For each frequency bin j (from hist_indices), if sample[k] > 0, then
    // sample[k] distinct molecules have frequency j in the bootstrap sample.
    let mut freq_of_freq: HashMap<u64, u64> = HashMap::new();
    let mut resample_total = 0u64;
    let mut resample_distinct = 0u64;
    for (k, &count) in sample.iter().enumerate() {
        if count > 0 {
            let freq = hist_indices[k];
            *freq_of_freq.entry(freq).or_insert(0) += count;
            resample_total += freq * count;
            resample_distinct += count;
        }
    }

    let mut hist: Vec<(u64, u64)> = freq_of_freq.into_iter().collect();
    hist.sort_by_key(|&(j, _)| j);
    (hist, resample_total, resample_distinct)
}

// ============================================================================
// Main estimation pipeline
// ============================================================================

/// Result of the library complexity estimation.
#[derive(Debug)]
pub struct PreseqResult {
    /// The complexity curve: (total_reads, expected_distinct, lower_ci, upper_ci).
    pub curve: Vec<(f64, f64, f64, f64)>,
}

/// Check the Good-Toulmin 2x extrapolation condition.
///
/// In preseq, this checks that `sum((-1)^(i+1) * hist[i])` is non-negative.
/// If negative, extrapolation is unreliable.
fn good_toulmin_2x_extrap(histogram: &[(u64, u64)]) -> f64 {
    let mut result = 0.0f64;
    for &(j, n_j) in histogram {
        if j == 0 {
            continue;
        }
        let sign = if j % 2 == 1 { 1.0 } else { -1.0 };
        result += sign * n_j as f64;
    }
    result
}

/// Run the full library complexity estimation pipeline.
///
/// Matches preseq's `lc_extrap` behavior: bootstrap replicates are used to
/// compute the median as the point estimate (not a separate non-bootstrap curve),
/// with quantile confidence intervals. Uses the original histogram's
/// `initial_distinct` for extrapolation, but each bootstrap's `vals_sum` for
/// the fold calculation.
///
/// # Arguments
/// * `histogram` - Frequency-of-frequencies from the accumulator.
/// * `total_reads` - Total number of reads/fragments observed.
/// * `n_distinct` - Total number of distinct fragments observed.
/// * `config` - Preseq configuration parameters.
///
/// # Returns
/// A `PreseqResult` containing the complexity curve with confidence intervals.
pub fn estimate_complexity(
    histogram: &[(u64, u64)],
    total_reads: u64,
    n_distinct: u64,
    config: &PreseqConfig,
) -> Result<PreseqResult> {
    if n_distinct < MIN_DISTINCT {
        bail!(
            "Too few distinct molecules ({}) for library complexity estimation (minimum: {})",
            n_distinct,
            MIN_DISTINCT
        );
    }

    info!(
        "Library complexity estimation: {} total reads, {} distinct",
        total_reads, n_distinct
    );

    // Good-Toulmin 2x extrapolation check
    let gt2x = good_toulmin_2x_extrap(histogram);
    if gt2x < 0.0 {
        bail!(
            "Good-Toulmin 2x extrapolation check failed (value={:.2}). \
             Library complexity estimation is unreliable for this sample.",
            gt2x
        );
    }

    let max_extrap = config.max_extrap;
    let step = config.step_size;

    // Build the evaluation grid: starts at step_size, goes to max_extrap
    let mut targets: Vec<f64> = Vec::new();
    let mut t = step;
    while t <= max_extrap {
        targets.push(t);
        t += step;
    }

    let n_bootstraps = config.n_bootstraps;

    if n_bootstraps == 0 {
        // No bootstrapping — return a single point estimate curve
        let point_curve = compute_curve(
            histogram,
            total_reads,
            n_distinct,
            &targets,
            config.max_terms,
            max_extrap,
            config.defects,
        )?;
        let curve: Vec<(f64, f64, f64, f64)> = point_curve
            .iter()
            .map(|&(t, e)| (t, e, f64::NAN, f64::NAN))
            .collect();
        return Ok(PreseqResult { curve });
    }

    // --- Bootstrap mode (matching preseq's extrap_bootstrap) ---
    // The point estimate is the MEDIAN of bootstrap replicates, NOT a separate
    // curve from the original histogram.
    // Uses original n_distinct (initial_distinct) for extrapolation, but each
    // bootstrap's vals_sum for the fold calculation.

    info!("Running {} bootstrap replicates...", n_bootstraps);

    // Generate deterministic per-replicate seeds from the master seed.
    // We generate 2x n_bootstraps seeds initially (enough for most cases
    // where failures are rare), and fall back to a sequential retry loop
    // if more are needed.
    let boot_max_terms = config.max_terms;
    let boot_defects = config.defects;
    let n_targets = targets.len();

    // Phase 1: generate seeds and run replicates in parallel.
    // Use 2x over-provisioning to handle occasional failures.
    let initial_batch = (2 * n_bootstraps as usize).max(n_bootstraps as usize + 10);
    let replicate_seeds: Vec<u64> = {
        let mut master_rng = ChaCha8Rng::seed_from_u64(config.seed);
        (0..initial_batch).map(|_| master_rng.next_u64()).collect()
    };

    let all_results: Vec<Option<Vec<f64>>> = replicate_seeds
        .par_iter()
        .map(|&seed| {
            let mut rng = ChaCha8Rng::seed_from_u64(seed);
            let (boot_hist, boot_total, _boot_distinct) = bootstrap_resample(histogram, &mut rng);

            if boot_total == 0 {
                return None;
            }

            // Use ORIGINAL n_distinct for extrapolation, but boot_total as vals_sum
            // This matches preseq's extrap_bootstrap behavior
            match compute_curve(
                &boot_hist,
                boot_total,
                n_distinct, // original initial_distinct
                &targets,
                boot_max_terms,
                max_extrap,
                boot_defects,
            ) {
                Ok(boot_curve) => {
                    let stable = boot_curve.iter().all(|&(_, e)| e.is_finite() && e >= 0.0);
                    if stable {
                        let vals: Vec<f64> =
                            boot_curve.iter().map(|&(_, expected)| expected).collect();
                        Some(vals)
                    } else {
                        None
                    }
                }
                Err(_) => None,
            }
        })
        .collect();

    // Collect successful results, stopping after n_bootstraps successes
    let mut bootstrap_curves: Vec<Vec<f64>> = vec![Vec::new(); n_targets];
    let mut successful_bootstraps = 0u32;

    for vals in all_results.into_iter().flatten() {
        for (i, &expected) in vals.iter().enumerate() {
            if i < bootstrap_curves.len() {
                bootstrap_curves[i].push(expected);
            }
        }
        successful_bootstraps += 1;
        if successful_bootstraps >= n_bootstraps {
            break;
        }
    }

    // Phase 2: if we still need more replicates, run sequentially with
    // fresh seeds (rare — only when failure rate > 50%).
    if successful_bootstraps < n_bootstraps {
        let max_extra = 100 * n_bootstraps as u64;
        let mut extra_rng = ChaCha8Rng::seed_from_u64(config.seed.wrapping_add(1));
        let mut extra_iter = 0u64;
        while successful_bootstraps < n_bootstraps && extra_iter < max_extra {
            extra_iter += 1;
            let seed = extra_rng.next_u64();
            let mut rng = ChaCha8Rng::seed_from_u64(seed);
            let (boot_hist, boot_total, _) = bootstrap_resample(histogram, &mut rng);
            if boot_total == 0 {
                continue;
            }
            if let Ok(boot_curve) = compute_curve(
                &boot_hist,
                boot_total,
                n_distinct,
                &targets,
                boot_max_terms,
                max_extrap,
                boot_defects,
            ) {
                let stable = boot_curve.iter().all(|&(_, e)| e.is_finite() && e >= 0.0);
                if stable {
                    for (i, &(_, expected)) in boot_curve.iter().enumerate() {
                        if i < bootstrap_curves.len() {
                            bootstrap_curves[i].push(expected);
                        }
                    }
                    successful_bootstraps += 1;
                }
            }
        }
    }

    // Compute median and quantile confidence intervals (matching preseq's
    // vector_median_and_ci)
    let alpha = 1.0 - config.confidence_level;
    let lower_q = alpha / 2.0;
    let upper_q = 1.0 - alpha / 2.0;

    let curve: Vec<(f64, f64, f64, f64)> = targets
        .iter()
        .enumerate()
        .map(|(i, &target)| {
            if i < bootstrap_curves.len() && bootstrap_curves[i].len() >= 2 {
                let mut vals = bootstrap_curves[i].clone();
                vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                let median = median_sorted(&vals);
                let lower = quantile(&vals, lower_q);
                let upper = quantile(&vals, upper_q);
                (target, median, lower, upper)
            } else if i < bootstrap_curves.len() && bootstrap_curves[i].len() == 1 {
                let val = bootstrap_curves[i][0];
                (target, val, f64::NAN, f64::NAN)
            } else {
                (target, f64::NAN, f64::NAN, f64::NAN)
            }
        })
        .collect();

    Ok(PreseqResult { curve })
}

/// Compute the median of a sorted slice.
fn median_sorted(sorted: &[f64]) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return f64::NAN;
    }
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    }
}

/// Compute a single complexity curve (point estimates only).
#[allow(clippy::too_many_arguments)]
fn compute_curve(
    histogram: &[(u64, u64)],
    total_reads: u64,
    n_distinct: u64,
    targets: &[f64],
    max_terms: usize,
    max_extrap: f64,
    use_defects: bool,
) -> Result<Vec<(f64, f64)>> {
    let n = total_reads as f64;

    // Cap max_terms at first_zero - 1 (matching preseq's load_histogram behavior).
    // first_zero is the first multiplicity index j where histogram has no entry.
    // Then make it even so the QD algorithm produces balanced CF approximants.
    let max_j = histogram.iter().map(|&(j, _)| j).max().unwrap_or(0) as usize;
    let first_zero = {
        let mut fz = max_j + 1; // default: no gap found
        for j in 1..=max_j {
            if !histogram.iter().any(|&(mult, _)| mult == j as u64) {
                fz = j;
                break;
            }
        }
        fz
    };
    let adjusted_max_terms = if first_zero > 1 {
        let capped = max_terms.min(first_zero - 1);
        // Make even
        capped - (capped % 2)
    } else {
        // No usable histogram entries
        0
    };
    let max_terms = adjusted_max_terms;

    // Compute PS coefficients
    let ps_coeffs = if use_defects {
        power_series_coeffs_defects(histogram, total_reads, max_terms)
    } else {
        power_series_coeffs(histogram, max_terms)
    };

    // Convert to continued fraction
    let cf_coeffs = qd_algorithm(&ps_coeffs)
        .context("Failed to compute continued fraction coefficients via QD algorithm")?;

    // Select optimal degree
    let degree =
        select_degree(&cf_coeffs, total_reads, max_terms, max_extrap).unwrap_or_else(|| {
            warn!("Could not find optimal CF degree, using max available");
            cf_coeffs.len().min(max_terms)
        });

    // Evaluate the curve
    let mut curve = Vec::with_capacity(targets.len());
    for &target in targets {
        let expected = if target <= n {
            // Interpolation region
            interpolate(histogram, total_reads, n_distinct, target)
        } else {
            // Extrapolation region
            match extrapolate(&cf_coeffs, degree, total_reads, n_distinct, target) {
                Some(v) => v,
                None => {
                    // If extrapolation fails, stop the curve here
                    break;
                }
            }
        };

        // Sanity check: expected should be non-negative
        let expected = expected.max(0.0);
        curve.push((target, expected));
    }

    Ok(curve)
}

/// Compute a quantile from a sorted slice.
fn quantile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    if sorted.len() == 1 {
        return sorted[0];
    }
    let n = sorted.len() as f64;
    let index = q * (n - 1.0);
    let lo = index.floor() as usize;
    let hi = index.ceil() as usize;
    let frac = index - lo as f64;
    if hi >= sorted.len() {
        sorted[sorted.len() - 1]
    } else {
        sorted[lo] * (1.0 - frac) + sorted[hi] * frac
    }
}

// ============================================================================
// Output
// ============================================================================

/// Write the preseq lc_extrap output in TSV format.
///
/// Output format matches preseq's lc_extrap output:
/// ```text
/// TOTAL_READS  EXPECTED_DISTINCT  LOWER_0.95CI  UPPER_0.95CI
/// 0.0  0.0  0.0  0.0
/// 1000000.0  500000.0  480000.0  520000.0
/// ...
/// ```
///
/// # Arguments
/// * `result` - The preseq estimation result.
/// * `output_path` - Path to write the output file.
/// * `confidence_level` - Confidence level used (for column header).
pub fn write_output(
    result: &PreseqResult,
    output_path: &Path,
    confidence_level: f64,
) -> Result<()> {
    let mut f = std::fs::File::create(output_path)
        .with_context(|| format!("Failed to create preseq output: {}", output_path.display()))?;

    // Write header
    writeln!(
        f,
        "TOTAL_READS\tEXPECTED_DISTINCT\tLOWER_{:.2}CI\tUPPER_{:.2}CI",
        confidence_level, confidence_level
    )?;

    // Write first row of zeros
    writeln!(f, "0.0\t0.0\t0.0\t0.0")?;

    // Write curve
    for &(total, expected, lower, upper) in &result.curve {
        if lower.is_nan() || upper.is_nan() {
            writeln!(f, "{:.1}\t{:.1}\t-\t-", total, expected)?;
        } else {
            writeln!(
                f,
                "{:.1}\t{:.1}\t{:.1}\t{:.1}",
                total, expected, lower, upper
            )?;
        }
    }

    info!("Wrote preseq output to {}", output_path.display());
    Ok(())
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ln_gamma_known_values() {
        // ln(Γ(1)) = 0
        assert!((ln_gamma(1.0) - 0.0).abs() < 1e-10);
        // ln(Γ(2)) = 0
        assert!((ln_gamma(2.0) - 0.0).abs() < 1e-10);
        // Γ(5) = 24, ln(24) ≈ 3.178054
        assert!((ln_gamma(5.0) - (24.0f64).ln()).abs() < 1e-10);
        // Γ(0.5) = √π, ln(√π) ≈ 0.572365
        let expected = std::f64::consts::PI.sqrt().ln();
        assert!(
            (ln_gamma(0.5) - expected).abs() < 1e-10,
            "ln_gamma(0.5) = {}, expected {}",
            ln_gamma(0.5),
            expected
        );
    }

    #[test]
    fn test_ln_binomial() {
        // C(10, 0) = 1
        assert!((ln_binomial(10, 0) - 0.0).abs() < 1e-10);
        // C(10, 10) = 1
        assert!((ln_binomial(10, 10) - 0.0).abs() < 1e-10);
        // C(10, 5) = 252
        assert!(
            (ln_binomial(10, 5) - (252.0f64).ln()).abs() < 1e-8,
            "ln_binomial(10,5) = {}, expected {}",
            ln_binomial(10, 5),
            (252.0f64).ln()
        );
        // C(20, 10) = 184756
        assert!(
            (ln_binomial(20, 10) - (184756.0f64).ln()).abs() < 1e-6,
            "ln_binomial(20,10) = {}, expected {}",
            ln_binomial(20, 10),
            (184756.0f64).ln()
        );
        // k > n should give -∞
        assert!(ln_binomial(5, 10).is_infinite());
    }

    #[test]
    fn test_compute_hash_deterministic() {
        let h1 = compute_hash(1, 100, 200, 0);
        let h2 = compute_hash(1, 100, 200, 0);
        assert_eq!(h1, h2);

        // Different inputs should (almost certainly) give different hashes
        let h3 = compute_hash(1, 100, 200, 1);
        assert_ne!(h1, h3);
    }

    #[test]
    fn test_histogram_construction() {
        let mut accum = PreseqAccum::new();

        // Simulate fragments: 3 seen once, 2 seen twice, 1 seen three times
        // Key is (tid, frag_start, frag_end)
        // Fragment A: 1 time
        accum.process_read(0, 100, 100);
        // Fragment B: 1 time
        accum.process_read(0, 200, 200);
        // Fragment C: 1 time
        accum.process_read(0, 300, 300);
        // Fragment D: 2 times
        accum.process_read(0, 400, 400);
        accum.process_read(0, 400, 400);
        // Fragment E: 2 times
        accum.process_read(0, 500, 500);
        accum.process_read(0, 500, 500);
        // Fragment F: 3 times
        accum.process_read(0, 600, 600);
        accum.process_read(0, 600, 600);
        accum.process_read(0, 600, 600);

        assert_eq!(accum.total_fragments, 10);
        assert_eq!(accum.n_distinct(), 6);

        let hist = accum.into_histogram();
        assert_eq!(hist.len(), 3);
        assert_eq!(hist[0], (1, 3)); // 3 fragments seen 1 time
        assert_eq!(hist[1], (2, 2)); // 2 fragments seen 2 times
        assert_eq!(hist[2], (3, 1)); // 1 fragment seen 3 times
    }

    #[test]
    fn test_accumulator_filters() {
        let mut accum = PreseqAccum::new();

        // Mapped read — should be counted
        accum.process_read(0, 100, 100);
        assert_eq!(accum.total_fragments, 1);

        // Unmapped (tid == -1) — should be skipped
        accum.process_read(-1, 200, 200);
        assert_eq!(accum.total_fragments, 1);

        // Another mapped read at different position — should be counted
        accum.process_read(0, 300, 300);
        assert_eq!(accum.total_fragments, 2);

        // Duplicate at same position — should be counted (whole point!)
        accum.process_read(0, 300, 300);
        assert_eq!(accum.total_fragments, 3);
    }

    #[test]
    fn test_pe_accumulator() {
        let mut accum = PreseqAccum::new();

        // PE fragment spanning 100-200 (callers compute frag_start/frag_end from BAM)
        accum.process_read(0, 100, 200);
        assert_eq!(accum.total_fragments, 1);

        // Same fragment again (PCR duplicate) — same key
        accum.process_read(0, 100, 200);
        assert_eq!(accum.total_fragments, 2);

        // Different fragment (same start, different end = different insert size)
        accum.process_read(0, 100, 300);
        assert_eq!(accum.total_fragments, 3);
        assert_eq!(accum.n_distinct(), 2); // Two distinct fragments
    }

    #[test]
    fn test_interleaved_fragments_grouped_correctly() {
        // Simulate fragments at the same start but different ends (different insert sizes):
        //   (0, 100, 200) — fragment A
        //   (0, 100, 300) — fragment B (different end)
        //   (0, 100, 200) — fragment A duplicate
        //
        // Hash-based counting correctly groups A (count 2), B (count 1) → 2 distinct

        let mut accum = PreseqAccum::new();
        accum.process_read(0, 100, 200); // A
        accum.process_read(0, 100, 300); // B (different end)
        accum.process_read(0, 100, 200); // A duplicate

        assert_eq!(accum.total_fragments, 3);
        assert_eq!(accum.n_distinct(), 2, "2 distinct (A and B)");

        let hist = accum.into_histogram();
        // A has count 2, B has count 1 → histogram: {1: 1, 2: 1}
        assert_eq!(hist.len(), 2);
        assert!(hist.contains(&(1, 1)), "One singleton (B)");
        assert!(hist.contains(&(2, 1)), "One duplicate pair (A)");
    }

    #[test]
    fn test_accumulator_merge() {
        let mut accum1 = PreseqAccum::new();
        accum1.process_read(0, 100, 100);
        accum1.process_read(0, 200, 200);

        let mut accum2 = PreseqAccum::new();
        accum2.process_read(0, 100, 100); // Same fragment as accum1
        accum2.process_read(0, 300, 300);

        accum1.merge(accum2);
        assert_eq!(accum1.total_fragments, 4);
        assert_eq!(accum1.n_distinct(), 3); // 100, 200, 300
    }

    #[test]
    fn test_interpolation_identity() {
        // At target = total_reads, interpolation should return n_distinct
        let hist = vec![(1, 50), (2, 30), (3, 20)];
        let total_reads = 150; // 50*1 + 30*2 + 20*3
        let n_distinct = 100; // 50 + 30 + 20
        let result = interpolate(&hist, total_reads, n_distinct, total_reads as f64);
        assert!(
            (result - n_distinct as f64).abs() < 1.0,
            "At N, interpolation should ≈ S, got {}",
            result
        );
    }

    #[test]
    fn test_interpolation_monotone() {
        // Interpolation should be monotonically increasing
        let hist = vec![(1, 100), (2, 50), (3, 25), (4, 10)];
        let total_reads = 335; // 100 + 100 + 75 + 40
        let n_distinct = 185;

        let mut prev = 0.0;
        for target in (10..=total_reads).step_by(10) {
            let val = interpolate(&hist, total_reads, n_distinct, target as f64);
            assert!(
                val >= prev - 1e-6,
                "Interpolation should be monotone: {} < {} at target={}",
                val,
                prev,
                target
            );
            prev = val;
        }
    }

    #[test]
    fn test_power_series_coeffs() {
        let hist = vec![(1, 100), (2, 50), (3, 25)];
        let coeffs = power_series_coeffs(&hist, 10);

        // coeffs[0] = (-1)^1 * n_1 = +100 (sign for j=0: (-1)^0 = +1, times n_1)
        assert!((coeffs[0] - 100.0).abs() < 1e-10);
        // coeffs[1] = (-1)^1 * n_2 = -50
        assert!((coeffs[1] - (-50.0)).abs() < 1e-10);
        // coeffs[2] = (-1)^0 * n_3 = +25
        assert!((coeffs[2] - 25.0).abs() < 1e-10);
    }

    #[test]
    fn test_qd_algorithm_produces_cf() {
        let hist = vec![(1, 1000), (2, 500), (3, 250), (4, 125), (5, 60)];
        let ps = power_series_coeffs(&hist, 10);
        let cf = qd_algorithm(&ps);
        assert!(
            cf.is_some(),
            "QD algorithm should succeed for reasonable input"
        );
        let cf = cf.unwrap();
        assert!(!cf.is_empty(), "CF coefficients should not be empty");
    }

    #[test]
    fn test_cf_evaluation() {
        // Single-term CF: evaluate_cf([c], t, 1) = c (Euler forward recurrence: h0=c, k0=1)
        // The caller multiplies by t to get the power series value: t * CF(t)
        let cf_coeffs = vec![100.0];
        let result = evaluate_cf(&cf_coeffs, 0.5, 1);
        assert!(result.is_some());
        assert!(
            (result.unwrap() - 100.0).abs() < 1e-6,
            "evaluate_cf with single coeff should give c=100, got {}",
            result.unwrap()
        );

        // With two terms [a, b]: CF = a / (1 + b*t) via Euler recurrence
        // Note: [100.0, -2.0] at t=0.5 causes k1=0 (degenerate), so we test
        // with coefficients that don't cause division by zero.
        let cf_coeffs3 = vec![100.0, 0.5];
        let result3 = evaluate_cf(&cf_coeffs3, 1.0, 2);
        assert!(result3.is_some());
        // h0=100, k0=1, h1=100+0.5*1.0*0=100, k1=1+0.5*1.0*1=1.5
        // CF = 100/1.5 = 66.667
        assert!(
            (result3.unwrap() - 100.0 / 1.5).abs() < 1e-6,
            "Two-term CF should give 100/1.5={}, got {}",
            100.0 / 1.5,
            result3.unwrap()
        );
    }

    #[test]
    fn test_quantile_computation() {
        let vals = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert!((quantile(&vals, 0.0) - 1.0).abs() < 1e-10);
        assert!((quantile(&vals, 0.5) - 3.0).abs() < 1e-10);
        assert!((quantile(&vals, 1.0) - 5.0).abs() < 1e-10);
        assert!((quantile(&vals, 0.25) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_bootstrap_reproducibility() {
        let hist = vec![(1, 100), (2, 50), (3, 25)];
        let total_reads = 225;

        let mut rng1 = ChaCha8Rng::seed_from_u64(42);
        let (h1, t1, d1) = bootstrap_resample(&hist, &mut rng1);

        let mut rng2 = ChaCha8Rng::seed_from_u64(42);
        let (h2, t2, d2) = bootstrap_resample(&hist, &mut rng2);

        assert_eq!(t1, t2, "Same seed should give same total");
        assert_eq!(d1, d2, "Same seed should give same distinct");
        assert_eq!(h1, h2, "Same seed should give same histogram");
    }

    #[test]
    fn test_curve_monotone_increasing() {
        // Create a realistic-ish histogram
        let hist: Vec<(u64, u64)> = vec![
            (1, 5000),
            (2, 2500),
            (3, 1250),
            (4, 625),
            (5, 300),
            (6, 150),
            (7, 75),
            (8, 35),
            (9, 18),
            (10, 9),
        ];
        let total_reads: u64 = hist.iter().map(|&(j, n_j)| j * n_j).sum();
        let n_distinct: u64 = hist.iter().map(|&(_, n_j)| n_j).sum();

        let config = PreseqConfig {
            enabled: true,
            max_extrap: total_reads as f64 * 10.0,
            step_size: total_reads as f64 / 10.0,
            n_bootstraps: 0, // No bootstrap for this test
            confidence_level: 0.95,
            seed: 1,
            max_terms: 100,
            defects: false,
        };

        let result = estimate_complexity(&hist, total_reads, n_distinct, &config);
        assert!(result.is_ok(), "Estimation should succeed: {:?}", result);

        let result = result.unwrap();
        let mut prev = 0.0;
        for &(_, expected, _, _) in &result.curve {
            assert!(
                expected >= prev - 1e-6,
                "Curve should be monotonically increasing: {} < {}",
                expected,
                prev
            );
            prev = expected;
        }
    }

    #[test]
    fn test_full_pipeline_with_bootstrap() {
        let hist: Vec<(u64, u64)> = vec![(1, 5000), (2, 2500), (3, 1250), (4, 625), (5, 300)];
        let total_reads: u64 = hist.iter().map(|&(j, n_j)| j * n_j).sum();
        let n_distinct: u64 = hist.iter().map(|&(_, n_j)| n_j).sum();

        let config = PreseqConfig {
            enabled: true,
            max_extrap: total_reads as f64 * 5.0,
            step_size: total_reads as f64,
            n_bootstraps: 10,
            confidence_level: 0.95,
            seed: 1,
            max_terms: 50,
            defects: false,
        };

        let result = estimate_complexity(&hist, total_reads, n_distinct, &config);
        assert!(result.is_ok(), "Full pipeline should succeed: {:?}", result);

        let result = result.unwrap();
        assert!(!result.curve.is_empty(), "Curve should not be empty");

        // Check that confidence intervals are ordered (lower ≤ upper)
        // Note: with small bootstrap counts, CIs may not bracket the point estimate
        for &(_, _expected, lower, upper) in &result.curve {
            if !lower.is_nan() && !upper.is_nan() {
                assert!(
                    lower <= upper + 1e-6,
                    "Lower CI ({}) should be ≤ upper CI ({})",
                    lower,
                    upper
                );
            }
        }
    }
}
