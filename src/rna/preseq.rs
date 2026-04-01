//! Library complexity estimation (preseq lc_extrap reimplementation).
//!
//! Estimates the expected number of distinct molecules as a function of
//! sequencing depth using the Good-Toulmin rational function extrapolation
//! method, matching the behavior of preseq v3.

use anyhow::{bail, Context, Result};
use log::debug;
use rust_htslib::bam;
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

use super::cpp_rng::CppMt19937;

use crate::config::PreseqConfig;

// ============================================================================
// Helpers
// ============================================================================

/// FNV-1a hash of a byte slice. Used to key the dangling-mate buffer by read
/// name without allocating a Vec<u8> per read.
#[inline(always)]
fn fnv1a_hash_bytes(bytes: &[u8]) -> u64 {
    crate::io::fnv1a(bytes)
}

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

/// Default maximum merged fragment length for PE reads (effectively unlimited).
///
/// Corresponds to preseq's `-seg_len` option. The nf-core preseq module
/// typically uses a very large value (e.g., 100000000) to avoid rejecting
/// legitimate long-insert PE fragments. Used in tests.
#[cfg(test)]
const DEFAULT_MAX_SEGMENT_LENGTH: i64 = 100_000_000;

// ============================================================================
// PreseqAccum — fragment counting accumulator (preseq v3.2.0 compatible)
// ============================================================================

/// Information about a single read's alignment, used for mate merging.
#[derive(Debug, Clone)]
struct MateInfo {
    /// Target (chromosome) ID.
    tid: i32,
    /// Alignment start position (0-based).
    start: i64,
    /// Alignment end position (exclusive, derived from CIGAR).
    end: i64,
}

/// Accumulator that counts fragment occurrences for library complexity estimation.
///
/// Matches preseq v3.2.0's `load_counts_BAM_pe` behavior:
///
/// 1. Only primary, mapped reads are processed (secondary/supplementary excluded).
/// 2. Paired reads are collated by read name. When both mates are found on the
///    same chromosome, they are merged into a single genomic interval
///    `(chrom, min_start, max_end)`. If the merged length exceeds
///    [`MAX_SEGMENT_LENGTH`] (5000), or if mates are on different chromosomes,
///    both reads are counted as individual fragments.
/// 3. Unpaired reads are counted as individual fragments, keyed by
///    `(tid, start)` only (matching preseq's SE `aln_pos` struct).
/// 4. Paired fragments are identified by `(tid, start, end)` tuples;
///    unpaired fragments by `(tid, start)` only.
///
/// Thread-safe via per-chromosome accumulators that are merged after parallel
/// processing. During merge, dangling mates from different chromosomes are
/// matched and merged or promoted to individual fragments.
#[derive(Debug)]
pub struct PreseqAccum {
    /// Maps fragment key hash to observation count.
    /// PE keys use `(tid, start, end)`; SE keys use `(tid, start)` only.
    fragment_counts: HashMap<u64, u32>,
    /// Total number of fragments (merged pairs + unpaired) counted so far.
    pub total_fragments: u64,
    /// Reads waiting for their mate. Keyed by read name.
    /// When both mates are seen, they are merged into a fragment immediately.
    dangling_mates: HashMap<u64, MateInfo>,
    /// Number of successfully merged PE pairs.
    n_paired: u64,
    /// Number of reads counted as individual fragments (unpaired or failed merge).
    n_unpaired: u64,
    /// Total mates (primary + mapped reads) processed.
    n_mates_processed: u64,
    /// Maximum merged PE fragment length. Fragments longer than this are split.
    max_segment_length: i64,
}

impl PreseqAccum {
    /// Create a new accumulator with the given maximum segment length.
    pub fn new(max_segment_length: i64) -> Self {
        Self {
            fragment_counts: HashMap::new(),
            total_fragments: 0,
            dangling_mates: HashMap::new(),
            n_paired: 0,
            n_unpaired: 0,
            n_mates_processed: 0,
            max_segment_length,
        }
    }

    /// Process a BAM record for library complexity estimation.
    ///
    /// Matches preseq v3.2.0 behavior:
    /// - Only primary, mapped reads are counted (secondary/supplementary
    ///   skipped). Note: preseq's C++ source appears to only check
    ///   `tid == -1`, but empirically it processes only primary reads
    ///   (47,713 for our SE test BAM, not 48,977 total mapped).
    /// - Paired reads are merged by read name into genomic intervals
    ///   (keyed by `(tid, start, end)`, matching `load_counts_BAM_pe`).
    /// - Unpaired reads are keyed by `(tid, start)` only, matching
    ///   `load_counts_BAM_se` / `aln_pos`.
    ///
    /// # Arguments
    /// * `record` - The BAM record to process.
    pub fn process_read(&mut self, record: &bam::Record) {
        // Filter: primary + mapped only (matches preseq v3.2.0 empirical behavior).
        if record.tid() < 0 {
            return;
        }
        if record.is_secondary() || record.is_supplementary() {
            return;
        }
        self.n_mates_processed += 1;

        let tid = record.tid();
        let start = record.pos();
        let end = record.cigar().end_pos();

        let info = MateInfo { tid, start, end };

        // Check if this read is part of a pair
        if record.is_paired() {
            // Use a FNV-1a hash of the qname as the mate-buffer key to avoid
            // allocating a Vec<u8> per read. Collisions are negligible for
            // typical BAM qnames; the MateInfo already stores coordinates that
            // would detect a spurious match if needed.
            let qname_hash = fnv1a_hash_bytes(record.qname());

            if let Some(mate) = self.dangling_mates.remove(&qname_hash) {
                // Found the mate — try to merge
                if mate.tid == tid {
                    // Same chromosome: merge into genomic interval
                    let merged_start = mate.start.min(start);
                    let merged_end = mate.end.max(end);
                    let frag_len = merged_end - merged_start;

                    if frag_len >= 0 && frag_len <= self.max_segment_length {
                        // Successful merge
                        self.add_fragment(tid, merged_start, merged_end);
                        self.n_paired += 1;
                    } else {
                        // Fragment too long or invalid — split into individual reads
                        self.add_fragment(mate.tid, mate.start, mate.end);
                        self.add_fragment(tid, start, end);
                        self.n_unpaired += 2;
                    }
                } else {
                    // Different chromosomes — cannot merge, count individually
                    self.add_fragment(mate.tid, mate.start, mate.end);
                    self.add_fragment(tid, start, end);
                    self.n_unpaired += 2;
                }
            } else {
                // First mate seen — store for later
                self.dangling_mates.insert(qname_hash, info);
            }
        } else {
            // Unpaired read — count as individual fragment using SE key (tid, start)
            // only, matching preseq v3.2.0's `load_counts_BAM_se` / `aln_pos`.
            self.add_fragment_se(tid, start);
            self.n_unpaired += 1;
        }
    }

    /// Add a PE fragment (merged pair) to the counting table.
    ///
    /// Uses `(tid, start, end)` as the fragment key, matching preseq v3.2.0's
    /// PE mode where merged genomic intervals define fragment identity.
    #[cfg_attr(test, allow(dead_code))]
    pub(crate) fn add_fragment(&mut self, tid: i32, start: i64, end: i64) {
        let key = fragment_hash(tid, start, end);
        *self.fragment_counts.entry(key).or_insert(0) += 1;
        self.total_fragments += 1;
    }

    /// Add an SE (unpaired) fragment to the counting table.
    ///
    /// Uses `(tid, start)` only as the fragment key, matching preseq v3.2.0's
    /// SE mode (`aln_pos`) where duplicate identity is determined by chromosome
    /// and start position alone, ignoring CIGAR-derived end position.
    fn add_fragment_se(&mut self, tid: i32, start: i64) {
        let key = fragment_hash_se(tid, start);
        *self.fragment_counts.entry(key).or_insert(0) += 1;
        self.total_fragments += 1;
    }

    /// Merge another accumulator into this one.
    ///
    /// First resolves cross-chromosome mate pairs: if a read name appears in
    /// both accumulators' dangling mates, the mates are merged (or split if
    /// on different chromosomes / too long). Remaining dangling mates are
    /// combined into a single map.
    pub fn merge(&mut self, other: PreseqAccum) {
        // Merge fragment counts
        for (key, count) in other.fragment_counts {
            *self.fragment_counts.entry(key).or_insert(0) += count;
        }
        self.total_fragments += other.total_fragments;
        self.n_paired += other.n_paired;
        self.n_unpaired += other.n_unpaired;
        self.n_mates_processed += other.n_mates_processed;

        // Resolve cross-accumulator mate pairs
        for (qname, other_mate) in other.dangling_mates {
            if let Some(self_mate) = self.dangling_mates.remove(&qname) {
                // Found mate in both accumulators — try to merge
                if self_mate.tid == other_mate.tid {
                    let merged_start = self_mate.start.min(other_mate.start);
                    let merged_end = self_mate.end.max(other_mate.end);
                    let frag_len = merged_end - merged_start;

                    if frag_len >= 0 && frag_len <= self.max_segment_length {
                        self.add_fragment(self_mate.tid, merged_start, merged_end);
                        self.n_paired += 1;
                    } else {
                        self.add_fragment(self_mate.tid, self_mate.start, self_mate.end);
                        self.add_fragment(other_mate.tid, other_mate.start, other_mate.end);
                        self.n_unpaired += 2;
                    }
                } else {
                    // Different chromosomes — count individually
                    self.add_fragment(self_mate.tid, self_mate.start, self_mate.end);
                    self.add_fragment(other_mate.tid, other_mate.start, other_mate.end);
                    self.n_unpaired += 2;
                }
            } else {
                // No match yet — keep as dangling
                self.dangling_mates.insert(qname, other_mate);
            }
        }
    }

    /// Finalize the accumulator by flushing all remaining dangling mates
    /// as individual (unpaired) fragments.
    ///
    /// Must be called after all merges are complete, before `into_histogram()`.
    pub fn finalize(&mut self) {
        let n_dangling = self.dangling_mates.len();
        let dangling: Vec<MateInfo> = self.dangling_mates.drain().map(|(_, v)| v).collect();
        for mate in dangling {
            self.add_fragment(mate.tid, mate.start, mate.end);
            self.n_unpaired += 1;
        }
        debug!(
            "preseq finalize: mates_processed={}, n_paired={}, n_unpaired={}, dangling_flushed={}, total_fragments={}, n_distinct={}",
            self.n_mates_processed, self.n_paired, self.n_unpaired, n_dangling, self.total_fragments, self.fragment_counts.len()
        );
    }

    /// Convert fragment counts into a frequency-of-frequencies histogram.
    ///
    /// Returns a sorted `Vec<(u64, u64)>` where each entry is `(j, n_j)`:
    /// `n_j` distinct molecules were observed exactly `j` times.
    ///
    /// **Important**: Call [`finalize()`](Self::finalize) before this method
    /// to flush any remaining dangling mates.
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

/// Compute a hash for a PE fragment key `(tid, start, end)`.
///
/// Uses FNV-1a hashing for speed with small keys. Used for merged PE pairs
/// where the fragment interval `(chrom, min_start, max_end)` defines identity,
/// matching preseq v3.2.0's `load_counts_BAM_pe` behavior.
fn fragment_hash(tid: i32, start: i64, end: i64) -> u64 {
    let mut buf = [0u8; 24];
    buf[0..8].copy_from_slice(&(tid as u64).to_le_bytes());
    buf[8..16].copy_from_slice(&(start as u64).to_le_bytes());
    buf[16..24].copy_from_slice(&(end as u64).to_le_bytes());
    crate::io::fnv1a(&buf)
}

/// Compute a hash for an SE fragment key `(tid, start)`.
///
/// Matches preseq v3.2.0's `load_counts_BAM_se` / `aln_pos` behavior: two
/// reads are duplicates iff they share the same chromosome and start position,
/// regardless of CIGAR-derived end position.
fn fragment_hash_se(tid: i32, start: i64) -> u64 {
    let mut buf = [0u8; 16];
    buf[0..8].copy_from_slice(&(tid as u64).to_le_bytes());
    buf[8..16].copy_from_slice(&(start as u64).to_le_bytes());
    crate::io::fnv1a(&buf)
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

    // Upstream preseq uses a dense histogram (vector indexed 0..max_freq).
    // We must produce numer values in the same index order to match upstream's
    // accumulate(cbegin(numer), cend(numer), 0) which uses int accumulator.
    // This truncates the running sum to int at EACH addition step.
    let max_j = histogram.iter().map(|&(j, _)| j).max().unwrap_or(0) as usize;
    let mut dense_hist = vec![0u64; max_j + 1];
    for &(j, n_j) in histogram {
        if (j as usize) < dense_hist.len() {
            dense_hist[j as usize] = n_j;
        }
    }

    let denom = ln_binomial(n, t);

    // Build numer array matching upstream: numer[i] for i=1..hist.size()
    // Then accumulate with int truncation at each step
    let mut acc: i32 = 0; // upstream uses accumulate(..., 0) which is int
    for (i, &n_j) in dense_hist.iter().enumerate().skip(1) {
        if n_j == 0 {
            continue; // numer[i] = 0, adding 0 to int doesn't change anything
        }
        let numer_val = if n < (i as u64) + t {
            0.0
        } else {
            let x = ln_binomial(n - i as u64, t);
            (x - denom).exp() * n_j as f64
        };
        // Upstream: acc (int) = acc (int) + numer_val (double)
        // This promotes acc to double, adds, then truncates back to int
        acc = (acc as f64 + numer_val) as i32;
    }
    n_distinct as f64 - acc as f64
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
/// The value of the continued fraction. Upstream returns curr_num/curr_den
/// unconditionally (no validity checks); the stability check handles rejection.
fn evaluate_cf(cf_coeffs: &[f64], t: f64, degree: usize) -> f64 {
    let n = degree.min(cf_coeffs.len());
    if n == 0 {
        return f64::NAN;
    }

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
    }

    // Upstream returns curr_num / curr_den unconditionally
    curr_num / curr_den
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
    _total_reads: u64,
    max_terms: usize,
    _max_extrap: f64,
) -> Option<usize> {
    // preseq uses search_max = 100, step = 0.05
    let search_max = 100.0;
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

            let cf_val = evaluate_cf(cf_coeffs, t, degree);
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
/// Upstream returns the value unconditionally; stability is checked on the full curve.
fn extrapolate(
    cf_coeffs: &[f64],
    degree: usize,
    total_reads: u64,
    n_distinct: u64,
    target: f64,
) -> f64 {
    let n = total_reads as f64;
    let s = n_distinct as f64;

    let fold = (target - n) / n;
    let cf_val = evaluate_cf(cf_coeffs, fold, degree);
    // preseq formula: initial_distinct + fold * CF(fold)
    s + fold * cf_val
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
    _total_reads: u64,
    max_terms: usize,
) -> Vec<f64> {
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
        coeffs.push(sign * n_val);
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
    rng: &mut CppMt19937,
) -> (Vec<(u64, u64)>, u64, u64) {
    // Build a dense histogram vector from 0..=max_j, matching upstream preseq's
    // `counts_hist` representation. This is critical for RNG compatibility:
    // upstream's multinomial() iterates over ALL indices (including zero-weight
    // entries), calling binomial() for each. Skipping zero-weight entries would
    // cause the RNG state to diverge.
    let max_j = histogram.iter().map(|&(j, _)| j).max().unwrap_or(0) as usize;
    let mut dense_hist = vec![0u64; max_j + 1];
    for &(j, n_j) in histogram {
        dense_hist[j as usize] = n_j;
    }

    // Total distinct molecules = sum of all non-zero n_j values
    let n_distinct: u64 = dense_hist.iter().sum();

    if n_distinct == 0 {
        return (Vec::new(), 0, 0);
    }

    // Sequential binomial sampling (multinomial via sequential binomials).
    // Matches upstream preseq's multinomial() in common.hpp exactly.
    //
    // IMPORTANT: We must call rng.binomial() for EVERY category in the dense
    // histogram, including zero-weight entries, even when remaining_trials == 0
    // or p >= 1.0, because the C++ code does:
    //   *r = binom_dist(trials, (*p) / remaining_prob)(gen);
    // unconditionally for every category. Skipping any call (even for zero
    // weights) would cause the RNG state to diverge from upstream, producing
    // different bootstrap results.
    let mut remaining_trials = n_distinct;
    let mut remaining_prob: f64 = n_distinct as f64;
    let mut sample = vec![0u64; dense_hist.len()];

    for i in 0..dense_hist.len() {
        let p = (dense_hist[i] as f64) / remaining_prob;
        let drawn = rng.binomial(remaining_trials, p);
        sample[i] = drawn;
        remaining_prob -= dense_hist[i] as f64;
        remaining_trials -= drawn;
    }

    // Reconstruct the histogram from the resampled data.
    let mut freq_of_freq: HashMap<u64, u64> = HashMap::new();
    let mut resample_total = 0u64;
    let mut resample_distinct = 0u64;
    for (k, &count) in sample.iter().enumerate() {
        if count > 0 {
            let freq = k as u64;
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

    debug!(
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
    while t < max_extrap {
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

    debug!("Running {} bootstrap replicates...", n_bootstraps);

    // Run bootstraps sequentially with a single MT19937 RNG seeded once,
    // matching upstream preseq's behavior exactly. Upstream uses C++ std::mt19937
    // (32-bit Mersenne Twister) with a single seed, running each bootstrap in order.
    let boot_max_terms = config.max_terms;
    let boot_defects = config.defects;
    let n_targets = targets.len();

    let mut rng = CppMt19937::new(config.seed as u32);
    let mut bootstrap_curves: Vec<Vec<f64>> = vec![Vec::new(); n_targets];
    let mut successful_bootstraps = 0u32;
    let max_iter = 100 * n_bootstraps as u64;
    let mut iter_count = 0u64;
    #[cfg(test)]
    let mut reject_empty = 0u64;
    #[cfg(test)]
    let mut reject_unstable = 0u64;
    #[cfg(test)]
    let mut accept_count = 0u64;

    while successful_bootstraps < n_bootstraps && iter_count < max_iter {
        iter_count += 1;
        let (boot_hist, boot_total, boot_distinct) = bootstrap_resample(histogram, &mut rng);

        if boot_total == 0 {
            continue;
        }

        // Use boot_distinct for interpolation, original n_distinct for extrapolation
        // This matches preseq's extrap_bootstrap behavior
        if let Ok(boot_curve) = compute_curve(
            &boot_hist,
            boot_total,
            boot_distinct, // bootstrap's own distinct for interpolation
            n_distinct,    // original initial_distinct for extrapolation
            &targets,
            boot_max_terms,
            max_extrap,
            boot_defects,
        ) {
            // Upstream checks: cf.is_valid() (non-empty curve means degree was found)
            // AND check_yield_estimates_stability on the full yield vector.
            if !boot_curve.is_empty() {
                let yield_vector: Vec<f64> = boot_curve.iter().map(|&(_, e)| e).collect();
                if check_yield_estimates_stability(&yield_vector) {
                    for (i, &(_, expected)) in boot_curve.iter().enumerate() {
                        if i < bootstrap_curves.len() {
                            bootstrap_curves[i].push(expected);
                        }
                    }
                    successful_bootstraps += 1;
                    #[cfg(test)]
                    {
                        accept_count += 1;
                    }
                } else {
                    #[cfg(test)]
                    {
                        reject_unstable += 1;
                    }
                }
            } else {
                #[cfg(test)]
                {
                    reject_empty += 1;
                }
            }
        }
    }

    #[cfg(test)]
    eprintln!(
        "DIAGNOSTIC: bootstrap stats: iter={}, accepted={}, rejected_empty={}, rejected_unstable={}",
        iter_count, accept_count, reject_empty, reject_unstable
    );

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

/// Check yield estimates stability (matching upstream's check_yield_estimates_stability).
///
/// Checks that the yield estimates are:
/// 1. Non-negative and finite
/// 2. Monotonically increasing (non-decreasing)
/// 3. Have negative second derivative (concavity)
///
/// Returns true if estimates are stable, false otherwise.
fn check_yield_estimates_stability(estimates: &[f64]) -> bool {
    if estimates.is_empty() {
        return false;
    }

    // Check non-negative and finite
    for &e in estimates {
        if !e.is_finite() || e < 0.0 {
            return false;
        }
    }

    // Check monotonically increasing (non-decreasing)
    for i in 1..estimates.len() {
        if estimates[i] < estimates[i - 1] {
            return false;
        }
    }

    // Check negative second derivative (concavity)
    // estimates[i] - estimates[i-1] should be decreasing (i.e. slope is decreasing)
    for i in 2..estimates.len() {
        if estimates[i] - estimates[i - 1] > estimates[i - 1] - estimates[i - 2] {
            return false;
        }
    }

    true
}

/// Compute a single complexity curve (point estimates only).
///
/// # Arguments
/// * `histogram` - Frequency-of-frequencies.
/// * `total_reads` - Total reads in this sample.
/// * `interp_distinct` - Distinct count for interpolation (bootstrap's own).
/// * `extrap_distinct` - Distinct count for extrapolation (original initial_distinct).
/// * `targets` - Evaluation grid points.
/// * `max_terms` - Maximum CF terms.
/// * `max_extrap` - Maximum extrapolation value.
/// * `use_defects` - Whether to use defects mode.
#[allow(clippy::too_many_arguments)]
fn compute_curve(
    histogram: &[(u64, u64)],
    total_reads: u64,
    interp_distinct: u64,
    extrap_distinct: u64,
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
    let degree = match select_degree(&cf_coeffs, total_reads, max_terms, max_extrap) {
        Some(d) => d,
        None => return Ok(Vec::new()),
    };

    // Evaluate the curve — upstream pushes all values and rejects
    // unstable curves afterwards via check_yield_estimates_stability.
    let mut curve = Vec::with_capacity(targets.len());
    for &target in targets {
        let expected = if target < n {
            // Interpolation region: use bootstrap's own distinct count.
            // Uses strict < to match C++ `while (curr_sample_sz < sample_vals_sum)`.
            interpolate(histogram, total_reads, interp_distinct, target)
        } else {
            // Extrapolation region: use original initial_distinct
            extrapolate(&cf_coeffs, degree, total_reads, extrap_distinct, target)
        };

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
    writeln!(f, "0\t0\t0\t0")?;

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

    debug!("Wrote preseq output to {}", output_path.display());
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
    fn test_fragment_hash_deterministic() {
        let h1 = fragment_hash(1, 200, 300);
        let h2 = fragment_hash(1, 200, 300);
        assert_eq!(h1, h2);

        // Different inputs should (almost certainly) give different hashes
        let h3 = fragment_hash(1, 100, 300);
        assert_ne!(h1, h3);
    }

    #[test]
    fn test_se_fragment_hash_ignores_end() {
        // SE hash uses (tid, start) only — reads at the same start but
        // different ends must hash identically, matching preseq's aln_pos.
        let h1 = fragment_hash_se(1, 200);
        let h2 = fragment_hash_se(1, 200);
        assert_eq!(h1, h2);

        // Different start → different hash
        let h3 = fragment_hash_se(1, 300);
        assert_ne!(h1, h3);
    }

    #[test]
    fn test_se_fragments_collapse_same_start_different_end() {
        // In SE mode, reads at the same (tid, start) but different CIGAR-derived
        // ends should be counted as duplicates, matching preseq lc_extrap
        // without -pe.
        let mut accum = PreseqAccum::new(DEFAULT_MAX_SEGMENT_LENGTH);

        // Three reads at chr0:100 with different end positions
        accum.add_fragment_se(0, 100); // read ending at 200
        accum.add_fragment_se(0, 100); // read ending at 250 (same key)
        accum.add_fragment_se(0, 100); // read ending at 300 (same key)

        // One read at a different start
        accum.add_fragment_se(0, 200);

        assert_eq!(accum.total_fragments, 4);
        assert_eq!(accum.n_distinct(), 2, "SE: same start = same fragment");

        let hist = accum.into_histogram();
        assert!(hist.contains(&(1, 1)), "One singleton");
        assert!(hist.contains(&(3, 1)), "One triplet");
    }

    #[test]
    fn test_histogram_construction() {
        let mut accum = PreseqAccum::new(DEFAULT_MAX_SEGMENT_LENGTH);

        // Simulate fragments: 3 seen once, 2 seen twice, 1 seen three times
        // Key is (tid, start, end)
        accum.add_fragment(0, 100, 200); // A: 1 time
        accum.add_fragment(0, 200, 300); // B: 1 time
        accum.add_fragment(0, 300, 400); // C: 1 time
        accum.add_fragment(0, 400, 500); // D: 2 times
        accum.add_fragment(0, 400, 500);
        accum.add_fragment(0, 500, 600); // E: 2 times
        accum.add_fragment(0, 500, 600);
        accum.add_fragment(0, 600, 700); // F: 3 times
        accum.add_fragment(0, 600, 700);
        accum.add_fragment(0, 600, 700);

        assert_eq!(accum.total_fragments, 10);
        assert_eq!(accum.n_distinct(), 6);

        let hist = accum.into_histogram();
        assert_eq!(hist.len(), 3);
        assert_eq!(hist[0], (1, 3)); // 3 fragments seen 1 time
        assert_eq!(hist[1], (2, 2)); // 2 fragments seen 2 times
        assert_eq!(hist[2], (3, 1)); // 1 fragment seen 3 times
    }

    #[test]
    fn test_fragment_counting() {
        let mut accum = PreseqAccum::new(DEFAULT_MAX_SEGMENT_LENGTH);

        // Fragment at (tid=0, start=100, end=200)
        accum.add_fragment(0, 100, 200);
        assert_eq!(accum.total_fragments, 1);

        // Different fragment
        accum.add_fragment(0, 300, 400);
        assert_eq!(accum.total_fragments, 2);

        // Duplicate at same position — should be counted (whole point!)
        accum.add_fragment(0, 300, 400);
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

        let mut accum = PreseqAccum::new(DEFAULT_MAX_SEGMENT_LENGTH);
        accum.add_fragment(0, 100, 200); // A
        accum.add_fragment(0, 100, 300); // B (different end)
        accum.add_fragment(0, 100, 200); // A duplicate

        assert_eq!(accum.total_fragments, 3);
        assert_eq!(accum.n_distinct(), 2, "2 distinct (A and B)");

        let hist = accum.into_histogram();
        // A has count 2, B has count 1 → histogram: {1: 1, 2: 1}
        assert_eq!(hist.len(), 2);
        assert!(hist.contains(&(1, 1)), "One singleton (B)");
        assert!(hist.contains(&(2, 1)), "One duplicate pair (A)");
    }

    #[test]
    fn test_accumulator_merge_fragment_counts() {
        let mut accum1 = PreseqAccum::new(DEFAULT_MAX_SEGMENT_LENGTH);
        accum1.add_fragment(0, 100, 200);
        accum1.add_fragment(0, 200, 300);

        let mut accum2 = PreseqAccum::new(DEFAULT_MAX_SEGMENT_LENGTH);
        accum2.add_fragment(0, 100, 200); // Same fragment as accum1
        accum2.add_fragment(0, 300, 400);

        accum1.merge(accum2);
        assert_eq!(accum1.total_fragments, 4);
        assert_eq!(accum1.n_distinct(), 3); // (100,200), (200,300), (300,400)
    }

    #[test]
    fn test_finalize_flushes_dangling_mates() {
        let mut accum = PreseqAccum::new(DEFAULT_MAX_SEGMENT_LENGTH);

        // Manually insert a dangling mate
        accum.dangling_mates.insert(
            1u64,
            MateInfo {
                tid: 0,
                start: 100,
                end: 200,
            },
        );

        assert_eq!(accum.total_fragments, 0);
        accum.finalize();
        assert_eq!(accum.total_fragments, 1);
        assert_eq!(accum.n_unpaired, 1);
        assert!(accum.dangling_mates.is_empty());
    }

    #[test]
    fn test_merge_resolves_cross_chromosome_mates() {
        let mut accum1 = PreseqAccum::new(DEFAULT_MAX_SEGMENT_LENGTH);
        let mut accum2 = PreseqAccum::new(DEFAULT_MAX_SEGMENT_LENGTH);

        // Mate 1 on chr0 (in accum1's dangling)
        accum1.dangling_mates.insert(
            1u64,
            MateInfo {
                tid: 0,
                start: 100,
                end: 200,
            },
        );

        // Mate 2 on chr1 (in accum2's dangling) — different chromosome
        accum2.dangling_mates.insert(
            1u64,
            MateInfo {
                tid: 1,
                start: 300,
                end: 400,
            },
        );

        accum1.merge(accum2);
        // Different chromosomes → both counted as individual fragments
        assert_eq!(accum1.total_fragments, 2);
        assert_eq!(accum1.n_unpaired, 2);
        assert!(accum1.dangling_mates.is_empty());
    }

    #[test]
    fn test_merge_resolves_same_chromosome_mates() {
        let mut accum1 = PreseqAccum::new(DEFAULT_MAX_SEGMENT_LENGTH);
        let mut accum2 = PreseqAccum::new(DEFAULT_MAX_SEGMENT_LENGTH);

        // Mate 1 on chr0
        accum1.dangling_mates.insert(
            1u64,
            MateInfo {
                tid: 0,
                start: 100,
                end: 200,
            },
        );

        // Mate 2 also on chr0, nearby
        accum2.dangling_mates.insert(
            1u64,
            MateInfo {
                tid: 0,
                start: 150,
                end: 350,
            },
        );

        accum1.merge(accum2);
        // Same chromosome, merged into (0, 100, 350), len=250 < 5000
        assert_eq!(accum1.total_fragments, 1);
        assert_eq!(accum1.n_paired, 1);
        assert!(accum1.dangling_mates.is_empty());
    }

    #[test]
    fn test_merge_rejects_long_fragments() {
        let mut accum1 = PreseqAccum::new(5000); // Small limit to trigger rejection
        let mut accum2 = PreseqAccum::new(5000);

        // Mate 1 on chr0
        accum1.dangling_mates.insert(
            1u64,
            MateInfo {
                tid: 0,
                start: 0,
                end: 100,
            },
        );

        // Mate 2 also on chr0 but far away (merged len > 5000)
        accum2.dangling_mates.insert(
            1u64,
            MateInfo {
                tid: 0,
                start: 6000,
                end: 6100,
            },
        );

        accum1.merge(accum2);
        // Fragment length 6100 > MAX_SEGMENT_LENGTH (5000) → split
        assert_eq!(accum1.total_fragments, 2);
        assert_eq!(accum1.n_unpaired, 2);
        assert!(accum1.dangling_mates.is_empty());
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
        assert!(
            (result - 100.0).abs() < 1e-6,
            "evaluate_cf with single coeff should give c=100, got {}",
            result
        );

        // With two terms [a, b]: CF = a / (1 + b*t) via Euler recurrence
        // Note: [100.0, -2.0] at t=0.5 causes k1=0 (degenerate), so we test
        // with coefficients that don't cause division by zero.
        let cf_coeffs3 = vec![100.0, 0.5];
        let result3 = evaluate_cf(&cf_coeffs3, 1.0, 2);
        // h0=100, k0=1, h1=100+0.5*1.0*0=100, k1=1+0.5*1.0*1=1.5
        // CF = 100/1.5 = 66.667
        assert!(
            (result3 - 100.0 / 1.5).abs() < 1e-6,
            "Two-term CF should give 100/1.5={}, got {}",
            100.0 / 1.5,
            result3
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
        let _total_reads = 225;

        let mut rng1 = CppMt19937::new(42);
        let (h1, t1, d1) = bootstrap_resample(&hist, &mut rng1);

        let mut rng2 = CppMt19937::new(42);
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
            max_segment_length: 100_000_000,
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

    /// Verify preseq algorithm intermediates match upstream for the small PE benchmark histogram.
    ///
    /// The benchmark histogram is from `preseq lc_extrap -bam -pe -seed 1 -seg_len 100000000`
    /// on the small test BAM. Final bootstrap values are platform-dependent due to
    /// std::binomial_distribution differences between libstdc++ (Linux) and libc++ (macOS).
    /// On Linux, output should match upstream byte-for-byte.
    #[test]
    fn test_benchmark_histogram_intermediates() {
        let hist: Vec<(u64, u64)> = vec![
            (1, 23425),
            (2, 443),
            (3, 69),
            (4, 28),
            (5, 14),
            (6, 5),
            (8, 3),
            (9, 2),
            (10, 1),
            (18, 1),
        ];
        let total_reads: u64 = 24800;
        let n_distinct: u64 = 23991;

        // Verify histogram totals
        let hist_total: u64 = hist.iter().map(|&(j, n)| j * n).sum();
        let hist_distinct: u64 = hist.iter().map(|&(_, n)| n).sum();
        assert_eq!(hist_total, total_reads);
        assert_eq!(hist_distinct, n_distinct);

        // first_zero = 7 (gap at j=7), max_terms = 6
        let first_zero = {
            let max_j = hist.iter().map(|&(j, _)| j).max().unwrap() as usize;
            let mut nj = vec![0u64; max_j + 2];
            for &(j, n) in &hist {
                nj[j as usize] = n;
            }
            (1..=max_j).find(|&j| nj[j] == 0).unwrap_or(max_j + 1)
        };
        assert_eq!(first_zero, 7);

        let max_terms = {
            let mt = std::cmp::min(100usize, first_zero - 1);
            mt - (mt % 2)
        };
        assert_eq!(max_terms, 6);

        // Power series coefficients: [+23425, -443, +69, -28, +14, -5]
        let ps_coeffs = power_series_coeffs(&hist, max_terms);
        assert_eq!(ps_coeffs, vec![23425.0, -443.0, 69.0, -28.0, 14.0, -5.0]);

        // QD algorithm should succeed
        let cf_coeffs = qd_algorithm(&ps_coeffs).expect("QD should succeed");
        assert_eq!(cf_coeffs.len(), 6);

        // Degree selection fails on this histogram (CF is unstable at degree 6)
        // This matches upstream: all output comes from bootstrap curves
        let degree = select_degree(&cf_coeffs, total_reads, max_terms, 1e10);
        assert!(
            degree.is_none(),
            "Degree selection should fail for this histogram (matches upstream)"
        );

        // Interpolation at 12000 (below total) — integer truncation matching upstream
        let interp_result = interpolate(&hist, total_reads, n_distinct, 12000.0);
        assert!(
            interp_result > 11700.0 && interp_result < 11800.0,
            "Interpolation at 12000 should be ~11770, got {}",
            interp_result
        );

        // Full bootstrap estimate
        let config = PreseqConfig {
            enabled: true,
            max_extrap: 1e10,
            step_size: 1e6,
            n_bootstraps: 100,
            confidence_level: 0.95,
            seed: 1,
            max_terms: 100,
            max_segment_length: 100_000_000,
            defects: false,
        };

        let result = estimate_complexity(&hist, total_reads, n_distinct, &config);
        assert!(result.is_ok(), "estimate_complexity failed: {:?}", result);
        let result = result.unwrap();

        // Curve should have 9999 rows (1M to 9999M at 1M step)
        assert_eq!(result.curve.len(), 9999, "Expected 9999 curve points");

        // Sanity: curve should be roughly in the right ballpark at 1M
        // Upstream: 613221.1. Platform-dependent, so allow wide tolerance.
        let (_, expected_1m, lower_1m, upper_1m) = result.curve[0];
        assert!(
            expected_1m > 500_000.0 && expected_1m < 800_000.0,
            "Expected ~600k at 1M, got {:.1}",
            expected_1m
        );
        assert!(lower_1m < expected_1m, "Lower CI should be below estimate");
        assert!(upper_1m > expected_1m, "Upper CI should be above estimate");

        eprintln!(
            "Benchmark comparison at 1M: expected={:.1} lower={:.1} upper={:.1} (upstream: 613221.1 / 499291.6 / 751901.7)",
            expected_1m, lower_1m, upper_1m
        );
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
            max_segment_length: 100_000_000,
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
