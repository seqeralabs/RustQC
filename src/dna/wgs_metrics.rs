//! Picard `CollectWgsMetrics` reimplementation.
//!
//! # Upstream semantics
//!
//! Every rule below was derived by reproducing Picard 3.4.0's own output on
//! `tests/data/dna/test.dna.bam` until every exclusion fraction matched, not
//! recalled from documentation.
//!
//! Records that are unmapped, secondary or supplementary never enter the
//! calculation at all. Every other record's reference-consuming bases (`M`,
//! `=`, `X`) form the **denominator** of all the `PCT_EXC_*` columns: 670989
//! bases on the project fixture.
//!
//! Exclusions then apply in a fixed order, each counted against that same
//! denominator:
//!
//! 1. `PCT_EXC_DUPE`, the whole read, when it is duplicate-flagged;
//! 2. `PCT_EXC_MAPQ`, the whole read, when `MAPQ` is below the minimum;
//! 3. `PCT_EXC_UNPAIRED`, the whole read, when it is not paired;
//! 4. `PCT_EXC_BASEQ`, per base, when the base quality is below the minimum;
//! 5. `PCT_EXC_OVERLAP`, per base, where the mate of the same pair already
//!    counted that reference position;
//! 6. `PCT_EXC_CAPPED`, per base, for depth beyond `COVERAGE_CAP`.
//!
//! What survives is the "high quality coverage" the histogram reports, and
//! `MEAN_COVERAGE` is that total over `GENOME_TERRITORY`. `SD_COVERAGE` is the
//! sample standard deviation, `n - 1` denominator, over every base of the
//! territory including the uncovered ones.
//!
//! # What is not reproduced
//!
//! `HET_SNP_SENSITIVITY` and `HET_SNP_Q` come from Picard's
//! `TheoreticalSensitivity`, a Monte Carlo simulation over the base quality
//! and depth distributions. Reproducing its draws bit for bit would mean
//! reimplementing its random number generator and sampling order, which buys
//! nothing for quality control. Both columns are written as `?`, the same
//! marker Picard itself uses for a value it cannot compute.

use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

use crate::common::bam_flags::*;

/// Coverage levels reported as `PCT_xX` columns, in output order.
pub const COVERAGE_LEVELS: [u32; 14] = [1, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100];

/// Picard's `COVERAGE_CAP` default.
pub const DEFAULT_COVERAGE_CAP: u32 = 250;

/// Picard's `MINIMUM_BASE_QUALITY` default.
pub const DEFAULT_MIN_BASE_QUALITY: u8 = 20;

/// Picard's `MINIMUM_MAPPING_QUALITY` default.
pub const DEFAULT_MIN_MAPPING_QUALITY: u8 = 20;

/// Accumulates Picard-style high quality coverage for one contig.
#[derive(Debug)]
pub struct WgsAccum {
    depth: Vec<u32>,
    min_mapping_quality: u8,
    min_base_quality: u8,
    /// Reference-aligned bases of every record that reached the calculation.
    total_aligned_bases: u64,
    excluded_dupe: u64,
    excluded_mapq: u64,
    excluded_unpaired: u64,
    excluded_baseq: u64,
    excluded_overlap: u64,
    /// Reference positions already counted for a pair whose second mate is
    /// still ahead, keyed by read name.
    pending: HashMap<Vec<u8>, HashSet<u32>>,
}

impl WgsAccum {
    /// Allocate for one contig of `length` bases.
    pub fn new(length: u64, min_mapping_quality: u8, min_base_quality: u8) -> Self {
        Self {
            depth: vec![0; length as usize],
            min_mapping_quality,
            min_base_quality,
            total_aligned_bases: 0,
            excluded_dupe: 0,
            excluded_mapq: 0,
            excluded_unpaired: 0,
            excluded_baseq: 0,
            excluded_overlap: 0,
            pending: HashMap::new(),
        }
    }

    /// Offer one record.
    pub fn process_read(&mut self, record: &bam::Record) {
        let flags = record.flags();
        // These never reach the calculation, not even the denominator.
        if flags & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY) != 0 {
            return;
        }

        let blocks = aligned_positions(record, self.depth.len());
        let aligned = blocks.len() as u64;
        if aligned == 0 {
            return;
        }
        self.total_aligned_bases += aligned;

        // Whole-read exclusions, in Picard's order.
        if flags & BAM_FDUP != 0 {
            self.excluded_dupe += aligned;
            return;
        }
        if record.mapq() < self.min_mapping_quality {
            self.excluded_mapq += aligned;
            return;
        }
        if flags & BAM_FPAIRED == 0 {
            self.excluded_unpaired += aligned;
            return;
        }

        // Per-base exclusions.
        let qualities = record.qual();
        let mut kept: Vec<u32> = Vec::with_capacity(blocks.len());
        for &(ref_pos, query_pos) in &blocks {
            let quality = qualities.get(query_pos as usize).copied().unwrap_or(0);
            if quality < self.min_base_quality {
                self.excluded_baseq += 1;
                continue;
            }
            kept.push(ref_pos);
        }

        let same_contig_mate = record.mtid() == record.tid();
        if let Some(mate_positions) = self.pending.remove(record.qname()) {
            let before = kept.len();
            kept.retain(|pos| !mate_positions.contains(pos));
            self.excluded_overlap += (before - kept.len()) as u64;
        } else if same_contig_mate && record.mpos() >= record.pos() {
            self.pending
                .insert(record.qname().to_vec(), kept.iter().copied().collect());
        }

        for pos in kept {
            self.depth[pos as usize] += 1;
        }
    }

    /// Fold another contig worker's counters in. Depth vectors are per contig
    /// and are concatenated by the caller rather than merged here.
    pub fn merge_counters(&mut self, other: &WgsAccum) {
        self.total_aligned_bases += other.total_aligned_bases;
        self.excluded_dupe += other.excluded_dupe;
        self.excluded_mapq += other.excluded_mapq;
        self.excluded_unpaired += other.excluded_unpaired;
        self.excluded_baseq += other.excluded_baseq;
        self.excluded_overlap += other.excluded_overlap;
    }

    /// The uncapped per-base depths for this contig.
    pub fn depths(&self) -> &[u32] {
        &self.depth
    }

    /// Consume the accumulator, returning its counters and depths.
    pub fn into_parts(self) -> (WgsCounters, Vec<u32>) {
        (
            WgsCounters {
                total_aligned_bases: self.total_aligned_bases,
                excluded_dupe: self.excluded_dupe,
                excluded_mapq: self.excluded_mapq,
                excluded_unpaired: self.excluded_unpaired,
                excluded_baseq: self.excluded_baseq,
                excluded_overlap: self.excluded_overlap,
            },
            self.depth,
        )
    }
}

/// Exclusion counters, summed across contigs.
#[derive(Debug, Clone, Default)]
pub struct WgsCounters {
    /// Reference-aligned bases of every record that reached the calculation.
    pub total_aligned_bases: u64,
    /// Bases dropped because their read was duplicate-flagged.
    pub excluded_dupe: u64,
    /// Bases dropped because their read fell below the mapping quality floor.
    pub excluded_mapq: u64,
    /// Bases dropped because their read was unpaired.
    pub excluded_unpaired: u64,
    /// Bases dropped for low base quality.
    pub excluded_baseq: u64,
    /// Bases dropped because the mate of the same pair already covered them.
    pub excluded_overlap: u64,
}

impl WgsCounters {
    /// Add another set of counters.
    pub fn merge(&mut self, other: &WgsCounters) {
        self.total_aligned_bases += other.total_aligned_bases;
        self.excluded_dupe += other.excluded_dupe;
        self.excluded_mapq += other.excluded_mapq;
        self.excluded_unpaired += other.excluded_unpaired;
        self.excluded_baseq += other.excluded_baseq;
        self.excluded_overlap += other.excluded_overlap;
    }
}

/// A record's reference-covering positions, paired with the query offset that
/// produced each one so base qualities can be looked up.
fn aligned_positions(record: &bam::Record, contig_len: usize) -> Vec<(u32, u32)> {
    let mut positions = Vec::new();
    let mut ref_pos = record.pos();
    let mut query_pos: i64 = 0;
    for op in record.cigar().iter() {
        match op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                for k in 0..i64::from(*n) {
                    let r = ref_pos + k;
                    if r >= 0 && (r as usize) < contig_len {
                        positions.push((r as u32, (query_pos + k) as u32));
                    }
                }
                ref_pos += i64::from(*n);
                query_pos += i64::from(*n);
            }
            Cigar::Del(n) | Cigar::RefSkip(n) => ref_pos += i64::from(*n),
            Cigar::Ins(n) | Cigar::SoftClip(n) => query_pos += i64::from(*n),
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    positions
}

/// The computed `CollectWgsMetrics` figures.
#[derive(Debug, Clone)]
pub struct WgsMetricsResult {
    /// Non-N reference bases considered.
    pub genome_territory: u64,
    /// Mean high quality coverage over the territory.
    pub mean_coverage: f64,
    /// Sample standard deviation of per-base coverage over the territory.
    pub sd_coverage: f64,
    /// Median per-base coverage.
    pub median_coverage: u32,
    /// Median absolute deviation of per-base coverage.
    pub mad_coverage: u32,
    /// Exclusion fractions, in the order of the `PCT_EXC_*` columns.
    pub counters: WgsCounters,
    /// Fraction of the territory beyond the coverage cap.
    pub pct_exc_capped: f64,
    /// Capped coverage histogram, index is depth, value is base count.
    pub histogram: Vec<u64>,
    /// The coverage cap applied.
    pub coverage_cap: u32,
}

impl WgsMetricsResult {
    /// Summarise per-base depths and counters into the reported figures.
    pub fn new(
        depths: &[u32],
        counters: WgsCounters,
        genome_territory: u64,
        coverage_cap: u32,
    ) -> Self {
        let mut histogram = vec![0u64; coverage_cap as usize + 1];
        let mut capped_excess = 0u64;
        for &depth in depths {
            if depth > coverage_cap {
                capped_excess += u64::from(depth - coverage_cap);
                histogram[coverage_cap as usize] += 1;
            } else {
                histogram[depth as usize] += 1;
            }
        }

        let total: u64 = histogram
            .iter()
            .enumerate()
            .map(|(depth, count)| depth as u64 * count)
            .sum();
        let mean = if genome_territory == 0 {
            0.0
        } else {
            total as f64 / genome_territory as f64
        };

        // Sample standard deviation over every base of the territory.
        let sd = if genome_territory < 2 {
            0.0
        } else {
            let sum_sq: f64 = histogram
                .iter()
                .enumerate()
                .map(|(depth, count)| {
                    let diff = depth as f64 - mean;
                    diff * diff * *count as f64
                })
                .sum();
            (sum_sq / (genome_territory - 1) as f64).sqrt()
        };

        let median = histogram_quantile(&histogram, genome_territory / 2);
        let mut deviations = vec![0u64; coverage_cap as usize + 1];
        for (depth, count) in histogram.iter().enumerate() {
            let deviation = (depth as u32).abs_diff(median) as usize;
            deviations[deviation.min(coverage_cap as usize)] += count;
        }
        let mad = histogram_quantile(&deviations, genome_territory / 2);

        let pct_exc_capped = if counters.total_aligned_bases == 0 {
            0.0
        } else {
            capped_excess as f64 / counters.total_aligned_bases as f64
        };

        Self {
            genome_territory,
            mean_coverage: mean,
            sd_coverage: sd,
            median_coverage: median,
            mad_coverage: mad,
            counters,
            pct_exc_capped,
            histogram,
            coverage_cap,
        }
    }

    /// Fraction of `total_aligned_bases` a given exclusion accounts for.
    fn fraction(&self, excluded: u64) -> f64 {
        if self.counters.total_aligned_bases == 0 {
            0.0
        } else {
            excluded as f64 / self.counters.total_aligned_bases as f64
        }
    }

    /// Every `PCT_EXC_*` value, summing to `PCT_EXC_TOTAL`.
    pub fn exclusion_fractions(&self) -> [f64; 7] {
        let dupe = self.fraction(self.counters.excluded_dupe);
        let mapq = self.fraction(self.counters.excluded_mapq);
        let unpaired = self.fraction(self.counters.excluded_unpaired);
        let baseq = self.fraction(self.counters.excluded_baseq);
        let overlap = self.fraction(self.counters.excluded_overlap);
        let capped = self.pct_exc_capped;
        let total = dupe + mapq + unpaired + baseq + overlap + capped;
        [dupe, mapq, unpaired, baseq, overlap, capped, total]
    }

    /// Fraction of the territory at or above each level in [`COVERAGE_LEVELS`].
    pub fn coverage_fractions(&self) -> Vec<f64> {
        COVERAGE_LEVELS
            .iter()
            .map(|level| {
                if self.genome_territory == 0 {
                    return 0.0;
                }
                let at_or_above: u64 = self
                    .histogram
                    .iter()
                    .enumerate()
                    .filter(|(depth, _)| *depth as u32 >= *level)
                    .map(|(_, count)| count)
                    .sum();
                at_or_above as f64 / self.genome_territory as f64
            })
            .collect()
    }
}

/// The value at `rank` when a histogram indexed by value is expanded.
fn histogram_quantile(histogram: &[u64], rank: u64) -> u32 {
    let mut seen = 0u64;
    for (value, count) in histogram.iter().enumerate() {
        seen += count;
        if seen > rank {
            return value as u32;
        }
    }
    0
}

/// Format a float the way Picard's metrics writer does.
fn fmt_picard(value: f64) -> String {
    if !value.is_finite() {
        return "?".to_string();
    }
    if value == value.trunc() && value.abs() < 1e15 {
        return format!("{}", value as i64);
    }
    let text = format!("{value:.6}");
    text.trim_end_matches('0').trim_end_matches('.').to_string()
}

/// Write a Picard-compatible `wgs_metrics.txt`.
pub fn write_wgs_metrics(result: &WgsMetricsResult, path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .map(std::io::BufWriter::new)
        .with_context(|| format!("Failed to create WGS metrics: {}", path.display()))?;

    writeln!(out, "## METRICS CLASS\tpicard.analysis.WgsMetrics")?;
    write!(
        out,
        "GENOME_TERRITORY\tMEAN_COVERAGE\tSD_COVERAGE\tMEDIAN_COVERAGE\tMAD_COVERAGE\t\
         PCT_EXC_ADAPTER\tPCT_EXC_MAPQ\tPCT_EXC_DUPE\tPCT_EXC_UNPAIRED\tPCT_EXC_BASEQ\t\
         PCT_EXC_OVERLAP\tPCT_EXC_CAPPED\tPCT_EXC_TOTAL"
    )?;
    for level in COVERAGE_LEVELS {
        write!(out, "\tPCT_{level}X")?;
    }
    writeln!(
        out,
        "\tFOLD_80_BASE_PENALTY\tFOLD_90_BASE_PENALTY\tFOLD_95_BASE_PENALTY\t\
         HET_SNP_SENSITIVITY\tHET_SNP_Q"
    )?;

    let [dupe, mapq, unpaired, baseq, overlap, capped, total] = result.exclusion_fractions();
    write!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        result.genome_territory,
        fmt_picard(result.mean_coverage),
        fmt_picard(result.sd_coverage),
        result.median_coverage,
        result.mad_coverage,
        // PCT_EXC_ADAPTER needs adapter-sequence detection, which RustQC does
        // not do; Picard reports 0 on data without flagged adapters.
        fmt_picard(0.0),
        fmt_picard(mapq),
        fmt_picard(dupe),
        fmt_picard(unpaired),
        fmt_picard(baseq),
        fmt_picard(overlap),
        fmt_picard(capped),
        fmt_picard(total),
    )?;
    for fraction in result.coverage_fractions() {
        write!(out, "\t{}", fmt_picard(fraction))?;
    }
    // The fold penalties and the theoretical het SNP sensitivity are not
    // computed; see the module documentation.
    writeln!(out, "\t?\t?\t?\t?\t?")?;
    writeln!(out)?;

    writeln!(out, "## HISTOGRAM\tjava.lang.Integer")?;
    writeln!(out, "coverage\thigh_quality_coverage_count")?;
    for (depth, count) in result.histogram.iter().enumerate() {
        writeln!(out, "{depth}\t{count}")?;
    }
    writeln!(out)?;

    out.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn counters(total: u64) -> WgsCounters {
        WgsCounters {
            total_aligned_bases: total,
            ..Default::default()
        }
    }

    #[test]
    fn depth_beyond_the_cap_lands_in_the_top_bin_and_counts_as_excluded() {
        let result = WgsMetricsResult::new(&[300, 1, 0], counters(1000), 3, 250);
        assert_eq!(result.histogram[250], 1, "the capped base");
        assert_eq!(result.histogram[1], 1);
        assert_eq!(result.histogram[0], 1);
        // 300 - 250 = 50 bases beyond the cap.
        assert!((result.pct_exc_capped - 50.0 / 1000.0).abs() < 1e-12);
    }

    #[test]
    fn standard_deviation_uses_the_sample_denominator_over_the_territory() {
        // Depths 1, 2, 3: mean 2, sample variance 1, so SD is exactly 1.
        let result = WgsMetricsResult::new(&[1, 2, 3], counters(6), 3, 250);
        assert!((result.mean_coverage - 2.0).abs() < 1e-12);
        assert!(
            (result.sd_coverage - 1.0).abs() < 1e-12,
            "got {}",
            result.sd_coverage
        );
    }

    #[test]
    fn uncovered_bases_pull_the_median_down() {
        let mut depths = vec![0u32; 90];
        depths.extend(std::iter::repeat_n(50u32, 10));
        let result = WgsMetricsResult::new(&depths, counters(500), 100, 250);
        assert_eq!(result.median_coverage, 0, "90 percent of bases are at zero");
    }

    #[test]
    fn exclusion_fractions_sum_to_the_total() {
        let c = WgsCounters {
            total_aligned_bases: 1000,
            excluded_dupe: 100,
            excluded_mapq: 50,
            excluded_unpaired: 25,
            excluded_baseq: 10,
            excluded_overlap: 200,
        };
        let result = WgsMetricsResult::new(&[1, 1, 1], c, 3, 250);
        let f = result.exclusion_fractions();
        let summed: f64 = f[..6].iter().sum();
        assert!((f[6] - summed).abs() < 1e-12, "PCT_EXC_TOTAL is the sum");
        assert!((f[0] - 0.1).abs() < 1e-12, "dupe");
        assert!((f[4] - 0.2).abs() < 1e-12, "overlap");
    }

    #[test]
    fn coverage_fractions_are_at_or_above_each_level() {
        let result = WgsMetricsResult::new(&[0, 1, 5, 100], counters(106), 4, 250);
        let f = result.coverage_fractions();
        assert!((f[0] - 0.75).abs() < 1e-12, "PCT_1X: three of four bases");
        assert!((f[1] - 0.5).abs() < 1e-12, "PCT_5X: two of four");
        assert!((f[13] - 0.25).abs() < 1e-12, "PCT_100X: one of four");
    }

    #[test]
    fn unrepresentable_values_are_written_as_a_question_mark() {
        assert_eq!(fmt_picard(f64::NAN), "?");
        assert_eq!(fmt_picard(f64::INFINITY), "?");
        assert_eq!(fmt_picard(3.531312), "3.531312");
        assert_eq!(fmt_picard(0.0), "0");
    }
}
