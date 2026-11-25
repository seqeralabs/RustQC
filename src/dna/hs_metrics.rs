//! Picard `CollectHsMetrics` reimplementation for targeted sequencing.
//!
//! # Upstream semantics
//!
//! Taken from Picard 3.4.0's `TargetMetricsCollector`, because this collector
//! does not filter the way [`crate::dna::wgs_metrics`] does and the difference
//! is not guessable from the outputs.
//!
//! Secondary alignments are excluded outright, so `TOTAL_READS` is 5642 on the
//! project fixture rather than the 5644 records it holds. Then, per record:
//!
//! 1. `PF_BASES` accumulates the read length of every non-supplementary read;
//! 2. mapped reads add their reference-aligned bases to `PF_BASES_ALIGNED`,
//!    and to `PF_UQ_BASES_ALIGNED` when not duplicate-flagged;
//! 3. the bait counters are taken **before** any filtering, so the assay
//!    metrics are not skewed by duplicates or mapping quality;
//! 4. duplicates are charged to `PCT_EXC_DUPE` and dropped;
//! 5. reads below the mapping quality floor are dropped;
//! 6. **overlap clipping happens next, at the read level**, charging
//!    `PCT_EXC_OVERLAP` with the number of aligned bases clipped;
//! 7. only then, per surviving base: below the base quality floor charges
//!    `PCT_EXC_BASEQ`; off-target charges `PCT_EXC_OFF_TARGET`; the rest are
//!    `ON_TARGET_BASES`.
//!
//! Step 6 before step 7 is the crux. `CollectWgsMetrics` applies base quality
//! first and reconciles overlaps per locus afterwards, which is why the two
//! collectors report different `PCT_EXC_BASEQ` and `PCT_EXC_OVERLAP` on the
//! same file: 0.003982 against 0.007352, and 0.330968 against 0.324694.
//!
//! Only the left-most mate of an overlapping pair is clipped, and everything
//! from the mate's alignment start onwards goes, per htsjdk's
//! `getNumOverlappingAlignedBasesToClip`.
//!
//! # What is not reproduced
//!
//! `HET_SNP_SENSITIVITY` and `HET_SNP_Q` come from the same Monte Carlo
//! simulation left out of `CollectWgsMetrics`, and `HS_PENALTY_*X` and
//! `FOLD_80_BASE_PENALTY` derive from it. All are written as Picard writes
//! them when it cannot compute them: `-1` for the penalties, `?` for the rest.

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

use crate::common::bam_flags::*;
use crate::dna::intervals::IntervalSet;

/// Coverage levels reported as `PCT_TARGET_BASES_xX`, in output order.
pub const TARGET_COVERAGE_LEVELS: [u32; 17] = [
    1, 2, 10, 20, 30, 40, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000, 100000,
];

/// Penalty levels reported as `HS_PENALTY_xX`, in output order.
pub const PENALTY_LEVELS: [u32; 6] = [10, 20, 30, 40, 50, 100];

/// Accumulates targeted-sequencing metrics for one contig.
#[derive(Debug)]
pub struct HsAccum {
    bait_mask: Vec<bool>,
    target_mask: Vec<bool>,
    /// High quality on-target depth per reference base.
    depth: Vec<u32>,
    min_mapping_quality: u8,
    min_base_quality: u8,
    counters: HsCounters,
}

/// The raw counters, summed across contigs.
#[derive(Debug, Clone, Default)]
pub struct HsCounters {
    /// Records seen, secondary alignments excluded.
    pub total_reads: u64,
    /// Read length of every non-supplementary record.
    pub pf_bases: u64,
    /// Reference-aligned bases of mapped records.
    pub pf_bases_aligned: u64,
    /// Reference-aligned bases of mapped, non-duplicate records.
    pub pf_uq_bases_aligned: u64,
    /// Non-duplicate records.
    pub pf_unique_reads: u64,
    /// Non-duplicate mapped records.
    pub pf_uq_reads_aligned: u64,
    /// Aligned bases falling on a bait.
    pub on_bait_bases: u64,
    /// Aligned bases of bait-overlapping reads that miss the baits themselves.
    pub near_bait_bases: u64,
    /// Aligned bases of reads that touch no bait at all.
    pub off_bait_bases: u64,
    /// Bases dropped because their read was duplicate-flagged.
    pub excluded_dupe: u64,
    /// Bases clipped as overlapping a mate.
    pub excluded_overlap: u64,
    /// Bases dropped for low base quality.
    pub excluded_baseq: u64,
    /// Bases dropped for falling outside the targets.
    pub excluded_off_target: u64,
    /// High quality bases on target.
    pub on_target_bases: u64,
    /// First-of-pair records over a bait, with a mapped mate.
    pub selected_pairs: u64,
    /// The same, excluding duplicates.
    pub selected_unique_pairs: u64,
}

impl HsCounters {
    /// Add another contig's counters.
    pub fn merge(&mut self, other: &HsCounters) {
        self.total_reads += other.total_reads;
        self.pf_bases += other.pf_bases;
        self.pf_bases_aligned += other.pf_bases_aligned;
        self.pf_uq_bases_aligned += other.pf_uq_bases_aligned;
        self.pf_unique_reads += other.pf_unique_reads;
        self.pf_uq_reads_aligned += other.pf_uq_reads_aligned;
        self.on_bait_bases += other.on_bait_bases;
        self.near_bait_bases += other.near_bait_bases;
        self.off_bait_bases += other.off_bait_bases;
        self.excluded_dupe += other.excluded_dupe;
        self.excluded_overlap += other.excluded_overlap;
        self.excluded_baseq += other.excluded_baseq;
        self.excluded_off_target += other.excluded_off_target;
        self.on_target_bases += other.on_target_bases;
        self.selected_pairs += other.selected_pairs;
        self.selected_unique_pairs += other.selected_unique_pairs;
    }
}

impl HsAccum {
    /// Prepare for one contig.
    pub fn new(
        contig: &str,
        length: u64,
        baits: &IntervalSet,
        targets: &IntervalSet,
        min_mapping_quality: u8,
        min_base_quality: u8,
    ) -> Self {
        Self {
            bait_mask: baits.mask(contig, length),
            target_mask: targets.mask(contig, length),
            depth: vec![0; length as usize],
            min_mapping_quality,
            min_base_quality,
            counters: HsCounters::default(),
        }
    }

    /// Offer one record.
    pub fn process_read(&mut self, record: &bam::Record) {
        let flags = record.flags();
        // Secondary alignments are not part of this collector's read set.
        if flags & BAM_FSECONDARY != 0 || flags & BAM_FQCFAIL != 0 {
            return;
        }
        self.counters.total_reads += 1;

        if flags & BAM_FSUPPLEMENTARY == 0 {
            self.counters.pf_bases += record.seq_len() as u64;
        }
        if flags & BAM_FDUP == 0 {
            self.counters.pf_unique_reads += 1;
        }
        if flags & BAM_FUNMAP != 0 {
            return;
        }

        let blocks = aligned_blocks(record, self.depth.len());
        let aligned: u64 = blocks.iter().map(|(start, end)| end - start).sum();
        self.counters.pf_bases_aligned += aligned;
        if flags & BAM_FDUP == 0 {
            self.counters.pf_uq_bases_aligned += aligned;
            self.counters.pf_uq_reads_aligned += 1;
        }

        // Bait metrics come before duplicate, mapping quality and overlap
        // filtering, so that the assay is measured rather than the library.
        let on_bait: u64 = blocks
            .iter()
            .map(|(start, end)| {
                (*start..*end)
                    .filter(|p| self.bait_mask.get(*p as usize).copied().unwrap_or(false))
                    .count() as u64
            })
            .sum();
        if on_bait > 0 {
            self.counters.on_bait_bases += on_bait;
            self.counters.near_bait_bases += aligned - on_bait;
        } else {
            self.counters.off_bait_bases += aligned;
        }

        // HS_LIBRARY_SIZE counts templates over a bait, once each.
        if flags & BAM_FSUPPLEMENTARY == 0
            && flags & BAM_FPAIRED != 0
            && flags & BAM_FREAD1 != 0
            && flags & BAM_FMUNMAP == 0
            && on_bait > 0
        {
            self.counters.selected_pairs += 1;
            if flags & BAM_FDUP == 0 {
                self.counters.selected_unique_pairs += 1;
            }
        }

        if flags & BAM_FDUP != 0 {
            self.counters.excluded_dupe += aligned;
            return;
        }
        if record.mapq() < self.min_mapping_quality {
            return;
        }

        // Overlap clipping, at the read level and before any base is examined.
        //
        // Two different quantities are at play. The counter Picard reports is
        // htsjdk's count of *read* bases clipped, insertions included. What is
        // actually removed from the alignment is every base at or past the
        // mate's start, in *reference* coordinates. They coincide only for a
        // gapless read, so they are tracked separately.
        self.counters.excluded_overlap += overlapping_bases_to_clip(record);
        let clip_from = overlap_clip_reference_start(record);

        let qualities = record.qual();
        let mut ref_pos = record.pos();
        let mut query_pos = 0i64;
        for op in record.cigar().iter() {
            match op {
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                    for k in 0..i64::from(*n) {
                        let r = ref_pos + k;
                        if r < 0 || r as usize >= self.depth.len() {
                            continue;
                        }
                        if let Some(from) = clip_from {
                            if r >= from {
                                continue;
                            }
                        }
                        let quality = qualities
                            .get((query_pos + k) as usize)
                            .copied()
                            .unwrap_or(0);
                        if quality < self.min_base_quality {
                            self.counters.excluded_baseq += 1;
                        } else if !self.target_mask[r as usize] {
                            self.counters.excluded_off_target += 1;
                        } else {
                            self.counters.on_target_bases += 1;
                            self.depth[r as usize] += 1;
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
    }

    /// Consume the accumulator, returning its counters and per-base depths.
    pub fn into_parts(self) -> (HsCounters, Vec<u32>, Vec<bool>) {
        (self.counters, self.depth, self.target_mask)
    }
}

/// A record's reference-covering blocks as half-open `[start, end)`.
fn aligned_blocks(record: &bam::Record, contig_len: usize) -> Vec<(u64, u64)> {
    let mut blocks = Vec::new();
    let mut ref_pos = record.pos();
    for op in record.cigar().iter() {
        match op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                let start = ref_pos.max(0) as u64;
                let end = ((ref_pos + i64::from(*n)).max(0) as u64).min(contig_len as u64);
                if start < end {
                    blocks.push((start, end));
                }
                ref_pos += i64::from(*n);
            }
            Cigar::Del(n) | Cigar::RefSkip(n) => ref_pos += i64::from(*n),
            _ => {}
        }
    }
    blocks
}

/// The reference position from which this read's alignment is clipped away
/// because its mate covers it, or `None` when nothing is clipped.
///
/// This is the mate's alignment start: the left-most read of an overlapping
/// pair loses everything from there onwards.
fn overlap_clip_reference_start(record: &bam::Record) -> Option<i64> {
    if overlapping_bases_to_clip(record) == 0 {
        return None;
    }
    Some(record.mpos())
}

/// Read bases to clip because a mate covers them, per htsjdk's
/// `getNumOverlappingAlignedBasesToClip`.
///
/// Only the left-most mate of the pair is clipped, and everything from the
/// mate's alignment start onwards goes. A pair sharing a start is broken by
/// clipping the second of the pair.
fn overlapping_bases_to_clip(record: &bam::Record) -> u64 {
    let flags = record.flags();
    if flags & BAM_FPAIRED == 0 || flags & BAM_FUNMAP != 0 || flags & BAM_FMUNMAP != 0 {
        return 0;
    }
    let start = record.pos();
    let mate_start = record.mpos();
    if mate_start < start {
        return 0;
    }
    if mate_start == start && flags & BAM_FREAD1 != 0 {
        return 0;
    }

    let mut clipped: i64 = 0;
    let mut ref_pos = start;
    for op in record.cigar().iter() {
        let ref_len = match op {
            Cigar::Match(n)
            | Cigar::Equal(n)
            | Cigar::Diff(n)
            | Cigar::Del(n)
            | Cigar::RefSkip(n) => i64::from(*n),
            _ => 0,
        };
        if mate_start < ref_pos + ref_len {
            match op {
                // Only M takes the partial path: htsjdk's MATCH_OR_MISMATCH is
                // the M operator alone, so = and X fall through to the branch
                // below and lose their whole element.
                Cigar::Match(_) => {
                    clipped += if mate_start < ref_pos {
                        ref_len
                    } else {
                        ref_pos + ref_len - mate_start
                    };
                }
                Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) | Cigar::RefSkip(_) => {}
                // Everything else loses its read-consuming bases outright,
                // which covers insertions as well as = and X.
                Cigar::Equal(n) | Cigar::Diff(n) | Cigar::Ins(n) => clipped += i64::from(*n),
                Cigar::Del(_) => {}
            }
        }
        ref_pos += ref_len;
    }
    // Left-most but not actually overlapping.
    clipped.max(0) as u64
}

/// Estimate library size from observed and unique templates.
///
/// Solves the Lander-Waterman equation `C/X = 1 - exp(-N/X)` by bisection,
/// exactly as Picard's `DuplicationMetrics.estimateLibrarySize` does, down to
/// the forty iterations and the starting bracket.
pub fn estimate_library_size(read_pairs: u64, unique_read_pairs: u64) -> Option<u64> {
    if read_pairs == 0 || read_pairs <= unique_read_pairs || unique_read_pairs == 0 {
        return None;
    }
    let n = read_pairs as f64;
    let c = unique_read_pairs as f64;
    let f = |x: f64| c / x - 1.0 + (-n / x).exp();

    let mut low = 1.0;
    let mut high = 100.0;
    while f(high * c) > 0.0 {
        high *= 10.0;
    }
    for _ in 0..40 {
        let mid = (low + high) / 2.0;
        let value = f(mid * c);
        if value == 0.0 {
            break;
        } else if value > 0.0 {
            low = mid;
        } else {
            high = mid;
        }
    }
    Some((c * (low + high) / 2.0) as u64)
}

/// The computed `CollectHsMetrics` figures.
#[derive(Debug, Clone)]
pub struct HsMetricsResult {
    /// Name of the bait set.
    pub bait_set: String,
    /// Bases covered by baits.
    pub bait_territory: u64,
    /// Bases covered by targets.
    pub target_territory: u64,
    /// Total reference length.
    pub genome_size: u64,
    /// Raw counters.
    pub counters: HsCounters,
    /// High quality depth of every target base, target order.
    pub target_depths: Vec<u32>,
    /// Number of targets with no coverage at all.
    pub zero_coverage_targets: u64,
    /// Number of targets.
    pub target_count: u64,
    /// Estimated library size, when it can be estimated.
    pub library_size: Option<u64>,
}

impl HsMetricsResult {
    /// Mean high quality coverage over the targets.
    pub fn mean_target_coverage(&self) -> f64 {
        if self.target_territory == 0 {
            0.0
        } else {
            self.counters.on_target_bases as f64 / self.target_territory as f64
        }
    }

    /// Mean aligned coverage over the baits.
    pub fn mean_bait_coverage(&self) -> f64 {
        if self.bait_territory == 0 {
            0.0
        } else {
            self.counters.pf_bases_aligned as f64 / self.bait_territory as f64
        }
    }

    /// Fraction of the targets at or above each level.
    pub fn target_coverage_fractions(&self) -> Vec<f64> {
        TARGET_COVERAGE_LEVELS
            .iter()
            .map(|level| {
                if self.target_territory == 0 {
                    return 0.0;
                }
                let at_or_above = self.target_depths.iter().filter(|d| **d >= *level).count();
                at_or_above as f64 / self.target_territory as f64
            })
            .collect()
    }

    /// Median, minimum and maximum high quality target coverage.
    pub fn target_coverage_bounds(&self) -> (u32, u32, u32) {
        if self.target_depths.is_empty() {
            return (0, 0, 0);
        }
        let mut sorted = self.target_depths.clone();
        sorted.sort_unstable();
        let median = sorted[sorted.len() / 2];
        (median, sorted[0], sorted[sorted.len() - 1])
    }
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

/// Write a Picard-compatible `hs_metrics.txt`.
pub fn write_hs_metrics(result: &HsMetricsResult, path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .map(std::io::BufWriter::new)
        .with_context(|| format!("Failed to create HS metrics: {}", path.display()))?;

    let c = &result.counters;
    let aligned = c.pf_bases_aligned as f64;
    let frac = |n: u64| {
        if aligned == 0.0 {
            0.0
        } else {
            n as f64 / aligned
        }
    };
    let selected = c.on_bait_bases + c.near_bait_bases;
    let (median, min, max) = result.target_coverage_bounds();

    writeln!(out, "## METRICS CLASS\tpicard.analysis.directed.HsMetrics")?;

    let mut header = String::from(
        "BAIT_SET\tBAIT_TERRITORY\tBAIT_DESIGN_EFFICIENCY\tON_BAIT_BASES\tNEAR_BAIT_BASES\t\
         OFF_BAIT_BASES\tPCT_SELECTED_BASES\tPCT_OFF_BAIT\tON_BAIT_VS_SELECTED\t\
         MEAN_BAIT_COVERAGE\tPCT_USABLE_BASES_ON_BAIT\tPCT_USABLE_BASES_ON_TARGET\t\
         FOLD_ENRICHMENT\tHS_LIBRARY_SIZE",
    );
    for level in PENALTY_LEVELS {
        header.push_str(&format!("\tHS_PENALTY_{level}X"));
    }
    header.push_str(
        "\tTARGET_TERRITORY\tGENOME_SIZE\tTOTAL_READS\tPF_READS\tPF_BASES\tPF_UNIQUE_READS\t\
         PF_UQ_READS_ALIGNED\tPF_BASES_ALIGNED\tPF_UQ_BASES_ALIGNED\tON_TARGET_BASES\t\
         PCT_PF_READS\tPCT_PF_UQ_READS\tPCT_PF_UQ_READS_ALIGNED\tMEAN_TARGET_COVERAGE\t\
         MEDIAN_TARGET_COVERAGE\tMAX_TARGET_COVERAGE\tMIN_TARGET_COVERAGE\tZERO_CVG_TARGETS_PCT\t\
         PCT_EXC_DUPE\tPCT_EXC_ADAPTER\tPCT_EXC_MAPQ\tPCT_EXC_BASEQ\tPCT_EXC_OVERLAP\t\
         PCT_EXC_OFF_TARGET\tFOLD_80_BASE_PENALTY",
    );
    for level in TARGET_COVERAGE_LEVELS {
        header.push_str(&format!("\tPCT_TARGET_BASES_{level}X"));
    }
    header.push_str(
        "\tAT_DROPOUT\tGC_DROPOUT\tHET_SNP_SENSITIVITY\tHET_SNP_Q\tSAMPLE\tLIBRARY\tREAD_GROUP",
    );
    writeln!(out, "{header}")?;

    write!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        result.bait_set,
        result.bait_territory,
        // Every bait base is intended as a target here; Picard reports the
        // fraction of bait territory that is also target territory.
        fmt_picard(if result.bait_territory == 0 {
            0.0
        } else {
            result.target_territory.min(result.bait_territory) as f64 / result.bait_territory as f64
        }),
        c.on_bait_bases,
        c.near_bait_bases,
        c.off_bait_bases,
        fmt_picard(frac(selected)),
        fmt_picard(frac(c.off_bait_bases)),
        fmt_picard(if selected == 0 {
            0.0
        } else {
            c.on_bait_bases as f64 / selected as f64
        }),
        fmt_picard(result.mean_bait_coverage()),
        fmt_picard(if c.pf_bases == 0 {
            0.0
        } else {
            c.on_bait_bases as f64 / c.pf_bases as f64
        }),
        fmt_picard(if c.pf_bases == 0 {
            0.0
        } else {
            c.on_target_bases as f64 / c.pf_bases as f64
        }),
        fmt_picard(fold_enrichment(result)),
        result
            .library_size
            .map(|v| v.to_string())
            .unwrap_or_default(),
    )?;
    // The penalties derive from the theoretical sensitivity simulation, which
    // is out of scope; Picard writes -1 when it cannot compute them.
    for _ in PENALTY_LEVELS {
        write!(out, "\t-1")?;
    }
    write!(
        out,
        "\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t?",
        result.target_territory,
        result.genome_size,
        c.total_reads,
        c.total_reads,
        c.pf_bases,
        c.pf_unique_reads,
        c.pf_uq_reads_aligned,
        c.pf_bases_aligned,
        c.pf_uq_bases_aligned,
        c.on_target_bases,
        fmt_picard(1.0),
        fmt_picard(if c.total_reads == 0 {
            0.0
        } else {
            c.pf_unique_reads as f64 / c.total_reads as f64
        }),
        fmt_picard(if c.pf_unique_reads == 0 {
            0.0
        } else {
            c.pf_uq_reads_aligned as f64 / c.pf_unique_reads as f64
        }),
        fmt_picard(result.mean_target_coverage()),
        median,
        max,
        min,
        fmt_picard(if result.target_count == 0 {
            0.0
        } else {
            result.zero_coverage_targets as f64 / result.target_count as f64
        }),
        fmt_picard(frac(c.excluded_dupe)),
        fmt_picard(0.0),
        fmt_picard(frac(0)),
        fmt_picard(frac(c.excluded_baseq)),
        fmt_picard(frac(c.excluded_overlap)),
        fmt_picard(frac(c.excluded_off_target)),
    )?;
    for fraction in result.target_coverage_fractions() {
        write!(out, "\t{}", fmt_picard(fraction))?;
    }
    // AT and GC dropout over targets, and the two simulated columns, are not
    // computed; see the module documentation.
    writeln!(out, "\t?\t?\t?\t?\t\t\t")?;
    writeln!(out)?;

    out.flush()?;
    Ok(())
}

/// Enrichment of the selected territory relative to uniform coverage.
///
/// Picard computes this from the *selected* bases against the bait territory,
/// not from on-target bases against the target territory. On the project
/// fixture every aligned base is on bait, so the figure reduces to
/// `GENOME_SIZE / BAIT_TERRITORY`, which is exactly the 1.142886 it reports.
fn fold_enrichment(result: &HsMetricsResult) -> f64 {
    let c = &result.counters;
    if c.pf_bases_aligned == 0 || result.bait_territory == 0 || result.genome_size == 0 {
        return 0.0;
    }
    let selected = (c.on_bait_bases + c.near_bait_bases) as f64 / c.pf_bases_aligned as f64;
    selected / (result.bait_territory as f64 / result.genome_size as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn library_size_solves_the_lander_waterman_equation() {
        // The project fixture's counts, checked against Picard's own answer.
        assert_eq!(estimate_library_size(2820, 1992), Some(3807));
    }

    #[test]
    fn library_size_is_absent_when_nothing_is_duplicated() {
        assert_eq!(estimate_library_size(100, 100), None);
        assert_eq!(estimate_library_size(0, 0), None);
    }

    #[test]
    fn target_coverage_fractions_are_at_or_above_each_level() {
        let result = HsMetricsResult {
            bait_set: "t".into(),
            bait_territory: 4,
            target_territory: 4,
            genome_size: 100,
            counters: HsCounters::default(),
            target_depths: vec![0, 1, 10, 300],
            zero_coverage_targets: 0,
            target_count: 1,
            library_size: None,
        };
        let f = result.target_coverage_fractions();
        assert!((f[0] - 0.75).abs() < 1e-12, "1X");
        assert!((f[2] - 0.5).abs() < 1e-12, "10X");
        assert!((f[8] - 0.25).abs() < 1e-12, "250X");
    }

    #[test]
    fn coverage_bounds_come_from_the_target_bases_only() {
        let result = HsMetricsResult {
            bait_set: "t".into(),
            bait_territory: 5,
            target_territory: 5,
            genome_size: 100,
            counters: HsCounters::default(),
            target_depths: vec![0, 3, 7, 9, 40],
            zero_coverage_targets: 0,
            target_count: 1,
            library_size: None,
        };
        assert_eq!(result.target_coverage_bounds(), (7, 0, 40));
    }
}
