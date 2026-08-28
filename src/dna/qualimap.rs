//! Qualimap `bamqc` reimplementation.
//!
//! # Upstream semantics
//!
//! Derived by reproducing Qualimap 2.3's own output on the project fixture
//! until each figure matched. Several rules are surprising and none of them
//! are guessable, so they are recorded here.
//!
//! **Windows.** The reference is split into `ceil(len / ceil(len / 400))`
//! windows, which is 397 windows of 101 bases on the 40001 base fixture, not
//! the round 400 the option name suggests.
//!
//! **Coverage.** Every primary mapped record contributes, with no duplicate,
//! mapping quality or base quality filtering and **no mate-overlap
//! correction**. Deletions count as covered. That is why Qualimap reports 16.77
//! mean coverage where mosdepth reports 6.20 on the same file: they are
//! measuring different things, and neither is wrong.
//!
//! **Mapping quality.** The global figure is the mean of the per-window means,
//! where a window with no reads contributes zero. That is why it reads 2.4178
//! rather than about 60. The per-position histogram truncates the mean rather
//! than rounding it.
//!
//! **Base composition.** Bases are counted in reference orientation, so
//! reverse-strand reads are reverse-complemented, but the clipped span that
//! selects which positions count is taken in *sequencing* orientation. Mixing
//! the two orientations is what Qualimap does; matching it means doing the
//! same.
//!
//! **Mismatches** are the `NM` tag less inserted bases only. Deleted bases are
//! not subtracted, which is what puts the fixture at 1350 rather than 1340.

use rust_htslib::bam;
use rust_htslib::bam::record::{Aux, Cigar};
use std::collections::BTreeMap;

use crate::common::bam_flags::*;

/// Qualimap's default target number of windows.
pub const DEFAULT_NUM_WINDOWS: usize = 400;

/// Highest coverage level reported in the genome fraction table.
const MAX_FRACTION_LEVEL: u32 = 51;

/// Per-contig accumulation for one alignment file.
#[derive(Debug)]
pub struct QualimapAccum {
    contig: String,
    length: u64,
    window_size: u64,
    /// Coverage per reference base, counting `M`, `=`, `X` and `D`.
    coverage: Vec<u32>,
    /// Sum of mapping quality over the reads covering each base.
    mapq_sum: Vec<u64>,
    /// Per-window sum of insert sizes and the number of reads contributing.
    insert_window_sum: Vec<i64>,
    insert_window_count: Vec<u64>,
    counters: QualimapCounters,
}

/// Read-level counters, summed across contigs.
#[derive(Debug, Clone, Default)]
pub struct QualimapCounters {
    /// Records seen, secondary alignments excluded and counted separately.
    pub reads: u64,
    /// Secondary alignments.
    pub secondary: u64,
    /// Mapped records.
    pub mapped: u64,
    /// Duplicate-flagged records.
    pub duplicates: u64,
    /// Mapped first-in-pair records with a mapped mate.
    pub paired_first: u64,
    /// Mapped second-in-pair records with a mapped mate.
    pub paired_second: u64,
    /// Mapped paired records whose mate is also mapped.
    pub paired_both: u64,
    /// Mapped paired records whose mate is not mapped.
    pub singletons: u64,
    /// Reference-consuming aligned bases, `M`, `=` and `X`.
    pub sequenced_bases: u64,
    /// Those plus deleted bases.
    pub mapped_bases: u64,
    /// Sum of the `NM` tag over mapped records.
    pub edit_distance: u64,
    /// Inserted bases.
    pub insertions: u64,
    /// Deleted bases.
    pub deletions: u64,
    /// Records carrying at least one insertion.
    pub reads_with_insertion: u64,
    /// Records carrying at least one deletion.
    pub reads_with_deletion: u64,
    /// Base composition in reference orientation, indexed by [`base_index`].
    pub base_counts: [u64; 5],
    /// Insert size histogram over positive `TLEN` values.
    pub insert_sizes: BTreeMap<u64, u64>,
    /// Per read position base composition, in reference orientation.
    pub nucleotide_by_position: Vec<[u64; 5]>,
    /// Per read position count of clipped bases.
    pub clipping_by_position: Vec<u64>,
    /// Total clipped bases, the denominator of the clipping profile.
    pub clipped_bases: u64,
    /// Homopolymer indel counts, indexed by [`base_index`], plus non-polymer.
    pub homopolymer_indels: [u64; 5],
    /// Indels not adjacent to a homopolymer run.
    pub non_polymer_indels: u64,
}

/// Index of a base in the fixed `A, C, G, T, N` ordering.
fn base_index(base: u8) -> usize {
    match base.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

/// The complement of a base, leaving anything unrecognised alone.
fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

impl QualimapCounters {
    /// Add another contig's counters.
    pub fn merge(&mut self, other: &QualimapCounters) {
        self.reads += other.reads;
        self.secondary += other.secondary;
        self.mapped += other.mapped;
        self.duplicates += other.duplicates;
        self.paired_first += other.paired_first;
        self.paired_second += other.paired_second;
        self.paired_both += other.paired_both;
        self.singletons += other.singletons;
        self.sequenced_bases += other.sequenced_bases;
        self.mapped_bases += other.mapped_bases;
        self.edit_distance += other.edit_distance;
        self.insertions += other.insertions;
        self.deletions += other.deletions;
        self.reads_with_insertion += other.reads_with_insertion;
        self.reads_with_deletion += other.reads_with_deletion;
        self.clipped_bases += other.clipped_bases;
        self.non_polymer_indels += other.non_polymer_indels;
        for (target, source) in self.base_counts.iter_mut().zip(&other.base_counts) {
            *target += source;
        }
        for (target, source) in self
            .homopolymer_indels
            .iter_mut()
            .zip(&other.homopolymer_indels)
        {
            *target += source;
        }
        for (size, count) in &other.insert_sizes {
            *self.insert_sizes.entry(*size).or_insert(0) += count;
        }
        if self.nucleotide_by_position.len() < other.nucleotide_by_position.len() {
            self.nucleotide_by_position
                .resize(other.nucleotide_by_position.len(), [0; 5]);
        }
        for (position, counts) in other.nucleotide_by_position.iter().enumerate() {
            for (target, source) in self.nucleotide_by_position[position].iter_mut().zip(counts) {
                *target += source;
            }
        }
        if self.clipping_by_position.len() < other.clipping_by_position.len() {
            self.clipping_by_position
                .resize(other.clipping_by_position.len(), 0);
        }
        for (position, count) in other.clipping_by_position.iter().enumerate() {
            self.clipping_by_position[position] += count;
        }
    }

    /// Mismatches, which Qualimap takes as `NM` less inserted bases only.
    pub fn mismatches(&self) -> u64 {
        self.edit_distance.saturating_sub(self.insertions)
    }

    /// Mismatches, insertions and deletions over sequenced bases.
    pub fn general_error_rate(&self) -> f64 {
        if self.sequenced_bases == 0 {
            return 0.0;
        }
        (self.mismatches() + self.insertions + self.deletions) as f64 / self.sequenced_bases as f64
    }

    /// Fraction of indels adjacent to a homopolymer run.
    pub fn homopolymer_fraction(&self) -> f64 {
        let poly: u64 = self.homopolymer_indels.iter().sum();
        let total = poly + self.non_polymer_indels;
        if total == 0 {
            0.0
        } else {
            poly as f64 / total as f64
        }
    }

    /// Mean, population standard deviation and median insert size.
    pub fn insert_size_stats(&self) -> (f64, f64, u64) {
        let n: u64 = self.insert_sizes.values().sum();
        if n == 0 {
            return (0.0, 0.0, 0);
        }
        let mean = self
            .insert_sizes
            .iter()
            .map(|(size, count)| *size as f64 * *count as f64)
            .sum::<f64>()
            / n as f64;
        let variance = self
            .insert_sizes
            .iter()
            .map(|(size, count)| {
                let diff = *size as f64 - mean;
                diff * diff * *count as f64
            })
            .sum::<f64>()
            / n as f64;
        let mut seen = 0u64;
        let mut median = 0u64;
        for (size, count) in &self.insert_sizes {
            seen += count;
            if seen > n / 2 {
                median = *size;
                break;
            }
        }
        (mean, variance.sqrt(), median)
    }
}

impl QualimapAccum {
    /// Prepare for one contig, splitting it into Qualimap's window grid.
    pub fn new(contig: &str, length: u64, num_windows: usize) -> Self {
        let window_size = length.div_ceil(num_windows as u64).max(1);
        let windows = length.div_ceil(window_size) as usize;
        Self {
            contig: contig.to_string(),
            length,
            window_size,
            coverage: vec![0; length as usize],
            mapq_sum: vec![0; length as usize],
            insert_window_sum: vec![0; windows],
            insert_window_count: vec![0; windows],
            counters: QualimapCounters::default(),
        }
    }

    /// Number of windows this contig is split into.
    pub fn window_count(&self) -> usize {
        self.insert_window_sum.len()
    }

    /// Width of each window; the last one may be shorter.
    pub fn window_size(&self) -> u64 {
        self.window_size
    }

    /// Offer one record.
    pub fn process_read(&mut self, record: &bam::Record) {
        let flags = record.flags();
        if flags & BAM_FSECONDARY != 0 {
            self.counters.secondary += 1;
            return;
        }
        self.counters.reads += 1;
        if flags & BAM_FUNMAP != 0 {
            return;
        }
        self.counters.mapped += 1;
        if flags & BAM_FDUP != 0 {
            self.counters.duplicates += 1;
        }

        if flags & BAM_FPAIRED != 0 {
            if flags & BAM_FMUNMAP != 0 {
                self.counters.singletons += 1;
            } else {
                self.counters.paired_both += 1;
                if flags & BAM_FREAD1 != 0 {
                    self.counters.paired_first += 1;
                }
                if flags & BAM_FREAD2 != 0 {
                    self.counters.paired_second += 1;
                }
            }
        }

        let mapq = u64::from(record.mapq());
        let sequence = record.seq().as_bytes();
        let reverse = flags & BAM_FREVERSE != 0;

        // Bases in reference orientation: reverse-complemented for a
        // reverse-strand read.
        let oriented: Vec<u8> = if reverse {
            sequence.iter().rev().map(|b| complement(*b)).collect()
        } else {
            sequence.clone()
        };

        let cigar = record.cigar();
        let ops: Vec<Cigar> = cigar.iter().copied().collect();

        // The clipped span is taken in sequencing orientation, unlike the
        // bases. That asymmetry is Qualimap's, and reproducing it is the only
        // way the composition figures agree.
        let leading_clip = match ops.first() {
            Some(Cigar::SoftClip(n)) | Some(Cigar::HardClip(n)) => *n as usize,
            _ => 0,
        };
        let trailing_clip = match ops.last() {
            Some(Cigar::SoftClip(n)) | Some(Cigar::HardClip(n)) => *n as usize,
            _ => 0,
        };

        let read_len = sequence.len();
        if self.counters.nucleotide_by_position.len() < read_len {
            self.counters
                .nucleotide_by_position
                .resize(read_len, [0; 5]);
            self.counters.clipping_by_position.resize(read_len, 0);
        }
        for position in 0..leading_clip.min(read_len) {
            self.counters.clipping_by_position[position] += 1;
            self.counters.clipped_bases += 1;
        }
        for offset in 0..trailing_clip.min(read_len) {
            let position = read_len - 1 - offset;
            self.counters.clipping_by_position[position] += 1;
            self.counters.clipped_bases += 1;
        }
        for position in leading_clip..read_len.saturating_sub(trailing_clip) {
            let base = oriented.get(position).copied().unwrap_or(b'N');
            self.counters.nucleotide_by_position[position][base_index(base)] += 1;
        }

        if let Ok(Aux::U8(nm)) = record.aux(b"NM") {
            self.counters.edit_distance += u64::from(nm);
        } else if let Ok(Aux::U16(nm)) = record.aux(b"NM") {
            self.counters.edit_distance += u64::from(nm);
        } else if let Ok(Aux::U32(nm)) = record.aux(b"NM") {
            self.counters.edit_distance += u64::from(nm);
        } else if let Ok(Aux::I32(nm)) = record.aux(b"NM") {
            self.counters.edit_distance += nm.max(0) as u64;
        }

        let mut reference_position = record.pos();
        let mut query_position = 0usize;
        let mut had_insertion = false;
        let mut had_deletion = false;

        for op in &ops {
            match op {
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                    let n = *n as usize;
                    for k in 0..n {
                        let position = reference_position + k as i64;
                        if position >= 0 && (position as usize) < self.coverage.len() {
                            self.coverage[position as usize] += 1;
                            self.mapq_sum[position as usize] += mapq;
                        }
                        let base = oriented.get(query_position + k).copied().unwrap_or(b'N');
                        self.counters.base_counts[base_index(base)] += 1;
                    }
                    self.counters.sequenced_bases += n as u64;
                    self.counters.mapped_bases += n as u64;
                    reference_position += n as i64;
                    query_position += n;
                }
                Cigar::Del(n) => {
                    let n = *n as usize;
                    for k in 0..n {
                        let position = reference_position + k as i64;
                        if position >= 0 && (position as usize) < self.coverage.len() {
                            self.coverage[position as usize] += 1;
                            self.mapq_sum[position as usize] += mapq;
                        }
                    }
                    self.counters.mapped_bases += n as u64;
                    self.counters.deletions += n as u64;
                    had_deletion = true;
                    self.classify_indel(&oriented, query_position);
                    reference_position += n as i64;
                }
                Cigar::Ins(n) => {
                    self.counters.insertions += u64::from(*n);
                    had_insertion = true;
                    self.classify_indel(&oriented, query_position);
                    query_position += *n as usize;
                }
                Cigar::RefSkip(n) => reference_position += i64::from(*n),
                Cigar::SoftClip(n) => query_position += *n as usize,
                Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }
        }
        if had_insertion {
            self.counters.reads_with_insertion += 1;
        }
        if had_deletion {
            self.counters.reads_with_deletion += 1;
        }

        let insert_size = record.insert_size();
        if insert_size > 0 {
            *self
                .counters
                .insert_sizes
                .entry(insert_size as u64)
                .or_insert(0) += 1;
            let window = (record.pos().max(0) as u64 / self.window_size) as usize;
            if window < self.insert_window_sum.len() {
                self.insert_window_sum[window] += insert_size;
                self.insert_window_count[window] += 1;
            }
        }
    }

    /// Charge an indel to a homopolymer bucket when the bases either side of
    /// it repeat, and to the non-polymer bucket otherwise.
    fn classify_indel(&mut self, oriented: &[u8], query_position: usize) {
        const RUN: usize = 4;
        let start = query_position.saturating_sub(RUN);
        let window = &oriented[start..query_position.min(oriented.len())];
        if window.len() == RUN && window.iter().all(|b| *b == window[0]) {
            self.counters.homopolymer_indels[base_index(window[0])] += 1;
        } else {
            self.counters.non_polymer_indels += 1;
        }
    }

    /// Consume the accumulator into its per-contig result.
    pub fn into_result(self) -> ContigQualimap {
        let window_size = self.window_size;
        let windows = self.insert_window_sum.len();
        let mut window_coverage = Vec::with_capacity(windows);
        let mut window_coverage_sd = Vec::with_capacity(windows);
        let mut window_mapq = Vec::with_capacity(windows);
        let mut window_insert = Vec::with_capacity(windows);
        let mut midpoints = Vec::with_capacity(windows);

        for window in 0..windows {
            let start = window as u64 * window_size;
            let end = ((window as u64 + 1) * window_size).min(self.length);
            let span = &self.coverage[start as usize..end as usize];
            let mapq_span = &self.mapq_sum[start as usize..end as usize];

            let mean = span.iter().map(|c| f64::from(*c)).sum::<f64>() / span.len() as f64;
            let variance = span
                .iter()
                .map(|c| {
                    let diff = f64::from(*c) - mean;
                    diff * diff
                })
                .sum::<f64>()
                / span.len() as f64;
            let covered: u64 = span.iter().map(|c| u64::from(*c)).sum();
            let mapq_total: u64 = mapq_span.iter().sum();

            window_coverage.push(mean);
            window_coverage_sd.push(variance.sqrt());
            window_mapq.push(if covered == 0 {
                0.0
            } else {
                mapq_total as f64 / covered as f64
            });
            window_insert.push(if self.insert_window_count[window] == 0 {
                0.0
            } else {
                self.insert_window_sum[window] as f64 / self.insert_window_count[window] as f64
            });
            midpoints.push((start + end + 1) as f64 / 2.0);
        }

        let mut coverage_histogram: BTreeMap<u32, u64> = BTreeMap::new();
        let mut mapq_histogram: BTreeMap<u32, u64> = BTreeMap::new();
        for (position, depth) in self.coverage.iter().enumerate() {
            *coverage_histogram.entry(*depth).or_insert(0) += 1;
            if *depth > 0 {
                // Truncated, not rounded: this is what Qualimap does.
                let mean = self.mapq_sum[position] / u64::from(*depth);
                *mapq_histogram.entry(mean as u32).or_insert(0) += 1;
            }
        }

        ContigQualimap {
            name: self.contig,
            length: self.length,
            coverage: self.coverage,
            window_size,
            midpoints,
            window_coverage,
            window_coverage_sd,
            window_mapq,
            window_insert,
            coverage_histogram,
            mapq_histogram,
            counters: self.counters,
        }
    }
}

/// One contig's Qualimap result.
#[derive(Debug, Clone)]
pub struct ContigQualimap {
    /// Contig name.
    pub name: String,
    /// Contig length.
    pub length: u64,
    /// Per-base coverage.
    pub coverage: Vec<u32>,
    /// Window width.
    pub window_size: u64,
    /// Window midpoints, as Qualimap reports positions.
    pub midpoints: Vec<f64>,
    /// Mean coverage per window.
    pub window_coverage: Vec<f64>,
    /// Coverage standard deviation per window.
    pub window_coverage_sd: Vec<f64>,
    /// Mean mapping quality per window, zero where uncovered.
    pub window_mapq: Vec<f64>,
    /// Mean insert size per window.
    pub window_insert: Vec<f64>,
    /// Bases at each exact coverage.
    pub coverage_histogram: BTreeMap<u32, u64>,
    /// Covered bases at each truncated mean mapping quality.
    pub mapq_histogram: BTreeMap<u32, u64>,
    /// Read-level counters gathered on this contig.
    pub counters: QualimapCounters,
}

impl ContigQualimap {
    /// Mean coverage over the contig.
    pub fn mean_coverage(&self) -> f64 {
        if self.length == 0 {
            0.0
        } else {
            self.coverage.iter().map(|c| f64::from(*c)).sum::<f64>() / self.length as f64
        }
    }

    /// Population standard deviation of per-base coverage.
    pub fn coverage_sd(&self) -> f64 {
        if self.length == 0 {
            return 0.0;
        }
        let mean = self.mean_coverage();
        let variance = self
            .coverage
            .iter()
            .map(|c| {
                let diff = f64::from(*c) - mean;
                diff * diff
            })
            .sum::<f64>()
            / self.length as f64;
        variance.sqrt()
    }

    /// Mean of the per-window mapping qualities, uncovered windows included.
    pub fn mean_mapping_quality(&self) -> f64 {
        if self.window_mapq.is_empty() {
            0.0
        } else {
            self.window_mapq.iter().sum::<f64>() / self.window_mapq.len() as f64
        }
    }

    /// Percentage of the contig at or above each coverage level.
    pub fn genome_fraction(&self) -> Vec<(u32, f64)> {
        (1..=MAX_FRACTION_LEVEL)
            .map(|level| {
                let at_or_above = self.coverage.iter().filter(|c| **c >= level).count();
                (level, 100.0 * at_or_above as f64 / self.length as f64)
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn window_grid_matches_qualimaps_arithmetic() {
        // 40001 bases into 400 windows: 101 bases each, and 397 of them.
        let accum = QualimapAccum::new("chr22", 40001, DEFAULT_NUM_WINDOWS);
        assert_eq!(accum.window_size(), 101);
        assert_eq!(accum.window_count(), 397);
    }

    #[test]
    fn a_short_contig_still_gets_one_window() {
        let accum = QualimapAccum::new("small", 10, DEFAULT_NUM_WINDOWS);
        assert_eq!(accum.window_size(), 1);
        assert_eq!(accum.window_count(), 10);
    }

    #[test]
    fn mismatches_subtract_insertions_but_not_deletions() {
        let mut counters = QualimapCounters {
            edit_distance: 1352,
            insertions: 2,
            deletions: 10,
            ..Default::default()
        };
        assert_eq!(counters.mismatches(), 1350, "deletions are not subtracted");
        counters.deletions = 0;
        assert_eq!(counters.mismatches(), 1350);
    }

    #[test]
    fn insert_size_statistics_use_the_population_denominator() {
        let mut counters = QualimapCounters::default();
        for size in [1u64, 2, 3] {
            counters.insert_sizes.insert(size, 1);
        }
        let (mean, sd, median) = counters.insert_size_stats();
        assert!((mean - 2.0).abs() < 1e-12);
        // Population variance of 1, 2, 3 is 2/3.
        assert!((sd - (2.0f64 / 3.0).sqrt()).abs() < 1e-12, "got {sd}");
        assert_eq!(median, 2);
    }

    #[test]
    fn base_indexing_folds_anything_unknown_into_n() {
        assert_eq!(base_index(b'A'), 0);
        assert_eq!(base_index(b'c'), 1);
        assert_eq!(base_index(b'N'), 4);
        assert_eq!(base_index(b'R'), 4, "ambiguity codes are counted as N");
    }

    #[test]
    fn complement_leaves_unknown_bases_alone() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'g'), b'C');
        assert_eq!(complement(b'N'), b'N');
        assert_eq!(complement(b'R'), b'R');
    }
}
