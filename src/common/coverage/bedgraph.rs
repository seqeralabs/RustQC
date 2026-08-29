//! Per-base coverage in `bedtools genomecov` semantics.
//!
//! # Upstream semantics
//!
//! Matches `bedtools genomecov -ibam -bg -split`. Two things about it are
//! worth stating, because neither matches the depth engines already in this
//! crate:
//!
//! `-split` counts only the reference blocks a read actually aligns to, so an
//! `N` in the CIGAR breaks the read into separate intervals rather than
//! covering the intron. Deletions break it too.
//!
//! There is no filtering. Duplicates, secondary alignments and low mapping
//! quality all contribute, and overlapping mates of a pair each count. That is
//! the opposite of what mosdepth does by default and of what
//! `CollectWgsMetrics` does, and it is why this has its own accumulator rather
//! than reusing either.

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

use crate::common::bam_flags::*;

/// One bedGraph interval: a half-open span at a constant depth.
#[derive(Debug, Clone, PartialEq)]
pub struct Interval {
    /// Contig name.
    pub chrom: String,
    /// Zero-based start.
    pub start: u32,
    /// Half-open end.
    pub end: u32,
    /// Depth over the span, after any scaling.
    pub value: f32,
}

/// Accumulates per-base coverage for one contig.
#[derive(Debug)]
pub struct CoverageAccum {
    /// Delta array of length `contig_len + 1`.
    deltas: Vec<i32>,
    len: usize,
    /// Only reads on this strand contribute; `None` counts every read.
    strand: Option<char>,
    /// Whether the array grows to fit the reads rather than being clamped to a
    /// known contig length. This is a field rather than something inferred
    /// from the array, because after its first read a growable accumulator is
    /// indistinguishable from a fixed one of that size and would silently stop
    /// growing.
    growable: bool,
}

impl CoverageAccum {
    /// Allocate for one contig of known length.
    ///
    /// `strand` restricts the reads counted, which is how the forward and
    /// reverse tracks are produced. It follows `bedtools genomecov -strand`:
    /// the read's own strand, not the fragment's.
    pub fn new(length: u64, strand: Option<char>) -> Self {
        let len = length as usize;
        Self {
            deltas: vec![0i32; len + 1],
            len,
            strand,
            growable: false,
        }
    }

    /// Allocate without knowing the contig length, growing as reads arrive.
    ///
    /// This is what lets the accumulator sit in the single alignment pass
    /// alongside the others, which are constructed before any contig is known.
    /// The length is only needed when the track is written, and it comes from
    /// the alignment header there.
    pub fn growable(strand: Option<char>) -> Self {
        Self {
            deltas: Vec::new(),
            len: 0,
            strand,
            growable: true,
        }
    }

    /// Offer one record.
    pub fn process_read(&mut self, record: &bam::Record) {
        if record.flags() & BAM_FUNMAP != 0 {
            return;
        }
        if let Some(wanted) = self.strand {
            let reverse = record.flags() & BAM_FREVERSE != 0;
            let actual = if reverse { '-' } else { '+' };
            if actual != wanted {
                return;
            }
        }

        let mut position = record.pos();
        for op in record.cigar().iter() {
            match op {
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                    let n = i64::from(*n);
                    self.add(position, position + n);
                    position += n;
                }
                // Both break the read into separate intervals under `-split`.
                Cigar::Del(n) | Cigar::RefSkip(n) => position += i64::from(*n),
                Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }
        }
    }

    /// Record a half-open block, clamped to the contig when its length is
    /// known and growing the array when it is not.
    fn add(&mut self, start: i64, end: i64) {
        let start = start.max(0) as usize;
        let mut end = end.max(0) as usize;
        if self.growable {
            if end + 1 > self.deltas.len() {
                self.deltas.resize(end + 1, 0);
                self.len = end;
            }
        } else {
            end = end.min(self.len);
        }
        if start >= end {
            return;
        }
        self.deltas[start] += 1;
        self.deltas[end] -= 1;
    }

    /// Fold another accumulator for the same contig in.
    pub fn merge(&mut self, other: CoverageAccum) {
        if other.deltas.len() > self.deltas.len() {
            self.deltas.resize(other.deltas.len(), 0);
            self.len = other.len.max(self.len);
        }
        for (index, delta) in other.deltas.iter().enumerate() {
            self.deltas[index] += delta;
        }
    }

    /// Collapse into bedGraph intervals, dropping the zero-depth spans that
    /// `bedtools genomecov -bg` omits.
    ///
    /// `scale` multiplies every depth, which is how normalised tracks are
    /// produced; a scale of 1.0 leaves the raw counts.
    pub fn into_intervals(self, chrom: &str, scale: f32) -> Vec<Interval> {
        let mut intervals = Vec::new();
        let mut running = 0i32;
        let mut run_start = 0usize;
        let mut run_depth = 0i32;

        for position in 0..self.len {
            running += self.deltas[position];
            if position == 0 {
                run_depth = running;
                run_start = 0;
                continue;
            }
            if running != run_depth {
                if run_depth > 0 {
                    intervals.push(Interval {
                        chrom: chrom.to_string(),
                        start: run_start as u32,
                        end: position as u32,
                        value: run_depth as f32 * scale,
                    });
                }
                run_depth = running;
                run_start = position;
            }
        }
        if run_depth > 0 && run_start < self.len {
            intervals.push(Interval {
                chrom: chrom.to_string(),
                start: run_start as u32,
                end: self.len as u32,
                value: run_depth as f32 * scale,
            });
        }
        intervals
    }
}

/// Write intervals as a bedGraph.
pub fn write_bedgraph(intervals: &[Interval], path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .map(std::io::BufWriter::new)
        .with_context(|| format!("Failed to create bedGraph: {}", path.display()))?;
    for interval in intervals {
        // Whole numbers print without a fractional part, as bedtools does for
        // unscaled counts; a scaled track keeps its decimals.
        if interval.value.fract() == 0.0 {
            writeln!(
                out,
                "{}\t{}\t{}\t{}",
                interval.chrom, interval.start, interval.end, interval.value as i64
            )?;
        } else {
            writeln!(
                out,
                "{}\t{}\t{}\t{}",
                interval.chrom, interval.start, interval.end, interval.value
            )?;
        }
    }
    out.flush()?;
    Ok(())
}

/// Coverage tracks across every contig, one accumulator per strand.
///
/// Positions are contig-local, so accumulators from different contigs must
/// never be summed. This keeps them keyed by contig, which is what makes the
/// per-chromosome workers safe to merge.
#[derive(Debug, Default)]
pub struct CoverageTracks {
    /// Strands to track, in output order. A single `None` means one combined
    /// track counting every read.
    strands: Vec<Option<char>>,
    /// Per contig, one accumulator per entry in `strands`.
    per_chrom: std::collections::HashMap<String, Vec<CoverageAccum>>,
}

impl CoverageTracks {
    /// Track the given strands. An empty list disables the tracks entirely.
    pub fn new(strands: Vec<Option<char>>) -> Self {
        Self {
            strands,
            per_chrom: std::collections::HashMap::new(),
        }
    }

    /// Whether anything is being tracked.
    pub fn is_enabled(&self) -> bool {
        !self.strands.is_empty()
    }

    /// Offer one record, which must belong to `chrom`.
    pub fn process_read(&mut self, record: &bam::Record, chrom: &str) {
        if self.strands.is_empty() {
            return;
        }
        let strands = &self.strands;
        let accums = self.per_chrom.entry(chrom.to_string()).or_insert_with(|| {
            strands
                .iter()
                .map(|strand| CoverageAccum::growable(*strand))
                .collect()
        });
        for accum in accums.iter_mut() {
            accum.process_read(record);
        }
    }

    /// Fold another set in, contig by contig.
    pub fn merge(&mut self, other: CoverageTracks) {
        for (chrom, theirs) in other.per_chrom {
            match self.per_chrom.get_mut(&chrom) {
                Some(mine) => {
                    for (mine, theirs) in mine.iter_mut().zip(theirs) {
                        mine.merge(theirs);
                    }
                }
                None => {
                    self.per_chrom.insert(chrom, theirs);
                }
            }
        }
    }

    /// Collapse into one interval list per tracked strand, contigs in the
    /// order given, which is the alignment header's order.
    pub fn into_intervals(
        mut self,
        chrom_order: &[String],
        scale: f32,
    ) -> Vec<(Option<char>, Vec<Interval>)> {
        let mut out: Vec<(Option<char>, Vec<Interval>)> =
            self.strands.iter().map(|s| (*s, Vec::new())).collect();
        for chrom in chrom_order {
            if let Some(accums) = self.per_chrom.remove(chrom) {
                for (index, accum) in accums.into_iter().enumerate() {
                    out[index].1.extend(accum.into_intervals(chrom, scale));
                }
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{CigarString, Record};

    fn record(pos: i64, cigar: Vec<Cigar>, flags: u16) -> Record {
        let query: usize = cigar
            .iter()
            .map(|op| match op {
                Cigar::Match(n) | Cigar::Ins(n) | Cigar::SoftClip(n) => *n as usize,
                Cigar::Equal(n) | Cigar::Diff(n) => *n as usize,
                _ => 0,
            })
            .sum();
        let mut r = Record::new();
        r.set(
            b"q",
            Some(&CigarString(cigar)),
            &vec![b'A'; query],
            &vec![30u8; query],
        );
        r.set_pos(pos);
        r.set_flags(flags);
        r
    }

    #[test]
    fn a_spliced_read_does_not_cover_its_intron() {
        let mut accum = CoverageAccum::new(20, None);
        accum.process_read(&record(
            0,
            vec![Cigar::Match(3), Cigar::RefSkip(5), Cigar::Match(3)],
            0,
        ));
        let intervals = accum.into_intervals("chr1", 1.0);
        assert_eq!(
            intervals,
            vec![
                Interval {
                    chrom: "chr1".into(),
                    start: 0,
                    end: 3,
                    value: 1.0
                },
                Interval {
                    chrom: "chr1".into(),
                    start: 8,
                    end: 11,
                    value: 1.0
                },
            ],
            "the intron is absent rather than at depth zero"
        );
    }

    #[test]
    fn a_deletion_also_splits_the_interval() {
        let mut accum = CoverageAccum::new(20, None);
        accum.process_read(&record(
            0,
            vec![Cigar::Match(2), Cigar::Del(2), Cigar::Match(2)],
            0,
        ));
        let intervals = accum.into_intervals("chr1", 1.0);
        assert_eq!(intervals.len(), 2);
        assert_eq!((intervals[0].start, intervals[0].end), (0, 2));
        assert_eq!((intervals[1].start, intervals[1].end), (4, 6));
    }

    #[test]
    fn runs_of_equal_depth_collapse_into_one_interval() {
        let mut accum = CoverageAccum::new(20, None);
        accum.process_read(&record(0, vec![Cigar::Match(10)], 0));
        accum.process_read(&record(0, vec![Cigar::Match(10)], 0));
        let intervals = accum.into_intervals("chr1", 1.0);
        assert_eq!(intervals.len(), 1, "one span at depth 2");
        assert_eq!(intervals[0].value, 2.0);
    }

    #[test]
    fn zero_depth_spans_are_omitted() {
        let mut accum = CoverageAccum::new(20, None);
        accum.process_read(&record(5, vec![Cigar::Match(3)], 0));
        let intervals = accum.into_intervals("chr1", 1.0);
        assert_eq!(intervals.len(), 1);
        assert_eq!((intervals[0].start, intervals[0].end), (5, 8));
    }

    #[test]
    fn nothing_is_filtered_out() {
        // Duplicates and secondary alignments contribute, unlike every other
        // depth engine in this crate.
        let mut accum = CoverageAccum::new(20, None);
        accum.process_read(&record(0, vec![Cigar::Match(4)], BAM_FDUP));
        accum.process_read(&record(0, vec![Cigar::Match(4)], BAM_FSECONDARY));
        assert_eq!(accum.into_intervals("chr1", 1.0)[0].value, 2.0);
    }

    #[test]
    fn an_unmapped_read_contributes_nothing() {
        let mut accum = CoverageAccum::new(20, None);
        accum.process_read(&record(0, vec![Cigar::Match(4)], BAM_FUNMAP));
        assert!(accum.into_intervals("chr1", 1.0).is_empty());
    }

    #[test]
    fn the_strand_filter_follows_the_reads_own_strand() {
        let mut forward = CoverageAccum::new(20, Some('+'));
        forward.process_read(&record(0, vec![Cigar::Match(4)], 0));
        forward.process_read(&record(0, vec![Cigar::Match(4)], BAM_FREVERSE));
        assert_eq!(forward.into_intervals("chr1", 1.0)[0].value, 1.0);

        let mut reverse = CoverageAccum::new(20, Some('-'));
        reverse.process_read(&record(0, vec![Cigar::Match(4)], 0));
        reverse.process_read(&record(0, vec![Cigar::Match(4)], BAM_FREVERSE));
        assert_eq!(reverse.into_intervals("chr1", 1.0)[0].value, 1.0);
    }

    #[test]
    fn scaling_multiplies_every_depth() {
        let mut accum = CoverageAccum::new(20, None);
        accum.process_read(&record(0, vec![Cigar::Match(4)], 0));
        assert_eq!(accum.into_intervals("chr1", 2.5)[0].value, 2.5);
    }

    #[test]
    fn a_growable_accumulator_needs_no_length_up_front() {
        let mut accum = CoverageAccum::growable(None);
        accum.process_read(&record(100, vec![Cigar::Match(5)], 0));
        let intervals = accum.into_intervals("chr1", 1.0);
        assert_eq!(intervals.len(), 1);
        assert_eq!((intervals[0].start, intervals[0].end), (100, 105));
    }

    #[test]
    fn merging_adds_two_workers_coverage() {
        let mut a = CoverageAccum::growable(None);
        a.process_read(&record(0, vec![Cigar::Match(4)], 0));
        let mut b = CoverageAccum::growable(None);
        b.process_read(&record(0, vec![Cigar::Match(4)], 0));
        b.process_read(&record(10, vec![Cigar::Match(2)], 0));
        a.merge(b);
        let intervals = a.into_intervals("chr1", 1.0);
        assert_eq!(intervals[0].value, 2.0, "the shared span doubles");
        assert_eq!((intervals[1].start, intervals[1].end), (10, 12));
    }

    #[test]
    fn a_read_running_past_the_contig_end_is_clipped() {
        let mut accum = CoverageAccum::new(6, None);
        accum.process_read(&record(4, vec![Cigar::Match(10)], 0));
        let intervals = accum.into_intervals("chr1", 1.0);
        assert_eq!((intervals[0].start, intervals[0].end), (4, 6));
    }
}
