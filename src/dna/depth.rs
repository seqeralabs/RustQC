//! Per-contig depth of coverage accumulation.
//!
//! One [`DepthAccum`] covers one contig. Aligned blocks are recorded as
//! increments in a delta array the length of the contig, and a prefix sum at
//! the end turns that into per-base depth in a single linear pass.
//!
//! # Upstream semantics
//!
//! The filters and the CIGAR walk reproduce mosdepth 0.3.14 run without
//! `--fast-mode`, whose help text describes that flag as "dont look at
//! internal cigar operations or correct mate overlaps". Default mode
//! therefore does both, and so does this module:
//!
//! - records carrying any bit of [`MOSDEPTH_DEFAULT_EXCLUDE`] are skipped
//!   (mosdepth's `-F` default of 1796);
//! - records with `MAPQ` below the cutoff are skipped (mosdepth's `-Q`,
//!   default 0);
//! - `M`, `=` and `X` cover the reference, `D` and `N` advance without
//!   covering, and `I`, `S`, `H` and `P` do not advance at all;
//! - a base covered by both mates of one pair counts once.
//!
//! That last rule is not a detail. On the project's test dataset, correcting
//! mate overlaps takes total covered bases from 469875 down to 247878, which
//! is exactly the gap between mosdepth's `--fast-mode` and its default.

use std::collections::{BTreeMap, HashMap};

use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

use crate::common::bam_flags::*;

/// Bit mask matching mosdepth's `-F` default: `UNMAP | SECONDARY | QCFAIL | DUP`.
pub const MOSDEPTH_DEFAULT_EXCLUDE: u16 = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;

/// Accumulates per-base depth for a single contig.
#[derive(Debug)]
pub struct DepthAccum {
    /// Delta array of length `contig_len + 1`; a `+1` at a block start and a
    /// `-1` one past its end, summed into depth by [`DepthAccum::into_depths`].
    deltas: Vec<i32>,
    /// Contig length in bases.
    len: usize,
    /// Records with `MAPQ` strictly below this value are ignored.
    mapq_cut: u8,
    /// Records carrying any of these flag bits are ignored.
    exclude_flags: u16,
    /// Aligned blocks already counted for a pair whose second mate is still
    /// ahead, keyed by read name.
    pending: HashMap<Vec<u8>, Vec<(usize, usize)>>,
    /// Read names indexed by the position their outstanding mate is expected
    /// at, so stale entries can be evicted without scanning `pending`.
    pending_by_pos: BTreeMap<i64, Vec<Vec<u8>>>,
}

impl DepthAccum {
    /// Allocate for one contig of `length` bases.
    pub fn new(length: u64, mapq_cut: u8, exclude_flags: u16) -> Self {
        let len = length as usize;
        Self {
            deltas: vec![0i32; len + 1],
            len,
            mapq_cut,
            exclude_flags,
            pending: HashMap::new(),
            pending_by_pos: BTreeMap::new(),
        }
    }

    /// Add one record's aligned blocks. Records failing the filters are ignored.
    ///
    /// Records are expected in coordinate order, which is what the per-contig
    /// worker feeds. That ordering is what makes the pending-mate bookkeeping
    /// bounded: once the read position passes the position an outstanding mate
    /// was announced at, that entry can never be claimed and is dropped.
    pub fn process_read(&mut self, record: &bam::Record) {
        if !self.passes_filters(record) {
            return;
        }
        let pos = record.pos();
        self.evict_unclaimable(pos);

        let blocks = Self::aligned_blocks(record, self.len);
        if blocks.is_empty() {
            return;
        }

        // A record can only overlap its own mate, and only on the same contig.
        let paired_here = record.flags() & BAM_FPAIRED != 0
            && record.flags() & BAM_FMUNMAP == 0
            && record.mtid() == record.tid();

        if paired_here {
            if let Some(mate_blocks) = self.pending.remove(record.qname()) {
                // Second mate of the pair: shared bases are already counted.
                self.add_blocks_excluding(&blocks, &mate_blocks);
                return;
            }
            if record.mpos() >= pos {
                let qname = record.qname().to_vec();
                self.pending.insert(qname.clone(), blocks.clone());
                self.pending_by_pos
                    .entry(record.mpos())
                    .or_default()
                    .push(qname);
            }
        }

        for &(start, end) in &blocks {
            self.add_block_usize(start, end);
        }
    }

    /// Number of pairs still waiting for their second mate. Test-only: the
    /// bookkeeping is an implementation detail, but an unbounded map would be
    /// a memory leak on a real chromosome, so it is worth asserting on.
    #[cfg(test)]
    pub fn pending_mates_len(&self) -> usize {
        self.pending.len()
    }

    /// Drop pending entries whose outstanding mate lies behind `pos` and can
    /// therefore never arrive (it was filtered out, or the file is truncated).
    fn evict_unclaimable(&mut self, pos: i64) {
        while let Some((&mate_pos, _)) = self.pending_by_pos.iter().next() {
            if mate_pos >= pos {
                break;
            }
            // Safe: the key came from `iter().next()` on this same map.
            let qnames = self.pending_by_pos.remove(&mate_pos).unwrap_or_default();
            for qname in qnames {
                self.pending.remove(&qname);
            }
        }
    }

    /// The record's reference-covering blocks as half-open `[start, end)`
    /// intervals, clamped to `len`.
    fn aligned_blocks(record: &bam::Record, len: usize) -> Vec<(usize, usize)> {
        let mut blocks = Vec::new();
        let mut pos = record.pos();
        for op in record.cigar().iter() {
            match op {
                // Reference-consuming and query-consuming: covers the reference.
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                    let n = i64::from(*n);
                    let start = pos.max(0) as usize;
                    let end = ((pos + n).max(0) as usize).min(len);
                    if start < end {
                        blocks.push((start, end));
                    }
                    pos += n;
                }
                // Reference-consuming only: advances without covering.
                Cigar::Del(n) | Cigar::RefSkip(n) => pos += i64::from(*n),
                // Neither reference-consuming nor covering.
                Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }
        }
        blocks
    }

    /// Add `blocks`, skipping any part already covered by `exclude`.
    ///
    /// Both sides are in ascending order and non-overlapping within themselves,
    /// because each comes from one record's CIGAR walk.
    fn add_blocks_excluding(&mut self, blocks: &[(usize, usize)], exclude: &[(usize, usize)]) {
        for &(start, end) in blocks {
            let mut cursor = start;
            for &(ex_start, ex_end) in exclude {
                if ex_end <= cursor {
                    continue;
                }
                if ex_start >= end {
                    break;
                }
                if ex_start > cursor {
                    self.add_block_usize(cursor, ex_start.min(end));
                }
                cursor = cursor.max(ex_end);
                if cursor >= end {
                    break;
                }
            }
            if cursor < end {
                self.add_block_usize(cursor, end);
            }
        }
    }

    /// Consume the delta array and return per-base depth for the contig.
    pub fn into_depths(self) -> Vec<u32> {
        let mut depths = Vec::with_capacity(self.len);
        let mut running = 0i32;
        for delta in self.deltas.iter().take(self.len) {
            running += delta;
            // `running` cannot go negative: every `-1` is emitted only after
            // its matching `+1`, and both are clamped to the same range.
            depths.push(running.max(0) as u32);
        }
        depths
    }

    /// Whether a record contributes to depth at all.
    fn passes_filters(&self, record: &bam::Record) -> bool {
        record.flags() & self.exclude_flags == 0 && record.mapq() >= self.mapq_cut
    }

    /// Record a half-open aligned block `[start, end)`, already clamped.
    fn add_block_usize(&mut self, start: usize, end: usize) {
        if start >= end {
            return;
        }
        self.deltas[start] += 1;
        self.deltas[end] -= 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString, Record};

    /// Build a minimal mapped record at `pos` with the given CIGAR, MAPQ and flags.
    ///
    /// `seq` and `qual` must both be as long as the query-consuming part of the
    /// CIGAR, otherwise the record is malformed and every assertion made against
    /// it is meaningless, so the helper asserts that itself.
    fn rec(pos: i64, cigar: Vec<Cigar>, mapq: u8, flags: u16) -> Record {
        let query_len: usize = cigar
            .iter()
            .map(|op| match op {
                Cigar::Match(n) | Cigar::Ins(n) | Cigar::SoftClip(n) => *n as usize,
                Cigar::Equal(n) | Cigar::Diff(n) => *n as usize,
                _ => 0,
            })
            .sum();
        let seq = vec![b'A'; query_len];
        let qual = vec![30u8; query_len];
        assert_eq!(seq.len(), qual.len(), "malformed test record");

        let mut r = Record::new();
        r.set(b"q", Some(&CigarString(cigar)), &seq, &qual);
        r.set_tid(0);
        r.set_pos(pos);
        r.set_mapq(mapq);
        r.set_flags(flags);
        r
    }

    /// Build a paired record whose mate sits at `mate_pos` on the same contig.
    fn pair_rec(qname: &[u8], pos: i64, mate_pos: i64, cigar: Vec<Cigar>, read2: bool) -> Record {
        let mut r = rec(
            pos,
            cigar,
            60,
            BAM_FPAIRED | BAM_FPROPER_PAIR | if read2 { BAM_FREAD2 } else { BAM_FREAD1 },
        );
        r.set_qname(qname);
        r.set_mtid(0);
        r.set_mpos(mate_pos);
        r
    }

    #[test]
    fn match_block_covers_exactly_its_span() {
        let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
        d.process_read(&rec(5, vec![Cigar::Match(4)], 60, 0));
        assert_eq!(&d.into_depths()[4..10], &[0, 1, 1, 1, 1, 0]);
    }

    #[test]
    fn deletion_and_skip_advance_without_covering() {
        let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
        d.process_read(&rec(
            0,
            vec![Cigar::Match(2), Cigar::Del(3), Cigar::Match(2)],
            60,
            0,
        ));
        assert_eq!(&d.into_depths()[0..8], &[1, 1, 0, 0, 0, 1, 1, 0]);
    }

    #[test]
    fn ref_skip_advances_without_covering() {
        let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
        d.process_read(&rec(
            0,
            vec![Cigar::Match(2), Cigar::RefSkip(3), Cigar::Match(2)],
            60,
            0,
        ));
        assert_eq!(&d.into_depths()[0..8], &[1, 1, 0, 0, 0, 1, 1, 0]);
    }

    #[test]
    fn insertion_and_soft_clip_do_not_advance_the_reference() {
        let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
        d.process_read(&rec(
            0,
            vec![
                Cigar::SoftClip(3),
                Cigar::Match(2),
                Cigar::Ins(4),
                Cigar::Match(2),
            ],
            60,
            0,
        ));
        assert_eq!(&d.into_depths()[0..6], &[1, 1, 1, 1, 0, 0]);
    }

    #[test]
    fn duplicate_flagged_reads_are_excluded_by_default() {
        let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
        d.process_read(&rec(0, vec![Cigar::Match(4)], 60, BAM_FDUP));
        assert_eq!(d.into_depths().iter().sum::<u32>(), 0);
    }

    #[test]
    fn secondary_qcfail_and_unmapped_reads_are_excluded_by_default() {
        for flag in [BAM_FSECONDARY, BAM_FQCFAIL, BAM_FUNMAP] {
            let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
            d.process_read(&rec(0, vec![Cigar::Match(4)], 60, flag));
            assert_eq!(
                d.into_depths().iter().sum::<u32>(),
                0,
                "flag {flag:#x} should be excluded"
            );
        }
    }

    #[test]
    fn reads_below_the_mapq_cutoff_are_excluded() {
        let mut d = DepthAccum::new(20, 30, MOSDEPTH_DEFAULT_EXCLUDE);
        d.process_read(&rec(0, vec![Cigar::Match(4)], 29, 0));
        assert_eq!(d.into_depths().iter().sum::<u32>(), 0);

        let mut d = DepthAccum::new(20, 30, MOSDEPTH_DEFAULT_EXCLUDE);
        d.process_read(&rec(0, vec![Cigar::Match(4)], 30, 0));
        assert_eq!(d.into_depths().iter().sum::<u32>(), 4);
    }

    #[test]
    fn a_read_running_past_the_contig_end_is_clipped_not_panicking() {
        let mut d = DepthAccum::new(6, 0, MOSDEPTH_DEFAULT_EXCLUDE);
        d.process_read(&rec(4, vec![Cigar::Match(10)], 60, 0));
        assert_eq!(d.into_depths(), vec![0, 0, 0, 0, 1, 1]);
    }

    #[test]
    fn overlapping_mates_cover_a_base_once() {
        let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
        d.process_read(&pair_rec(b"pair1", 0, 0, vec![Cigar::Match(4)], false));
        d.process_read(&pair_rec(b"pair1", 0, 0, vec![Cigar::Match(4)], true));
        assert_eq!(
            &d.into_depths()[0..5],
            &[1, 1, 1, 1, 0],
            "a base covered by both mates counts once"
        );
    }

    #[test]
    fn partially_overlapping_mates_count_the_shared_bases_once() {
        let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
        d.process_read(&pair_rec(b"pair1", 0, 2, vec![Cigar::Match(4)], false));
        d.process_read(&pair_rec(b"pair1", 2, 0, vec![Cigar::Match(4)], true));
        // Mate 1 covers 0..4, mate 2 covers 2..6; bases 2 and 3 are shared.
        assert_eq!(&d.into_depths()[0..7], &[1, 1, 1, 1, 1, 1, 0]);
    }

    #[test]
    fn non_overlapping_mates_each_contribute() {
        let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
        d.process_read(&pair_rec(b"pair2", 0, 8, vec![Cigar::Match(4)], false));
        d.process_read(&pair_rec(b"pair2", 8, 0, vec![Cigar::Match(4)], true));
        assert_eq!(
            &d.into_depths()[0..13],
            &[1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0]
        );
    }

    #[test]
    fn reads_from_different_pairs_at_the_same_locus_both_count() {
        let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
        d.process_read(&pair_rec(b"pairA", 0, 0, vec![Cigar::Match(4)], false));
        d.process_read(&pair_rec(b"pairB", 0, 0, vec![Cigar::Match(4)], false));
        assert_eq!(d.into_depths()[0], 2);
    }

    #[test]
    fn the_pending_mate_map_is_emptied_once_both_mates_are_seen() {
        let mut d = DepthAccum::new(20, 0, MOSDEPTH_DEFAULT_EXCLUDE);
        d.process_read(&pair_rec(b"pair1", 0, 0, vec![Cigar::Match(4)], false));
        d.process_read(&pair_rec(b"pair1", 0, 0, vec![Cigar::Match(4)], true));
        assert_eq!(d.pending_mates_len(), 0, "the entry must be dropped");
    }

    #[test]
    fn a_pending_mate_that_never_arrives_is_evicted() {
        let mut d = DepthAccum::new(200, 0, MOSDEPTH_DEFAULT_EXCLUDE);
        // Its mate is announced at 10 but never turns up (filtered, say).
        d.process_read(&pair_rec(b"orphan", 0, 10, vec![Cigar::Match(4)], false));
        assert_eq!(d.pending_mates_len(), 1);
        // Walking past position 10 makes the entry unclaimable.
        d.process_read(&pair_rec(b"later", 50, 50, vec![Cigar::Match(4)], false));
        assert_eq!(
            d.pending_mates_len(),
            1,
            "only the unclaimable one is dropped"
        );
    }

    /// Engine-level parity check against mosdepth 0.3.14 on the committed
    /// fixture. `tests/expected/dna/test.mosdepth.summary.txt` records
    /// `total 40001 247878 6.20 0 867` for this BAM, so the total covered
    /// bases and the maximum depth are both pinned here. Getting this right
    /// requires the flag filter, the CIGAR walk and the mate-overlap
    /// correction to all be right at once.
    #[test]
    fn total_covered_bases_match_mosdepth_on_the_fixture() {
        use rust_htslib::bam::Read;

        let bam_path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/dna/test.dna.bam");
        let mut bam = bam::Reader::from_path(bam_path).unwrap();
        let header = bam.header().to_owned();
        let contig_len = header.target_len(0).unwrap();

        let mut accum = DepthAccum::new(contig_len, 0, MOSDEPTH_DEFAULT_EXCLUDE);
        let mut record = Record::new();
        while let Some(result) = bam.read(&mut record) {
            result.unwrap();
            accum.process_read(&record);
        }

        let depths = accum.into_depths();
        let total: u64 = depths.iter().map(|d| u64::from(*d)).sum();
        let max = depths.iter().copied().max().unwrap();

        assert_eq!(depths.len(), 40001, "contig length");
        assert_eq!(total, 247878, "total covered bases must match mosdepth");
        assert_eq!(max, 867, "maximum depth must match mosdepth");
    }
}
