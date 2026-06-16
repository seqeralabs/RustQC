//! Genome coverage accumulation matching `bedtools genomecov`.
//!
//! Replicates nf-core/rnaseq bigWig generation parameters:
//! - Combined track: `-split -bg`
//! - Per-strand tracks (stranded libraries): `-split -du -strand +/- -bg`
//!
//! # Streaming, memory-bounded design
//!
//! A naive difference array allocates one counter per genome base, which for a
//! mammalian genome is tens of GB per track (and the parallel path keeps every
//! chromosome resident at once). Instead, because the input BAM is
//! coordinate-sorted, each chromosome's reads arrive in ascending start order.
//! We keep only a small sparse map of pending coverage deltas and a sweep
//! cursor: once the cursor passes a position, no future read can change
//! coverage there, so we emit the finished bedGraph intervals and drop the
//! deltas below the cursor. Peak memory is therefore bounded by the widest
//! single-read genomic span (introns included), not by chromosome length or
//! sequencing depth.
//!
//! The emitted bedGraph is bit-for-bit identical to the equivalent dense
//! difference-array sweep (see the brute-force equivalence tests below), which
//! in turn matches `bedtools genomecov`.

use std::collections::BTreeMap;

use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

use crate::rna::bam_flags::BAM_FREAD2;
use crate::Strandedness;

/// Strand filter for per-strand coverage tracks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrandFilter {
    /// Forward reference strand (`-strand +`).
    Plus,
    /// Reverse reference strand (`-strand -`).
    Minus,
}

/// Settings for a single coverage track.
///
/// All nf-core/rnaseq tracks use `-split` (reads are always split at splicing
/// blocks), so splitting is unconditional rather than a toggle.
#[derive(Debug, Clone, Copy, Default)]
struct TrackSettings {
    /// Flip second-mate strand for strand-specific libraries (`-du`).
    deduplicate_strand: bool,
    /// Optional strand filter (`-strand`).
    strand: Option<StrandFilter>,
}

/// A single bedGraph interval: `(chrom, start, end, depth)` with `depth > 0`.
type BedGraphEntry = (String, u32, u32, u32);

/// Streaming coverage accumulator for one track (combined, forward, or reverse).
///
/// Holds finalized bedGraph intervals plus the sweep state for the chromosome
/// currently being processed. Reads must be supplied in coordinate-sorted order
/// within each chromosome (guaranteed by the BAM read loop).
#[derive(Debug, Default, Clone)]
pub struct TrackAccum {
    settings: TrackSettings,
    /// Finalized bedGraph intervals for all completed chromosomes.
    bedgraph: Vec<BedGraphEntry>,

    // --- sweep state for the active chromosome ---
    /// Name of the chromosome currently being accumulated (empty if none).
    chrom_name: String,
    /// Length of the active chromosome, used to clamp interval ends.
    chrom_size: u32,
    /// Pending coverage deltas at positions at or ahead of the sweep cursor.
    /// `+1` at a block start, `-1` at a block's half-open end.
    deltas: BTreeMap<u32, i32>,
    /// Start of the currently open (constant-depth) run.
    cur_start: u32,
    /// Coverage depth of the currently open run.
    cur_depth: i64,
}

impl TrackAccum {
    fn new(settings: TrackSettings) -> Self {
        Self {
            settings,
            ..Default::default()
        }
    }

    /// Begin accumulating a new chromosome, finalizing any previous one.
    fn start_chrom(&mut self, name: &str, size: u64) {
        self.finish_chrom();
        self.chrom_name.clear();
        self.chrom_name.push_str(name);
        self.chrom_size = size.min(u32::MAX as u64) as u32;
        self.cur_start = 0;
        self.cur_depth = 0;
        debug_assert!(self.deltas.is_empty());
    }

    /// Record a half-open alignment block `[start, end)`, clamped to the
    /// chromosome, as a pair of coverage deltas.
    fn add_block(&mut self, start: u64, end: u64) {
        let size = self.chrom_size as u64;
        if size == 0 || start >= size {
            return;
        }
        let end = end.min(size);
        if end <= start {
            return;
        }
        *self.deltas.entry(start as u32).or_insert(0) += 1;
        *self.deltas.entry(end as u32).or_insert(0) -= 1;
    }

    /// Apply and emit all deltas strictly below `pos`. Safe to call once the
    /// read cursor has advanced past `pos`, since no later read can add
    /// coverage below it.
    fn flush_below(&mut self, pos: u32) {
        let keep = self.deltas.split_off(&pos);
        let flush = std::mem::replace(&mut self.deltas, keep);
        self.consume(flush);
    }

    /// Finalize the active chromosome, emitting all remaining intervals.
    fn finish_chrom(&mut self) {
        if self.chrom_name.is_empty() {
            return;
        }
        let flush = std::mem::take(&mut self.deltas);
        self.consume(flush);
        // A balanced difference array always returns to depth 0.
        debug_assert_eq!(self.cur_depth, 0);
        self.cur_depth = 0;
        self.chrom_name.clear();
    }

    /// Drive the emit state machine over a run of (position, delta) pairs in
    /// ascending position order. Only depth *changes* create interval
    /// boundaries, so net-zero positions merge automatically (matching
    /// bedtools' equal-depth interval merging).
    fn consume(&mut self, deltas: BTreeMap<u32, i32>) {
        for (pos, delta) in deltas {
            if delta == 0 {
                continue;
            }
            if self.cur_depth > 0 {
                self.bedgraph.push((
                    self.chrom_name.clone(),
                    self.cur_start,
                    pos,
                    self.cur_depth as u32,
                ));
            }
            self.cur_start = pos;
            self.cur_depth += delta as i64;
        }
    }

    fn merge(&mut self, mut other: Self) {
        debug_assert!(other.chrom_name.is_empty(), "merge of unfinished track");
        self.bedgraph.append(&mut other.bedgraph);
    }
}

/// Bundle of coverage accumulators for all nf-core/rnaseq bigWig tracks.
#[derive(Debug, Clone)]
pub struct GenomeCovAccum {
    stranded: Strandedness,
    /// tid of the chromosome currently being accumulated (-1 = none).
    current_tid: i32,
    /// Strand-agnostic combined coverage (`-split -bg`).
    pub combined: TrackAccum,
    /// Forward-genome-strand coverage (stranded libraries only).
    pub forward: Option<TrackAccum>,
    /// Reverse-genome-strand coverage (stranded libraries only).
    pub reverse: Option<TrackAccum>,
    /// Scratch buffer for split-read alignment blocks, reused across reads to
    /// avoid a per-record heap allocation in the hot loop.
    blocks_buf: Vec<(u64, u64)>,
}

impl GenomeCovAccum {
    /// Create accumulators for the given strandedness.
    ///
    /// Per-strand tracks are only created for stranded libraries.
    pub fn new(stranded: Strandedness) -> Self {
        let combined = TrackAccum::new(TrackSettings {
            deduplicate_strand: false,
            strand: None,
        });
        let (forward, reverse) = if stranded == Strandedness::Unstranded {
            (None, None)
        } else {
            let (fwd_strand, rev_strand) = match stranded {
                Strandedness::Forward => (StrandFilter::Plus, StrandFilter::Minus),
                Strandedness::Reverse => (StrandFilter::Minus, StrandFilter::Plus),
                Strandedness::Unstranded => unreachable!(),
            };
            (
                Some(TrackAccum::new(TrackSettings {
                    deduplicate_strand: true,
                    strand: Some(fwd_strand),
                })),
                Some(TrackAccum::new(TrackSettings {
                    deduplicate_strand: true,
                    strand: Some(rev_strand),
                })),
            )
        };
        Self {
            stranded,
            current_tid: -1,
            combined,
            forward,
            reverse,
            blocks_buf: Vec::new(),
        }
    }

    /// Merge another worker's (finished) accumulators into this one.
    pub fn merge(&mut self, other: Self) {
        self.combined.merge(other.combined);
        if let (Some(ref mut f), Some(o)) = (&mut self.forward, other.forward) {
            f.merge(o);
        }
        if let (Some(ref mut r), Some(o)) = (&mut self.reverse, other.reverse) {
            r.merge(o);
        }
    }

    /// Finalize the active chromosome on every track. Must be called once all
    /// reads have been processed, before merging or reading results.
    pub fn finish(&mut self) {
        self.combined.finish_chrom();
        if let Some(ref mut f) = self.forward {
            f.finish_chrom();
        }
        if let Some(ref mut r) = self.reverse {
            r.finish_chrom();
        }
        self.current_tid = -1;
    }

    /// Process a mapped BAM record, updating all enabled tracks.
    ///
    /// `tid`/`name`/`size` identify the record's reference. Records must arrive
    /// in coordinate-sorted order, grouped by chromosome.
    pub fn process_read(&mut self, record: &bam::Record, tid: u32, name: &str, size: u64) {
        if record.is_unmapped() {
            return;
        }

        // Chromosome transition: finalize the previous chromosome and reset
        // the sweep state on every track.
        if self.current_tid != tid as i32 {
            self.combined.start_chrom(name, size);
            if let Some(ref mut f) = self.forward {
                f.start_chrom(name, size);
            }
            if let Some(ref mut r) = self.reverse {
                r.start_chrom(name, size);
            }
            self.current_tid = tid as i32;
        }

        // Extract alignment blocks once; strand filtering differs per track but
        // the blocks themselves do not.
        genomecov_blocks(record, &mut self.blocks_buf);
        let pos = record.pos().max(0) as u32;

        add_read_to_track(&mut self.combined, record, &self.blocks_buf);
        if self.stranded != Strandedness::Unstranded {
            if let Some(ref mut f) = self.forward {
                add_read_to_track(f, record, &self.blocks_buf);
            }
            if let Some(ref mut r) = self.reverse {
                add_read_to_track(r, record, &self.blocks_buf);
            }
        }

        // Advance the sweep cursor: positions below the current read's start are
        // now final and can be emitted, releasing their deltas.
        self.combined.flush_below(pos);
        if let Some(ref mut f) = self.forward {
            f.flush_below(pos);
        }
        if let Some(ref mut r) = self.reverse {
            r.flush_below(pos);
        }
    }
}

/// Apply a read's blocks to one track if it passes the track's strand filter.
fn add_read_to_track(track: &mut TrackAccum, record: &bam::Record, blocks: &[(u64, u64)]) {
    let mut is_reverse = record.is_reverse();
    // bedtools `-du`: flip second mate strand so both mates appear on the
    // same fragment strand (dUTP/strand-specific PE libraries).
    if track.settings.deduplicate_strand
        && record.is_paired()
        && !record.is_mate_unmapped()
        && record.flags() & BAM_FREAD2 != 0
    {
        is_reverse = !is_reverse;
    }

    if let Some(filter) = track.settings.strand {
        let want_reverse = filter == StrandFilter::Minus;
        if is_reverse != want_reverse {
            return;
        }
    }

    for &(start, end) in blocks {
        track.add_block(start, end);
    }
}

/// Extract alignment blocks for `-split` mode, matching bedtools `GetBamBlocks`
/// with `breakOnDeletionOps=true` (default, no `-ignoreD`).
fn genomecov_blocks(record: &bam::Record, blocks: &mut Vec<(u64, u64)>) {
    blocks.clear();
    let mut ref_pos = record.pos().max(0) as u64;

    for op in record.cigar().iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len = *len as u64;
                blocks.push((ref_pos, ref_pos + len));
                ref_pos += len;
            }
            Cigar::RefSkip(len) | Cigar::Del(len) => {
                ref_pos += *len as u64;
            }
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
}

/// Merged genome coverage ready for bedGraph / bigWig output.
#[derive(Debug)]
pub struct GenomeCovResult {
    /// Combined strand-agnostic track.
    pub combined: TrackAccum,
    /// Forward-genome-strand track (None for unstranded libraries).
    pub forward: Option<TrackAccum>,
    /// Reverse-genome-strand track (None for unstranded libraries).
    pub reverse: Option<TrackAccum>,
    /// Chromosome names and lengths in BAM header order, carried through so
    /// output does not need to re-open the BAM to recover them.
    pub chrom_sizes: Vec<(String, u64)>,
}

impl GenomeCovAccum {
    /// Finalize the accumulator into a result, attaching the BAM header-ordered
    /// chromosome sizes used for bedGraph ordering and bigWig writing.
    ///
    /// [`finish`](Self::finish) must have been called first.
    pub fn into_result(self, chrom_sizes: Vec<(String, u64)>) -> GenomeCovResult {
        debug_assert_eq!(self.current_tid, -1, "into_result before finish()");
        GenomeCovResult {
            combined: self.combined,
            forward: self.forward,
            reverse: self.reverse,
            chrom_sizes,
        }
    }
}

/// Clip bedGraph intervals to chromosome sizes (UCSC `bedClip` behaviour).
///
/// The streaming accumulator already clamps interval ends to each chromosome's
/// length, so in the normal pipeline this is defensive: it guarantees `bedClip`
/// parity and drops any entry whose chromosome is absent from `chrom_sizes`.
pub fn clip_bedgraph(
    entries: Vec<BedGraphEntry>,
    chrom_sizes: &[(String, u64)],
) -> Vec<BedGraphEntry> {
    let sizes: std::collections::HashMap<&str, u64> =
        chrom_sizes.iter().map(|(n, s)| (n.as_str(), *s)).collect();

    entries
        .into_iter()
        .filter_map(|(chrom, start, end, val)| {
            let size = *sizes.get(chrom.as_str())? as u32;
            if start >= size {
                return None;
            }
            let clipped_end = end.min(size);
            if start >= clipped_end {
                return None;
            }
            Some((chrom, start, clipped_end, val))
        })
        .collect()
}

/// Return a track's bedGraph intervals grouped in `chrom_order` (BAM header)
/// order. Intervals within each chromosome are already in ascending position
/// order from the streaming sweep, so a stable sort by chromosome suffices.
pub fn track_to_bedgraph(track: &TrackAccum, chrom_order: &[(String, u64)]) -> Vec<BedGraphEntry> {
    let order: std::collections::HashMap<&str, usize> = chrom_order
        .iter()
        .enumerate()
        .map(|(i, (n, _))| (n.as_str(), i))
        .collect();

    let mut entries = track.bedgraph.clone();
    entries.sort_by_key(|(chrom, _, _, _)| *order.get(chrom.as_str()).unwrap_or(&usize::MAX));
    entries
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Brute-force reference: per-base coverage from half-open blocks, clamped
    /// to `size`, converted to merged bedGraph intervals. Obviously correct and
    /// independent of the streaming implementation.
    fn brute_force_bedgraph(
        chrom: &str,
        size: u32,
        blocks: &[(u64, u64)],
    ) -> Vec<(String, u32, u32, u32)> {
        let mut depth = vec![0u32; size as usize];
        for &(start, end) in blocks {
            let end = end.min(size as u64);
            for p in start..end {
                if p < size as u64 {
                    depth[p as usize] += 1;
                }
            }
        }
        let mut out = Vec::new();
        let mut i = 0usize;
        while i < depth.len() {
            if depth[i] == 0 {
                i += 1;
                continue;
            }
            let d = depth[i];
            let start = i;
            while i < depth.len() && depth[i] == d {
                i += 1;
            }
            out.push((chrom.to_string(), start as u32, i as u32, d));
        }
        out
    }

    /// Feed blocks (as one-block "reads", sorted by start) through the streaming
    /// accumulator and return the combined bedGraph.
    fn stream_bedgraph(chrom: &str, size: u64, blocks: &[(u64, u64)]) -> Vec<BedGraphEntry> {
        let mut track = TrackAccum::new(TrackSettings::default());
        track.start_chrom(chrom, size);
        // Blocks must be applied in ascending start order with the sweep cursor.
        let mut sorted = blocks.to_vec();
        sorted.sort_by_key(|&(s, _)| s);
        for &(s, e) in &sorted {
            track.add_block(s, e);
            track.flush_below(s as u32);
        }
        track.finish_chrom();
        track.bedgraph
    }

    #[test]
    fn test_simple_overlap() {
        let blocks = [(2u64, 5u64), (2, 5)];
        assert_eq!(
            stream_bedgraph("chr1", 10, &blocks),
            brute_force_bedgraph("chr1", 10, &blocks)
        );
    }

    #[test]
    fn test_adjacent_and_nested() {
        let blocks = [(0u64, 10u64), (3, 6), (6, 9), (9, 12)];
        assert_eq!(
            stream_bedgraph("chr1", 20, &blocks),
            brute_force_bedgraph("chr1", 20, &blocks)
        );
    }

    #[test]
    fn test_clamping_past_chrom_end() {
        let blocks = [(8u64, 15u64), (5, 12)];
        assert_eq!(
            stream_bedgraph("chr1", 10, &blocks),
            brute_force_bedgraph("chr1", 10, &blocks)
        );
    }

    #[test]
    fn test_net_zero_position_merges() {
        // A block ending exactly where another begins keeps depth constant
        // across the boundary, so the intervals must merge.
        let blocks = [(0u64, 5u64), (5, 10)];
        let got = stream_bedgraph("chr1", 10, &blocks);
        assert_eq!(got, vec![("chr1".to_string(), 0, 10, 1)]);
        assert_eq!(got, brute_force_bedgraph("chr1", 10, &blocks));
    }

    #[test]
    fn test_randomised_against_brute_force() {
        // Deterministic LCG so the test is reproducible without rand.
        let mut state = 0x1234_5678u64;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as u32
        };
        for trial in 0..200 {
            let size = 50 + (next() % 150);
            let n = next() % 60;
            let mut blocks = Vec::new();
            for _ in 0..n {
                let start = next() % size;
                let len = 1 + next() % 20;
                blocks.push((start as u64, (start + len) as u64));
            }
            let stream = stream_bedgraph("c", size as u64, &blocks);
            let brute = brute_force_bedgraph("c", size, &blocks);
            assert_eq!(stream, brute, "mismatch on trial {trial}");
        }
    }

    #[test]
    fn test_clip_bedgraph() {
        let entries = vec![("chr1".to_string(), 100, 200, 5)];
        let sizes = vec![("chr1".to_string(), 150u64)];
        let clipped = clip_bedgraph(entries, &sizes);
        assert_eq!(clipped, vec![("chr1".to_string(), 100, 150, 5)]);
    }

    #[test]
    fn test_track_to_bedgraph_orders_by_header() {
        let mut track = TrackAccum::new(TrackSettings::default());
        // Emit chr2 before chr1 to exercise reordering.
        track.bedgraph = vec![("chr2".to_string(), 0, 5, 1), ("chr1".to_string(), 0, 5, 1)];
        let order = vec![("chr1".to_string(), 100u64), ("chr2".to_string(), 100u64)];
        let out = track_to_bedgraph(&track, &order);
        assert_eq!(
            out,
            vec![("chr1".to_string(), 0, 5, 1), ("chr2".to_string(), 0, 5, 1),]
        );
    }
}
