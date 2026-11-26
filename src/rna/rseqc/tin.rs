//! TIN (Transcript Integrity Number) analysis.
//!
//! Measures transcript integrity via Shannon entropy of read coverage
//! uniformity across sampled exonic positions. Reimplementation of
//! RSeQC's `tin.py` tool.

use rust_htslib::bam;
use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap, HashMap, HashSet};
use std::hash::{BuildHasher, Hasher};
use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use indexmap::IndexMap;
use log::debug;
use rand::Rng;

use crate::gtf::{Gene, Transcript};

// ===================================================================
// Deterministic hash state for reproducible TIN
// ===================================================================

/// A [`BuildHasher`] that produces deterministic hashes seeded from a `u64`.
///
/// When `--tin-seed` is provided, the seed is used directly, making the
/// [`HashSet`] layout (and therefore threshold-crossing behaviour during
/// parallel accumulation) identical across runs. Without a seed, a random
/// value is drawn so the default non-deterministic behaviour is preserved.
#[derive(Clone, Debug)]
pub(crate) struct TinHashState(u64);

impl BuildHasher for TinHashState {
    type Hasher = std::collections::hash_map::DefaultHasher;

    fn build_hasher(&self) -> Self::Hasher {
        let mut h = std::collections::hash_map::DefaultHasher::new();
        h.write_u64(self.0);
        h
    }
}

impl TinHashState {
    fn new(seed: Option<u64>) -> Self {
        TinHashState(seed.unwrap_or_else(|| rand::rng().next_u64()))
    }
}

// ===================================================================
// Data structures
// ===================================================================

/// A single transcript's sampling metadata.
#[derive(Debug, Clone)]
pub struct TranscriptSampling {
    /// Gene identifier (gene_id from GTF).
    pub gene_id: String,
    /// Chromosome name (original case, for output).
    pub chrom: String,
    /// Chromosome name uppercased (for lookup matching).
    pub chrom_upper: String,
    /// Transcript start position (0-based).
    pub tx_start: u64,
    /// Transcript end position (0-based, exclusive).
    pub tx_end: u64,
    /// Sorted exon blocks as (start, end) in 0-based half-open coords.
    #[allow(dead_code)] // stored for API completeness; values used during construction
    pub exon_regions: Vec<(u64, u64)>,
    /// Total exonic bases across all exon blocks.
    #[allow(dead_code)] // stored for API completeness; values used during construction
    pub exon_length: u64,
    /// Unique genomic positions sampled within exonic regions (sorted, deduped).
    /// Used for coverage counting — each position gets one coverage slot.
    pub sampled_positions: Vec<u64>,
    /// Total position count including duplicates from upstream's `chose_bases`.
    /// Used as the denominator `l` in the TIN formula to match upstream behavior.
    /// For the small mRNA case, upstream does NOT dedup chose_bases, so
    /// `total_position_count > len(sampled_positions)` when boundary positions
    /// overlap with exonic positions.
    pub total_position_count: usize,
}

/// Index for fast lookup of sampled positions during BAM processing.
///
/// Maps chromosome → sorted list of (genomic_pos, transcript_idx, slot_idx)
/// for binary search during per-read dispatch.
#[derive(Debug, Default)]
pub struct TinIndex {
    /// All transcripts with their sampling metadata.
    pub transcripts: Vec<TranscriptSampling>,
    /// Per-chromosome sorted position entries for binary search.
    /// Each entry: (genomic_position, transcript_index, slot_index).
    pub chrom_positions: HashMap<String, Vec<(u64, u32, u32)>>,
    /// Per-chromosome transcript spans for overlap pre-filtering.
    /// Each entry: (tx_start, tx_end, transcript_index), sorted by tx_start.
    pub chrom_spans: HashMap<String, Vec<(u64, u64, u32)>>,
}

/// Result of TIN computation for a single transcript.
#[derive(Debug)]
pub struct TinResult {
    /// Gene identifier.
    pub gene_id: String,
    /// Chromosome name.
    pub chrom: String,
    /// Transcript start (0-based).
    pub tx_start: u64,
    /// Transcript end (0-based, exclusive).
    pub tx_end: u64,
    /// Computed TIN score (0-100), or 0.0 if below coverage threshold.
    pub tin: f64,
    /// Whether the transcript passed the minimum coverage threshold.
    /// Only transcripts with `passed_threshold == true` contribute to
    /// summary statistics, matching upstream RSeQC's `sample_TINs` list.
    pub passed_threshold: bool,
}

/// Collection of all TIN results for a BAM file.
#[derive(Debug)]
pub struct TinResults {
    /// Per-transcript TIN scores.
    pub transcripts: Vec<TinResult>,
}

impl TinResults {
    /// Number of transcripts with computed TIN scores.
    pub fn len(&self) -> usize {
        self.transcripts.len()
    }

    /// Whether there are no transcripts with computed TIN scores.
    pub fn is_empty(&self) -> bool {
        self.transcripts.is_empty()
    }
}

// ===================================================================
// TIN index construction
// ===================================================================

impl TinIndex {
    /// Build TIN index from GTF gene annotations.
    ///
    /// Emits every transcript from every gene, matching RSeQC's
    /// per-transcript behaviour.
    ///
    /// tin.py writes one row per BED line, in BED order. nf-core/rnaseq
    /// (3.25.0 and later) builds that BED with the `eautils/gtf2bed` module,
    /// which sorts transcripts by chromosome name (string comparison), then by
    /// the start of each transcript's *first exon line* in the GTF (for
    /// minus-strand transcripts that is usually the highest exon, not the
    /// lowest), then by transcript ID (string comparison). The transcripts are
    /// indexed in that order so the output rows follow it. (Older pipelines
    /// used a `gtf2bed` that broke ties by Perl hash order, which is random.)
    ///
    /// `gtf2bed` skips exons without a `transcript_id` (and, being Perl, an
    /// ID or chromosome named `0`), so tin.py has no row for them; neither
    /// does RustQC.
    pub fn from_genes(genes: &IndexMap<String, Gene>, sample_size: usize) -> Self {
        let mut index = TinIndex::default();

        let mut bed_order: Vec<(&Gene, &Transcript)> = genes
            .values()
            .flat_map(|gene| gene.transcripts.iter().map(move |tx| (gene, tx)))
            .filter(|(_, tx)| in_gtf2bed(tx))
            .collect();
        bed_order.sort_by(|(_, a), (_, b)| {
            a.chrom
                .cmp(&b.chrom)
                .then(a.first_exon_start.cmp(&b.first_exon_start))
                .then(a.transcript_id.cmp(&b.transcript_id))
        });

        for (gene, tx) in bed_order {
            let exon_regions: Vec<(u64, u64)> = tx
                .exons
                .iter()
                .map(|&(s, e)| (s.saturating_sub(1), e)) // GTF 1-based to 0-based half-open
                .collect();
            let exon_length: u64 = exon_regions.iter().map(|(s, e)| e - s).sum();

            if exon_length == 0 {
                continue;
            }

            let chrom = gene.chrom.clone();
            let chrom_upper = chrom.to_uppercase();
            let tx_start = exon_regions.first().map(|r| r.0).unwrap_or(0);
            let tx_end = exon_regions.last().map(|r| r.1).unwrap_or(0);

            let (sampled, total_count) =
                sample_exonic_positions(&exon_regions, sample_size, tx_start, tx_end);
            if sampled.is_empty() {
                continue;
            }

            index.add_transcript(TranscriptSampling {
                gene_id: tx.transcript_id.clone(),
                chrom,
                chrom_upper,
                tx_start,
                tx_end,
                exon_regions,
                exon_length,
                sampled_positions: sampled,
                total_position_count: total_count,
            });
        }

        index.build();
        index
    }

    /// Add a transcript to the index (before calling `build()`).
    fn add_transcript(&mut self, tx: TranscriptSampling) {
        let tx_idx = self.transcripts.len() as u32;
        let chrom_upper = tx.chrom_upper.clone();
        let tx_start = tx.tx_start;
        let tx_end = tx.tx_end;

        // Add sampled positions to per-chromosome position list
        for (slot_idx, &pos) in tx.sampled_positions.iter().enumerate() {
            self.chrom_positions
                .entry(chrom_upper.clone())
                .or_default()
                .push((pos, tx_idx, slot_idx as u32));
        }

        // Add transcript span
        self.chrom_spans
            .entry(chrom_upper)
            .or_default()
            .push((tx_start, tx_end, tx_idx));

        self.transcripts.push(tx);
    }

    /// Sort per-chromosome data for binary search.
    fn build(&mut self) {
        for positions in self.chrom_positions.values_mut() {
            positions.sort_unstable();
        }
        for spans in self.chrom_spans.values_mut() {
            spans.sort_unstable();
        }
        debug!(
            "TIN index: {} transcripts, {} chromosomes",
            self.transcripts.len(),
            self.chrom_positions.len()
        );
    }
}

/// Whether nf-core's `gtf2bed` writes a BED row for this transcript: its
/// `next unless $id && $f[0]` skips features without a `transcript_id`
/// attribute, and Perl also treats an empty string or `"0"` as false.
fn in_gtf2bed(tx: &Transcript) -> bool {
    let truthy = |s: &str| !s.is_empty() && s != "0";
    tx.has_transcript_id && truthy(&tx.transcript_id) && truthy(&tx.chrom)
}

/// Sample positions within exonic regions, matching upstream RSeQC `tin.py`'s
/// `genomic_positions()` algorithm.
///
/// Upstream uses 1-based genomic coordinates internally. All positions returned
/// here are 1-based to match pileup semantics (pysam's `pileupcolumn.pos + 1`).
///
/// When `mRNA_size <= sample_size`: returns ALL exonic bases (1-based) plus
/// `tx_start+1` and `tx_end` as boundary markers. The upstream does NOT dedup
/// these, so duplicates inflate `len(chose_bases)` which becomes the TIN
/// denominator `l`. We reproduce this by returning non-unique positions.
/// When `mRNA_size > sample_size`: returns every `step_size`-th exonic base
/// plus exon boundary positions (first and last base of each exon), deduplicated
/// and sorted.
/// Returns `(unique_positions, total_count)` where:
/// - `unique_positions`: sorted, deduplicated genomic positions for coverage counting
/// - `total_count`: length including duplicates, used as TIN denominator `l`
fn sample_exonic_positions(
    exon_regions: &[(u64, u64)],
    n: usize,
    tx_start: u64,
    tx_end: u64,
) -> (Vec<u64>, usize) {
    let mrna_size: u64 = exon_regions.iter().map(|(s, e)| e - s).sum();
    if mrna_size == 0 {
        return (Vec::new(), 0);
    }

    // Build 1-based list of all exonic positions (upstream: range(st+1, end+1))
    // exon_regions are 0-based half-open: (start, end) means bases start..end-1
    // 1-based equivalents: start+1 .. end (inclusive)
    if mrna_size <= n as u64 {
        // Small transcript: chose_bases = [tx_start+1, tx_end] + all exonic bases
        // Upstream does NOT deduplicate — duplicates inflate len(chose_bases)
        // which is used as the TIN denominator `l`.
        let mut positions: Vec<u64> = Vec::with_capacity(mrna_size as usize + 2);
        positions.push(tx_start + 1); // transcript boundary (1-based)
        positions.push(tx_end); // transcript boundary (1-based, 0-based exclusive = 1-based inclusive)
        for &(start, end) in exon_regions {
            for pos in (start + 1)..=end {
                positions.push(pos);
            }
        }
        positions.sort_unstable();
        let total_count = positions.len(); // includes duplicates

        // Dedup for coverage counting (pileup only returns unique positions)
        positions.dedup();
        (positions, total_count)
    } else {
        // Large transcript: strided sampling + exon boundaries
        // Build all 1-based exonic positions
        let mut gene_all_base: Vec<u64> = Vec::with_capacity(mrna_size as usize);
        let mut exon_bounds: Vec<u64> = Vec::new();
        for &(start, end) in exon_regions {
            for pos in (start + 1)..=end {
                gene_all_base.push(pos);
            }
            // Exon boundary positions (1-based)
            exon_bounds.push(start + 1); // first base of exon
            exon_bounds.push(end); // last base of exon
        }

        // Strided sampling: step_size = int(mRNA_size / sample_size)
        let step_size = (mrna_size as usize) / n;
        let step_size = step_size.max(1);

        let mut chose_bases: Vec<u64> = (0..gene_all_base.len())
            .step_by(step_size)
            .map(|i| gene_all_base[i])
            .collect();

        // Merge exon boundaries with sampled positions (upstream: uniqify(exon_bounds + chose_bases))
        // uniqify preserves first occurrence order and removes duplicates
        let mut all_positions = exon_bounds;
        all_positions.append(&mut chose_bases);

        let mut seen = HashSet::new();
        all_positions.retain(|x| seen.insert(*x));
        all_positions.sort_unstable();

        let total_count = all_positions.len(); // already deduped by uniqify
        (all_positions, total_count)
    }
}

// ===================================================================
// TIN accumulator
// ===================================================================

/// Minimum base quality for a base to count towards pileup coverage.
///
/// tin.py calls `samfile.pileup()` with pysam's defaults, and `min_base_quality=13`
/// makes `PileupColumn.pileups` skip bases below this quality.
const PILEUP_MIN_BASE_QUALITY: u8 = 13;

/// One sampled position covered by an aligned base of a read.
#[derive(Debug, Clone, Copy)]
struct CoveredBase {
    /// Start of the range of entries in the chromosome's position list at this position
    first: u32,
    /// End (exclusive) of that range
    last: u32,
    /// Query offset of the read base aligned to this position
    qpos: u32,
}

/// pysam's default `max_depth` for `pileup()`, which it passes to htslib's
/// `bam_mplp_set_maxcnt`.
const PILEUP_MAX_DEPTH: usize = 8000;

/// What one transcript's pileup does with a read (htslib `bam_plp_push`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Push {
    /// The read is outside the transcript's region, so `fetch()` never returns it
    Absent,
    /// Dropped by the depth cap (`max_depth`)
    Dropped,
    /// Not buffered: no reference length, and the pileup is already at its position
    Silent,
    /// Buffered and piled up
    Kept,
}

/// A properly paired pileup read, as htslib's overlap hash sees it.
///
/// The hash (`overlap_push` / `overlap_remove`) is keyed by read name, so
/// supplementary copies and their mates compete for one slot. Each tin.py
/// transcript gets its own pileup over its own region, so the slot's history
/// is replayed per transcript from the records of one name (see
/// [`overlap_partner`]).
#[derive(Debug)]
struct MateRec {
    /// Alignment start (0-based)
    pos: u64,
    /// Reference end: start plus the CIGAR's reference length (may equal `pos`)
    end: u64,
    /// End as `fetch()` sees it (at least `pos + 1`); the read is in a
    /// transcript's pileup iff this is above the transcript start
    fetch_end: u64,
    /// Pileup position when the read was pushed: reads ending before it had
    /// already left the buffer (and the hash)
    frontier: u64,
    /// Passes `overlap_push`'s checks (mate mapped on the same chromosome, not
    /// a "wild" CIGAR); others never look up or enter the hash
    candidate: bool,
    /// Stored in the hash when no entry exists (its mate is still to come)
    insertable: bool,
    /// For a read without reference length: transcripts starting below this
    /// value already have a read at this position, so they do not buffer it
    silent_below: u64,
    /// Depth-capped transcripts whose pileup dropped this read (sorted)
    dropped_in: Vec<u32>,
    /// What a tweak needs, kept only for reads that may be stored in the hash
    /// and span a sampled position
    data: Option<Box<MateData>>,
}

/// Alignment of a stored mate, for [`tweak_overlap_quality`].
#[derive(Debug)]
struct MateData {
    /// Raw BAM CIGAR
    cigar: Vec<u32>,
    /// 4-bit packed read sequence, as stored in the BAM
    seq: Vec<u8>,
    /// Base qualities as stored in the BAM (what the read was counted with)
    qual: Vec<u8>,
    /// The read's covered sampled positions
    bases: Vec<CoveredBase>,
}

impl MateRec {
    /// How the pileup of a transcript starting at `tx_start` handles this read.
    /// `tx` names the transcript when its own depth-cap drops must be honoured.
    fn push_in(&self, tx_start: u64, tx: Option<u32>) -> Push {
        if self.fetch_end <= tx_start {
            Push::Absent
        } else if tx.is_some_and(|t| self.dropped_in.binary_search(&t).is_ok()) {
            Push::Dropped
        } else if tx_start < self.silent_below {
            Push::Silent
        } else {
            Push::Kept
        }
    }
}

/// Records of one read name whose hash slot may still be in use.
#[derive(Debug, Default)]
struct MateGroup {
    /// Records in BAM order
    recs: Vec<MateRec>,
    /// Largest `end` of the records: once the pileup has moved past it, every
    /// record has left the buffer and the slot is empty in every transcript
    max_end: u64,
}

/// Replay htslib's overlap hash for one read name in the pileup of the
/// transcript starting at `tx_start`, and return the earlier record that the
/// last record of `recs` is paired with (`tweak_overlap_quality(a, b)`), if any.
///
/// The pileup only sees reads overlapping the transcript. A record leaves the
/// buffer, and empties the slot (`overlap_remove` is by name), once the pileup
/// has moved past its end; a later record is pushed after that iff its
/// frontier is beyond that end. A record dropped by the depth cap also empties
/// the slot.
fn overlap_partner(
    recs: &[MateRec],
    tx_start: u64,
    tx: Option<u32>,
    gone: &mut Vec<bool>,
) -> Option<usize> {
    let last = recs.len().checked_sub(1)?;
    gone.clear();
    gone.resize(recs.len(), false);
    let mut slot: Option<usize> = None;
    for (i, x) in recs.iter().enumerate() {
        for (h, r) in recs[..i].iter().enumerate() {
            if !gone[h] && r.end < x.frontier && r.push_in(tx_start, tx) == Push::Kept {
                gone[h] = true;
                slot = None;
            }
        }
        match x.push_in(tx_start, tx) {
            Push::Kept if x.candidate => {
                if let Some(a) = slot {
                    if i == last {
                        return Some(a);
                    }
                    slot = None;
                } else if x.insertable {
                    slot = Some(i);
                }
            }
            Push::Dropped => slot = None,
            _ => {}
        }
    }
    None
}

/// Which earlier record the newest record of a group is paired with, for
/// every transcript that can contain both.
struct Pairing {
    /// Transcript-start thresholds where the outcome may change (sorted)
    thresholds: Vec<u64>,
    /// Partner for transcripts starting in `[thresholds[j-1], thresholds[j])`
    by_class: Vec<Option<usize>>,
    /// Partner for depth-capped transcripts that dropped one of the records
    overrides: Vec<(u32, Option<usize>)>,
}

impl Pairing {
    fn new(recs: &[MateRec], gone: &mut Vec<bool>) -> Self {
        let mut thresholds: Vec<u64> = recs
            .iter()
            .flat_map(|r| [r.fetch_end, r.silent_below])
            .filter(|&v| v > 0)
            .collect();
        thresholds.sort_unstable();
        thresholds.dedup();
        let by_class = (0..=thresholds.len())
            .map(|j| {
                let tx_start = if j == 0 { 0 } else { thresholds[j - 1] };
                overlap_partner(recs, tx_start, None, gone)
            })
            .collect();
        let mut capped: Vec<u32> = recs
            .iter()
            .flat_map(|r| r.dropped_in.iter().copied())
            .collect();
        capped.sort_unstable();
        capped.dedup();
        let overrides = capped.into_iter().map(|t| (t, None)).collect();
        Pairing {
            thresholds,
            by_class,
            overrides,
        }
    }

    /// Fill in the partners of the depth-capped transcripts.
    fn resolve_overrides(&mut self, recs: &[MateRec], index: &TinIndex, gone: &mut Vec<bool>) {
        for (t, partner) in &mut self.overrides {
            let tx_start = index.transcripts[*t as usize].tx_start;
            *partner = overlap_partner(recs, tx_start, Some(*t), gone);
        }
    }

    /// Partner of the newest record in the pileup of transcript `tx`.
    fn partner(&self, tx: u32, tx_start: u64) -> Option<usize> {
        if let Ok(k) = self.overrides.binary_search_by_key(&tx, |o| o.0) {
            return self.overrides[k].1;
        }
        self.by_class[self.thresholds.partition_point(|&v| v <= tx_start)]
    }

    /// Distinct partners over all transcripts.
    fn partners(&self) -> Vec<usize> {
        let mut out: Vec<usize> = self
            .by_class
            .iter()
            .chain(self.overrides.iter().map(|o| &o.1))
            .flatten()
            .copied()
            .collect();
        out.sort_unstable();
        out.dedup();
        out
    }
}

/// Pileup state of one transcript whose depth may reach `max_depth`.
///
/// htslib drops a read when the pileup is already at its start position and
/// holds more than `max_depth` reads (`bam_plp_push`: `mp->cnt > maxcnt`,
/// where the count includes one spare node). Which reads that drops depends on
/// the transcript's own region, so such transcripts track their buffer here.
#[derive(Debug)]
struct CappedPileup {
    /// Transcript end (reads starting here or later are not fetched)
    tx_end: u64,
    /// Position of the last read pushed, if the pileup is at it
    last_pos: Option<u64>,
    /// Ends of the buffered reads
    ends: BinaryHeap<Reverse<u64>>,
    /// Largest end of a dropped read (the pileup differs while it would still be buffered)
    dropped_until: Option<u64>,
}

impl CappedPileup {
    /// `bam_plp_push` for a read of this transcript's region.
    fn push(&mut self, pos: u64, end: u64) -> Push {
        if self.last_pos != Some(pos) {
            // First read at this position: never dropped. Reads ending before
            // it leave the buffer once the pileup moves here.
            self.last_pos = Some(pos);
            while self.ends.peek().is_some_and(|e| e.0 < pos) {
                self.ends.pop();
            }
            self.ends.push(Reverse(end));
            Push::Kept
        } else if self.ends.len() + 1 > PILEUP_MAX_DEPTH {
            self.dropped_until = Some(self.dropped_until.map_or(end, |d| d.max(end)));
            Push::Dropped
        } else if end > pos {
            self.ends.push(Reverse(end));
            Push::Kept
        } else {
            Push::Silent
        }
    }
}

/// One mate as seen by [`tweak_overlap_quality`].
struct MateView<'a> {
    /// Alignment start (0-based)
    pos: i64,
    /// Raw BAM CIGAR
    cigar: &'a [u32],
    /// 4-bit packed read sequence
    seq: &'a [u8],
    /// Base qualities, rewritten in place
    qual: &'a mut [u8],
}

/// BAM CIGAR operation codes used by the htslib port below.
const CIGAR_MATCH: u32 = 0;
const CIGAR_INS: u32 = 1;
const CIGAR_DEL: u32 = 2;
const CIGAR_REF_SKIP: u32 = 3;
const CIGAR_SOFT_CLIP: u32 = 4;
const CIGAR_HARD_CLIP: u32 = 5;
const CIGAR_PAD: u32 = 6;
const CIGAR_EQUAL: u32 = 7;
const CIGAR_DIFF: u32 = 8;

/// Outcome of a CIGAR cursor step (htslib returns BAM_CMATCH, -1 or -2).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Step {
    /// Positioned on an aligned (M/=/X) base
    Match,
    /// End of the CIGAR, or the requested position is not covered
    End,
    /// Unexpected CIGAR operation
    Error,
}

/// Position within a CIGAR, ported from htslib's `cigar_iref2iseq_set` /
/// `cigar_iref2iseq_next` (htslib 1.24, as bundled with pysam 0.24).
struct CigarCursor<'a> {
    /// Raw BAM CIGAR operations
    ops: &'a [u32],
    /// Index of the current operation (htslib's `*cigar` pointer)
    k: usize,
    /// Offset within the current operation
    icig: i64,
    /// Query offset of the current base
    iseq: i64,
    /// Reference offset of the current base relative to the alignment start
    iref: i64,
}

impl<'a> CigarCursor<'a> {
    /// Position on the first aligned base at or after reference offset `target`
    /// (`cigar_iref2iseq_set`).
    fn set(ops: &'a [u32], target: i64) -> (Self, Step) {
        let mut c = CigarCursor {
            ops,
            k: 0,
            icig: 0,
            iseq: 0,
            iref: 0,
        };
        let mut pos = target;
        if pos < 0 {
            return (c, Step::End);
        }
        while let Some(&op) = c.ops.get(c.k) {
            let ncig = (op >> 4) as i64;
            match op & 0xf {
                CIGAR_SOFT_CLIP | CIGAR_INS => {
                    c.k += 1;
                    c.iseq += ncig;
                    c.icig = 0;
                }
                CIGAR_HARD_CLIP | CIGAR_PAD => {
                    c.k += 1;
                    c.icig = 0;
                }
                CIGAR_MATCH | CIGAR_EQUAL | CIGAR_DIFF => {
                    pos -= ncig;
                    if pos < 0 {
                        c.icig = ncig + pos;
                        c.iseq += c.icig;
                        c.iref += c.icig;
                        return (c, Step::Match);
                    }
                    c.k += 1;
                    c.iseq += ncig;
                    c.icig = 0;
                    c.iref += ncig;
                }
                CIGAR_DEL | CIGAR_REF_SKIP => {
                    pos = (pos - ncig).max(0);
                    c.k += 1;
                    c.icig = 0;
                    c.iref += ncig;
                }
                _ => return (c, Step::Error),
            }
        }
        c.iseq = -1;
        (c, Step::End)
    }

    /// Step to the next aligned base (`cigar_iref2iseq_next`).
    fn next(&mut self) -> Step {
        while let Some(&op) = self.ops.get(self.k) {
            let ncig = (op >> 4) as i64;
            match op & 0xf {
                CIGAR_MATCH | CIGAR_EQUAL | CIGAR_DIFF => {
                    if self.icig >= ncig - 1 {
                        self.icig = -1;
                        self.k += 1;
                        continue;
                    }
                    self.iseq += 1;
                    self.icig += 1;
                    self.iref += 1;
                    return Step::Match;
                }
                CIGAR_DEL | CIGAR_REF_SKIP => {
                    self.k += 1;
                    self.iref += ncig;
                    self.icig = -1;
                }
                CIGAR_INS | CIGAR_SOFT_CLIP => {
                    self.k += 1;
                    self.iseq += ncig;
                    self.icig = -1;
                }
                CIGAR_HARD_CLIP | CIGAR_PAD => {
                    self.k += 1;
                    self.icig = -1;
                }
                _ => return Step::Error,
            }
        }
        self.iseq = -1;
        self.iref = -1;
        Step::End
    }

    /// Whether the operation before the current one is a deletion
    /// (htslib: `cigar > bam_get_cigar(b) && bam_cigar_op(cigar[-1]) == BAM_CDEL`).
    fn after_deletion(&self) -> bool {
        self.k > 0
            && self
                .ops
                .get(self.k - 1)
                .is_some_and(|op| op & 0xf == CIGAR_DEL)
    }
}

/// htslib's `__ac_X31_hash_string` followed by `__ac_Wang_hash` (khash.h).
fn htslib_name_hash(qname: &[u8]) -> u32 {
    let mut h: u32 = 0;
    if let Some((&first, rest)) = qname.split_first() {
        h = first as i8 as u32;
        if h != 0 {
            for &c in rest {
                if c == 0 {
                    break;
                }
                h = (h << 5).wrapping_sub(h).wrapping_add(c as i8 as u32);
            }
        }
    }
    let mut key = h;
    key = key.wrapping_add(!(key << 15));
    key ^= key >> 10;
    key = key.wrapping_add(key << 3);
    key ^= key >> 6;
    key = key.wrapping_add(!(key << 11));
    key ^= key >> 16;
    key
}

/// 4-bit base code at query offset `i` of a packed BAM sequence (`bam_seqi`).
fn seqi(seq: &[u8], i: i64) -> u8 {
    let byte = seq.get((i >> 1) as usize).copied().unwrap_or(0);
    if i & 1 == 0 {
        byte >> 4
    } else {
        byte & 0xf
    }
}

/// `x * 0.8` truncated to a quality byte, as C's implicit double-to-uint8 conversion.
fn scale_qual(q: u8) -> u8 {
    (q as f64 * 0.8) as u8
}

/// Port of htslib 1.24's `tweak_overlap_quality`, which the pileup engine runs
/// when the second mate of an overlapping pair arrives (pysam's default
/// `ignore_overlaps=True`). Qualities are rewritten in place.
///
/// Where both mates align to the same reference base, one mate (chosen by a
/// hash of the read name for matching bases, else the higher quality) keeps the
/// summed or 0.8-scaled quality and the other gets 0. Inside a deletion of one
/// mate the other mate's bases are scaled or zeroed. The walk skips positions
/// after a reference skip (`N`), so those keep both bases, exactly as htslib.
///
/// # Arguments
/// * `a` - The first (left) mate
/// * `b` - The second mate
/// * `qname` - The read name (selects which mate keeps matching bases)
fn tweak_overlap_quality(a: &mut MateView, b: &mut MateView, qname: &[u8]) {
    let mut iref = b.pos;
    let (mut ac, mut a_ret) = CigarCursor::set(a.cigar, iref - a.pos);
    if a_ret != Step::Match {
        return;
    }
    let (mut bc, mut b_ret) = CigarCursor::set(b.cigar, iref - b.pos);
    if b_ret != Step::Match {
        return;
    }

    let (amul, bmul): (u8, u8) = if htslib_name_hash(qname) & 1 != 0 {
        (1, 0)
    } else {
        (0, 1)
    };
    let a_len = a.qual.len() as i64;
    let b_len = b.qual.len() as i64;
    // Quality index guard: htslib only bounds-checks `iseq > l_qseq`
    let idx = |i: i64, len: i64| -> Option<usize> { (i >= 0 && i < len).then_some(i as usize) };

    loop {
        while a_ret == Step::Match && ac.iref >= 0 && ac.iref < iref - a.pos {
            a_ret = ac.next();
        }
        if a_ret != Step::Match {
            return;
        }
        while b_ret == Step::Match && bc.iref >= 0 && bc.iref < iref - b.pos {
            b_ret = bc.next();
        }
        if b_ret != Step::Match {
            return;
        }
        iref = iref.max(ac.iref + a.pos).max(bc.iref + b.pos);
        iref += 1;

        if ac.iref + a.pos != bc.iref + b.pos {
            if ac.iref + a.pos < bc.iref + b.pos && bc.after_deletion() {
                // Deletion in b: catch a up, scaling or zeroing its qualities
                loop {
                    if let Some(i) = idx(ac.iseq, a_len) {
                        a.qual[i] = if amul != 0 { scale_qual(a.qual[i]) } else { 0 };
                    }
                    if ac.next() != Step::Match {
                        return;
                    }
                    if ac.iref + a.pos >= bc.iref + b.pos {
                        break;
                    }
                }
            } else if ac.after_deletion() {
                // Deletion in a: catch b up
                loop {
                    if let Some(i) = idx(bc.iseq, b_len) {
                        b.qual[i] = if bmul != 0 { scale_qual(b.qual[i]) } else { 0 };
                    }
                    if bc.next() != Step::Match {
                        return;
                    }
                    if bc.iref + b.pos >= ac.iref + a.pos {
                        break;
                    }
                }
            } else {
                // Anything else (e.g. reference skip) is not handled by htslib
                continue;
            }
        }

        if ac.iseq > a_len || bc.iseq > b_len {
            return;
        }
        let (Some(ai), Some(bi)) = (idx(ac.iseq, a_len), idx(bc.iseq, b_len)) else {
            return;
        };
        let (qa, qb) = (a.qual[ai], b.qual[bi]);
        if seqi(a.seq, ac.iseq) == seqi(b.seq, bc.iseq) {
            let sum = (qa as u32 + qb as u32).min(200) as u8;
            a.qual[ai] = amul * sum;
            b.qual[bi] = bmul * sum;
        } else if qa > qb {
            a.qual[ai] = scale_qual(qa);
            b.qual[bi] = 0;
        } else if qa < qb {
            b.qual[bi] = scale_qual(qb);
            a.qual[ai] = 0;
        } else {
            a.qual[ai] = if amul != 0 { scale_qual(qa) } else { 0 };
            b.qual[bi] = if bmul != 0 { scale_qual(qb) } else { 0 };
        }
    }
}

/// Per-read accumulator for TIN computation.
///
/// Tracks coverage at sampled positions and unique read start positions
/// per transcript during the BAM pass. Unique start tracking matches
/// RSeQC's `check_min_reads()`, which uses a `set()` of read start
/// positions.
///
/// Coverage reproduces what tin.py's `genebody_coverage()` gets from
/// `pysam.AlignmentFile.pileup(chrom, tx_start, tx_end, truncate=True)` with
/// its default settings, one pileup per transcript: reads with an unmapped,
/// secondary, QC-fail or duplicate flag are skipped, and so are orphans
/// (`ignore_orphans`: paired reads without the proper-pair flag); overlapping
/// mates are counted once (`ignore_overlaps`, htslib's overlap hash keyed by
/// read name); bases below [`PILEUP_MIN_BASE_QUALITY`] are not counted; and
/// where more than `max_depth` reads pile up, htslib drops reads
/// ([`PILEUP_MAX_DEPTH`]). The last two depend on which reads the
/// transcript's region holds, so they are tracked per transcript where they
/// can differ, still in one pass over the BAM.
///
/// To reduce memory pressure, each transcript's `unique_starts` set is
/// capped at `min_cov + 1` entries. Once the threshold is exceeded the
/// set is cleared and `exceeded_threshold` is set, since only the
/// pass/fail decision matters for the TIN score.
#[derive(Debug)]
pub struct TinAccum {
    /// Per-transcript per-slot coverage counts.
    /// Indexed as `coverage[tx_idx][slot_idx]`.
    pub coverage: Vec<Vec<u32>>,
    /// Per-transcript unique read start positions, capped at `min_cov + 1`.
    /// Once exceeded, the set is drained and `exceeded_threshold[tx_idx]`
    /// is set instead.
    pub(crate) unique_starts: Vec<HashSet<u64, TinHashState>>,
    /// Per-transcript flag: true once unique start count exceeded `min_cov`.
    /// Avoids further HashSet inserts for high-coverage transcripts.
    pub exceeded_threshold: Vec<bool>,
    /// Per-transcript flag: true if any pileup read spans a sampled position
    /// (aligned, deleted or skipped). pysam then yields a pileup column there,
    /// so tin.py's coverage list is non-empty even when every count is zero.
    has_column: Vec<bool>,
    /// Number of sampled slots per transcript.
    #[allow(dead_code)]
    pub n_samples: Vec<u32>,
    /// Minimum MAPQ threshold (kept for API compatibility but not used;
    /// upstream tin.py does not filter by MAPQ).
    #[allow(dead_code)]
    pub mapq_cut: u8,
    /// Minimum coverage threshold (unique read start positions).
    pub min_cov: u32,
    /// Hash state used for unique_starts sets (deterministic when seeded).
    hash_state: TinHashState,
    /// Reusable scratch buffer for aligned blocks, avoids per-read allocation.
    blocks_buf: Vec<AlignedBlock>,
    /// Reusable scratch buffer for the current read's covered sampled bases.
    covered_buf: Vec<CoveredBase>,
    /// Reusable scratch buffer for [`overlap_partner`].
    gone_buf: Vec<bool>,
    /// Last chromosome seen -- used to detect chromosome transitions and
    /// reset the per-chromosome state below.
    cursor_chrom: String,
    /// Cursor into `chrom_spans` for the current chromosome: spans before it
    /// start at or before the last read start.
    cursor_span: usize,
    /// Transcripts whose span contains the last read start: (tx_end, tx_idx).
    active: Vec<(u64, u32)>,
    /// Cursor into `chrom_positions` for the current chromosome. Each read's
    /// alignment start is >= the previous read's, so we search forward from here.
    cursor_pos: usize,
    /// Start of the last pileup read (the position the pileups are at).
    plp_pos: Option<u64>,
    /// Largest fetch end of the pileup reads at `plp_pos` so far.
    plp_pos_reach: u64,
    /// Every pileup read that may still be buffered by some transcript's
    /// pileup: (end, fetch end, silent_below). An upper bound on any
    /// transcript's buffer, and the source for [`CappedPileup`]s.
    plp_buf: BinaryHeap<Reverse<(u64, u64, u64)>>,
    /// Transcripts whose pileup may hit the depth cap, by index.
    capped: BTreeMap<u32, CappedPileup>,
    /// Properly paired pileup reads whose hash slot may be in use, by name.
    mates: HashMap<Vec<u8>, MateGroup>,
    /// (max_end, name) of `mates` entries, to forget a group once the pileup
    /// has moved past it.
    mates_by_end: BinaryHeap<Reverse<(u64, Vec<u8>)>>,
}

impl TinAccum {
    /// Create a new TIN accumulator for the given index.
    ///
    /// When `seed` is `Some(n)`, the internal hash sets use a deterministic
    /// hash function seeded from `n`, making TIN results reproducible across
    /// runs. When `None`, the default random hash state is used.
    pub fn new(index: &TinIndex, mapq_cut: u8, min_cov: u32, seed: Option<u64>) -> Self {
        let n_transcripts = index.transcripts.len();
        let mut coverage = Vec::with_capacity(n_transcripts);
        let mut n_samples = Vec::with_capacity(n_transcripts);

        let hash_state = TinHashState::new(seed);

        for tx in &index.transcripts {
            let n = tx.sampled_positions.len();
            coverage.push(vec![0u32; n]);
            n_samples.push(n as u32);
        }

        let unique_starts = (0..n_transcripts)
            .map(|_| HashSet::with_hasher(hash_state.clone()))
            .collect();

        TinAccum {
            coverage,
            unique_starts,
            exceeded_threshold: vec![false; n_transcripts],
            has_column: vec![false; n_transcripts],
            n_samples,
            mapq_cut,
            min_cov,
            hash_state,
            blocks_buf: Vec::with_capacity(8),
            covered_buf: Vec::new(),
            gone_buf: Vec::new(),
            cursor_chrom: String::new(),
            cursor_span: 0,
            active: Vec::new(),
            cursor_pos: 0,
            plp_pos: None,
            plp_pos_reach: 0,
            plp_buf: BinaryHeap::new(),
            capped: BTreeMap::new(),
            mates: HashMap::new(),
            mates_by_end: BinaryHeap::new(),
        }
    }

    /// Process a single BAM record for TIN accumulation.
    ///
    /// Updates coverage at sampled positions that overlap the read's
    /// aligned blocks, and tracks read start positions.
    ///
    /// Matches upstream RSeQC `tin.py` filtering:
    /// - `check_min_reads()`: skips qcfail, unmapped, secondary only
    /// - `genebody_coverage()` pileup: additionally skips duplicates and
    ///   orphans, de-duplicates overlapping mates, skips bases with quality
    ///   below 13 and caps the depth (pysam `pileup()` defaults); deleted and
    ///   skipped (`N`) bases are not counted
    /// - Neither function filters by MAPQ or supplementary flags
    pub fn process_read(&mut self, record: &bam::Record, chrom_upper: &str, index: &TinIndex) {
        let flags = record.flags();

        // Skip unmapped, QC-fail, secondary only (matching upstream tin.py)
        if flags & 0x4 != 0 || flags & 0x200 != 0 || flags & 0x100 != 0 {
            return;
        }

        // Reset the per-chromosome state on chromosome transition.
        if chrom_upper != self.cursor_chrom {
            self.cursor_chrom.clear();
            self.cursor_chrom.push_str(chrom_upper);
            self.cursor_span = 0;
            self.active.clear();
            self.cursor_pos = 0;
            self.plp_pos = None;
            self.plp_pos_reach = 0;
            self.plp_buf.clear();
            self.capped.clear();
            self.mates.clear();
            self.mates_by_end.clear();
        }

        // Get chromosome positions for this chromosome
        let chrom_positions = match index.chrom_positions.get(chrom_upper) {
            Some(p) => p,
            None => return,
        };

        // Get transcript spans for this chromosome
        let chrom_spans = match index.chrom_spans.get(chrom_upper) {
            Some(s) => s,
            None => return,
        };

        // Upstream uses `aligned_read.pos` (the leftmost reference position).
        let read_start = record.pos().max(0) as u64;

        // Transcripts whose span contains the read start (tx_start <=
        // read_start < tx_end). Reads are coordinate-sorted, so spans only
        // enter (in tx_start order) and leave.
        let span_end = self.cursor_span
            + chrom_spans[self.cursor_span..].partition_point(|s| s.0 <= read_start);
        for &(_, tx_end, tx_idx) in &chrom_spans[self.cursor_span..span_end] {
            self.active.push((tx_end, tx_idx));
        }
        self.cursor_span = span_end;
        self.active.retain(|&(tx_end, _)| tx_end > read_start);

        for &(_, tx_idx) in &self.active {
            let idx = tx_idx as usize;
            // Skip HashSet insert if this transcript already exceeded the
            // min_cov threshold -- we only need the pass/fail decision.
            if !self.exceeded_threshold[idx] {
                self.unique_starts[idx].insert(read_start);
                if self.unique_starts[idx].len() > self.min_cov as usize {
                    self.exceeded_threshold[idx] = true;
                    // Free the HashSet memory; the flag is sufficient.
                    self.unique_starts[idx] = HashSet::with_hasher(self.hash_state.clone());
                }
            }
        }

        // Advance the position cursor to the first sampled position past the
        // read start (1-based P is in the read's span iff P > read_start).
        // Reads are coordinate-sorted, so this only moves forward.
        let base = self.cursor_pos;
        self.cursor_pos = base + chrom_positions[base..].partition_point(|p| p.0 <= read_start);

        // Reads that pysam's pileup never sees (default `stepper="samtools"`
        // flag_filter and `ignore_orphans=True`). check_min_reads() uses fetch()
        // and does count them, so this comes after the unique-start tracking.
        let is_paired = flags & 0x1 != 0;
        if flags & 0x400 != 0 || (is_paired && flags & 0x2 == 0) {
            return;
        }

        // Reference end as htslib's pileup computes it (raw CIGAR reference
        // length, may be 0), and as fetch() does (at least one base).
        let read_end = fill_aligned_blocks(record, &mut self.blocks_buf);
        let fetch_end = read_end.max(read_start + 1);

        // --- Pileup position: what every transcript's pileup has done before
        // this read is pushed.
        let first_at_pos = self.plp_pos != Some(read_start);
        let frontier = if first_at_pos {
            self.plp_pos.unwrap_or(0)
        } else {
            read_start
        };
        if first_at_pos {
            self.plp_pos = Some(read_start);
            self.plp_pos_reach = 0;
            while self.plp_buf.peek().is_some_and(|r| r.0 .0 < read_start) {
                self.plp_buf.pop();
            }
            let hot = self.plp_buf.len() + 1 > PILEUP_MAX_DEPTH;
            self.capped.retain(|_, c| {
                c.tx_end > read_start && (hot || c.dropped_until.is_some_and(|d| d >= read_start))
            });
        }
        while self.mates_by_end.peek().is_some_and(|e| e.0 .0 < frontier) {
            if let Some(Reverse((end, qname))) = self.mates_by_end.pop() {
                if self.mates.get(&qname).is_some_and(|g| g.max_end == end) {
                    self.mates.remove(&qname);
                }
            }
        }

        // --- Depth cap. A transcript's pileup holds a subset of `plp_buf`, so
        // no pileup can drop this read unless `plp_buf` is over the cap. Those
        // that could start tracking their own buffer here.
        if !first_at_pos && self.plp_buf.len() + 1 > PILEUP_MAX_DEPTH {
            let lo = chrom_spans.partition_point(|s| s.0 <= read_start);
            let hi = lo + chrom_spans[lo..].partition_point(|s| s.0 < fetch_end);
            let fresh: Vec<(u64, u64, u32)> = self
                .active
                .iter()
                .map(|&(tx_end, tx)| (index.transcripts[tx as usize].tx_start, tx_end, tx))
                .chain(chrom_spans[lo..hi].iter().copied())
                .filter(|&(_, _, tx)| !self.capped.contains_key(&tx))
                .collect();
            for (tx_start, tx_end, tx) in fresh {
                let ends = self
                    .plp_buf
                    .iter()
                    .filter(|r| r.0 .1 > tx_start && tx_start >= r.0 .2)
                    .map(|r| Reverse(r.0 .0))
                    .collect();
                let last_pos = (self.plp_pos_reach > tx_start).then_some(read_start);
                self.capped.insert(
                    tx,
                    CappedPileup {
                        tx_end,
                        last_pos,
                        ends,
                        dropped_until: None,
                    },
                );
            }
        }
        let mut dropped_in: Vec<u32> = Vec::new();
        for (&tx, c) in self.capped.iter_mut() {
            if index.transcripts[tx as usize].tx_start < fetch_end
                && c.push(read_start, read_end) == Push::Dropped
            {
                dropped_in.push(tx);
            }
        }
        let silent_below = if read_end == read_start {
            self.plp_pos_reach
        } else {
            0
        };
        self.plp_buf
            .push(Reverse((read_end, fetch_end, silent_below)));
        self.plp_pos_reach = self.plp_pos_reach.max(fetch_end);

        // --- Coverage. Walk the sampled positions inside the read's reference
        // span (read_start, read_end]. Every one gets a pileup column; those
        // inside an aligned block are covered by one of the read's bases.
        let kept_in = |tx: u32| dropped_in.is_empty() || dropped_in.binary_search(&tx).is_err();
        self.covered_buf.clear();
        let mut in_span = false;
        let mut block_idx = 0;
        let mut i = self.cursor_pos;
        while i < chrom_positions.len() && chrom_positions[i].0 <= read_end {
            in_span = true;
            let pos = chrom_positions[i].0;
            let mut j = i;
            while j < chrom_positions.len() && chrom_positions[j].0 == pos {
                if kept_in(chrom_positions[j].1) {
                    self.has_column[chrom_positions[j].1 as usize] = true;
                }
                j += 1;
            }

            // Aligned blocks are sorted; skip blocks that end before pos.
            let ref0 = pos - 1;
            while block_idx < self.blocks_buf.len() && self.blocks_buf[block_idx].ref_end <= ref0 {
                block_idx += 1;
            }
            if let Some(block) = self.blocks_buf.get(block_idx) {
                if block.ref_start <= ref0 {
                    self.covered_buf.push(CoveredBase {
                        first: i as u32,
                        last: j as u32,
                        qpos: (block.query_start + (ref0 - block.ref_start)) as u32,
                    });
                }
            }
            i = j;
        }

        // Count this read's bases with its own qualities; overlap tweaks
        // below correct them.
        let quals = record.qual();
        for b in &self.covered_buf {
            if quals.get(b.qpos as usize).copied().unwrap_or(0) >= PILEUP_MIN_BASE_QUALITY {
                for &(_, tx_idx, slot_idx) in &chrom_positions[b.first as usize..b.last as usize] {
                    if kept_in(tx_idx) {
                        self.coverage[tx_idx as usize][slot_idx as usize] += 1;
                    }
                }
            }
        }

        // --- Overlapping mates (`ignore_overlaps=True`). Only properly
        // paired reads enter htslib's overlap hash.
        if flags & 0x2 == 0 {
            return;
        }
        let mtid = record.mtid();
        let mpos = record.mpos();
        let rec = MateRec {
            pos: read_start,
            end: read_end,
            fetch_end,
            frontier,
            candidate: flags & 0x8 == 0
                && !(mtid >= 0 && mtid != record.tid())
                && !(record.insert_size().unsigned_abs() >= 2 * record.seq_len() as u64
                    && mpos >= read_end as i64),
            insertable: mpos >= record.pos() || (is_paired && mpos == -1),
            silent_below,
            dropped_in,
            data: None,
        };
        let qname = record.qname();
        if !self.mates.contains_key(qname) {
            self.mates.insert(qname.to_vec(), MateGroup::default());
        }
        let Some(group) = self.mates.get_mut(qname) else {
            return;
        };
        let could_pair = rec.candidate && group.recs.iter().any(|r| r.candidate && r.insertable);
        let store = rec.candidate && rec.insertable && in_span;
        if rec.end > group.max_end || group.recs.is_empty() {
            group.max_end = group.max_end.max(rec.end);
            self.mates_by_end
                .push(Reverse((group.max_end, qname.to_vec())));
        }
        group.recs.push(rec);

        if could_pair && in_span {
            let recs = &group.recs;
            let mut pairing = Pairing::new(recs, &mut self.gone_buf);
            pairing.resolve_overrides(recs, index, &mut self.gone_buf);
            for a in pairing.partners() {
                let Some(first) = recs[a].data.as_deref() else {
                    continue;
                };
                // Rewrite both reads' qualities as htslib does, then move each
                // affected base's count in the transcripts where they pair.
                let mut a_qual = first.qual.clone();
                let mut b_qual = quals.to_vec();
                tweak_overlap_quality(
                    &mut MateView {
                        pos: recs[a].pos as i64,
                        cigar: &first.cigar,
                        seq: &first.seq,
                        qual: &mut a_qual,
                    },
                    &mut MateView {
                        pos: record.pos(),
                        cigar: record.raw_cigar(),
                        seq: record.seq().encoded,
                        qual: &mut b_qual,
                    },
                    qname,
                );
                let passes = |q: &[u8], b: &CoveredBase| {
                    q.get(b.qpos as usize).copied().unwrap_or(0) >= PILEUP_MIN_BASE_QUALITY
                };
                for (bases, before, after) in [
                    (&first.bases, &first.qual[..], &a_qual[..]),
                    (&self.covered_buf, quals, &b_qual[..]),
                ] {
                    for b in bases {
                        let (was, now) = (passes(before, b), passes(after, b));
                        if was == now {
                            continue;
                        }
                        for &(_, tx_idx, slot_idx) in
                            &chrom_positions[b.first as usize..b.last as usize]
                        {
                            let tx_start = index.transcripts[tx_idx as usize].tx_start;
                            if pairing.partner(tx_idx, tx_start) != Some(a) {
                                continue;
                            }
                            let c = &mut self.coverage[tx_idx as usize][slot_idx as usize];
                            if now {
                                *c += 1;
                            } else {
                                *c -= 1;
                            }
                        }
                    }
                }
            }
        }

        // Keep what a later read of this name needs to tweak this one.
        if store {
            if let Some(last) = group.recs.last_mut() {
                last.data = Some(Box::new(MateData {
                    cigar: record.raw_cigar().to_vec(),
                    seq: record.seq().encoded.to_vec(),
                    qual: quals.to_vec(),
                    bases: self.covered_buf.clone(),
                }));
            }
        }
    }
    /// Merge another accumulator into this one.
    pub fn merge(&mut self, other: TinAccum) {
        for (i, other_cov) in other.coverage.into_iter().enumerate() {
            for (j, count) in other_cov.into_iter().enumerate() {
                self.coverage[i][j] += count;
            }
            self.has_column[i] |= other.has_column[i];
            // If either side already exceeded the threshold, the merged
            // result also exceeds it -- no need to union the sets.
            if self.exceeded_threshold[i] || other.exceeded_threshold[i] {
                self.exceeded_threshold[i] = true;
                self.unique_starts[i] = HashSet::with_hasher(self.hash_state.clone());
            } else {
                self.unique_starts[i].extend(other.unique_starts[i].iter());
                if self.unique_starts[i].len() > self.min_cov as usize {
                    self.exceeded_threshold[i] = true;
                    self.unique_starts[i] = HashSet::with_hasher(self.hash_state.clone());
                }
            }
        }
    }

    /// Compute TIN scores for all transcripts.
    pub fn into_result(self, index: &TinIndex) -> TinResults {
        let mut transcripts = Vec::with_capacity(index.transcripts.len());

        for (i, tx) in index.transcripts.iter().enumerate() {
            // A transcript passes if the exceeded_threshold flag was set
            // (unique starts > min_cov during accumulation or merge), or
            // if the remaining set still exceeds the threshold.
            let passed =
                self.exceeded_threshold[i] || self.unique_starts[i].len() as u32 > self.min_cov;
            let coverage = &self.coverage[i];

            // Upstream: `if len(read_count) > cutoff` -- strictly greater
            if !passed {
                transcripts.push(TinResult {
                    gene_id: tx.gene_id.clone(),
                    chrom: tx.chrom.clone(),
                    tx_start: tx.tx_start,
                    tx_end: tx.tx_end,
                    tin: 0.0,
                    passed_threshold: false,
                });
                continue;
            }

            // Use total_position_count (which includes duplicates from upstream's
            // chose_bases list) as the denominator, matching upstream's l = len(pick_positions)
            let n_total = tx.total_position_count;
            let tin = if self.has_column[i] {
                compute_tin(coverage, n_total)
            } else {
                // No pileup column at any sampled position: tin.py's coverage
                // list is empty and tin_score() returns 0.
                0.0
            };

            transcripts.push(TinResult {
                gene_id: tx.gene_id.clone(),
                chrom: tx.chrom.clone(),
                tx_start: tx.tx_start,
                tx_end: tx.tx_end,
                tin,
                passed_threshold: true,
            });
        }

        TinResults { transcripts }
    }
}

/// Compute TIN score from coverage vector using Shannon entropy.
///
/// Matches upstream RSeQC's `tin_score()` for a non-empty coverage list:
///   cvg_eff = [c for c in cvg if c > 0]
///   entropy = shannon_entropy(cvg_eff)
///   tin = 100 * exp(entropy) / l
///
/// where `l` is the total number of sampled positions (including zeros),
/// NOT just the nonzero count. A transcript with many uncovered positions
/// gets a lower TIN because the denominator is larger. When no position has
/// coverage, the entropy is 0 and the TIN is `100 / l`, as upstream. (The
/// caller handles the empty-list case, where upstream returns 0.)
fn compute_tin(coverage: &[u32], n_total_positions: usize) -> f64 {
    if n_total_positions == 0 {
        return 0.0;
    }

    let nonzero: Vec<f64> = coverage
        .iter()
        .filter(|&&c| c > 0)
        .map(|&c| c as f64)
        .collect();

    let total: f64 = nonzero.iter().sum();

    // Shannon entropy over non-zero positions only (matching upstream)
    let entropy: f64 = if total == 0.0 {
        0.0
    } else {
        nonzero
            .iter()
            .map(|&c| {
                let p = c / total;
                -p * p.ln()
            })
            .sum()
    };

    // Upstream: tin = 100 * exp(entropy) / l
    // where l = total number of sampled positions (not just nonzero)
    100.0 * entropy.exp() / n_total_positions as f64
}

/// An aligned (M/=/X) block of a read.
#[derive(Debug, Clone, Copy)]
struct AlignedBlock {
    /// Reference start (0-based, inclusive)
    ref_start: u64,
    /// Reference end (0-based, exclusive)
    ref_end: u64,
    /// Query offset of the block's first base
    query_start: u64,
}

/// Fill `buf` with the aligned blocks of a BAM record's CIGAR, reusing the
/// existing Vec capacity to avoid per-read heap allocation.
///
/// # Returns
/// The reference end of the alignment (0-based, exclusive), including
/// deletions and skipped regions.
fn fill_aligned_blocks(record: &bam::Record, buf: &mut Vec<AlignedBlock>) -> u64 {
    use rust_htslib::bam::record::Cigar;

    buf.clear();
    let mut pos = record.pos().max(0) as u64;
    let mut qpos = 0u64;

    for op in record.cigar().iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                buf.push(AlignedBlock {
                    ref_start: pos,
                    ref_end: pos + *len as u64,
                    query_start: qpos,
                });
                pos += *len as u64;
                qpos += *len as u64;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                pos += *len as u64;
            }
            Cigar::Ins(len) | Cigar::SoftClip(len) => {
                qpos += *len as u64;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    pos
}

// ===================================================================
// Output writers
// ===================================================================

/// Format a finite TIN value as Python formats a float.
fn format_tin_float(value: f64) -> String {
    let mut text = format!("{value:?}");
    // Python pads a single-digit negative exponent, for example 1e-05.
    if text
        .split_once("e-")
        .is_some_and(|(_, exponent)| exponent.len() == 1)
    {
        text.insert(text.len() - 1, '0');
    }
    text
}

/// Write per-transcript TIN scores in RSeQC format.
///
/// Format: TSV with header geneID\tchrom\ttx_start\ttx_end\tTIN
pub fn write_tin(results: &TinResults, output_path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(output_path)
        .with_context(|| format!("Failed to create TIN output: {}", output_path.display()))?;

    writeln!(f, "geneID\tchrom\ttx_start\ttx_end\tTIN")?;

    for r in &results.transcripts {
        // Match upstream formatting:
        // - Failed check_min_reads: explicitly writes "0.0" (Python float literal)
        // - Passed but tin_score returned int 0 (empty coverage): writes "0"
        // - Non-zero TIN: writes float value
        if !r.passed_threshold {
            // Below coverage threshold → upstream writes explicit 0.0
            writeln!(
                f,
                "{}\t{}\t{}\t{}\t0.0",
                r.gene_id, r.chrom, r.tx_start, r.tx_end
            )?;
        } else if r.tin == 0.0 {
            // Passed threshold but zero TIN → upstream tin_score returns int 0
            writeln!(
                f,
                "{}\t{}\t{}\t{}\t0",
                r.gene_id, r.chrom, r.tx_start, r.tx_end
            )?;
        } else {
            writeln!(
                f,
                "{}\t{}\t{}\t{}\t{}",
                r.gene_id,
                r.chrom,
                r.tx_start,
                r.tx_end,
                format_tin_float(r.tin)
            )?;
        }
    }

    Ok(())
}

/// Sum `values` the way NumPy's `add.reduce` does for a contiguous float64
/// array (`pairwise_sum` in NumPy's `loops_utils.h.src`).
///
/// tin.py reports `numpy.mean` and `numpy.std` of its TIN list, so the
/// summary must use the same summation order to print the same digits:
/// blocks of up to 128 values are summed with eight interleaved accumulators,
/// and longer runs are split in two (at a multiple of 8) and summed recursively.
fn numpy_pairwise_sum(values: &[f64]) -> f64 {
    /// NumPy's `PW_BLOCKSIZE`.
    const BLOCK: usize = 128;
    let n = values.len();
    if n < 8 {
        values.iter().fold(-0.0, |acc, &v| acc + v)
    } else if n <= BLOCK {
        let mut r = [0.0f64; 8];
        r.copy_from_slice(&values[..8]);
        let unrolled_end = n - n % 8;
        let (chunks, _) = values[8..unrolled_end].as_chunks::<8>();
        for chunk in chunks {
            for (acc, &v) in r.iter_mut().zip(chunk) {
                *acc += v;
            }
        }
        let mut res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for &v in &values[unrolled_end..] {
            res += v;
        }
        res
    } else {
        let mut n2 = n / 2;
        n2 -= n2 % 8;
        numpy_pairwise_sum(&values[..n2]) + numpy_pairwise_sum(&values[n2..])
    }
}

/// Write TIN summary statistics.
///
/// Format: Bam_file\tTIN(mean)\tTIN(median)\tTIN(stdev)
///
/// Matches tin.py: the BAM is reported by its file name, and the mean and
/// (population) standard deviation are NumPy's, over the TINs of the
/// transcripts that passed the read check, in output row order.
pub fn write_tin_summary(results: &TinResults, bam_path: &str, output_path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(output_path)
        .with_context(|| format!("Failed to create TIN summary: {}", output_path.display()))?;

    let scores: Vec<f64> = results
        .transcripts
        .iter()
        .filter(|r| r.passed_threshold)
        .map(|r| r.tin)
        .collect();

    // tin.py prints os.path.basename() of the BAM path.
    let bam_name = Path::new(bam_path).file_name().map_or_else(
        || bam_path.to_string(),
        |n| n.to_string_lossy().into_owned(),
    );

    writeln!(f, "Bam_file\tTIN(mean)\tTIN(median)\tTIN(stdev)")?;
    // NumPy reports nan when no transcript passed the read check.
    if scores.is_empty() {
        writeln!(f, "{bam_name}\tnan\tnan\tnan")?;
        return Ok(());
    }

    let n = scores.len() as f64;
    // numpy.mean: add.reduce (initial 0.0) divided by the count.
    let mean = (0.0 + numpy_pairwise_sum(&scores)) / n;
    let median = crate::io::median(&scores);

    // numpy.std: sqrt of the mean squared deviation from numpy.mean.
    let squared: Vec<f64> = scores.iter().map(|x| (x - mean) * (x - mean)).collect();
    let stdev = ((0.0 + numpy_pairwise_sum(&squared)) / n).sqrt();

    writeln!(
        f,
        "{}\t{}\t{}\t{}",
        bam_name,
        format_tin_float(mean),
        format_tin_float(median),
        format_tin_float(stdev)
    )?;

    Ok(())
}

// ===================================================================
// Unit tests
// ===================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn gene_with(gene_id: &str, chrom: &str, txs: Vec<(&str, Vec<(u64, u64)>, u64)>) -> Gene {
        let transcripts: Vec<Transcript> = txs
            .into_iter()
            .map(|(id, exons, first_exon_start)| Transcript {
                transcript_id: id.to_string(),
                has_transcript_id: true,
                chrom: chrom.to_string(),
                start: exons[0].0,
                end: exons[exons.len() - 1].1,
                strand: '-',
                exons,
                first_exon_start,
                cds_start: None,
                cds_end: None,
            })
            .collect();
        Gene {
            gene_id: gene_id.to_string(),
            chrom: chrom.to_string(),
            start: transcripts.iter().map(|t| t.start).min().unwrap(),
            end: transcripts.iter().map(|t| t.end).max().unwrap(),
            strand: '-',
            exons: Vec::new(),
            effective_length: 0,
            attributes: HashMap::new(),
            transcripts,
        }
    }

    #[test]
    fn test_tin_summary_with_no_passing_transcripts() {
        let results = TinResults {
            transcripts: vec![TinResult {
                gene_id: "tx".to_owned(),
                chrom: "chr1".to_owned(),
                tx_start: 0,
                tx_end: 100,
                tin: 0.0,
                passed_threshold: false,
            }],
        };
        let file =
            std::env::temp_dir().join(format!("rustqc-tin-summary-{}.txt", std::process::id()));
        write_tin_summary(&results, "staged/test.bam", &file).unwrap();
        let summary = std::fs::read_to_string(&file).unwrap();
        std::fs::remove_file(file).unwrap();
        assert_eq!(
            summary,
            "Bam_file\tTIN(mean)\tTIN(median)\tTIN(stdev)\ntest.bam\tnan\tnan\tnan\n"
        );
    }

    #[test]
    fn test_tin_float_formatting() {
        assert_eq!(format_tin_float(20.0), "20.0");
        assert_eq!(format_tin_float(0.0), "0.0");
        assert_eq!(format_tin_float(1e-5), "1e-05");
    }

    #[test]
    fn test_numpy_pairwise_sum() {
        // Short inputs are summed left to right.
        assert_eq!(numpy_pairwise_sum(&[]), -0.0);
        assert_eq!(numpy_pairwise_sum(&[1.0, 2.0, 3.0]), 6.0);
        // The pairwise order shows in the last bits. Reference values are
        // numpy.add.reduce (NumPy 2.5.5) of 1/k.
        let v: Vec<f64> = (1..=300).map(|k| 1.0 / k as f64).collect();
        assert_eq!(numpy_pairwise_sum(&v[..20]), 3.597739657143682);
        assert_eq!(numpy_pairwise_sum(&v), 6.282663880299504);
        assert_eq!(v.iter().sum::<f64>(), 6.282663880299502);
    }

    #[test]
    fn test_index_follows_gtf2bed_order() {
        // gtf2bed sorts by chromosome name, then by the first exon line's
        // start, then by transcript ID.
        let mut genes = IndexMap::new();
        genes.insert(
            "G1".to_string(),
            gene_with(
                "G1",
                "chr2",
                vec![
                    ("Z", vec![(100, 200), (900, 1000)], 900),
                    ("B", vec![(300, 400)], 300),
                ],
            ),
        );
        genes.insert(
            "G2".to_string(),
            gene_with(
                "G2",
                "chr10",
                vec![
                    ("C", vec![(5000, 5100)], 5000),
                    ("D", vec![(500, 600)], 500),
                ],
            ),
        );
        genes.insert(
            "G3".to_string(),
            gene_with("G3", "chr2", vec![("E", vec![(900, 950)], 900)]),
        );
        // No transcript_id: gtf2bed writes no BED row, so tin.py has no row.
        let mut no_tx = gene_with("G4", "chr2", vec![("G4__no_tx", vec![(50, 60)], 50)]);
        no_tx.transcripts[0].has_transcript_id = false;
        genes.insert("G4".to_string(), no_tx);
        let index = TinIndex::from_genes(&genes, 100);
        let order: Vec<&str> = index
            .transcripts
            .iter()
            .map(|t| t.gene_id.as_str())
            .collect();
        assert_eq!(order, vec!["D", "C", "B", "E", "Z"]);
    }

    #[test]
    fn test_sample_exonic_positions() {
        // Single exon of 10 bases: 0-based [100, 110)
        // tx_start=100, tx_end=110 (matching BED format)
        // mRNA_size = 10, sample_size = 5, so mRNA_size > sample_size
        // step_size = 10 / 5 = 2
        // gene_all_base (1-based) = [101, 102, ..., 110]
        // exon_bounds = [101, 110]
        // strided: indices 0,2,4,6,8 -> [101, 103, 105, 107, 109]
        // merged with bounds [101, 110] -> [101, 103, 105, 107, 109, 110]
        // Large mRNA case: deduped by uniqify, total_count = unique count
        let exons = vec![(100, 110)];
        let (positions, total_count) = sample_exonic_positions(&exons, 5, 100, 110);
        assert_eq!(positions.len(), 6); // 5 strided + exon end boundary 110
        assert_eq!(total_count, 6); // no duplicates in large mRNA case
                                    // Positions should be 1-based
        assert_eq!(positions[0], 101);
        assert_eq!(positions[5], 110);
    }

    #[test]
    fn test_sample_multi_exon() {
        // Two exons: 0-based [100,105) and [200,205) → 10 bases total
        // tx_start=100, tx_end=205
        // mRNA_size = 10, sample_size = 10, so mRNA_size <= sample_size
        // chose_bases = [tx_start+1, tx_end] + all exonic = [101, 205, 101..=105, 201..=205]
        // Sorted with dupes: [101, 101, 102, 103, 104, 105, 201, 202, 203, 204, 205, 205]
        // total_count = 12 (with duplicates), unique positions = 10
        let exons = vec![(100, 105), (200, 205)];
        let (positions, total_count) = sample_exonic_positions(&exons, 10, 100, 205);
        assert_eq!(positions.len(), 10); // 10 unique positions
        assert_eq!(total_count, 12); // 10 + 2 boundary duplicates
        assert_eq!(positions[0], 101);
        assert_eq!(positions[9], 205);
    }

    #[test]
    fn test_sample_tiny_exon() {
        // Exon smaller than sample size: 0-based [100, 103) → 3 bases
        // tx_start=100, tx_end=103
        // mRNA_size = 3 <= sample_size = 100
        // chose_bases = [101, 103, 101, 102, 103] -> sorted: [101, 101, 102, 103, 103]
        // total_count = 5 (with duplicates), unique positions = 3
        let exons = vec![(100, 103)];
        let (positions, total_count) = sample_exonic_positions(&exons, 100, 100, 103);
        assert_eq!(positions.len(), 3); // 3 unique exonic positions
        assert_eq!(total_count, 5); // 3 + 2 boundary duplicates
        assert_eq!(positions[0], 101);
        assert_eq!(positions[2], 103);
    }

    #[test]
    fn test_compute_tin_uniform() {
        // Uniform coverage → max TIN
        let cov = vec![10, 10, 10, 10, 10, 10, 10, 10, 10, 10];
        let tin = compute_tin(&cov, cov.len());
        assert!(
            (tin - 100.0).abs() < 0.01,
            "Uniform coverage TIN should be ~100: {tin}"
        );
    }

    #[test]
    fn test_compute_tin_uniform_with_zeros() {
        // Uniform coverage but only 5 of 10 positions covered
        // TIN = 100 * exp(H) / 10; H entropy over 5 equal values -> exp(H)=5
        // TIN = 100 * 5 / 10 = 50
        let cov = vec![10, 10, 10, 10, 10, 0, 0, 0, 0, 0];
        let tin = compute_tin(&cov, cov.len());
        assert!(
            (tin - 50.0).abs() < 0.01,
            "Half-covered uniform TIN should be ~50: {tin}"
        );
    }

    #[test]
    fn test_compute_tin_degraded() {
        // 5' degradation: high at start, low at end
        let cov = vec![100, 80, 60, 40, 20, 10, 5, 2, 1, 1];
        let tin = compute_tin(&cov, cov.len());
        assert!(
            tin > 0.0 && tin < 100.0,
            "Degraded TIN should be between 0 and 100: {tin}"
        );
    }

    /// Pack an ACGT string into BAM 4-bit encoding.
    fn pack_seq(seq: &str) -> Vec<u8> {
        let code = |c: u8| match c {
            b'A' => 1u8,
            b'C' => 2,
            b'G' => 4,
            b'T' => 8,
            _ => 15,
        };
        seq.as_bytes()
            .chunks(2)
            .map(|p| (code(p[0]) << 4) | p.get(1).map_or(0, |&c| code(c)))
            .collect()
    }

    /// Raw BAM CIGAR operation.
    fn cig(op: u32, len: u32) -> u32 {
        (len << 4) | op
    }

    /// Run `tweak_overlap_quality` on two mates and return their new qualities.
    #[allow(clippy::too_many_arguments)]
    fn tweak(
        a_pos: i64,
        a_cigar: &[u32],
        a_seq: &str,
        a_q: u8,
        b_pos: i64,
        b_cigar: &[u32],
        b_seq: &str,
        b_q: u8,
        qname: &[u8],
    ) -> (Vec<u8>, Vec<u8>) {
        let (a_packed, b_packed) = (pack_seq(a_seq), pack_seq(b_seq));
        let mut a_qual = vec![a_q; a_seq.len()];
        let mut b_qual = vec![b_q; b_seq.len()];
        tweak_overlap_quality(
            &mut MateView {
                pos: a_pos,
                cigar: a_cigar,
                seq: &a_packed,
                qual: &mut a_qual,
            },
            &mut MateView {
                pos: b_pos,
                cigar: b_cigar,
                seq: &b_packed,
                qual: &mut b_qual,
            },
            qname,
        );
        (a_qual, b_qual)
    }

    #[test]
    fn test_htslib_name_hash_parity() {
        // Which mate keeps the summed quality of matching overlapping bases, as
        // observed with pysam 0.24.1 (htslib 1.24): odd hash = first mate.
        for (name, first_keeps) in [
            ("pair1", false),
            ("pair2", true),
            ("pair3", false),
            ("pair4", false),
            ("SRR5665260.26358182", true),
            ("read/1", true),
        ] {
            assert_eq!(
                htslib_name_hash(name.as_bytes()) & 1 == 1,
                first_keeps,
                "hash parity for {name}"
            );
        }
    }

    #[test]
    fn test_tweak_overlap_matching_and_mismatching_bases() {
        let m10 = [cig(CIGAR_MATCH, 10)];
        // Matching bases (Q10 + Q10): one mate gets 20, the other 0, at the 5
        // overlapping positions; bases outside the overlap are untouched.
        let (a, b) = tweak(
            200,
            &m10,
            "CCCCCCCCCC",
            10,
            205,
            &m10,
            "CCCCCCCCCC",
            10,
            b"x",
        );
        for i in 0..5 {
            assert_eq!(a[i], 10, "a outside overlap");
            assert_eq!(b[5 + i], 10, "b outside overlap");
            assert_eq!(a[5 + i].max(b[i]), 20);
            assert_eq!(a[5 + i].min(b[i]), 0);
        }
        // Mismatching bases: the higher quality keeps 0.8 * q (truncated)
        let (a, b) = tweak(
            300,
            &m10,
            "GGGGGGGGGG",
            17,
            305,
            &m10,
            "TTTTTTTTTT",
            16,
            b"x",
        );
        assert_eq!((a[5], b[0]), (13, 0));
        let (a, b) = tweak(
            300,
            &m10,
            "GGGGGGGGGG",
            16,
            305,
            &m10,
            "TTTTTTTTTT",
            16,
            b"x",
        );
        assert_eq!(a[5].max(b[0]), 12);
        assert_eq!(a[5].min(b[0]), 0);
    }

    #[test]
    fn test_tweak_overlap_after_ref_skip_and_deletion() {
        // First mate 5M100xN/D50M at 1000 (second block starts at 1105), second
        // mate 20M at 1104. Observed with pysam 0.24.1: after N the first base of
        // the block keeps both qualities; after D the walk catches up instead.
        let a_seq = "A".repeat(55);
        let b_seq = "A".repeat(20);
        let b_cigar = [cig(CIGAR_MATCH, 20)];

        let a_cigar = [
            cig(CIGAR_MATCH, 5),
            cig(CIGAR_REF_SKIP, 100),
            cig(CIGAR_MATCH, 50),
        ];
        let (a, b) = tweak(1000, &a_cigar, &a_seq, 30, 1104, &b_cigar, &b_seq, 30, b"x");
        // Ref 1105 = a[5], b[1]: untouched
        assert_eq!((a[5], b[1]), (30, 30));
        // Ref 1106 = a[6], b[2]: de-duplicated
        assert_eq!(a[6].min(b[2]), 0);
        assert_eq!(a[6].max(b[2]), 60);

        let a_cigar = [
            cig(CIGAR_MATCH, 5),
            cig(CIGAR_DEL, 100),
            cig(CIGAR_MATCH, 50),
        ];
        let (a, b) = tweak(1000, &a_cigar, &a_seq, 30, 1104, &b_cigar, &b_seq, 30, b"x");
        // Ref 1104 is inside a's deletion: b's base is scaled or zeroed
        assert!(b[0] == 24 || b[0] == 0, "b in deletion: {}", b[0]);
        // Ref 1105 is de-duplicated
        assert_eq!(a[5].min(b[1]), 0);
        assert_eq!(a[5].max(b[1]), 60);
    }

    #[test]
    fn test_compute_tin_all_zero() {
        // tin.py: a non-empty coverage list with no positive value has entropy 0,
        // so TIN = 100 * exp(0) / l = 100 / 5
        let cov = vec![0, 0, 0, 0, 0];
        let tin = compute_tin(&cov, cov.len());
        assert!((tin - 20.0).abs() < 1e-12, "All-zero coverage TIN: {tin}");
    }

    #[test]
    fn test_compute_tin_single_nonzero() {
        // Only one position covered out of 5 -> entropy=0 -> exp(0)=1
        // TIN = 100 * 1 / 5 = 20
        let cov = vec![0, 0, 100, 0, 0];
        let tin = compute_tin(&cov, cov.len());
        assert!(
            (tin - 20.0).abs() < 0.01,
            "Single position coverage should give TIN=20: {tin}"
        );
    }
}
