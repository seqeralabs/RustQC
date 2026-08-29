//! Per-read accumulator structs for single-pass RSeQC integration.
//!
//! Each RSeQC tool has an accumulator that collects per-read data during the
//! main BAM counting loop. Accumulators are created per chromosome worker and
//! merged after parallel processing, just like `ChromResult`.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::hash::{Hash, Hasher};

use anyhow::Result;
use indexmap::IndexMap;
use rust_htslib::bam;

// BamStatAccum is read-level and assay-agnostic; it lives in `crate::common`
// and is shared with the dna pipeline. Re-exported so existing paths resolve.
use super::common::{self, KnownJunctionSet, ReferenceJunctions};
use super::infer_experiment::{GeneModel, InferExperimentResult};
use super::inner_distance::{
    build_histogram, ExonBitset, InnerDistanceResult, PairRecord, TranscriptTree,
};
use super::junction_annotation::{Junction, JunctionClass, JunctionResults};
use super::junction_saturation::SaturationResult;
use super::read_distribution::{ChromIntervals, ReadDistributionResult, RegionSets};
use super::read_duplication::ReadDuplicationResult;
use super::tin::TinAccum;
pub use crate::common::bam_stat_accum::BamStatAccum;
use crate::rna::preseq::PreseqAccum;

use crate::rna::bam_flags::*;

// ===================================================================
// Shared references to annotation data
// ===================================================================

/// Read-only annotation data shared across all chromosome workers.
///
/// Each field is `Option` — `None` when the corresponding tool is disabled.
pub struct RseqcAnnotations<'a> {
    /// Gene model for infer_experiment.
    pub gene_model: Option<&'a GeneModel>,
    /// Reference junctions for junction_annotation.
    pub ref_junctions: Option<&'a ReferenceJunctions>,
    /// Genomic region sets for read_distribution.
    pub rd_regions: Option<&'a RegionSets>,
    /// Exon bitset for inner_distance.
    pub exon_bitset: Option<&'a ExonBitset>,
    /// Transcript tree for inner_distance.
    pub transcript_tree: Option<&'a TranscriptTree>,

    /// TIN index for transcript integrity number calculation.
    pub tin_index: Option<&'a super::tin::TinIndex>,
}

/// Per-tool configuration parameters.
///
/// Some fields are only used at accumulator construction time or during
/// result conversion, not during per-read processing.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct RseqcConfig {
    /// Strands to build coverage track accumulators for. Empty disables them;
    /// a single `None` entry counts every read into one combined track.
    pub coverage_strands: Vec<Option<char>>,
    /// MAPQ cutoff for read quality filtering.
    pub mapq_cut: u8,
    /// Maximum reads to sample for infer_experiment.
    pub infer_experiment_sample_size: u64,
    /// Minimum intron size for junction filtering.
    pub min_intron: u64,
    /// Minimum coverage for junction saturation.
    pub junction_saturation_min_coverage: u32,
    /// Sampling start percentage for junction saturation.
    pub junction_saturation_sample_start: u32,
    /// Sampling end percentage for junction saturation.
    pub junction_saturation_sample_end: u32,
    /// Sampling step percentage for junction saturation.
    pub junction_saturation_sample_step: u32,
    /// Maximum read pairs to sample for inner distance.
    pub inner_distance_sample_size: u64,
    /// Lower bound of inner distance histogram.
    pub inner_distance_lower_bound: i64,
    /// Upper bound of inner distance histogram.
    pub inner_distance_upper_bound: i64,
    /// Bin width for inner distance histogram.
    pub inner_distance_step: i64,
    /// Whether bam_stat analysis is enabled.
    pub bam_stat_enabled: bool,
    /// Whether infer_experiment analysis is enabled.
    pub infer_experiment_enabled: bool,
    /// Whether read_duplication analysis is enabled.
    pub read_duplication_enabled: bool,
    /// Whether read_distribution analysis is enabled.
    pub read_distribution_enabled: bool,
    /// Whether junction_annotation analysis is enabled.
    pub junction_annotation_enabled: bool,
    /// Whether junction_saturation analysis is enabled.
    pub junction_saturation_enabled: bool,
    /// Whether inner_distance analysis is enabled.
    pub inner_distance_enabled: bool,
    /// Whether TIN analysis is enabled.
    pub tin_enabled: bool,
    /// Number of equally-spaced sampling positions per transcript for TIN.
    pub tin_sample_size: usize,
    /// Minimum number of read starts for a transcript to compute TIN.
    pub tin_min_coverage: u32,
    /// Random seed for reproducible TIN results (deterministic hash state).
    pub tin_seed: Option<u64>,
    /// Random seed for junction_saturation observation shuffle.
    /// Defaults to 42 when not set via `--seed`.
    pub junction_saturation_seed: u64,
    /// Whether preseq library complexity estimation is enabled.
    pub preseq_enabled: bool,
    /// Maximum merged PE fragment length for preseq (preseq's `-seg_len`).
    pub preseq_max_segment_length: i64,
}

// ===================================================================
// Per-tool accumulators
// ===================================================================

// -------------------------------------------------------------------
// infer_experiment accumulator
// -------------------------------------------------------------------

/// infer_experiment accumulator — strand protocol inference.
///
/// Processes ALL reads (no sampling cap). Upstream RSeQC `infer_experiment.py`
/// samples 200K reads, but since RustQC processes in a single pass, processing
/// all reads gives equivalent-or-better results on small data and avoids
/// sampling divergence on large data.
#[derive(Debug, Default)]
pub struct InferExpAccum {
    /// Paired-end strand class counts (keys: "1++", "1--", "2+-", "2-+", etc.)
    pub p_strandness: HashMap<String, u64>,
    /// Single-end strand class counts (keys: "++", "--", "+-", "-+", etc.)
    pub s_strandness: HashMap<String, u64>,
    /// Reusable key buffer to avoid per-read format!() allocations.
    key_buf: String,
}

impl InferExpAccum {
    /// Process a single BAM record.
    pub fn process_read(
        &mut self,
        record: &bam::Record,
        chrom: &str,
        model: &GeneModel,
        mapq_cut: u8,
    ) {
        let flags = record.flags();

        // Skip QC-fail, dup, secondary, unmapped.
        // Upstream RSeQC does not filter supplementary (0x800), so we don't
        // either — this matches the upstream filter set exactly.
        if flags & BAM_FQCFAIL != 0
            || flags & BAM_FDUP != 0
            || flags & BAM_FSECONDARY != 0
            || flags & BAM_FUNMAP != 0
        {
            return;
        }

        if record.mapq() < mapq_cut {
            return;
        }

        let map_strand = if record.is_reverse() { '-' } else { '+' };

        // Compute query alignment length (M+I+=+X, NO soft clips) to match
        // upstream RSeQC's `qlen` which is pysam's `query_alignment_length`.
        // pos + qlen gives an approximate reference end (slightly overestimated
        // due to insertions, matching upstream's approximation). Soft clips must
        // be excluded because they don't consume the reference — including them
        // would widen the query interval, causing more spurious overlaps with
        // genes on both strands and inflating the "failed to determine" fraction.
        let read_start = record.pos() as u64;
        let qalen: u64 = record
            .cigar()
            .iter()
            .filter_map(|op| {
                use rust_htslib::bam::record::Cigar::*;
                match op {
                    Match(len) | Ins(len) | Equal(len) | Diff(len) => Some(*len as u64),
                    _ => None,
                }
            })
            .sum();
        let read_end = read_start + qalen;

        let strands = model.find_strands(chrom, read_start, read_end);
        if strands.is_empty() {
            return;
        }

        // Build key in reusable buffer to avoid per-read allocation.
        // Keys are small (e.g. "1++", "2+-", "++:-") with ~12 distinct values.
        self.key_buf.clear();
        let map = if record.is_paired() {
            self.key_buf.push(if record.is_first_in_template() {
                '1'
            } else {
                '2'
            });
            &mut self.p_strandness
        } else {
            &mut self.s_strandness
        };
        self.key_buf.push(map_strand);
        for (i, &s) in strands.iter().enumerate() {
            if i > 0 {
                self.key_buf.push(':');
            }
            self.key_buf.push(s as char);
        }
        // Only allocate a new String for first-seen keys
        match map.get_mut(self.key_buf.as_str()) {
            Some(count) => *count += 1,
            None => {
                map.insert(self.key_buf.clone(), 1);
            }
        }
    }

    /// Merge another accumulator into this one.
    ///
    /// Raw counts are accumulated without scaling. `into_result()` computes
    /// fractions from raw counts — the absolute count doesn't matter, only
    /// the proportions.
    pub fn merge(&mut self, other: InferExpAccum) {
        for (key, count) in other.p_strandness {
            *self.p_strandness.entry(key).or_insert(0) += count;
        }
        for (key, count) in other.s_strandness {
            *self.s_strandness.entry(key).or_insert(0) += count;
        }
    }
}

// -------------------------------------------------------------------
// read_duplication accumulator
// -------------------------------------------------------------------

/// read_duplication accumulator — sequence-based and position-based dedup.
#[derive(Debug, Default)]
pub struct ReadDupAccum {
    /// Sequence hash → occurrence count (hash-based dedup to save memory).
    pub seq_dup: HashMap<u128, u64>,
    /// Position key hash → occurrence count (hash-based dedup to save memory).
    pub pos_dup: HashMap<u64, u64>,
}

impl ReadDupAccum {
    /// Process a single BAM record.
    pub fn process_read(&mut self, record: &bam::Record, chrom: &str, mapq_cut: u8) {
        let flags = record.flags();

        // Filter: unmapped, QC-fail. Does NOT skip dup or secondary (intentional).
        if flags & BAM_FUNMAP != 0 || flags & BAM_FQCFAIL != 0 {
            return;
        }
        if record.mapq() < mapq_cut {
            return;
        }

        // Sequence-based: hash directly from BAM 4-bit encoding (no allocation)
        let seq_hash = hash_sequence_encoded(&record.seq());
        *self.seq_dup.entry(seq_hash).or_insert(0) += 1;

        // Position-based: hash key from CIGAR (avoids string allocation)
        let pos = record.pos();
        let cigar = record.cigar();
        let key = hash_position_key(chrom, pos, &cigar);
        *self.pos_dup.entry(key).or_insert(0) += 1;
    }

    /// Merge another accumulator into this one.
    pub fn merge(&mut self, other: ReadDupAccum) {
        for (hash, count) in other.seq_dup {
            *self.seq_dup.entry(hash).or_insert(0) += count;
        }
        for (key, count) in other.pos_dup {
            *self.pos_dup.entry(key).or_insert(0) += count;
        }
    }
}

/// Hash a read sequence to u128 for deduplication.
///
/// Reads directly from the BAM 4-bit encoded nucleotide representation,
/// avoiding the allocation that `seq.as_bytes()` would require. Each base
/// is already case-insensitive in 4-bit encoding, so no uppercasing is needed.
///
/// Uses two rounds of SipHash-1-3 (via `DefaultHasher`) to produce a 128-bit
/// fingerprint. Effective collision resistance is ~64 bits (birthday bound ~2^32),
/// more than sufficient for typical RNA-seq datasets (< 1B distinct reads).
fn hash_sequence_encoded(seq: &bam::record::Seq<'_>) -> u128 {
    let len = seq.len();
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    for i in 0..len {
        // encoded_base returns a 4-bit IUPAC code (0-15), inherently case-insensitive
        seq.encoded_base(i).hash(&mut hasher);
    }
    // DefaultHasher produces u64; extend to u128 by double-hashing with length
    let h1 = hasher.finish();
    let mut hasher2 = std::collections::hash_map::DefaultHasher::new();
    len.hash(&mut hasher2);
    h1.hash(&mut hasher2);
    let h2 = hasher2.finish();
    (h1 as u128) << 64 | (h2 as u128)
}

/// Hash position key matching RSeQC's `fetch_exon` + position key logic.
/// Uses FNV-1a hashing to avoid string allocation per read.
fn hash_position_key(chrom: &str, pos: i64, cigar: &bam::record::CigarStringView) -> u64 {
    use rust_htslib::bam::record::Cigar;

    let mut h = crate::io::FNV1A_OFFSET;
    crate::io::fnv1a_update(&mut h, chrom.as_bytes());
    crate::io::fnv1a_update(&mut h, &pos.to_le_bytes());

    let mut ref_pos = pos;
    for op in cigar.iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let end = ref_pos + *len as i64;
                crate::io::fnv1a_update(&mut h, &ref_pos.to_le_bytes());
                crate::io::fnv1a_update(&mut h, &end.to_le_bytes());
                ref_pos = end;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                ref_pos += *len as i64;
            }
            Cigar::SoftClip(len) => {
                // RSeQC bug: S advances reference position
                ref_pos += *len as i64;
            }
            Cigar::Ins(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    h
}

// -------------------------------------------------------------------
// read_distribution accumulator
// -------------------------------------------------------------------

/// read_distribution accumulator — region classification counters.
#[derive(Debug, Default)]
pub struct ReadDistAccum {
    /// Total reads processed.
    pub total_reads: u64,
    /// Total assigned tags (read fragments).
    pub total_tags: u64,
    /// Tags overlapping CDS exons.
    pub cds_tags: u64,
    /// Tags overlapping 5' UTR regions.
    pub utr5_tags: u64,
    /// Tags overlapping 3' UTR regions.
    pub utr3_tags: u64,
    /// Tags overlapping intron regions.
    pub intron_tags: u64,
    /// Tags within 1 kb upstream of TSS.
    pub tss_1k_tags: u64,
    /// Tags within 5 kb upstream of TSS.
    pub tss_5k_tags: u64,
    /// Tags within 10 kb upstream of TSS.
    pub tss_10k_tags: u64,
    /// Tags within 1 kb downstream of TES.
    pub tes_1k_tags: u64,
    /// Tags within 5 kb downstream of TES.
    pub tes_5k_tags: u64,
    /// Tags within 10 kb downstream of TES.
    pub tes_10k_tags: u64,
    /// Tags not overlapping any annotated region.
    pub unassigned: u64,
}

impl ReadDistAccum {
    /// Process a single BAM record.
    pub fn process_read(&mut self, record: &bam::Record, chrom_upper: &str, regions: &RegionSets) {
        let flags = record.flags();

        // Filter: QC-fail, dup, secondary, unmapped. No MAPQ filter.
        if flags & BAM_FQCFAIL != 0
            || flags & BAM_FDUP != 0
            || flags & BAM_FSECONDARY != 0
            || flags & BAM_FUNMAP != 0
        {
            return;
        }

        self.total_reads += 1;

        // Extract exon blocks (RSeQC-compatible: M-only, S advances)
        let exon_blocks = fetch_exon_blocks_rseqc(record);

        for (block_start, block_end) in exon_blocks {
            self.total_tags += 1;
            let midpoint = block_start + (block_end - block_start) / 2;

            // Priority cascade matching RSeQC read_distribution.py:
            // CDS > UTR (5'/3', ambiguous if both) > Intron > Intergenic (TSS/TES cumulative)
            if point_in(&regions.cds_exon, chrom_upper, midpoint) {
                self.cds_tags += 1;
            } else {
                let in_utr5 = point_in(&regions.utr_5, chrom_upper, midpoint);
                let in_utr3 = point_in(&regions.utr_3, chrom_upper, midpoint);
                if in_utr5 && in_utr3 {
                    // Ambiguous UTR — RSeQC counts as unassigned
                    self.unassigned += 1;
                } else if in_utr5 {
                    self.utr5_tags += 1;
                } else if in_utr3 {
                    self.utr3_tags += 1;
                } else if point_in(&regions.intron, chrom_upper, midpoint) {
                    self.intron_tags += 1;
                } else {
                    // Intergenic — TSS/TES with cumulative counting
                    let in_tss_10k = point_in(&regions.tss_up_10kb, chrom_upper, midpoint);
                    let in_tes_10k = point_in(&regions.tes_down_10kb, chrom_upper, midpoint);
                    if in_tss_10k && in_tes_10k {
                        // Ambiguous intergenic — RSeQC counts as unassigned
                        self.unassigned += 1;
                    } else if in_tss_10k {
                        // Cumulative: 1kb ⊂ 5kb ⊂ 10kb
                        self.tss_10k_tags += 1;
                        if point_in(&regions.tss_up_5kb, chrom_upper, midpoint) {
                            self.tss_5k_tags += 1;
                            if point_in(&regions.tss_up_1kb, chrom_upper, midpoint) {
                                self.tss_1k_tags += 1;
                            }
                        }
                    } else if in_tes_10k {
                        self.tes_10k_tags += 1;
                        if point_in(&regions.tes_down_5kb, chrom_upper, midpoint) {
                            self.tes_5k_tags += 1;
                            if point_in(&regions.tes_down_1kb, chrom_upper, midpoint) {
                                self.tes_1k_tags += 1;
                            }
                        }
                    } else {
                        self.unassigned += 1;
                    }
                }
            }
        }
    }

    /// Merge another accumulator into this one.
    pub fn merge(&mut self, other: ReadDistAccum) {
        self.total_reads += other.total_reads;
        self.total_tags += other.total_tags;
        self.cds_tags += other.cds_tags;
        self.utr5_tags += other.utr5_tags;
        self.utr3_tags += other.utr3_tags;
        self.intron_tags += other.intron_tags;
        self.tss_1k_tags += other.tss_1k_tags;
        self.tss_5k_tags += other.tss_5k_tags;
        self.tss_10k_tags += other.tss_10k_tags;
        self.tes_1k_tags += other.tes_1k_tags;
        self.tes_5k_tags += other.tes_5k_tags;
        self.tes_10k_tags += other.tes_10k_tags;
        self.unassigned += other.unassigned;
    }
}

// -------------------------------------------------------------------
// junction_annotation accumulator
// -------------------------------------------------------------------

/// junction_annotation accumulator — junction classification.
#[derive(Debug, Default)]
pub struct JuncAnnotAccum {
    /// Per-junction read counts and classification.
    pub junction_counts: IndexMap<Junction, (u64, JunctionClass)>,
    /// Total splicing events observed.
    pub total_events: u64,
    /// Events matching known (annotated) junctions.
    pub known_events: u64,
    /// Events with one known and one novel splice site.
    pub partial_novel_events: u64,
    /// Events with both splice sites novel.
    pub complete_novel_events: u64,
    /// Events filtered out (below minimum intron length).
    pub filtered_events: u64,
}

impl JuncAnnotAccum {
    /// Process a single BAM record.
    pub fn process_read(
        &mut self,
        record: &bam::Record,
        chrom_upper: &str,
        ref_junctions: &ReferenceJunctions,
        min_intron: u64,
        mapq_cut: u8,
    ) {
        let flags = record.flags();

        // Filter: QC-fail, dup, secondary, unmapped, MAPQ
        if flags & BAM_FQCFAIL != 0
            || flags & BAM_FDUP != 0
            || flags & BAM_FSECONDARY != 0
            || flags & BAM_FUNMAP != 0
        {
            return;
        }
        if record.mapq() < mapq_cut {
            return;
        }

        let start_pos = record.pos() as u64;
        let cigar = record.cigar();
        let introns = common::fetch_introns(start_pos, cigar.as_ref());

        for (intron_start, intron_end) in introns {
            self.total_events += 1;

            let intron_size = intron_end.saturating_sub(intron_start);
            if intron_size < min_intron {
                self.filtered_events += 1;
                continue;
            }

            let class = classify_junction(chrom_upper, intron_start, intron_end, ref_junctions);

            match class {
                JunctionClass::Annotated => self.known_events += 1,
                JunctionClass::PartialNovel => self.partial_novel_events += 1,
                JunctionClass::CompleteNovel => self.complete_novel_events += 1,
            }

            let junction = Junction {
                chrom: chrom_upper.to_string(),
                intron_start,
                intron_end,
            };
            self.junction_counts
                .entry(junction)
                .and_modify(|(c, _)| *c += 1)
                .or_insert((1, class));
        }
    }

    /// Merge another accumulator into this one.
    pub fn merge(&mut self, other: JuncAnnotAccum) {
        for (junction, (count, class)) in other.junction_counts {
            self.junction_counts
                .entry(junction)
                .and_modify(|(c, _)| *c += count)
                .or_insert((count, class));
        }
        self.total_events += other.total_events;
        self.known_events += other.known_events;
        self.partial_novel_events += other.partial_novel_events;
        self.complete_novel_events += other.complete_novel_events;
        self.filtered_events += other.filtered_events;
    }
}

/// Classify a junction — delegates to the shared implementation in common.rs.
fn classify_junction(
    chrom: &str,
    intron_start: u64,
    intron_end: u64,
    reference: &ReferenceJunctions,
) -> JunctionClass {
    super::common::classify_junction(chrom, intron_start, intron_end, reference)
}

// -------------------------------------------------------------------
// junction_saturation accumulator
// -------------------------------------------------------------------

/// junction_saturation accumulator — collects all junction observations.
///
/// Uses hashed `u64` keys instead of heap-allocated Strings to reduce memory
/// from ~25 bytes/observation to 8 bytes/observation. The hash function is
/// the same `SipHash` used for read_duplication sequence dedup.
#[derive(Debug, Default)]
pub struct JuncSatAccum {
    /// All junction observation keys (hashed `"CHROM:start-end"`).
    pub observations: Vec<u64>,
    /// String keys corresponding to each hash, for the known_junctions lookup
    /// during into_result(). Stored as `(hash, string)` only for unique junctions.
    unique_keys: HashMap<u64, String>,
}

/// Hash a junction key string into a u64 for compact storage.
fn hash_junction_key(chrom: &str, start: u64, end: u64) -> u64 {
    use std::hash::Hasher;
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    chrom.hash(&mut hasher);
    start.hash(&mut hasher);
    end.hash(&mut hasher);
    hasher.finish()
}

impl JuncSatAccum {
    /// Process a single BAM record.
    pub fn process_read(
        &mut self,
        record: &bam::Record,
        chrom_upper: &str,
        min_intron: u64,
        mapq_cut: u8,
    ) {
        let flags = record.flags();

        // Filter: QC-fail, dup, secondary, unmapped, MAPQ
        if flags & BAM_FQCFAIL != 0
            || flags & BAM_FDUP != 0
            || flags & BAM_FSECONDARY != 0
            || flags & BAM_FUNMAP != 0
        {
            return;
        }
        if record.mapq() < mapq_cut {
            return;
        }

        let start = record.pos() as u64;
        let cigar = record.cigar();
        let introns = common::fetch_introns(start, cigar.as_ref());

        for (istart, iend) in introns {
            if iend - istart < min_intron {
                continue;
            }
            let h = hash_junction_key(chrom_upper, istart, iend);
            self.observations.push(h);
            // Track the string key only for the first occurrence of each hash
            self.unique_keys
                .entry(h)
                .or_insert_with(|| format!("{}:{}-{}", chrom_upper, istart, iend));
        }
    }

    /// Merge another accumulator into this one.
    pub fn merge(&mut self, other: JuncSatAccum) {
        self.observations.extend(other.observations);
        // Merge unique key maps (first writer wins — all map to the same string)
        for (h, key) in other.unique_keys {
            self.unique_keys.entry(h).or_insert(key);
        }
    }
}

// -------------------------------------------------------------------
// inner_distance accumulator
// -------------------------------------------------------------------

/// A single read pair's inner distance record (same as inner_distance::PairRecord).
#[derive(Debug)]
pub struct InnerDistPair {
    /// Read name (raw bytes from BAM, avoids per-read UTF-8 validation).
    pub name: Vec<u8>,
    /// Inner distance (None if different chromosomes).
    pub distance: Option<i64>,
    /// Classification label (always a static string literal).
    pub classification: &'static str,
}

/// inner_distance accumulator — paired-end inner distance sampling.
#[derive(Debug, Default)]
pub struct InnerDistAccum {
    /// Per-pair detail records.
    pub pairs: Vec<InnerDistPair>,
    /// Distances for histogram building.
    pub distances: Vec<i64>,
    /// Number of pairs processed so far.
    pub pair_num: u64,
    /// Maximum pairs to sample.
    pub sample_size: u64,
}

impl InnerDistAccum {
    /// Create a new accumulator with the given sample size.
    pub fn new(sample_size: u64) -> Self {
        InnerDistAccum {
            sample_size,
            ..Default::default()
        }
    }

    /// Process a single BAM record.
    ///
    /// Matches upstream RSeQC inner_distance.py exactly:
    /// - `read1_end = pos + qlen + splice_intron_size` (not cigar.end_pos())
    /// - Overlap uses 1-based exon positions from CIGAR M-blocks
    /// - Transcript membership checked at `read1_end - 1` and `read2_start`
    /// - mRNA distance = exonic bases between `read1_end` and `read2_start`
    pub fn process_read(
        &mut self,
        record: &bam::Record,
        chrom_upper: &str,
        exon_bitset: &ExonBitset,
        transcript_tree: &TranscriptTree,
        mapq_cut: u8,
    ) {
        // Already reached sample size
        if self.pair_num >= self.sample_size {
            return;
        }

        let flags = record.flags();

        // Filter: QC-fail, dup, secondary, unmapped, unpaired, mate unmapped, MAPQ
        if flags & BAM_FQCFAIL != 0
            || flags & BAM_FDUP != 0
            || flags & BAM_FSECONDARY != 0
            || flags & BAM_FUNMAP != 0
        {
            return;
        }
        if flags & BAM_FPAIRED == 0 || record.is_mate_unmapped() {
            return;
        }
        if record.mapq() < mapq_cut {
            return;
        }

        let read1_start = record.pos() as u64;
        let read2_start = record.mpos() as u64;

        // Different chromosomes: only process from the lower-tid side to avoid
        // double-counting in parallel mode (each chromosome is a separate worker).
        if record.tid() != record.mtid() {
            if record.tid() > record.mtid() {
                return;
            }

            self.pair_num += 1;

            self.pairs.push(InnerDistPair {
                name: record.qname().to_vec(),
                distance: None,
                classification: "sameChrom=No",
            });
            return;
        }

        // Same chromosome: mate dedup by position (both mates in same worker)
        if read2_start < read1_start {
            return;
        }
        // Same position: skip if this is read1 (upstream sets inner_distance=0 and continues)
        if read2_start == read1_start && record.is_first_in_template() {
            return;
        }

        self.pair_num += 1;

        let read_name = record.qname().to_vec();

        // Compute read1_end matching upstream RSeQC:
        //   read1_len = aligned_read.qlen  (= query_alignment_length: M+I+=/X)
        //   splice_intron_size = sum of N operations
        //   read1_end = read1_start + read1_len + splice_intron_size
        //
        // This differs from cigar.end_pos() which uses M+D+N (not I, includes D).
        // Upstream includes I but excludes D in the query length component.
        let (qalen, splice_intron_size) = compute_qalen_and_intron_size(record);
        let read1_end = read1_start + qalen + splice_intron_size;

        // Compute inner distance matching upstream exactly
        let inner_dist: i64 = if read2_start >= read1_end {
            (read2_start - read1_end) as i64
        } else {
            // Overlap: upstream enumerates 1-based exon positions from CIGAR M-blocks,
            // then counts those in range (read2_start, read1_end] (1-based).
            // This equals counting 0-based exonic positions in [read2_start, read1_end).
            let exon_blocks = fetch_exon_blocks_rseqc(record);
            let mut overlap_count: i64 = 0;
            for (ex_start, ex_end) in &exon_blocks {
                // Exon block is [ex_start, ex_end) in 0-based
                // Overlap region is [read2_start, read1_end) in 0-based
                let ov_start = (*ex_start).max(read2_start);
                let ov_end = (*ex_end).min(read1_end);
                if ov_start < ov_end {
                    overlap_count += (ov_end - ov_start) as i64;
                }
            }
            -overlap_count
        };

        // Check transcript membership using read1_end (upstream formula)
        let read1_genes =
            transcript_tree.find_overlapping(chrom_upper, read1_end.saturating_sub(1));
        let read2_genes = transcript_tree.find_overlapping(chrom_upper, read2_start);
        let common_genes: HashSet<_> = read1_genes.intersection(&read2_genes).collect();

        let classification: &'static str;

        if common_genes.is_empty() {
            classification = "sameTranscript=No,dist=genomic";
        } else if inner_dist > 0 {
            if !exon_bitset.has_chrom(chrom_upper) {
                classification = "unknownChromosome,dist=genomic";
            } else {
                let exonic_bases =
                    exon_bitset.count_exonic_bases(chrom_upper, read1_end, read2_start);

                if exonic_bases as i64 == inner_dist {
                    // sameExon: all bases between reads are exonic
                    // Upstream reports `size` which equals `inner_distance` here
                    classification = "sameTranscript=Yes,sameExon=Yes,dist=mRNA";
                } else if exonic_bases > 0 {
                    // Different exon: report mRNA distance (exonic bases only)
                    let mrna_dist = exonic_bases as i64;
                    self.pairs.push(InnerDistPair {
                        name: read_name,
                        distance: Some(mrna_dist),
                        classification: "sameTranscript=Yes,sameExon=No,dist=mRNA",
                    });
                    self.distances.push(mrna_dist);
                    return;
                } else {
                    classification = "sameTranscript=Yes,nonExonic=Yes,dist=genomic";
                }
            }
        } else {
            classification = "readPairOverlap";
        }

        self.pairs.push(InnerDistPair {
            name: read_name,
            distance: Some(inner_dist),
            classification,
        });
        self.distances.push(inner_dist);
    }

    /// Merge another accumulator into this one, respecting the sample size limit.
    ///
    /// In parallel mode each worker accumulates pairs independently.  After
    /// merging all workers the total may exceed `sample_size`, so we truncate
    /// to the limit.  This matches the upstream RSeQC behaviour of sampling
    /// at most `sample_size` pairs regardless of processing order.
    pub fn merge(&mut self, other: InnerDistAccum) {
        self.pairs.extend(other.pairs);
        self.distances.extend(other.distances);
        self.pair_num += other.pair_num;

        // Enforce the sampling limit after merging
        if self.pair_num > self.sample_size {
            let limit = self.sample_size as usize;
            self.pairs.truncate(limit);
            self.distances.truncate(limit);
            self.pair_num = self.sample_size;
        }
    }
}

// ===================================================================
// Shared CIGAR helpers
// ===================================================================

/// Compute query alignment length (M+I+=/X) and total splice intron size (N)
/// from a BAM record's CIGAR, matching upstream RSeQC's inner_distance.py.
///
/// Upstream computes `read1_end = pos + qlen + splice_intron_size` where:
/// - `qlen` = pysam `query_alignment_length` = sum of M+I+=/X CIGAR ops
/// - `splice_intron_size` = sum of N (RefSkip) CIGAR ops
///
/// This differs from `cigar.end_pos()` which uses M+D+N (includes D, excludes I).
fn compute_qalen_and_intron_size(record: &bam::Record) -> (u64, u64) {
    let mut qalen: u64 = 0;
    let mut intron_size: u64 = 0;
    for op in record.cigar().iter() {
        use rust_htslib::bam::record::Cigar::*;
        match op {
            Match(len) | Equal(len) | Diff(len) => qalen += *len as u64,
            Ins(len) => qalen += *len as u64,
            RefSkip(len) => intron_size += *len as u64,
            _ => {}
        }
    }
    (qalen, intron_size)
}

/// Extract exon blocks from CIGAR matching RSeQC's `bam_cigar.fetch_exon()`.
///
/// Only M (Match) creates blocks. D/N/S advance reference position.
/// =/X/I/H/P are ignored. The S-advances behavior is an RSeQC bug we replicate.
fn fetch_exon_blocks_rseqc(record: &bam::Record) -> Vec<(u64, u64)> {
    let mut exons = Vec::new();
    let mut chrom_st = record.pos() as u64;

    for op in record.cigar().iter() {
        use rust_htslib::bam::record::Cigar::*;
        match op {
            Match(len) => {
                let start = chrom_st;
                chrom_st += *len as u64;
                exons.push((start, chrom_st));
            }
            Del(len) | RefSkip(len) => chrom_st += *len as u64,
            SoftClip(len) => chrom_st += *len as u64, // RSeQC bug
            _ => {}
        }
    }

    exons
}

// ===================================================================
// Top-level accumulator bundle
// ===================================================================

/// Bundle of all RSeQC accumulators. Each field is `Option` — `None` when
/// that tool is disabled.
#[derive(Debug, Default)]
pub struct RseqcAccumulators {
    /// bam_stat accumulator (`None` when disabled).
    pub bam_stat: Option<BamStatAccum>,
    /// infer_experiment accumulator (`None` when disabled).
    pub infer_exp: Option<InferExpAccum>,
    /// read_duplication accumulator (`None` when disabled).
    pub read_dup: Option<ReadDupAccum>,
    /// read_distribution accumulator (`None` when disabled).
    pub read_dist: Option<ReadDistAccum>,
    /// junction_annotation accumulator (`None` when disabled).
    pub junc_annot: Option<JuncAnnotAccum>,
    /// junction_saturation accumulator (`None` when disabled).
    pub junc_sat: Option<JuncSatAccum>,
    /// inner_distance accumulator (`None` when disabled).
    pub inner_dist: Option<InnerDistAccum>,
    /// TIN accumulator (`None` when disabled).
    pub tin: Option<TinAccum>,
    /// preseq library complexity accumulator (`None` when disabled).
    pub preseq: Option<PreseqAccum>,
    /// Coverage tracks, kept per contig so workers merge safely.
    pub coverage: crate::common::coverage::bedgraph::CoverageTracks,
}

impl RseqcAccumulators {
    /// Create an empty set with no accumulators enabled.
    pub fn empty() -> Self {
        Self {
            bam_stat: None,
            infer_exp: None,
            read_dup: None,
            read_dist: None,
            junc_annot: None,
            junc_sat: None,
            inner_dist: None,
            tin: None,
            preseq: None,
            coverage: crate::common::coverage::bedgraph::CoverageTracks::default(),
        }
    }

    /// Create accumulators for all enabled tools.
    pub fn new(config: &RseqcConfig, annotations: Option<&RseqcAnnotations>) -> Self {
        RseqcAccumulators {
            bam_stat: if config.bam_stat_enabled {
                Some(BamStatAccum::default())
            } else {
                None
            },
            infer_exp: if config.infer_experiment_enabled {
                Some(InferExpAccum::default())
            } else {
                None
            },
            read_dup: if config.read_duplication_enabled {
                Some(ReadDupAccum::default())
            } else {
                None
            },
            read_dist: if config.read_distribution_enabled {
                Some(ReadDistAccum::default())
            } else {
                None
            },
            junc_annot: if config.junction_annotation_enabled {
                Some(JuncAnnotAccum::default())
            } else {
                None
            },
            junc_sat: if config.junction_saturation_enabled {
                Some(JuncSatAccum::default())
            } else {
                None
            },
            inner_dist: if config.inner_distance_enabled {
                Some(InnerDistAccum::new(config.inner_distance_sample_size))
            } else {
                None
            },
            tin: if config.tin_enabled {
                annotations.and_then(|a| a.tin_index).map(|idx| {
                    TinAccum::new(
                        idx,
                        config.mapq_cut,
                        config.tin_min_coverage,
                        config.tin_seed,
                    )
                })
            } else {
                None
            },
            preseq: if config.preseq_enabled {
                Some(PreseqAccum::new(config.preseq_max_segment_length))
            } else {
                None
            },
            coverage: crate::common::coverage::bedgraph::CoverageTracks::new(
                config.coverage_strands.clone(),
            ),
        }
    }

    /// Dispatch a BAM record to all enabled tool accumulators.
    ///
    /// This is called for EVERY record before the counting.rs filter cascade,
    /// so each tool applies its own filters internally.
    #[allow(clippy::too_many_arguments)]
    pub fn process_read(
        &mut self,
        record: &bam::Record,
        chrom: &str,
        chrom_upper: &str,
        annotations: &RseqcAnnotations,
        config: &RseqcConfig,
    ) {
        // bam_stat: sees all records, applies its own filters
        if let Some(ref mut accum) = self.bam_stat {
            accum.process_read(record, config.mapq_cut);
        }

        // read_duplication: needs chrom for position key, applies its own filters
        if let Some(ref mut accum) = self.read_dup {
            accum.process_read(record, chrom, config.mapq_cut);
        }

        // infer_experiment: needs gene model overlap
        if let (Some(ref mut accum), Some(model)) = (&mut self.infer_exp, annotations.gene_model) {
            accum.process_read(record, chrom, model, config.mapq_cut);
        }

        // Coverage tracks: no filtering at all, matching bedtools genomecov.
        self.coverage.process_read(record, chrom);

        // read_distribution: needs region sets, uses uppercased chrom
        if let (Some(ref mut accum), Some(regions)) = (&mut self.read_dist, annotations.rd_regions)
        {
            accum.process_read(record, chrom_upper, regions);
        }

        // junction_annotation: needs reference junctions, uses uppercased chrom
        if let (Some(ref mut accum), Some(ref_junctions)) =
            (&mut self.junc_annot, annotations.ref_junctions)
        {
            accum.process_read(
                record,
                chrom_upper,
                ref_junctions,
                config.min_intron,
                config.mapq_cut,
            );
        }

        // junction_saturation: uses uppercased chrom
        if let Some(ref mut accum) = &mut self.junc_sat {
            accum.process_read(record, chrom_upper, config.min_intron, config.mapq_cut);
        }

        // inner_distance: needs exon bitset + transcript tree, uses uppercased chrom
        if let (Some(ref mut accum), Some(exon_bitset), Some(transcript_tree)) = (
            &mut self.inner_dist,
            annotations.exon_bitset,
            annotations.transcript_tree,
        ) {
            accum.process_read(
                record,
                chrom_upper,
                exon_bitset,
                transcript_tree,
                config.mapq_cut,
            );
        }

        // tin: needs TinIndex, uses uppercased chrom for position lookup
        if let (Some(ref mut accum), Some(tin_index)) = (&mut self.tin, annotations.tin_index) {
            accum.process_read(record, chrom_upper, tin_index);
        }

        // preseq: counts unique fragments for library complexity estimation.
        //
        // Matches preseq v3.2.0's load_counts_BAM_pe behavior:
        //   - Only primary, mapped reads are processed.
        //   - Paired reads are merged by read name into genomic intervals.
        //   - Unpaired reads counted as individual fragments.
        //
        // Filtering (primary + mapped, no secondary/supplementary) is handled
        // inside PreseqAccum::process_read().
        if let Some(ref mut accum) = self.preseq {
            accum.process_read(record);
        }
    }

    /// Merge another set of accumulators into this one.
    pub fn merge(&mut self, other: RseqcAccumulators) {
        if let (Some(ref mut a), Some(b)) = (&mut self.bam_stat, other.bam_stat) {
            a.merge(b);
        }
        if let (Some(ref mut a), Some(b)) = (&mut self.infer_exp, other.infer_exp) {
            a.merge(b);
        }
        if let (Some(ref mut a), Some(b)) = (&mut self.read_dup, other.read_dup) {
            a.merge(b);
        }
        if let (Some(ref mut a), Some(b)) = (&mut self.read_dist, other.read_dist) {
            a.merge(b);
        }
        if let (Some(ref mut a), Some(b)) = (&mut self.junc_annot, other.junc_annot) {
            a.merge(b);
        }
        if let (Some(ref mut a), Some(b)) = (&mut self.junc_sat, other.junc_sat) {
            a.merge(b);
        }
        if let (Some(ref mut a), Some(b)) = (&mut self.inner_dist, other.inner_dist) {
            a.merge(b);
        }
        if let (Some(ref mut a), Some(b)) = (&mut self.tin, other.tin) {
            a.merge(b);
        }
        self.coverage.merge(other.coverage);
        if let (Some(ref mut a), Some(b)) = (&mut self.preseq, other.preseq) {
            a.merge(b);
        }
    }
}

// ===================================================================
// RegionSets point-query helper
// ===================================================================

/// Check if a point falls within any interval for the given chromosome
/// in a region map (HashMap<String, ChromIntervals>).
fn point_in(region_map: &HashMap<String, ChromIntervals>, chrom: &str, point: u64) -> bool {
    region_map.get(chrom).is_some_and(|ci| ci.contains(point))
}

// ===================================================================
// Converter methods: accumulator → result types for output functions
// ===================================================================

impl InferExpAccum {
    /// Convert accumulated strand counts into an `InferExperimentResult`.
    pub fn into_result(self) -> InferExperimentResult {
        let p_total: u64 = self.p_strandness.values().sum();
        let s_total: u64 = self.s_strandness.values().sum();
        let total = p_total + s_total;

        if total == 0 {
            return InferExperimentResult {
                total_sampled: 0,
                library_type: String::from("Undetermined"),
                frac_failed: 0.0,
                frac_protocol1: 0.0,
                frac_protocol2: 0.0,
            };
        }

        // PE keys for spec1: "1++", "1--", "2+-", "2-+"
        // PE keys for spec2: "1+-", "1-+", "2++", "2--"
        // SE keys for spec1: "++", "--"
        // SE keys for spec2: "+-", "-+"
        // Sum individual strand-specific keys to get protocol fractions.
        // PE spec1: 1++,1--,2+-,2-+  (fr-secondstrand / stranded)
        // PE spec2: 1+-,1-+,2++,2--  (fr-firststrand / reversely stranded)
        // SE spec1: ++,--            (sense)
        // SE spec2: +-,-+            (antisense)
        let pe_spec1 = *self.p_strandness.get("1++").unwrap_or(&0)
            + *self.p_strandness.get("1--").unwrap_or(&0)
            + *self.p_strandness.get("2+-").unwrap_or(&0)
            + *self.p_strandness.get("2-+").unwrap_or(&0);
        let pe_spec2 = *self.p_strandness.get("1+-").unwrap_or(&0)
            + *self.p_strandness.get("1-+").unwrap_or(&0)
            + *self.p_strandness.get("2++").unwrap_or(&0)
            + *self.p_strandness.get("2--").unwrap_or(&0);
        let se_spec1 =
            *self.s_strandness.get("++").unwrap_or(&0) + *self.s_strandness.get("--").unwrap_or(&0);
        let se_spec2 =
            *self.s_strandness.get("+-").unwrap_or(&0) + *self.s_strandness.get("-+").unwrap_or(&0);

        let (library_type, spec1, spec2) = if p_total > 0 && s_total > 0 {
            (
                "Mixture".to_string(),
                pe_spec1 + se_spec1,
                pe_spec2 + se_spec2,
            )
        } else if p_total > 0 {
            ("PairEnd".to_string(), pe_spec1, pe_spec2)
        } else {
            ("SingleEnd".to_string(), se_spec1, se_spec2)
        };

        let determined = spec1 + spec2;
        let failed = total - determined;
        let total_f = total as f64;

        InferExperimentResult {
            total_sampled: total,
            library_type,
            frac_failed: failed as f64 / total_f,
            frac_protocol1: spec1 as f64 / total_f,
            frac_protocol2: spec2 as f64 / total_f,
        }
    }
}

impl ReadDupAccum {
    /// Convert accumulated hash maps into a `ReadDuplicationResult`.
    ///
    /// The accumulators use u128 hash keys for sequence dedup (memory-efficient),
    /// so this just builds the duplication-level histograms from the raw counts.
    pub fn into_result(self) -> ReadDuplicationResult {
        let pos_histogram = build_dup_histogram(&self.pos_dup);
        let seq_histogram = build_dup_histogram(&self.seq_dup);
        ReadDuplicationResult {
            pos_histogram,
            seq_histogram,
        }
    }
}

/// Build duplication-level histogram from a count map.
/// Key = duplication level, Value = number of positions/sequences at that level.
fn build_dup_histogram<K: Eq + std::hash::Hash>(counts: &HashMap<K, u64>) -> BTreeMap<u64, u64> {
    let mut histogram = BTreeMap::new();
    for &count in counts.values() {
        *histogram.entry(count).or_insert(0) += 1;
    }
    histogram
}

impl ReadDistAccum {
    /// Convert accumulated tag counters into a `ReadDistributionResult`.
    pub fn into_result(self, regions: &RegionSets) -> ReadDistributionResult {
        fn sum_bases(map: &HashMap<String, ChromIntervals>) -> u64 {
            map.values().map(|ci| ci.total_bases()).sum()
        }

        let region_data = vec![
            (
                "CDS_Exons".to_string(),
                sum_bases(&regions.cds_exon),
                self.cds_tags,
            ),
            (
                "5'UTR_Exons".to_string(),
                sum_bases(&regions.utr_5),
                self.utr5_tags,
            ),
            (
                "3'UTR_Exons".to_string(),
                sum_bases(&regions.utr_3),
                self.utr3_tags,
            ),
            (
                "Introns".to_string(),
                sum_bases(&regions.intron),
                self.intron_tags,
            ),
            (
                "TSS_up_1kb".to_string(),
                sum_bases(&regions.tss_up_1kb),
                self.tss_1k_tags,
            ),
            (
                "TSS_up_5kb".to_string(),
                sum_bases(&regions.tss_up_5kb),
                self.tss_5k_tags,
            ),
            (
                "TSS_up_10kb".to_string(),
                sum_bases(&regions.tss_up_10kb),
                self.tss_10k_tags,
            ),
            (
                "TES_down_1kb".to_string(),
                sum_bases(&regions.tes_down_1kb),
                self.tes_1k_tags,
            ),
            (
                "TES_down_5kb".to_string(),
                sum_bases(&regions.tes_down_5kb),
                self.tes_5k_tags,
            ),
            (
                "TES_down_10kb".to_string(),
                sum_bases(&regions.tes_down_10kb),
                self.tes_10k_tags,
            ),
        ];

        ReadDistributionResult {
            total_reads: self.total_reads,
            total_tags: self.total_tags,
            regions: region_data,
            unassigned_tags: self.unassigned,
        }
    }
}

impl JuncAnnotAccum {
    /// Convert accumulated junction data into `JunctionResults`.
    ///
    /// Junctions are reordered to match the upstream RSeQC output: grouped
    /// by BAM header chromosome order, with within-chromosome insertion
    /// order preserved.  Because each worker thread processes entire
    /// chromosomes sequentially, the per-chromosome insertion order already
    /// matches what the upstream single-threaded Python produces.  We just
    /// need to interleave chromosome groups in BAM header order.
    pub fn into_result(self, bam_header_refs: &[(String, u64)]) -> JunctionResults {
        // Build chrom → index map (uppercased to match Junction.chrom).
        let chrom_order: std::collections::HashMap<String, usize> = bam_header_refs
            .iter()
            .enumerate()
            .map(|(i, (name, _))| (name.to_uppercase(), i))
            .collect();
        let sentinel = bam_header_refs.len();

        // Group junctions by chromosome, preserving within-group insertion order.
        let mut per_chrom: std::collections::BTreeMap<usize, Vec<_>> =
            std::collections::BTreeMap::new();
        for (junction, (count, class)) in self.junction_counts {
            let idx = chrom_order
                .get(&junction.chrom)
                .copied()
                .unwrap_or(sentinel);
            per_chrom
                .entry(idx)
                .or_default()
                .push((junction, (count, class)));
        }

        // Flatten groups in BAM header chromosome order.
        let junctions: IndexMap<_, _> = per_chrom
            .into_values()
            .flat_map(|v| v.into_iter())
            .collect();

        JunctionResults {
            junctions,
            total_events: self.total_events,
            known_events: self.known_events,
            partial_novel_events: self.partial_novel_events,
            complete_novel_events: self.complete_novel_events,
            filtered_events: self.filtered_events,
        }
    }
}

impl JuncSatAccum {
    /// Post-process accumulated observations into a `SaturationResult`.
    ///
    /// Performs the shuffle (Phase 2) and incremental subsampling (Phase 3)
    /// that were previously done in the standalone `junction_saturation()` function.
    ///
    /// Uses incremental known/novel counting: instead of iterating all unique
    /// junctions at each percentage step (O(U*P)), maintains running counters
    /// updated only when new unique junctions first appear or cross the
    /// min_coverage threshold.
    pub fn into_result(
        mut self,
        known_junctions: &KnownJunctionSet,
        sample_start: u32,
        sample_end: u32,
        sample_step: u32,
        min_coverage: u32,
        seed: u64,
    ) -> SaturationResult {
        use rand::seq::SliceRandom;
        use rand::SeedableRng;
        use rand_chacha::ChaCha8Rng;

        // Pre-build a HashSet<u64> from known_junctions for O(1) hash-based lookup
        let known_hashes: HashSet<u64> = self
            .unique_keys
            .iter()
            .filter(|(_, key)| known_junctions.junctions.contains(*key))
            .map(|(&h, _)| h)
            .collect();

        // Phase 2: deterministic shuffle
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        self.observations.shuffle(&mut rng);

        // Build percentage series
        let mut percentages: Vec<u32> = (sample_start..=sample_end)
            .step_by(sample_step as usize)
            .collect();
        if *percentages.last().unwrap_or(&0) != 100 {
            percentages.push(100);
        }

        // Phase 3: incremental sampling with running counters
        let total = self.observations.len();
        let mut junction_counts: HashMap<u64, u32> = HashMap::new();
        let mut prev_end = 0;
        let mut known_counts = Vec::with_capacity(percentages.len());
        let mut novel_counts = Vec::with_capacity(percentages.len());
        let mut all_counts = Vec::with_capacity(percentages.len());

        // Running counters — updated incrementally as new observations arrive
        let mut running_known: usize = 0;
        let mut running_novel: usize = 0;

        for &pct in &percentages {
            let index_end = total * pct as usize / 100;
            for &obs_hash in &self.observations[prev_end..index_end] {
                let count = junction_counts.entry(obs_hash).or_insert(0);
                *count += 1;

                let is_known = known_hashes.contains(&obs_hash);
                if *count == 1 {
                    // First time seeing this junction
                    if is_known {
                        if min_coverage <= 1 {
                            running_known += 1;
                        }
                        // else: known but below threshold, don't count yet
                    } else {
                        running_novel += 1;
                    }
                } else if is_known && *count == min_coverage && min_coverage > 1 {
                    // Junction just crossed the min_coverage threshold
                    running_known += 1;
                }
            }
            prev_end = index_end;

            known_counts.push(running_known);
            novel_counts.push(running_novel);
            all_counts.push(junction_counts.len());
        }

        SaturationResult {
            percentages,
            known_counts,
            novel_counts,
            all_counts,
        }
    }
}

impl InnerDistAccum {
    /// Convert accumulated pair data into an `InnerDistanceResult`.
    pub fn into_result(
        self,
        lower_bound: i64,
        upper_bound: i64,
        step: i64,
    ) -> Result<InnerDistanceResult> {
        let pairs: Vec<PairRecord> = self
            .pairs
            .into_iter()
            .map(|p| PairRecord {
                name: String::from_utf8_lossy(&p.name).into_owned(),
                distance: p.distance,
                classification: p.classification.to_owned(),
            })
            .collect();
        let histogram = build_histogram(&self.distances, lower_bound, upper_bound, step)?;
        Ok(InnerDistanceResult {
            pairs,
            histogram,
            total_pairs: self.pair_num,
        })
    }
}
