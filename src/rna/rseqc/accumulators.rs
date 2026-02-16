//! Per-read accumulator structs for single-pass RSeQC integration.
//!
//! Each RSeQC tool has an accumulator that collects per-read data during the
//! main BAM counting loop. Accumulators are created per chromosome worker and
//! merged after parallel processing, just like `ChromResult`.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::hash::{Hash, Hasher};

use rust_htslib::bam;

use super::bam_stat::BamStatResult;
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
use crate::rna::preseq::PreseqAccum;

// ===================================================================
// BAM flag constants
// ===================================================================

const BAM_FPAIRED: u16 = 0x1;
const BAM_FPROPER_PAIR: u16 = 0x2;
const BAM_FUNMAP: u16 = 0x4;
const BAM_FREVERSE: u16 = 0x10;
const BAM_FREAD1: u16 = 0x40;
const BAM_FREAD2: u16 = 0x80;
const BAM_FSECONDARY: u16 = 0x100;
const BAM_FQCFAIL: u16 = 0x200;
const BAM_FDUP: u16 = 0x400;
const BAM_FSUPPLEMENTARY: u16 = 0x800;

// ===================================================================
// Shared references to annotation data
// ===================================================================

/// Read-only annotation data shared across all chromosome workers.
///
/// Each field is `Option` — `None` when the corresponding tool is disabled
/// or the annotation source doesn't support it (e.g., BED mode has no GTF).
pub struct RseqcAnnotations<'a> {
    /// Gene model for infer_experiment.
    pub gene_model: Option<&'a GeneModel>,
    /// Reference junctions for junction_annotation.
    pub ref_junctions: Option<&'a ReferenceJunctions>,
    /// Known junction set for junction_saturation (used at result-conversion time).
    #[allow(dead_code)]
    pub known_junctions: Option<&'a KnownJunctionSet>,
    /// Genomic region sets for read_distribution.
    pub rd_regions: Option<&'a RegionSets>,
    /// Exon bitset for inner_distance.
    pub exon_bitset: Option<&'a ExonBitset>,
    /// Transcript tree for inner_distance.
    pub transcript_tree: Option<&'a TranscriptTree>,
    /// Chromosomes present in the known junction set (precomputed for fast lookup).
    pub ref_chroms: Option<&'a HashSet<String>>,
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
    /// Which tools are enabled.
    pub bam_stat_enabled: bool,
    pub infer_experiment_enabled: bool,
    pub read_duplication_enabled: bool,
    pub read_distribution_enabled: bool,
    pub junction_annotation_enabled: bool,
    pub junction_saturation_enabled: bool,
    pub inner_distance_enabled: bool,
    /// Whether TIN analysis is enabled.
    pub tin_enabled: bool,
    /// Number of equally-spaced sampling positions per transcript for TIN.
    pub tin_sample_size: usize,
    /// Minimum number of read starts for a transcript to compute TIN.
    pub tin_min_coverage: u32,
    /// Whether preseq library complexity estimation is enabled.
    pub preseq_enabled: bool,
}

// ===================================================================
// Per-tool accumulators
// ===================================================================

/// bam_stat accumulator — simple flag/MAPQ counting.
///
/// Also collects the additional counters needed for samtools-compatible
/// flagstat, idxstats, and stats output.
#[derive(Debug, Default)]
pub struct BamStatAccum {
    // --- RSeQC bam_stat fields (original) ---
    /// Total BAM records seen (primary + secondary + supplementary + unmapped).
    pub total_records: u64,
    /// Records with QC-fail flag (0x200).
    pub qc_failed: u64,
    /// Records with duplicate flag (0x400).
    pub duplicates: u64,
    /// Secondary alignment records (0x100). RSeQC calls these "non-primary".
    pub non_primary: u64,
    /// Unmapped reads (0x4).
    pub unmapped: u64,
    /// Mapped reads with MAPQ < cutoff.
    pub non_unique: u64,
    /// Mapped reads with MAPQ >= cutoff (uniquely mapped).
    pub unique: u64,
    /// Among unique reads: read1 in a pair.
    pub read_1: u64,
    /// Among unique reads: read2 in a pair.
    pub read_2: u64,
    /// Among unique reads: forward strand.
    pub forward: u64,
    /// Among unique reads: reverse strand.
    pub reverse: u64,
    /// Among unique reads: has splice junction (CIGAR N).
    pub splice: u64,
    /// Among unique reads: no splice junctions.
    pub non_splice: u64,
    /// Among unique reads: in proper pairs (0x2).
    pub proper_pairs: u64,
    /// Among proper-paired unique reads: mates on different chromosomes.
    pub proper_pair_diff_chrom: u64,
    /// MAPQ distribution for primary, non-QC-fail, non-dup, mapped reads.
    pub mapq_distribution: BTreeMap<u8, u64>,

    // --- samtools flagstat additional fields ---
    /// Secondary alignments (0x100) — counted independently of QC/dup.
    pub secondary: u64,
    /// Supplementary alignments (0x800) — counted independently of QC/dup.
    pub supplementary: u64,
    /// All mapped records (not 0x4), regardless of QC/dup.
    pub mapped: u64,
    /// Paired reads (0x1), regardless of QC/dup.
    pub paired_flagstat: u64,
    /// Read1 in pair (0x40), regardless of QC/dup — for flagstat.
    pub read1_flagstat: u64,
    /// Read2 in pair (0x80), regardless of QC/dup — for flagstat.
    pub read2_flagstat: u64,
    /// First fragments for samtools stats: primary reads that are not "last fragments".
    pub first_fragments: u64,
    /// Last fragments for samtools stats: primary reads with 0x80 flag.
    pub last_fragments: u64,
    /// Properly paired reads (0x1 + 0x2), regardless of QC/dup.
    pub properly_paired: u64,
    /// Both mates mapped (paired + both !unmapped).
    pub both_mapped: u64,
    /// Singletons (paired, this mapped, mate unmapped).
    pub singletons: u64,
    /// Paired, both mapped, different reference.
    pub mate_diff_chr: u64,
    /// Paired, both mapped, different reference, MAPQ >= 5.
    pub mate_diff_chr_mapq5: u64,

    // --- samtools idxstats additional fields ---
    /// Per-reference (tid) mapped and unmapped counts.
    pub chrom_counts: HashMap<i32, (u64, u64)>,
    /// Unmapped reads with no reference (tid < 0).
    pub unplaced_unmapped: u64,

    // --- samtools stats SN additional fields ---
    /// Sum of query sequence lengths for all primary reads (non-secondary, non-supplementary).
    pub total_len: u64,
    /// Sum of first fragment (read1 or unpaired) sequence lengths.
    pub total_first_fragment_len: u64,
    /// Sum of last fragment (read2) sequence lengths.
    pub total_last_fragment_len: u64,
    /// Sum of query lengths for mapped primary reads.
    pub bases_mapped: u64,
    /// Sum of M/=/X CIGAR operations for mapped primary reads.
    pub bases_mapped_cigar: u64,
    /// Sum of query lengths for duplicate-flagged primary reads.
    pub bases_duplicated: u64,
    /// Maximum query sequence length (among primary reads).
    pub max_len: u64,
    /// Maximum first-fragment sequence length.
    pub max_first_fragment_len: u64,
    /// Maximum last-fragment sequence length.
    pub max_last_fragment_len: u64,
    /// Sum of average per-read base qualities (for average-of-averages).
    pub quality_sum: f64,
    /// Number of reads contributing to quality_sum (primary, non-QC-fail).
    pub quality_count: u64,
    /// Sum of NM tag values across mapped primary reads.
    pub mismatches: u64,
    /// Sum of absolute TLEN for properly paired primary reads (for mean insert size).
    pub insert_size_sum: f64,
    /// Sum of squared TLEN for insert size standard deviation.
    pub insert_size_sq_sum: f64,
    /// Number of reads contributing to insert size stats.
    pub insert_size_count: u64,
    /// Inward-oriented pairs (FR).
    pub inward_pairs: u64,
    /// Outward-oriented pairs (RF).
    pub outward_pairs: u64,
    /// Other orientation pairs (FF, RR).
    pub other_orientation: u64,
    /// Primary reads that are paired (for "raw total sequences" = paired means /2 fragments).
    pub primary_paired: u64,
    /// Total primary reads (non-secondary, non-supplementary).
    pub primary_count: u64,
    /// Primary mapped reads count (non-secondary, non-supplementary, !unmapped).
    pub primary_mapped: u64,
    /// Primary duplicate reads.
    pub primary_duplicates: u64,
    /// All mapped reads with MAPQ = 0 (including secondary/supplementary).
    pub reads_mq0: u64,
    /// Primary non-QC-fail mapped paired reads where mate is also mapped.
    pub reads_mapped_and_paired: u64,
}

/// Mate unmapped flag (0x8).
const BAM_FMUNMAP: u16 = 0x8;
/// Mate reverse strand flag (0x20).
const BAM_FMREVERSE: u16 = 0x20;

impl BamStatAccum {
    /// Process a single BAM record. Called for EVERY record (before counting filters).
    ///
    /// Collects counters for:
    /// - RSeQC bam_stat (original cascade with early returns)
    /// - samtools flagstat (counts all records independently)
    /// - samtools idxstats (per-reference mapped/unmapped counts)
    /// - samtools stats SN section (sequence lengths, quality, insert size, etc.)
    pub fn process_read(&mut self, record: &bam::Record, mapq_cut: u8) {
        let flags = record.flags();
        self.total_records += 1;

        let is_secondary = flags & BAM_FSECONDARY != 0;
        let is_supplementary = flags & BAM_FSUPPLEMENTARY != 0;
        let is_unmapped = flags & BAM_FUNMAP != 0;
        let is_paired = flags & BAM_FPAIRED != 0;
        let is_dup = flags & BAM_FDUP != 0;
        let is_qcfail = flags & BAM_FQCFAIL != 0;
        let is_primary = !is_secondary && !is_supplementary;
        let is_mapped = !is_unmapped;
        let tid = record.tid();
        let mapq = record.mapq();

        // =================================================================
        // samtools flagstat counters (count ALL records, no early returns)
        // =================================================================
        if is_secondary {
            self.secondary += 1;
        }
        if is_supplementary {
            self.supplementary += 1;
        }
        if is_mapped {
            self.mapped += 1;
            // samtools stats: reads MQ0 counts all mapped reads (including
            // secondary/supplementary) with MAPQ=0
            if record.mapq() == 0 {
                self.reads_mq0 += 1;
            }
        }
        // samtools stats: "1st fragments" / "last fragments" count primary reads only
        // For paired reads: read2 flag -> last, everything else -> 1st
        // For SE reads (no PAIRED flag): all counted as 1st fragments
        if is_primary {
            if flags & BAM_FREAD2 != 0 {
                self.last_fragments += 1;
            } else {
                self.first_fragments += 1;
            }
        }
        // samtools flagstat: paired-read metrics count PRIMARY reads only
        // (secondary/supplementary are excluded from paired/read1/read2/properly-paired counts)
        if is_paired && is_primary {
            self.paired_flagstat += 1;
            if flags & BAM_FREAD1 != 0 {
                self.read1_flagstat += 1;
            }
            if flags & BAM_FREAD2 != 0 {
                self.read2_flagstat += 1;
            }
            if flags & BAM_FPROPER_PAIR != 0 {
                self.properly_paired += 1;
            }
            let mate_unmapped = flags & BAM_FMUNMAP != 0;
            if is_mapped && !mate_unmapped {
                self.both_mapped += 1;
                if tid != record.mtid() {
                    self.mate_diff_chr += 1;
                    if mapq >= 5 {
                        self.mate_diff_chr_mapq5 += 1;
                    }
                }
            }
            if is_mapped && mate_unmapped {
                self.singletons += 1;
            }
        }

        // =================================================================
        // samtools idxstats counters (per-reference)
        // =================================================================
        if is_unmapped {
            if tid >= 0 {
                // Unmapped read placed on a reference (has tid)
                self.chrom_counts.entry(tid).or_insert((0, 0)).1 += 1;
            } else {
                self.unplaced_unmapped += 1;
            }
        } else if tid >= 0 {
            // Mapped read
            self.chrom_counts.entry(tid).or_insert((0, 0)).0 += 1;
        }

        // =================================================================
        // samtools stats SN counters (primary reads only)
        // =================================================================
        if is_primary {
            self.primary_count += 1;
            let seq_len = record.seq_len() as u64;
            let mate_unmapped = flags & BAM_FMUNMAP != 0;

            self.total_len += seq_len;
            if is_paired && flags & BAM_FREAD2 != 0 {
                self.total_last_fragment_len += seq_len;
                if seq_len > self.max_last_fragment_len {
                    self.max_last_fragment_len = seq_len;
                }
            } else {
                self.total_first_fragment_len += seq_len;
                if seq_len > self.max_first_fragment_len {
                    self.max_first_fragment_len = seq_len;
                }
            }
            if seq_len > self.max_len {
                self.max_len = seq_len;
            }

            if is_paired {
                self.primary_paired += 1;
            }
            if is_dup {
                self.primary_duplicates += 1;
                self.bases_duplicated += seq_len;
            }
            // "reads mapped and paired" for samtools stats: primary, non-QC-fail,
            // mapped, paired, mate also mapped
            if is_mapped && is_paired && !is_qcfail && !mate_unmapped {
                self.reads_mapped_and_paired += 1;
            }
            if is_mapped {
                self.primary_mapped += 1;
                self.bases_mapped += seq_len;

                // CIGAR-based mapped bases (M/=/X operations)
                use rust_htslib::bam::record::Cigar as C;
                let cigar_mapped: u64 = record
                    .cigar()
                    .iter()
                    .map(|op| match op {
                        C::Match(n) | C::Equal(n) | C::Diff(n) => u64::from(*n),
                        _ => 0,
                    })
                    .sum();
                self.bases_mapped_cigar += cigar_mapped;

                // NM tag (edit distance)
                if let Ok(rust_htslib::bam::record::Aux::U8(nm)) = record.aux(b"NM") {
                    self.mismatches += u64::from(nm);
                } else if let Ok(rust_htslib::bam::record::Aux::U16(nm)) = record.aux(b"NM") {
                    self.mismatches += u64::from(nm);
                } else if let Ok(rust_htslib::bam::record::Aux::U32(nm)) = record.aux(b"NM") {
                    self.mismatches += u64::from(nm);
                } else if let Ok(rust_htslib::bam::record::Aux::I8(nm)) = record.aux(b"NM") {
                    if nm > 0 {
                        self.mismatches += nm as u64;
                    }
                } else if let Ok(rust_htslib::bam::record::Aux::I16(nm)) = record.aux(b"NM") {
                    if nm > 0 {
                        self.mismatches += nm as u64;
                    }
                } else if let Ok(rust_htslib::bam::record::Aux::I32(nm)) = record.aux(b"NM") {
                    if nm > 0 {
                        self.mismatches += nm as u64;
                    }
                }

                // Insert size for properly paired primary reads on the same
                // chromosome. samtools stats counts insert size only for the read
                // with positive TLEN to avoid double-counting, and excludes
                // reads with TLEN > 65536 (outliers / structural variants).
                if is_paired && flags & BAM_FPROPER_PAIR != 0 {
                    let tid = record.tid();
                    let mtid = record.mtid();
                    if tid == mtid {
                        let tlen = record.insert_size();
                        if tlen > 0 && tlen <= 65536 {
                            let abs_tlen = tlen as f64;
                            self.insert_size_sum += abs_tlen;
                            self.insert_size_sq_sum += abs_tlen * abs_tlen;
                            self.insert_size_count += 1;
                        }
                    }
                }

                // Pair orientation for mapped paired reads where mate is also mapped.
                // samtools stats determines orientation based on leftmost position
                // and only counts from the read with the smaller position (to avoid
                // double-counting).
                if is_paired && !mate_unmapped {
                    let pos = record.pos();
                    let mpos = record.mpos();
                    let tid = record.tid();
                    let mtid = record.mtid();

                    // Only count when this read is the leftmost (or when same
                    // position, use read1 as tiebreaker)
                    let is_leftmost = if tid != mtid {
                        false // different chromosomes - skip orientation
                    } else if pos != mpos {
                        pos < mpos
                    } else {
                        flags & BAM_FREAD1 != 0 // tiebreaker: read1
                    };

                    if is_leftmost {
                        let this_rev = flags & BAM_FREVERSE != 0;
                        let mate_rev = flags & BAM_FMREVERSE != 0;

                        // Orientation based on strand of left and right read
                        if !this_rev && mate_rev {
                            self.inward_pairs += 1; // F...R -> inward
                        } else if this_rev && !mate_rev {
                            self.outward_pairs += 1; // R...F -> outward
                        } else {
                            self.other_orientation += 1; // F...F or R...R
                        }
                    }
                }
            }

            // Average quality for primary non-QC-fail reads
            if !is_qcfail {
                let quals = record.qual();
                if !quals.is_empty() {
                    let avg_q: f64 =
                        quals.iter().map(|&q| f64::from(q)).sum::<f64>() / quals.len() as f64;
                    self.quality_sum += avg_q;
                    self.quality_count += 1;
                }
            }
        }

        // =================================================================
        // RSeQC bam_stat cascade (original logic, with early returns)
        // =================================================================

        // 1. QC-failed
        if is_qcfail {
            self.qc_failed += 1;
            return;
        }

        // 2. Duplicate
        if is_dup {
            self.duplicates += 1;
            return;
        }

        // 3. Secondary (non-primary) — NOT supplementary
        if is_secondary {
            self.non_primary += 1;
            return;
        }

        // 4. Unmapped
        if is_unmapped {
            self.unmapped += 1;
            return;
        }

        // Mapped primary non-QC-fail non-dup: record MAPQ distribution
        *self.mapq_distribution.entry(mapq).or_insert(0) += 1;

        // 5. MAPQ classification
        if mapq < mapq_cut {
            self.non_unique += 1;
            return;
        }

        // Uniquely mapped
        self.unique += 1;

        if flags & BAM_FREAD1 != 0 {
            self.read_1 += 1;
        }
        if flags & BAM_FREAD2 != 0 {
            self.read_2 += 1;
        }
        if flags & BAM_FREVERSE != 0 {
            self.reverse += 1;
        } else {
            self.forward += 1;
        }

        // Splice detection: CIGAR N operation
        let has_splice = record
            .cigar()
            .iter()
            .any(|op| matches!(op, rust_htslib::bam::record::Cigar::RefSkip(_)));
        if has_splice {
            self.splice += 1;
        } else {
            self.non_splice += 1;
        }

        // Proper pair analysis
        if is_paired && flags & BAM_FPROPER_PAIR != 0 {
            self.proper_pairs += 1;
            if tid != record.mtid() {
                self.proper_pair_diff_chrom += 1;
            }
        }
    }

    /// Merge another accumulator into this one.
    pub fn merge(&mut self, other: BamStatAccum) {
        // RSeQC bam_stat fields
        self.total_records += other.total_records;
        self.qc_failed += other.qc_failed;
        self.duplicates += other.duplicates;
        self.non_primary += other.non_primary;
        self.unmapped += other.unmapped;
        self.non_unique += other.non_unique;
        self.unique += other.unique;
        self.read_1 += other.read_1;
        self.read_2 += other.read_2;
        self.forward += other.forward;
        self.reverse += other.reverse;
        self.splice += other.splice;
        self.non_splice += other.non_splice;
        self.proper_pairs += other.proper_pairs;
        self.proper_pair_diff_chrom += other.proper_pair_diff_chrom;
        for (mapq, count) in other.mapq_distribution {
            *self.mapq_distribution.entry(mapq).or_insert(0) += count;
        }

        // samtools flagstat fields
        self.secondary += other.secondary;
        self.supplementary += other.supplementary;
        self.mapped += other.mapped;
        self.paired_flagstat += other.paired_flagstat;
        self.read1_flagstat += other.read1_flagstat;
        self.read2_flagstat += other.read2_flagstat;
        self.first_fragments += other.first_fragments;
        self.last_fragments += other.last_fragments;
        self.properly_paired += other.properly_paired;
        self.both_mapped += other.both_mapped;
        self.singletons += other.singletons;
        self.mate_diff_chr += other.mate_diff_chr;
        self.mate_diff_chr_mapq5 += other.mate_diff_chr_mapq5;

        // samtools idxstats fields
        for (tid, (m, u)) in other.chrom_counts {
            let entry = self.chrom_counts.entry(tid).or_insert((0, 0));
            entry.0 += m;
            entry.1 += u;
        }
        self.unplaced_unmapped += other.unplaced_unmapped;

        // samtools stats SN fields
        self.total_len += other.total_len;
        self.total_first_fragment_len += other.total_first_fragment_len;
        self.total_last_fragment_len += other.total_last_fragment_len;
        self.bases_mapped += other.bases_mapped;
        self.bases_mapped_cigar += other.bases_mapped_cigar;
        self.bases_duplicated += other.bases_duplicated;
        if other.max_len > self.max_len {
            self.max_len = other.max_len;
        }
        if other.max_first_fragment_len > self.max_first_fragment_len {
            self.max_first_fragment_len = other.max_first_fragment_len;
        }
        if other.max_last_fragment_len > self.max_last_fragment_len {
            self.max_last_fragment_len = other.max_last_fragment_len;
        }
        self.quality_sum += other.quality_sum;
        self.quality_count += other.quality_count;
        self.mismatches += other.mismatches;
        self.insert_size_sum += other.insert_size_sum;
        self.insert_size_sq_sum += other.insert_size_sq_sum;
        self.insert_size_count += other.insert_size_count;
        self.inward_pairs += other.inward_pairs;
        self.outward_pairs += other.outward_pairs;
        self.other_orientation += other.other_orientation;
        self.primary_paired += other.primary_paired;
        self.primary_count += other.primary_count;
        self.primary_mapped += other.primary_mapped;
        self.primary_duplicates += other.primary_duplicates;
        self.reads_mq0 += other.reads_mq0;
        self.reads_mapped_and_paired += other.reads_mapped_and_paired;
    }
}

// -------------------------------------------------------------------
// infer_experiment accumulator
// -------------------------------------------------------------------

/// infer_experiment accumulator — strand protocol inference via sampling.
#[derive(Debug, Default)]
pub struct InferExpAccum {
    /// Paired-end strand class counts (keys: "1++", "1--", "2+-", "2-+", etc.)
    pub p_strandness: HashMap<String, u64>,
    /// Single-end strand class counts (keys: "++", "--", "+-", "-+", etc.)
    pub s_strandness: HashMap<String, u64>,
    /// Number of usable reads sampled so far.
    pub count: u64,
    /// Maximum reads to sample.
    pub sample_size: u64,
}

impl InferExpAccum {
    /// Create a new accumulator with the given sample size.
    pub fn new(sample_size: u64) -> Self {
        InferExpAccum {
            sample_size,
            ..Default::default()
        }
    }

    /// Process a single BAM record.
    pub fn process_read(
        &mut self,
        record: &bam::Record,
        chrom: &str,
        model: &GeneModel,
        mapq_cut: u8,
    ) {
        // Already reached sample size
        if self.count >= self.sample_size {
            return;
        }

        let flags = record.flags();

        // Skip QC-fail, dup, secondary, unmapped, supplementary
        if flags & BAM_FQCFAIL != 0
            || flags & BAM_FDUP != 0
            || flags & BAM_FSECONDARY != 0
            || flags & BAM_FUNMAP != 0
            || flags & BAM_FSUPPLEMENTARY != 0
        {
            return;
        }

        if record.mapq() < mapq_cut {
            return;
        }

        let map_strand = if record.is_reverse() { '-' } else { '+' };

        // Compute query alignment length (M+I+=+X) to match RSeQC's qlen
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

        let strand_str: String = strands
            .iter()
            .map(|&s| (s as char).to_string())
            .collect::<Vec<String>>()
            .join(":");

        if record.is_paired() {
            let read_id = if record.is_first_in_template() {
                "1"
            } else {
                "2"
            };
            let key = format!("{}{}{}", read_id, map_strand, strand_str);
            *self.p_strandness.entry(key).or_insert(0) += 1;
        } else {
            let key = format!("{}{}", map_strand, strand_str);
            *self.s_strandness.entry(key).or_insert(0) += 1;
        }

        self.count += 1;
    }

    /// Merge another accumulator into this one.
    pub fn merge(&mut self, other: InferExpAccum) {
        for (key, count) in other.p_strandness {
            *self.p_strandness.entry(key).or_insert(0) += count;
        }
        for (key, count) in other.s_strandness {
            *self.s_strandness.entry(key).or_insert(0) += count;
        }
        self.count += other.count;
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

    // FNV-1a hash of the position key components
    let mut h: u64 = 0xcbf29ce484222325;
    let mix = |h: &mut u64, bytes: &[u8]| {
        for &b in bytes {
            *h ^= b as u64;
            *h = h.wrapping_mul(0x100000001b3);
        }
    };

    mix(&mut h, chrom.as_bytes());
    mix(&mut h, &pos.to_le_bytes());

    let mut ref_pos = pos;
    for op in cigar.iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let end = ref_pos + *len as i64;
                mix(&mut h, &ref_pos.to_le_bytes());
                mix(&mut h, &end.to_le_bytes());
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
    pub total_reads: u64,
    pub total_tags: u64,
    pub cds_tags: u64,
    pub utr5_tags: u64,
    pub utr3_tags: u64,
    pub intron_tags: u64,
    pub tss_1k_tags: u64,
    pub tss_5k_tags: u64,
    pub tss_10k_tags: u64,
    pub tes_1k_tags: u64,
    pub tes_5k_tags: u64,
    pub tes_10k_tags: u64,
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
    pub junction_counts: HashMap<Junction, (u64, JunctionClass)>,
    pub total_events: u64,
    pub known_events: u64,
    pub partial_novel_events: u64,
    pub complete_novel_events: u64,
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

/// Classify a junction (duplicated from junction_annotation to avoid circular deps).
fn classify_junction(
    chrom: &str,
    intron_start: u64,
    intron_end: u64,
    reference: &ReferenceJunctions,
) -> JunctionClass {
    let start_known = reference
        .intron_starts
        .get(chrom)
        .is_some_and(|s| s.contains(&intron_start));
    let end_known = reference
        .intron_ends
        .get(chrom)
        .is_some_and(|s| s.contains(&intron_end));

    match (start_known, end_known) {
        (true, true) => JunctionClass::Annotated,
        (false, false) => JunctionClass::CompleteNovel,
        _ => JunctionClass::PartialNovel,
    }
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
        ref_chroms: &HashSet<String>,
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

        // Skip chromosomes not in reference
        if !ref_chroms.contains(chrom_upper) {
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
    /// Read name.
    pub name: String,
    /// Inner distance (None if different chromosomes).
    pub distance: Option<i64>,
    /// Classification string.
    pub classification: String,
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

            let read_name = String::from_utf8_lossy(record.qname()).to_string();
            self.pairs.push(InnerDistPair {
                name: read_name,
                distance: None,
                classification: "sameChrom=No".to_string(),
            });
            return;
        }

        // Same chromosome: mate dedup by position (both mates in same worker)
        if read2_start < read1_start {
            return;
        }
        // Same position: skip if this is read1
        if read2_start == read1_start && record.is_first_in_template() {
            return;
        }

        self.pair_num += 1;

        let read_name = String::from_utf8_lossy(record.qname()).to_string();

        // Compute read1_end from CIGAR
        let read1_end = record.cigar().end_pos() as u64;

        // Compute inner distance
        let inner_dist: i64 = if read2_start >= read1_end {
            (read2_start - read1_end) as i64
        } else {
            // Overlap: count exonic positions of read1 in overlap region
            let exon_blocks = fetch_exon_blocks_rseqc(record);
            let mut overlap_count: i64 = 0;
            for (ex_start, ex_end) in &exon_blocks {
                let ov_start = (*ex_start).max(read2_start);
                let ov_end = (*ex_end).min(read1_end);
                if ov_start < ov_end {
                    overlap_count += (ov_end - ov_start) as i64;
                }
            }
            -overlap_count
        };

        // Check transcript membership
        let read1_genes =
            transcript_tree.find_overlapping(chrom_upper, read1_end.saturating_sub(1));
        let read2_genes = transcript_tree.find_overlapping(chrom_upper, read2_start);
        let common_genes: HashSet<_> = read1_genes.intersection(&read2_genes).collect();

        let classification: String;

        if common_genes.is_empty() {
            classification = "sameTranscript=No,dist=genomic".to_string();
        } else if inner_dist > 0 {
            if !exon_bitset.has_chrom(chrom_upper) {
                classification = "unknownChromosome,dist=genomic".to_string();
            } else {
                let exonic_bases =
                    exon_bitset.count_exonic_bases(chrom_upper, read1_end, read2_start);

                if exonic_bases as i64 == inner_dist {
                    classification = "sameTranscript=Yes,sameExon=Yes,dist=mRNA".to_string();
                } else if exonic_bases > 0 {
                    classification = "sameTranscript=Yes,sameExon=No,dist=mRNA".to_string();
                    let mrna_dist = exonic_bases as i64;
                    self.pairs.push(InnerDistPair {
                        name: read_name,
                        distance: Some(mrna_dist),
                        classification,
                    });
                    self.distances.push(mrna_dist);
                    return;
                } else {
                    classification = "sameTranscript=Yes,nonExonic=Yes,dist=genomic".to_string();
                }
            }
        } else {
            classification = "readPairOverlap".to_string();
        }

        self.pairs.push(InnerDistPair {
            name: read_name,
            distance: Some(inner_dist),
            classification,
        });
        self.distances.push(inner_dist);
    }

    /// Merge another accumulator into this one.
    pub fn merge(&mut self, other: InnerDistAccum) {
        self.pairs.extend(other.pairs);
        self.distances.extend(other.distances);
        self.pair_num += other.pair_num;
    }
}

// ===================================================================
// Shared CIGAR helper
// ===================================================================

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
    pub bam_stat: Option<BamStatAccum>,
    pub infer_exp: Option<InferExpAccum>,
    pub read_dup: Option<ReadDupAccum>,
    pub read_dist: Option<ReadDistAccum>,
    pub junc_annot: Option<JuncAnnotAccum>,
    pub junc_sat: Option<JuncSatAccum>,
    pub inner_dist: Option<InnerDistAccum>,
    pub tin: Option<TinAccum>,
    /// preseq library complexity accumulator.
    pub preseq: Option<PreseqAccum>,
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
                Some(InferExpAccum::new(config.infer_experiment_sample_size))
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
                annotations
                    .and_then(|a| a.tin_index)
                    .map(|idx| TinAccum::new(idx, config.mapq_cut, config.tin_min_coverage))
            } else {
                None
            },
            preseq: if config.preseq_enabled {
                Some(PreseqAccum::new())
            } else {
                None
            },
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

        // junction_saturation: needs known junction chroms, uses uppercased chrom
        if let (Some(ref mut accum), Some(ref_chroms)) =
            (&mut self.junc_sat, &annotations.ref_chroms)
        {
            accum.process_read(
                record,
                chrom_upper,
                ref_chroms,
                config.min_intron,
                config.mapq_cut,
            );
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

        // preseq: counts unique fragment fingerprints for library complexity.
        //
        // Matches preseq v3.2.0 PE BAM behavior (load_counts_BAM_pe):
        //   - Only skip SECONDARY reads (flag 0x100). Note: preseq v3.2.0's
        //     `is_primary` check is `!(flag & 0x100)` — it does NOT filter
        //     supplementary alignments (0x800), so we include them.
        //   - Must also be mapped (not flag 0x4).
        //   - For PROPER PAIRS (flag 0x2): only count read1 — one per fragment.
        //     Fragment key = (chrom, frag_start, frag_end), matching preseq's
        //     `merge_mates()` which merges pairs by read name.
        //   - For non-proper-pair reads (including paired reads without 0x2):
        //     treat each read individually with key = (tid, pos, cigar_end).
        //     preseq pushes both mates as individual GenomicRegions in this case.
        //   - We compute frag_start/frag_end from pos + insert_size (TLEN) to
        //     avoid needing a read-name lookup buffer.
        //
        // NOTE: preseq v3.2.0 (used by nf-core/rnaseq) differs significantly
        // from preseq master branch which uses a simple aln_pos_pair template.
        // We match v3.2.0 behavior here.
        //
        // Unmapped reads (tid < 0) are additionally skipped inside process_read().
        if let Some(ref mut accum) = self.preseq {
            // preseq's is_primary = !(flag & 0x100) — only filters secondary,
            // not supplementary. Must also be mapped.
            if !record.is_secondary() && !record.is_unmapped() {
                if record.is_proper_pair() {
                    // Proper pair (flag 0x2): merge mates into one fragment.
                    // Only count read1 to avoid double-counting the pair.
                    if record.is_first_in_template() {
                        let tlen = record.insert_size();
                        let (frag_start, frag_end) = if tlen > 0 {
                            // Read1 is leftmost: frag spans pos to pos+tlen
                            (record.pos(), record.pos() + tlen)
                        } else if tlen < 0 {
                            // Read1 is rightmost: fragment end = cigar_end,
                            // fragment start = cigar_end + tlen (tlen is negative,
                            // abs(tlen) = frag_end - frag_start per SAM spec)
                            let cigar_end = record.cigar().end_pos();
                            (cigar_end + tlen, cigar_end)
                        } else {
                            // TLEN == 0: mate unmapped or same position
                            (record.pos(), record.cigar().end_pos())
                        };
                        accum.process_read(record.tid(), frag_start, frag_end);
                    }
                } else {
                    // Non-proper-pair or unpaired: count each read individually.
                    // preseq pushes these as individual GenomicRegions with
                    // key = (chrom, start, end) from CIGAR.
                    accum.process_read(record.tid(), record.pos(), record.cigar().end_pos());
                }
            }
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

impl BamStatAccum {
    /// Convert accumulated counters into a `BamStatResult` for output.
    pub fn into_result(self) -> BamStatResult {
        BamStatResult {
            // RSeQC bam_stat fields
            total_records: self.total_records,
            qc_failed: self.qc_failed,
            duplicates: self.duplicates,
            non_primary: self.non_primary,
            unmapped: self.unmapped,
            non_unique: self.non_unique,
            unique: self.unique,
            read_1: self.read_1,
            read_2: self.read_2,
            forward: self.forward,
            reverse: self.reverse,
            splice: self.splice,
            non_splice: self.non_splice,
            proper_pairs: self.proper_pairs,
            proper_pair_diff_chrom: self.proper_pair_diff_chrom,
            mapq_distribution: self.mapq_distribution,
            // samtools flagstat fields
            secondary: self.secondary,
            supplementary: self.supplementary,
            mapped: self.mapped,
            paired_flagstat: self.paired_flagstat,
            read1_flagstat: self.read1_flagstat,
            read2_flagstat: self.read2_flagstat,
            first_fragments: self.first_fragments,
            last_fragments: self.last_fragments,
            properly_paired: self.properly_paired,
            both_mapped: self.both_mapped,
            singletons: self.singletons,
            mate_diff_chr: self.mate_diff_chr,
            mate_diff_chr_mapq5: self.mate_diff_chr_mapq5,
            // samtools idxstats fields
            chrom_counts: self.chrom_counts,
            unplaced_unmapped: self.unplaced_unmapped,
            // samtools stats SN fields
            total_len: self.total_len,
            total_first_fragment_len: self.total_first_fragment_len,
            total_last_fragment_len: self.total_last_fragment_len,
            bases_mapped: self.bases_mapped,
            bases_mapped_cigar: self.bases_mapped_cigar,
            bases_duplicated: self.bases_duplicated,
            max_len: self.max_len,
            max_first_fragment_len: self.max_first_fragment_len,
            max_last_fragment_len: self.max_last_fragment_len,
            quality_sum: self.quality_sum,
            quality_count: self.quality_count,
            mismatches: self.mismatches,
            insert_size_sum: self.insert_size_sum,
            insert_size_sq_sum: self.insert_size_sq_sum,
            insert_size_count: self.insert_size_count,
            inward_pairs: self.inward_pairs,
            outward_pairs: self.outward_pairs,
            other_orientation: self.other_orientation,
            primary_paired: self.primary_paired,
            primary_count: self.primary_count,
            primary_mapped: self.primary_mapped,
            primary_duplicates: self.primary_duplicates,
            reads_mq0: self.reads_mq0,
            reads_mapped_and_paired: self.reads_mapped_and_paired,
        }
    }
}

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
        let seq_histogram = build_dup_histogram_from_hash(&self.seq_dup);
        ReadDuplicationResult {
            pos_histogram,
            seq_histogram,
        }
    }
}

/// Build duplication-level histogram from a count map.
/// Key = duplication level, Value = number of positions/sequences at that level.
fn build_dup_histogram(counts: &HashMap<u64, u64>) -> BTreeMap<u64, u64> {
    let mut histogram = BTreeMap::new();
    for &count in counts.values() {
        *histogram.entry(count).or_insert(0) += 1;
    }
    histogram
}

/// Build duplication-level histogram from the hash-based sequence map.
fn build_dup_histogram_from_hash(counts: &HashMap<u128, u64>) -> BTreeMap<u64, u64> {
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
    pub fn into_result(self) -> JunctionResults {
        JunctionResults {
            junctions: self.junction_counts,
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
        let mut rng = ChaCha8Rng::seed_from_u64(42);
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
    pub fn into_result(self, lower_bound: i64, upper_bound: i64, step: i64) -> InnerDistanceResult {
        let pairs: Vec<PairRecord> = self
            .pairs
            .into_iter()
            .map(|p| PairRecord {
                name: p.name,
                distance: p.distance,
                classification: p.classification,
            })
            .collect();
        let histogram = build_histogram(&self.distances, lower_bound, upper_bound, step);
        InnerDistanceResult {
            pairs,
            histogram,
            total_pairs: self.pair_num,
        }
    }
}
