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

use super::bam_stat::{BamStatResult, GcDepthBin};
use super::flagstat::FlagstatCounts;
use super::stats_coverage::{resolve_carries, ContigCovSummary, ContigCoverage};

/// Default GC-depth bin size in base pairs (matches upstream samtools default).
const GCD_BIN_SIZE: u64 = 20_000;
/// samtools stats insert size cap (`-i`, default 8000).
const MAX_INSERT_SIZE: u64 = 8000;
/// Longest indel in samtools stats' ID section (`nindels`, fixed at 300).
const MAX_INDEL_LEN: i64 = 300;
/// samtools stats read order: first fragment (unpaired, or read1 only).
const READ_ORDER_FIRST: u8 = 1;
/// samtools stats read order: last fragment (read2 only).
const READ_ORDER_LAST: u8 = 2;

/// samtools stats read order of a record: 1 = first, 2 = last, 0 or 3 = other.
fn read_order(flags: u16) -> u8 {
    if flags & BAM_FPAIRED == 0 {
        return READ_ORDER_FIRST;
    }
    u8::from(flags & BAM_FREAD1 != 0) + 2 * u8::from(flags & BAM_FREAD2 != 0)
}
use super::common::{self, KnownJunctionSet, ReferenceJunctions};
use super::infer_experiment::{GeneModel, InferExperimentResult};
use super::inner_distance::{
    build_histogram, ExonBitset, InnerDistanceResult, PairRecord, TranscriptTree,
};
use super::junction_annotation::{Junction, JunctionClass, JunctionResults};
use super::junction_saturation::SaturationResult;
use super::read_distribution::{ChromIntervals, ReadDistributionResult, RegionSets};
use super::read_duplication::ReadDuplicationResult;
use super::split_bam::{SplitBamAccum, SplitBamLookup};
use super::tin::TinAccum;
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
    /// Per-target BED blocks for split_bam (built per alignment file).
    pub split_bam: Option<&'a SplitBamLookup<'a>>,
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
    /// Whether split_bam BED classification is enabled.
    pub split_bam_enabled: bool,
}

// ===================================================================
// Per-tool accumulators
// ===================================================================

/// bam_stat accumulator — simple flag/MAPQ counting.
///
/// Also collects the additional counters needed for samtools-compatible
/// flagstat, idxstats, and stats output.
#[derive(Debug)]
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

    // --- samtools flagstat additional fields ---
    /// Secondary alignments (0x100) — counted independently of QC/dup.
    pub secondary: u64,
    /// Supplementary alignments (0x800) that are not also secondary, as in
    /// samtools (SECONDARY takes priority) — independent of QC/dup.
    pub supplementary: u64,
    /// All mapped records (not 0x4), regardless of QC/dup.
    pub mapped: u64,
    /// Paired reads (0x1), regardless of QC/dup.
    pub paired_flagstat: u64,
    /// Read1 in pair (0x40), regardless of QC/dup — for flagstat.
    pub read1_flagstat: u64,
    /// Read2 in pair (0x80), regardless of QC/dup — for flagstat.
    pub read2_flagstat: u64,
    /// samtools stats first fragments: primary reads with a sequence that are
    /// unpaired, or paired with only the read1 flag.
    pub first_fragments: u64,
    /// samtools stats last fragments: primary paired reads with a sequence and
    /// only the read2 flag.
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
    /// Sum of first fragment sequence lengths.
    pub total_first_fragment_len: u64,
    /// Sum of last fragment sequence lengths.
    pub total_last_fragment_len: u64,
    /// Sum of query lengths for mapped primary reads.
    pub bases_mapped: u64,
    /// Sum of M/I/=/X CIGAR operations for mapped non-secondary reads.
    pub bases_mapped_cigar: u64,
    /// Sum of query lengths for duplicate-flagged primary reads.
    pub bases_duplicated: u64,
    /// Maximum read length including hard clips (non-secondary reads).
    pub max_len: u64,
    /// Maximum first-fragment read length including hard clips.
    pub max_first_fragment_len: u64,
    /// Maximum last-fragment read length including hard clips.
    pub max_last_fragment_len: u64,
    /// Sum of base qualities of primary first and last fragments.
    pub quality_sum: f64,
    /// Number of bases contributing to quality_sum.
    pub quality_count: u64,
    /// Sum of NM tag values across mapped non-secondary reads (wrapping).
    pub mismatches: u64,
    /// Insert size with orientation: abs_tlen → [total, inward, outward, other].
    /// Both mates contribute (halved at output), capped at 8000.
    pub is_hist: HashMap<u64, [u64; 4]>,
    /// Inward-oriented pairs (FR).
    pub inward_pairs: u64,
    /// Outward-oriented pairs (RF).
    pub outward_pairs: u64,
    /// Other orientation pairs (FF, RR).
    pub other_orientation: u64,
    /// Total primary reads (non-secondary, non-supplementary).
    pub primary_count: u64,
    /// Primary mapped reads count (non-secondary, non-supplementary, !unmapped).
    pub primary_mapped: u64,
    /// Primary duplicate reads.
    pub primary_duplicates: u64,
    /// Primary mapped reads with MAPQ = 0 (matching upstream samtools stats).
    pub reads_mq0: u64,
    /// Primary mapped paired reads with a sequence whose mate is also mapped.
    pub reads_mapped_and_paired: u64,

    // --- samtools stats histogram/distribution fields ---
    /// Read length histogram (primary reads, length including hard clips): length → count.
    pub rl_hist: HashMap<u64, u64>,
    /// First fragment read length histogram: length → count.
    pub frl_hist: HashMap<u64, u64>,
    /// Last fragment read length histogram: length → count.
    pub lrl_hist: HashMap<u64, u64>,
    /// MAPQ histogram: primary, mapped, !qcfail, !dup (quality 0-255).
    pub mapq_hist: [u64; 256],
    /// Per-cycle quality for primary first fragments (read-oriented cycles).
    /// Outer: cycle index. Inner: quality value (0-255) → count.
    pub ffq: Vec<[u64; 256]>,
    /// Per-cycle quality for primary last fragments.
    pub lfq: Vec<[u64; 256]>,
    /// GC content step-function for first fragments, 200 bins (matching samtools ngc=200).
    /// Each bin i stores the number of reads with gc_count * 199 / seq_len <= i.
    pub gcf: [u64; 200],
    /// GC content step-function for last fragments, 200 bins.
    pub gcl: [u64; 200],
    /// Per-cycle base composition for first fragments (primary, mapped, !qcfail, !dup).
    /// [A, C, G, T, N, Other] per cycle.
    pub fbc: Vec<[u64; 6]>,
    /// Per-cycle base composition for last fragments.
    pub lbc: Vec<[u64; 6]>,
    /// Per-cycle base composition (read-oriented) for first fragments.
    /// Reverse strand reads contribute in reversed cycle order.
    pub fbc_ro: Vec<[u64; 6]>,
    /// Per-cycle base composition (read-oriented) for last fragments.
    pub lbc_ro: Vec<[u64; 6]>,
    /// Per-cycle base composition (reverse-complemented for reverse-strand reads,
    /// combined first+last fragments). Used for GCT output. [A, C, G, T] only.
    pub gcc_rc: Vec<[u64; 4]>,
    /// Total base counters for first fragments: [A, C, G, T, N].
    pub ftc: [u64; 5],
    /// Total base counters for last fragments: [A, C, G, T, N].
    pub ltc: [u64; 5],
    /// Indel distribution by size: length → [insertions, deletions].
    pub id_hist: HashMap<u64, [u64; 2]>,
    /// Indels per cycle: cycle → [ins_fwd, ins_rev, del_fwd, del_rev].
    pub ic: Vec<[u64; 4]>,
    /// CRC32 checksum sums: [names, sequences, qualities].
    /// Each is the wrapping u32 sum of per-read CRC32 values.
    pub chk: [u32; 3],
    /// Coverage distribution: depth → number of reference positions at that depth.
    /// Populated from a round buffer pileup during sorted BAM processing.
    pub cov_hist: HashMap<u32, u64>,
    /// Coverage round buffer of the contig being read.
    cov_contig: Option<ContigCoverage>,
    /// Finished contigs, for the slot samtools carries between contigs.
    cov_done: Vec<ContigCovSummary>,

    // --- samtools stats counters without a bam_stat/flagstat equivalent ---
    /// Primary reads with a sequence that are neither first nor last fragments
    /// (paired with both or neither of the read1/read2 flags).
    pub other_fragments: u64,
    /// Primary reads with a sequence and the QC-fail flag.
    pub stats_qcfail: u64,
    /// Primary unmapped reads with a sequence.
    pub stats_unmapped: u64,
    /// Primary reads with a sequence and the paired flag.
    pub stats_paired_tech: u64,
    /// Primary reads with a sequence, both mates mapped and the proper-pair flag.
    pub stats_properly_paired: u64,
    /// Primary reads with a sequence, both mates mapped, mate on another contig.
    pub stats_anomalous: u64,
    /// Duplicate-flagged non-secondary reads with a sequence.
    pub stats_dup_reads: u64,
    /// Sequence length of `stats_dup_reads`.
    pub stats_dup_bases: u64,

    /// samtools flagstat counters, split by QC-fail.
    pub flagstat: FlagstatCounts,

    // --- GC-depth (GCD section) fields ---
    /// GC-depth bins (one per `GCD_BIN_SIZE`-bp window) with their contig and start.
    gcd_bins: Vec<(i32, i64, GcDepthBin)>,
    /// Index in `gcd_bins` of the bin being filled.
    gcd_cur: Option<usize>,
    /// Start position of the current GCD bin.
    gcd_pos: i64,
    /// Chromosome tid of the current GCD bin.
    gcd_tid: i32,
}

impl Default for BamStatAccum {
    fn default() -> Self {
        Self {
            total_records: 0,
            qc_failed: 0,
            duplicates: 0,
            non_primary: 0,
            unmapped: 0,
            non_unique: 0,
            unique: 0,
            read_1: 0,
            read_2: 0,
            forward: 0,
            reverse: 0,
            splice: 0,
            non_splice: 0,
            proper_pairs: 0,
            proper_pair_diff_chrom: 0,
            secondary: 0,
            supplementary: 0,
            mapped: 0,
            paired_flagstat: 0,
            read1_flagstat: 0,
            read2_flagstat: 0,
            first_fragments: 0,
            last_fragments: 0,
            properly_paired: 0,
            both_mapped: 0,
            singletons: 0,
            mate_diff_chr: 0,
            mate_diff_chr_mapq5: 0,
            chrom_counts: HashMap::new(),
            unplaced_unmapped: 0,
            total_len: 0,
            total_first_fragment_len: 0,
            total_last_fragment_len: 0,
            bases_mapped: 0,
            bases_mapped_cigar: 0,
            bases_duplicated: 0,
            max_len: 0,
            max_first_fragment_len: 0,
            max_last_fragment_len: 0,
            quality_sum: 0.0,
            quality_count: 0,
            mismatches: 0,
            is_hist: HashMap::new(),
            inward_pairs: 0,
            outward_pairs: 0,
            other_orientation: 0,
            primary_count: 0,
            primary_mapped: 0,
            primary_duplicates: 0,
            reads_mq0: 0,
            reads_mapped_and_paired: 0,
            rl_hist: HashMap::new(),
            frl_hist: HashMap::new(),
            lrl_hist: HashMap::new(),
            mapq_hist: [0u64; 256],
            ffq: Vec::new(),
            lfq: Vec::new(),
            gcf: [0u64; 200],
            gcl: [0u64; 200],
            fbc: Vec::new(),
            lbc: Vec::new(),
            fbc_ro: Vec::new(),
            lbc_ro: Vec::new(),
            gcc_rc: Vec::new(),
            ftc: [0u64; 5],
            ltc: [0u64; 5],
            id_hist: HashMap::new(),
            ic: Vec::new(),
            chk: [0u32; 3],
            cov_hist: HashMap::new(),
            cov_contig: None,
            cov_done: Vec::new(),
            other_fragments: 0,
            stats_qcfail: 0,
            stats_unmapped: 0,
            stats_paired_tech: 0,
            stats_properly_paired: 0,
            stats_anomalous: 0,
            stats_dup_reads: 0,
            stats_dup_bases: 0,
            flagstat: FlagstatCounts::default(),
            gcd_bins: Vec::new(),
            gcd_cur: None,
            gcd_pos: -1,
            gcd_tid: -1,
        }
    }
}

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
        // samtools flagstat counters (count ALL records, no early returns).
        // `flagstat` is what the flagstat output uses; the individual fields
        // are kept for the library API.
        // =================================================================
        // SECONDARY takes priority over SUPPLEMENTARY, as in samtools
        // `flagstat_loop` and `stats.c`: a record with both flags is secondary
        // only, so `primary + secondary + supplementary == total`.
        if is_secondary {
            self.secondary += 1;
        } else if is_supplementary {
            self.supplementary += 1;
        }
        if is_mapped {
            self.mapped += 1;
        }
        self.flagstat.add(flags, tid, record.mtid(), mapq);
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
        // CHK checksums: computed on ALL reads (including secondary and
        // supplementary). Matches samtools stats.c update_checksum() which
        // is called before the secondary-read early return.
        // =================================================================
        {
            let qname = record.qname();
            let name_crc = crc32fast::hash(qname);
            self.chk[0] = self.chk[0].wrapping_add(name_crc);

            let seq_len = record.seq_len();
            if seq_len > 0 {
                // SAFETY: We access the raw BAM record data to compute CRC32
                // checksums matching samtools' approach. The pointer arithmetic
                // replicates htslib's bam_get_seq() macro:
                //   data + l_qname + (n_cigar << 2)
                // The seq_len > 0 guard above ensures sequence data exists.
                // The slice length seq_len.div_ceil(2) matches the BAM spec's
                // 4-bit encoded sequence format: (seq_len+1)/2 bytes.
                let seq_bytes = unsafe {
                    let inner = record.inner();
                    let data = inner.data;
                    let seq_offset =
                        inner.core.l_qname as isize + ((inner.core.n_cigar as isize) << 2);
                    let seq_nbytes = seq_len.div_ceil(2);
                    std::slice::from_raw_parts(data.offset(seq_offset), seq_nbytes)
                };
                let seq_crc = crc32fast::hash(seq_bytes);
                self.chk[1] = self.chk[1].wrapping_add(seq_crc);

                let qual = record.qual();
                let qual_crc = crc32fast::hash(qual);
                self.chk[2] = self.chk[2].wrapping_add(qual_crc);
            }
        }

        // Primary-read totals kept for the library API (`BamStatResult`).
        let seq_len = record.seq_len();
        if is_primary {
            self.primary_count += 1;
            if is_dup {
                self.primary_duplicates += 1;
                self.bases_duplicated += seq_len as u64;
            }
            if is_mapped {
                self.primary_mapped += 1;
            }
        }

        // =================================================================
        // samtools stats. Secondary records only count towards
        // "non-primary alignments" and the CHK checksums above, and records
        // without a sequence stop here, as in samtools `collect_stats`.
        // =================================================================
        if !is_secondary && seq_len > 0 {
            self.process_stats_read(record, flags, seq_len);
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

    /// samtools stats counters for one non-secondary record with a sequence.
    ///
    /// Follows `collect_stats` in samtools `stats.c` (1.22.1) with default
    /// options: no flag filters, no reference, coordinate-sorted input.
    fn process_stats_read(&mut self, record: &bam::Record, flags: u16, seq_len: usize) {
        use rust_htslib::bam::record::Cigar as C;

        let is_primary = flags & BAM_FSUPPLEMENTARY == 0;
        let is_unmapped = flags & BAM_FUNMAP != 0;
        let is_paired = flags & BAM_FPAIRED != 0;
        let mate_unmapped = flags & BAM_FMUNMAP != 0;
        let tid = record.tid();
        let pos = record.pos();
        let cigar = record.cigar();

        if flags & BAM_FDUP != 0 {
            self.stats_dup_reads += 1;
            self.stats_dup_bases += seq_len as u64;
        }

        let order = read_order(flags);

        // samtools `unclipped_length`: the sequence plus hard clips.
        let read_len = seq_len as u64
            + cigar
                .iter()
                .map(|op| match op {
                    C::HardClip(n) => u64::from(*n),
                    _ => 0,
                })
                .sum::<u64>();
        self.max_len = self.max_len.max(read_len);
        if order == READ_ORDER_FIRST {
            self.max_first_fragment_len = self.max_first_fragment_len.max(read_len);
        }
        if order == READ_ORDER_LAST {
            self.max_last_fragment_len = self.max_last_fragment_len.max(read_len);
        }

        // MAPQ histogram: !(UNMAP|SECOND|SUPPL|QCFAIL|DUP).
        if is_primary && !is_unmapped && flags & (BAM_FQCFAIL | BAM_FDUP) == 0 {
            self.mapq_hist[record.mapq() as usize] += 1;
        }

        // GC of the read for the GC-depth bins: samtools only computes it for
        // primary first/last fragments; other reads add depth with GC 0.
        let mut gc_count: u32 = 0;
        if is_primary {
            *self.rl_hist.entry(read_len).or_insert(0) += 1;
            if order == READ_ORDER_FIRST {
                *self.frl_hist.entry(read_len).or_insert(0) += 1;
            }
            if order == READ_ORDER_LAST {
                *self.lrl_hist.entry(read_len).or_insert(0) += 1;
            }
            gc_count = self.collect_orig_read_stats(record, flags, order, seq_len);
        }

        if is_unmapped {
            return;
        }

        self.count_indels(&cigar, flags & BAM_FREVERSE != 0, order, seq_len);

        // Insert size and pair orientation: primary, paired, both mapped.
        // Both mates are counted; output halves the counts.
        if is_primary && is_paired && !mate_unmapped {
            let mtid = record.mtid();
            let isize = record.insert_size().unsigned_abs().min(MAX_INSERT_SIZE);
            if isize > 0 || tid == mtid {
                let pos_fst = record.mpos() - pos;
                let is_fst: i64 = if flags & BAM_FREAD1 != 0 { 1 } else { -1 };
                let is_fwd: i64 = if flags & BAM_FREVERSE != 0 { -1 } else { 1 };
                let is_mfwd: i64 = if flags & BAM_FMREVERSE != 0 { -1 } else { 1 };

                // orientation index: 1 = inward, 2 = outward, 3 = other
                let orientation_idx = if is_fwd * is_mfwd > 0 {
                    self.other_orientation += 1;
                    3usize
                } else if is_fst * pos_fst > 0 {
                    if is_fst * is_fwd > 0 {
                        self.inward_pairs += 1;
                        1
                    } else {
                        self.outward_pairs += 1;
                        2
                    }
                } else if is_fst * pos_fst < 0 {
                    if is_fst * is_fwd > 0 {
                        self.outward_pairs += 1;
                        2
                    } else {
                        self.inward_pairs += 1;
                        1
                    }
                } else {
                    // exactly overlapping reads count as inward
                    self.inward_pairs += 1;
                    1
                };
                let entry = self.is_hist.entry(isize).or_insert([0; 4]);
                entry[0] += 1;
                entry[orientation_idx] += 1;
            }
        }

        // Mismatches from NM (`bam_aux2i`: signed, non-integer types count 0).
        if let Some(nm) = get_aux_int(record, b"NM") {
            self.mismatches = self.mismatches.wrapping_add(nm as u64);
        }

        // Bases mapped (cigar): M, I, = and X operations.
        self.bases_mapped_cigar += cigar
            .iter()
            .map(|op| match op {
                C::Match(n) | C::Ins(n) | C::Equal(n) | C::Diff(n) => u64::from(*n),
                _ => 0,
            })
            .sum::<u64>();

        // GC-depth bins (no reference): a new bin on the first read, a contig
        // change, or a read starting more than GCD_BIN_SIZE past the bin start.
        let new_bin = self.gcd_cur.is_none()
            || tid != self.gcd_tid
            || pos - self.gcd_pos > GCD_BIN_SIZE as i64;
        if new_bin {
            self.gcd_tid = tid;
            self.gcd_pos = pos;
            self.gcd_cur = Some(self.gcd_bins.len());
            self.gcd_bins
                .push((tid, pos, GcDepthBin { gc: 0.0, depth: 0 }));
        }
        if let Some((_, _, bin)) = self.gcd_cur.and_then(|i| self.gcd_bins.get_mut(i)) {
            bin.depth += 1;
            bin.gc += gc_count as f32 / seq_len as f32;
        }

        // Coverage distribution: one round buffer per contig (see stats_coverage).
        if self.cov_contig.as_ref().is_none_or(|c| c.tid() != tid) {
            self.finish_cov_contig();
            self.cov_contig = Some(ContigCoverage::new(tid));
        }
        if let Some(cov) = self.cov_contig.as_mut() {
            cov.ensure_read_len(read_len as usize);
            cov.flush(pos, &mut self.cov_hist);
            let mut p = pos;
            for op in cigar.iter() {
                match op {
                    C::Match(n) | C::Equal(n) | C::Diff(n) => {
                        cov.insert(p, p + i64::from(*n));
                        p += i64::from(*n);
                    }
                    C::Del(n) | C::RefSkip(n) => p += i64::from(*n),
                    _ => {}
                }
            }
        }
    }

    /// samtools `collect_orig_read_stats`: per-read statistics of a primary
    /// record with a sequence. Returns the read's G+C count (0 unless it is a
    /// first or last fragment).
    fn collect_orig_read_stats(
        &mut self,
        record: &bam::Record,
        flags: u16,
        order: u8,
        seq_len: usize,
    ) -> u32 {
        self.total_len += seq_len as u64;
        if flags & BAM_FQCFAIL != 0 {
            self.stats_qcfail += 1;
        }
        if flags & BAM_FPAIRED != 0 {
            self.stats_paired_tech += 1;
        }

        let mut gc_count: u32 = 0;
        if order == READ_ORDER_FIRST || order == READ_ORDER_LAST {
            let is_reverse = flags & BAM_FREVERSE != 0;
            let seq = record.seq();
            let quals = record.qual();

            let (qual_arr, base_arr, base_ro_arr, gc_arr, tc_arr) = if order == READ_ORDER_LAST {
                (
                    &mut self.lfq,
                    &mut self.lbc,
                    &mut self.lbc_ro,
                    &mut self.gcl,
                    &mut self.ltc,
                )
            } else {
                (
                    &mut self.ffq,
                    &mut self.fbc,
                    &mut self.fbc_ro,
                    &mut self.gcf,
                    &mut self.ftc,
                )
            };

            if seq_len > qual_arr.len() {
                qual_arr.resize(seq_len, [0u64; 256]);
            }
            if seq_len > base_arr.len() {
                base_arr.resize(seq_len, [0u64; 6]);
            }
            if seq_len > base_ro_arr.len() {
                base_ro_arr.resize(seq_len, [0u64; 6]);
            }
            if seq_len > self.gcc_rc.len() {
                self.gcc_rc.resize(seq_len, [0u64; 4]);
            }

            // BAM 4-bit encoding: A=1, C=2, G=4, T=8, N=15, others=0,3,5..14
            // BASE_IDX[nibble] → 0=A, 1=C, 2=G, 3=T, 4=N, 5=Other
            const BASE_IDX: [u8; 16] = [5, 0, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 4];
            // RC_IDX[base_idx] → reverse-complement base_idx (A↔T, C↔G)
            const RC_IDX: [u8; 6] = [3, 2, 1, 0, 4, 5];

            let mut qual_sum: u64 = 0;
            // Hoist the strand branch out of the inner loop.
            if !is_reverse {
                for i in 0..seq_len {
                    let q = quals[i];
                    qual_sum += u64::from(q);
                    qual_arr[i][usize::from(q)] += 1;

                    let base_idx = BASE_IDX[seq.encoded_base(i) as usize] as usize;
                    base_arr[i][base_idx] += 1;
                    base_ro_arr[i][base_idx] += 1;
                    if base_idx < 4 {
                        self.gcc_rc[i][base_idx] += 1;
                    }
                    if base_idx == 1 || base_idx == 2 {
                        gc_count += 1;
                    }
                    if base_idx < 5 {
                        tc_arr[base_idx] += 1;
                    }
                }
            } else {
                for i in 0..seq_len {
                    let ro_cycle = seq_len - 1 - i;
                    let q = quals[i];
                    qual_sum += u64::from(q);
                    qual_arr[ro_cycle][usize::from(q)] += 1;

                    let base_idx = BASE_IDX[seq.encoded_base(i) as usize] as usize;
                    base_arr[i][base_idx] += 1;
                    base_ro_arr[ro_cycle][base_idx] += 1;
                    if base_idx < 4 {
                        self.gcc_rc[ro_cycle][RC_IDX[base_idx] as usize] += 1;
                    }
                    if base_idx == 1 || base_idx == 2 {
                        gc_count += 1;
                    }
                    if base_idx < 5 {
                        tc_arr[base_idx] += 1;
                    }
                }
            }
            self.quality_sum += qual_sum as f64;
            self.quality_count += seq_len as u64;

            // GC content: cumulative step function with ngc = 200 bins.
            let ngc: usize = 200;
            let gc_idx_min = gc_count as usize * (ngc - 1) / seq_len;
            let gc_idx_max = ((gc_count as usize + 1) * (ngc - 1) / seq_len).min(ngc - 1);
            for item in gc_arr.iter_mut().take(gc_idx_max).skip(gc_idx_min) {
                *item += 1;
            }
        }

        match order {
            READ_ORDER_FIRST => {
                self.first_fragments += 1;
                self.total_first_fragment_len += seq_len as u64;
            }
            READ_ORDER_LAST => {
                self.last_fragments += 1;
                self.total_last_fragment_len += seq_len as u64;
            }
            _ => self.other_fragments += 1,
        }

        if flags & BAM_FUNMAP != 0 {
            self.stats_unmapped += 1;
        } else {
            self.bases_mapped += seq_len as u64;
            if record.mapq() == 0 {
                self.reads_mq0 += 1;
            }
            if flags & BAM_FPAIRED != 0 && flags & BAM_FMUNMAP == 0 {
                self.reads_mapped_and_paired += 1;
                if flags & BAM_FPROPER_PAIR != 0 {
                    self.stats_properly_paired += 1;
                }
                if record.tid() != record.mtid() {
                    self.stats_anomalous += 1;
                }
            }
        }
        gc_count
    }

    /// samtools `count_indels`: indel size distribution (ID) and indels per
    /// read-oriented cycle (IC) of a mapped, non-secondary record.
    fn count_indels(
        &mut self,
        cigar: &rust_htslib::bam::record::CigarStringView,
        is_reverse: bool,
        order: u8,
        seq_len: usize,
    ) {
        use rust_htslib::bam::record::Cigar as C;

        let read_len = seq_len as i64;
        let mut icycle: i64 = 0;
        for op in cigar.iter() {
            let ncig = i64::from(op.len());
            if ncig == 0 {
                continue; // "curiously, this can happen: 0D"
            }
            match op {
                C::Ins(_) => {
                    let idx = if is_reverse {
                        read_len - icycle - ncig
                    } else {
                        icycle
                    };
                    if idx >= 0 {
                        self.add_indel_cycle(idx as usize, order, 0);
                    }
                    icycle += ncig;
                    if ncig <= MAX_INDEL_LEN {
                        self.id_hist.entry(ncig as u64).or_insert([0; 2])[0] += 1;
                    }
                }
                C::Del(_) => {
                    let idx = if is_reverse {
                        read_len - icycle - 1
                    } else {
                        icycle - 1
                    };
                    if idx < 0 {
                        continue; // samtools discards these "meaningless deletions"
                    }
                    self.add_indel_cycle(idx as usize, order, 2);
                    if ncig <= MAX_INDEL_LEN {
                        self.id_hist.entry(ncig as u64).or_insert([0; 2])[1] += 1;
                    }
                }
                C::RefSkip(_) | C::HardClip(_) | C::Pad(_) => {}
                _ => icycle += ncig,
            }
        }
    }

    /// Count an insertion (`kind` 0) or deletion (`kind` 2) at cycle `idx`
    /// for first (`+0`) or last (`+1`) fragments.
    fn add_indel_cycle(&mut self, idx: usize, order: u8, kind: usize) {
        let col = match order {
            READ_ORDER_FIRST => kind,
            READ_ORDER_LAST => kind + 1,
            _ => return,
        };
        if idx >= self.ic.len() {
            self.ic.resize(idx + 1, [0u64; 4]);
        }
        self.ic[idx][col] += 1;
    }

    /// End the coverage round buffer of the current contig, if any.
    fn finish_cov_contig(&mut self) {
        if let Some(cov) = self.cov_contig.take() {
            self.cov_done.push(cov.finish(&mut self.cov_hist));
        }
    }

    /// Merge another accumulator into this one.
    pub fn merge(&mut self, mut other: BamStatAccum) {
        // End the other worker's open contig so its coverage is complete.
        other.finish_cov_contig();

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
        self.mismatches = self.mismatches.wrapping_add(other.mismatches);
        for (isize_val, counts) in other.is_hist {
            let entry = self.is_hist.entry(isize_val).or_insert([0; 4]);
            for i in 0..4 {
                entry[i] += counts[i];
            }
        }
        self.inward_pairs += other.inward_pairs;
        self.outward_pairs += other.outward_pairs;
        self.other_orientation += other.other_orientation;
        self.primary_count += other.primary_count;
        self.primary_mapped += other.primary_mapped;
        self.primary_duplicates += other.primary_duplicates;
        self.reads_mq0 += other.reads_mq0;
        self.reads_mapped_and_paired += other.reads_mapped_and_paired;

        // Histogram/distribution fields
        for (len, count) in other.rl_hist {
            *self.rl_hist.entry(len).or_insert(0) += count;
        }
        for (len, count) in other.frl_hist {
            *self.frl_hist.entry(len).or_insert(0) += count;
        }
        for (len, count) in other.lrl_hist {
            *self.lrl_hist.entry(len).or_insert(0) += count;
        }
        for i in 0..256 {
            self.mapq_hist[i] += other.mapq_hist[i];
        }

        // Per-cycle quality arrays (FFQ/LFQ)
        merge_vec_arrays(&mut self.ffq, other.ffq);
        merge_vec_arrays(&mut self.lfq, other.lfq);

        // GC content distributions (200 bins)
        for i in 0..200 {
            self.gcf[i] += other.gcf[i];
            self.gcl[i] += other.gcl[i];
        }

        // Per-cycle base composition (FBC/LBC and read-oriented)
        merge_vec_arrays(&mut self.fbc, other.fbc);
        merge_vec_arrays(&mut self.lbc, other.lbc);
        merge_vec_arrays(&mut self.fbc_ro, other.fbc_ro);
        merge_vec_arrays(&mut self.lbc_ro, other.lbc_ro);
        merge_vec_arrays(&mut self.gcc_rc, other.gcc_rc);

        // Total base counters
        for i in 0..5 {
            self.ftc[i] += other.ftc[i];
            self.ltc[i] += other.ltc[i];
        }

        // Indel distribution
        for (len, counts) in other.id_hist {
            let entry = self.id_hist.entry(len).or_insert([0; 2]);
            entry[0] += counts[0];
            entry[1] += counts[1];
        }

        // Indels per cycle
        merge_vec_arrays(&mut self.ic, other.ic);

        // CHK checksums (wrapping u32 addition)
        for i in 0..3 {
            self.chk[i] = self.chk[i].wrapping_add(other.chk[i]);
        }

        // COV histogram (additive merge)
        for (depth, count) in other.cov_hist {
            *self.cov_hist.entry(depth).or_insert(0) += count;
        }

        // Finished contigs; the carried slots are resolved in into_result().
        self.cov_done.append(&mut other.cov_done);

        // samtools stats-only counters
        self.other_fragments += other.other_fragments;
        self.stats_qcfail += other.stats_qcfail;
        self.stats_unmapped += other.stats_unmapped;
        self.stats_paired_tech += other.stats_paired_tech;
        self.stats_properly_paired += other.stats_properly_paired;
        self.stats_anomalous += other.stats_anomalous;
        self.stats_dup_reads += other.stats_dup_reads;
        self.stats_dup_bases += other.stats_dup_bases;
        self.flagstat.merge(&other.flagstat);

        // GCD bins (whole contigs per worker; put back in BAM order in into_result()).
        self.gcd_bins.append(&mut other.gcd_bins);
    }
}

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
///
/// The observation list is shuffled before sampling, so its order must not
/// depend on the thread count. It is kept in BAM record order (then intron
/// order within each read), which is also RSeQC's order. Observations are
/// therefore kept in per-chromosome segments, in the order the records were
/// read. A worker reads each chromosome in one contiguous `fetch`, and
/// [`merge`](Self::merge) sorts the segments by reference ID, which is file
/// order for a coordinate-sorted (indexed) BAM, whatever the thread count.
/// The sequential path never merges, so it keeps file order for any BAM.
#[derive(Debug, Default)]
pub struct JuncSatAccum {
    /// Junction observation keys (hashed `"CHROM:start-end"`), grouped in
    /// runs of consecutive records with the same reference ID.
    segments: Vec<(i32, Vec<u64>)>,
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
            self.push_observation(record.tid(), chrom_upper, istart, iend);
        }
    }

    /// Append one junction observation from a record on reference `tid`.
    fn push_observation(&mut self, tid: i32, chrom_upper: &str, istart: u64, iend: u64) {
        let h = hash_junction_key(chrom_upper, istart, iend);
        match self.segments.last_mut() {
            Some((seg_tid, obs)) if *seg_tid == tid => obs.push(h),
            _ => self.segments.push((tid, vec![h])),
        }
        // Track the string key only for the first occurrence of each hash
        self.unique_keys
            .entry(h)
            .or_insert_with(|| format!("{}:{}-{}", chrom_upper, istart, iend));
    }

    /// Merge another accumulator into this one.
    ///
    /// Segments are stably sorted by reference ID, so the merged order is the
    /// BAM file order no matter how chromosomes were split between workers.
    pub fn merge(&mut self, other: JuncSatAccum) {
        self.segments.extend(other.segments);
        self.segments.sort_by_key(|(tid, _)| *tid);
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
    /// split_bam accumulator (`None` when disabled).
    pub split_bam: Option<SplitBamAccum>,
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
            split_bam: None,
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
            split_bam: if config.split_bam_enabled {
                Some(SplitBamAccum::default())
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

        // split_bam: classifies EVERY record (incl. unmapped, QC-fail,
        // secondary, supplementary), keyed by the alignment's own tid so BED
        // names match alignment names, not GTF-mapped names (as RSeQC does).
        if let (Some(ref mut accum), Some(lookup)) = (&mut self.split_bam, annotations.split_bam) {
            accum.process_read(record, lookup);
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
        if let (Some(ref mut a), Some(b)) = (&mut self.split_bam, other.split_bam) {
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
// Merge helpers for Vec<[u64; N]> per-cycle arrays
// ===================================================================

/// Merge two `Vec<[u64; N]>` arrays element-wise, extending target if shorter.
fn merge_vec_arrays<const N: usize>(target: &mut Vec<[u64; N]>, source: Vec<[u64; N]>) {
    if source.len() > target.len() {
        target.resize(source.len(), [0u64; N]);
    }
    for (i, arr) in source.into_iter().enumerate() {
        for j in 0..N {
            target[i][j] += arr[j];
        }
    }
}

// ===================================================================
// Converter methods: accumulator → result types for output functions
// ===================================================================

impl BamStatAccum {
    /// Convert accumulated counters into a `BamStatResult` for output.
    pub fn into_result(mut self) -> BamStatResult {
        // End the last contig and resolve the slots carried between contigs.
        self.finish_cov_contig();
        resolve_carries(std::mem::take(&mut self.cov_done), &mut self.cov_hist);

        // GC-depth bins in BAM order (contig, then position): samtools treats
        // the last bin specially, so the order must not depend on threads.
        self.gcd_bins.sort_by_key(|&(tid, pos, _)| (tid, pos));
        let gcd_bins = self.gcd_bins.into_iter().map(|(_, _, bin)| bin).collect();

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
            is_hist: self.is_hist,
            inward_pairs: self.inward_pairs,
            outward_pairs: self.outward_pairs,
            other_orientation: self.other_orientation,
            primary_count: self.primary_count,
            primary_mapped: self.primary_mapped,
            primary_duplicates: self.primary_duplicates,
            reads_mq0: self.reads_mq0,
            reads_mapped_and_paired: self.reads_mapped_and_paired,
            // Histogram/distribution fields
            rl_hist: self.rl_hist,
            frl_hist: self.frl_hist,
            lrl_hist: self.lrl_hist,
            mapq_hist: self.mapq_hist,
            ffq: self.ffq,
            lfq: self.lfq,
            gcf: self.gcf,
            gcl: self.gcl,
            fbc: self.fbc,
            lbc: self.lbc,
            fbc_ro: self.fbc_ro,
            lbc_ro: self.lbc_ro,
            gcc_rc: self.gcc_rc,
            ftc: self.ftc,
            ltc: self.ltc,
            id_hist: self.id_hist,
            ic: self.ic,
            chk: self.chk,
            cov_hist: self.cov_hist,
            gcd_bins,
            other_fragments: self.other_fragments,
            stats_qcfail: self.stats_qcfail,
            stats_unmapped: self.stats_unmapped,
            stats_paired_tech: self.stats_paired_tech,
            stats_properly_paired: self.stats_properly_paired,
            stats_anomalous: self.stats_anomalous,
            stats_dup_reads: self.stats_dup_reads,
            stats_dup_bases: self.stats_dup_bases,
            flagstat: self.flagstat,
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
    /// Follows RSeQC's `ParseBAM.saturation_junction()`, except that the shuffle
    /// uses RustQC's own seeded RNG, so the intermediate points are
    /// reproducible but do not match RSeQC's (unseeded) random draws. The
    /// x-axis values and the 100% point are deterministic and match RSeQC:
    ///
    /// 1. Drop observations on chromosomes without an annotated transcript
    ///    (RSeQC skips reads whose chromosome is not in the BED gene model).
    /// 2. Shuffle the observations (in BAM file order) with a seeded ChaCha8 RNG.
    /// 3. For each sampling percentage `p` in
    ///    `range(sample_start, sample_end, sample_step) + [100]`, add the
    ///    observations from `int(N * ((p - sample_step) / 100.0))` (clamped
    ///    to 0) to `int(N * (p / 100.0))` to the running junction counts,
    ///    using the same floating-point index arithmetic. The windows can
    ///    overlap or leave gaps for some parameter combinations, as in RSeQC.
    ///
    /// Known/novel/all junction counts are maintained incrementally instead
    /// of re-scanning all unique junctions at every sampling step.
    pub fn into_result(
        self,
        known_junctions: &KnownJunctionSet,
        sample_start: u32,
        sample_end: u32,
        sample_step: u32,
        min_coverage: u32,
        seed: u64,
    ) -> SaturationResult {
        // Hashes of known junctions, and of junctions on annotated chromosomes.
        let mut known_hashes: HashSet<u64> = HashSet::new();
        let mut kept_hashes: HashSet<u64> = HashSet::new();
        for (&h, key) in &self.unique_keys {
            if known_junctions.junctions.contains(key) {
                known_hashes.insert(h);
            }
            let chrom = key.rsplit_once(':').map_or(key.as_str(), |(c, _)| c);
            if known_junctions.chroms.contains(chrom) {
                kept_hashes.insert(h);
            }
        }

        // Phase 1: the observation list in RSeQC's order (BAM file order).
        let mut observations: Vec<u64> = self
            .segments
            .into_iter()
            .flat_map(|(_, obs)| obs)
            .filter(|h| kept_hashes.contains(h))
            .collect();

        // Phase 2: deterministic shuffle
        {
            use rand::seq::SliceRandom;
            use rand::SeedableRng;
            observations.shuffle(&mut rand_chacha::ChaCha8Rng::seed_from_u64(seed));
        }

        // Phase 3: sampling, with RSeQC's percentage series and index arithmetic.
        let percentages = sampling_percentages(sample_start, sample_end, sample_step);
        let total = observations.len();
        let index_at = |pct: i64| -> usize {
            // Python: int(SR_num * (pct / 100.0)), truncating toward zero.
            let idx = (total as f64 * (pct as f64 / 100.0)) as i64;
            idx.clamp(0, total as i64) as usize
        };
        // A junction counts as known once seen `max(min_coverage, 1)` times
        // (RSeQC: `uniqSpliceSites[sj] >= recur` over junctions seen at least once).
        let known_threshold = min_coverage.max(1);

        let mut junction_counts: HashMap<u64, u32> = HashMap::new();
        let mut known_counts = Vec::with_capacity(percentages.len());
        let mut novel_counts = Vec::with_capacity(percentages.len());
        let mut all_counts = Vec::with_capacity(percentages.len());
        let mut running_known: usize = 0;
        let mut running_novel: usize = 0;

        for &pct in &percentages {
            let index_st = index_at(i64::from(pct) - i64::from(sample_step));
            let index_end = index_at(i64::from(pct));
            for &obs_hash in observations.get(index_st..index_end).unwrap_or(&[]) {
                let count = junction_counts.entry(obs_hash).or_insert(0);
                *count += 1;
                if known_hashes.contains(&obs_hash) {
                    if *count == known_threshold {
                        running_known += 1;
                    }
                } else if *count == 1 {
                    running_novel += 1;
                }
            }

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

/// RSeQC's sampling percentages: `list(range(start, end, step))` followed by 100.
///
/// A zero step would make Python's `range()` raise; here it yields only `[100]`.
fn sampling_percentages(sample_start: u32, sample_end: u32, sample_step: u32) -> Vec<u32> {
    let mut percentages: Vec<u32> = if sample_step == 0 {
        Vec::new()
    } else {
        (sample_start..sample_end)
            .step_by(sample_step as usize)
            .collect()
    };
    percentages.push(100);
    percentages
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sampling_percentages_match_python_range() {
        // list(range(5, 100, 5)) + [100]
        let default: Vec<u32> = (1..=20).map(|i| i * 5).collect();
        assert_eq!(sampling_percentages(5, 100, 5), default);
        // range() excludes the end, then RSeQC appends 100.
        assert_eq!(sampling_percentages(10, 50, 10), vec![10, 20, 30, 40, 100]);
        assert_eq!(sampling_percentages(1, 100, 5).last(), Some(&100));
        assert_eq!(sampling_percentages(1, 100, 5)[19], 96);
        assert_eq!(sampling_percentages(50, 50, 5), vec![100]);
        assert_eq!(sampling_percentages(5, 100, 0), vec![100]);
    }

    /// Observations for the saturation tests: 40 on CHR1 (tid 0) interleaved
    /// with 4 on CHRX (tid 1), which is not annotated.
    fn test_observations() -> Vec<(i32, &'static str, u64)> {
        let mut out = Vec::new();
        for i in 0..40u64 {
            if i % 11 == 3 {
                out.push((1, "CHRX", 1));
            }
            out.push((0, "CHR1", (i * 7) % 9 * 100));
        }
        out
    }

    fn test_known() -> KnownJunctionSet {
        KnownJunctionSet {
            junctions: ["CHR1:0-50", "CHR1:100-150", "CHR1:200-250", "CHR1:700-750"]
                .iter()
                .map(|s| s.to_string())
                .collect(),
            chroms: ["CHR1".to_string()].into_iter().collect(),
        }
    }

    fn push(accum: &mut JuncSatAccum, tid: i32, chrom: &str, start: u64) {
        let end = if chrom == "CHRX" { 100 } else { start + 50 };
        accum.push_observation(tid, chrom, start, end);
    }

    fn filled(obs: &[(i32, &str, u64)]) -> JuncSatAccum {
        let mut accum = JuncSatAccum::default();
        for &(tid, chrom, st) in obs {
            push(&mut accum, tid, chrom, st);
        }
        accum
    }

    #[test]
    fn test_saturation_final_point_and_axis() {
        // The x axis and (when the windows cover the whole list) the 100%
        // point do not depend on the shuffle.
        let known = test_known();
        for seed in [1, 7, 42] {
            for (start, end, step) in [(5, 100, 5), (1, 100, 5), (10, 100, 10)] {
                let r = filled(&test_observations()).into_result(&known, start, end, step, 1, seed);
                assert_eq!(r.percentages, sampling_percentages(start, end, step));
                // All 40 CHR1 observations are used (CHRX is dropped): 9
                // distinct junctions, 4 of them known.
                assert_eq!(*r.all_counts.last().unwrap(), 9);
                assert_eq!(*r.known_counts.last().unwrap(), 4);
                assert_eq!(*r.novel_counts.last().unwrap(), 5);
            }
        }
    }

    #[test]
    fn test_saturation_deterministic_and_seeded() {
        let known = test_known();
        let run = |seed| filled(&test_observations()).into_result(&known, 5, 100, 5, 1, seed);
        assert_eq!(run(3).all_counts, run(3).all_counts);
    }

    #[test]
    fn test_saturation_merge_order_independent() {
        // A coordinate-sorted layout: tid 0 records, then tid 1, then tid 2.
        let obs: Vec<(i32, &str, u64)> = (0..30u64)
            .map(|i| ((i / 10) as i32, "CHR1", (i * 13) % 7 * 100))
            .collect();
        let known = test_known();
        let run = |order: &[i32]| {
            let mut merged = JuncSatAccum::default();
            for &tid in order {
                let part: Vec<_> = obs.iter().copied().filter(|o| o.0 == tid).collect();
                merged.merge(filled(&part));
            }
            merged.into_result(&known, 5, 100, 5, 1, 3)
        };
        let expected = filled(&obs).into_result(&known, 5, 100, 5, 1, 3);
        for order in [[0, 1, 2], [2, 0, 1], [1, 2, 0]] {
            let r = run(&order);
            assert_eq!(r.all_counts, expected.all_counts);
            assert_eq!(r.known_counts, expected.known_counts);
            assert_eq!(r.novel_counts, expected.novel_counts);
        }
    }
}
