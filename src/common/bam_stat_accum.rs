//! Read-level alignment statistics accumulator.
//!
//! [`BamStatAccum`] gathers, in a single pass over the records, every counter
//! consumed by RSeQC `bam_stat` and by the samtools-compatible `stats`,
//! `flagstat` and `idxstats` writers. It needs no annotation and no library
//! protocol, so both the `rna` and `dna` pipelines drive the same struct: each
//! parallel worker owns one, and they are merged before conversion.

use std::collections::HashMap;

use rust_htslib::bam;

use crate::common::bam_flags::*;
use crate::common::bam_stat::{BamStatResult, GcDepthBin};

/// Default GC-depth bin size in base pairs (matches upstream samtools default).
const GCD_BIN_SIZE: u64 = 20_000;

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
    /// Insert size with orientation: abs_tlen → [total, inward, outward, other].
    /// Only one mate per pair contributes (upstream mate), capped at 8000.
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
    /// Primary non-QC-fail mapped paired reads where mate is also mapped.
    pub reads_mapped_and_paired: u64,

    // --- samtools stats histogram/distribution fields ---
    /// Read length histogram (all primary reads): length → count.
    pub rl_hist: HashMap<u64, u64>,
    /// First fragment read length histogram: length → count.
    pub frl_hist: HashMap<u64, u64>,
    /// Last fragment read length histogram: length → count.
    pub lrl_hist: HashMap<u64, u64>,
    /// MAPQ histogram: primary, mapped, !qcfail, !dup (quality 0-255).
    pub mapq_hist: [u64; 256],
    /// Per-cycle quality for first fragments (primary, mapped, !qcfail, !dup).
    /// Outer: cycle index. Inner: quality value → count (64 buckets covers Q0-Q63).
    pub ffq: Vec<[u64; 64]>,
    /// Per-cycle quality for last fragments.
    pub lfq: Vec<[u64; 64]>,
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
    /// Circular buffer for coverage pileup, matching upstream samtools design.
    /// `cov_buf[cov_buf_idx]` corresponds to reference position `cov_buf_pos`.
    /// The buffer grows dynamically to accommodate `max_read_length * 5`.
    cov_buf: Vec<u32>,
    /// Index into `cov_buf` corresponding to `cov_buf_pos`.
    cov_buf_idx: usize,
    /// Reference position of the element at `cov_buf[cov_buf_idx]`.
    cov_buf_pos: i64,
    /// Current chromosome tid for round buffer tracking.
    cov_buf_tid: i32,

    // --- GC-depth (GCD section) fields ---
    /// Accumulated GC-depth bins (one per `GCD_BIN_SIZE`-bp genomic window).
    gcd_bins: Vec<GcDepthBin>,
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
            cov_buf: vec![0u32; 1500], // matches upstream samtools: nbases * 5 = 300 * 5
            cov_buf_idx: 0,
            cov_buf_pos: 0,
            cov_buf_tid: -1,
            gcd_bins: Vec::new(),
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

        // Track gc_count from the primary-read per-cycle loop so the GCD
        // section below can reuse it without re-scanning the sequence.
        let mut primary_gc_count: u64 = 0;

        // =================================================================
        // samtools stats SN counters (primary reads only)
        // =================================================================
        if is_primary {
            self.primary_count += 1;
            let seq_len = record.seq_len() as u64;
            let mate_unmapped = flags & BAM_FMUNMAP != 0;

            self.total_len += seq_len;
            let is_last_fragment = is_paired && flags & BAM_FREAD2 != 0;
            if is_last_fragment {
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

            // RL/FRL/LRL: read length histograms (all primary reads)
            *self.rl_hist.entry(seq_len).or_insert(0) += 1;
            if is_last_fragment {
                *self.lrl_hist.entry(seq_len).or_insert(0) += 1;
            } else {
                *self.frl_hist.entry(seq_len).or_insert(0) += 1;
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

                // samtools stats: reads MQ0 counts primary mapped reads with MAPQ=0
                // (upstream stats.c: MQ0 is counted inside collect_orig_read_stats,
                // which is only called for IS_ORIGINAL reads = non-secondary, non-supplementary)
                if record.mapq() == 0 {
                    self.reads_mq0 += 1;
                }

                // NOTE: bases_mapped_cigar is now computed in the IC/ID CIGAR
                // loop below (for all mapped non-secondary reads) to avoid a
                // separate full CIGAR traversal here.

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

                // Insert size + orientation for paired primary reads where both
                // mates are mapped. Matches samtools stats gate:
                //   IS_PAIRED_AND_MAPPED && IS_ORIGINAL
                //   if (isize > 0 || tid == mtid)
                // Both mates contribute; samtools divides by 2 at output.
                // We do the same in write_insert_size() and the SN section.
                if is_paired && !mate_unmapped {
                    let tid = record.tid();
                    let mtid = record.mtid();
                    let tlen = record.insert_size();
                    let abs_tlen = tlen.unsigned_abs();

                    if abs_tlen > 0 || tid == mtid {
                        let pos = record.pos();
                        let mpos = record.mpos();

                        // Compute orientation (only meaningful for same-chromosome)
                        let pos_fst = mpos - pos;
                        let is_fst: i64 = if flags & BAM_FREAD1 != 0 { 1 } else { -1 };
                        let is_fwd: i64 = if flags & BAM_FREVERSE != 0 { -1 } else { 1 };
                        let is_mfwd: i64 = if flags & BAM_FMREVERSE != 0 { -1 } else { 1 };

                        // orientation_idx: 1=inward, 2=outward, 3=other
                        let orientation_idx = if is_fwd * is_mfwd > 0 {
                            self.other_orientation += 1;
                            3usize
                        } else if is_fst * pos_fst > 0 {
                            if is_fst * is_fwd > 0 {
                                self.inward_pairs += 1;
                                1usize
                            } else {
                                self.outward_pairs += 1;
                                2usize
                            }
                        } else if is_fst * pos_fst < 0 {
                            if is_fst * is_fwd > 0 {
                                self.outward_pairs += 1;
                                2usize
                            } else {
                                self.inward_pairs += 1;
                                1usize
                            }
                        } else {
                            self.inward_pairs += 1;
                            1usize
                        };

                        if abs_tlen > 0 {
                            // Cap at MAX_INSERT_SIZE (8000), matching
                            // samtools stats which accumulates overflow
                            // into the cap bucket.
                            let capped = abs_tlen.min(8000);
                            let entry = self.is_hist.entry(capped).or_insert([0; 4]);
                            entry[0] += 1; // total
                            entry[orientation_idx] += 1;
                        }
                    }
                }
            }

            // Average quality for primary non-QC-fail reads.
            // Upstream samtools stats computes per-BASE quality average:
            // sum of all individual base qualities / total bases.
            // (Not a per-read average of averages.)
            if !is_qcfail {
                let quals = record.qual();
                if !quals.is_empty() {
                    let base_qual_sum: f64 = quals.iter().map(|&q| f64::from(q)).sum::<f64>();
                    self.quality_sum += base_qual_sum;
                    self.quality_count += quals.len() as u64;
                }
            }

            // =============================================================
            // MAPQ histogram: primary + mapped + !qcfail + !dup
            // (matches samtools stats.c:1239 five-flag exclusion)
            // =============================================================
            if is_mapped && !is_qcfail && !is_dup {
                self.mapq_hist[mapq as usize] += 1;
            }

            // =============================================================
            // Per-cycle quality & base composition histograms:
            // FFQ/LFQ, FBC/LBC, GCF/GCL, FTC/LTC, FBC_RO/LBC_RO
            //
            // Upstream samtools stats includes duplicates, unmapped, and
            // qcfail reads in these histograms (collect_orig_read_stats
            // has no such checks). Only secondary+supplementary are
            // excluded (via IS_ORIGINAL), which is already handled by
            // the outer is_primary guard.
            // =============================================================
            {
                let is_reverse = flags & BAM_FREVERSE != 0;

                let seq = record.seq();
                let quals = record.qual();
                let read_len = seq.len();

                // Determine which arrays to use (first vs last fragment)
                // If paired: read2 = last, read1 = first. If SE: all = first.
                let (qual_arr, base_arr, base_ro_arr, gc_arr, tc_arr) = if is_last_fragment {
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

                // Ensure per-cycle arrays are large enough
                if read_len > qual_arr.len() {
                    qual_arr.resize(read_len, [0u64; 64]);
                }
                if read_len > base_arr.len() {
                    base_arr.resize(read_len, [0u64; 6]);
                }
                if read_len > base_ro_arr.len() {
                    base_ro_arr.resize(read_len, [0u64; 6]);
                }
                if read_len > self.gcc_rc.len() {
                    self.gcc_rc.resize(read_len, [0u64; 4]);
                }

                let mut gc_count: u64 = 0;

                // Pre-built lookup tables for the per-cycle inner loop,
                // avoiding branches and match overhead on every base.
                //
                // BAM 4-bit encoding: A=1, C=2, G=4, T=8, N=15, others=0,3,5..14
                // BASE_IDX[nibble] → 0=A, 1=C, 2=G, 3=T, 4=N, 5=Other
                const BASE_IDX: [u8; 16] = [5, 0, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 4];
                // RC_IDX[base_idx] → reverse-complement base_idx (A↔T, C↔G)
                // Only meaningful for base_idx 0-3 (ACGT). Index 4/5 not used.
                const RC_IDX: [u8; 6] = [3, 2, 1, 0, 4, 5]; // A→T, C→G, G→C, T→A

                // Hoist the is_reverse branch outside the inner loop so the
                // compiler can version the loop and potentially auto-vectorize
                // each variant independently.
                if !is_reverse {
                    for i in 0..read_len {
                        let q = quals[i] as usize;
                        qual_arr[i][q.min(63)] += 1;

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
                    for i in 0..read_len {
                        let ro_cycle = read_len - 1 - i;
                        let q = quals[i] as usize;
                        qual_arr[ro_cycle][q.min(63)] += 1;

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

                // Save gc_count for GCD section below (avoids re-scanning the sequence).
                primary_gc_count = gc_count;

                // GC content: cumulative step function with ngc=200 bins.
                // Matches samtools stats.c:925-941. For a read with gc_count G/C
                // bases out of read_len total, increment bins gc_idx_min..gc_idx_max.
                let ngc: usize = 200;
                if let (Some(gc_idx_min), Some(gc_idx_max)) = (
                    (gc_count as usize * (ngc - 1)).checked_div(read_len),
                    ((gc_count as usize + 1) * (ngc - 1)).checked_div(read_len),
                ) {
                    let gc_idx_max = gc_idx_max.min(ngc - 1);
                    for item in gc_arr.iter_mut().take(gc_idx_max).skip(gc_idx_min) {
                        *item += 1;
                    }
                }
            }
        } // if is_primary

        // =============================================================
        // Indel distribution (ID) and indels per cycle (IC) from CIGAR.
        //
        // Upstream samtools stats calls count_indels() AFTER the
        // secondary-read early return (line 1206-1210) and the
        // IS_UNMAPPED return (line 1255), but OUTSIDE IS_ORIGINAL().
        // This means: all mapped, non-secondary reads are included
        // (supplementary, duplicate, qcfail all contribute).
        //
        // IC uses first-fragment/last-fragment read order (not
        // forward/reverse strand) and read-oriented cycle indices,
        // matching upstream count_indels().
        // =============================================================
        // =============================================================
        // Combined single-CIGAR-pass block for IC/ID (indel distribution),
        // bases_mapped_cigar, and COV (coverage ring-buffer pileup).
        //
        // Both IC/ID and COV apply to the same read set (mapped,
        // non-secondary). Merging them into one CIGAR traversal
        // eliminates two redundant record.cigar() calls per read.
        //
        // IC/ID: Upstream samtools stats calls count_indels() outside
        //   IS_ORIGINAL() — supplementary/dup/qcfail all contribute.
        //   IC uses first/last-fragment order and read-oriented cycles.
        //
        // COV: Circular-buffer pileup; buffer flushed up to read start
        //   before CIGAR walk; M/=/X blocks inserted as ranges.
        //   Buffer grown to max_read_len * 5 as needed.
        // =============================================================
        if is_mapped && !is_secondary {
            use rust_htslib::bam::record::Cigar as C;
            let is_reverse = flags & BAM_FREVERSE != 0;
            let read_len = record.seq_len();
            let tid = record.tid();
            let pos = record.pos(); // 0-based

            // Upstream order: paired ? (read1?FIRST:0)+(read2?LAST:0) : FIRST
            let order: u32 = if is_paired {
                (if flags & BAM_FREAD1 != 0 { 1 } else { 0 })
                    + (if flags & BAM_FREAD2 != 0 { 2 } else { 0 })
            } else {
                1 // unpaired → FIRST
            };

            // COV buffer setup (must happen before CIGAR walk).
            // Skip reads with no sequence (upstream samtools early-return).
            let do_cov = read_len > 0;
            let buf_size = if do_cov {
                // Grow buffer to max_read_len * 5 if needed.
                // When growing, linearise the circular data just like
                // upstream samtools: copy [idx..old_size] then [0..idx]
                // into a fresh buffer, and reset idx to 0.
                let need = read_len * 5;
                if need > self.cov_buf.len() {
                    let old_size = self.cov_buf.len();
                    let mut new_buf = vec![0u32; need];
                    let head = old_size - self.cov_buf_idx;
                    new_buf[..head].copy_from_slice(&self.cov_buf[self.cov_buf_idx..]);
                    new_buf[head..head + self.cov_buf_idx]
                        .copy_from_slice(&self.cov_buf[..self.cov_buf_idx]);
                    self.cov_buf = new_buf;
                    self.cov_buf_idx = 0;
                }
                let bs = self.cov_buf.len();
                // Flush entire buffer on chromosome change
                if tid != self.cov_buf_tid {
                    self.flush_cov_buf_all();
                    self.cov_buf_tid = tid;
                    self.cov_buf_pos = pos;
                    self.cov_buf_idx = 0;
                }
                // Flush positions from cov_buf_pos up to read start
                self.cov_buf_flush_to(pos, bs);
                bs
            } else {
                0
            };

            // Single CIGAR traversal serving IC/ID + bases_mapped_cigar + COV
            let cigar = record.cigar();
            let mut icycle: usize = 0;
            let mut cigar_mapped: u64 = 0;
            let mut ref_pos = pos;

            for op in cigar.iter() {
                match op {
                    C::Ins(n) => {
                        let ncig = *n as usize;
                        let len = *n as u64;
                        cigar_mapped += len; // I counts toward bases_mapped_cigar

                        // ID: indel size distribution
                        let id_entry = self.id_hist.entry(len).or_insert([0; 2]);
                        id_entry[0] += 1; // insertions

                        // IC: indels per cycle (read-oriented index)
                        let idx = if is_reverse {
                            read_len.saturating_sub(icycle + ncig)
                        } else {
                            icycle
                        };
                        if idx >= self.ic.len() {
                            self.ic.resize(idx + 1, [0u64; 4]);
                        }
                        if order == 1 {
                            self.ic[idx][0] += 1; // ins_1st
                        }
                        if order == 2 {
                            self.ic[idx][1] += 1; // ins_2nd
                        }

                        icycle += ncig; // I advances query cycle; ref unchanged
                                        // COV: I consumes no reference positions
                    }
                    C::Del(n) => {
                        let len = *n as u64;
                        // ID: indel size distribution
                        let id_entry = self.id_hist.entry(len).or_insert([0; 2]);
                        id_entry[1] += 1; // deletions

                        // IC: indels per cycle (read-oriented index)
                        let idx = if is_reverse {
                            if icycle == 0 {
                                // Discard meaningless deletions at cycle 0
                                // (upstream: "if (idx<0) continue;")
                                ref_pos += *n as i64; // still advance ref for COV
                                continue;
                            }
                            read_len.saturating_sub(icycle + 1)
                        } else {
                            if icycle == 0 {
                                ref_pos += *n as i64;
                                continue;
                            }
                            icycle - 1
                        };
                        if idx >= self.ic.len() {
                            self.ic.resize(idx + 1, [0u64; 4]);
                        }
                        if order == 1 {
                            self.ic[idx][2] += 1; // del_1st
                        }
                        if order == 2 {
                            self.ic[idx][3] += 1; // del_2nd
                        }
                        // D does NOT advance query cycle; does advance ref
                        ref_pos += *n as i64;
                    }
                    C::Match(n) | C::Equal(n) | C::Diff(n) => {
                        let len = *n as u64;
                        cigar_mapped += len; // M/=/X count toward bases_mapped_cigar
                        icycle += *n as usize;
                        // COV: M/=/X consumes reference positions
                        if do_cov {
                            let end = ref_pos + *n as i64;
                            self.cov_buf_insert(ref_pos, end, buf_size);
                            ref_pos = end;
                        } else {
                            ref_pos += *n as i64;
                        }
                    }
                    C::RefSkip(n) => {
                        ref_pos += *n as i64; // N advances ref (COV skips it)
                    }
                    C::SoftClip(n) => {
                        icycle += *n as usize; // S advances query cycle
                                               // COV: S consumes no reference positions
                    }
                    C::HardClip(_) | C::Pad(_) => {}
                }
            }
            self.bases_mapped_cigar += cigar_mapped;
        } // if is_mapped && !is_secondary (IC/ID + COV combined)

        // =============================================================
        // GCD: GC-depth accumulation (no-reference path).
        //
        // Matches upstream samtools stats without --ref-seq: bins of
        // GCD_BIN_SIZE bp, depth incremented for each read, GC fraction
        // accumulated from the read's sequence.
        //
        // Included reads: mapped, non-secondary (same as COV).
        //
        // NOTE: gc_count_for_gcd is set from the primary-read per-cycle
        // loop above (when is_primary is true), or computed here only for
        // non-primary mapped reads, avoiding a redundant full sequence scan.
        // =============================================================
        if is_mapped && !is_secondary {
            let tid = record.tid();
            let pos = record.pos();
            let seq_len = record.seq_len();

            if seq_len > 0 {
                // Start a new bin on: first read, chromosome change, or
                // read beyond current bin boundary.
                let new_bin = self.gcd_pos < 0
                    || tid != self.gcd_tid
                    || pos - self.gcd_pos > GCD_BIN_SIZE as i64;

                if new_bin {
                    self.gcd_bins.push(GcDepthBin { gc: 0.0, depth: 0 });
                    self.gcd_pos = pos;
                    self.gcd_tid = tid;
                }

                // Increment depth and accumulate GC fraction from read seq.
                if let Some(bin) = self.gcd_bins.last_mut() {
                    bin.depth += 1;
                    // For primary reads, gc_count was already computed in the
                    // per-cycle base loop above. For non-primary mapped reads
                    // (supplementary, etc.) compute it here from the sequence.
                    let gc_count: u32 = if is_primary {
                        primary_gc_count as u32
                    } else {
                        let seq = record.seq();
                        let mut count: u32 = 0;
                        for i in 0..seq_len {
                            let base = seq.encoded_base(i);
                            if base == 2 || base == 4 {
                                count += 1;
                            }
                        }
                        count
                    };
                    bin.gc += gc_count as f32 / seq_len as f32;
                }
            }
        } // if is_mapped && !is_secondary (GCD)

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

    /// Flush all remaining positions in the coverage round buffer into cov_hist.
    /// Must be called after processing all reads (or when switching chromosomes).
    /// Flush the circular buffer from `cov_buf_pos` up to (but not including) `pos`.
    /// Each slot's depth is recorded in `cov_hist` and the slot is zeroed.
    /// Matches upstream `round_buffer_flush` logic from samtools stats.c.
    fn cov_buf_flush_to(&mut self, pos: i64, buf_size: usize) {
        if pos - self.cov_buf_pos >= buf_size as i64 {
            // Gap exceeds buffer size.  Match upstream samtools exactly:
            // flush `size - 1` positions (from cov_buf_pos to
            // cov_buf_pos + size - 2), leaving the LAST slot untouched.
            // Then advance idx by `size - 1` and jump pos.
            //
            // Upstream (stats.c round_buffer_flush lines 334-366):
            //   pos = rbuf.pos + size - 1;          // cap at last slot
            //   ito = lidx2ridx(start, size, rbuf.pos, pos-1);
            //   // flush from start to ito (size-1 slots)
            //   rbuf.start = lidx2ridx(start, size, rbuf.pos, pos);
            //   rbuf.pos = new_pos;
            let flush_count = buf_size - 1; // flush all but the last slot
            for _ in 0..flush_count {
                let depth = self.cov_buf[self.cov_buf_idx];
                if depth > 0 {
                    *self.cov_hist.entry(depth).or_insert(0) += 1;
                    self.cov_buf[self.cov_buf_idx] = 0;
                }
                self.cov_buf_idx += 1;
                if self.cov_buf_idx >= buf_size {
                    self.cov_buf_idx = 0;
                }
            }
            // idx now points to the ONE unflushed slot (the last position
            // in the old window).  Jump pos to the new read position.
            self.cov_buf_pos = pos;
        } else {
            // Normal case: flush slot by slot.
            while self.cov_buf_pos < pos {
                let depth = self.cov_buf[self.cov_buf_idx];
                if depth > 0 {
                    *self.cov_hist.entry(depth).or_insert(0) += 1;
                    self.cov_buf[self.cov_buf_idx] = 0;
                }
                self.cov_buf_idx += 1;
                if self.cov_buf_idx >= buf_size {
                    self.cov_buf_idx = 0;
                }
                self.cov_buf_pos += 1;
            }
        }
    }

    /// Insert a contiguous reference range `[from, to)` into the circular buffer,
    /// incrementing depth for each position. The range must fit within `buf_size`.
    fn cov_buf_insert(&mut self, from: i64, to: i64, buf_size: usize) {
        for ref_pos in from..to {
            // Map ref_pos to buffer index: offset from cov_buf_idx by (ref_pos - cov_buf_pos)
            let offset = (ref_pos - self.cov_buf_pos) as usize;
            let idx = (self.cov_buf_idx + offset) % buf_size;
            self.cov_buf[idx] += 1;
        }
    }

    /// Flush the entire circular buffer and reset tracking state.
    pub fn flush_cov_buf_all(&mut self) {
        for slot in self.cov_buf.iter_mut() {
            if *slot > 0 {
                *self.cov_hist.entry(*slot).or_insert(0) += 1;
                *slot = 0;
            }
        }
        self.cov_buf_idx = 0;
        self.cov_buf_pos = 0;
        self.cov_buf_tid = -1;
    }

    /// Merge another accumulator into this one.
    pub fn merge(&mut self, mut other: BamStatAccum) {
        // Flush any remaining positions in the other's round buffer into its
        // cov_hist before merging.  Without this, positions still in the
        // round buffer would be silently lost during parallel merges.
        other.flush_cov_buf_all();

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
        self.mismatches += other.mismatches;
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

        // GCD bins (concatenate — bins from different chromosome workers
        // are independent and will be sorted during output).
        self.gcd_bins.append(&mut other.gcd_bins);
    }
}

impl BamStatAccum {
    /// Convert accumulated counters into a `BamStatResult` for output.
    pub fn into_result(mut self) -> BamStatResult {
        // Flush remaining positions in the coverage round buffer
        self.flush_cov_buf_all();
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
            gcd_bins: self.gcd_bins,
        }
    }
}
