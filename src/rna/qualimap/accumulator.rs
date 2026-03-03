//! Qualimap-compatible RNA-Seq QC accumulator.
//!
//! Processes BAM records using Qualimap's counting rules:
//! - M-only CIGAR blocks (no insertion/deletion/softclip contribution)
//! - Enclosure-based gene assignment (each M-block must be fully contained in exons of the same gene)
//! - PE reads: combine both mates' M-blocks, assign on combined set, count numReads=2
//! - Multi-mappers excluded entirely (NH > 1)
//! - Junction motif extraction from read sequence at N-operation boundaries

use std::collections::{HashMap, HashSet};

use coitrees::IntervalTree;
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

use super::coverage::{MergedGeneCoverage, TranscriptCoverage};
use super::index::QualimapIndex;

// ============================================================
// BAM flag constants
// ============================================================

const BAM_FUNMAP: u16 = 0x4;
const BAM_FSECONDARY: u16 = 0x100;
const BAM_FQCFAIL: u16 = 0x200;
const BAM_FSUPPLEMENTARY: u16 = 0x800;
const BAM_FPAIRED: u16 = 0x1;
const BAM_FPROPER_PAIR: u16 = 0x2;
const BAM_FREAD1: u16 = 0x40;
const BAM_FREAD2: u16 = 0x80;
const BAM_FREVERSE: u16 = 0x10;

// ============================================================
// Mate buffer for PE reads
// ============================================================

/// Key for mate buffering: (qname_hash, r1_tid, r1_pos, r2_tid, r2_pos).
type MateKey = (u64, i32, i64, i32, i64);

/// Buffered mate information for PE reconciliation.
#[derive(Debug, Clone)]
struct MateInfo {
    /// M-only aligned blocks (0-based half-open).
    m_blocks: Vec<(i32, i32)>,
    /// Genes whose exons enclose ALL M-blocks of this mate.
    enclosing_genes: HashSet<u32>,
    /// Read sequence bytes at junction boundaries (for motif extraction).
    junction_motifs: Vec<JunctionMotif>,
    /// Whether this mate is on the reverse strand (flag 0x10).
    is_reverse: bool,
    /// Whether this mate is the first read in pair (flag 0x40).
    is_first_of_pair: bool,
}

/// A splice junction motif extracted from the read sequence at an N-operation.
#[derive(Debug, Clone)]
#[allow(dead_code)]
struct JunctionMotif {
    /// Genomic start of the intron (0-based).
    intron_start: i32,
    /// Genomic end of the intron (0-based).
    intron_end: i32,
    /// 2bp donor motif (from read, e.g., "GT").
    donor: [u8; 2],
    /// 2bp acceptor motif (from read, e.g., "AG").
    acceptor: [u8; 2],
}

// ============================================================
// Qualimap counters
// ============================================================

/// Accumulated Qualimap-compatible QC counters.
///
/// These match the counters in Qualimap's `RNASeqQCAnalysis` and are output
/// in `rnaseq_qc_results.txt`.
#[derive(Debug, Clone, Default)]
pub struct QualimapCounters {
    /// Total primary alignments processed.
    pub primary_alignments: u64,
    /// Total secondary alignments seen (skipped).
    pub secondary_alignments: u64,
    /// Unmapped reads (flag 0x4).
    pub not_aligned: u64,
    /// Multi-mapped reads (NH > 1), excluded from counting.
    pub alignment_not_unique: u64,
    /// Reads assigned to exactly one gene (exonic).
    pub exonic_reads: u64,
    /// Reads not assigned to any gene and not intronic (intergenic).
    pub no_feature: u64,
    /// Reads assigned to more than one gene.
    pub ambiguous: u64,
    /// Reads with any M-block overlapping an intron.
    pub intronic_reads: u64,
    /// Reads not overlapping any exon or intron.
    pub intergenic_reads: u64,
    /// No-feature reads that overlap an exon region ("intronic/intergenic
    /// overlapping exon" in Qualimap terminology). Only counted once per
    /// read/fragment that has any M-block overlapping an exon.
    pub overlapping_exon: u64,
    /// Total reads counted (numReads: for PE, 2 per fragment).
    pub read_count: u64,
    /// Total fragments counted (PE: 1 per pair, SE: 1 per read).
    pub fragment_count: u64,
    /// Left mates in proper pairs (flag 0x2 + 0x40).
    pub left_proper_in_pair: u64,
    /// Right mates in proper pairs (flag 0x2 + 0x80).
    pub right_proper_in_pair: u64,
    /// Both mates with proper-pair flag set (for numberOfMappedPairs = count / 2).
    pub both_proper_in_pair: u64,
    /// Reads at splice junctions (with N-operations in CIGAR).
    pub reads_at_junctions: u64,
    /// Total supplementary alignments seen (skipped).
    pub supplementary: u64,
    /// SSP forward-strand estimation counter.
    pub ssp_fwd: u64,
    /// SSP reverse-strand estimation counter.
    pub ssp_rev: u64,
}

impl QualimapCounters {
    /// Merge counters from another accumulator.
    pub fn merge(&mut self, other: &QualimapCounters) {
        self.primary_alignments += other.primary_alignments;
        self.secondary_alignments += other.secondary_alignments;
        self.not_aligned += other.not_aligned;
        self.alignment_not_unique += other.alignment_not_unique;
        self.exonic_reads += other.exonic_reads;
        self.no_feature += other.no_feature;
        self.ambiguous += other.ambiguous;
        self.intronic_reads += other.intronic_reads;
        self.intergenic_reads += other.intergenic_reads;
        self.overlapping_exon += other.overlapping_exon;
        self.read_count += other.read_count;
        self.fragment_count += other.fragment_count;
        self.left_proper_in_pair += other.left_proper_in_pair;
        self.right_proper_in_pair += other.right_proper_in_pair;
        self.both_proper_in_pair += other.both_proper_in_pair;
        self.reads_at_junctions += other.reads_at_junctions;
        self.supplementary += other.supplementary;
        self.ssp_fwd += other.ssp_fwd;
        self.ssp_rev += other.ssp_rev;
    }
}

// ============================================================
// Junction motif accumulator
// ============================================================

/// Accumulated junction motif counts (donor-acceptor pairs).
#[derive(Debug, Clone, Default)]
pub struct JunctionMotifCounts {
    /// Counts keyed by (donor_2bp, acceptor_2bp) as uppercase ASCII.
    pub counts: HashMap<([u8; 2], [u8; 2]), u64>,
}

impl JunctionMotifCounts {
    /// Merge counts from another accumulator.
    pub fn merge(&mut self, other: &JunctionMotifCounts) {
        for (key, &count) in &other.counts {
            *self.counts.entry(*key).or_insert(0) += count;
        }
    }
}

// ============================================================
// QualimapAccum — the main accumulator
// ============================================================

/// Qualimap RNA-Seq QC accumulator.
///
/// Processes BAM records using Qualimap's counting rules and accumulates:
/// - Read/fragment classification counters
/// - Per-transcript coverage
/// - Junction motif counts
///
/// This accumulator is designed to integrate into the existing BAM pass via
/// `RseqcAccumulators`, receiving each BAM record through `process_read()`.
#[derive(Debug, Clone)]
pub struct QualimapAccum {
    /// QC counters matching Qualimap's output.
    pub counters: QualimapCounters,
    /// Per-transcript per-base coverage.
    pub coverage: TranscriptCoverage,
    /// Per-gene merged-model coverage.
    pub merged_gene_coverage: MergedGeneCoverage,
    /// Junction splice motif counts.
    pub junction_motifs: JunctionMotifCounts,
    /// Mate buffer for PE reconciliation.
    mate_buffer: HashMap<MateKey, MateInfo>,
    /// Strandedness protocol: 0=unstranded, 1=forward, 2=reverse.
    stranded: u8,
}

impl QualimapAccum {
    /// Create a new Qualimap accumulator.
    ///
    /// # Arguments
    /// * `stranded` - Strandedness: 0=unstranded, 1=forward, 2=reverse.
    pub fn new(stranded: u8) -> Self {
        Self {
            counters: QualimapCounters::default(),
            coverage: TranscriptCoverage::new(),
            merged_gene_coverage: MergedGeneCoverage::new(),
            junction_motifs: JunctionMotifCounts::default(),
            mate_buffer: HashMap::new(),
            stranded,
        }
    }

    /// Process a single BAM record.
    ///
    /// Applies Qualimap's filtering and counting logic:
    /// 1. Skip unmapped, secondary, QC-fail, supplementary reads
    /// 2. Skip multi-mappers (NH > 1)
    /// 3. Extract M-only CIGAR blocks
    /// 4. For SE: assign immediately via enclosure check
    /// 5. For PE: buffer mate, assign when pair is complete
    ///
    /// # Arguments
    /// * `record` - The BAM record to process.
    /// * `chrom` - Chromosome name for this record.
    /// * `index` - The Qualimap annotation index.
    pub fn process_read(&mut self, record: &bam::Record, chrom: &str, index: &QualimapIndex) {
        let flags = record.flags();
        log::trace!(
            "QM process_read chrom={} flags={} pos={}",
            chrom,
            flags,
            record.pos()
        );

        // Skip unmapped
        if flags & BAM_FUNMAP != 0 {
            self.counters.not_aligned += 1;
            return;
        }

        // Multi-mapper check on ALL mapped reads (including secondary).
        // Qualimap increments alignmentNotUnique before the skipSecondary check,
        // so secondary reads with NH>1 are counted in both secondary AND non-unique.
        let is_multi_mapped = get_nh_tag(record).is_some_and(|nh| nh > 1);
        if is_multi_mapped {
            self.counters.alignment_not_unique += 1;
        }

        // Track secondary vs primary
        let is_secondary = flags & BAM_FSECONDARY != 0;
        if is_secondary {
            self.counters.secondary_alignments += 1;
            // Qualimap skipSecondary=true by default in RNA-Seq mode.
            // Secondary reads are counted but excluded from all downstream
            // gene assignment, SSP, and junction analysis.
            return;
        }
        self.counters.primary_alignments += 1;

        // Count left/right reads for ALL primary paired reads (not just proper pair).
        // Qualimap: leftProperInPair counts first-of-pair, rightProperInPair counts
        // second-of-pair — despite the name, they don't check the proper-pair flag.
        // bothProperInPair counts reads that ARE proper pair (used for numberOfMappedPairs).
        if flags & BAM_FPAIRED != 0 {
            if flags & BAM_FREAD1 != 0 {
                self.counters.left_proper_in_pair += 1;
            }
            if flags & BAM_FREAD2 != 0 {
                self.counters.right_proper_in_pair += 1;
            }
            if flags & BAM_FPROPER_PAIR != 0 {
                self.counters.both_proper_in_pair += 1;
            }
        }

        // Skip supplementary and QC-fail (these never contribute)
        if flags & BAM_FSUPPLEMENTARY != 0 {
            self.counters.supplementary += 1;
            return;
        }
        if flags & BAM_FQCFAIL != 0 {
            return;
        }

        // readCount: primary reads that pass initial filters
        // (Qualimap: readCount counts all non-secondary, non-supplementary reads)
        self.counters.read_count += 1;

        // Extract M-only CIGAR blocks (replicating Qualimap's CIGAR bug where
        // offset += length for ALL operations, including I/S/H/P).
        let m_blocks = extract_m_blocks(record);

        // Skip multi-mappers from gene assignment, junction counting, and coverage.
        // Coverage is added only for uniquely-mapped reads, restricted to fully-enclosed
        // M-blocks, matching the correct Qualimap coverage model.
        if is_multi_mapped {
            return;
        }
        if m_blocks.is_empty() {
            return;
        }

        // Extract junction motifs from N-operations in CIGAR.
        // Only uniquely-mapped primary reads contribute junctions.
        let (n_op_count, motifs) = extract_junction_motifs(record);
        // Qualimap increments numReadsWithJunction for every N-operation,
        // even when the motif extraction guard fails. Match that behavior.
        self.counters.reads_at_junctions += n_op_count as u64;

        // Determine enclosing genes for this read's M-blocks (with strand filtering)
        let is_reverse = flags & BAM_FREVERSE != 0;
        let is_first_of_pair = flags & BAM_FREAD1 != 0;
        let enclosing_genes = find_enclosing_genes(
            &m_blocks,
            chrom,
            index,
            self.stranded,
            is_reverse,
            is_first_of_pair,
        );

        // --- Coverage: per-block, strand-independent (Qualimap model) ---
        // Qualimap adds coverage per M-block to each enclosing gene, without
        // strand filtering. This is separate from gene assignment counting.
        // For SE reads: add coverage now. For PE: deferred to reconcile_pair().
        if flags & BAM_FPAIRED == 0 || flags & BAM_FPROPER_PAIR == 0 {
            self.add_per_block_coverage(&m_blocks, chrom, index);
        }

        // PE path: check BAM flags directly (no library-level paired flag needed)
        if flags & BAM_FPAIRED != 0 && flags & BAM_FPROPER_PAIR != 0 {
            // PE proper pair: buffer and combine with mate (Qualimap name-sorts
            // and processes both mates' intervals together as a fragment).
            self.process_pe_read(record, chrom, m_blocks, enclosing_genes, motifs, index);
        } else {
            // SE mode or non-proper pair: assign immediately (numReads = 1)
            let is_reverse = flags & BAM_FREVERSE != 0;
            let is_first = flags & BAM_FREAD1 != 0;
            self.assign_se_read(
                &m_blocks,
                enclosing_genes,
                motifs,
                chrom,
                index,
                is_reverse,
                is_first,
            );
        }
    }

    /// Process a PE read — buffer and reconcile with mate.
    fn process_pe_read(
        &mut self,
        record: &bam::Record,
        chrom: &str,
        m_blocks: Vec<(i32, i32)>,
        enclosing_genes: HashSet<u32>,
        motifs: Vec<JunctionMotif>,
        index: &QualimapIndex,
    ) {
        let key = make_mate_key(record);

        let flags = record.flags();
        let is_reverse = flags & BAM_FREVERSE != 0;
        let is_first_of_pair = flags & BAM_FREAD1 != 0;

        if let Some(mate) = self.mate_buffer.remove(&key) {
            // Found the mate — reconcile the pair.
            // Current read is r1 (second-seen), buffered mate is r2 (first-seen).
            self.reconcile_pair(
                &m_blocks,
                &motifs,
                is_reverse,
                is_first_of_pair,
                &mate.m_blocks,
                &mate.junction_motifs,
                mate.is_reverse,
                mate.is_first_of_pair,
                chrom,
                index,
            );
        } else {
            // First mate seen — buffer it
            self.mate_buffer.insert(
                key,
                MateInfo {
                    m_blocks,
                    enclosing_genes,
                    junction_motifs: motifs,
                    is_reverse,
                    is_first_of_pair,
                },
            );
        }
    }

    /// Reconcile a PE pair: combine M-blocks and gene assignments.
    ///
    /// Qualimap intersects the enclosing gene sets of both mates.
    /// If exactly one gene encloses ALL M-blocks of BOTH mates → assigned.
    /// If >1 gene → ambiguous. If 0 genes → check intronic/intergenic.
    #[allow(clippy::too_many_arguments)]
    fn reconcile_pair(
        &mut self,
        r1_blocks: &[(i32, i32)],
        r1_motifs: &[JunctionMotif],
        r1_is_reverse: bool,
        r1_is_first_of_pair: bool,
        r2_blocks: &[(i32, i32)],
        r2_motifs: &[JunctionMotif],
        r2_is_reverse: bool,
        _r2_is_first_of_pair: bool,
        chrom: &str,
        index: &QualimapIndex,
    ) {
        // Combine M-blocks from both mates (used for coverage and overlapping-exon checks)
        let mut combined_blocks: Vec<(i32, i32)> =
            Vec::with_capacity(r1_blocks.len() + r2_blocks.len());
        combined_blocks.extend_from_slice(r1_blocks);
        combined_blocks.extend_from_slice(r2_blocks);

        // Qualimap's findIntersectingFeatures processes each read's intervals with its own
        // corrected strand, then a gene must cover intervals from BOTH mates. We replicate
        // this by running find_enclosing_genes per mate and intersecting the two gene sets.
        let r2_is_first_of_pair = !r1_is_first_of_pair;
        let r1_genes = find_enclosing_genes(
            r1_blocks,
            chrom,
            index,
            self.stranded,
            r1_is_reverse,
            r1_is_first_of_pair,
        );
        let r2_genes = find_enclosing_genes(
            r2_blocks,
            chrom,
            index,
            self.stranded,
            r2_is_reverse,
            r2_is_first_of_pair,
        );
        let assigned_genes: HashSet<u32> = r1_genes.intersection(&r2_genes).copied().collect();

        // --- Coverage: per-block, strand-independent (Qualimap model) ---
        // Qualimap adds coverage per M-block to each enclosing gene, without
        // strand filtering. Uses combined blocks from both mates.
        self.add_per_block_coverage(&combined_blocks, chrom, index);

        // Collect junction motifs from both mates
        let all_motifs: Vec<&JunctionMotif> = r1_motifs.iter().chain(r2_motifs.iter()).collect();
        for motif in &all_motifs {
            *self
                .junction_motifs
                .counts
                .entry((motif.donor, motif.acceptor))
                .or_insert(0) += 1;
        }

        // Fragment assignment
        self.counters.fragment_count += 1;

        // Check for overlapping_exon (runs for ALL reads, not just no-feature).
        // Qualimap increments this in findIntersectingFeatures() for any read
        // where an M-block overlaps an exon without being fully enclosed by it.
        self.check_overlapping_exon(&combined_blocks, chrom, index);

        match assigned_genes.len() {
            0 => {
                // No gene encloses all blocks — classify as intronic or intergenic
                self.classify_no_feature(&combined_blocks, chrom, index, 2);
                // Qualimap: noFeature += numReads (2 for PE fragment)
                self.counters.no_feature += 2;
            }
            1 => {
                // Exactly one gene — assigned (exonic)
                // Qualimap: exonicReads += numReads (2 for PE fragment)
                self.counters.exonic_reads += 2;
            }
            _ => {
                // Multiple genes — ambiguous
                // Qualimap: ambiguous += numReads (2 for PE fragment)
                self.counters.ambiguous += 2;
            }
        }

        // SSP estimation: counted per (M-block, gene) pair for ALL fragments —
        // including ambiguous — because Qualimap counts SSP in findIntersectingFeatures()
        // before the exonic/ambiguous/no-feature decision. Count each mate separately
        // with its own strand and first-of-pair flags.
        self.count_ssp_for_blocks(r1_blocks, chrom, index, r1_is_reverse, r1_is_first_of_pair);
        let r2_is_first = !r1_is_first_of_pair;
        self.count_ssp_for_blocks(r2_blocks, chrom, index, r2_is_reverse, r2_is_first);
    }

    /// Assign an SE read immediately.
    #[allow(clippy::too_many_arguments)]
    fn assign_se_read(
        &mut self,
        m_blocks: &[(i32, i32)],
        enclosing_genes: HashSet<u32>,
        motifs: Vec<JunctionMotif>,
        chrom: &str,
        index: &QualimapIndex,
        is_reverse: bool,
        is_first_of_pair: bool,
    ) {
        // Record junction motifs
        for motif in &motifs {
            *self
                .junction_motifs
                .counts
                .entry((motif.donor, motif.acceptor))
                .or_insert(0) += 1;
        }

        self.counters.fragment_count += 1;

        // Check for overlapping_exon (runs for ALL reads, not just no-feature).
        // Qualimap increments this in findIntersectingFeatures() for any read
        // where an M-block overlaps an exon without being fully enclosed by it.
        self.check_overlapping_exon(m_blocks, chrom, index);

        match enclosing_genes.len() {
            0 => {
                self.classify_no_feature(m_blocks, chrom, index, 1);
                self.counters.no_feature += 1;
            }
            1 => {
                self.counters.exonic_reads += 1;
            }
            _ => {
                self.counters.ambiguous += 1;
            }
        }

        // SSP estimation: counted per (M-block, gene) pair for ALL reads —
        // including ambiguous — because Qualimap counts SSP in findIntersectingFeatures()
        // before the exonic/ambiguous/no-feature decision.
        self.count_ssp_for_blocks(m_blocks, chrom, index, is_reverse, is_first_of_pair);
    }

    /// Check if any M-block overlaps an exon region without being fully enclosed.
    ///
    /// Qualimap Java increments `overlapping_exon` in `findIntersectingFeatures()`
    /// for ANY read (exonic, ambiguous, or no-feature) where at least one M-block
    /// overlaps an exon but is NOT fully enclosed by that exon. The `readOverlaps`
    /// flag is set once per read/fragment.
    fn check_overlapping_exon(
        &mut self,
        m_blocks: &[(i32, i32)],
        chrom: &str,
        index: &QualimapIndex,
    ) -> bool {
        if let Some(merged_tree) = index.merged_exon_tree(chrom) {
            let mut overlaps_not_enclosed = false;
            for &(start, end) in m_blocks {
                // COITree uses closed [first, last] intervals internally, but our
                // coordinates are 0-based half-open. The tree may return nodes
                // that are merely abutting (not truly overlapping) in half-open
                // semantics. We verify actual overlap explicitly.
                merged_tree.query(start, end, |iv| {
                    // Verify overlap in half-open coordinates
                    let exon_start = iv.metadata.exon_start;
                    let exon_end = iv.metadata.exon_end;
                    if exon_start >= end || start >= exon_end {
                        return; // Abutting or non-overlapping — skip
                    }
                    // Qualimap checks per-node: if this specific exon does NOT
                    // enclose the M-block, it sets readOverlaps = true. Even if
                    // a different exon DOES enclose the block, the non-enclosing
                    // overlap still counts.
                    let enclosed = exon_start <= start && end <= exon_end;
                    if !enclosed {
                        overlaps_not_enclosed = true;
                    }
                });
                if overlaps_not_enclosed {
                    break;
                }
            }
            if overlaps_not_enclosed {
                // Qualimap counts overlapping_exon once per read/fragment.
                self.counters.overlapping_exon += 1;
                return true;
            }
        }
        false
    }

    /// Add per-transcript and merged-gene coverage, per-block, strand-independent.
    ///
    /// Matches Qualimap's coverage model exactly: in `findIntersectingFeatures`,
    /// `addCoverage` is called inside the per-interval loop for every feature
    /// whose merged-exon node encloses that individual M-block. Coverage is:
    /// - **Per-block**: each M-block is independently checked for enclosure
    /// - **Strand-independent**: no strand protocol check for coverage
    /// - **All enclosing genes**: every gene enclosing a block gets coverage
    ///
    /// # Arguments
    /// * `m_blocks` - M-only aligned blocks (0-based half-open).
    /// * `chrom`    - Chromosome name for tree lookup.
    /// * `index`    - Qualimap annotation index.
    fn add_per_block_coverage(
        &mut self,
        m_blocks: &[(i32, i32)],
        chrom: &str,
        index: &QualimapIndex,
    ) {
        if m_blocks.is_empty() {
            return;
        }
        let tree = match index.merged_exon_tree(chrom) {
            Some(t) => t,
            None => return,
        };

        for &(block_start, block_end) in m_blocks {
            let mut block_genes: HashSet<u32> = HashSet::new();
            tree.query(block_start, block_end, |iv| {
                if block_start >= iv.metadata.exon_start && block_end <= iv.metadata.exon_end {
                    block_genes.insert(iv.metadata.gene_idx);
                }
            });

            let single_block = &[(block_start, block_end)];
            for &gidx in &block_genes {
                let (tx_start, _) = index.gene_transcript_ranges[gidx as usize];
                let transcripts = index.gene_transcripts(gidx);
                self.coverage
                    .add_coverage(single_block, transcripts, tx_start);
                if let Some(model) = index.merged_gene_models.get(gidx as usize) {
                    self.merged_gene_coverage
                        .add_coverage(gidx, single_block, model);
                }
            }
        }
    }

    /// Classify a read that wasn't assigned to any gene as intronic or intergenic.
    fn classify_no_feature(
        &mut self,
        m_blocks: &[(i32, i32)],
        chrom: &str,
        index: &QualimapIndex,
        num_reads: u64,
    ) {
        // Classify as intronic or intergenic
        if let Some(intron_tree) = index.intron_tree(chrom) {
            for &(start, end) in m_blocks {
                let mut found_intron = false;
                intron_tree.query(start, end, |_iv| {
                    found_intron = true;
                });
                if found_intron {
                    self.counters.intronic_reads += num_reads;
                    return;
                }
            }
        }
        self.counters.intergenic_reads += num_reads;
    }

    /// Count SSP votes for a set of M-blocks across all their enclosing genes.
    ///
    /// Qualimap's `findIntersectingFeatures()` counts SSP per (M-block, gene) pair —
    /// i.e., once per gene per block, for every gene whose merged exon encloses that
    /// block. This happens before the exonic/ambiguous/no-feature decision, so
    /// ambiguous reads (and even no-feature reads that happen to be enclosed) still
    /// contribute to the SSP estimation.
    ///
    /// `is_first_of_pair` maps to Qualimap's interval name 'f'/'r':
    ///   - 'f' (first): same strand → fwd++, else rev++
    ///   - 'r' (second): same strand → rev++, else fwd++
    fn count_ssp_for_blocks(
        &mut self,
        m_blocks: &[(i32, i32)],
        chrom: &str,
        index: &QualimapIndex,
        is_reverse: bool,
        is_first_of_pair: bool,
    ) {
        let tree = match index.merged_exon_tree(chrom) {
            Some(t) => t,
            None => return,
        };

        for &(block_start, block_end) in m_blocks {
            // Collect the set of genes that enclose this block, along with
            // their strand. A block may be enclosed by multiple genes simultaneously
            // (Qualimap counts SSP once per gene per block, hence the HashSet of
            // (gene_idx, strand) pairs to avoid double-counting from multiple nodes
            // of the same gene).
            let mut gene_strands: HashSet<(u32, char)> = HashSet::new();
            tree.query(block_start, block_end, |iv| {
                if block_start >= iv.metadata.exon_start && block_end <= iv.metadata.exon_end {
                    let gidx = iv.metadata.gene_idx;
                    let (tx_start, _) = index.gene_transcript_ranges[gidx as usize];
                    let strand = index.transcripts[tx_start as usize].strand;
                    gene_strands.insert((gidx, strand));
                }
            });

            for (_gidx, gene_strand) in &gene_strands {
                let same_strand =
                    (is_reverse && *gene_strand == '-') || (!is_reverse && *gene_strand == '+');
                if is_first_of_pair {
                    if same_strand {
                        self.counters.ssp_fwd += 1;
                    } else {
                        self.counters.ssp_rev += 1;
                    }
                } else if same_strand {
                    self.counters.ssp_rev += 1;
                } else {
                    self.counters.ssp_fwd += 1;
                }
            }
        }
    }

    /// Flush any remaining unpaired mates at the end of a chromosome batch.
    ///
    /// PE mates that were never reconciled are treated as singletons and
    /// Qualimap skips proper-pair fragments where only one mate survived
    /// filters (numSingleReadOnly). We simply discard orphaned mates.
    #[allow(dead_code)]
    pub fn flush_unpaired(&mut self, _index: &QualimapIndex) {
        // Qualimap behaviour: if a proper-pair read's mate was filtered out
        // (by NH>1, unmapped, etc.), the fragment is skipped entirely —
        // neither exonic nor noFeature is counted. We match this by simply
        // draining the mate buffer without processing.
        self.mate_buffer.clear();
    }

    /// Merge another `QualimapAccum` into this one.
    ///
    /// Used when combining results from parallel chromosome batches.
    pub fn merge(&mut self, other: QualimapAccum) {
        self.counters.merge(&other.counters);
        self.coverage.merge(other.coverage);
        self.junction_motifs.merge(&other.junction_motifs);

        // Merge mate buffers (cross-chromosome mates)
        for (key, info) in other.mate_buffer {
            // If the same key exists in both, that means both mates were
            // buffered in different batches — reconcile them
            if let Some(existing) = self.mate_buffer.remove(&key) {
                // We don't have the chrom handy for intronic classification,
                // so treat orphans conservatively
                self.counters.fragment_count += 1;

                let common_genes: HashSet<u32> = existing
                    .enclosing_genes
                    .intersection(&info.enclosing_genes)
                    .copied()
                    .collect();

                let assigned = if existing.enclosing_genes.is_empty()
                    && !info.enclosing_genes.is_empty()
                {
                    info.enclosing_genes.clone()
                } else if info.enclosing_genes.is_empty() && !existing.enclosing_genes.is_empty() {
                    existing.enclosing_genes.clone()
                } else {
                    common_genes
                };

                match assigned.len() {
                    0 => {
                        self.counters.no_feature += 1;
                        self.counters.intergenic_reads += 1;
                    }
                    1 => {
                        self.counters.exonic_reads += 1;
                    }
                    _ => {
                        self.counters.ambiguous += 1;
                    }
                }

                // Junction motifs
                for motif in existing
                    .junction_motifs
                    .iter()
                    .chain(info.junction_motifs.iter())
                {
                    *self
                        .junction_motifs
                        .counts
                        .entry((motif.donor, motif.acceptor))
                        .or_insert(0) += 1;
                }
            } else {
                self.mate_buffer.insert(key, info);
            }
        }
    }

    /// Get the number of unmatched mates still in the buffer.
    #[allow(dead_code)]
    pub fn unmatched_count(&self) -> usize {
        self.mate_buffer.len()
    }

    /// Convert this accumulator into a `QualimapResult`.
    ///
    /// This should be called after all reads have been processed and
    /// `flush_unpaired()` has been called.
    pub fn into_result(self, index: &QualimapIndex) -> super::QualimapResult {
        // Convert junction motif byte arrays to readable strings
        let mut junction_motifs_str = HashMap::new();
        for ((donor, acceptor), count) in &self.junction_motifs.counts {
            let motif = format!(
                "{}{}{}{}",
                donor[0] as char, donor[1] as char, acceptor[0] as char, acceptor[1] as char,
            );
            *junction_motifs_str.entry(motif).or_insert(0u64) += count;
        }

        // Clone the raw coverage before converting to string-keyed map
        let transcript_coverage_raw = self.coverage.clone();

        // Convert per-transcript coverage from flat indices to transcript_id strings
        let mut transcript_coverage = HashMap::new();
        for (flat_idx, depth) in self.coverage.iter() {
            if let Some(tx_info) = index.transcripts.get(flat_idx as usize) {
                transcript_coverage.insert(tx_info.transcript_id.clone(), depth.to_vec());
            }
        }

        super::QualimapResult {
            primary_alignments: self.counters.primary_alignments,
            secondary_alignments: self.counters.secondary_alignments,
            not_aligned: self.counters.not_aligned,
            alignment_not_unique: self.counters.alignment_not_unique,
            exonic_reads: self.counters.exonic_reads,
            ambiguous_reads: self.counters.ambiguous,
            no_feature: self.counters.no_feature,
            intronic_reads: self.counters.intronic_reads,
            intergenic_reads: self.counters.intergenic_reads,
            overlapping_exon_reads: self.counters.overlapping_exon,
            read_count: self.counters.read_count,
            fragment_count: self.counters.fragment_count,
            left_proper_in_pair: self.counters.left_proper_in_pair,
            right_proper_in_pair: self.counters.right_proper_in_pair,
            both_proper_in_pair: self.counters.both_proper_in_pair,
            reads_at_junctions: self.counters.reads_at_junctions,
            junction_motifs: junction_motifs_str,
            transcript_coverage,
            transcript_coverage_raw,
            merged_gene_coverage: self.merged_gene_coverage,
            ssp_fwd: self.counters.ssp_fwd,
            ssp_rev: self.counters.ssp_rev,
        }
    }
}

// ============================================================
// Helper functions
// ============================================================

/// Extract NH tag (number of reported alignments) from a BAM record.
fn get_nh_tag(record: &bam::Record) -> Option<u32> {
    match record.aux(b"NH") {
        Ok(rust_htslib::bam::record::Aux::U8(v)) => Some(v as u32),
        Ok(rust_htslib::bam::record::Aux::U16(v)) => Some(v as u32),
        Ok(rust_htslib::bam::record::Aux::U32(v)) => Some(v),
        Ok(rust_htslib::bam::record::Aux::I8(v)) => Some(v as u32),
        Ok(rust_htslib::bam::record::Aux::I16(v)) => Some(v as u32),
        Ok(rust_htslib::bam::record::Aux::I32(v)) => Some(v as u32),
        _ => None,
    }
}

/// Extract M-only aligned blocks from a BAM record's CIGAR.
///
/// Returns a vector of (start, end) in 0-based half-open coordinates.
/// Only `M` (alignment match) operations are included.
///
/// Extract correct alignment blocks from a BAM record's CIGAR string.
///
/// Unlike [`extract_m_blocks`], this function uses **correct** CIGAR semantics:
/// only M/D/N/=/X operations advance the reference position, while I/S/H/P
/// **Qualimap compatibility note:** Qualimap's Java code advances the reference
/// offset for **all** CIGAR operations, including insertions and soft-clips.
/// This is technically incorrect (I/S consume only query, not reference), but
/// we replicate it here to produce identical enclosure results. The effect is
/// that M-blocks following an I or S operation are shifted rightward by the
/// length of the I/S, which can push them past exon boundaries and change
/// the exonic/intronic classification.
fn extract_m_blocks(record: &bam::Record) -> Vec<(i32, i32)> {
    let mut blocks = Vec::new();
    let mut ref_pos = record.pos() as i32; // 0-based

    for op in record.cigar().iter() {
        match op {
            // Qualimap only creates alignment blocks for M (match/mismatch),
            // not for = (sequence match) or X (sequence mismatch).
            // In Qualimap's getReadIntervals, only CigarOperator.M creates blocks.
            Cigar::Match(len) => {
                let len = *len as i32;
                blocks.push((ref_pos, ref_pos + len));
                ref_pos += len;
            }
            // Qualimap bug: EQ and X advance offset but don't create blocks.
            Cigar::Equal(len) | Cigar::Diff(len) => {
                ref_pos += *len as i32;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                ref_pos += *len as i32;
            }
            // Qualimap bug: insertions and soft-clips advance the reference
            // position even though they only consume query bases.
            Cigar::Ins(len) | Cigar::SoftClip(len) => {
                ref_pos += *len as i32;
            }
            // Qualimap bug: hard clips and pads also advance the reference
            // offset even though they shouldn't consume reference bases.
            Cigar::HardClip(len) | Cigar::Pad(len) => {
                ref_pos += *len as i32;
            }
        }
    }

    blocks
}

/// Extract junction motifs from N-operations in CIGAR.
///
/// For each N (intron/RefSkip) operation, extracts the 2bp donor and 2bp acceptor
/// motifs from the read sequence at the junction boundaries.
/// Returns `(n_op_count, motifs)` where `n_op_count` is the total number of
/// N-operations in the CIGAR (used for `reads_at_junctions`), and `motifs`
/// contains the successfully extracted junction motifs. Qualimap counts
/// `numReadsWithJunction` for every N-operation regardless of whether the
/// motif could be extracted.
fn extract_junction_motifs(record: &bam::Record) -> (usize, Vec<JunctionMotif>) {
    let cigar = record.cigar();
    let seq = record.seq();
    let seq_len = seq.len();

    let mut motifs = Vec::new();
    let mut n_op_count: usize = 0;
    let mut ref_pos = record.pos() as i32;
    let mut seq_pos: usize = 0;

    for op in cigar.iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len = *len as usize;
                ref_pos += len as i32;
                seq_pos += len;
            }
            Cigar::RefSkip(len) => {
                let intron_start = ref_pos;
                let intron_end = ref_pos + *len as i32;
                // Qualimap increments numReadsWithJunction for every N-op,
                // regardless of whether the motif can be extracted.
                n_op_count += 1;

                // Extract donor motif: 2bp at end of preceding exon in read
                // Extract acceptor motif: 2bp at start of next exon in read
                if seq_pos >= 2 && seq_pos + 1 < seq_len {
                    let donor = [seq.encoded_base(seq_pos - 2), seq.encoded_base(seq_pos - 1)];
                    let acceptor = [seq.encoded_base(seq_pos), seq.encoded_base(seq_pos + 1)];

                    // Convert 4-bit encoding to ASCII
                    let donor_ascii = [decode_base(donor[0]), decode_base(donor[1])];
                    let acceptor_ascii = [decode_base(acceptor[0]), decode_base(acceptor[1])];

                    motifs.push(JunctionMotif {
                        intron_start,
                        intron_end,
                        donor: donor_ascii,
                        acceptor: acceptor_ascii,
                    });
                }

                ref_pos = intron_end;
            }
            Cigar::Ins(len) | Cigar::SoftClip(len) => {
                seq_pos += *len as usize;
            }
            Cigar::Del(len) => {
                ref_pos += *len as i32;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    (n_op_count, motifs)
}

/// Decode a 4-bit encoded base to ASCII uppercase.
fn decode_base(encoded: u8) -> u8 {
    match encoded {
        1 => b'A',
        2 => b'C',
        4 => b'G',
        8 => b'T',
        15 => b'N',
        _ => b'N',
    }
}

/// Find genes whose merged exons enclose ALL M-blocks of a read.
///
/// Qualimap Java merges overlapping/abutting exons per gene before building
/// its interval tree. A gene "encloses" a read if EVERY M-block is enclosed
/// by at least one merged exon interval of that gene. This matches Qualimap's
/// `processAlignmentIntervals()` → `findIntersectingFeatures()` → BitSet
/// cardinality check.
///
/// When `stranded > 0`, only genes whose strand matches the read's effective
/// strand are considered. Effective strand follows Qualimap's `getReadIntervals()`
/// logic: for `stranded=2` (reverse), read1 strand is flipped; for `stranded=1`
/// (forward), read2 strand is flipped.
fn find_enclosing_genes(
    m_blocks: &[(i32, i32)],
    chrom: &str,
    index: &QualimapIndex,
    stranded: u8,
    is_reverse: bool,
    is_first_of_pair: bool,
) -> HashSet<u32> {
    let tree = match index.merged_exon_tree(chrom) {
        Some(t) => t,
        None => return HashSet::new(),
    };

    if m_blocks.is_empty() {
        return HashSet::new();
    }

    // Compute the effective strand of this read after protocol correction.
    // Qualimap's getReadIntervals() flips the strand of:
    //   - read1 for strand-specific-reverse (stranded=2)
    //   - read2 for strand-specific-forward (stranded=1)
    // Result: read_on_plus = true means this read maps to the + gene strand.
    let read_on_plus = if stranded == 0 {
        true // no filtering; handled by skipping strand check below
    } else {
        let flip = (stranded == 2 && is_first_of_pair) || (stranded == 1 && !is_first_of_pair);
        if flip {
            is_reverse
        } else {
            !is_reverse
        }
    };

    // Closure: does an exon's strand match this read?
    let strand_ok = |exon_strand: u8| -> bool {
        if stranded == 0 {
            return true;
        }
        (read_on_plus && exon_strand == b'+') || (!read_on_plus && exon_strand == b'-')
    };

    // For the first M-block, find all strand-compatible genes with an enclosing merged exon
    let mut candidate_genes: HashSet<u32> = HashSet::new();
    let (first_start, first_end) = m_blocks[0];
    tree.query(first_start, first_end, |iv| {
        if first_start >= iv.metadata.exon_start
            && first_end <= iv.metadata.exon_end
            && strand_ok(iv.metadata.strand)
        {
            candidate_genes.insert(iv.metadata.gene_idx);
        }
    });

    if candidate_genes.is_empty() {
        return candidate_genes;
    }

    // For subsequent M-blocks, intersect with candidates
    for &(block_start, block_end) in &m_blocks[1..] {
        let mut block_genes: HashSet<u32> = HashSet::new();
        tree.query(block_start, block_end, |iv| {
            if block_start >= iv.metadata.exon_start
                && block_end <= iv.metadata.exon_end
                && candidate_genes.contains(&iv.metadata.gene_idx)
                && strand_ok(iv.metadata.strand)
            {
                block_genes.insert(iv.metadata.gene_idx);
            }
        });
        candidate_genes = block_genes;
        if candidate_genes.is_empty() {
            break;
        }
    }

    candidate_genes
}

/// Build a mate buffer key from a BAM record.
///
/// Uses FNV-1a hash of qname + r1-centric ordering of positions.
fn make_mate_key(record: &bam::Record) -> MateKey {
    let qname_hash = hash_qname(record.qname());
    let tid = record.tid();
    let pos = record.pos();
    let mate_tid = record.mtid();
    let mate_pos = record.mpos();

    // Always store with read1 position first for consistent lookup
    if tid < mate_tid || (tid == mate_tid && pos <= mate_pos) {
        (qname_hash, tid, pos, mate_tid, mate_pos)
    } else {
        (qname_hash, mate_tid, mate_pos, tid, pos)
    }
}

/// FNV-1a hash of a byte slice (for qname hashing).
fn hash_qname(qname: &[u8]) -> u64 {
    let mut hash: u64 = 0xcbf29ce484222325;
    for &byte in qname {
        hash ^= byte as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

// ============================================================
// Unit tests
// ============================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hash_qname_deterministic() {
        let h1 = hash_qname(b"read1");
        let h2 = hash_qname(b"read1");
        let h3 = hash_qname(b"read2");
        assert_eq!(h1, h2);
        assert_ne!(h1, h3);
    }

    #[test]
    fn test_decode_base() {
        assert_eq!(decode_base(1), b'A');
        assert_eq!(decode_base(2), b'C');
        assert_eq!(decode_base(4), b'G');
        assert_eq!(decode_base(8), b'T');
        assert_eq!(decode_base(15), b'N');
        assert_eq!(decode_base(0), b'N');
    }

    #[test]
    fn test_find_enclosing_genes_empty_chrom() {
        let genes = indexmap::IndexMap::new();
        let index = QualimapIndex::from_genes(&genes);
        let result = find_enclosing_genes(&[(100, 200)], "chr1", &index, 0, false, true);
        assert!(result.is_empty());
    }

    #[test]
    fn test_counters_merge() {
        let mut c1 = QualimapCounters::default();
        c1.primary_alignments = 10;
        c1.exonic_reads = 5;
        c1.intronic_reads = 3;

        let c2 = QualimapCounters {
            primary_alignments: 20,
            exonic_reads: 8,
            intronic_reads: 4,
            ..Default::default()
        };

        c1.merge(&c2);
        assert_eq!(c1.primary_alignments, 30);
        assert_eq!(c1.exonic_reads, 13);
        assert_eq!(c1.intronic_reads, 7);
    }

    #[test]
    fn test_junction_motif_counts_merge() {
        let mut j1 = JunctionMotifCounts::default();
        j1.counts.insert(([b'G', b'T'], [b'A', b'G']), 5);

        let mut j2 = JunctionMotifCounts::default();
        j2.counts.insert(([b'G', b'T'], [b'A', b'G']), 3);
        j2.counts.insert(([b'G', b'C'], [b'A', b'G']), 1);

        j1.merge(&j2);
        assert_eq!(j1.counts[&([b'G', b'T'], [b'A', b'G'])], 8);
        assert_eq!(j1.counts[&([b'G', b'C'], [b'A', b'G'])], 1);
    }
}
