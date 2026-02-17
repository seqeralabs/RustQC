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

use super::coverage::TranscriptCoverage;
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
    /// Junction splice motif counts.
    pub junction_motifs: JunctionMotifCounts,
    /// Mate buffer for PE reconciliation.
    mate_buffer: HashMap<MateKey, MateInfo>,
}

impl QualimapAccum {
    /// Create a new Qualimap accumulator.
    pub fn new() -> Self {
        Self {
            counters: QualimapCounters::default(),
            coverage: TranscriptCoverage::new(),
            junction_motifs: JunctionMotifCounts::default(),
            mate_buffer: HashMap::new(),
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

        // Skip multi-mappers from gene assignment (already counted above)
        if is_multi_mapped {
            return;
        }

        // Extract M-only CIGAR blocks
        let m_blocks = extract_m_blocks(record);
        if m_blocks.is_empty() {
            return;
        }

        // Extract junction motifs from N-operations in CIGAR
        let motifs = extract_junction_motifs(record);
        if !motifs.is_empty() {
            // Qualimap counts per-junction (per N operation), not per-read
            self.counters.reads_at_junctions += motifs.len() as u64;
        }

        // Determine enclosing genes for this read's M-blocks
        let enclosing_genes = find_enclosing_genes(&m_blocks, chrom, index);

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
        // Combine M-blocks from both mates
        let mut combined_blocks: Vec<(i32, i32)> =
            Vec::with_capacity(r1_blocks.len() + r2_blocks.len());
        combined_blocks.extend_from_slice(r1_blocks);
        combined_blocks.extend_from_slice(r2_blocks);

        // Re-run enclosure check on all combined blocks (Qualimap's approach:
        // processAlignmentIntervals passes ALL intervals to findIntersectingFeatures)
        let assigned_genes = find_enclosing_genes(&combined_blocks, chrom, index);

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

        match assigned_genes.len() {
            0 => {
                // No gene encloses all blocks — classify as intronic or intergenic
                self.classify_no_feature(&combined_blocks, chrom, index, 2);
                // Qualimap: noFeature += numReads (2 for PE fragment)
                self.counters.no_feature += 2;
            }
            1 => {
                // Exactly one gene — assigned (exonic)
                let gene_idx = *assigned_genes.iter().next().unwrap();
                // Qualimap: exonicReads += numReads (2 for PE fragment)
                self.counters.exonic_reads += 2;

                // SSP estimation: for each enclosed M-block, compare read strand to gene strand.
                // Qualimap counts SSP per M-block, not per read.
                let (tx_start, _tx_end) = index.gene_transcript_ranges[gene_idx as usize];
                let gene_strand = index.transcripts[tx_start as usize].strand;

                // R1 blocks: use r1_is_reverse, r1_is_first_of_pair
                for &(block_start, block_end) in r1_blocks {
                    // Check if this block is enclosed by any exon of the gene
                    if Self::block_enclosed_by_gene(block_start, block_end, gene_idx, chrom, index)
                    {
                        let same_strand = (r1_is_reverse && gene_strand == '-')
                            || (!r1_is_reverse && gene_strand == '+');
                        if r1_is_first_of_pair {
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

                // R2 blocks: mate has opposite first-of-pair, and its own strand
                for &(block_start, block_end) in r2_blocks {
                    if Self::block_enclosed_by_gene(block_start, block_end, gene_idx, chrom, index)
                    {
                        let same_strand = (r2_is_reverse && gene_strand == '-')
                            || (!r2_is_reverse && gene_strand == '+');
                        // r2 is the mate, so first-of-pair is !r1_is_first_of_pair
                        let r2_is_first = !r1_is_first_of_pair;
                        if r2_is_first {
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

                // Add coverage for all transcripts of this gene
                let transcripts = index.gene_transcripts(gene_idx);
                self.coverage
                    .add_coverage(&combined_blocks, transcripts, tx_start);
            }
            _ => {
                // Multiple genes — ambiguous
                // Qualimap: ambiguous += numReads (2 for PE fragment)
                self.counters.ambiguous += 2;
            }
        }
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

        match enclosing_genes.len() {
            0 => {
                self.classify_no_feature(m_blocks, chrom, index, 1);
                self.counters.no_feature += 1;
            }
            1 => {
                let gene_idx = *enclosing_genes.iter().next().unwrap();
                self.counters.exonic_reads += 1;

                // SSP estimation per enclosed M-block: for each M-block that is
                // enclosed by an exon of this gene, compare read strand to gene
                // strand. Qualimap counts SSP per M-block, not per read.
                let (tx_start, _tx_end) = index.gene_transcript_ranges[gene_idx as usize];
                let gene_strand = index.transcripts[tx_start as usize].strand;
                for &(block_start, block_end) in m_blocks {
                    if Self::block_enclosed_by_gene(block_start, block_end, gene_idx, chrom, index)
                    {
                        let same_strand = (is_reverse && gene_strand == '-')
                            || (!is_reverse && gene_strand == '+');
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

                let transcripts = index.gene_transcripts(gene_idx);
                self.coverage.add_coverage(m_blocks, transcripts, tx_start);
            }
            _ => {
                self.counters.ambiguous += 1;
            }
        }
    }

    /// Classify a read that wasn't assigned to any gene as intronic or intergenic.
    ///
    /// Also checks if any M-block overlaps an exon region — if so, increments
    /// `overlapping_exon`. In Qualimap terminology this is "intronic/intergenic
    /// overlapping exon": a read that is NOT enclosed by exons of any gene, but
    /// still physically overlaps an exon region.
    fn classify_no_feature(
        &mut self,
        m_blocks: &[(i32, i32)],
        chrom: &str,
        index: &QualimapIndex,
        num_reads: u64,
    ) {
        // Check if any M-block overlaps an exon region (for overlapping_exon count).
        // Qualimap counts this once per read/fragment via readOverlaps flag.
        if let Some(exon_tree) = index.exon_tree(chrom) {
            let mut overlaps_exon = false;
            for &(start, end) in m_blocks {
                exon_tree.query(start, end, |_iv| {
                    overlaps_exon = true;
                });
                if overlaps_exon {
                    break;
                }
            }
            if overlaps_exon {
                // Qualimap counts overlapping_exon once per read/fragment, not per-read-in-pair.
                self.counters.overlapping_exon += 1;
            }
        }

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

    /// Check if a single M-block is enclosed by any exon of the specified gene.
    ///
    /// Used for per-block SSP estimation where we only count SSP for blocks
    /// that are actually enclosed by an exon of the assigned gene.
    fn block_enclosed_by_gene(
        block_start: i32,
        block_end: i32,
        gene_idx: u32,
        chrom: &str,
        index: &QualimapIndex,
    ) -> bool {
        if let Some(tree) = index.exon_tree(chrom) {
            let mut enclosed = false;
            tree.query(block_start, block_end, |iv| {
                if iv.metadata.gene_idx == gene_idx
                    && block_start >= iv.metadata.exon_start
                    && block_end <= iv.metadata.exon_end
                {
                    enclosed = true;
                }
            });
            enclosed
        } else {
            false
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
        // DEBUG: Print counting diagnostics
        eprintln!("=== QM DEBUG ===");
        eprintln!("  read_count: {}", self.counters.read_count);
        eprintln!("  primary: {}", self.counters.primary_alignments);
        eprintln!("  secondary: {}", self.counters.secondary_alignments);
        eprintln!("  non-unique: {}", self.counters.alignment_not_unique);
        eprintln!("  exonic: {}", self.counters.exonic_reads);
        eprintln!("  ambiguous: {}", self.counters.ambiguous);
        eprintln!("  no_feature: {}", self.counters.no_feature);
        eprintln!("  intronic: {}", self.counters.intronic_reads);
        eprintln!("  intergenic: {}", self.counters.intergenic_reads);
        eprintln!("  overlapping_exon: {}", self.counters.overlapping_exon);
        // Convert junction motif byte arrays to readable strings
        // Convert junction motif byte arrays to readable strings
        let mut junction_motifs_str = HashMap::new();
        for ((donor, acceptor), count) in &self.junction_motifs.counts {
            let motif = format!(
                "{}{}{}{}",
                donor[0] as char, donor[1] as char, acceptor[0] as char, acceptor[1] as char,
            );
            *junction_motifs_str.entry(motif).or_insert(0u64) += count;
        }

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
fn extract_junction_motifs(record: &bam::Record) -> Vec<JunctionMotif> {
    let cigar = record.cigar();
    let seq = record.seq();
    let seq_len = seq.len();

    let mut motifs = Vec::new();
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

    motifs
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

/// Find genes whose exons enclose ALL M-blocks of a read.
///
/// Qualimap's enclosure check: for each M-block, find all exon nodes that
/// fully contain it (block_start >= exon_start AND block_end <= exon_end).
/// A gene "encloses" a read if EVERY M-block is enclosed by at least one exon
/// of that gene.
fn find_enclosing_genes(
    m_blocks: &[(i32, i32)],
    chrom: &str,
    index: &QualimapIndex,
) -> HashSet<u32> {
    let tree = match index.exon_tree(chrom) {
        Some(t) => t,
        None => return HashSet::new(),
    };

    if m_blocks.is_empty() {
        return HashSet::new();
    }

    // For the first M-block, find all genes with an enclosing exon
    let mut candidate_genes: HashSet<u32> = HashSet::new();
    let (first_start, first_end) = m_blocks[0];
    tree.query(first_start, first_end, |iv| {
        if first_start >= iv.metadata.exon_start && first_end <= iv.metadata.exon_end {
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
        let result = find_enclosing_genes(&[(100, 200)], "chr1", &index);
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
