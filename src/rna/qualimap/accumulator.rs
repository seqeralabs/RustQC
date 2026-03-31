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

use crate::cli::Strandedness;

use super::coverage::TranscriptCoverage;
use super::index::QualimapIndex;

// ============================================================
// Cached tree query result for a single M-block
// ============================================================

/// Pre-computed results from a single merged exon tree query for one M-block.
///
/// By querying the tree once per block and caching the results, we avoid
/// redundant COITree queries across `find_enclosing_genes`,
/// `add_per_block_coverage`, `check_overlapping_exon`, and
/// `count_ssp_for_blocks`.
#[derive(Debug, Clone)]
struct CachedBlockHits {
    /// Genes whose merged exon fully encloses this M-block, with their strand.
    /// Deduplicated by `gene_idx`; each gene appears at most once.
    enclosing_genes: Vec<(u32, u8)>,
    /// Whether any node overlaps this M-block without fully enclosing it.
    /// Used for the Qualimap `overlapping_exon` counter.
    has_overlap_not_enclosed: bool,
}

use crate::rna::bam_flags::*;

// ============================================================
// Mate buffer for PE reads
// ============================================================

/// Key for mate buffering: (qname_hash, r1_tid, r1_pos, r2_tid, r2_pos).
type MateKey = (u64, i32, i64, i32, i64);

/// A splice junction motif extracted from the read sequence at an N-operation.
#[derive(Debug, Clone)]
struct JunctionMotif {
    /// 2bp donor motif (from read, e.g., "GT").
    donor: [u8; 2],
    /// 2bp acceptor motif (from read, e.g., "AG").
    acceptor: [u8; 2],
}

/// Per-mate data passed through the PE reconciliation pipeline.
///
/// Groups the per-read data needed for `reconcile_pair()` and
/// `classify_and_count()`, reducing parameter counts.
#[derive(Debug, Clone)]
struct MateData {
    /// M-only aligned blocks (0-based half-open).
    m_blocks: Vec<(i32, i32)>,
    /// Genes whose exons enclose ALL M-blocks of this mate.
    enclosing_genes: HashSet<u32>,
    /// Cached per-block tree query results for this mate.
    cached_hits: Vec<CachedBlockHits>,
    /// Whether this mate is on the reverse strand (flag 0x10).
    is_reverse: bool,
    /// Whether this mate is the first read in pair (flag 0x40).
    is_first_of_pair: bool,
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
    /// Left mates in proper pairs (flag 0x2 + 0x40).
    pub left_proper_in_pair: u64,
    /// Right mates in proper pairs (flag 0x2 + 0x80).
    pub right_proper_in_pair: u64,
    /// Both mates with proper-pair flag set (for numberOfMappedPairs = count / 2).
    pub both_proper_in_pair: u64,
    /// Reads at splice junctions (with N-operations in CIGAR).
    pub reads_at_junctions: u64,
    /// Junction events classified as "known" (both endpoints match known exon boundaries).
    ///
    /// Matches Qualimap Java's `knownJunctions` counter: for each N-operation, if both
    /// the 5' and 3' junction positions overlap with known exon boundary positions
    /// (within ±1 tolerance), the event is "known".
    pub known_junction_events: u64,
    /// Junction events classified as "partly known" (exactly one endpoint matches).
    ///
    /// Matches Qualimap Java's `partlyKnownJunctions` counter.
    pub partly_known_junction_events: u64,
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
        self.left_proper_in_pair += other.left_proper_in_pair;
        self.right_proper_in_pair += other.right_proper_in_pair;
        self.both_proper_in_pair += other.both_proper_in_pair;
        self.reads_at_junctions += other.reads_at_junctions;
        self.known_junction_events += other.known_junction_events;
        self.partly_known_junction_events += other.partly_known_junction_events;
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
    mate_buffer: HashMap<MateKey, MateData>,
    /// Strandedness protocol.
    stranded: Strandedness,
}

impl QualimapAccum {
    /// Create a new Qualimap accumulator.
    pub fn new(stranded: Strandedness) -> Self {
        Self {
            counters: QualimapCounters::default(),
            coverage: TranscriptCoverage::new(),
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

        // Skip multi-mappers from gene assignment, junction counting, and coverage.
        // Coverage is added only for uniquely-mapped reads, restricted to fully-enclosed
        // M-blocks, matching the correct Qualimap coverage model.
        if is_multi_mapped {
            return;
        }

        // Extract M-only CIGAR blocks (replicating Qualimap's CIGAR bug where
        // offset += length for ALL operations, including I/S/H/P) AND junction
        // motifs in a single CIGAR pass, avoiding a second record.cigar() call.
        let (m_blocks, n_op_count, motifs, junction_ref_positions) =
            extract_m_blocks_and_junctions(record);

        if m_blocks.is_empty() {
            return;
        }
        // Qualimap increments numReadsWithJunction for every N-operation,
        // even when the motif extraction guard fails. Match that behavior.
        self.counters.reads_at_junctions += n_op_count as u64;

        // Count junction motifs per-read, matching Qualimap which
        // calls collectJunctionInfo independently for each SAMRecord.
        for motif in &motifs {
            *self
                .junction_motifs
                .counts
                .entry((motif.donor, motif.acceptor))
                .or_insert(0) += 1;
        }

        // Classify junction events as known/partly known/novel using the
        // Qualimap junction location map (exon boundary positions ±1 tolerance).
        for &(pos_ref1, pos_ref2) in &junction_ref_positions {
            let j1 = index.has_junction_overlap(chrom, pos_ref1);
            let j2 = index.has_junction_overlap(chrom, pos_ref2);
            if j1 && j2 {
                self.counters.known_junction_events += 1;
            } else if j1 || j2 {
                self.counters.partly_known_junction_events += 1;
            }
            // Novel junctions (neither matches) are computed as:
            // reads_at_junctions - known - partly_known
        }

        // --- Cache tree query results for all M-blocks ---
        // Query the merged exon tree ONCE per block set and reuse the results
        // across enclosure check, coverage, overlapping-exon, and SSP.
        let cached_hits = cache_block_hits(&m_blocks, chrom, index);

        // Determine enclosing genes using cached results (with strand filtering)
        let is_reverse = flags & BAM_FREVERSE != 0;
        // SE reads lack BAM_FREAD1 (0x40), but should be treated as "first of
        // pair" for strand flip logic, matching Qualimap's getReadIntervals():
        //   boolean firstOfPair = true;
        //   if (pairedRead) { firstOfPair = read.getFirstOfPairFlag(); ... }
        //   else { if (protocol == STRAND_SPECIFIC_REVERSE) strand = !strand; }
        // See: bitbucket.org/kokonech/qualimap ComputeCountsTask.java
        let is_first_of_pair = if flags & BAM_FPAIRED != 0 {
            flags & BAM_FREAD1 != 0
        } else {
            true
        };
        let enclosing_genes =
            find_enclosing_genes_cached(&cached_hits, self.stranded, is_reverse, is_first_of_pair);

        // --- Coverage: per-block, strand-independent (Qualimap model) ---
        // Qualimap adds coverage per M-block to each enclosing gene, without
        // strand filtering. This is separate from gene assignment counting.
        // For SE reads: add coverage now. For PE: deferred to reconcile_pair().
        if flags & BAM_FPAIRED == 0 || flags & BAM_FPROPER_PAIR == 0 {
            self.add_per_block_coverage_cached(&m_blocks, &cached_hits, index);
        }

        // Bundle per-read data for downstream use
        let data = MateData {
            m_blocks,
            enclosing_genes,
            cached_hits,
            is_reverse,
            is_first_of_pair,
        };

        // PE path: check BAM flags directly (no library-level paired flag needed)
        if flags & BAM_FPAIRED != 0 && flags & BAM_FPROPER_PAIR != 0 {
            // PE proper pair: buffer and combine with mate (Qualimap name-sorts
            // and processes both mates' intervals together as a fragment).
            self.process_pe_read(record, data, chrom, index);
        } else {
            // SE mode or non-proper pair: assign immediately (numReads = 1)
            self.classify_and_count(&data, chrom, index, 1);
        }
    }

    /// Process a PE read — buffer first mate, reconcile when second arrives.
    fn process_pe_read(
        &mut self,
        record: &bam::Record,
        data: MateData,
        chrom: &str,
        index: &QualimapIndex,
    ) {
        let key = make_mate_key(record);

        if let Some(mate) = self.mate_buffer.remove(&key) {
            // Found the mate — reconcile the pair.
            // Current read is r1 (second-seen), buffered mate is r2 (first-seen).
            self.reconcile_pair(&data, &mate, chrom, index);
        } else {
            // First mate seen — buffer it
            self.mate_buffer.insert(key, data);
        }
    }

    /// Reconcile a PE pair: combine M-blocks and gene assignments.
    ///
    /// Qualimap intersects the enclosing gene sets of both mates.
    /// If exactly one gene encloses ALL M-blocks of BOTH mates → assigned.
    /// If >1 gene → ambiguous. If 0 genes → check intronic/intergenic.
    fn reconcile_pair(&mut self, r1: &MateData, r2: &MateData, chrom: &str, index: &QualimapIndex) {
        // Intersect the enclosing gene sets from both mates (already computed).
        let assigned_genes: HashSet<u32> = r1
            .enclosing_genes
            .intersection(&r2.enclosing_genes)
            .copied()
            .collect();

        // --- Coverage: per-block, strand-independent (Qualimap model) ---
        // Use cached hits from both mates to add coverage without re-querying the tree.
        self.add_per_block_coverage_cached(&r1.m_blocks, &r1.cached_hits, index);
        self.add_per_block_coverage_cached(&r2.m_blocks, &r2.cached_hits, index);

        // Check for overlapping_exon using cached hits from both mates.
        // At most one increment per fragment: check r2 only if r1 didn't trigger.
        if !check_has_overlap_not_enclosed(&r1.cached_hits) {
            self.check_overlapping_exon_cached(&r2.cached_hits);
        } else {
            self.counters.overlapping_exon += 1;
        }

        // Combine blocks for intronic/intergenic classification only when needed
        match assigned_genes.len() {
            0 => {
                // No gene encloses all blocks — classify as intronic or intergenic
                let mut combined_blocks: Vec<(i32, i32)> =
                    Vec::with_capacity(r1.m_blocks.len() + r2.m_blocks.len());
                combined_blocks.extend_from_slice(&r1.m_blocks);
                combined_blocks.extend_from_slice(&r2.m_blocks);
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

        // SSP estimation using cached hits — avoids re-querying the tree.
        self.count_ssp_for_blocks_cached(&r1.cached_hits, r1.is_reverse, r1.is_first_of_pair);
        let r2_is_first = !r1.is_first_of_pair;
        self.count_ssp_for_blocks_cached(&r2.cached_hits, r2.is_reverse, r2_is_first);
    }

    /// Classify and count a fragment (SE read or single-mate data).
    ///
    /// Shared logic for SE reads, non-proper pairs, and cross-chromosome
    /// mate reconciliation. `num_reads` is 1 for SE / non-proper and 2 for
    /// cross-chromosome PE fragments.
    fn classify_and_count(
        &mut self,
        data: &MateData,
        chrom: &str,
        index: &QualimapIndex,
        num_reads: u64,
    ) {
        // Check for overlapping_exon using cached results.
        self.check_overlapping_exon_cached(&data.cached_hits);

        match data.enclosing_genes.len() {
            0 => {
                self.classify_no_feature(&data.m_blocks, chrom, index, num_reads);
                self.counters.no_feature += num_reads;
            }
            1 => {
                self.counters.exonic_reads += num_reads;
            }
            _ => {
                self.counters.ambiguous += num_reads;
            }
        }

        // SSP estimation using cached hits.
        self.count_ssp_for_blocks_cached(&data.cached_hits, data.is_reverse, data.is_first_of_pair);
    }

    /// Check if any M-block overlaps an exon region without being fully enclosed,
    /// using pre-computed cached block hits. Increments the counter at most once.
    fn check_overlapping_exon_cached(&mut self, cached_hits: &[CachedBlockHits]) {
        if check_has_overlap_not_enclosed(cached_hits) {
            self.counters.overlapping_exon += 1;
        }
    }

    /// Add per-transcript and merged-gene coverage, per-block, strand-independent.
    ///
    /// Uses pre-computed cached block hits to avoid redundant tree queries.
    /// Matches Qualimap's coverage model exactly: coverage is added per-block
    /// for every gene whose merged exon encloses that individual M-block.
    ///
    /// # Arguments
    /// * `m_blocks`    - M-only aligned blocks (0-based half-open).
    /// * `cached_hits` - Pre-computed tree query results per block.
    /// * `index`       - Qualimap annotation index.
    fn add_per_block_coverage_cached(
        &mut self,
        m_blocks: &[(i32, i32)],
        cached_hits: &[CachedBlockHits],
        index: &QualimapIndex,
    ) {
        for (i, &(block_start, block_end)) in m_blocks.iter().enumerate() {
            let hit = &cached_hits[i];
            let single_block = &[(block_start, block_end)];
            for &(gidx, _strand) in &hit.enclosing_genes {
                let (tx_start, _) = index.gene_transcript_ranges[gidx as usize];
                let transcripts = index.gene_transcripts(gidx);
                self.coverage
                    .add_coverage(single_block, transcripts, tx_start);
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

    /// Count SSP votes using pre-computed cached block hits.
    ///
    /// Uses the `(gene_idx, strand)` pairs from cached hits directly — each gene
    /// appears at most once per block (deduplicated in `cache_block_hits`), so no
    /// additional deduplication is needed.
    fn count_ssp_for_blocks_cached(
        &mut self,
        cached_hits: &[CachedBlockHits],
        is_reverse: bool,
        is_first_of_pair: bool,
    ) {
        for hit in cached_hits {
            for &(_gidx, strand) in &hit.enclosing_genes {
                let same_strand = (is_reverse && strand == b'-') || (!is_reverse && strand == b'+');
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

        // Merge mate buffers (cross-chromosome mates).
        // When both mates were buffered in different per-chromosome batches,
        // reconcile them here. These are PE fragments, so num_reads = 2.
        for (key, info) in other.mate_buffer {
            if let Some(existing) = self.mate_buffer.remove(&key) {
                // We don't have the chrom handy for intronic classification,
                // so treat no-feature conservatively as intergenic.
                let common_genes: HashSet<u32> = existing
                    .enclosing_genes
                    .intersection(&info.enclosing_genes)
                    .copied()
                    .collect();

                // Union fallback: if one mate has no enclosing genes, use
                // the other mate's set (cross-chrom pairs often have one
                // mate in an unannotated region).
                let assigned = if existing.enclosing_genes.is_empty()
                    && !info.enclosing_genes.is_empty()
                {
                    info.enclosing_genes.clone()
                } else if info.enclosing_genes.is_empty() && !existing.enclosing_genes.is_empty() {
                    existing.enclosing_genes.clone()
                } else {
                    common_genes
                };

                // PE fragment: count both reads (num_reads = 2)
                match assigned.len() {
                    0 => {
                        self.counters.no_feature += 2;
                        self.counters.intergenic_reads += 2;
                    }
                    1 => {
                        self.counters.exonic_reads += 2;
                    }
                    _ => {
                        self.counters.ambiguous += 2;
                    }
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
    pub fn into_result(self) -> super::QualimapResult {
        // Convert junction motif byte arrays to readable strings
        let mut junction_motifs_str = HashMap::new();
        for ((donor, acceptor), count) in &self.junction_motifs.counts {
            let motif = format!(
                "{}{}{}{}",
                donor[0] as char, donor[1] as char, acceptor[0] as char, acceptor[1] as char,
            );
            *junction_motifs_str.entry(motif).or_insert(0u64) += count;
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
            left_proper_in_pair: self.counters.left_proper_in_pair,
            right_proper_in_pair: self.counters.right_proper_in_pair,
            both_proper_in_pair: self.counters.both_proper_in_pair,
            reads_at_junctions: self.counters.reads_at_junctions,
            known_junction_events: self.counters.known_junction_events,
            partly_known_junction_events: self.counters.partly_known_junction_events,
            junction_motifs: junction_motifs_str,
            transcript_coverage_raw: self.coverage,
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
    crate::rna::bam_flags::get_aux_int(record, b"NH").map(|v| v as u32)
}

/// Extract M-only aligned blocks from a BAM record's CIGAR.
///
/// Returns a vector of (start, end) in 0-based half-open coordinates.
/// Only `M` (alignment match) operations are included.
///
/// Extract correct alignment blocks from a BAM record's CIGAR string.
///
/// Unlike `extract_m_blocks_and_junctions`, this function uses **correct** CIGAR semantics:
/// only M/D/N/=/X operations advance the reference position, while I/S/H/P
/// **Qualimap compatibility note:** Qualimap's Java code advances the reference
/// offset for **all** CIGAR operations, including insertions and soft-clips.
/// This is technically incorrect (I/S consume only query, not reference), but
/// we replicate it here to produce identical enclosure results. The effect is
/// that M-blocks following an I or S operation are shifted rightward by the
/// length of the I/S, which can push them past exon boundaries and change
/// the exonic/intronic classification.
///
/// Also extracts junction motifs, n_op_count, and junction reference positions
/// in the same CIGAR pass.
///
/// Returns `(m_blocks, n_op_count, motifs, junction_ref_positions)`.
/// Each junction position pair `(pos_ref1, pos_ref2)` matches Qualimap Java's
/// `posInRef1 = alignmentStart + posInRead` and `posInRef2 = posInRef1 + clippedLength`.
#[allow(
    clippy::needless_range_loop,
    clippy::type_complexity,
    unused_assignments,
    unused_variables
)]
fn extract_m_blocks_and_junctions(
    record: &bam::Record,
) -> (Vec<(i32, i32)>, usize, Vec<JunctionMotif>, Vec<(i32, i32)>) {
    let cigar = record.cigar();
    let seq = record.seq();
    let seq_len = seq.len();

    let mut blocks = Vec::new();
    let mut motifs = Vec::new();
    let mut junction_ref_positions = Vec::new();
    let mut n_op_count: usize = 0;
    let mut ref_pos = record.pos() as i32; // 0-based
    let alignment_start_1based = record.pos() as i32 + 1;
    // seq_pos tracks query position for junction motif extraction and
    // junction position calculation (Qualimap's posInRead).
    // ref_pos for extract_m_blocks uses its own tracking with the Qualimap
    // bug (I/S/H/P advance ref_pos too).
    let mut seq_pos: usize = 0;
    // Qualimap-bug ref_pos tracks the shifted position for M-block extraction.
    // Junction motif extraction uses a separate ref tracking (standard).
    // We run both in one pass by tracking seq_pos alongside ref_pos.

    for op in cigar.iter() {
        match op {
            // Qualimap only creates alignment blocks for M (match/mismatch),
            // not for = (sequence match) or X (sequence mismatch).
            Cigar::Match(len) => {
                let len = *len as i32;
                blocks.push((ref_pos, ref_pos + len));
                ref_pos += len;
                seq_pos += len as usize;
            }
            // Qualimap bug: EQ and X advance offset but don't create blocks.
            Cigar::Equal(len) | Cigar::Diff(len) => {
                ref_pos += *len as i32;
                seq_pos += *len as usize;
            }
            Cigar::RefSkip(len) => {
                // Junction motif extraction (N-operations).
                // Qualimap increments numReadsWithJunction for every N-op,
                // regardless of whether the motif can be extracted.
                n_op_count += 1;

                // Record junction reference positions for known/novel classification.
                // Matches Qualimap Java: posInRef1 = alignmentStart + posInRead,
                // posInRef2 = posInRef1 + clippedLength.
                let pos_ref1 = alignment_start_1based + seq_pos as i32;
                let pos_ref2 = pos_ref1 + *len as i32;
                junction_ref_positions.push((pos_ref1, pos_ref2));

                // Extract donor motif: 2bp at end of preceding exon in read
                // Extract acceptor motif: 2bp at start of next exon in read
                // Guard matches Qualimap: posInRead - 2 > 0, i.e. posInRead >= 3.
                // This requires at least 3 query bases before the junction.
                if seq_pos > 2 && seq_pos + 1 < seq_len {
                    let donor = [seq.encoded_base(seq_pos - 2), seq.encoded_base(seq_pos - 1)];
                    let acceptor = [seq.encoded_base(seq_pos), seq.encoded_base(seq_pos + 1)];

                    let donor_ascii = [decode_base(donor[0]), decode_base(donor[1])];
                    let acceptor_ascii = [decode_base(acceptor[0]), decode_base(acceptor[1])];

                    motifs.push(JunctionMotif {
                        donor: donor_ascii,
                        acceptor: acceptor_ascii,
                    });
                }

                ref_pos += *len as i32;
                // N does not advance query (seq_pos unchanged)
            }
            // Qualimap bug: insertions and soft-clips advance the reference
            // position even though they only consume query bases.
            Cigar::Ins(len) | Cigar::SoftClip(len) => {
                ref_pos += *len as i32;
                seq_pos += *len as usize;
            }
            Cigar::Del(len) => {
                ref_pos += *len as i32;
                // D does not advance query
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    (blocks, n_op_count, motifs, junction_ref_positions)
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

/// Check whether any cached block hit has an overlap that is not fully enclosed.
///
/// Pure predicate — does not modify any counters. Used in `reconcile_pair` to
/// avoid double-incrementing `overlapping_exon` when checking both mates.
fn check_has_overlap_not_enclosed(cached_hits: &[CachedBlockHits]) -> bool {
    cached_hits.iter().any(|h| h.has_overlap_not_enclosed)
}

/// Query the merged exon tree ONCE per M-block and cache the results.
///
/// This single tree query per block replaces the 4+ separate queries that were
/// previously done in `find_enclosing_genes`, `add_per_block_coverage`,
/// `check_overlapping_exon`, and `count_ssp_for_blocks`.
fn cache_block_hits(
    m_blocks: &[(i32, i32)],
    chrom: &str,
    index: &QualimapIndex,
) -> Vec<CachedBlockHits> {
    let tree = match index.merged_exon_tree(chrom) {
        Some(t) => t,
        None => {
            return m_blocks
                .iter()
                .map(|_| CachedBlockHits {
                    enclosing_genes: Vec::new(),
                    has_overlap_not_enclosed: false,
                })
                .collect();
        }
    };

    m_blocks
        .iter()
        .map(|&(block_start, block_end)| {
            let mut enclosing_genes: Vec<(u32, u8)> = Vec::new();
            let mut has_overlap_not_enclosed = false;

            tree.query(block_start, block_end, |iv| {
                let exon_start = iv.metadata.exon_start;
                let exon_end = iv.metadata.exon_end;

                // Verify actual overlap in half-open coordinates
                if exon_start >= block_end || block_start >= exon_end {
                    return; // Abutting or non-overlapping — skip
                }

                // Check enclosure
                if block_start >= exon_start && block_end <= exon_end {
                    let gene_idx = iv.metadata.gene_idx;
                    let strand = iv.metadata.strand;
                    // Deduplicate by gene_idx (a gene can have multiple
                    // merged exons enclosing the same block)
                    if !enclosing_genes.iter().any(|(g, _)| *g == gene_idx) {
                        enclosing_genes.push((gene_idx, strand));
                    }
                } else {
                    // Overlaps but does NOT enclose
                    has_overlap_not_enclosed = true;
                }
            });

            CachedBlockHits {
                enclosing_genes,
                has_overlap_not_enclosed,
            }
        })
        .collect()
}

/// Find genes whose merged exons enclose ALL M-blocks of a read, using cached hits.
///
/// Operates on pre-computed `CachedBlockHits` (which store `(gene_idx, strand)`
/// tuples) instead of querying the tree directly. A gene must enclose every
/// M-block and pass the strand filter to be included.
fn find_enclosing_genes_cached(
    cached_hits: &[CachedBlockHits],
    stranded: Strandedness,
    is_reverse: bool,
    is_first_of_pair: bool,
) -> HashSet<u32> {
    if cached_hits.is_empty() {
        return HashSet::new();
    }

    // Compute the effective strand of this read after protocol correction.
    // Qualimap's getReadIntervals() flips the strand of:
    //   - read1 for strand-specific-reverse (Reverse)
    //   - read2 for strand-specific-forward (Forward)
    let read_on_plus = if stranded == Strandedness::Unstranded {
        true // no filtering; handled by skipping strand check below
    } else {
        let flip = (stranded == Strandedness::Reverse && is_first_of_pair)
            || (stranded == Strandedness::Forward && !is_first_of_pair);
        if flip {
            is_reverse
        } else {
            !is_reverse
        }
    };

    // Closure: does a gene's strand match this read?
    let strand_ok = |gene_strand: u8| -> bool {
        if stranded == Strandedness::Unstranded {
            return true;
        }
        (read_on_plus && gene_strand == b'+') || (!read_on_plus && gene_strand == b'-')
    };

    // First M-block: seed candidates from strand-compatible enclosing genes
    let mut candidate_genes: HashSet<u32> = cached_hits[0]
        .enclosing_genes
        .iter()
        .filter(|(_, strand)| strand_ok(*strand))
        .map(|(gidx, _)| *gidx)
        .collect();

    if candidate_genes.is_empty() {
        return candidate_genes;
    }

    // Subsequent M-blocks: intersect — keep only genes that enclose every block
    for hit in &cached_hits[1..] {
        candidate_genes.retain(|gidx| {
            hit.enclosing_genes
                .iter()
                .any(|(g, strand)| g == gidx && strand_ok(*strand))
        });
        if candidate_genes.is_empty() {
            break;
        }
    }

    candidate_genes
}

/// Find genes whose merged exons enclose ALL M-blocks of a read (original version).
///
/// This non-cached version is kept for use in unit tests and the
/// `find_enclosing_genes_empty_chrom` test.
#[cfg(test)]
fn find_enclosing_genes(
    m_blocks: &[(i32, i32)],
    chrom: &str,
    index: &QualimapIndex,
    stranded: Strandedness,
    is_reverse: bool,
    is_first_of_pair: bool,
) -> HashSet<u32> {
    let cached_hits = cache_block_hits(m_blocks, chrom, index);
    find_enclosing_genes_cached(&cached_hits, stranded, is_reverse, is_first_of_pair)
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
    crate::io::fnv1a(qname)
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
        let result = find_enclosing_genes(
            &[(100, 200)],
            "chr1",
            &index,
            Strandedness::Unstranded,
            false,
            true,
        );
        assert!(result.is_empty());
    }

    #[test]
    fn test_counters_merge() {
        let mut c1 = QualimapCounters {
            primary_alignments: 10,
            exonic_reads: 5,
            intronic_reads: 3,
            ..Default::default()
        };

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
