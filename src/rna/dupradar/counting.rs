//! Read counting engine for alignment files (SAM/BAM/CRAM).
//!
//! Assigns reads from an alignment file to genes based on GTF annotation,
//! producing four count vectors matching the dupRadar approach:
//! 1. All reads including multimappers (with duplicates)
//! 2. All reads including multimappers (without duplicates)
//! 3. Uniquely mapped reads only (with duplicates)
//! 4. Uniquely mapped reads only (without duplicates)
//!
//! This implements a simplified featureCounts-compatible counting strategy.

use crate::gtf::Gene;
use crate::rna::qualimap::QualimapAccum;
use crate::rna::rseqc::accumulators::{RseqcAccumulators, RseqcAnnotations, RseqcConfig};
use anyhow::{Context, Result};
use coitrees::{COITree, Interval, IntervalTree};
use indexmap::IndexMap;
use log::{debug, info, warn};
use rayon::prelude::*;
use rust_htslib::bam::{self, FetchDefinition, Read as BamRead};
use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, Ordering};

// ===================================================================
// Gene ID interning
// ===================================================================

/// Interned gene ID index. Avoids cloning gene ID strings in the hot loop.
type GeneIdx = u32;

/// Maps gene ID strings to compact integer indices for fast, allocation-free
/// lookups in the counting hot path.
#[derive(Debug)]
struct GeneIdInterner {
    /// Maps gene ID string to its index
    id_to_idx: HashMap<String, GeneIdx>,
    /// Maps index back to gene ID string (for output)
    idx_to_id: Vec<String>,
}

impl GeneIdInterner {
    /// Build an interner from the gene map, preserving insertion order.
    fn from_genes(genes: &IndexMap<String, Gene>) -> Self {
        let mut id_to_idx = HashMap::with_capacity(genes.len());
        let mut idx_to_id = Vec::with_capacity(genes.len());
        for (i, gene_id) in genes.keys().enumerate() {
            id_to_idx.insert(gene_id.clone(), i as GeneIdx);
            idx_to_id.push(gene_id.clone());
        }
        GeneIdInterner {
            id_to_idx,
            idx_to_id,
        }
    }

    /// Look up the index for a gene ID. Returns None if the gene is unknown.
    fn get(&self, gene_id: &str) -> Option<GeneIdx> {
        self.id_to_idx.get(gene_id).copied()
    }

    /// Get the gene ID string for an index.
    fn name(&self, idx: GeneIdx) -> &str {
        &self.idx_to_id[idx as usize]
    }

    /// Number of interned gene IDs.
    fn len(&self) -> usize {
        self.idx_to_id.len()
    }
}

// ===================================================================
// BAM header validation
// ===================================================================

/// Known duplicate-marking tool identifiers.
///
/// These strings are matched (case-insensitively) against the `ID:` and `PN:`
/// (program name) fields of `@PG` header entries. Matching is restricted to these
/// fields to avoid false positives from command-line strings or unrelated tools
/// (e.g., `picard SortSam` or `sambamba sort`).
const KNOWN_DUP_MARKERS: &[&str] = &[
    "markduplicates",
    "samblaster",
    "sambamba markdup",
    "sambamba_markdup",
    "biobambam",
    "estreamer",
    "fgbio",
    "umis",
    "umi_tools",
    "umi-tools",
    "gencore",
    "gatk markduplicates",
    "sentieon dedup",
];

/// Extracts the value of a specific tag from a SAM header line.
///
/// For example, given the line `@PG\tID:MarkDuplicates\tPN:MarkDuplicates` and
/// the tag `"ID"`, returns `Some("MarkDuplicates")`.
fn extract_header_tag<'a>(line: &'a str, tag: &str) -> Option<&'a str> {
    let prefix = format!("{}:", tag);
    line.split('\t')
        .find(|field| field.starts_with(&prefix))
        .map(|field| &field[prefix.len()..])
}

/// Checks whether a SAM header text contains evidence of a duplicate-marking tool.
///
/// Only the `ID:` and `PN:` fields of `@PG` lines are inspected to avoid false
/// positives from command-line arguments or unrelated tools that happen to share
/// a name (e.g., `picard SortSam`).
///
/// Returns `true` if a known duplicate-marking tool is found.
fn header_has_dup_marker(header_text: &str) -> bool {
    for line in header_text.lines() {
        // Only inspect @PG (program) header lines
        if !line.starts_with("@PG") {
            continue;
        }

        // Extract only the ID and PN fields for matching
        let id = extract_header_tag(line, "ID").unwrap_or("");
        let pn = extract_header_tag(line, "PN").unwrap_or("");
        let id_lower = id.to_lowercase();
        let pn_lower = pn.to_lowercase();

        if KNOWN_DUP_MARKERS
            .iter()
            .any(|marker| id_lower.contains(marker) || pn_lower.contains(marker))
        {
            debug!("Found duplicate-marking tool in @PG header: {}", line);
            return true;
        }
    }
    false
}

/// Checks the BAM `@PG` header lines for evidence that a duplicate-marking tool has been run.
///
/// Returns `Ok(())` if a known duplicate-marking tool is found, or an error with
/// a descriptive message if none is detected.
///
/// # Arguments
///
/// * `header` - The BAM header view to inspect
/// * `bam_path` - The BAM file path (used in the error message)
fn verify_duplicates_marked(header: &bam::HeaderView, bam_path: &str) -> Result<()> {
    let header_text = String::from_utf8_lossy(header.as_bytes());
    if header_has_dup_marker(&header_text) {
        return Ok(());
    }

    anyhow::bail!(
        "No duplicate-marking tool found in BAM header of '{}'.\n\
         \n\
         RustQC requires that BAM files have duplicates marked (SAM flag 0x400)\n\
         but NOT removed. The BAM @PG header lines do not contain evidence of a\n\
         known duplicate-marking tool.\n\
         \n\
         Please run one of the following tools before using RustQC:\n\
         \n\
           - Picard MarkDuplicates: picard MarkDuplicates I=input.bam O=marked.bam M=metrics.txt\n\
           - samblaster:            samtools view -h input.bam | samblaster | samtools view -bS - > marked.bam\n\
           - sambamba markdup:      sambamba markdup input.bam marked.bam\n\
         \n\
         If you are certain that duplicates are already marked, use --skip-dup-check to bypass this check.",
        bam_path
    )
}

// ===================================================================
// BAM flag constants
// ===================================================================

use crate::rna::bam_flags::*;

/// Counts for a single gene across the four counting modes.
#[derive(Debug, Clone, Default)]
pub struct GeneCounts {
    /// Count with multimappers, with duplicates
    pub all_multi: u64,
    /// Count with multimappers, without duplicates (dups excluded)
    pub nodup_multi: u64,
    /// Count without multimappers, with duplicates
    pub all_unique: u64,
    /// Count without multimappers, without duplicates (dups excluded)
    pub nodup_unique: u64,

    // --- featureCounts per-read counting ---
    // featureCounts with `-p` (no `--countReadPairs`) processes each read
    // independently in single-end mode. These counters track per-read
    // assignments for featureCounts-compatible output.
    /// Per-read count: uniquely-mapped reads assigned to this gene (all reads
    /// including duplicates). Matches featureCounts' `-p` behaviour where each
    /// read is independently assigned.
    pub fc_reads: u64,
}

impl GeneCounts {
    /// Increment the appropriate counters for a fragment assigned to this gene.
    fn count_fragment(&mut self, is_dup: bool, is_multi: bool) {
        self.all_multi += 1;
        if !is_dup {
            self.nodup_multi += 1;
        }
        if !is_multi {
            self.all_unique += 1;
            if !is_dup {
                self.nodup_unique += 1;
            }
        }
    }
}

/// Merge gene hits from two mates into a combined list of best-scoring genes.
///
/// Each input is a sorted, deduplicated slice of gene indices that a read
/// overlaps. Each gene gets +1 per mate that overlaps it. Genes present in
/// both mates get score 2, genes in only one get score 1. Returns the gene(s)
/// with the highest combined score. This replaces a per-fragment `HashMap`
/// with a sorted-merge algorithm optimized for the common case of 0-2 gene
/// hits per mate.
fn merge_gene_hits(a: &[GeneIdx], b: &[GeneIdx]) -> Vec<GeneIdx> {
    // Fast paths for the common cases
    if a.is_empty() {
        return b.to_vec();
    }
    if b.is_empty() {
        return a.to_vec();
    }

    // Both non-empty: sorted merge to find genes in both mates (score 2)
    // vs genes in only one mate (score 1). Prefer genes in both.
    let mut both = Vec::new();
    let mut ai = 0;
    let mut bi = 0;
    while ai < a.len() && bi < b.len() {
        match a[ai].cmp(&b[bi]) {
            std::cmp::Ordering::Equal => {
                both.push(a[ai]);
                ai += 1;
                bi += 1;
            }
            std::cmp::Ordering::Less => ai += 1,
            std::cmp::Ordering::Greater => bi += 1,
        }
    }

    // If any gene was hit by both mates, return only those (score 2 > score 1)
    if !both.is_empty() {
        return both;
    }

    // No overlap: all genes have score 1, return the union
    let mut union = Vec::with_capacity(a.len() + b.len());
    ai = 0;
    bi = 0;
    while ai < a.len() && bi < b.len() {
        match a[ai].cmp(&b[bi]) {
            std::cmp::Ordering::Less => {
                union.push(a[ai]);
                ai += 1;
            }
            std::cmp::Ordering::Greater => {
                union.push(b[bi]);
                bi += 1;
            }
            std::cmp::Ordering::Equal => {
                // Should not happen (handled above), but be safe
                union.push(a[ai]);
                ai += 1;
                bi += 1;
            }
        }
    }
    union.extend_from_slice(&a[ai..]);
    union.extend_from_slice(&b[bi..]);
    union
}

/// Assign a fragment to a gene if exactly one gene was hit.
/// Returns true if the fragment was assigned to a gene.
fn assign_fragment_to_gene(
    gene_hits: &[GeneIdx],
    gene_counts: &mut [GeneCounts],
    is_dup: bool,
    is_multi: bool,
) -> bool {
    if gene_hits.len() == 1 {
        let idx = gene_hits[0] as usize;
        if idx < gene_counts.len() {
            gene_counts[idx].count_fragment(is_dup, is_multi);
            return true;
        }
    }
    false
}

/// Result of the counting step, including per-gene counts and assignment statistics.
///
/// Contains everything needed to produce both dupRadar outputs and
/// featureCounts-compatible output files, plus optional RSeQC accumulator results
/// collected during the same BAM pass.
#[derive(Debug)]
pub struct CountResult {
    /// Per-gene counts indexed by gene_id
    pub gene_counts: IndexMap<String, GeneCounts>,
    /// Total mapped reads (including multimappers, including duplicates)
    pub total_reads_multi_dup: u64,
    /// Total mapped reads (unique only, including duplicates)
    pub total_reads_unique_dup: u64,

    // --- dupRadar fragment-level assignment statistics ---
    /// Total reads/records seen in the BAM file
    pub stat_total_reads: u64,
    /// Fragments successfully assigned to exactly one gene
    pub stat_assigned: u64,
    /// Fragments overlapping multiple genes (ambiguous)
    pub stat_ambiguous: u64,
    /// Fragments overlapping no annotated gene
    pub stat_no_features: u64,
    /// Total fragments (single-end reads or paired-end pairs)
    pub stat_total_fragments: u64,

    // --- featureCounts per-read statistics ---
    // featureCounts with `-p` (no `--countReadPairs`) counts each read
    // independently. These stats match that behaviour for compatible output.
    /// Per-read: reads assigned to exactly one gene (non-multimapped)
    pub fc_assigned: u64,
    /// Per-read: reads overlapping multiple genes (non-multimapped)
    pub fc_ambiguous: u64,
    /// Per-read: reads overlapping no gene (non-multimapped)
    pub fc_no_features: u64,
    /// Per-read: multimapping reads (NH > 1)
    pub fc_multimapping: u64,
    /// Per-read: unmapped reads
    pub fc_unmapped: u64,

    // --- RSeQC results (collected during the same BAM pass) ---
    /// Accumulated RSeQC tool results, if RSeQC tools were enabled
    pub rseqc: Option<RseqcAccumulators>,
    /// Qualimap RNA-Seq QC results (if enabled).
    #[allow(dead_code)]
    pub qualimap: Option<crate::rna::qualimap::QualimapResult>,
}

/// Metadata stored with each interval in the cache-oblivious interval tree.
#[derive(Debug, Clone, Copy, Default)]
struct IvMeta {
    /// Interned gene ID index (avoids String cloning in hot path)
    gene_idx: GeneIdx,
    /// Strand of the exon ('+', '-', or '.')
    strand: char,
}

/// Per-chromosome interval index backed by a cache-oblivious interval tree
/// (coitrees). Provides very fast overlap queries on sorted genomic intervals.
type ChromIndex = COITree<IvMeta, u32>;

/// Build a spatial index from gene annotations for fast overlap queries.
///
/// Constructs a `COITree` per chromosome from exon intervals. Uses the
/// interner to store compact gene ID indices as interval metadata,
/// avoiding per-exon String cloning.
fn build_index(
    genes: &IndexMap<String, Gene>,
    interner: &GeneIdInterner,
) -> HashMap<String, ChromIndex> {
    let mut chrom_intervals: HashMap<String, Vec<Interval<IvMeta>>> = HashMap::new();

    for gene in genes.values() {
        // Look up the interned gene index once per gene
        let gene_idx = interner
            .get(&gene.gene_id)
            .expect("gene must be in interner"); // safe: interner built from same gene map

        // Add each exon as an interval (featureCounts assigns reads at exon level)
        for exon in &gene.exons {
            // Convert from 1-based inclusive GTF to 0-based half-open, then
            // to coitrees' end-inclusive coordinates: [start-1, end-1].
            let iv = Interval::new(
                (exon.start - 1) as i32,
                (exon.end - 1) as i32,
                IvMeta {
                    gene_idx,
                    strand: exon.strand,
                },
            );
            chrom_intervals
                .entry(exon.chrom.clone())
                .or_default()
                .push(iv);
        }
    }

    chrom_intervals
        .into_iter()
        .map(|(chrom, intervals)| (chrom, COITree::new(&intervals)))
        .collect()
}

/// Determine if a read's strand matches the expected gene strand,
/// given the library strandedness.
///
/// # Arguments
/// * `read_reverse` - Whether the read is mapped to the reverse strand
/// * `is_read1` - Whether this is the first read in a pair (for paired-end)
/// * `paired` - Whether the library is paired-end
/// * `gene_strand` - The strand of the gene ('+' or '-')
/// * `stranded` - Library strandedness (0=unstranded, 1=forward, 2=reverse)
fn strand_matches(
    read_reverse: bool,
    is_read1: bool,
    paired: bool,
    gene_strand: char,
    stranded: u8,
) -> bool {
    if stranded == 0 || gene_strand == '.' {
        return true; // Unstranded: everything matches
    }

    // Determine the "effective" strand of the read
    // For paired-end: read1 maps to the transcript strand, read2 maps to the opposite
    // For single-end: the read maps directly
    let read_on_plus = !read_reverse;

    let effective_plus = if paired {
        if is_read1 {
            read_on_plus
        } else {
            !read_on_plus // Read2 is the complement
        }
    } else {
        read_on_plus
    };

    match stranded {
        1 => {
            // Forward stranded: read (or read1) aligns to the same strand as the gene
            (effective_plus && gene_strand == '+') || (!effective_plus && gene_strand == '-')
        }
        2 => {
            // Reverse stranded: read (or read1) aligns to the opposite strand of the gene
            (!effective_plus && gene_strand == '+') || (effective_plus && gene_strand == '-')
        }
        _ => true,
    }
}

/// Extract aligned (match) blocks from a CIGAR string into a reusable buffer.
///
/// Walks the CIGAR operations and appends genomic intervals for each
/// M/=/X operation to the provided buffer (which is cleared first).
/// N (intron skip) and D (deletion) operations advance the reference
/// position without producing an aligned block.
/// This ensures spliced reads don't falsely overlap with genes in introns.
///
/// # Arguments
/// * `start` - The reference start position of the read
/// * `cigar` - The CIGAR string view from the alignment record
/// * `blocks` - Reusable buffer for aligned blocks (cleared before use)
fn cigar_to_aligned_blocks(
    start: u64,
    cigar: &rust_htslib::bam::record::CigarStringView,
    blocks: &mut Vec<(u64, u64)>,
) {
    use rust_htslib::bam::record::Cigar;
    blocks.clear();
    let mut ref_pos = start;

    for op in cigar.iter() {
        match op {
            // Alignment match (M), sequence match (=), sequence mismatch (X):
            // These consume both reference and query bases. They represent
            // actual aligned segments.
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len = *len as u64;
                blocks.push((ref_pos, ref_pos + len));
                ref_pos += len;
            }
            // Reference skip (N): intron in RNA-seq. Advances reference only.
            // Deletion (D): deletion from reference. Advances reference only.
            Cigar::RefSkip(len) | Cigar::Del(len) => {
                ref_pos += *len as u64;
            }
            // Insertion (I), soft clip (S), hard clip (H), padding (P):
            // These do not advance the reference position.
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
}

/// Metadata collected per read for paired-end mate buffering.
/// When counting paired-end fragments, we need to see both mates to
/// properly determine gene overlap (featureCounts checks both mates
/// independently). This struct stores the information we need from each mate.
#[derive(Debug)]
struct MateInfo {
    /// Interned gene ID indices this mate overlaps (after strand filtering, deduplicated)
    gene_hits: Vec<GeneIdx>,
    /// Whether this read is a PCR/optical duplicate
    is_dup: bool,
    /// Whether this read is a multimapper (NH > 1)
    is_multi: bool,
}

/// Key for the mate buffer, matching featureCounts' `SAM_pairer_get_read_full_name()`.
///
/// featureCounts pairs mates using: (read_name, r1_refID, r1_pos, r2_refID, r2_pos, HI_tag).
/// R1/R2 roles are determined by the FLAG 0x40 bit (BAM_FREAD1), so both mates of the
/// same alignment pair always compute an identical key. The HI (Hit Index) tag disambiguates
/// multiple alignment pairs of the same multi-mapped read.
///
/// We use a u64 hash of the read name instead of the full `Vec<u8>` to avoid per-read
/// heap allocations. Collisions are astronomically unlikely given the compound key
/// includes tid, pos, and hi_tag alongside the hash.
type MateBufferKey = (u64, i32, i64, i32, i64, i32);

/// Hash a read name (byte slice) into a u64 using FNV-1a.
/// This avoids allocating a `Vec<u8>` for every read in paired-end mode.
#[inline(always)]
fn hash_qname(qname: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325; // FNV offset basis
    for &b in qname {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3); // FNV prime
    }
    h
}

/// Accumulated results from processing a batch of chromosomes.
/// Each parallel worker produces one of these, which are then merged.
#[derive(Debug)]
struct ChromResult {
    /// Per-gene counts (flat vec indexed by GeneIdx)
    gene_counts: Vec<GeneCounts>,
    /// Unmatched mates left over from this chromosome batch
    unmatched_mates: HashMap<MateBufferKey, MateInfo>,
    /// Total reads seen (all reads, including skipped)
    total_reads: u64,
    /// Total mapped reads (after filtering)
    total_mapped: u64,
    /// Total duplicate reads
    total_dup: u64,
    /// Total multimapping reads
    total_multi: u64,
    /// Total fragments (single-end reads or paired-end pairs)
    total_fragments: u64,
    /// Fragments assigned to exactly one gene
    stat_assigned: u64,
    /// Fragments overlapping multiple genes
    stat_ambiguous: u64,
    /// Fragments overlapping no annotated gene
    stat_no_features: u64,
    /// Mapped reads/fragments per counting mode (for N calculation)
    n_multi_dup: u64,
    n_multi_nodup: u64,
    n_unique_dup: u64,
    n_unique_nodup: u64,

    // --- featureCounts per-read statistics ---
    // Each read is independently assessed (no mate-pair merging),
    // matching featureCounts' behaviour with `-p` (no `--countReadPairs`).
    fc_assigned: u64,
    fc_ambiguous: u64,
    fc_no_features: u64,
    fc_multimapping: u64,
    fc_unmapped: u64,
    /// Qualimap RNA-Seq QC accumulator (if enabled).
    qualimap: Option<QualimapAccum>,
}

impl ChromResult {
    /// Create a new empty result with the given number of genes.
    fn new(num_genes: usize) -> Self {
        ChromResult {
            gene_counts: vec![GeneCounts::default(); num_genes],
            unmatched_mates: HashMap::new(),
            total_reads: 0,
            total_mapped: 0,
            total_dup: 0,
            total_multi: 0,
            total_fragments: 0,
            stat_assigned: 0,
            stat_ambiguous: 0,
            stat_no_features: 0,
            n_multi_dup: 0,
            n_multi_nodup: 0,
            n_unique_dup: 0,
            n_unique_nodup: 0,
            fc_assigned: 0,
            fc_ambiguous: 0,
            fc_no_features: 0,
            fc_multimapping: 0,
            fc_unmapped: 0,
            qualimap: None,
        }
    }

    /// Merge another ChromResult into this one (additive).
    fn merge(&mut self, other: ChromResult) {
        for (i, counts) in other.gene_counts.into_iter().enumerate() {
            self.gene_counts[i].all_multi += counts.all_multi;
            self.gene_counts[i].nodup_multi += counts.nodup_multi;
            self.gene_counts[i].all_unique += counts.all_unique;
            self.gene_counts[i].nodup_unique += counts.nodup_unique;
            self.gene_counts[i].fc_reads += counts.fc_reads;
        }
        // Merge unmatched mates — these will be reconciled in a separate step
        self.unmatched_mates.extend(other.unmatched_mates);
        self.total_reads += other.total_reads;
        self.total_mapped += other.total_mapped;
        self.total_dup += other.total_dup;
        self.total_multi += other.total_multi;
        self.total_fragments += other.total_fragments;
        self.stat_assigned += other.stat_assigned;
        self.stat_ambiguous += other.stat_ambiguous;
        self.stat_no_features += other.stat_no_features;
        self.n_multi_dup += other.n_multi_dup;
        self.n_multi_nodup += other.n_multi_nodup;
        self.n_unique_dup += other.n_unique_dup;
        self.n_unique_nodup += other.n_unique_nodup;
        self.fc_assigned += other.fc_assigned;
        self.fc_ambiguous += other.fc_ambiguous;
        self.fc_no_features += other.fc_no_features;
        self.fc_multimapping += other.fc_multimapping;
        self.fc_unmapped += other.fc_unmapped;
        // Merge Qualimap accumulator
        if let Some(other_qm) = other.qualimap {
            if let Some(ref mut self_qm) = self.qualimap {
                self_qm.merge(other_qm);
            } else {
                self.qualimap = Some(other_qm);
            }
        }
    }
}

/// Per-read featureCounts classification.
///
/// featureCounts with `-p` (no `--countReadPairs`) processes each read
/// independently. Multi-mapped reads (NH > 1) are categorised as
/// Unassigned_MultiMapping and excluded from gene assignment.
/// When multiple genes are hit, the read is Unassigned_Ambiguity
/// (matching default featureCounts `-g gene_id` behaviour where each
/// gene is its own meta-feature).
fn classify_read_fc(is_multi: bool, gene_hits: &[GeneIdx], result: &mut ChromResult) {
    if is_multi {
        result.fc_multimapping += 1;
    } else if gene_hits.is_empty() {
        result.fc_no_features += 1;
    } else if gene_hits.len() == 1 {
        result.fc_assigned += 1;
        let idx = gene_hits[0] as usize;
        if idx < result.gene_counts.len() {
            result.gene_counts[idx].fc_reads += 1;
        }
    } else {
        // Multiple gene hits → Ambiguous (default featureCounts behaviour)
        result.fc_ambiguous += 1;
    }
}

/// Process a batch of chromosomes from an alignment file, counting reads against gene annotations.
///
/// Opens its own indexed reader and seeks to each chromosome in the batch.
/// Returns accumulated results for all chromosomes in the batch, plus optional
/// RSeQC accumulator results collected during the same pass.
#[allow(clippy::too_many_arguments)]
fn process_chromosome_batch(
    bam_path: &str,
    tids: &[u32],
    tid_to_name: &[String],
    index: &HashMap<String, ChromIndex>,
    num_genes: usize,
    stranded: u8,
    paired: bool,
    chrom_mapping: &HashMap<String, String>,
    chrom_prefix: Option<&str>,
    global_read_counter: &AtomicU64,
    reference: Option<&str>,
    rseqc_config: Option<&RseqcConfig>,
    rseqc_annotations: Option<&RseqcAnnotations>,
    htslib_threads: usize,
    qualimap_index: Option<&crate::rna::qualimap::QualimapIndex>,
) -> Result<(ChromResult, Option<RseqcAccumulators>)> {
    let mut result = ChromResult::new(num_genes);
    if qualimap_index.is_some() {
        result.qualimap = Some(crate::rna::qualimap::QualimapAccum::new(stranded));
    }
    let mut rseqc_accums = rseqc_config.map(|cfg| RseqcAccumulators::new(cfg, rseqc_annotations));

    // Pre-compute resolved chromosome names per TID (apply prefix/mapping once,
    // not on every read). Used by both featureCounts counting and RSeQC tools.
    let tid_to_gtf_chrom: Vec<String> = tid_to_name
        .iter()
        .map(|bam_name| {
            if let Some(mapped) = chrom_mapping.get(bam_name.as_str()) {
                mapped.clone()
            } else if let Some(prefix) = chrom_prefix {
                format!("{}{}", prefix, bam_name)
            } else {
                bam_name.clone()
            }
        })
        .collect();
    let tid_to_rseqc_chrom = &tid_to_gtf_chrom;
    let tid_to_rseqc_chrom_upper: Vec<String> =
        tid_to_gtf_chrom.iter().map(|s| s.to_uppercase()).collect();

    // Open an indexed reader for this thread (supports BAM with .bai/.csi and CRAM with .crai)
    let mut bam = bam::IndexedReader::from_path(bam_path)
        .with_context(|| format!("Failed to open indexed alignment file: {}", bam_path))?;
    if let Some(ref_path) = reference {
        bam.set_reference(ref_path)
            .with_context(|| format!("Failed to set reference FASTA: {}", ref_path))?;
    }
    if htslib_threads > 0 {
        bam.set_threads(htslib_threads)
            .context("Failed to set htslib decompression threads")?;
    }

    // Reusable buffers
    let mut aligned_blocks_buf: Vec<(u64, u64)> = Vec::new();
    let mut gene_hits: Vec<GeneIdx> = Vec::new();
    let mut mate_buffer: HashMap<MateBufferKey, MateInfo> = HashMap::new();
    let mut record = bam::Record::new();

    for &tid in tids {
        // Seek to this chromosome
        bam.fetch(tid)
            .with_context(|| format!("Failed to seek to tid {}", tid))?;

        while let Some(read_result) = bam.read(&mut record) {
            read_result.context("Error reading alignment record")?;
            result.total_reads += 1;

            // Periodic progress logging (approximate, using atomic counter)
            let prev = global_read_counter.fetch_add(1, Ordering::Relaxed);
            if (prev + 1).is_multiple_of(5_000_000) {
                debug!("Processed ~{} reads...", prev + 1);
            }

            // --- RSeQC per-read dispatch (before counting filters) ---
            // RSeQC tools apply their own filter logic internally.
            // bam_stat and read_duplication need to see records that
            // counting.rs filters out (QC-fail, secondary, etc.).
            if let (Some(ref mut accums), Some(annots), Some(cfg)) =
                (&mut rseqc_accums, rseqc_annotations, rseqc_config)
            {
                let rec_tid = record.tid();
                let (chrom, chrom_upper) =
                    if rec_tid >= 0 && (rec_tid as usize) < tid_to_rseqc_chrom.len() {
                        (
                            tid_to_rseqc_chrom[rec_tid as usize].as_str(),
                            tid_to_rseqc_chrom_upper[rec_tid as usize].as_str(),
                        )
                    } else {
                        ("", "")
                    };
                accums.process_read(&record, chrom, chrom_upper, annots, cfg);
            }

            // --- Qualimap per-read dispatch (before counting filters) ---
            // Qualimap accumulator handles its own filtering (unmapped, secondary,
            // QC-fail, supplementary, NH>1) and uses enclosure-based gene assignment.
            if let Some(ref mut qm) = result.qualimap {
                let tid = record.tid();
                if tid >= 0 && (tid as usize) < tid_to_gtf_chrom.len() {
                    let qm_chrom = &tid_to_gtf_chrom[tid as usize];
                    qm.process_read(&record, qm_chrom, qualimap_index.unwrap());
                }
            }

            let flags = record.flags();

            // Skip unmapped reads (count for featureCounts summary)
            if flags & BAM_FUNMAP != 0 {
                result.fc_unmapped += 1;
                continue;
            }

            // Skip supplementary alignments (but NOT secondary - featureCounts processes
            // secondary alignments as separate counting events for multi-mapped reads)
            if flags & BAM_FSUPPLEMENTARY != 0 {
                continue;
            }

            // Skip QC-failed reads
            if flags & BAM_FQCFAIL != 0 {
                continue;
            }

            // For paired-end data, must actually be a paired read
            if paired && flags & BAM_FPAIRED == 0 {
                continue;
            }

            result.total_mapped += 1;

            let is_dup = flags & BAM_FDUP != 0;
            if is_dup {
                result.total_dup += 1;
            }

            // Determine if the read is a multimapper (NH tag)
            let is_multi = match record.aux(b"NH") {
                Ok(rust_htslib::bam::record::Aux::U8(nh)) => nh > 1,
                Ok(rust_htslib::bam::record::Aux::U16(nh)) => nh > 1,
                Ok(rust_htslib::bam::record::Aux::U32(nh)) => nh > 1,
                Ok(rust_htslib::bam::record::Aux::I8(nh)) => nh > 1,
                Ok(rust_htslib::bam::record::Aux::I16(nh)) => nh > 1,
                Ok(rust_htslib::bam::record::Aux::I32(nh)) => nh > 1,
                _ => false,
            };
            if is_multi {
                result.total_multi += 1;
            }

            let is_reverse = flags & BAM_FREVERSE != 0;
            let is_read1 = flags & BAM_FREAD1 != 0;

            // Get the chromosome name
            let rec_tid = record.tid();
            if rec_tid < 0 || rec_tid as usize >= tid_to_name.len() {
                continue;
            }
            // Use pre-computed GTF chromosome name (prefix/mapping applied once at init)
            let chrom = tid_to_gtf_chrom[rec_tid as usize].as_str();

            // Find gene overlaps using CIGAR-aware aligned blocks
            gene_hits.clear();
            if let Some(chrom_idx) = index.get(chrom) {
                cigar_to_aligned_blocks(
                    record.pos() as u64,
                    &record.cigar(),
                    &mut aligned_blocks_buf,
                );

                for &(block_start, block_end) in &aligned_blocks_buf {
                    chrom_idx.query(block_start as i32, (block_end - 1) as i32, |node| {
                        let meta = node.metadata;
                        if strand_matches(is_reverse, is_read1, paired, meta.strand, stranded) {
                            gene_hits.push(meta.gene_idx);
                        }
                    });
                }
                gene_hits.sort_unstable();
                gene_hits.dedup();
            }

            // --- Per-read featureCounts counting (independent of mate pairing) ---
            classify_read_fc(is_multi, &gene_hits, &mut result);

            // --- Single-end counting ---
            if !paired {
                result.n_multi_dup += 1;
                result.n_unique_dup += 1;
                if !is_dup {
                    result.n_multi_nodup += 1;
                    result.n_unique_nodup += 1;
                }
                result.total_fragments += 1;

                if gene_hits.is_empty() {
                    result.stat_no_features += 1;
                } else if gene_hits.len() > 1 {
                    result.stat_ambiguous += 1;
                } else if assign_fragment_to_gene(
                    &gene_hits,
                    &mut result.gene_counts,
                    is_dup,
                    is_multi,
                ) {
                    result.stat_assigned += 1;
                }
                continue;
            }

            // --- Paired-end counting: buffer mates and combine ---
            let qname_hash = hash_qname(record.qname());
            let own_pos = record.pos();
            let own_tid = record.tid();
            let mate_pos_val = record.mpos();
            let mate_tid = record.mtid();

            let hi_tag: i32 = match record.aux(b"HI") {
                Ok(rust_htslib::bam::record::Aux::U8(v)) => v as i32,
                Ok(rust_htslib::bam::record::Aux::U16(v)) => v as i32,
                Ok(rust_htslib::bam::record::Aux::U32(v)) => v as i32,
                Ok(rust_htslib::bam::record::Aux::I8(v)) => v as i32,
                Ok(rust_htslib::bam::record::Aux::I16(v)) => v as i32,
                Ok(rust_htslib::bam::record::Aux::I32(v)) => v,
                _ => -1,
            };

            let buffer_key: MateBufferKey = if is_read1 {
                (qname_hash, own_tid, own_pos, mate_tid, mate_pos_val, hi_tag)
            } else {
                (qname_hash, mate_tid, mate_pos_val, own_tid, own_pos, hi_tag)
            };

            if let Some(mate_info) = mate_buffer.remove(&buffer_key) {
                result.total_fragments += 1;

                let frag_is_dup = if is_read1 { is_dup } else { mate_info.is_dup };
                let frag_is_multi = if is_read1 {
                    is_multi
                } else {
                    mate_info.is_multi
                };

                result.n_multi_dup += 1;
                result.n_unique_dup += 1;
                if !frag_is_dup {
                    result.n_multi_nodup += 1;
                    result.n_unique_nodup += 1;
                }

                // featureCounts union-with-both-end-preference voting: each gene
                // gets +1 per mate that overlaps it (max 2).  Genes hit by both
                // mates (score 2) beat genes hit by only one (score 1).  If
                // multiple genes tie at the highest score the fragment is
                // ambiguous.  This matches Rsubread::featureCounts behaviour.
                let combined_genes: Vec<GeneIdx> =
                    merge_gene_hits(&mate_info.gene_hits, &gene_hits);

                if combined_genes.is_empty() {
                    result.stat_no_features += 1;
                } else if combined_genes.len() > 1 {
                    result.stat_ambiguous += 1;
                } else if assign_fragment_to_gene(
                    &combined_genes,
                    &mut result.gene_counts,
                    frag_is_dup,
                    frag_is_multi,
                ) {
                    result.stat_assigned += 1;
                }
            } else {
                // Take ownership of gene_hits instead of cloning (it's cleared at loop start)
                mate_buffer.insert(
                    buffer_key,
                    MateInfo {
                        gene_hits: std::mem::take(&mut gene_hits),
                        is_dup,
                        is_multi,
                    },
                );
            }
        }
    }

    // Move unmatched mates into the result for cross-chromosome reconciliation
    result.unmatched_mates = mate_buffer;
    Ok((result, rseqc_accums))
}

/// Partition chromosome indices across workers using greedy bin-packing
/// (largest-first scheduling). Assigns each chromosome to the worker with the
/// smallest current total length, producing a more balanced distribution than
/// round-robin when chromosome sizes vary widely.
fn partition_chromosomes(lengths: &[u64], num_workers: usize) -> Vec<Vec<u32>> {
    let mut batches: Vec<Vec<u32>> = vec![Vec::new(); num_workers];
    let mut loads: Vec<u64> = vec![0; num_workers];

    // Sort chromosome indices by length descending
    let mut order: Vec<u32> = (0..lengths.len() as u32).collect();
    order.sort_unstable_by(|&a, &b| lengths[b as usize].cmp(&lengths[a as usize]));

    // Assign each chromosome to the least-loaded worker
    for tid in order {
        let min_worker = loads
            .iter()
            .enumerate()
            .min_by_key(|&(_, &load)| load)
            .map(|(i, _)| i)
            .unwrap_or(0);
        batches[min_worker].push(tid);
        loads[min_worker] += lengths[tid as usize];
    }

    batches
}

///   an index for `threads > 1`; SAM files always use single-threaded mode)
/// * `genes` - Gene annotation map from GTF parsing
/// * `stranded` - Library strandedness (0, 1, or 2)
/// * `paired` - Whether the library is paired-end
/// * `threads` - Number of threads for alignment processing
/// * `reference` - Optional path to reference FASTA (required for CRAM files)
/// * `skip_dup_check` - If true, skip the BAM header check for duplicate-marking tools
#[allow(clippy::too_many_arguments)]
pub fn count_reads(
    bam_path: &str,
    genes: &IndexMap<String, Gene>,
    stranded: u8,
    paired: bool,
    threads: usize,
    chrom_mapping: &HashMap<String, String>,
    chrom_prefix: Option<&str>,
    reference: Option<&str>,
    skip_dup_check: bool,
    rseqc_config: Option<&RseqcConfig>,
    rseqc_annotations: Option<&RseqcAnnotations>,
    qualimap_index: Option<&crate::rna::qualimap::QualimapIndex>,
) -> Result<CountResult> {
    // Build gene ID interner for allocation-free lookups in the hot loop
    let interner = GeneIdInterner::from_genes(genes);

    // Build spatial index (uses interned gene indices)
    let index = build_index(genes, &interner);

    // Get chromosome names from header using a temporary reader,
    // and verify that duplicates have been marked in the BAM file.
    let (tid_to_name, tid_to_len): (Vec<String>, Vec<u64>) = {
        let mut bam = bam::Reader::from_path(bam_path)
            .with_context(|| format!("Failed to open alignment file: {}", bam_path))?;
        if let Some(ref_path) = reference {
            bam.set_reference(ref_path)
                .with_context(|| format!("Failed to set reference FASTA: {}", ref_path))?;
        }
        let header = bam.header().clone();

        // Check for evidence of duplicate-marking in @PG header lines
        if skip_dup_check {
            info!("Skipping duplicate-marking verification (--skip-dup-check)");
        } else {
            verify_duplicates_marked(&header, bam_path)?;
        }

        let names = (0..header.target_count())
            .map(|tid| String::from_utf8_lossy(header.tid2name(tid)).to_string())
            .collect();
        let lengths = (0..header.target_count())
            .map(|tid| header.target_len(tid).unwrap_or(0))
            .collect();
        (names, lengths)
    };

    // Check if an alignment index is available for parallel processing
    // (BAM uses .bai/.csi, CRAM uses .crai; SAM has no index format)
    let use_parallel = threads > 1 && {
        let idx_check = bam::IndexedReader::from_path(bam_path);
        if let Ok(mut reader) = idx_check {
            if let Some(ref_path) = reference {
                reader.set_reference(ref_path).is_ok()
            } else {
                true
            }
        } else {
            false
        }
    };
    if threads > 1 && !use_parallel {
        warn!(
            "Alignment index not found for {}; falling back to single-threaded processing. \
             Create an index with 'samtools index' to enable parallel processing.",
            bam_path
        );
    }

    // Build gene_idx → biotype_id lookup for per-read featureCounts counting.
    // When featureCounts runs with `-g gene_biotype`, exons are grouped by their
    let (mut merged, mut merged_rseqc) = if use_parallel {
        // --- Parallel chromosome processing ---
        //
        // Divide chromosomes into batches (one per thread) and process each
        // batch on its own thread with its own IndexedReader. The interval
        // index is shared read-only across all threads (COITree is Send+Sync).
        let num_chroms = tid_to_name.len();
        let num_workers = threads.min(num_chroms).max(1);

        // Distribute chromosomes using greedy bin-packing (largest-first) for
        // balanced load — assigns each chromosome to the worker with the smallest
        // current total, preventing large chromosomes from clustering on one thread.
        let batches = partition_chromosomes(&tid_to_len, num_workers);
        if log::log_enabled!(log::Level::Debug) {
            let worker_loads: Vec<u64> = batches
                .iter()
                .map(|b| b.iter().map(|&tid| tid_to_len[tid as usize]).sum())
                .collect();
            debug!(
                "Chromosome load distribution across {} workers: {:?}",
                num_workers, worker_loads
            );
        }

        // Configure rayon thread pool
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_workers)
            .build()
            .context("Failed to build rayon thread pool")?;

        // Global read counter for progress logging across threads
        let global_read_counter = AtomicU64::new(0);

        info!(
            "Processing {} chromosomes across {} threads",
            num_chroms, num_workers
        );

        // Calculate htslib decompression threads per worker
        let htslib_threads = if num_workers > 0 {
            (threads.saturating_sub(num_workers)) / num_workers
        } else {
            0
        };

        // Process chromosome batches in parallel
        let results: Vec<Result<(ChromResult, Option<RseqcAccumulators>)>> = pool.install(|| {
            batches
                .par_iter()
                .map(|batch| {
                    process_chromosome_batch(
                        bam_path,
                        batch,
                        &tid_to_name,
                        &index,
                        interner.len(),
                        stranded,
                        paired,
                        chrom_mapping,
                        chrom_prefix,
                        &global_read_counter,
                        reference,
                        rseqc_config,
                        rseqc_annotations,
                        htslib_threads,
                        qualimap_index,
                    )
                })
                .collect()
        });

        // Merge all chromosome results (both dupRadar and RSeQC)
        let mut merged = ChromResult::new(interner.len());
        let mut merged_rseqc: Option<RseqcAccumulators> =
            rseqc_config.map(|cfg| RseqcAccumulators::new(cfg, rseqc_annotations));
        for result in results {
            let (chrom_result, rseqc_result) = result?;
            merged.merge(chrom_result);
            if let (Some(ref mut merged_acc), Some(worker_acc)) = (&mut merged_rseqc, rseqc_result)
            {
                merged_acc.merge(worker_acc);
            }
        }
        (merged, merged_rseqc)
    } else {
        // --- Single-threaded fallback (no index or threads=1) ---
        //
        // Process all reads sequentially through a single pass over the alignment file.
        // This avoids the need for an index file.
        let global_read_counter = AtomicU64::new(0);
        let mut result = ChromResult::new(interner.len());
        let mut rseqc_accums: Option<RseqcAccumulators> =
            rseqc_config.map(|cfg| RseqcAccumulators::new(cfg, rseqc_annotations));
        let mut qualimap_accum: Option<crate::rna::qualimap::QualimapAccum> =
            qualimap_index.map(|_| crate::rna::qualimap::QualimapAccum::new(stranded));

        // Pre-compute resolved chromosome names for RSeQC tools
        // (apply chromosome prefix and mapping, same as parallel path)
        let tid_to_rseqc_chrom: Vec<String> = tid_to_name
            .iter()
            .map(|bam_name| {
                if let Some(mapped) = chrom_mapping.get(bam_name.as_str()) {
                    mapped.clone()
                } else if let Some(prefix) = chrom_prefix {
                    format!("{}{}", prefix, bam_name)
                } else {
                    bam_name.clone()
                }
            })
            .collect();
        let tid_to_rseqc_chrom_upper: Vec<String> = tid_to_rseqc_chrom
            .iter()
            .map(|s| s.to_uppercase())
            .collect();

        let mut bam = bam::Reader::from_path(bam_path)
            .with_context(|| format!("Failed to open alignment file: {}", bam_path))?;
        if let Some(ref_path) = reference {
            bam.set_reference(ref_path)
                .with_context(|| format!("Failed to set reference FASTA: {}", ref_path))?;
        }
        // Enable htslib multi-threaded decompression for the single-threaded path
        if threads > 1 {
            bam.set_threads(threads.saturating_sub(1))
                .context("Failed to set htslib decompression threads")?;
        }

        // Reusable buffers
        let mut aligned_blocks_buf: Vec<(u64, u64)> = Vec::new();
        let mut gene_hits: Vec<GeneIdx> = Vec::new();
        let mut mate_buffer: HashMap<MateBufferKey, MateInfo> = HashMap::new();
        let mut record = bam::Record::new();

        while let Some(read_result) = bam.read(&mut record) {
            read_result.context("Error reading alignment record")?;
            result.total_reads += 1;

            let prev = global_read_counter.fetch_add(1, Ordering::Relaxed);
            if (prev + 1).is_multiple_of(5_000_000) {
                debug!("Processed {} reads...", prev + 1);
            }

            let flags = record.flags();

            // --- RSeQC per-read dispatch (before counting filters) ---
            // bam_stat needs to see ALL records including QC-fail/dup/secondary
            if let (Some(ref mut accums), Some(annots), Some(cfg)) =
                (&mut rseqc_accums, rseqc_annotations, rseqc_config)
            {
                let tid = record.tid();
                let (chrom, chrom_upper) = if tid >= 0 && (tid as usize) < tid_to_rseqc_chrom.len()
                {
                    (
                        tid_to_rseqc_chrom[tid as usize].as_str(),
                        tid_to_rseqc_chrom_upper[tid as usize].as_str(),
                    )
                } else {
                    ("", "")
                };
                accums.process_read(&record, chrom, chrom_upper, annots, cfg);
            }

            // --- Qualimap per-read dispatch (before counting filters) ---
            // Qualimap has its own filtering (M-only CIGAR, NH>1, enclosure-based)
            if let (Some(ref mut qm_accum), Some(qm_index)) = (&mut qualimap_accum, qualimap_index)
            {
                let tid = record.tid();
                if tid >= 0 && (tid as usize) < tid_to_rseqc_chrom.len() {
                    let chrom = &tid_to_rseqc_chrom[tid as usize];
                    qm_accum.process_read(&record, chrom, qm_index);
                }
            }

            if flags & BAM_FUNMAP != 0 {
                result.fc_unmapped += 1;
                continue;
            }
            if flags & BAM_FSUPPLEMENTARY != 0 {
                continue;
            }
            if flags & BAM_FQCFAIL != 0 {
                continue;
            }
            if paired && flags & BAM_FPAIRED == 0 {
                continue;
            }

            result.total_mapped += 1;

            let is_dup = flags & BAM_FDUP != 0;
            if is_dup {
                result.total_dup += 1;
            }

            let is_multi = match record.aux(b"NH") {
                Ok(rust_htslib::bam::record::Aux::U8(nh)) => nh > 1,
                Ok(rust_htslib::bam::record::Aux::U16(nh)) => nh > 1,
                Ok(rust_htslib::bam::record::Aux::U32(nh)) => nh > 1,
                Ok(rust_htslib::bam::record::Aux::I8(nh)) => nh > 1,
                Ok(rust_htslib::bam::record::Aux::I16(nh)) => nh > 1,
                Ok(rust_htslib::bam::record::Aux::I32(nh)) => nh > 1,
                _ => false,
            };
            if is_multi {
                result.total_multi += 1;
            }

            let is_reverse = flags & BAM_FREVERSE != 0;
            let is_read1 = flags & BAM_FREAD1 != 0;

            let tid = record.tid();
            if tid < 0 || tid as usize >= tid_to_name.len() {
                continue;
            }
            let chrom = &tid_to_rseqc_chrom[tid as usize];

            gene_hits.clear();
            if let Some(chrom_idx) = index.get(chrom.as_str()) {
                cigar_to_aligned_blocks(
                    record.pos() as u64,
                    &record.cigar(),
                    &mut aligned_blocks_buf,
                );
                for &(block_start, block_end) in &aligned_blocks_buf {
                    chrom_idx.query(block_start as i32, (block_end - 1) as i32, |node| {
                        let meta = node.metadata;
                        if strand_matches(is_reverse, is_read1, paired, meta.strand, stranded) {
                            gene_hits.push(meta.gene_idx);
                        }
                    });
                }
                gene_hits.sort_unstable();
                gene_hits.dedup();
            }

            // --- Per-read featureCounts counting (independent of mate pairing) ---
            classify_read_fc(is_multi, &gene_hits, &mut result);

            if !paired {
                result.n_multi_dup += 1;
                result.n_unique_dup += 1;
                if !is_dup {
                    result.n_multi_nodup += 1;
                    result.n_unique_nodup += 1;
                }
                result.total_fragments += 1;

                if gene_hits.is_empty() {
                    result.stat_no_features += 1;
                } else if gene_hits.len() > 1 {
                    result.stat_ambiguous += 1;
                } else if assign_fragment_to_gene(
                    &gene_hits,
                    &mut result.gene_counts,
                    is_dup,
                    is_multi,
                ) {
                    result.stat_assigned += 1;
                }
                continue;
            }

            // Paired-end mate buffering (same logic as process_chromosome_batch)
            let qname_hash = hash_qname(record.qname());
            let own_pos = record.pos();
            let own_tid = record.tid();
            let mate_pos_val = record.mpos();
            let mate_tid_val = record.mtid();

            let hi_tag: i32 = match record.aux(b"HI") {
                Ok(rust_htslib::bam::record::Aux::U8(v)) => v as i32,
                Ok(rust_htslib::bam::record::Aux::U16(v)) => v as i32,
                Ok(rust_htslib::bam::record::Aux::U32(v)) => v as i32,
                Ok(rust_htslib::bam::record::Aux::I8(v)) => v as i32,
                Ok(rust_htslib::bam::record::Aux::I16(v)) => v as i32,
                Ok(rust_htslib::bam::record::Aux::I32(v)) => v,
                _ => -1,
            };

            let buffer_key: MateBufferKey = if is_read1 {
                (
                    qname_hash,
                    own_tid,
                    own_pos,
                    mate_tid_val,
                    mate_pos_val,
                    hi_tag,
                )
            } else {
                (
                    qname_hash,
                    mate_tid_val,
                    mate_pos_val,
                    own_tid,
                    own_pos,
                    hi_tag,
                )
            };

            if let Some(mate_info) = mate_buffer.remove(&buffer_key) {
                result.total_fragments += 1;
                let frag_is_dup = if is_read1 { is_dup } else { mate_info.is_dup };
                let frag_is_multi = if is_read1 {
                    is_multi
                } else {
                    mate_info.is_multi
                };

                result.n_multi_dup += 1;
                result.n_unique_dup += 1;
                if !frag_is_dup {
                    result.n_multi_nodup += 1;
                    result.n_unique_nodup += 1;
                }

                // featureCounts union-with-both-end-preference voting
                // (same logic as parallel path — see merge_gene_hits)
                let combined_genes: Vec<GeneIdx> =
                    merge_gene_hits(&mate_info.gene_hits, &gene_hits);

                if combined_genes.is_empty() {
                    result.stat_no_features += 1;
                } else if combined_genes.len() > 1 {
                    result.stat_ambiguous += 1;
                } else if assign_fragment_to_gene(
                    &combined_genes,
                    &mut result.gene_counts,
                    frag_is_dup,
                    frag_is_multi,
                ) {
                    result.stat_assigned += 1;
                }
            } else {
                // Take ownership of gene_hits instead of cloning (it's cleared at loop start)
                mate_buffer.insert(
                    buffer_key,
                    MateInfo {
                        gene_hits: std::mem::take(&mut gene_hits),
                        is_dup,
                        is_multi,
                    },
                );
            }
        }

        result.unmatched_mates = mate_buffer;
        result.qualimap = qualimap_accum;
        (result, rseqc_accums)
    };

    // --- Sweep unmapped reads (parallel path only) ---
    //
    // The parallel path uses IndexedReader.fetch(tid) per chromosome, which
    // never visits reads in the unmapped segment (tid=-1) at the end of the
    // BAM file. We sweep them here so that bam_stat totals, fc_unmapped,
    // and RSeQC accumulators see every record.
    if use_parallel {
        let mut bam = bam::IndexedReader::from_path(bam_path).with_context(|| {
            format!(
                "Failed to open indexed BAM for unmapped sweep: {}",
                bam_path
            )
        })?;
        if let Some(ref_path) = reference {
            bam.set_reference(ref_path)
                .with_context(|| format!("Failed to set reference FASTA: {}", ref_path))?;
        }
        bam.fetch(FetchDefinition::Unmapped)
            .context("Failed to seek to unmapped segment")?;

        let mut record = bam::Record::new();
        let mut unmapped_count: u64 = 0;

        while let Some(read_result) = bam.read(&mut record) {
            read_result.context("Error reading unmapped record")?;
            unmapped_count += 1;
            merged.total_reads += 1;

            // RSeQC dispatch — empty chromosome strings for unplaced reads
            if let (Some(ref mut accums), Some(annots), Some(cfg)) =
                (&mut merged_rseqc, rseqc_annotations, rseqc_config)
            {
                accums.process_read(&record, "", "", annots, cfg);
            }

            // Qualimap dispatch — qualimap will skip internally since tid < 0
            if let Some(ref mut qm) = merged.qualimap {
                // tid < 0 for unmapped reads, so qualimap cannot resolve a
                // chromosome name. Skip rather than panic on out-of-bounds.
                let tid = record.tid();
                if tid >= 0 {
                    // Shouldn't happen in the unmapped segment, but be safe
                    if let Some(qm_index) = qualimap_index {
                        if (tid as usize) < tid_to_name.len() {
                            qm.process_read(&record, &tid_to_name[tid as usize], qm_index);
                        }
                    }
                }
            }

            let flags = record.flags();
            if flags & BAM_FUNMAP != 0 {
                merged.fc_unmapped += 1;
            }
        }

        if unmapped_count > 0 {
            info!(
                "Parallel unmapped sweep: processed {} unmapped reads",
                unmapped_count
            );
        }
    }

    // --- Cross-chromosome mate reconciliation ---
    //
    // In parallel mode, mates on different chromosomes end up in different
    // workers' mate buffers. Reconcile them here: try to match each unmatched
    // mate against the others. Any that still don't match are treated as
    // singletons (same as the sequential path).
    if paired && !merged.unmatched_mates.is_empty() {
        // Drain all unmatched mates and try to pair them
        let unmatched: Vec<(MateBufferKey, MateInfo)> = merged.unmatched_mates.drain().collect();
        let mut still_unmatched: HashMap<MateBufferKey, MateInfo> = HashMap::new();

        for (key, info) in unmatched {
            if let Some(mate_info) = still_unmatched.remove(&key) {
                // Found a cross-chromosome mate pair!
                merged.total_fragments += 1;

                // Determine which is read1 — we don't track is_read1 in MateInfo,
                // but both mates compute the same key. Use read1's dup/multi status.
                // Since we can't distinguish which is read1 from MateInfo alone,
                // use the first mate seen (info) as the "read1" for dup/multi.
                let frag_is_dup = info.is_dup;
                let frag_is_multi = info.is_multi;

                merged.n_multi_dup += 1;
                merged.n_unique_dup += 1;
                if !frag_is_dup {
                    merged.n_multi_nodup += 1;
                    merged.n_unique_nodup += 1;
                }

                // Combine gene hits with featureCounts scoring
                let combined_genes: Vec<GeneIdx> =
                    merge_gene_hits(&mate_info.gene_hits, &info.gene_hits);

                if combined_genes.is_empty() {
                    merged.stat_no_features += 1;
                } else if combined_genes.len() > 1 {
                    merged.stat_ambiguous += 1;
                } else if assign_fragment_to_gene(
                    &combined_genes,
                    &mut merged.gene_counts,
                    frag_is_dup,
                    frag_is_multi,
                ) {
                    merged.stat_assigned += 1;
                    // Coverage recording skipped — no aligned blocks available for
                    // cross-chromosome mate reconciliation.
                }
            } else {
                still_unmatched.insert(key, info);
            }
        }

        // Handle remaining singletons (mates whose partner was filtered out)
        for (_key, mate_info) in still_unmatched.drain() {
            merged.total_fragments += 1;

            merged.n_multi_dup += 1;
            merged.n_unique_dup += 1;
            if !mate_info.is_dup {
                merged.n_multi_nodup += 1;
                merged.n_unique_nodup += 1;
            }

            if mate_info.gene_hits.is_empty() {
                merged.stat_no_features += 1;
            } else if mate_info.gene_hits.len() > 1 {
                merged.stat_ambiguous += 1;
            } else if assign_fragment_to_gene(
                &mate_info.gene_hits,
                &mut merged.gene_counts,
                mate_info.is_dup,
                mate_info.is_multi,
            ) {
                merged.stat_assigned += 1;
                // Coverage recording skipped here — no aligned blocks available for
                // cross-chromosome singleton mates (rare edge case).
            }
        }
    }

    info!(
        "Read {} total reads, {} mapped, {} fragments, {} duplicates, {} multimappers",
        merged.total_reads,
        merged.total_mapped,
        merged.total_fragments,
        merged.total_dup,
        merged.total_multi
    );
    info!(
        "Assignment stats: {} assigned, {} ambiguous, {} no_features (total fragments: {})",
        merged.stat_assigned,
        merged.stat_ambiguous,
        merged.stat_no_features,
        merged.total_fragments
    );

    // Verify that at least some reads were flagged as duplicates.
    // Even if the @PG header check passed (or was skipped), it's possible that
    // duplicates were not actually marked. Warn the user in that case.
    if merged.total_dup == 0 && merged.total_mapped > 0 && !skip_dup_check {
        anyhow::bail!(
            "No duplicate-flagged reads found among {} mapped reads in '{}'.\n\
             \n\
             Although the BAM header suggests a duplicate-marking tool was run,\n\
             no reads have the duplicate flag (0x400) set. This likely means\n\
             duplicates were removed rather than marked, or the tool did not\n\
             flag any duplicates.\n\
             \n\
             RustQC requires duplicates to be marked (flagged) but NOT removed.\n\
             Please re-run your duplicate-marking tool without removing duplicates.\n\
             \n\
             If this is expected, use --skip-dup-check to bypass this check.",
            merged.total_mapped,
            bam_path
        );
    }

    // Detect chromosome name mismatch
    let genes_with_reads = merged
        .gene_counts
        .iter()
        .filter(|c| c.all_multi > 0)
        .count();
    if merged.total_mapped > 0 && genes_with_reads == 0 {
        let bam_chroms: Vec<&str> = tid_to_name.iter().take(5).map(|s| s.as_str()).collect();
        let gtf_chroms: Vec<&str> = index.keys().take(5).map(|s| s.as_str()).collect();
        anyhow::bail!(
            "Chromosome name mismatch: no reads could be assigned to any gene.\n\
             \n\
             Alignment chromosomes (first 5): {}\n\
             GTF chromosomes (first 5): {}\n\
             \n\
             The alignment and GTF files appear to use different chromosome naming conventions.\n\
             To fix this, create a YAML config file with a chromosome_mapping section and pass it via --config:\n\
             \n\
             Example config.yaml:\n\
             \n\
             chromosome_mapping:\n\
             {}",
            bam_chroms.join(", "),
            gtf_chroms.join(", "),
            bam_chroms
                .iter()
                .zip(gtf_chroms.iter())
                .map(|(b, g)| format!("  {}: {}", g, b))
                .collect::<Vec<_>>()
                .join("\n")
        );
    }

    // Convert flat Vec<GeneCounts> back to IndexMap<String, GeneCounts> for output
    let gene_counts_map: IndexMap<String, GeneCounts> = (0..interner.len())
        .map(|i| {
            let name = interner.name(i as GeneIdx).to_string();
            let counts = std::mem::take(&mut merged.gene_counts[i]);
            (name, counts)
        })
        .collect();

    Ok(CountResult {
        gene_counts: gene_counts_map,
        total_reads_multi_dup: merged.n_multi_dup,
        total_reads_unique_dup: merged.n_unique_dup,
        stat_total_reads: merged.total_reads,
        stat_assigned: merged.stat_assigned,
        stat_ambiguous: merged.stat_ambiguous,
        stat_no_features: merged.stat_no_features,
        stat_total_fragments: merged.total_fragments,
        fc_assigned: merged.fc_assigned,
        fc_ambiguous: merged.fc_ambiguous,
        fc_no_features: merged.fc_no_features,
        fc_multimapping: merged.fc_multimapping,
        fc_unmapped: merged.fc_unmapped,
        rseqc: merged_rseqc,
        qualimap: match (merged.qualimap, qualimap_index) {
            (Some(mut a), Some(idx)) => {
                a.flush_unpaired(idx);
                Some(a.into_result())
            }
            _ => None,
        },
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strand_matching_unstranded() {
        // Unstranded: everything matches
        assert!(strand_matches(false, true, false, '+', 0));
        assert!(strand_matches(true, true, false, '+', 0));
        assert!(strand_matches(false, true, false, '-', 0));
        assert!(strand_matches(true, true, false, '-', 0));
    }

    #[test]
    fn test_strand_matching_forward_single() {
        // Forward stranded, single-end
        // Read on + strand should match + gene
        assert!(strand_matches(false, true, false, '+', 1));
        // Read on - strand should match - gene
        assert!(strand_matches(true, true, false, '-', 1));
        // Read on + strand should NOT match - gene
        assert!(!strand_matches(false, true, false, '-', 1));
        // Read on - strand should NOT match + gene
        assert!(!strand_matches(true, true, false, '+', 1));
    }

    #[test]
    fn test_strand_matching_reverse_single() {
        // Reverse stranded, single-end
        // Read on + strand should match - gene (reversed)
        assert!(strand_matches(false, true, false, '-', 2));
        // Read on - strand should match + gene (reversed)
        assert!(strand_matches(true, true, false, '+', 2));
    }

    #[test]
    fn test_strand_matching_forward_paired() {
        // Forward stranded, paired-end
        // Read1 on + strand: effective + -> matches + gene
        assert!(strand_matches(false, true, true, '+', 1));
        // Read2 on + strand: effective - -> matches - gene
        assert!(strand_matches(false, false, true, '-', 1));
        // Read1 on - strand: effective - -> matches - gene
        assert!(strand_matches(true, true, true, '-', 1));
        // Read2 on - strand: effective + -> matches + gene
        assert!(strand_matches(true, false, true, '+', 1));
    }

    // ---------------------------------------------------------------
    // Duplicate marking validation tests
    // ---------------------------------------------------------------

    #[test]
    fn test_extract_header_tag() {
        let line = "@PG\tID:MarkDuplicates\tPN:MarkDuplicates\tVN:2.27.4\tCL:picard MarkDuplicates I=in.bam O=out.bam";
        assert_eq!(extract_header_tag(line, "ID"), Some("MarkDuplicates"));
        assert_eq!(extract_header_tag(line, "PN"), Some("MarkDuplicates"));
        assert_eq!(extract_header_tag(line, "VN"), Some("2.27.4"));
        assert_eq!(
            extract_header_tag(line, "CL"),
            Some("picard MarkDuplicates I=in.bam O=out.bam")
        );
        assert_eq!(extract_header_tag(line, "XX"), None);
    }

    #[test]
    fn test_dup_check_picard_markduplicates() {
        let header = "@HD\tVN:1.6\tSO:coordinate\n\
                       @PG\tID:MarkDuplicates\tPN:MarkDuplicates\tVN:2.27.4";
        assert!(header_has_dup_marker(header));
    }

    #[test]
    fn test_dup_check_samblaster() {
        let header = "@HD\tVN:1.6\n\
                       @PG\tID:samblaster\tPN:samblaster\tVN:0.1.26";
        assert!(header_has_dup_marker(header));
    }

    #[test]
    fn test_dup_check_sambamba_markdup() {
        let header = "@HD\tVN:1.6\n\
                       @PG\tID:sambamba_markdup\tPN:sambamba markdup";
        assert!(header_has_dup_marker(header));
    }

    #[test]
    fn test_dup_check_biobambam() {
        let header = "@HD\tVN:1.6\n\
                       @PG\tID:biobambam2\tPN:biobambam";
        assert!(header_has_dup_marker(header));
    }

    #[test]
    fn test_dup_check_case_insensitive() {
        // MarkDuplicates with different casing
        let header = "@PG\tID:MARKDUPLICATES\tPN:markduplicates";
        assert!(header_has_dup_marker(header));
    }

    #[test]
    fn test_dup_check_no_dup_marker() {
        // Header with only alignment tools, no dup marker
        let header = "@HD\tVN:1.6\tSO:coordinate\n\
                       @PG\tID:bwa\tPN:bwa\tVN:0.7.17\n\
                       @PG\tID:samtools\tPN:samtools\tVN:1.17";
        assert!(!header_has_dup_marker(header));
    }

    #[test]
    fn test_dup_check_empty_header() {
        assert!(!header_has_dup_marker(""));
        assert!(!header_has_dup_marker("@HD\tVN:1.6\tSO:coordinate"));
    }

    #[test]
    fn test_dup_check_picard_sortsam_no_false_positive() {
        // Picard SortSam should NOT match — only the CL field mentions "picard"
        let header = "@PG\tID:SortSam\tPN:SortSam\tCL:picard SortSam I=in.bam O=out.bam";
        assert!(!header_has_dup_marker(header));
    }

    #[test]
    fn test_dup_check_picard_collect_metrics_no_false_positive() {
        // Another Picard tool that should not match
        let header =
            "@PG\tID:CollectInsertSizeMetrics\tPN:CollectInsertSizeMetrics\tCL:picard CollectInsertSizeMetrics";
        assert!(!header_has_dup_marker(header));
    }

    #[test]
    fn test_dup_check_sambamba_sort_no_false_positive() {
        // sambamba sort should NOT match
        let header = "@PG\tID:sambamba\tPN:sambamba sort";
        assert!(!header_has_dup_marker(header));
    }

    #[test]
    fn test_dup_check_multiple_pg_lines() {
        // MarkDuplicates appears as a later @PG entry in the chain
        let header = "@HD\tVN:1.6\tSO:coordinate\n\
                       @PG\tID:bwa\tPN:bwa\tVN:0.7.17\n\
                       @PG\tID:samtools\tPN:samtools\tVN:1.17\n\
                       @PG\tID:MarkDuplicates\tPN:MarkDuplicates\tPP:samtools";
        assert!(header_has_dup_marker(header));
    }

    #[test]
    fn test_dup_check_dup_marker_only_in_cl_no_match() {
        let header_text = "@PG\tID:STAR\tPN:STAR\tVN:2.7.10a\tCL:STAR --readFilesIn sample.fastq --outSAMtype BAM SortedByCoordinate\n";
        assert!(!header_has_dup_marker(header_text));
    }

    // --- Per-read featureCounts classification tests ---

    /// Helper to create a ChromResult for testing classify_read_fc
    fn make_test_chrom_result(num_genes: usize) -> ChromResult {
        ChromResult::new(num_genes)
    }

    #[test]
    fn test_fc_classify_multimapped_read() {
        let mut result = make_test_chrom_result(3);
        let gene_hits: Vec<GeneIdx> = vec![0]; // has a hit, but is_multi=true

        classify_read_fc(true, &gene_hits, &mut result);

        assert_eq!(
            result.fc_multimapping, 1,
            "Multi-mapped read should be counted as fc_multimapping"
        );
        assert_eq!(result.fc_assigned, 0);
        assert_eq!(result.fc_ambiguous, 0);
        assert_eq!(result.fc_no_features, 0);
        assert_eq!(
            result.gene_counts[0].fc_reads, 0,
            "Multi-mapped read should not credit any gene"
        );
    }

    #[test]
    fn test_fc_classify_no_gene_hits() {
        let mut result = make_test_chrom_result(3);
        let gene_hits: Vec<GeneIdx> = vec![];

        classify_read_fc(false, &gene_hits, &mut result);

        assert_eq!(
            result.fc_no_features, 1,
            "Read with no hits should be fc_no_features"
        );
        assert_eq!(result.fc_assigned, 0);
        assert_eq!(result.fc_ambiguous, 0);
        assert_eq!(result.fc_multimapping, 0);
    }

    #[test]
    fn test_fc_classify_single_gene_hit() {
        let mut result = make_test_chrom_result(3);
        let gene_hits: Vec<GeneIdx> = vec![1];

        classify_read_fc(false, &gene_hits, &mut result);

        assert_eq!(
            result.fc_assigned, 1,
            "Single gene hit should be fc_assigned"
        );
        assert_eq!(
            result.gene_counts[1].fc_reads, 1,
            "Single gene hit should credit that gene"
        );
        assert_eq!(result.gene_counts[0].fc_reads, 0);
        assert_eq!(result.gene_counts[2].fc_reads, 0);
        assert_eq!(result.fc_ambiguous, 0);
        assert_eq!(result.fc_no_features, 0);
    }

    #[test]
    fn test_fc_classify_multi_gene_hits_ambiguous() {
        // Two genes hit → always Ambiguous (matching featureCounts -g gene_id)
        let mut result = make_test_chrom_result(4);
        let gene_hits: Vec<GeneIdx> = vec![0, 1];

        classify_read_fc(false, &gene_hits, &mut result);

        assert_eq!(
            result.fc_ambiguous, 1,
            "Multi-gene hits should be fc_ambiguous"
        );
        assert_eq!(result.fc_assigned, 0);
        assert_eq!(
            result.gene_counts[0].fc_reads, 0,
            "Ambiguous read should not credit any gene"
        );
        assert_eq!(result.gene_counts[1].fc_reads, 0);
    }

    #[test]
    fn test_fc_classify_three_gene_hits_ambiguous() {
        // Three genes hit → Ambiguous
        let mut result = make_test_chrom_result(5);
        let gene_hits: Vec<GeneIdx> = vec![0, 1, 2];

        classify_read_fc(false, &gene_hits, &mut result);

        assert_eq!(
            result.fc_ambiguous, 1,
            "Three-gene hits should be fc_ambiguous"
        );
        assert_eq!(result.fc_assigned, 0);
        assert_eq!(result.gene_counts[0].fc_reads, 0);
        assert_eq!(result.gene_counts[1].fc_reads, 0);
        assert_eq!(result.gene_counts[2].fc_reads, 0);
    }

    #[test]
    fn test_fc_classify_cumulative_counts() {
        // Verify that calling classify_read_fc multiple times accumulates counts
        let mut result = make_test_chrom_result(3);

        // First call: single hit to gene 0
        classify_read_fc(false, &[0], &mut result);
        // Second call: single hit to gene 0 again
        classify_read_fc(false, &[0], &mut result);
        // Third call: multi-mapped
        classify_read_fc(true, &[0], &mut result);
        // Fourth call: no features
        classify_read_fc(false, &[], &mut result);
        // Fifth call: ambiguous (multiple genes)
        classify_read_fc(false, &[0, 1], &mut result);

        assert_eq!(result.fc_assigned, 2);
        assert_eq!(result.fc_multimapping, 1);
        assert_eq!(result.fc_no_features, 1);
        assert_eq!(result.fc_ambiguous, 1);
        assert_eq!(result.gene_counts[0].fc_reads, 2);
        assert_eq!(result.gene_counts[1].fc_reads, 0);
    }

    // --- Chromosome partitioning tests ---

    #[test]
    fn test_partition_chromosomes_balanced() {
        // Simulate human-like chromosome sizes (large variation)
        let lengths = vec![
            250, 243, 198, 191, 182, 171, 159, 146, 138, 133, 135, 130, 57,
        ];
        let batches = partition_chromosomes(&lengths, 4);

        assert_eq!(batches.len(), 4);

        // Every chromosome should appear exactly once
        let mut all_tids: Vec<u32> = batches.iter().flatten().copied().collect();
        all_tids.sort();
        assert_eq!(all_tids, (0..13).collect::<Vec<u32>>());

        // Check that loads are reasonably balanced (no worker > 1.3× average)
        let total: u64 = lengths.iter().sum();
        let avg = total as f64 / 4.0;
        for batch in &batches {
            let load: u64 = batch.iter().map(|&tid| lengths[tid as usize]).sum();
            assert!(
                (load as f64) < avg * 1.3,
                "Worker load {} exceeds 1.3× average {:.0}",
                load,
                avg
            );
        }
    }

    #[test]
    fn test_partition_chromosomes_single_worker() {
        let lengths = vec![100, 200, 300];
        let batches = partition_chromosomes(&lengths, 1);
        assert_eq!(batches.len(), 1);
        assert_eq!(batches[0].len(), 3);
    }

    #[test]
    fn test_partition_chromosomes_more_workers_than_chroms() {
        let lengths = vec![100, 200];
        let batches = partition_chromosomes(&lengths, 4);
        assert_eq!(batches.len(), 4);
        let non_empty: usize = batches.iter().filter(|b| !b.is_empty()).count();
        assert_eq!(non_empty, 2);
    }
}
