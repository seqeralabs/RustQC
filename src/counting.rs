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
use anyhow::{Context, Result};
use coitrees::{COITree, Interval, IntervalTree};
use indexmap::IndexMap;
use log::{debug, info, warn};
use rayon::prelude::*;
use rust_htslib::bam::{self, Read as BamRead};
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

/// Flag indicating the read is a PCR or optical duplicate (0x400).
const BAM_FDUP: u16 = 0x400;
/// Flag indicating the read is unmapped (0x4).
const BAM_FUNMAP: u16 = 0x4;
/// Flag indicating the read failed quality checks (0x200).
const BAM_FQCFAIL: u16 = 0x200;
/// Flag indicating a secondary alignment (0x100).
#[allow(dead_code)]
const BAM_FSECONDARY: u16 = 0x100;
/// Flag indicating a supplementary alignment (0x800).
const BAM_FSUPPLEMENTARY: u16 = 0x800;
/// Flag indicating the read is paired (0x1).
const BAM_FPAIRED: u16 = 0x1;
/// Flag indicating it is the first read in a pair (0x40).
const BAM_FREAD1: u16 = 0x40;
/// Flag indicating the read is mapped in a proper pair (0x2).
#[allow(dead_code)]
const BAM_FPROPER_PAIR: u16 = 0x2;
/// Flag indicating the mate is unmapped (0x8).
#[allow(dead_code)]
const BAM_FMUNMAP: u16 = 0x8;
/// Flag indicating the read is reverse-complemented (0x10).
const BAM_FREVERSE: u16 = 0x10;

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
/// featureCounts-compatible output files.
#[derive(Debug)]
pub struct CountResult {
    /// Per-gene counts indexed by gene_id
    pub gene_counts: IndexMap<String, GeneCounts>,
    /// Total mapped reads (including multimappers, including duplicates)
    pub total_reads_multi_dup: u64,
    /// Total mapped reads (including multimappers, excluding duplicates)
    #[allow(dead_code)]
    pub total_reads_multi_nodup: u64,
    /// Total mapped reads (unique only, including duplicates)
    pub total_reads_unique_dup: u64,
    /// Total mapped reads (unique only, excluding duplicates)
    #[allow(dead_code)]
    pub total_reads_unique_nodup: u64,

    // --- featureCounts summary statistics ---
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
    /// Total duplicate reads
    #[allow(dead_code)]
    pub stat_total_dup: u64,
    /// Total multimapping reads
    #[allow(dead_code)]
    pub stat_total_multi: u64,

    // --- BAM stat counters (RSeQC bam_stat.py equivalent) ---
    /// Unmapped reads (BAM flag 0x4)
    pub stat_unmapped: u64,
    /// QC-failed reads (BAM flag 0x200)
    pub stat_qc_failed: u64,
    /// Secondary alignments (BAM flag 0x100)
    pub stat_secondary: u64,
    /// Supplementary alignments (BAM flag 0x800)
    pub stat_supplementary: u64,
    /// Non-primary hits (secondary + supplementary)
    pub stat_non_primary_hits: u64,
    /// Reads mapped in proper pairs (BAM flag 0x2)
    pub stat_proper_pairs: u64,
    /// Read-1 reads among mapped (BAM flag 0x40)
    pub stat_read1: u64,
    /// Read-2 reads among mapped (paired but not read-1)
    pub stat_read2: u64,
    /// Mapped reads on the forward (+) strand
    pub stat_plus_strand: u64,
    /// Mapped reads on the reverse (-) strand
    pub stat_minus_strand: u64,
    /// Singletons: mapped reads whose mate is unmapped
    pub stat_singletons: u64,
    /// Reads whose mate maps to a different chromosome
    pub stat_mate_mapped_diff_chr: u64,
    /// Spliced reads (CIGAR contains 'N' / RefSkip operation)
    pub stat_splice_reads: u64,
    /// Non-spliced mapped reads (no 'N' in CIGAR)
    pub stat_non_splice_reads: u64,
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
type MateBufferKey = (Vec<u8>, i32, i64, i32, i64, i32);

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

    // BAM stat counters
    stat_unmapped: u64,
    stat_qc_failed: u64,
    stat_secondary: u64,
    stat_supplementary: u64,
    stat_proper_pairs: u64,
    stat_read1: u64,
    stat_read2: u64,
    stat_plus_strand: u64,
    stat_minus_strand: u64,
    stat_singletons: u64,
    stat_mate_mapped_diff_chr: u64,
    stat_splice_reads: u64,
    stat_non_splice_reads: u64,
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
            stat_unmapped: 0,
            stat_qc_failed: 0,
            stat_secondary: 0,
            stat_supplementary: 0,
            stat_proper_pairs: 0,
            stat_read1: 0,
            stat_read2: 0,
            stat_plus_strand: 0,
            stat_minus_strand: 0,
            stat_singletons: 0,
            stat_mate_mapped_diff_chr: 0,
            stat_splice_reads: 0,
            stat_non_splice_reads: 0,
        }
    }

    /// Merge another ChromResult into this one (additive).
    fn merge(&mut self, other: ChromResult) {
        for (i, counts) in other.gene_counts.into_iter().enumerate() {
            self.gene_counts[i].all_multi += counts.all_multi;
            self.gene_counts[i].nodup_multi += counts.nodup_multi;
            self.gene_counts[i].all_unique += counts.all_unique;
            self.gene_counts[i].nodup_unique += counts.nodup_unique;
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
        self.stat_unmapped += other.stat_unmapped;
        self.stat_qc_failed += other.stat_qc_failed;
        self.stat_secondary += other.stat_secondary;
        self.stat_supplementary += other.stat_supplementary;
        self.stat_proper_pairs += other.stat_proper_pairs;
        self.stat_read1 += other.stat_read1;
        self.stat_read2 += other.stat_read2;
        self.stat_plus_strand += other.stat_plus_strand;
        self.stat_minus_strand += other.stat_minus_strand;
        self.stat_singletons += other.stat_singletons;
        self.stat_mate_mapped_diff_chr += other.stat_mate_mapped_diff_chr;
        self.stat_splice_reads += other.stat_splice_reads;
        self.stat_non_splice_reads += other.stat_non_splice_reads;
    }
}

/// Process a batch of chromosomes from an alignment file, counting reads against gene annotations.
///
/// Opens its own indexed reader and seeks to each chromosome in the batch.
/// Returns accumulated results for all chromosomes in the batch.
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
) -> Result<ChromResult> {
    let mut result = ChromResult::new(num_genes);

    // Open an indexed reader for this thread (supports BAM with .bai/.csi and CRAM with .crai)
    let mut bam = bam::IndexedReader::from_path(bam_path)
        .with_context(|| format!("Failed to open indexed alignment file: {}", bam_path))?;
    if let Some(ref_path) = reference {
        bam.set_reference(ref_path)
            .with_context(|| format!("Failed to set reference FASTA: {}", ref_path))?;
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

            let flags = record.flags();

            // Count BAM stat categories before skip filters
            if flags & BAM_FUNMAP != 0 {
                result.stat_unmapped += 1;
            }
            if flags & BAM_FQCFAIL != 0 {
                result.stat_qc_failed += 1;
            }
            if flags & BAM_FSECONDARY != 0 {
                result.stat_secondary += 1;
            }
            if flags & BAM_FSUPPLEMENTARY != 0 {
                result.stat_supplementary += 1;
            }

            // Skip unmapped reads
            if flags & BAM_FUNMAP != 0 {
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

            // BAM stat counters for mapped reads
            if is_reverse {
                result.stat_minus_strand += 1;
            } else {
                result.stat_plus_strand += 1;
            }
            if flags & BAM_FPAIRED != 0 {
                if is_read1 {
                    result.stat_read1 += 1;
                } else {
                    result.stat_read2 += 1;
                }
                if flags & BAM_FPROPER_PAIR != 0 {
                    result.stat_proper_pairs += 1;
                }
                if flags & BAM_FMUNMAP != 0 {
                    result.stat_singletons += 1;
                } else if record.mtid() != record.tid() {
                    result.stat_mate_mapped_diff_chr += 1;
                }
            }

            // Splice detection: check CIGAR for RefSkip ('N') operations
            {
                use rust_htslib::bam::record::Cigar;
                let has_splice = record
                    .cigar()
                    .iter()
                    .any(|op| matches!(op, Cigar::RefSkip(_)));
                if has_splice {
                    result.stat_splice_reads += 1;
                } else {
                    result.stat_non_splice_reads += 1;
                }
            }

            // Get the chromosome name
            let rec_tid = record.tid();
            if rec_tid < 0 || rec_tid as usize >= tid_to_name.len() {
                continue;
            }
            let bam_chrom = &tid_to_name[rec_tid as usize];
            // Apply chromosome name prefix and/or mapping (alignment name -> GTF name)
            let prefixed_chrom;
            let chrom = if let Some(mapped) = chrom_mapping.get(bam_chrom.as_str()) {
                mapped.as_str()
            } else if let Some(prefix) = chrom_prefix {
                prefixed_chrom = format!("{}{}", prefix, bam_chrom);
                prefixed_chrom.as_str()
            } else {
                bam_chrom.as_str()
            };

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
            let read_name = record.qname().to_vec();
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
                (
                    read_name.clone(),
                    own_tid,
                    own_pos,
                    mate_tid,
                    mate_pos_val,
                    hi_tag,
                )
            } else {
                (
                    read_name.clone(),
                    mate_tid,
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

                let combined_genes: Vec<GeneIdx> =
                    if mate_info.gene_hits.is_empty() && gene_hits.is_empty() {
                        Vec::new()
                    } else {
                        let mut scores: HashMap<GeneIdx, u8> = HashMap::new();
                        for &g in &mate_info.gene_hits {
                            *scores.entry(g).or_insert(0) += 1;
                        }
                        for &g in &gene_hits {
                            *scores.entry(g).or_insert(0) += 1;
                        }
                        let max_score = scores.values().copied().max().unwrap_or(0);
                        let mut best: Vec<GeneIdx> = scores
                            .iter()
                            .filter(|(_, &s)| s == max_score)
                            .map(|(&g, _)| g)
                            .collect();
                        best.sort_unstable();
                        best
                    };

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
                mate_buffer.insert(
                    buffer_key,
                    MateInfo {
                        gene_hits: gene_hits.clone(),
                        is_dup,
                        is_multi,
                    },
                );
            }
        }
    }

    // Move unmatched mates into the result for cross-chromosome reconciliation
    result.unmatched_mates = mate_buffer;
    Ok(result)
}

/// Count reads from an alignment file (SAM/BAM/CRAM) and assign them to genes.
///
/// Performs four simultaneous counting modes matching dupRadar's approach:
/// - With/without multimappers
/// - With/without PCR duplicates
///
/// When `threads > 1`, processing is parallelized by chromosome: the alignment
/// index is used to divide chromosomes among threads, each opening its own reader
/// and processing independently. Per-chromosome results are merged, and any
/// unmatched paired-end mates (from cross-chromosome pairs) are reconciled in a
/// final pass.
///
/// For paired-end data, this function buffers mates by a composite key matching
/// featureCounts' `SAM_pairer_get_read_full_name()` (read name, R1/R2 refIDs,
/// R1/R2 positions, HI tag). Gene assignment uses featureCounts' scoring strategy:
/// genes overlapped by both mates score higher than genes overlapped by only one.
///
/// # Arguments
/// * `bam_path` - Path to the duplicate-marked alignment file (BAM/CRAM must have
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
) -> Result<CountResult> {
    // Build gene ID interner for allocation-free lookups in the hot loop
    let interner = GeneIdInterner::from_genes(genes);

    // Build spatial index (uses interned gene indices)
    let index = build_index(genes, &interner);

    // Get chromosome names from header using a temporary reader,
    // and verify that duplicates have been marked in the BAM file.
    let tid_to_name: Vec<String> = {
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

        (0..header.target_count())
            .map(|tid| String::from_utf8_lossy(header.tid2name(tid)).to_string())
            .collect()
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

    let mut merged = if use_parallel {
        // --- Parallel chromosome processing ---
        //
        // Divide chromosomes into batches (one per thread) and process each
        // batch on its own thread with its own IndexedReader. The interval
        // index is shared read-only across all threads (COITree is Send+Sync).
        let num_chroms = tid_to_name.len();
        let num_workers = threads.min(num_chroms).max(1);

        // Distribute chromosomes round-robin across workers for balanced load
        // (chromosomes vary widely in size, so round-robin spreads large ones)
        let mut batches: Vec<Vec<u32>> = vec![Vec::new(); num_workers];
        for (i, _name) in tid_to_name.iter().enumerate() {
            batches[i % num_workers].push(i as u32);
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

        // Process chromosome batches in parallel
        let results: Vec<Result<ChromResult>> = pool.install(|| {
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
                    )
                })
                .collect()
        });

        // Merge all chromosome results
        let mut merged = ChromResult::new(interner.len());
        for result in results {
            merged.merge(result?);
        }
        merged
    } else {
        // --- Single-threaded fallback (no index or threads=1) ---
        //
        // Process all reads sequentially through a single pass over the alignment file.
        // This avoids the need for an index file.
        let global_read_counter = AtomicU64::new(0);
        let mut result = ChromResult::new(interner.len());

        let mut bam = bam::Reader::from_path(bam_path)
            .with_context(|| format!("Failed to open alignment file: {}", bam_path))?;
        if let Some(ref_path) = reference {
            bam.set_reference(ref_path)
                .with_context(|| format!("Failed to set reference FASTA: {}", ref_path))?;
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

            // Count BAM stat categories before skip filters
            if flags & BAM_FUNMAP != 0 {
                result.stat_unmapped += 1;
            }
            if flags & BAM_FQCFAIL != 0 {
                result.stat_qc_failed += 1;
            }
            if flags & BAM_FSECONDARY != 0 {
                result.stat_secondary += 1;
            }
            if flags & BAM_FSUPPLEMENTARY != 0 {
                result.stat_supplementary += 1;
            }

            if flags & BAM_FUNMAP != 0 {
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

            if is_reverse {
                result.stat_minus_strand += 1;
            } else {
                result.stat_plus_strand += 1;
            }
            if flags & BAM_FPAIRED != 0 {
                if is_read1 {
                    result.stat_read1 += 1;
                } else {
                    result.stat_read2 += 1;
                }
                if flags & BAM_FPROPER_PAIR != 0 {
                    result.stat_proper_pairs += 1;
                }
                if flags & BAM_FMUNMAP != 0 {
                    result.stat_singletons += 1;
                } else if record.mtid() != record.tid() {
                    result.stat_mate_mapped_diff_chr += 1;
                }
            }

            {
                use rust_htslib::bam::record::Cigar;
                let has_splice = record
                    .cigar()
                    .iter()
                    .any(|op| matches!(op, Cigar::RefSkip(_)));
                if has_splice {
                    result.stat_splice_reads += 1;
                } else {
                    result.stat_non_splice_reads += 1;
                }
            }

            let tid = record.tid();
            if tid < 0 || tid as usize >= tid_to_name.len() {
                continue;
            }
            let bam_chrom = &tid_to_name[tid as usize];
            let prefixed_chrom;
            let chrom = if let Some(mapped) = chrom_mapping.get(bam_chrom.as_str()) {
                mapped.as_str()
            } else if let Some(prefix) = chrom_prefix {
                prefixed_chrom = format!("{}{}", prefix, bam_chrom);
                prefixed_chrom.as_str()
            } else {
                bam_chrom.as_str()
            };

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
            let read_name = record.qname().to_vec();
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
                    read_name.clone(),
                    own_tid,
                    own_pos,
                    mate_tid_val,
                    mate_pos_val,
                    hi_tag,
                )
            } else {
                (
                    read_name.clone(),
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

                let combined_genes: Vec<GeneIdx> =
                    if mate_info.gene_hits.is_empty() && gene_hits.is_empty() {
                        Vec::new()
                    } else {
                        let mut scores: HashMap<GeneIdx, u8> = HashMap::new();
                        for &g in &mate_info.gene_hits {
                            *scores.entry(g).or_insert(0) += 1;
                        }
                        for &g in &gene_hits {
                            *scores.entry(g).or_insert(0) += 1;
                        }
                        let max_score = scores.values().copied().max().unwrap_or(0);
                        let mut best: Vec<GeneIdx> = scores
                            .iter()
                            .filter(|(_, &s)| s == max_score)
                            .map(|(&g, _)| g)
                            .collect();
                        best.sort_unstable();
                        best
                    };

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
                mate_buffer.insert(
                    buffer_key,
                    MateInfo {
                        gene_hits: gene_hits.clone(),
                        is_dup,
                        is_multi,
                    },
                );
            }
        }

        result.unmatched_mates = mate_buffer;
        result
    };

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
                    if mate_info.gene_hits.is_empty() && info.gene_hits.is_empty() {
                        Vec::new()
                    } else {
                        let mut scores: HashMap<GeneIdx, u8> = HashMap::new();
                        for &g in &mate_info.gene_hits {
                            *scores.entry(g).or_insert(0) += 1;
                        }
                        for &g in &info.gene_hits {
                            *scores.entry(g).or_insert(0) += 1;
                        }
                        let max_score = scores.values().copied().max().unwrap_or(0);
                        let mut best: Vec<GeneIdx> = scores
                            .iter()
                            .filter(|(_, &s)| s == max_score)
                            .map(|(&g, _)| g)
                            .collect();
                        best.sort_unstable();
                        best
                    };

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
        total_reads_multi_nodup: merged.n_multi_nodup,
        total_reads_unique_dup: merged.n_unique_dup,
        total_reads_unique_nodup: merged.n_unique_nodup,
        stat_total_reads: merged.total_reads,
        stat_assigned: merged.stat_assigned,
        stat_ambiguous: merged.stat_ambiguous,
        stat_no_features: merged.stat_no_features,
        stat_total_fragments: merged.total_fragments,
        stat_total_dup: merged.total_dup,
        stat_total_multi: merged.total_multi,
        stat_unmapped: merged.stat_unmapped,
        stat_qc_failed: merged.stat_qc_failed,
        stat_secondary: merged.stat_secondary,
        stat_supplementary: merged.stat_supplementary,
        stat_proper_pairs: merged.stat_proper_pairs,
        stat_read1: merged.stat_read1,
        stat_read2: merged.stat_read2,
        stat_plus_strand: merged.stat_plus_strand,
        stat_minus_strand: merged.stat_minus_strand,
        stat_singletons: merged.stat_singletons,
        stat_mate_mapped_diff_chr: merged.stat_mate_mapped_diff_chr,
        stat_splice_reads: merged.stat_splice_reads,
        stat_non_splice_reads: merged.stat_non_splice_reads,
        stat_non_primary_hits: merged.stat_secondary + merged.stat_supplementary,
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
    fn test_dup_check_gatk_markduplicates() {
        let header = "@PG\tID:GATK MarkDuplicates\tPN:GATK MarkDuplicates\tVN:4.4.0";
        assert!(header_has_dup_marker(header));
    }

    #[test]
    fn test_dup_check_dup_marker_only_in_cl_no_match() {
        // A tool that mentions MarkDuplicates only in the CL field should not match
        // (e.g., a wrapper script that calls MarkDuplicates but has its own ID/PN)
        let header = "@PG\tID:my_pipeline\tPN:my_pipeline\tCL:java -jar picard.jar MarkDuplicates";
        assert!(!header_has_dup_marker(header));
    }
}
