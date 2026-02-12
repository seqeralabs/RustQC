//! Read counting engine for BAM files.
//!
//! Assigns reads from a BAM file to genes based on GTF annotation,
//! producing four count vectors matching the dupRadar approach:
//! 1. All reads including multimappers (with duplicates)
//! 2. All reads including multimappers (without duplicates)
//! 3. Uniquely mapped reads only (with duplicates)
//! 4. Uniquely mapped reads only (without duplicates)
//!
//! This implements a simplified featureCounts-compatible counting strategy.

use crate::gtf::Gene;
use anyhow::{Context, Result};
use indexmap::IndexMap;
use log::{debug, info};
use rust_htslib::bam::{self, Read as BamRead};
use std::collections::HashMap;

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
    gene_hits: &[String],
    gene_counts: &mut IndexMap<String, GeneCounts>,
    is_dup: bool,
    is_multi: bool,
) -> bool {
    if gene_hits.len() == 1 {
        if let Some(counts) = gene_counts.get_mut(gene_hits[0].as_str()) {
            counts.count_fragment(is_dup, is_multi);
            return true;
        }
    }
    false
}

/// Result of the counting step, including per-gene counts and total mapped reads.
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
}

/// An interval tree node for efficient overlap queries.
/// Uses a simple sorted-interval approach with binary search.
#[derive(Debug, Clone)]
struct GeneInterval {
    /// Start position (0-based, half-open for internal use)
    start: u64,
    /// End position (0-based, half-open)
    end: u64,
    /// Gene ID index in the genes map
    gene_id: String,
    /// Strand
    strand: char,
}

/// A simple interval index for a single chromosome.
/// Intervals are sorted by start position for binary search.
/// Augmented with max_end for efficient overlap queries.
#[derive(Debug)]
struct ChromIndex {
    intervals: Vec<GeneInterval>,
    /// For each position i, the maximum end value among intervals[0..=i].
    /// This allows pruning the scan: if max_end[i] <= query_start, all
    /// intervals from 0..=i can be skipped.
    max_end: Vec<u64>,
}

impl ChromIndex {
    fn new(mut intervals: Vec<GeneInterval>) -> Self {
        intervals.sort_by_key(|iv| (iv.start, iv.end));

        // Build max_end array
        let mut max_end = Vec::with_capacity(intervals.len());
        let mut running_max = 0u64;
        for iv in &intervals {
            running_max = running_max.max(iv.end);
            max_end.push(running_max);
        }

        ChromIndex { intervals, max_end }
    }

    /// Find all gene intervals overlapping the query range [start, end).
    /// An interval overlaps if interval.start < end AND interval.end > start.
    fn query(&self, start: u64, end: u64) -> Vec<&GeneInterval> {
        let mut results = Vec::new();
        if self.intervals.is_empty() {
            return results;
        }

        // Quick check: if the max end of all intervals <= start, no overlaps possible
        if *self.max_end.last().unwrap() <= start {
            return results;
        }

        // Find the first interval where iv.start >= end using binary search.
        // All intervals at or after this index have start >= end, so they can't overlap.
        let upper = self.intervals.partition_point(|iv| iv.start < end);

        // Now scan intervals [0..upper) and check if iv.end > start.
        // Use max_end to find a lower bound: find the leftmost position where
        // max_end > start. All intervals before that can't overlap.
        // Binary search on max_end (which is monotonically non-decreasing).
        let lower = self.max_end[..upper].partition_point(|&me| me <= start);

        for iv in &self.intervals[lower..upper] {
            if iv.end > start {
                results.push(iv);
            }
        }

        results
    }
}

/// Build a spatial index from gene annotations for fast overlap queries.
fn build_index(genes: &IndexMap<String, Gene>) -> HashMap<String, ChromIndex> {
    let mut chrom_intervals: HashMap<String, Vec<GeneInterval>> = HashMap::new();

    for gene in genes.values() {
        // Add each exon as an interval (featureCounts assigns reads at exon level)
        for exon in &gene.exons {
            let interval = GeneInterval {
                // Convert from 1-based inclusive GTF to 0-based half-open
                start: exon.start - 1,
                end: exon.end,
                gene_id: gene.gene_id.clone(),
                strand: exon.strand,
            };
            chrom_intervals
                .entry(exon.chrom.clone())
                .or_default()
                .push(interval);
        }
    }

    chrom_intervals
        .into_iter()
        .map(|(chrom, intervals)| (chrom, ChromIndex::new(intervals)))
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

/// Extract aligned (match) blocks from a CIGAR string.
///
/// Walks the CIGAR operations and returns genomic intervals for each
/// M/=/X operation. N (intron skip) and D (deletion) operations advance
/// the reference position without producing an aligned block.
/// This ensures spliced reads don't falsely overlap with genes in introns.
fn cigar_to_aligned_blocks(
    start: u64,
    cigar: &rust_htslib::bam::record::CigarStringView,
) -> Vec<(u64, u64)> {
    use rust_htslib::bam::record::Cigar;
    let mut blocks = Vec::new();
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

    blocks
}

/// Metadata collected per read for paired-end mate buffering.
/// When counting paired-end fragments, we need to see both mates to
/// properly determine gene overlap (featureCounts checks both mates
/// independently). This struct stores the information we need from each mate.
struct MateInfo {
    /// Gene IDs this mate overlaps (after strand filtering, deduplicated)
    gene_hits: Vec<String>,
    /// Whether this read is a PCR/optical duplicate
    is_dup: bool,
    /// Whether this read is a multimapper (NH > 1)
    is_multi: bool,
}

/// Key for the mate buffer: (read_name, expected_pair_key)
/// For a read at position P with mate at position M, the pair key is (min(P,M), max(P,M))
/// so both mates of the same alignment pair produce the same key.
type MateBufferKey = (Vec<u8>, (i64, i64));

/// Count reads from a BAM file and assign them to genes.
///
/// Performs four simultaneous counting modes matching dupRadar's approach:
/// - With/without multimappers
/// - With/without PCR duplicates
///
/// For paired-end data, this function buffers mates by read name and combines
/// their gene overlaps to match featureCounts' fragment-level counting behavior.
/// A fragment is assigned to a gene if either mate overlaps that gene's exons.
///
/// # Arguments
/// * `bam_path` - Path to the duplicate-marked BAM file
/// * `genes` - Gene annotation map from GTF parsing
/// * `stranded` - Library strandedness (0, 1, or 2)
/// * `paired` - Whether the library is paired-end
/// * `threads` - Number of threads for BAM reading
pub fn count_reads(
    bam_path: &str,
    genes: &IndexMap<String, Gene>,
    stranded: u8,
    paired: bool,
    threads: usize,
    chrom_mapping: &HashMap<String, String>,
    chrom_prefix: Option<&str>,
) -> Result<CountResult> {
    // Build spatial index
    let index = build_index(genes);

    // Initialize gene counts
    let mut gene_counts: IndexMap<String, GeneCounts> = IndexMap::new();
    for gene_id in genes.keys() {
        gene_counts.insert(gene_id.clone(), GeneCounts::default());
    }

    // Open BAM file
    let mut bam = bam::Reader::from_path(bam_path)
        .with_context(|| format!("Failed to open BAM file: {}", bam_path))?;

    if threads > 1 {
        bam.set_threads(threads - 1)
            .context("Failed to set BAM reader threads")?;
    }

    // Get chromosome names from header
    let header = bam.header().clone();
    let tid_to_name: Vec<String> = (0..header.target_count())
        .map(|tid| String::from_utf8_lossy(header.tid2name(tid)).to_string())
        .collect();

    // Track statistics
    let mut total_reads: u64 = 0;
    let mut total_mapped: u64 = 0;
    let mut total_dup: u64 = 0;
    let mut total_multi: u64 = 0;
    let mut total_fragments: u64 = 0;

    // featureCounts-style assignment statistics
    let mut stat_assigned: u64 = 0;
    let mut stat_ambiguous: u64 = 0;
    let mut stat_no_features: u64 = 0;
    // Totals for N calculation (mapped reads per mode)
    let mut n_multi_dup: u64 = 0;
    let mut n_multi_nodup: u64 = 0;
    let mut n_unique_dup: u64 = 0;
    let mut n_unique_nodup: u64 = 0;

    // For paired-end: buffer mates by (read_name, position_pair) to combine
    // their overlaps. Using position pairs ensures secondary alignments at
    // different locations form separate fragments (matching featureCounts).
    let mut mate_buffer: HashMap<MateBufferKey, MateInfo> = HashMap::new();

    // Iterate over all records
    let mut record = bam::Record::new();
    while let Some(result) = bam.read(&mut record) {
        result.context("Error reading BAM record")?;
        total_reads += 1;

        if total_reads.is_multiple_of(5_000_000) {
            debug!("Processed {} reads...", total_reads);
        }

        let flags = record.flags();

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

        total_mapped += 1;

        let is_dup = flags & BAM_FDUP != 0;
        if is_dup {
            total_dup += 1;
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
            total_multi += 1;
        }

        let is_reverse = flags & BAM_FREVERSE != 0;
        let is_read1 = flags & BAM_FREAD1 != 0;

        // Get the chromosome name
        let tid = record.tid();
        if tid < 0 || tid as usize >= tid_to_name.len() {
            continue;
        }
        let bam_chrom = &tid_to_name[tid as usize];
        // Apply chromosome name prefix and/or mapping (BAM name -> GTF name)
        let prefixed_chrom;
        let chrom = if let Some(mapped) = chrom_mapping.get(bam_chrom.as_str()) {
            // Explicit mapping takes priority
            mapped.as_str()
        } else if let Some(prefix) = chrom_prefix {
            // Apply prefix: e.g. "chr" + "1" -> "chr1"
            prefixed_chrom = format!("{}{}", prefix, bam_chrom);
            prefixed_chrom.as_str()
        } else {
            bam_chrom.as_str()
        };

        // Find gene overlaps using CIGAR-aware aligned blocks
        let gene_hits = if let Some(chrom_idx) = index.get(chrom) {
            // Extract aligned blocks from CIGAR (M/=/X operations only).
            // This avoids false overlaps with genes in introns of spliced reads.
            let aligned_blocks = cigar_to_aligned_blocks(record.pos() as u64, &record.cigar());

            let mut overlaps = Vec::new();
            for (block_start, block_end) in &aligned_blocks {
                overlaps.extend(chrom_idx.query(*block_start, *block_end));
            }

            // Filter by strand and deduplicate gene IDs
            let mut genes_hit: Vec<String> = overlaps
                .iter()
                .filter(|iv| strand_matches(is_reverse, is_read1, paired, iv.strand, stranded))
                .map(|iv| iv.gene_id.clone())
                .collect();
            genes_hit.sort_unstable();
            genes_hit.dedup();
            genes_hit
        } else {
            // Chromosome not in annotation
            Vec::new()
        };

        // --- Single-end counting: assign directly ---
        if !paired {
            // Update N totals
            n_multi_dup += 1;
            n_unique_dup += 1;
            if !is_dup {
                n_multi_nodup += 1;
                n_unique_nodup += 1;
            }
            total_fragments += 1;

            // Assign to gene if unambiguous
            if gene_hits.is_empty() {
                stat_no_features += 1;
            } else if gene_hits.len() > 1 {
                stat_ambiguous += 1;
            } else if assign_fragment_to_gene(&gene_hits, &mut gene_counts, is_dup, is_multi) {
                stat_assigned += 1;
            }
            continue;
        }

        // --- Paired-end counting: buffer mates and combine ---
        //
        // featureCounts with isPairedEnd=TRUE and countReadPairs=TRUE assigns
        // a fragment to a gene if EITHER mate overlaps that gene's exons.
        // We buffer the first mate we see and combine overlaps when the second
        // mate arrives. The duplicate flag from read1 is used for the fragment.
        //
        // For multi-mapped reads with secondary alignments, each alignment
        // location pair is counted as a separate fragment. We use
        // (read_name, position_pair) as the buffer key to distinguish them.
        let read_name = record.qname().to_vec();
        let own_pos = record.pos();
        let mate_pos_val = record.mpos();

        // Create a position pair key that both mates will compute identically
        let pos_pair = if own_pos <= mate_pos_val {
            (own_pos, mate_pos_val)
        } else {
            (mate_pos_val, own_pos)
        };
        let buffer_key: MateBufferKey = (read_name.clone(), pos_pair);

        if let Some(mate_info) = mate_buffer.remove(&buffer_key) {
            // We've seen the other mate - combine overlaps and assign
            total_fragments += 1;

            // For the fragment, use read1's dup/multi status (featureCounts
            // considers a fragment as duplicate if read1 is flagged as duplicate)
            let frag_is_dup = if is_read1 { is_dup } else { mate_info.is_dup };
            let frag_is_multi = if is_read1 {
                is_multi
            } else {
                mate_info.is_multi
            };

            // Update N totals (once per fragment)
            n_multi_dup += 1;
            n_unique_dup += 1;
            if !frag_is_dup {
                n_multi_nodup += 1;
                n_unique_nodup += 1;
            }

            // Combine gene hits from both mates using INTERSECTION strategy
            // (matching featureCounts behaviour).
            //
            // If both mates overlap genes, we use the INTERSECTION of their
            // gene sets.  This resolves spurious ambiguity when one mate
            // barely clips a neighbouring gene's exon while the other mate
            // clearly overlaps only the intended gene.
            //
            // If only one mate has gene hits (the other overlaps nothing),
            // we use that mate's gene set directly.
            let combined_genes: Vec<String> =
                if mate_info.gene_hits.is_empty() && gene_hits.is_empty() {
                    Vec::new()
                } else if mate_info.gene_hits.is_empty() {
                    gene_hits
                } else if gene_hits.is_empty() {
                    mate_info.gene_hits
                } else {
                    // Both mates have gene hits - use INTERSECTION
                    let set_a: std::collections::HashSet<&String> =
                        mate_info.gene_hits.iter().collect();
                    let intersection: Vec<String> = gene_hits
                        .iter()
                        .filter(|g| set_a.contains(g))
                        .cloned()
                        .collect();
                    if intersection.is_empty() {
                        // Mates disagree on gene assignment - treat as ambiguous
                        // by returning the union (which will have len > 1)
                        let mut union = mate_info.gene_hits;
                        union.extend(gene_hits);
                        union.sort_unstable();
                        union.dedup();
                        union
                    } else {
                        intersection
                    }
                };

            // Assign to gene if unambiguous (exactly one gene from combined overlaps)
            if combined_genes.is_empty() {
                stat_no_features += 1;
            } else if combined_genes.len() > 1 {
                stat_ambiguous += 1;
            } else if assign_fragment_to_gene(
                &combined_genes,
                &mut gene_counts,
                frag_is_dup,
                frag_is_multi,
            ) {
                stat_assigned += 1;
            }
        } else {
            // First mate seen - buffer it and wait for the other mate
            mate_buffer.insert(
                buffer_key,
                MateInfo {
                    gene_hits,
                    is_dup,
                    is_multi,
                },
            );
        }
    }

    // Handle any unpaired mates left in the buffer (singletons).
    // These are reads whose mate was filtered out (e.g., mate unmapped,
    // mate on a different chromosome and filtered by featureCounts, etc.)
    // featureCounts by default still counts these as fragments.
    for (_key, mate_info) in mate_buffer.drain() {
        total_fragments += 1;

        n_multi_dup += 1;
        n_unique_dup += 1;
        if !mate_info.is_dup {
            n_multi_nodup += 1;
            n_unique_nodup += 1;
        }

        if mate_info.gene_hits.is_empty() {
            stat_no_features += 1;
        } else if mate_info.gene_hits.len() > 1 {
            stat_ambiguous += 1;
        } else if assign_fragment_to_gene(
            &mate_info.gene_hits,
            &mut gene_counts,
            mate_info.is_dup,
            mate_info.is_multi,
        ) {
            stat_assigned += 1;
        }
    }

    info!(
        "Read {} total reads, {} mapped, {} fragments, {} duplicates, {} multimappers",
        total_reads, total_mapped, total_fragments, total_dup, total_multi
    );
    info!(
        "Assignment stats: {} assigned, {} ambiguous, {} no_features (total fragments: {})",
        stat_assigned, stat_ambiguous, stat_no_features, total_fragments
    );

    // Detect chromosome name mismatch: if we processed mapped reads but
    // no genes received any counts, the GTF and BAM likely use different
    // chromosome naming conventions (e.g., "chr1" vs "1").
    let genes_with_reads = gene_counts.values().filter(|c| c.all_multi > 0).count();
    if total_mapped > 0 && genes_with_reads == 0 {
        // Collect example chromosome names from each source for the error message
        let bam_chroms: Vec<&str> = tid_to_name.iter().take(5).map(|s| s.as_str()).collect();
        let gtf_chroms: Vec<&str> = index.keys().take(5).map(|s| s.as_str()).collect();
        anyhow::bail!(
            "Chromosome name mismatch: no BAM reads could be assigned to any gene.\n\
             \n\
             BAM chromosomes (first 5): {}\n\
             GTF chromosomes (first 5): {}\n\
             \n\
             The BAM and GTF files appear to use different chromosome naming conventions.\n\
             To fix this, create a YAML config file with a chromosome_mapping section and pass it via --config:\n\
             \n\
             Example config.yaml:\n\
             \n\
             chromosome_mapping:\n\
             {}",
            bam_chroms.join(", "),
            gtf_chroms.join(", "),
            bam_chroms.iter().zip(gtf_chroms.iter())
                .map(|(b, g)| format!("  {}: {}", g, b))
                .collect::<Vec<_>>()
                .join("\n")
        );
    }

    Ok(CountResult {
        gene_counts,
        total_reads_multi_dup: n_multi_dup,
        total_reads_multi_nodup: n_multi_nodup,
        total_reads_unique_dup: n_unique_dup,
        total_reads_unique_nodup: n_unique_nodup,
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
}
