//! Qualimap-compatible annotation index for RNA-Seq QC.
//!
//! Builds per-chromosome COITrees for exon overlap queries and intron interval
//! trees. Each exon node stores its parent gene and transcript, enabling
//! Qualimap's enclosure-based gene assignment (every M-block must be fully
//! enclosed by exons of the same gene).

use coitrees::{COITree, Interval, IntervalTree};
use indexmap::IndexMap;
use log::debug;
use std::collections::HashMap;

use crate::gtf::Gene;

// ============================================================
// Exon node metadata stored in the COITree
// ============================================================

/// Metadata for a single exon interval in the COITree.
///
/// Each node represents one exon of one transcript of one gene.
/// Qualimap assigns reads by checking that all M-blocks are *enclosed*
/// (fully contained) within exons of the same gene.
#[derive(Debug, Clone, Copy)]
pub struct ExonMeta {
    /// Index into the gene list (position in the IndexMap).
    #[allow(dead_code)]
    pub gene_idx: u32,
    /// Index into the transcript list within the gene.
    #[allow(dead_code)]
    pub transcript_idx: u16,
    /// Exon start in 0-based half-open coordinates.
    #[allow(dead_code)]
    pub exon_start: i32,
    /// Exon end in 0-based half-open coordinates.
    #[allow(dead_code)]
    pub exon_end: i32,
    /// Strand ('+', '-', or '.').
    #[allow(dead_code)]
    pub strand: u8,
}

impl Default for ExonMeta {
    fn default() -> Self {
        Self {
            gene_idx: 0,
            transcript_idx: 0,
            exon_start: 0,
            exon_end: 0,
            strand: b'.',
        }
    }
}

/// Per-chromosome exon COITree for overlap queries.
pub type ExonTree = COITree<ExonMeta, u32>;

// ============================================================
// Intron node metadata
// ============================================================

/// Metadata for an intron interval (gap between consecutive exons of a transcript).
#[derive(Debug, Clone, Copy, Default)]
#[allow(dead_code)]
pub struct IntronMeta {
    /// Index into the gene list.
    pub gene_idx: u32,
}

/// Per-chromosome intron COITree for intronic read classification.
pub type IntronTree = COITree<IntronMeta, u32>;

// ============================================================
// Merged exon node metadata
// ============================================================

/// Metadata for a merged exon interval in the COITree.
///
/// Qualimap Java merges overlapping/abutting exons per gene before building
/// its IntervalTree. This merged representation is used for the enclosure
/// check (gene assignment). A read block enclosed by a merged interval is
/// considered "exonic" even if no single transcript exon fully contains it.
#[derive(Debug, Clone, Copy)]
pub struct MergedExonMeta {
    /// Index into the gene list (position in the IndexMap).
    pub gene_idx: u32,
    /// Exon start in 0-based half-open coordinates.
    pub exon_start: i32,
    /// Exon end in 0-based half-open coordinates.
    pub exon_end: i32,
    /// Strand ('+', '-', or '.').
    #[allow(dead_code)]
    pub strand: u8,
}

impl Default for MergedExonMeta {
    fn default() -> Self {
        Self {
            gene_idx: 0,
            exon_start: 0,
            exon_end: 0,
            strand: b'.',
        }
    }
}

/// Per-chromosome merged exon COITree for enclosure-based gene assignment.
///
/// Uses merged per-gene exon intervals (union of all transcript exons) to match
/// Qualimap Java's interval tree behavior.
pub type MergedExonTree = COITree<MergedExonMeta, u32>;

// ============================================================
// Per-transcript metadata
// ============================================================

/// Pre-computed metadata for a single transcript, used for coverage tracking.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct TranscriptInfo {
    /// Gene index (position in the IndexMap).
    pub gene_idx: u32,
    /// Transcript index within the gene.
    pub transcript_idx: u16,
    /// Gene ID string (for output).
    pub gene_id: String,
    /// Transcript ID string.
    pub transcript_id: String,
    /// Chromosome name.
    pub chrom: String,
    /// Transcript start (0-based).
    pub start: i32,
    /// Transcript end (0-based half-open).
    pub end: i32,
    /// Strand.
    pub strand: char,
    /// Total transcript length (end - start).
    pub length: u32,
    /// Sorted exon intervals in 0-based half-open coordinates.
    pub exons: Vec<(i32, i32)>,
    /// Total exonic bases (sum of exon lengths).
    pub exonic_length: u32,
}

// ============================================================
// Merged gene model — union of all transcript exons per gene
// ============================================================

/// Merged exon model for a gene: union of all exons across all transcripts.
///
/// Qualimap (via Picard's `Gene.Transcript`) merges all exons from all isoforms
/// into a single non-redundant exon set per gene. Coverage is tracked on this
/// merged model, which gives the correct bias and coverage profile values.
#[derive(Debug, Clone)]
pub struct MergedGeneModel {
    /// Gene index (position in the IndexMap).
    #[allow(dead_code)]
    pub gene_idx: u32,
    /// Strand ('+', '-', or '.').
    #[allow(dead_code)]
    pub strand: char,
    /// Merged non-overlapping exon intervals in 0-based half-open coordinates, sorted by start.
    pub exons: Vec<(i32, i32)>,
    /// Total merged exonic length (sum of merged exon lengths).
    pub exonic_length: u32,
}

impl MergedGeneModel {
    /// Build a merged gene model from all transcripts of a gene.
    ///
    /// Takes all exon intervals (0-based half-open) from all transcripts,
    /// sorts them, and merges overlapping/adjacent intervals into a single
    /// non-redundant set.
    pub fn from_transcripts(gene_idx: u32, strand: char, all_exons: &[(i32, i32)]) -> Self {
        let exons = merge_intervals(all_exons);
        let exonic_length: u32 = exons.iter().map(|(s, e)| (e - s) as u32).sum();
        Self {
            gene_idx,
            strand,
            exons,
            exonic_length,
        }
    }
}

/// Merge overlapping and adjacent intervals into non-redundant sorted intervals.
///
/// Input: arbitrary list of (start, end) 0-based half-open intervals.
/// Output: sorted, non-overlapping intervals where overlapping/touching intervals
/// have been merged.
fn merge_intervals(intervals: &[(i32, i32)]) -> Vec<(i32, i32)> {
    if intervals.is_empty() {
        return Vec::new();
    }
    let mut sorted: Vec<(i32, i32)> = intervals.to_vec();
    sorted.sort_unstable_by_key(|(s, _)| *s);

    let mut merged = Vec::with_capacity(sorted.len());
    let (mut cur_start, mut cur_end) = sorted[0];

    for &(start, end) in &sorted[1..] {
        if start < cur_end {
            // Strictly overlapping — extend (abutting exons are kept separate)
            cur_end = cur_end.max(end);
        } else {
            // Gap — emit current and start new
            merged.push((cur_start, cur_end));
            cur_start = start;
            cur_end = end;
        }
    }
    merged.push((cur_start, cur_end));
    merged
}

// ============================================================
// QualimapIndex — the main index structure
// ============================================================

/// Complete Qualimap annotation index built from GTF genes.
///
/// Contains per-chromosome exon and intron trees, plus transcript metadata
/// and merged gene models for coverage tracking.
/// This index is *separate* from the existing dupRadar/featureCounts index
/// because Qualimap uses fundamentally different assignment logic.
pub struct QualimapIndex {
    /// Per-chromosome exon COITrees (key = chromosome name).
    /// Individual per-transcript exons — used for coverage tracking.
    pub exon_trees: HashMap<String, ExonTree>,
    /// Per-chromosome merged exon COITrees (key = chromosome name).
    /// Merged per-gene exons — used for enclosure-based gene assignment.
    /// Matches Qualimap Java's `GenomicRegionSet.addRegion()` which merges
    /// overlapping/abutting exons per gene via `concatenateIntervals()`.
    pub merged_exon_trees: HashMap<String, MergedExonTree>,
    /// Per-chromosome intron COITrees (key = chromosome name).
    pub intron_trees: HashMap<String, IntronTree>,
    /// All transcript metadata, indexed by a flat transcript ID.
    /// The index is `(gene_idx, transcript_idx)` mapped to a flat position.
    pub transcripts: Vec<TranscriptInfo>,
    /// Lookup: gene_idx -> range of transcript indices in `transcripts` vec.
    /// `gene_transcript_ranges[gene_idx] = (start, end)` half-open into `transcripts`.
    pub gene_transcript_ranges: Vec<(u32, u32)>,
    /// Merged gene models: one per gene, indexed by gene_idx.
    /// Each merges all exons from all transcripts into non-redundant intervals.
    pub merged_gene_models: Vec<MergedGeneModel>,
    /// Total number of genes.
    #[allow(dead_code)]
    pub num_genes: u32,
}

impl QualimapIndex {
    /// Build the Qualimap index from GTF gene annotations.
    ///
    /// Creates per-chromosome exon COITrees (with enclosure metadata) and intron
    /// interval trees. Also builds per-transcript metadata for coverage tracking.
    ///
    /// # Arguments
    /// * `genes` - Gene annotations from GTF parsing (insertion-order preserved).
    pub fn from_genes(genes: &IndexMap<String, Gene>) -> Self {
        // Collect exon intervals and intron intervals per chromosome
        let mut exon_intervals: HashMap<String, Vec<Interval<ExonMeta>>> = HashMap::new();
        let mut intron_intervals: HashMap<String, Vec<Interval<IntronMeta>>> = HashMap::new();
        let mut transcripts = Vec::new();
        let mut gene_transcript_ranges = Vec::with_capacity(genes.len());

        let mut merged_gene_models = Vec::with_capacity(genes.len());

        for (gene_idx, gene) in genes.values().enumerate() {
            let gene_idx = gene_idx as u32;
            let tx_start_idx = transcripts.len() as u32;

            // Collect all exons from all transcripts for merged gene model
            let mut all_gene_exons: Vec<(i32, i32)> = Vec::new();

            for (tx_idx, tx) in gene.transcripts.iter().enumerate() {
                let tx_idx = tx_idx as u16;

                // Convert transcript exons from 1-based inclusive (GTF) to 0-based half-open
                let mut exons_0based: Vec<(i32, i32)> = tx
                    .exons
                    .iter()
                    .map(|(s, e)| ((*s as i32) - 1, *e as i32))
                    .collect();
                exons_0based.sort_unstable_by_key(|(s, _)| *s);

                let exonic_length: u32 = exons_0based.iter().map(|(s, e)| (e - s) as u32).sum();

                let tx_start = exons_0based.first().map(|(s, _)| *s).unwrap_or(0);
                let tx_end = exons_0based.last().map(|(_, e)| *e).unwrap_or(0);

                // Collect for merged gene model
                all_gene_exons.extend_from_slice(&exons_0based);

                // Add exon intervals to the per-chrom collection
                let chrom_exons = exon_intervals.entry(tx.chrom.clone()).or_default();
                for &(start, end) in &exons_0based {
                    chrom_exons.push(Interval::new(
                        start,
                        end,
                        ExonMeta {
                            gene_idx,
                            transcript_idx: tx_idx,
                            exon_start: start,
                            exon_end: end,
                            strand: tx.strand as u8,
                        },
                    ));
                }

                // Add intron intervals (gaps between consecutive exons)
                let chrom_introns = intron_intervals.entry(tx.chrom.clone()).or_default();
                for window in exons_0based.windows(2) {
                    let intron_start = window[0].1;
                    let intron_end = window[1].0;
                    if intron_end > intron_start {
                        chrom_introns.push(Interval::new(
                            intron_start,
                            intron_end,
                            IntronMeta { gene_idx },
                        ));
                    }
                }

                // Store transcript metadata
                transcripts.push(TranscriptInfo {
                    gene_idx,
                    transcript_idx: tx_idx,
                    gene_id: gene.gene_id.clone(),
                    transcript_id: tx.transcript_id.clone(),
                    chrom: tx.chrom.clone(),
                    start: tx_start,
                    end: tx_end,
                    strand: tx.strand,
                    length: (tx_end - tx_start) as u32,
                    exons: exons_0based,
                    exonic_length,
                });
            }

            let tx_end_idx = transcripts.len() as u32;
            gene_transcript_ranges.push((tx_start_idx, tx_end_idx));

            // Build merged gene model from all exons across all transcripts
            merged_gene_models.push(MergedGeneModel::from_transcripts(
                gene_idx,
                gene.strand,
                &all_gene_exons,
            ));
        }

        // Build merged exon intervals from merged gene models.
        // Qualimap Java merges overlapping/abutting exons per gene before inserting
        // them into its IntervalTree. We build a separate COITree from these merged
        // intervals for the enclosure-based gene assignment check.
        let mut merged_exon_intervals: HashMap<String, Vec<Interval<MergedExonMeta>>> =
            HashMap::new();
        for (gene_idx_usize, gene) in genes.values().enumerate() {
            let gene_idx = gene_idx_usize as u32;
            let model = &merged_gene_models[gene_idx_usize];
            // Use the first transcript's chromosome (all transcripts of a gene share the same chrom)
            if let Some(tx) = gene.transcripts.first() {
                let chrom_merged = merged_exon_intervals.entry(tx.chrom.clone()).or_default();
                for &(start, end) in &model.exons {
                    chrom_merged.push(Interval::new(
                        start,
                        end,
                        MergedExonMeta {
                            gene_idx,
                            exon_start: start,
                            exon_end: end,
                            strand: model.strand as u8,
                        },
                    ));
                }
            }
        }

        // Build COITrees from collected intervals
        let total_exon_ivs: usize = exon_intervals.values().map(|v| v.len()).sum();
        let total_merged_ivs: usize = merged_exon_intervals.values().map(|v| v.len()).sum();
        eprintln!(
            "QM_INDEX: {} exon intervals ({} merged), {} transcripts, {} genes",
            total_exon_ivs,
            total_merged_ivs,
            transcripts.len(),
            gene_transcript_ranges.len()
        );
        let exon_trees: HashMap<String, ExonTree> = exon_intervals
            .into_iter()
            .map(|(chrom, intervals)| (chrom, COITree::new(&intervals)))
            .collect();

        let merged_exon_trees: HashMap<String, MergedExonTree> = merged_exon_intervals
            .into_iter()
            .map(|(chrom, intervals)| (chrom, COITree::new(&intervals)))
            .collect();

        let intron_trees: HashMap<String, IntronTree> = intron_intervals
            .into_iter()
            .map(|(chrom, intervals)| (chrom, COITree::new(&intervals)))
            .collect();

        let num_genes = genes.len() as u32;

        debug!(
            "Built Qualimap index: {} genes, {} transcripts, {} chromosomes with exons",
            num_genes,
            transcripts.len(),
            exon_trees.len(),
        );

        Self {
            exon_trees,
            merged_exon_trees,
            intron_trees,
            transcripts,
            gene_transcript_ranges,
            merged_gene_models,
            num_genes,
        }
    }

    /// Get the per-transcript exon tree for a chromosome (used for coverage tracking).
    #[allow(dead_code)]
    pub fn exon_tree(&self, chrom: &str) -> Option<&ExonTree> {
        self.exon_trees.get(chrom)
    }

    /// Get the merged exon tree for a chromosome (used for enclosure-based gene assignment).
    ///
    /// Qualimap Java merges overlapping/abutting exons per gene before checking
    /// enclosure. This tree contains those merged intervals, ensuring reads that
    /// span two overlapping transcript exons of the same gene are correctly
    /// classified as "exonic".
    pub fn merged_exon_tree(&self, chrom: &str) -> Option<&MergedExonTree> {
        self.merged_exon_trees.get(chrom)
    }

    /// Get the intron tree for a chromosome.
    pub fn intron_tree(&self, chrom: &str) -> Option<&IntronTree> {
        self.intron_trees.get(chrom)
    }

    /// Get transcript info entries for a gene by gene index.
    ///
    /// Returns a slice of `TranscriptInfo` for all transcripts of the gene.
    pub fn gene_transcripts(&self, gene_idx: u32) -> &[TranscriptInfo] {
        let (start, end) = self.gene_transcript_ranges[gene_idx as usize];
        &self.transcripts[start as usize..end as usize]
    }
}

// ============================================================
// Unit tests
// ============================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gtf::{Exon, Gene, Transcript};

    /// Helper to build a minimal gene for testing.
    fn make_gene(gene_id: &str, chrom: &str, strand: char, transcripts: Vec<Transcript>) -> Gene {
        let start = transcripts
            .iter()
            .flat_map(|t| t.exons.iter().map(|(s, _)| *s))
            .min()
            .unwrap_or(1);
        let end = transcripts
            .iter()
            .flat_map(|t| t.exons.iter().map(|(_, e)| *e))
            .max()
            .unwrap_or(1);
        let exons = transcripts
            .iter()
            .flat_map(|t| {
                t.exons.iter().map(|(s, e)| Exon {
                    chrom: chrom.to_string(),
                    start: *s,
                    end: *e,
                    strand,
                })
            })
            .collect();
        Gene {
            gene_id: gene_id.to_string(),
            chrom: chrom.to_string(),
            start,
            end,
            strand,
            exons,
            effective_length: 0,
            attributes: HashMap::new(),
            transcripts,
        }
    }

    fn make_transcript(
        tx_id: &str,
        chrom: &str,
        strand: char,
        exons: Vec<(u64, u64)>,
    ) -> Transcript {
        let start = exons.iter().map(|(s, _)| *s).min().unwrap_or(1);
        let end = exons.iter().map(|(_, e)| *e).max().unwrap_or(1);
        Transcript {
            transcript_id: tx_id.to_string(),
            chrom: chrom.to_string(),
            start,
            end,
            strand,
            exons,
            cds_start: None,
            cds_end: None,
        }
    }

    #[test]
    fn test_index_from_single_gene() {
        let mut genes = IndexMap::new();
        let tx = make_transcript("tx1", "chr1", '+', vec![(101, 200), (301, 400)]);
        let gene = make_gene("gene1", "chr1", '+', vec![tx]);
        genes.insert("gene1".to_string(), gene);

        let index = QualimapIndex::from_genes(&genes);

        assert_eq!(index.num_genes, 1);
        assert_eq!(index.transcripts.len(), 1);
        assert!(index.exon_tree("chr1").is_some());
        assert!(index.intron_tree("chr1").is_some());
        assert!(index.exon_tree("chr2").is_none());

        // Check transcript metadata
        let ti = &index.transcripts[0];
        assert_eq!(ti.gene_id, "gene1");
        assert_eq!(ti.transcript_id, "tx1");
        assert_eq!(ti.exons.len(), 2);
        // GTF 1-based inclusive (101, 200) -> 0-based half-open (100, 200)
        assert_eq!(ti.exons[0], (100, 200));
        // GTF 1-based inclusive (301, 400) -> 0-based half-open (300, 400)
        assert_eq!(ti.exons[1], (300, 400));
        assert_eq!(ti.exonic_length, 200); // 100 + 100

        // Check gene transcript ranges
        let txs = index.gene_transcripts(0);
        assert_eq!(txs.len(), 1);
        assert_eq!(txs[0].transcript_id, "tx1");

        // Check that exon tree has entries
        let tree = index.exon_tree("chr1").unwrap();
        let mut hits = Vec::new();
        tree.query(150, 160, |iv| hits.push(iv.metadata.clone()));
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].gene_idx, 0);

        // Check intron tree: intron from 200..300
        let itree = index.intron_tree("chr1").unwrap();
        let mut ihits = Vec::new();
        itree.query(250, 260, |iv| ihits.push(iv.metadata.clone()));
        assert_eq!(ihits.len(), 1);
        assert_eq!(ihits[0].gene_idx, 0);
    }

    #[test]
    fn test_index_multiple_transcripts() {
        let mut genes = IndexMap::new();
        let tx1 = make_transcript("tx1", "chr1", '+', vec![(101, 200), (301, 400)]);
        let tx2 = make_transcript("tx2", "chr1", '+', vec![(101, 250), (351, 400)]);
        let gene = make_gene("gene1", "chr1", '+', vec![tx1, tx2]);
        genes.insert("gene1".to_string(), gene);

        let index = QualimapIndex::from_genes(&genes);

        assert_eq!(index.transcripts.len(), 2);
        let txs = index.gene_transcripts(0);
        assert_eq!(txs.len(), 2);
        assert_eq!(txs[0].transcript_id, "tx1");
        assert_eq!(txs[1].transcript_id, "tx2");
    }

    #[test]
    fn test_enclosure_check() {
        // Simulate what the accumulator will do: check if an M-block
        // is fully enclosed within any single exon of a gene.
        let mut genes = IndexMap::new();
        let tx = make_transcript("tx1", "chr1", '+', vec![(101, 200), (301, 400)]);
        let gene = make_gene("gene1", "chr1", '+', vec![tx]);
        genes.insert("gene1".to_string(), gene);

        let index = QualimapIndex::from_genes(&genes);
        let tree = index.exon_tree("chr1").unwrap();

        // M-block 120..180 (0-based) should be enclosed by exon (100, 200)
        let mut enclosed_genes = Vec::new();
        tree.query(120, 180, |iv| {
            if 120 >= iv.metadata.exon_start && 180 <= iv.metadata.exon_end {
                enclosed_genes.push(iv.metadata.gene_idx);
            }
        });
        assert_eq!(enclosed_genes, vec![0]);

        // M-block 190..310 (0-based) overlaps both exons but is NOT enclosed by either
        let mut enclosed_genes2 = Vec::new();
        tree.query(190, 310, |iv| {
            if 190 >= iv.metadata.exon_start && 310 <= iv.metadata.exon_end {
                enclosed_genes2.push(iv.metadata.gene_idx);
            }
        });
        assert!(enclosed_genes2.is_empty());
    }

    #[test]
    fn test_merged_exon_enclosure() {
        // Qualimap Java merges overlapping/abutting exons per gene before
        // building the interval tree. A read block that spans two overlapping
        // transcript exons is enclosed by the merged interval but NOT by
        // either individual exon.
        let mut genes = IndexMap::new();
        // tx1: exon [100, 200), tx2: exon [150, 300) → merged: [100, 300)
        let tx1 = make_transcript("tx1", "chr1", '+', vec![(101, 200)]);
        let tx2 = make_transcript("tx2", "chr1", '+', vec![(151, 300)]);
        let gene = make_gene("gene1", "chr1", '+', vec![tx1, tx2]);
        genes.insert("gene1".to_string(), gene);

        let index = QualimapIndex::from_genes(&genes);

        // Block [120, 250) is NOT enclosed by either individual exon
        let per_tx_tree = index.exon_tree("chr1").unwrap();
        let mut per_tx_enclosed = Vec::new();
        per_tx_tree.query(120, 250, |iv| {
            if 120 >= iv.metadata.exon_start && 250 <= iv.metadata.exon_end {
                per_tx_enclosed.push(iv.metadata.gene_idx);
            }
        });
        assert!(
            per_tx_enclosed.is_empty(),
            "Block should NOT be enclosed by any individual transcript exon"
        );

        // But it IS enclosed by the merged exon [100, 300)
        let merged_tree = index.merged_exon_tree("chr1").unwrap();
        let mut merged_enclosed = Vec::new();
        merged_tree.query(120, 250, |iv| {
            if 120 >= iv.metadata.exon_start && 250 <= iv.metadata.exon_end {
                merged_enclosed.push(iv.metadata.gene_idx);
            }
        });
        assert_eq!(
            merged_enclosed,
            vec![0],
            "Block should be enclosed by merged exon interval"
        );
    }
}
