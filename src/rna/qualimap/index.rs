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
    pub gene_idx: u32,
    /// Index into the transcript list within the gene.
    #[allow(dead_code)]
    pub transcript_idx: u16,
    /// Exon start in 0-based half-open coordinates.
    pub exon_start: i32,
    /// Exon end in 0-based half-open coordinates.
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
// QualimapIndex — the main index structure
// ============================================================

/// Complete Qualimap annotation index built from GTF genes.
///
/// Contains per-chromosome exon and intron trees, plus transcript metadata.
/// This index is *separate* from the existing dupRadar/featureCounts index
/// because Qualimap uses fundamentally different assignment logic.
pub struct QualimapIndex {
    /// Per-chromosome exon COITrees (key = chromosome name).
    pub exon_trees: HashMap<String, ExonTree>,
    /// Per-chromosome intron COITrees (key = chromosome name).
    pub intron_trees: HashMap<String, IntronTree>,
    /// All transcript metadata, indexed by a flat transcript ID.
    /// The index is `(gene_idx, transcript_idx)` mapped to a flat position.
    pub transcripts: Vec<TranscriptInfo>,
    /// Lookup: gene_idx -> range of transcript indices in `transcripts` vec.
    /// `gene_transcript_ranges[gene_idx] = (start, end)` half-open into `transcripts`.
    pub gene_transcript_ranges: Vec<(u32, u32)>,
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

        for (gene_idx, gene) in genes.values().enumerate() {
            let gene_idx = gene_idx as u32;
            let tx_start_idx = transcripts.len() as u32;

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
        }

        // Build COITrees from collected intervals
        let exon_trees: HashMap<String, ExonTree> = exon_intervals
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
            intron_trees,
            transcripts,
            gene_transcript_ranges,
            num_genes,
        }
    }

    /// Get the exon tree for a chromosome.
    pub fn exon_tree(&self, chrom: &str) -> Option<&ExonTree> {
        self.exon_trees.get(chrom)
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
}
