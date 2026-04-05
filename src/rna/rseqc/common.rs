//! Shared utilities for RSeQC tool reimplementations.
//!
//! Common functions used by multiple RSeQC tools, including CIGAR-based intron
//! extraction and GTF-based junction building.

use std::collections::{HashMap, HashSet};

use indexmap::IndexMap;
use log::debug;
use noodles_sam::alignment::record::cigar::op::Kind;

use crate::gtf::Gene;

// ===================================================================
// CIGAR intron extraction
// ===================================================================

/// Extract intron intervals from a CIGAR string.
///
/// Matches RSeQC's `fetch_intron()` behavior: only CIGAR `N` operations produce
/// intron intervals. Soft clips do NOT advance position (unlike `fetch_exon()`).
/// `=` and `X` operations are ignored (do not advance position) — matching the
/// original Python bug.
///
/// # Arguments
/// * `start_pos` - Alignment start position (0-based, from BAM record).
/// * `cigar_ops` - CIGAR operations from the BAM record.
///
/// # Returns
/// Vector of `(intron_start, intron_end)` tuples (0-based coordinates).
pub fn fetch_introns(
    start_pos: u64,
    cigar_ops: &dyn noodles_sam::alignment::record::Cigar,
) -> std::io::Result<Vec<(u64, u64)>> {
    let mut pos = start_pos;
    let mut introns = Vec::new();

    for op_result in cigar_ops.iter() {
        let op = op_result?;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                pos += op.len() as u64;
            }
            Kind::Insertion | Kind::SoftClip => {}
            Kind::Deletion => {
                pos += op.len() as u64;
            }
            Kind::Skip => {
                // N: intron!
                let intron_start = pos;
                let intron_end = pos + op.len() as u64;
                introns.push((intron_start, intron_end));
                pos = intron_end;
            }
            Kind::HardClip => {}
            Kind::Pad => {} // Padding - no position change
        }
    }

    Ok(introns)
}

// ===================================================================
// Reference junction data structures
// ===================================================================

/// Reference annotation: sets of known intron start and end positions per chromosome.
///
/// Used by `junction_annotation` for classifying junctions as annotated,
/// partial novel, or complete novel.
#[derive(Debug, Default)]
pub struct ReferenceJunctions {
    /// Known intron start positions per chromosome (uppercased).
    pub intron_starts: HashMap<String, HashSet<u64>>,
    /// Known intron end positions per chromosome (uppercased).
    pub intron_ends: HashMap<String, HashSet<u64>>,
}

/// Classify a junction as annotated, partial novel, or complete novel
/// based on whether its start and end positions appear in the reference.
pub fn classify_junction(
    chrom: &str,
    intron_start: u64,
    intron_end: u64,
    reference: &ReferenceJunctions,
) -> super::junction_annotation::JunctionClass {
    let start_known = reference
        .intron_starts
        .get(chrom)
        .is_some_and(|s| s.contains(&intron_start));
    let end_known = reference
        .intron_ends
        .get(chrom)
        .is_some_and(|s| s.contains(&intron_end));

    use super::junction_annotation::JunctionClass;
    match (start_known, end_known) {
        (true, true) => JunctionClass::Annotated,
        (false, false) => JunctionClass::CompleteNovel,
        _ => JunctionClass::PartialNovel,
    }
}

/// Set of known junction keys in `"CHROM:start-end"` format.
///
/// Used by `junction_saturation` for classifying junctions.
#[derive(Debug, Default)]
pub struct KnownJunctionSet {
    /// Junction keys in `"CHROM:start-end"` format (uppercased chromosome).
    pub junctions: HashSet<String>,
}

// ===================================================================
// GTF-based junction extraction
// ===================================================================

/// Build reference junctions from GTF gene annotations.
///
/// Extracts intron boundaries from transcript-level exon structures. For each
/// transcript, introns are the gaps between sorted exon blocks. Exon coordinates
/// are converted from GTF 1-based inclusive to 0-based half-open.
///
/// # Arguments
/// * `genes` - Parsed GTF gene annotations with transcript-level data.
///
/// # Returns
/// A `ReferenceJunctions` with intron start and end position sets per chromosome.
pub fn build_reference_junctions_from_genes(genes: &IndexMap<String, Gene>) -> ReferenceJunctions {
    let mut result = ReferenceJunctions::default();
    let mut transcript_count = 0u64;

    for gene in genes.values() {
        for tx in &gene.transcripts {
            if tx.exons.len() <= 1 {
                continue;
            }

            let chrom = tx.chrom.to_uppercase();
            let starts_set = result.intron_starts.entry(chrom.clone()).or_default();
            let ends_set = result.intron_ends.entry(chrom).or_default();

            // Exons are already sorted by start in the Transcript struct (GTF 1-based inclusive).
            // Convert to 0-based half-open: start_0based = start - 1, end_0based = end
            // Introns are gaps between consecutive exons.
            for i in 0..tx.exons.len() - 1 {
                let intron_start = tx.exons[i].1; // exon end in GTF is inclusive, so 0-based exclusive = end
                let intron_end = tx.exons[i + 1].0 - 1; // next exon start in GTF is 1-based, so 0-based = start - 1
                starts_set.insert(intron_start);
                ends_set.insert(intron_end);
            }

            transcript_count += 1;
        }
    }

    debug!(
        "Extracted reference junctions from {} multi-exon transcripts (GTF)",
        transcript_count
    );

    result
}

/// Build known junction set from GTF gene annotations.
///
/// Returns junction keys in `"CHROM:start-end"` format, matching the BED-based
/// format for compatibility with junction saturation analysis.
///
/// # Arguments
/// * `genes` - Parsed GTF gene annotations with transcript-level data.
pub fn build_known_junctions_from_genes(genes: &IndexMap<String, Gene>) -> KnownJunctionSet {
    let mut result = KnownJunctionSet::default();

    for gene in genes.values() {
        for tx in &gene.transcripts {
            if tx.exons.len() <= 1 {
                continue;
            }

            let chrom = tx.chrom.to_uppercase();

            // Same coordinate conversion as build_reference_junctions_from_genes
            for i in 0..tx.exons.len() - 1 {
                let intron_start = tx.exons[i].1;
                let intron_end = tx.exons[i + 1].0 - 1;
                let key = format!("{chrom}:{intron_start}-{intron_end}");
                result.junctions.insert(key);
            }
        }
    }

    debug!(
        "Extracted {} known splice junctions from GTF",
        result.junctions.len()
    );

    result
}

// ===================================================================
// Unit tests
// ===================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gtf::{Exon, Transcript};
    use std::collections::HashMap;

    #[test]
    fn test_fetch_introns_simple() {
        use noodles_sam::alignment::{record::cigar::op::Kind, record_buf::Cigar};
        // 50M500N50M — one intron at position 100+50=150 to 150+500=650
        let ops = vec![
            sam::alignment::record::cigar::Op::new(Kind::Match, 50),
            sam::alignment::record::cigar::Op::new(Kind::Skip, 500),
            sam::alignment::record::cigar::Op::new(Kind::Match, 50),
        ];
        let cigar = Cigar::from(ops);
        let introns = fetch_introns(100, &cigar).unwrap();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], (150, 650));
    }

    #[test]
    fn test_fetch_introns_multiple() {
        use noodles_sam::alignment::{record::cigar::op::Kind, record_buf::Cigar};
        // 10M500N20M300N10M — two introns
        let ops = vec![
            sam::alignment::record::cigar::Op::new(Kind::Match, 10),
            sam::alignment::record::cigar::Op::new(Kind::Skip, 500),
            sam::alignment::record::cigar::Op::new(Kind::Match, 20),
            sam::alignment::record::cigar::Op::new(Kind::Skip, 300),
            sam::alignment::record::cigar::Op::new(Kind::Match, 10),
        ];
        let cigar = Cigar::from(ops);
        let introns = fetch_introns(100, &cigar).unwrap();
        assert_eq!(introns.len(), 2);
        assert_eq!(introns[0], (110, 610));
        assert_eq!(introns[1], (630, 930));
    }

    #[test]
    fn test_fetch_introns_with_deletions() {
        use noodles_sam::alignment::{record::cigar::op::Kind, record_buf::Cigar};
        // 10M5D10M500N10M
        let ops = vec![
            sam::alignment::record::cigar::Op::new(Kind::Match, 10),
            sam::alignment::record::cigar::Op::new(Kind::Deletion, 5),
            sam::alignment::record::cigar::Op::new(Kind::Match, 10),
            sam::alignment::record::cigar::Op::new(Kind::Skip, 500),
            sam::alignment::record::cigar::Op::new(Kind::Match, 10),
        ];
        let cigar = Cigar::from(ops);
        let introns = fetch_introns(100, &cigar).unwrap();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], (125, 625)); // 100+10+5+10=125
    }

    #[test]
    fn test_fetch_introns_no_introns() {
        use noodles_sam::alignment::{record::cigar::op::Kind, record_buf::Cigar};
        let ops = vec![sam::alignment::record::cigar::Op::new(Kind::Match, 100)];
        let cigar = Cigar::from(ops);
        let introns = fetch_introns(100, &cigar).unwrap();
        assert!(introns.is_empty());
    }

    #[test]
    fn test_fetch_introns_soft_clip_no_advance() {
        use noodles_sam::alignment::{record::cigar::op::Kind, record_buf::Cigar};
        // 5S50M500N50M — soft clip should NOT advance position
        let ops = vec![
            sam::alignment::record::cigar::Op::new(Kind::SoftClip, 5),
            sam::alignment::record::cigar::Op::new(Kind::Match, 50),
            sam::alignment::record::cigar::Op::new(Kind::Skip, 500),
            sam::alignment::record::cigar::Op::new(Kind::Match, 50),
        ];
        let cigar = Cigar::from(ops);
        let introns = fetch_introns(100, &cigar).unwrap();
        assert_eq!(introns.len(), 1);
        assert_eq!(introns[0], (150, 650));
    }

    /// Helper to create a Gene with transcripts for testing.
    fn make_gene_with_transcripts(
        gene_id: &str,
        chrom: &str,
        transcripts: Vec<Transcript>,
    ) -> Gene {
        let exons: Vec<Exon> = transcripts
            .iter()
            .flat_map(|tx| {
                tx.exons.iter().map(|&(start, end)| Exon {
                    chrom: chrom.to_string(),
                    start,
                    end,
                    strand: '+',
                })
            })
            .collect();

        let start = exons.iter().map(|e| e.start).min().unwrap_or(0);
        let end = exons.iter().map(|e| e.end).max().unwrap_or(0);

        Gene {
            gene_id: gene_id.to_string(),
            chrom: chrom.to_string(),
            start,
            end,
            strand: '+',
            exons,
            effective_length: 0,
            attributes: HashMap::new(),
            transcripts,
        }
    }

    #[test]
    fn test_build_reference_junctions_from_genes() {
        // Gene with one transcript: 3 exons -> 2 introns
        // GTF coords (1-based inclusive): exon1=[101,200], exon2=[301,400], exon3=[501,600]
        // 0-based half-open: exon1=[100,200), exon2=[300,400), exon3=[500,600)
        // Introns: [200, 300) and [400, 500)
        let tx = Transcript {
            transcript_id: "TX1".to_string(),

            chrom: "chr1".to_string(),
            start: 101,
            end: 600,
            strand: '+',
            exons: vec![(101, 200), (301, 400), (501, 600)],
            cds_start: None,
            cds_end: None,
        };

        let gene = make_gene_with_transcripts("G1", "chr1", vec![tx]);
        let mut genes = IndexMap::new();
        genes.insert("G1".to_string(), gene);

        let ref_junctions = build_reference_junctions_from_genes(&genes);

        let chr1_starts = ref_junctions.intron_starts.get("CHR1").unwrap();
        let chr1_ends = ref_junctions.intron_ends.get("CHR1").unwrap();

        assert_eq!(chr1_starts.len(), 2);
        assert!(chr1_starts.contains(&200)); // exon1 end (GTF inclusive) = 0-based exclusive
        assert!(chr1_starts.contains(&400)); // exon2 end

        assert_eq!(chr1_ends.len(), 2);
        assert!(chr1_ends.contains(&300)); // exon2 start (GTF 1-based) - 1 = 300
        assert!(chr1_ends.contains(&500)); // exon3 start - 1
    }

    #[test]
    fn test_build_known_junctions_from_genes() {
        let tx = Transcript {
            transcript_id: "TX1".to_string(),

            chrom: "chr1".to_string(),
            start: 101,
            end: 600,
            strand: '+',
            exons: vec![(101, 200), (301, 400), (501, 600)],
            cds_start: None,
            cds_end: None,
        };

        let gene = make_gene_with_transcripts("G1", "chr1", vec![tx]);
        let mut genes = IndexMap::new();
        genes.insert("G1".to_string(), gene);

        let known = build_known_junctions_from_genes(&genes);

        assert_eq!(known.junctions.len(), 2);
        assert!(known.junctions.contains("CHR1:200-300"));
        assert!(known.junctions.contains("CHR1:400-500"));
    }

    #[test]
    fn test_build_junctions_skips_single_exon() {
        let tx = Transcript {
            transcript_id: "TX1".to_string(),

            chrom: "chr1".to_string(),
            start: 101,
            end: 200,
            strand: '+',
            exons: vec![(101, 200)],
            cds_start: None,
            cds_end: None,
        };

        let gene = make_gene_with_transcripts("G1", "chr1", vec![tx]);
        let mut genes = IndexMap::new();
        genes.insert("G1".to_string(), gene);

        let ref_junctions = build_reference_junctions_from_genes(&genes);
        assert!(ref_junctions.intron_starts.is_empty());

        let known = build_known_junctions_from_genes(&genes);
        assert!(known.junctions.is_empty());
    }

    #[test]
    fn test_build_junctions_multiple_transcripts() {
        // Two transcripts with different splicing patterns
        let tx1 = Transcript {
            transcript_id: "TX1".to_string(),

            chrom: "chr1".to_string(),
            start: 101,
            end: 600,
            strand: '+',
            exons: vec![(101, 200), (301, 600)],
            cds_start: None,
            cds_end: None,
        };

        let tx2 = Transcript {
            transcript_id: "TX2".to_string(),

            chrom: "chr1".to_string(),
            start: 101,
            end: 600,
            strand: '+',
            exons: vec![(101, 250), (401, 600)],
            cds_start: None,
            cds_end: None,
        };

        let gene = make_gene_with_transcripts("G1", "chr1", vec![tx1, tx2]);
        let mut genes = IndexMap::new();
        genes.insert("G1".to_string(), gene);

        let ref_junctions = build_reference_junctions_from_genes(&genes);
        let chr1_starts = ref_junctions.intron_starts.get("CHR1").unwrap();

        // TX1 intron: start=200, TX2 intron: start=250
        assert!(chr1_starts.contains(&200));
        assert!(chr1_starts.contains(&250));
    }
}
