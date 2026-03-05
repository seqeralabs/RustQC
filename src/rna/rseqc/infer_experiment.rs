//! Infer library strandedness from RNA-Seq alignments.
//!
//! Reimplementation of RSeQC's `infer_experiment.py`. Samples reads that overlap
//! gene models (BED12) and determines the fraction consistent with each strand
//! protocol.

use crate::gtf::Gene;
use anyhow::{Context, Result};
use indexmap::IndexMap;
use log::info;
use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, Write};
use std::path::Path;

// ============================================================================
// BED12 gene model
// ============================================================================

/// A genomic interval with strand information, parsed from a BED12 file.
#[derive(Debug, Clone)]
struct TranscriptInterval {
    /// 0-based start position (BED start)
    start: u64,
    /// 0-based exclusive end position (BED end)
    end: u64,
    /// Strand: `b'+'` or `b'-'`
    strand: u8,
}

/// Per-chromosome collection of transcript intervals, sorted by start position.
#[derive(Debug, Default)]
pub struct GeneModel {
    /// Map from chromosome name to sorted list of transcript intervals.
    intervals: HashMap<String, Vec<TranscriptInterval>>,
}

impl GeneModel {
    /// Build a gene model from parsed GTF gene annotations.
    ///
    /// Creates one interval per **transcript** (not per gene) to match
    /// upstream RSeQC's BED12-based model. Using gene-level spans would
    /// include introns and inter-transcript gaps, causing reads in
    /// non-transcript regions to dilute the strand signal.
    ///
    /// GTF coordinates are 1-based inclusive; BED coordinates are 0-based
    /// half-open. Conversion: BED start = GTF start - 1, BED end = GTF end
    /// (since GTF end is inclusive and BED end is exclusive, the numeric value
    /// is the same).
    pub fn from_genes(genes: &IndexMap<String, Gene>) -> Self {
        let mut model = GeneModel::default();
        let mut count: u64 = 0;

        for gene in genes.values() {
            // Use transcript-level intervals to match RSeQC's BED12 approach
            if gene.transcripts.is_empty() {
                // Fallback: use gene span if no transcripts were parsed
                let strand = match gene.strand {
                    '+' => b'+',
                    '-' => b'-',
                    _ => continue,
                };
                let start = gene.start.saturating_sub(1);
                let end = gene.end;
                model
                    .intervals
                    .entry(gene.chrom.clone())
                    .or_default()
                    .push(TranscriptInterval { start, end, strand });
                count += 1;
            } else {
                for tx in &gene.transcripts {
                    let strand = match tx.strand {
                        '+' => b'+',
                        '-' => b'-',
                        _ => continue,
                    };
                    // Convert 1-based inclusive GTF to 0-based half-open BED
                    let start = tx.start.saturating_sub(1);
                    let end = tx.end;
                    model
                        .intervals
                        .entry(tx.chrom.clone())
                        .or_default()
                        .push(TranscriptInterval { start, end, strand });
                    count += 1;
                }
            }
        }

        // Sort intervals by start position for binary search
        for intervals in model.intervals.values_mut() {
            intervals.sort_by_key(|iv| iv.start);
        }

        info!("Loaded {} transcript intervals from GTF annotation", count);
        model
    }

    /// Load a BED12 (or BED6+) file into the gene model.
    ///
    /// Extracts columns: chrom (0), start (1), end (2), strand (5).
    pub fn from_bed(path: &str) -> Result<Self> {
        let reader = crate::io::open_reader(path)
            .with_context(|| format!("Failed to open BED file: {}", path))?;

        let mut model = GeneModel::default();
        let mut count: u64 = 0;

        for line in reader.lines() {
            let line = line.context("Failed to read BED line")?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') || line.starts_with("track") {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 6 {
                continue; // Need at least 6 columns for strand
            }

            let chrom = fields[0].to_string();
            let start: u64 = fields[1]
                .parse()
                .with_context(|| format!("Invalid BED start: {}", fields[1]))?;
            let end: u64 = fields[2]
                .parse()
                .with_context(|| format!("Invalid BED end: {}", fields[2]))?;
            let strand = match fields[5] {
                "+" => b'+',
                "-" => b'-',
                _ => continue, // Skip unknown strand
            };

            model
                .intervals
                .entry(chrom)
                .or_default()
                .push(TranscriptInterval { start, end, strand });
            count += 1;
        }

        // Sort intervals by start position for binary search
        for intervals in model.intervals.values_mut() {
            intervals.sort_by_key(|iv| iv.start);
        }

        info!("Loaded {} transcripts from BED file", count);
        Ok(model)
    }

    /// Find all strands of transcripts overlapping the query interval [qstart, qend).
    ///
    /// Returns a set of strand bytes found among overlapping transcripts.
    pub fn find_strands(&self, chrom: &str, qstart: u64, qend: u64) -> Vec<u8> {
        let mut strands = Vec::new();
        if let Some(intervals) = self.intervals.get(chrom) {
            // Binary search to find the first interval that could overlap
            // An interval overlaps if interval.start < qend && interval.end > qstart
            let idx = intervals.partition_point(|iv| iv.start < qend);
            // Check all intervals from 0..idx that could overlap
            // (their start < qend, but we also need end > qstart)
            // We scan backwards from idx since intervals are sorted by start
            // Actually, we need to scan from the beginning because intervals with
            // small start could still overlap if their end > qstart
            for iv in intervals.iter().take(idx) {
                if iv.end > qstart && !strands.contains(&iv.strand) {
                    strands.push(iv.strand);
                }
            }
        }
        strands
    }
}

// ============================================================================
// Strandedness inference
// ============================================================================

/// Result of strandedness inference.
#[derive(Debug)]
pub struct InferExperimentResult {
    /// Number of usable reads sampled.
    pub total_sampled: u64,
    /// Library type: "PairEnd", "SingleEnd", or "Mixture".
    pub library_type: String,
    /// Fraction of reads that failed to determine strandedness.
    pub frac_failed: f64,
    /// Fraction explained by protocol 1 (PE: 1++,1--,2+-,2-+; SE: ++,--)
    pub frac_protocol1: f64,
    /// Fraction explained by protocol 2 (PE: 1+-,1-+,2++,2--; SE: +-,-+)
    pub frac_protocol2: f64,
}

// ============================================================================
// Output
// ============================================================================

/// Write infer_experiment results to a file in RSeQC-compatible format.
///
/// # Arguments
///
/// * `result` - The inference result.
/// * `output_path` - Path to write the output file.
pub fn write_infer_experiment<P: AsRef<Path>>(
    result: &InferExperimentResult,
    output_path: P,
) -> Result<()> {
    let output_path = output_path.as_ref();
    let mut writer = fs::File::create(output_path)
        .with_context(|| format!("Failed to create output file: {}", output_path.display()))?;

    match result.library_type.as_str() {
        "PairEnd" => {
            writeln!(writer)?;
            writeln!(writer)?;
            writeln!(writer, "This is PairEnd Data")?;
            writeln!(
                writer,
                "Fraction of reads failed to determine: {:.4}",
                result.frac_failed
            )?;
            writeln!(
                writer,
                "Fraction of reads explained by \"1++,1--,2+-,2-+\": {:.4}",
                result.frac_protocol1
            )?;
            writeln!(
                writer,
                "Fraction of reads explained by \"1+-,1-+,2++,2--\": {:.4}",
                result.frac_protocol2
            )?;
        }
        "SingleEnd" => {
            writeln!(writer)?;
            writeln!(writer)?;
            writeln!(writer, "This is SingleEnd Data")?;
            writeln!(
                writer,
                "Fraction of reads failed to determine: {:.4}",
                result.frac_failed
            )?;
            writeln!(
                writer,
                "Fraction of reads explained by \"++,--\": {:.4}",
                result.frac_protocol1
            )?;
            writeln!(
                writer,
                "Fraction of reads explained by \"+-,-+\": {:.4}",
                result.frac_protocol2
            )?;
        }
        _ => {
            writeln!(writer)?;
            writeln!(writer)?;
            writeln!(writer, "Unknown Data type")?;
        }
    }

    Ok(())
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gtf::{Exon, Gene};
    use indexmap::IndexMap;

    /// Helper to create a simple Gene for testing.
    fn make_gene(gene_id: &str, chrom: &str, start: u64, end: u64, strand: char) -> Gene {
        Gene {
            gene_id: gene_id.to_string(),
            chrom: chrom.to_string(),
            start,
            end,
            strand,
            exons: vec![Exon {
                chrom: chrom.to_string(),
                start,
                end,
                strand,
            }],
            effective_length: end - start + 1,
            attributes: HashMap::new(),
            transcripts: Vec::new(),
        }
    }

    #[test]
    fn test_from_genes_basic() {
        let mut genes = IndexMap::new();
        genes.insert(
            "GENE1".to_string(),
            make_gene("GENE1", "chr1", 100, 500, '+'),
        );
        genes.insert(
            "GENE2".to_string(),
            make_gene("GENE2", "chr1", 1000, 2000, '-'),
        );

        let model = GeneModel::from_genes(&genes);

        // Should have chr1 with 2 intervals
        assert_eq!(model.intervals.len(), 1);
        let chr1 = model.intervals.get("chr1").unwrap();
        assert_eq!(chr1.len(), 2);

        // First interval: GTF 100-500 -> BED 99-500
        assert_eq!(chr1[0].start, 99);
        assert_eq!(chr1[0].end, 500);
        assert_eq!(chr1[0].strand, b'+');

        // Second interval: GTF 1000-2000 -> BED 999-2000
        assert_eq!(chr1[1].start, 999);
        assert_eq!(chr1[1].end, 2000);
        assert_eq!(chr1[1].strand, b'-');
    }

    #[test]
    fn test_from_genes_skips_unknown_strand() {
        let mut genes = IndexMap::new();
        genes.insert(
            "GENE1".to_string(),
            make_gene("GENE1", "chr1", 100, 500, '.'),
        );
        genes.insert(
            "GENE2".to_string(),
            make_gene("GENE2", "chr1", 600, 800, '+'),
        );

        let model = GeneModel::from_genes(&genes);
        let chr1 = model.intervals.get("chr1").unwrap();
        assert_eq!(chr1.len(), 1); // Only GENE2, GENE1 skipped
        assert_eq!(chr1[0].strand, b'+');
    }

    #[test]
    fn test_from_genes_find_strands() {
        let mut genes = IndexMap::new();
        genes.insert(
            "GENE1".to_string(),
            make_gene("GENE1", "chr1", 100, 500, '+'),
        );
        genes.insert(
            "GENE2".to_string(),
            make_gene("GENE2", "chr1", 300, 800, '-'),
        );

        let model = GeneModel::from_genes(&genes);

        // Query overlapping both genes (BED coords: 350-400)
        let strands = model.find_strands("chr1", 350, 400);
        assert_eq!(strands.len(), 2);
        assert!(strands.contains(&b'+'));
        assert!(strands.contains(&b'-'));

        // Query overlapping only GENE2 (BED coords: 550-600)
        let strands = model.find_strands("chr1", 550, 600);
        assert_eq!(strands, vec![b'-']);
    }

    #[test]
    fn test_bed_loading() {
        // Test that BED loading works with the small test file if available
        let bed_path = "benchmark/input/small/chr6.bed";
        if Path::new(bed_path).exists() {
            let model = GeneModel::from_bed(bed_path).unwrap();
            assert!(!model.intervals.is_empty(), "Should have loaded intervals");
            // chr6.bed should have 8443 transcripts
            let total: usize = model.intervals.values().map(|v| v.len()).sum();
            assert!(total > 0, "Should have transcripts");
        }
    }

    #[test]
    fn test_find_strands_no_overlap() {
        let model = GeneModel::default();
        let strands = model.find_strands("chr1", 100, 200);
        assert!(strands.is_empty());
    }

    #[test]
    fn test_find_strands_single_overlap() {
        let mut model = GeneModel::default();
        model.intervals.insert(
            "chr1".to_string(),
            vec![TranscriptInterval {
                start: 100,
                end: 500,
                strand: b'+',
            }],
        );
        let strands = model.find_strands("chr1", 200, 300);
        assert_eq!(strands, vec![b'+']);
    }

    #[test]
    fn test_find_strands_both_strands() {
        let mut model = GeneModel::default();
        model.intervals.insert(
            "chr1".to_string(),
            vec![
                TranscriptInterval {
                    start: 100,
                    end: 500,
                    strand: b'+',
                },
                TranscriptInterval {
                    start: 200,
                    end: 600,
                    strand: b'-',
                },
            ],
        );
        let strands = model.find_strands("chr1", 250, 350);
        assert_eq!(strands.len(), 2);
        assert!(strands.contains(&b'+'));
        assert!(strands.contains(&b'-'));
    }
}
