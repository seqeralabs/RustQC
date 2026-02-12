//! GTF annotation file parser.
//!
//! Parses GTF/GFF2 format files to extract gene and exon information.
//! Computes effective gene lengths as the total number of non-overlapping
//! exon bases (matching featureCounts behavior).

use anyhow::{Context, Result};
use indexmap::IndexMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Represents a single exon interval.
#[derive(Debug, Clone)]
pub struct Exon {
    /// Chromosome/contig name
    pub chrom: String,
    /// Start position (1-based, inclusive, as in GTF)
    pub start: u64,
    /// End position (1-based, inclusive, as in GTF)
    pub end: u64,
    /// Strand: '+', '-', or '.'
    pub strand: char,
}

/// Represents a gene with all its exons.
#[derive(Debug, Clone)]
pub struct Gene {
    /// Gene identifier (from gene_id attribute)
    pub gene_id: String,
    /// Chromosome/contig name
    #[allow(dead_code)]
    pub chrom: String,
    /// Gene start (minimum of all exon starts, 1-based)
    pub start: u64,
    /// Gene end (maximum of all exon ends, 1-based)
    pub end: u64,
    /// Strand: '+', '-', or '.'
    #[allow(dead_code)]
    pub strand: char,
    /// All exons belonging to this gene
    pub exons: Vec<Exon>,
    /// Effective length: total non-overlapping exon bases
    pub effective_length: u64,
}

/// Parse a GTF attribute string to extract a specific attribute value.
///
/// GTF attributes are semicolon-separated key-value pairs like:
/// `gene_id "ENSG00000223972"; gene_name "DDX11L1";`
fn get_attribute(attributes: &str, key: &str) -> Option<String> {
    for attr in attributes.split(';') {
        let attr = attr.trim();
        if attr.is_empty() {
            continue;
        }
        // Split on first whitespace
        if let Some(pos) = attr.find(|c: char| c.is_whitespace()) {
            let (k, v) = attr.split_at(pos);
            if k.trim() == key {
                // Remove surrounding quotes and whitespace
                let v = v.trim().trim_matches('"');
                return Some(v.to_string());
            }
        }
    }
    None
}

/// Compute the total number of non-overlapping bases across a set of intervals.
///
/// This matches featureCounts behavior for computing effective gene length:
/// merge overlapping exons and sum their lengths.
fn compute_non_overlapping_length(exons: &[Exon]) -> u64 {
    if exons.is_empty() {
        return 0;
    }

    // Collect all intervals as (start, end) and sort by start
    let mut intervals: Vec<(u64, u64)> = exons.iter().map(|e| (e.start, e.end)).collect();
    intervals.sort_unstable();

    // Merge overlapping intervals
    let mut total_bases: u64 = 0;
    let mut current_start = intervals[0].0;
    let mut current_end = intervals[0].1;

    for &(start, end) in &intervals[1..] {
        if start <= current_end + 1 {
            // Overlapping or adjacent, extend
            current_end = current_end.max(end);
        } else {
            // Non-overlapping, count previous interval
            total_bases += current_end - current_start + 1;
            current_start = start;
            current_end = end;
        }
    }
    // Count the last interval
    total_bases += current_end - current_start + 1;

    total_bases
}

/// Parse a GTF file and return a map of gene_id -> Gene.
///
/// Extracts all exon features and groups them by gene_id.
/// Computes effective gene length from non-overlapping exon bases.
/// Returns genes in the order they are first encountered in the GTF.
///
/// # Arguments
/// * `path` - Path to the GTF annotation file
///
/// # Returns
/// An IndexMap preserving insertion order of gene_id -> Gene
pub fn parse_gtf(path: &str) -> Result<IndexMap<String, Gene>> {
    let file = File::open(Path::new(path))
        .with_context(|| format!("Failed to open GTF file: {}", path))?;
    let reader = BufReader::new(file);

    let mut genes: IndexMap<String, Gene> = IndexMap::new();

    for line in reader.lines() {
        let line = line?;

        // Skip comments and empty lines
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        // We only care about exon features (matching featureCounts default)
        if fields[2] != "exon" {
            continue;
        }

        let chrom = fields[0].to_string();
        let start: u64 = fields[3]
            .parse()
            .with_context(|| format!("Invalid start position: {}", fields[3]))?;
        let end: u64 = fields[4]
            .parse()
            .with_context(|| format!("Invalid end position: {}", fields[4]))?;
        let strand = fields[6].chars().next().unwrap_or('.');

        let attributes = fields[8];
        let gene_id = match get_attribute(attributes, "gene_id") {
            Some(id) => id,
            None => continue, // Skip exons without gene_id
        };

        let exon = Exon {
            chrom: chrom.clone(),
            start,
            end,
            strand,
        };

        genes
            .entry(gene_id.clone())
            .and_modify(|gene| {
                gene.start = gene.start.min(start);
                gene.end = gene.end.max(end);
                gene.exons.push(exon.clone());
            })
            .or_insert_with(|| Gene {
                gene_id: gene_id.clone(),
                chrom,
                start,
                end,
                strand,
                exons: vec![exon],
                effective_length: 0, // computed later
            });
    }

    // Compute effective lengths for all genes
    for gene in genes.values_mut() {
        gene.effective_length = compute_non_overlapping_length(&gene.exons);
    }

    Ok(genes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_attribute() {
        let attrs =
            r#"gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; gene_name "DDX11L1";"#;
        assert_eq!(
            get_attribute(attrs, "gene_id"),
            Some("ENSG00000223972".to_string())
        );
        assert_eq!(
            get_attribute(attrs, "gene_name"),
            Some("DDX11L1".to_string())
        );
        assert_eq!(get_attribute(attrs, "missing"), None);
    }

    #[test]
    fn test_non_overlapping_length() {
        // Single exon: 100-200 = 101 bases
        let exons = vec![Exon {
            chrom: "chr1".to_string(),
            start: 100,
            end: 200,
            strand: '+',
        }];
        assert_eq!(compute_non_overlapping_length(&exons), 101);

        // Two non-overlapping exons
        let exons = vec![
            Exon {
                chrom: "chr1".to_string(),
                start: 100,
                end: 200,
                strand: '+',
            },
            Exon {
                chrom: "chr1".to_string(),
                start: 300,
                end: 400,
                strand: '+',
            },
        ];
        assert_eq!(compute_non_overlapping_length(&exons), 202);

        // Two overlapping exons: 100-200, 150-250 -> merged 100-250 = 151
        let exons = vec![
            Exon {
                chrom: "chr1".to_string(),
                start: 100,
                end: 200,
                strand: '+',
            },
            Exon {
                chrom: "chr1".to_string(),
                start: 150,
                end: 250,
                strand: '+',
            },
        ];
        assert_eq!(compute_non_overlapping_length(&exons), 151);

        // Empty
        let exons: Vec<Exon> = vec![];
        assert_eq!(compute_non_overlapping_length(&exons), 0);
    }
}
