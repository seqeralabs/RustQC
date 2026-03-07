//! GTF annotation file parser.
//!
//! Parses GTF/GFF2 format files to extract gene, transcript, and exon
//! information. Computes effective gene lengths as the total number of
//! non-overlapping exon bases (matching featureCounts behavior).
//! Also extracts transcript-level structure (exon blocks, CDS ranges) for
//! use by RSeQC-equivalent tools.

use anyhow::{Context, Result};
use indexmap::IndexMap;
use std::collections::HashMap;
use std::io::BufRead;

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

/// Represents a single transcript with its exon blocks and optional CDS range.
///
/// Transcript-level structure is needed by RSeQC-equivalent tools that require
/// per-transcript exon/intron/CDS information (e.g. read_distribution,
/// junction_annotation, inner_distance).
#[derive(Debug, Clone)]
pub struct Transcript {
    /// Transcript identifier (from transcript_id attribute)
    pub transcript_id: String,
    /// Chromosome/contig name
    pub chrom: String,
    /// Transcript start (minimum of all exon starts, 1-based inclusive)
    pub start: u64,
    /// Transcript end (maximum of all exon ends, 1-based inclusive)
    pub end: u64,
    /// Strand: '+', '-', or '.'
    pub strand: char,
    /// Sorted exon blocks as (start, end) pairs (1-based inclusive)
    pub exons: Vec<(u64, u64)>,
    /// CDS start (minimum of all CDS feature starts, 1-based inclusive).
    /// `None` for non-coding transcripts.
    pub cds_start: Option<u64>,
    /// CDS end (maximum of all CDS feature ends, 1-based inclusive).
    /// `None` for non-coding transcripts.
    pub cds_end: Option<u64>,
}

/// Represents a gene with all its exons and transcripts.
#[derive(Debug, Clone)]
pub struct Gene {
    /// Gene identifier (from gene_id attribute)
    pub gene_id: String,
    /// Chromosome/contig name (populated from GTF, used by downstream tools)
    pub chrom: String,
    /// Gene start (minimum of all exon starts, 1-based)
    pub start: u64,
    /// Gene end (maximum of all exon ends, 1-based)
    pub end: u64,
    /// Strand: '+', '-', or '.' (populated from GTF, used by downstream tools)
    pub strand: char,
    /// All exons belonging to this gene (flattened across transcripts)
    pub exons: Vec<Exon>,
    /// Effective length: total non-overlapping exon bases
    pub effective_length: u64,
    /// Additional extracted attributes (e.g. gene_biotype, gene_name).
    ///
    /// Keys are attribute names, values are the attribute values from the GTF.
    /// Only attributes requested via `extra_attributes` in [`parse_gtf`] are stored.
    pub attributes: HashMap<String, String>,
    /// Per-transcript exon/CDS structure, built from transcript_id grouping.
    ///
    /// This is populated during GTF parsing and provides the transcript-level
    /// detail needed by RSeQC-equivalent tools. May be empty if the GTF lacks
    /// transcript_id attributes.
    pub transcripts: Vec<Transcript>,
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

/// Intermediate accumulator for building a transcript from parsed GTF features.
struct TranscriptBuilder {
    /// Transcript identifier
    transcript_id: String,
    /// Chromosome (from first exon)
    chrom: String,
    /// Strand (from first exon)
    strand: char,
    /// Exon intervals: (start, end), 1-based inclusive
    exons: Vec<(u64, u64)>,
    /// CDS intervals: (start, end), 1-based inclusive
    cds: Vec<(u64, u64)>,
    /// start_codon intervals: (start, end), 1-based inclusive
    start_codons: Vec<(u64, u64)>,
    /// stop_codon intervals: (start, end), 1-based inclusive
    stop_codons: Vec<(u64, u64)>,
}

/// Parse a GTF file and return a map of gene_id -> Gene.
///
/// Extracts all exon and CDS features, groups them by gene_id, and builds
/// transcript-level structures by grouping exons/CDS by transcript_id.
/// Computes effective gene length from non-overlapping exon bases.
/// Returns genes in the order they are first encountered in the GTF.
///
/// # Arguments
/// * `path` - Path to the GTF annotation file
/// * `extra_attributes` - Additional GTF attribute names to extract (e.g. `["gene_biotype", "gene_name"]`)
///
/// # Returns
/// An IndexMap preserving insertion order of gene_id -> Gene
pub fn parse_gtf(path: &str, extra_attributes: &[String]) -> Result<IndexMap<String, Gene>> {
    let reader = crate::io::open_reader(path)
        .with_context(|| format!("Failed to open GTF file: {}", path))?;

    let mut genes: IndexMap<String, Gene> = IndexMap::new();

    // Accumulate transcript-level data: (gene_id, transcript_id) -> TranscriptBuilder
    // Use IndexMap to preserve insertion order within each gene
    let mut tx_builders: IndexMap<(String, String), TranscriptBuilder> = IndexMap::new();

    for line in reader.lines() {
        let line = line.context("Failed to read line from GTF file")?;

        // Skip comments and empty lines
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let feature_type = fields[2];

        // We care about exon, CDS, start_codon, and stop_codon features
        if feature_type != "exon"
            && feature_type != "CDS"
            && feature_type != "start_codon"
            && feature_type != "stop_codon"
        {
            continue;
        }

        let chrom = fields[0].to_string();
        let start: u64 = fields[3]
            .parse()
            .with_context(|| format!("Invalid start position: {}", fields[3]))?;
        let end: u64 = fields[4]
            .parse()
            .with_context(|| format!("Invalid end position: {}", fields[4]))?;
        anyhow::ensure!(
            start <= end,
            "Malformed GTF: start ({}) > end ({}) for {} on {}",
            start,
            end,
            feature_type,
            fields[0]
        );
        let strand = fields[6].chars().next().unwrap_or('.');

        let attr_str = fields[8];
        let gene_id = match get_attribute(attr_str, "gene_id") {
            Some(id) => id,
            None => continue, // Skip features without gene_id
        };

        // For exon features, update the gene-level data (unchanged from original)
        if feature_type == "exon" {
            let exon = Exon {
                chrom: chrom.clone(),
                start,
                end,
                strand,
            };

            // Extract extra attributes from the first exon encountered for this gene
            genes
                .entry(gene_id.clone())
                .and_modify(|gene| {
                    gene.start = gene.start.min(start);
                    gene.end = gene.end.max(end);
                    gene.exons.push(exon.clone());
                })
                .or_insert_with(|| {
                    let mut attrs = HashMap::new();
                    for attr_name in extra_attributes {
                        if let Some(val) = get_attribute(attr_str, attr_name) {
                            attrs.insert(attr_name.clone(), val);
                        }
                    }
                    Gene {
                        gene_id: gene_id.clone(),
                        chrom: chrom.clone(),
                        start,
                        end,
                        strand,
                        exons: vec![exon],
                        effective_length: 0, // computed later
                        attributes: attrs,
                        transcripts: Vec::new(), // populated later
                    }
                });
        }

        // Accumulate transcript-level data for both exon and CDS features.
        // Features without transcript_id are grouped under a synthetic key.
        let transcript_id = get_attribute(attr_str, "transcript_id")
            .unwrap_or_else(|| format!("{}__no_tx", gene_id));

        let key = (gene_id.clone(), transcript_id.clone());
        let builder = tx_builders.entry(key).or_insert_with(|| TranscriptBuilder {
            transcript_id,
            chrom: chrom.clone(),
            strand,
            exons: Vec::new(),
            cds: Vec::new(),
            start_codons: Vec::new(),
            stop_codons: Vec::new(),
        });

        match feature_type {
            "exon" => builder.exons.push((start, end)),
            "CDS" => builder.cds.push((start, end)),
            "start_codon" => builder.start_codons.push((start, end)),
            "stop_codon" => builder.stop_codons.push((start, end)),
            _ => unreachable!(),
        }
    }

    // Warn if no genes were extracted — likely indicates a format problem
    if genes.is_empty() {
        log::warn!(
            "No genes extracted from GTF file '{}'. Check that the file is in GTF (not GFF3) \
             format and contains exon features with gene_id attributes.",
            path
        );
    }

    // Build Transcript structs from accumulated data and attach to genes
    for ((gene_id, _), mut builder) in tx_builders {
        // Skip transcript builders that have no exons (CDS-only entries are useless)
        if builder.exons.is_empty() {
            continue;
        }

        // Sort exon blocks by start position
        builder.exons.sort_unstable();

        let tx_start = builder.exons.first().map(|e| e.0).unwrap_or(0);
        let tx_end = builder.exons.last().map(|e| e.1).unwrap_or(0);

        // Compute CDS range (thickStart/thickEnd) using start_codon/stop_codon
        // positions with strand-aware logic matching the Perl gtf2bed script
        // used by upstream nf-core/rnaseq.
        //
        // The Perl script stores the start of the last start_codon line and the
        // end of the last stop_codon line (overwriting on each occurrence).
        // For split codons spanning exon boundaries, multiple GTF lines exist
        // for the same codon; the Perl script keeps the last line's values.
        //
        // For + strand: thick_start = last_sc_start, thick_end = last_stc_end
        // For - strand: the Perl script swaps and adjusts:
        //   thick_start = last_stc_end - 2, thick_end = last_sc_start + 2
        //   (this reconstructs the codon span from a single position,
        //    correct for non-split codons and matching upstream behavior)
        // Missing codons fall back to transcript boundaries (tx_start/tx_end).
        // Transcripts with CDS features but no codon features get tx_start/tx_end.
        let (cds_start, cds_end) = if builder.cds.is_empty() {
            (None, None)
        } else {
            let has_start_codon = !builder.start_codons.is_empty();
            let has_stop_codon = !builder.stop_codons.is_empty();

            let thick_start;
            let thick_end;

            if has_start_codon || has_stop_codon {
                // Use the last codon line seen (matching Perl overwrite semantics).
                // For split codons across exon boundaries, multiple GTF lines exist;
                // the Perl script simply overwrites with each line, keeping the last.
                let last_start_codon_start = builder.start_codons.last().map(|c| c.0);
                let last_stop_codon_end = builder.stop_codons.last().map(|c| c.1);

                match builder.strand {
                    '+' => {
                        // + strand: thick_start = last start_codon start
                        //           thick_end   = last stop_codon end
                        thick_start = last_start_codon_start.unwrap_or(tx_start);
                        thick_end = last_stop_codon_end.unwrap_or(tx_end);
                    }
                    '-' => {
                        // - strand: Perl swaps then adjusts:
                        //   thick_start = last_stop_codon_end - 2
                        //   thick_end   = last_start_codon_start + 2
                        thick_start = last_stop_codon_end
                            .map(|e| e.saturating_sub(2))
                            .unwrap_or(tx_start);
                        thick_end = last_start_codon_start.map(|s| s + 2).unwrap_or(tx_end);
                    }
                    _ => {
                        // Unknown strand: fall back to transcript boundaries
                        thick_start = tx_start;
                        thick_end = tx_end;
                    }
                }
            } else {
                // CDS features but no codon features: use transcript boundaries
                thick_start = tx_start;
                thick_end = tx_end;
            }

            (Some(thick_start), Some(thick_end))
        };

        let transcript = Transcript {
            transcript_id: builder.transcript_id,
            chrom: builder.chrom,
            start: tx_start,
            end: tx_end,
            strand: builder.strand,
            exons: builder.exons,
            cds_start,
            cds_end,
        };

        // Attach to the parent gene (which must exist since we only created
        // transcript builders for features that also had exon lines)
        if let Some(gene) = genes.get_mut(&gene_id) {
            gene.transcripts.push(transcript);
        }
    }

    // Compute effective lengths for all genes
    for gene in genes.values_mut() {
        gene.effective_length = compute_non_overlapping_length(&gene.exons);
    }

    Ok(genes)
}

/// Check whether a given attribute name exists in a GTF file.
///
/// Scans up to `max_lines` data lines (non-comment, non-empty) and returns
/// `true` if the attribute is found in at least one exon feature.
pub fn attribute_exists_in_gtf(path: &str, attribute_name: &str, max_lines: usize) -> bool {
    let reader = match crate::io::open_reader(path) {
        Ok(r) => r,
        Err(_) => return false,
    };
    let mut checked = 0;
    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(_) => break,
        };
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 || fields[2] != "exon" {
            continue;
        }
        if get_attribute(fields[8], attribute_name).is_some() {
            return true;
        }
        checked += 1;
        if checked >= max_lines {
            break;
        }
    }
    false
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

    /// Helper to create a temporary GTF file for testing.
    fn write_temp_gtf(content: &str) -> (std::path::PathBuf, std::fs::File) {
        use std::io::Write;
        use std::sync::atomic::{AtomicU64, Ordering};
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        let dir = std::env::temp_dir();
        let id = COUNTER.fetch_add(1, Ordering::Relaxed);
        let path = dir.join(format!(
            "rustqc_test_{:?}_{}.gtf",
            std::thread::current().id(),
            id
        ));
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f.flush().unwrap();
        (path, f)
    }

    #[test]
    fn test_parse_gtf_transcripts_basic() {
        // T1 has CDS + start_codon + stop_codon → precise CDS range
        // T2 has no CDS features → non-coding
        let (path, _f) = write_temp_gtf(
            "\
chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\texon\t300\t400\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\texon\t100\t250\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\";\n\
chr1\ttest\tCDS\t120\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tCDS\t300\t380\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tstart_codon\t120\t122\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tstop_codon\t381\t383\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
        );

        let genes = parse_gtf(path.to_str().unwrap(), &[]).unwrap();
        assert_eq!(genes.len(), 1);
        let gene = &genes["G1"];

        // Gene-level: 3 exons total, spanning 100-400
        assert_eq!(gene.exons.len(), 3);
        assert_eq!(gene.start, 100);
        assert_eq!(gene.end, 400);
        assert_eq!(gene.strand, '+');

        // Should have 2 transcripts
        assert_eq!(gene.transcripts.len(), 2);

        let t1 = gene
            .transcripts
            .iter()
            .find(|t| t.transcript_id == "T1")
            .unwrap();
        assert_eq!(t1.exons, vec![(100, 200), (300, 400)]);
        assert_eq!(t1.start, 100);
        assert_eq!(t1.end, 400);
        // + strand with start_codon and stop_codon:
        // thick_start = start_codon start = 120
        // thick_end = stop_codon end = 383
        assert_eq!(t1.cds_start, Some(120));
        assert_eq!(t1.cds_end, Some(383));

        let t2 = gene
            .transcripts
            .iter()
            .find(|t| t.transcript_id == "T2")
            .unwrap();
        assert_eq!(t2.exons, vec![(100, 250)]);
        assert_eq!(t2.start, 100);
        assert_eq!(t2.end, 250);
        assert_eq!(t2.cds_start, None); // non-coding
        assert_eq!(t2.cds_end, None);
    }

    #[test]
    fn test_parse_gtf_no_transcript_id() {
        // Exons without transcript_id should be grouped under a synthetic key
        let (path, _f) = write_temp_gtf(
            "\
chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id \"G1\";\n\
chr1\ttest\texon\t300\t400\t.\t+\t.\tgene_id \"G1\";\n",
        );

        let genes = parse_gtf(path.to_str().unwrap(), &[]).unwrap();
        let gene = &genes["G1"];
        // Should have exactly 1 transcript with synthetic ID
        assert_eq!(gene.transcripts.len(), 1);
        assert_eq!(gene.transcripts[0].transcript_id, "G1__no_tx");
        assert_eq!(gene.transcripts[0].exons, vec![(100, 200), (300, 400)]);
    }

    #[test]
    fn test_parse_gtf_cds_only_no_exon() {
        // CDS features without corresponding exons should not create transcripts
        let (path, _f) = write_temp_gtf(
            "\
chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tCDS\t500\t600\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T_orphan\";\n",
        );

        let genes = parse_gtf(path.to_str().unwrap(), &[]).unwrap();
        let gene = &genes["G1"];
        // T_orphan should be skipped because it has no exons
        assert_eq!(gene.transcripts.len(), 1);
        assert_eq!(gene.transcripts[0].transcript_id, "T1");
    }

    #[test]
    fn test_parse_gtf_multi_gene_transcripts() {
        // G1/T1: CDS but no codon features → falls back to tx_start/tx_end
        // G2/T2: non-coding (no CDS)
        let (path, _f) = write_temp_gtf(
            "\
chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tCDS\t120\t180\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr2\ttest\texon\t500\t600\t.\t-\t.\tgene_id \"G2\"; transcript_id \"T2\";\n\
chr2\ttest\texon\t700\t800\t.\t-\t.\tgene_id \"G2\"; transcript_id \"T2\";\n",
        );

        let genes = parse_gtf(path.to_str().unwrap(), &[]).unwrap();
        assert_eq!(genes.len(), 2);

        let g1 = &genes["G1"];
        assert_eq!(g1.transcripts.len(), 1);
        // CDS features but no codon features → tx_start=100, tx_end=200
        assert_eq!(g1.transcripts[0].cds_start, Some(100));
        assert_eq!(g1.transcripts[0].cds_end, Some(200));

        let g2 = &genes["G2"];
        assert_eq!(g2.transcripts.len(), 1);
        assert_eq!(g2.transcripts[0].exons, vec![(500, 600), (700, 800)]);
        assert_eq!(g2.transcripts[0].strand, '-');
        assert_eq!(g2.transcripts[0].cds_start, None);
    }

    #[test]
    fn test_stop_codon_extends_cds_range_forward_strand() {
        // On + strand with start_codon and stop_codon:
        // thick_start = start_codon start, thick_end = stop_codon end
        let (path, _f) = write_temp_gtf(
            "\
chr1\ttest\texon\t1000\t2000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\texon\t3000\t4000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tCDS\t1200\t2000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tCDS\t3000\t3500\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tstart_codon\t1200\t1202\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tstop_codon\t3501\t3503\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
        );

        let genes = parse_gtf(path.to_str().unwrap(), &[]).unwrap();
        let t1 = &genes["G1"].transcripts[0];
        // + strand: thick_start = start_codon start = 1200
        //           thick_end = stop_codon end = 3503
        assert_eq!(t1.cds_start, Some(1200));
        assert_eq!(t1.cds_end, Some(3503));
    }

    #[test]
    fn test_stop_codon_extends_cds_range_reverse_strand() {
        // On - strand with start_codon and stop_codon (non-split):
        // Perl gtf2bed: swap then adjust → thick_start = stc_end - 2, thick_end = sc_start + 2
        // For a contiguous 3-base stop_codon [1497,1499]: 1499 - 2 = 1497
        // For a contiguous 3-base start_codon [3798,3800]: 3798 + 2 = 3800
        let (path, _f) = write_temp_gtf(
            "\
chr1\ttest\texon\t1000\t2000\t.\t-\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\texon\t3000\t4000\t.\t-\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tCDS\t1500\t2000\t.\t-\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tCDS\t3000\t3800\t.\t-\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tstart_codon\t3798\t3800\t.\t-\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tstop_codon\t1497\t1499\t.\t-\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
        );

        let genes = parse_gtf(path.to_str().unwrap(), &[]).unwrap();
        let t1 = &genes["G1"].transcripts[0];
        // - strand: thick_start = last_stc_end - 2 = 1499 - 2 = 1497
        //           thick_end   = last_sc_start + 2 = 3798 + 2 = 3800
        assert_eq!(t1.cds_start, Some(1497));
        assert_eq!(t1.cds_end, Some(3800));
    }

    #[test]
    fn test_split_stop_codon_reverse_strand() {
        // On - strand with a stop codon split across two exons:
        //   stop_codon line 1: [524, 524] (1 base, exon 2)
        //   stop_codon line 2: [101, 102] (2 bases, exon 3)
        //
        // The Perl gtf2bed script overwrites $sc{id}[1] with each stop_codon
        // line, keeping the last end seen (102). Then for - strand:
        //   thick_start = 102 - 2 = 100 (1-based)
        //
        // Using min/max would incorrectly give thick_start = min(524,101) - 1
        // which is wrong. We must match the Perl overwrite behavior.
        let (path, _f) = write_temp_gtf(
            "\
chr1\ttest\texon\t100\t200\t.\t-\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\texon\t400\t600\t.\t-\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\texon\t800\t1000\t.\t-\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tCDS\t110\t200\t.\t-\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tCDS\t400\t600\t.\t-\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tCDS\t800\t900\t.\t-\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tstart_codon\t898\t900\t.\t-\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tstop_codon\t524\t524\t.\t-\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tstop_codon\t101\t102\t.\t-\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
        );

        let genes = parse_gtf(path.to_str().unwrap(), &[]).unwrap();
        let t1 = &genes["G1"].transcripts[0];
        // - strand: last stop_codon end = 102, last start_codon start = 898
        // thick_start = 102 - 2 = 100 (Perl: stc_end - 2)
        // thick_end = 898 + 2 = 900 (Perl: sc_start + 2)
        assert_eq!(t1.cds_start, Some(100));
        assert_eq!(t1.cds_end, Some(900));
    }

    #[test]
    fn test_split_start_codon_forward_strand() {
        // On + strand with a start codon split across two exons:
        //   start_codon line 1: [198, 199] (2 bases, exon 1)
        //   start_codon line 2: [401, 401] (1 base, exon 2)
        //
        // Perl keeps the last start_codon start seen (401).
        // For + strand: thick_start = last_sc_start = 401
        let (path, _f) = write_temp_gtf(
            "\
chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\texon\t400\t600\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tCDS\t198\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tCDS\t400\t550\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tstart_codon\t198\t199\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tstart_codon\t401\t401\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tstop_codon\t551\t553\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
        );

        let genes = parse_gtf(path.to_str().unwrap(), &[]).unwrap();
        let t1 = &genes["G1"].transcripts[0];
        // + strand: thick_start = last start_codon start = 401
        //           thick_end   = last stop_codon end = 553
        assert_eq!(t1.cds_start, Some(401));
        assert_eq!(t1.cds_end, Some(553));
    }

    #[test]
    fn test_stop_codon_without_cds_does_not_create_coding() {
        // stop_codon without CDS features should not make a non-coding transcript coding
        let (path, _f) = write_temp_gtf(
            "\
chr1\ttest\texon\t1000\t2000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
chr1\ttest\tstop_codon\t1500\t1502\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
        );

        let genes = parse_gtf(path.to_str().unwrap(), &[]).unwrap();
        let t1 = &genes["G1"].transcripts[0];
        assert_eq!(t1.cds_start, None);
        assert_eq!(t1.cds_end, None);
    }
}
