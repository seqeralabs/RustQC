//! featureCounts-compatible output file generation.
//!
//! Produces the same output files as Subread's featureCounts:
//! - Main counts file (`.featureCounts.tsv`): per-gene read counts with annotation
//! - Summary file (`.featureCounts.tsv.summary`): assignment statistics
//! - Biotype counts file: counts aggregated by a GTF attribute (e.g. `gene_biotype`)
//! - MultiQC biotype bargraph file: biotype counts formatted for MultiQC

use crate::gtf::Gene;
use crate::rna::dupradar::counting::CountResult;
use anyhow::Result;
use indexmap::IndexMap;
use std::io::Write;
use std::path::Path;

/// Write the featureCounts-compatible main counts file.
///
/// Produces a tab-separated file matching featureCounts output format:
/// - Line 1: comment with program version and command
/// - Line 2: header with column names
/// - Lines 3+: one row per gene with annotation columns and count
///
/// # Columns
/// `Geneid`, `Chr`, `Start`, `End`, `Strand`, `Length`, `<bam_filename>`
pub fn write_counts_file(
    path: &Path,
    genes: &IndexMap<String, Gene>,
    counts: &CountResult,
    bam_path: &str,
    command_line: &str,
) -> Result<()> {
    let file = std::fs::File::create(path)?;
    let mut w = std::io::BufWriter::new(file);

    // Line 1: comment header with program info (matches featureCounts format)
    writeln!(
        w,
        "# Program:RustQC v{}; Command:\"{}\"",
        env!("CARGO_PKG_VERSION"),
        command_line
    )?;

    // Line 2: column header
    let bam_name = Path::new(bam_path)
        .file_name()
        .map(|n| n.to_string_lossy().to_string())
        .unwrap_or_else(|| bam_path.to_string());
    writeln!(w, "Geneid\tChr\tStart\tEnd\tStrand\tLength\t{}", bam_name)?;

    // Data rows
    for (gene_id, gene) in genes.iter() {
        // Skip genes with no exons
        if gene.exons.is_empty() {
            continue;
        }

        // Build semicolon-separated annotation fields from exons
        let chrs: Vec<&str> = gene.exons.iter().map(|e| e.chrom.as_str()).collect();
        let starts: Vec<String> = gene.exons.iter().map(|e| e.start.to_string()).collect();
        let ends: Vec<String> = gene.exons.iter().map(|e| e.end.to_string()).collect();
        let strands: Vec<String> = gene.exons.iter().map(|e| e.strand.to_string()).collect();

        // Use per-read counts for featureCounts output (matches featureCounts -p behavior
        // where each read is counted independently, not as fragments)
        let count = counts
            .gene_counts
            .get(gene_id)
            .map(|gc| gc.fc_reads)
            .unwrap_or(0);

        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            gene_id,
            chrs.join(";"),
            starts.join(";"),
            ends.join(";"),
            strands.join(";"),
            gene.effective_length,
            count,
        )?;
    }

    Ok(())
}

/// Write the featureCounts-compatible summary file.
///
/// Produces a tab-separated file with read assignment statistics.
/// The status categories match featureCounts output.
pub fn write_summary_file(path: &Path, counts: &CountResult, bam_path: &str) -> Result<()> {
    let file = std::fs::File::create(path)?;
    let mut w = std::io::BufWriter::new(file);

    let bam_name = Path::new(bam_path)
        .file_name()
        .map(|n| n.to_string_lossy().to_string())
        .unwrap_or_else(|| bam_path.to_string());

    writeln!(w, "Status\t{}", bam_name)?;
    // Use per-read featureCounts stats (matching featureCounts -p behavior where
    // each read is counted independently, not as fragments). The fc_* fields track
    // per-read assignment with multi-mapping reads separated out.
    writeln!(w, "Assigned\t{}", counts.fc_assigned)?;
    writeln!(w, "Unassigned_Unmapped\t{}", counts.fc_unmapped)?;
    writeln!(w, "Unassigned_Read_Type\t0")?;
    writeln!(w, "Unassigned_Singleton\t0")?;
    writeln!(w, "Unassigned_MappingQuality\t0")?;
    writeln!(w, "Unassigned_Chimera\t0")?;
    writeln!(w, "Unassigned_FragmentLength\t0")?;
    writeln!(w, "Unassigned_Duplicate\t0")?;
    writeln!(w, "Unassigned_MultiMapping\t{}", counts.fc_multimapping)?;
    writeln!(w, "Unassigned_Secondary\t0")?;
    writeln!(w, "Unassigned_NonSplit\t0")?;
    writeln!(w, "Unassigned_NoFeatures\t{}", counts.fc_no_features)?;
    writeln!(w, "Unassigned_Overlapping_Length\t0")?;
    writeln!(w, "Unassigned_Ambiguity\t{}", counts.fc_ambiguous)?;

    Ok(())
}

/// Aggregate counts by a given GTF attribute (e.g. `gene_biotype`).
///
/// Returns an ordered map of attribute_value -> total count.
/// Genes missing the attribute are grouped under `"unknown"`.
pub fn aggregate_biotype_counts(
    genes: &IndexMap<String, Gene>,
    counts: &CountResult,
    attribute_name: &str,
) -> IndexMap<String, u64> {
    let mut biotype_counts: IndexMap<String, u64> = IndexMap::new();

    for (gene_id, gene) in genes.iter() {
        let biotype = gene
            .attributes
            .get(attribute_name)
            .cloned()
            .unwrap_or_else(|| "unknown".to_string());

        // Use per-read counts for featureCounts biotype output
        let count = counts
            .gene_counts
            .get(gene_id)
            .map(|gc| gc.fc_reads)
            .unwrap_or(0);

        *biotype_counts.entry(biotype).or_insert(0) += count;
    }

    // Sort by count descending for readability
    let mut pairs: Vec<(String, u64)> = biotype_counts.into_iter().collect();
    pairs.sort_by(|a, b| b.1.cmp(&a.1));
    pairs.into_iter().collect()
}

/// Write biotype counts to a plain TSV file.
///
/// Two columns: biotype name and count, tab-separated.
pub fn write_biotype_counts(path: &Path, biotype_counts: &IndexMap<String, u64>) -> Result<()> {
    let file = std::fs::File::create(path)?;
    let mut w = std::io::BufWriter::new(file);

    writeln!(w, "Biotype\tCount")?;
    for (biotype, count) in biotype_counts.iter() {
        writeln!(w, "{}\t{}", biotype, count)?;
    }

    Ok(())
}

/// Write biotype counts formatted for MultiQC custom content (bargraph).
///
/// This matches the format used by the nf-core/rnaseq pipeline's
/// `MULTIQC_CUSTOM_BIOTYPE` process.
pub fn write_biotype_counts_mqc(path: &Path, biotype_counts: &IndexMap<String, u64>) -> Result<()> {
    let file = std::fs::File::create(path)?;
    let mut w = std::io::BufWriter::new(file);

    // MultiQC YAML header
    writeln!(w, "# id: 'biotype_counts'")?;
    writeln!(w, "# section_name: 'Biotype Counts'")?;
    writeln!(
        w,
        "# description: \"shows reads overlapping genomic features of different biotypes,\""
    )?;
    writeln!(w, "# plot_type: 'bargraph'")?;
    writeln!(w, "# anchor: 'featurecounts_biotype'")?;
    writeln!(w, "# pconfig:")?;
    writeln!(w, "#     id: \"featurecounts_biotype_plot\"")?;
    writeln!(w, "#     title: \"featureCounts: Biotypes\"")?;
    writeln!(w, "#     xlab: \"# Reads\"")?;
    writeln!(w, "#     cpswitch_counts_label: \"Number of Reads\"")?;

    // Data: biotype\tcount (no header row — MultiQC custom content format)
    for (biotype, count) in biotype_counts.iter() {
        writeln!(w, "{}\t{}", biotype, count)?;
    }

    Ok(())
}

/// Write biotype rRNA percentage for MultiQC general stats.
///
/// Extracts the rRNA count and computes the percentage of total assigned reads.
pub fn write_biotype_rrna_mqc(
    path: &Path,
    biotype_counts: &IndexMap<String, u64>,
    total_assigned: u64,
    sample_name: &str,
) -> Result<()> {
    let file = std::fs::File::create(path)?;
    let mut w = std::io::BufWriter::new(file);

    let rrna_count = biotype_counts.get("rRNA").copied().unwrap_or(0);
    let rrna_pct = if total_assigned > 0 {
        rrna_count as f64 / total_assigned as f64 * 100.0
    } else {
        0.0
    };

    // MultiQC YAML header for general stats
    writeln!(w, "# id: 'biotype_counts_rrna'")?;
    writeln!(w, "# plot_type: 'generalstats'")?;
    writeln!(w, "# pconfig:")?;
    writeln!(w, "#     percent_rRNA:")?;
    writeln!(w, "#         title: '% rRNA'")?;
    writeln!(w, "#         namespace: 'Biotype Counts'")?;
    writeln!(
        w,
        "#         description: 'Percentage of reads assigned to rRNA genes'"
    )?;
    writeln!(w, "#         max: 100")?;
    writeln!(w, "#         min: 0")?;
    writeln!(w, "#         format: '{{:.2f}}'")?;
    writeln!(w, "Sample\tpercent_rRNA")?;
    writeln!(w, "{}\t{:.4}", sample_name, rrna_pct)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rna::dupradar::counting::GeneCounts;
    use std::collections::HashMap;

    fn make_test_gene(gene_id: &str, biotype: Option<&str>) -> Gene {
        let mut attributes = HashMap::new();
        if let Some(bt) = biotype {
            attributes.insert("gene_biotype".to_string(), bt.to_string());
        }
        Gene {
            gene_id: gene_id.to_string(),
            chrom: "chr1".to_string(),
            start: 100,
            end: 200,
            strand: '+',
            exons: vec![crate::gtf::Exon {
                chrom: "chr1".to_string(),
                start: 100,
                end: 200,
                strand: '+',
            }],
            effective_length: 101,
            attributes,
            transcripts: Vec::new(),
        }
    }

    #[test]
    fn test_aggregate_biotype_counts() {
        let mut genes = IndexMap::new();
        genes.insert(
            "gene1".to_string(),
            make_test_gene("gene1", Some("protein_coding")),
        );
        genes.insert(
            "gene2".to_string(),
            make_test_gene("gene2", Some("protein_coding")),
        );
        genes.insert("gene3".to_string(), make_test_gene("gene3", Some("rRNA")));
        genes.insert("gene4".to_string(), make_test_gene("gene4", None));

        let mut gene_counts = IndexMap::new();
        gene_counts.insert(
            "gene1".to_string(),
            GeneCounts {
                all_unique: 100,
                fc_reads: 100,
                ..Default::default()
            },
        );
        gene_counts.insert(
            "gene2".to_string(),
            GeneCounts {
                all_unique: 50,
                fc_reads: 50,
                ..Default::default()
            },
        );
        gene_counts.insert(
            "gene3".to_string(),
            GeneCounts {
                all_unique: 25,
                fc_reads: 25,
                ..Default::default()
            },
        );
        gene_counts.insert(
            "gene4".to_string(),
            GeneCounts {
                all_unique: 10,
                fc_reads: 10,
                ..Default::default()
            },
        );

        let count_result = CountResult {
            gene_counts,
            total_reads_multi_dup: 185,
            total_reads_unique_dup: 185,
            stat_total_reads: 200,
            stat_assigned: 185,
            stat_ambiguous: 5,
            stat_no_features: 10,
            stat_total_fragments: 200,
            fc_assigned: 185,
            fc_ambiguous: 5,
            fc_no_features: 10,
            fc_multimapping: 0,
            fc_unmapped: 0,
            rseqc: None,
            genebody: None,
        };

        let biotypes = aggregate_biotype_counts(&genes, &count_result, "gene_biotype");

        assert_eq!(biotypes.get("protein_coding"), Some(&150));
        assert_eq!(biotypes.get("rRNA"), Some(&25));
        assert_eq!(biotypes.get("unknown"), Some(&10));
    }
}
