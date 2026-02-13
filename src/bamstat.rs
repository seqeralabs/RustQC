//! BAM alignment statistics output.
//!
//! Produces summary statistics for BAM files, compatible with the output format
//! of RSeQC's `bam_stat.py`. All counters are collected during the main counting
//! pass — no additional BAM traversal is required.

use crate::counting::CountResult;
use anyhow::{Context, Result};
use std::io::Write;
use std::path::Path;

/// Write BAM alignment statistics in RSeQC bam_stat.py compatible format.
///
/// # Arguments
///
/// * `path` - Output file path
/// * `result` - Counting results containing all BAM statistics
/// * `bam_path` - Path to the source BAM file (included in header)
pub fn write_bam_stat(path: &Path, result: &CountResult, bam_path: &str) -> Result<()> {
    let file = std::fs::File::create(path)
        .with_context(|| format!("Failed to create BAM stat file: {}", path.display()))?;
    let mut writer = std::io::BufWriter::new(file);

    let non_primary = result.stat_non_primary_hits;
    let mapped = result.stat_total_reads
        - result.stat_unmapped
        - result.stat_qc_failed
        - result.stat_supplementary;
    let unique = mapped.saturating_sub(result.stat_total_multi);

    writeln!(writer, "#=== RustQC bam_stat (RSeQC-compatible) ===")?;
    writeln!(writer, "#")?;
    writeln!(writer, "# Input file: {}", bam_path)?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Total records:                          {}",
        result.stat_total_reads
    )?;
    writeln!(
        writer,
        "QC failed:                              {}",
        result.stat_qc_failed
    )?;
    writeln!(
        writer,
        "Optical/PCR duplicate:                  {}",
        result.stat_total_dup
    )?;
    writeln!(
        writer,
        "Non primary hits:                       {}",
        non_primary
    )?;
    writeln!(
        writer,
        "  Secondary alignments:                 {}",
        result.stat_secondary
    )?;
    writeln!(
        writer,
        "  Supplementary alignments:             {}",
        result.stat_supplementary
    )?;
    writeln!(
        writer,
        "Unmapped reads:                         {}",
        result.stat_unmapped
    )?;
    writeln!(writer)?;
    writeln!(
        writer,
        "mapq < mapq_cut (non-unique):           {}",
        result.stat_total_multi
    )?;
    writeln!(writer, "mapq >= mapq_cut (unique):              {}", unique)?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Read-1:                                 {}",
        result.stat_read1
    )?;
    writeln!(
        writer,
        "Read-2:                                 {}",
        result.stat_read2
    )?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Reads map to '+':                       {}",
        result.stat_plus_strand
    )?;
    writeln!(
        writer,
        "Reads map to '-':                       {}",
        result.stat_minus_strand
    )?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Non-splice reads:                       {}",
        result.stat_non_splice_reads
    )?;
    writeln!(
        writer,
        "Splice reads:                           {}",
        result.stat_splice_reads
    )?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Reads mapped in proper pairs:           {}",
        result.stat_proper_pairs
    )?;
    writeln!(
        writer,
        "Proper-paired reads map to different chrom: {}",
        result.stat_mate_mapped_diff_chr
    )?;
    writeln!(
        writer,
        "Singletons:                             {}",
        result.stat_singletons
    )?;

    writer.flush().context("Failed to flush BAM stat output")?;

    Ok(())
}

/// Write BAM statistics in MultiQC general stats format.
///
/// Produces a TSV file with a YAML header that MultiQC recognises as a
/// `generalstats` data source, providing key mapping metrics in the
/// MultiQC General Statistics table.
///
/// # Arguments
///
/// * `path` - Output file path
/// * `result` - Counting results containing all BAM statistics
/// * `sample_name` - Sample name for the MultiQC row
pub fn write_bam_stat_mqc(path: &Path, result: &CountResult, sample_name: &str) -> Result<()> {
    let file = std::fs::File::create(path)
        .with_context(|| format!("Failed to create BAM stat MultiQC file: {}", path.display()))?;
    let mut writer = std::io::BufWriter::new(file);

    let total = result.stat_total_reads as f64;
    let mapped = (result.stat_total_reads
        - result.stat_unmapped
        - result.stat_qc_failed
        - result.stat_supplementary) as f64;
    let mapped_pct = if total > 0.0 {
        mapped / total * 100.0
    } else {
        0.0
    };
    let dup_pct = if mapped > 0.0 {
        result.stat_total_dup as f64 / mapped * 100.0
    } else {
        0.0
    };
    let proper_pair_pct = if mapped > 0.0 {
        result.stat_proper_pairs as f64 / mapped * 100.0
    } else {
        0.0
    };

    writeln!(writer, "# id: 'rustqc_bamstat'")?;
    writeln!(writer, "# section_name: 'RustQC BAM Stats'")?;
    writeln!(
        writer,
        "# description: 'BAM alignment statistics from RustQC (RSeQC bam_stat replacement)'"
    )?;
    writeln!(writer, "# plot_type: 'generalstats'")?;
    writeln!(writer, "# pconfig:")?;
    writeln!(writer, "#     total_reads:")?;
    writeln!(writer, "#         title: 'Total Reads'")?;
    writeln!(
        writer,
        "#         description: 'Total number of records in the BAM file'"
    )?;
    writeln!(writer, "#         format: '{{:,.0f}}'")?;
    writeln!(writer, "#     mapped_pct:")?;
    writeln!(writer, "#         title: '% Mapped'")?;
    writeln!(
        writer,
        "#         description: 'Percentage of reads that are mapped'"
    )?;
    writeln!(writer, "#         max: 100")?;
    writeln!(writer, "#         min: 0")?;
    writeln!(writer, "#         suffix: '%'")?;
    writeln!(writer, "#         format: '{{:.1f}}'")?;
    writeln!(writer, "#     duplicate_pct:")?;
    writeln!(writer, "#         title: '% Dups'")?;
    writeln!(
        writer,
        "#         description: 'Percentage of mapped reads marked as duplicates'"
    )?;
    writeln!(writer, "#         max: 100")?;
    writeln!(writer, "#         min: 0")?;
    writeln!(writer, "#         suffix: '%'")?;
    writeln!(writer, "#         format: '{{:.1f}}'")?;
    writeln!(writer, "#     proper_pair_pct:")?;
    writeln!(writer, "#         title: '% Proper Pairs'")?;
    writeln!(
        writer,
        "#         description: 'Percentage of mapped reads in proper pairs'"
    )?;
    writeln!(writer, "#         max: 100")?;
    writeln!(writer, "#         min: 0")?;
    writeln!(writer, "#         suffix: '%'")?;
    writeln!(writer, "#         format: '{{:.1f}}'")?;
    writeln!(
        writer,
        "Sample\ttotal_reads\tmapped_pct\tduplicate_pct\tproper_pair_pct"
    )?;
    writeln!(
        writer,
        "{}\t{}\t{:.1}\t{:.1}\t{:.1}",
        sample_name, result.stat_total_reads, mapped_pct, dup_pct, proper_pair_pct
    )?;

    writer
        .flush()
        .context("Failed to flush BAM stat MultiQC output")?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::counting::CountResult;
    use indexmap::IndexMap;
    use std::fs;
    use std::path::PathBuf;

    fn test_output_dir() -> PathBuf {
        let dir = PathBuf::from("tests/output_bamstat_unit");
        let _ = fs::create_dir_all(&dir);
        dir
    }

    fn make_test_result() -> CountResult {
        CountResult {
            gene_counts: IndexMap::new(),
            total_reads_multi_dup: 10,
            total_reads_multi_nodup: 20,
            total_reads_unique_dup: 30,
            total_reads_unique_nodup: 40,
            stat_total_reads: 1000,
            stat_assigned: 500,
            stat_ambiguous: 50,
            stat_no_features: 200,
            stat_total_fragments: 800,
            stat_total_dup: 100,
            stat_total_multi: 150,
            stat_unmapped: 50,
            stat_qc_failed: 10,
            stat_secondary: 5,
            stat_supplementary: 15,
            stat_non_primary_hits: 20,
            stat_proper_pairs: 600,
            stat_read1: 400,
            stat_read2: 400,
            stat_plus_strand: 450,
            stat_minus_strand: 475,
            stat_singletons: 20,
            stat_mate_mapped_diff_chr: 10,
            stat_splice_reads: 300,
            stat_non_splice_reads: 625,
        }
    }

    #[test]
    fn test_write_bam_stat() {
        let dir = test_output_dir();
        let path = dir.join("test.bam_stat.txt");
        let result = make_test_result();

        write_bam_stat(&path, &result, "test.bam").unwrap();

        let content = fs::read_to_string(&path).unwrap();
        assert!(
            content.contains("Total records:"),
            "Missing 'Total records' header"
        );
        assert!(content.contains("1000"), "Missing total records count");
        assert!(content.contains("QC failed:"), "Missing 'QC failed' header");
        assert!(
            content.contains("Unmapped reads:"),
            "Missing 'Unmapped reads' header"
        );
        assert!(
            content.contains("Non primary hits:"),
            "Missing 'Non primary hits' header"
        );
        assert!(
            content.contains("Splice reads:"),
            "Missing 'Splice reads' header"
        );
        assert!(
            content.contains("Reads mapped in proper pairs:"),
            "Missing proper pairs header"
        );

        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_write_bam_stat_mqc() {
        let dir = test_output_dir();
        let path = dir.join("test.bam_stat_mqc.txt");
        let result = make_test_result();

        write_bam_stat_mqc(&path, &result, "sample1").unwrap();

        let content = fs::read_to_string(&path).unwrap();
        assert!(
            content.contains("rustqc_bamstat"),
            "Missing MultiQC section ID"
        );
        assert!(
            content.contains("generalstats"),
            "Missing plot_type generalstats"
        );
        assert!(content.contains("sample1"), "Missing sample name");
        assert!(
            content.contains("1000"),
            "Missing total reads in MultiQC output"
        );

        let _ = fs::remove_dir_all(&dir);
    }
}
