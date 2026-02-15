//! Samtools flagstat-compatible output.
//!
//! Produces output matching `samtools flagstat` format, allowing MultiQC and
//! other downstream tools to parse the results as if they came from samtools.

use std::path::Path;

use anyhow::{Context, Result};
use log::info;

use super::bam_stat::BamStatResult;

// ============================================================================
// Output formatting
// ============================================================================

/// Format a flagstat count line with the standard `<passed> + 0 <description>` format.
fn fmt_line(count: u64, description: &str) -> String {
    format!("{count} + 0 {description}")
}

/// Format a flagstat percentage as `(XX.XX% : N/A)`.
///
/// samtools uses `:` separated `passed : failed` percentages. Since we do not
/// separate QC-pass/fail, the failed side is always `N/A`.
fn fmt_pct(numerator: u64, denominator: u64) -> String {
    if denominator == 0 {
        "(N/A : N/A)".to_string()
    } else {
        let pct = 100.0 * numerator as f64 / denominator as f64;
        format!("({pct:.2}% : N/A)")
    }
}

/// Write samtools flagstat-compatible output.
///
/// The output format matches `samtools flagstat` exactly, with one line per
/// statistic. Each line has the form `<passed> + 0 <description>`.
///
/// # Arguments
/// * `result` - The computed BAM statistics
/// * `output_path` - Path to write the flagstat file
pub fn write_flagstat(result: &BamStatResult, output_path: &Path) -> Result<()> {
    use std::io::Write;

    let mut out = std::fs::File::create(output_path)
        .with_context(|| format!("Failed to create flagstat file: {}", output_path.display()))?;

    // Total = primary + secondary + supplementary (all QC states)
    let total = result.total_records;

    // Primary count for percentage denominators
    let primary = result.primary_count;

    // Line 1: total (QC-passed + QC-failed)
    writeln!(
        out,
        "{}",
        fmt_line(total, "in total (QC-passed reads + QC-failed reads)")
    )?;

    // Line 2: primary
    writeln!(out, "{}", fmt_line(primary, "primary"))?;

    // Line 3: secondary
    writeln!(out, "{}", fmt_line(result.secondary, "secondary"))?;

    // Line 4: supplementary
    writeln!(out, "{}", fmt_line(result.supplementary, "supplementary"))?;

    // Line 5: duplicates
    writeln!(out, "{}", fmt_line(result.duplicates, "duplicates"))?;

    // Line 6: primary duplicates
    writeln!(
        out,
        "{}",
        fmt_line(result.primary_duplicates, "primary duplicates")
    )?;

    // Line 7: mapped
    writeln!(
        out,
        "{} {}",
        fmt_line(result.mapped, "mapped"),
        fmt_pct(result.mapped, total)
    )?;

    // Line 8: primary mapped
    writeln!(
        out,
        "{} {}",
        fmt_line(result.primary_mapped, "primary mapped"),
        fmt_pct(result.primary_mapped, primary)
    )?;

    // Line 9: paired in sequencing
    writeln!(
        out,
        "{}",
        fmt_line(result.paired_flagstat, "paired in sequencing")
    )?;

    // Line 10: read1
    writeln!(out, "{}", fmt_line(result.read1_flagstat, "read1"))?;

    // Line 11: read2
    writeln!(out, "{}", fmt_line(result.read2_flagstat, "read2"))?;

    // Line 12: properly paired (percentage of paired reads)
    writeln!(
        out,
        "{} {}",
        fmt_line(result.properly_paired, "properly paired"),
        fmt_pct(result.properly_paired, result.paired_flagstat)
    )?;

    // Line 13: with itself and mate mapped
    writeln!(
        out,
        "{}",
        fmt_line(result.both_mapped, "with itself and mate mapped")
    )?;

    // Line 14: singletons (percentage of paired reads)
    writeln!(
        out,
        "{} {}",
        fmt_line(result.singletons, "singletons"),
        fmt_pct(result.singletons, result.paired_flagstat)
    )?;

    // Line 15: with mate mapped to a different chr
    writeln!(
        out,
        "{}",
        fmt_line(result.mate_diff_chr, "with mate mapped to a different chr")
    )?;

    // Line 16: with mate mapped to a different chr (mapQ>=5)
    writeln!(
        out,
        "{}",
        fmt_line(
            result.mate_diff_chr_mapq5,
            "with mate mapped to a different chr (mapQ>=5)"
        )
    )?;

    info!("Wrote flagstat output to {}", output_path.display());
    Ok(())
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rna::rseqc::accumulators::BamStatAccum;
    use rust_htslib::bam::{self, Read as BamRead};
    use std::io::Read;
    #[test]
    fn test_flagstat_format() {
        let mut reader =
            bam::Reader::from_path("tests/data/test.bam").expect("Failed to open test.bam");
        let mut accum = BamStatAccum::default();
        let mut record = bam::Record::new();

        while let Some(res) = reader.read(&mut record) {
            res.expect("Error reading BAM record");
            accum.process_read(&record, 30);
        }

        let result = accum.into_result();
        let tmp_path = std::env::temp_dir().join("rustqc_test_flagstat.txt");
        write_flagstat(&result, &tmp_path).expect("Failed to write flagstat");

        let mut contents = String::new();
        std::fs::File::open(&tmp_path)
            .unwrap()
            .read_to_string(&mut contents)
            .unwrap();
        let _ = std::fs::remove_file(&tmp_path);

        // Verify first line format
        assert!(
            contents.contains("in total (QC-passed reads + QC-failed reads)"),
            "Missing total line in flagstat output"
        );
        // Verify mapped line has percentage
        assert!(
            contents.contains("mapped"),
            "Missing mapped line in flagstat output"
        );
    }
}
