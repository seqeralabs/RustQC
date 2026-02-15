//! Samtools stats SN-section compatible output.
//!
//! Produces the Summary Numbers (SN) section from `samtools stats`, which is
//! the section parsed by MultiQC for key alignment statistics.

use std::path::Path;

use anyhow::{Context, Result};
use log::info;

use super::bam_stat::BamStatResult;

// ============================================================================
// Output formatting
// ============================================================================

/// Format a float value to 1 decimal place.
fn fmt_1dp(v: f64) -> String {
    format!("{v:.1}")
}

/// Format a float value in scientific notation matching samtools output.
///
/// samtools uses C's `%e` format which gives 6 decimal places with a
/// zero-padded two-digit exponent (e.g. `1.000000e+00`). Rust's `{:.6e}`
/// does not pad the exponent, so we reformat manually.
fn fmt_sci(v: f64) -> String {
    let s = format!("{v:.6e}");
    // Rust gives e.g. "0.000000e0" or "1.200000e-3" — we need "e+00" / "e-03"
    if let Some(pos) = s.find('e') {
        let (mantissa, exp_part) = s.split_at(pos);
        let exp_str = &exp_part[1..]; // skip 'e'
        let (sign, digits) = if let Some(rest) = exp_str.strip_prefix('-') {
            ("-", rest)
        } else if let Some(rest) = exp_str.strip_prefix('+') {
            ("+", rest)
        } else {
            ("+", exp_str)
        };
        let exp_num: i32 = digits.parse().unwrap_or(0);
        format!("{mantissa}e{sign}{exp_num:02}")
    } else {
        s
    }
}

/// Write samtools stats SN-section compatible output.
///
/// The output starts with a samtools-compatible comment header (important for
/// MultiQC detection) and then emits `SN\t<key>:\t<value>\t# <comment>` lines.
///
/// # Arguments
/// * `result` - The computed BAM statistics
/// * `output_path` - Path to write the stats file
pub fn write_stats(result: &BamStatResult, output_path: &Path) -> Result<()> {
    use std::io::Write;

    let mut out = std::fs::File::create(output_path)
        .with_context(|| format!("Failed to create stats file: {}", output_path.display()))?;

    // Header comment — MultiQC detects "This file was produced by samtools stats"
    writeln!(out, "# This file was produced by samtools stats and RustQC")?;
    writeln!(
        out,
        "# The command line was: rustqc rna (samtools stats compatible output)"
    )?;

    // Derived values
    let raw_total = result.primary_count; // primary = non-secondary, non-supplementary
    let filtered = result.qc_failed;

    // "sequences" = primary reads that passed QC
    let sequences = raw_total.saturating_sub(filtered);

    // Average lengths
    let avg_len = if raw_total > 0 {
        result.total_len as f64 / raw_total as f64
    } else {
        0.0
    };
    let avg_first = if result.first_fragments > 0 {
        result.total_first_fragment_len as f64 / result.first_fragments as f64
    } else {
        0.0
    };
    let avg_last = if result.last_fragments > 0 {
        result.total_last_fragment_len as f64 / result.last_fragments as f64
    } else {
        0.0
    };

    // Average quality
    let avg_quality = if result.quality_count > 0 {
        result.quality_sum / result.quality_count as f64
    } else {
        0.0
    };

    // Insert size stats
    let avg_insert = if result.insert_size_count > 0 {
        result.insert_size_sum / result.insert_size_count as f64
    } else {
        0.0
    };
    let insert_sd = if result.insert_size_count > 1 {
        let n = result.insert_size_count as f64;
        let mean = result.insert_size_sum / n;
        let variance = (result.insert_size_sq_sum / n) - (mean * mean);
        if variance > 0.0 {
            variance.sqrt()
        } else {
            0.0
        }
    } else {
        0.0
    };

    // Error rate = mismatches / bases_mapped_cigar
    let error_rate = if result.bases_mapped_cigar > 0 {
        result.mismatches as f64 / result.bases_mapped_cigar as f64
    } else {
        0.0
    };

    // Write SN lines
    sn(
        &mut out,
        "raw total sequences:",
        raw_total,
        "excluding supplementary and secondary reads",
    )?;
    sn_no_comment(&mut out, "filtered sequences:", filtered)?;
    sn_no_comment(&mut out, "sequences:", sequences)?;
    sn(&mut out, "is sorted:", 1_u64, "sorted by coordinate")?;
    sn_no_comment(&mut out, "1st fragments:", result.first_fragments)?;
    sn_no_comment(&mut out, "last fragments:", result.last_fragments)?;

    // Reads mapped
    sn_no_comment(&mut out, "reads mapped:", result.primary_mapped)?;
    sn(
        &mut out,
        "reads mapped and paired:",
        result.reads_mapped_and_paired,
        "paired-end technology bit set + both mates mapped",
    )?;
    sn_no_comment(&mut out, "reads unmapped:", result.unmapped)?;
    sn(
        &mut out,
        "reads properly paired:",
        result.properly_paired,
        "proper-pair bit set",
    )?;
    sn(
        &mut out,
        "reads paired:",
        result.paired_flagstat,
        "paired-end technology bit set",
    )?;
    sn(
        &mut out,
        "reads duplicated:",
        result.primary_duplicates,
        "PCR or optical duplicate bit set",
    )?;
    sn(&mut out, "reads MQ0:", result.reads_mq0, "mapped and MQ=0")?;

    // Quality and length stats
    sn_no_comment(&mut out, "reads QC failed:", filtered)?;
    sn_no_comment(
        &mut out,
        "non-primary alignments:",
        result.secondary + result.supplementary,
    )?;
    sn_no_comment(&mut out, "supplementary alignments:", result.supplementary)?;
    sn(
        &mut out,
        "total length:",
        result.total_len,
        "ignores clipping",
    )?;
    sn(
        &mut out,
        "total first fragment length:",
        result.total_first_fragment_len,
        "ignores clipping",
    )?;
    sn(
        &mut out,
        "total last fragment length:",
        result.total_last_fragment_len,
        "ignores clipping",
    )?;
    sn(
        &mut out,
        "bases mapped:",
        result.bases_mapped,
        "ignores clipping",
    )?;
    sn(
        &mut out,
        "bases mapped (cigar):",
        result.bases_mapped_cigar,
        "more accurate",
    )?;
    sn_no_comment(&mut out, "bases trimmed:", 0_u64)?;
    sn_no_comment(&mut out, "bases duplicated:", result.bases_duplicated)?;
    sn(&mut out, "mismatches:", result.mismatches, "from NM fields")?;
    sn_str(
        &mut out,
        "error rate:",
        &fmt_sci(error_rate),
        "mismatches / bases mapped (cigar)",
    )?;
    sn_no_comment(&mut out, "average length:", avg_len as u64)?;
    sn_no_comment(&mut out, "average first fragment length:", avg_first as u64)?;
    sn_no_comment(&mut out, "average last fragment length:", avg_last as u64)?;
    sn_no_comment(&mut out, "maximum length:", result.max_len)?;
    sn_no_comment(
        &mut out,
        "maximum first fragment length:",
        result.max_first_fragment_len,
    )?;
    sn_no_comment(
        &mut out,
        "maximum last fragment length:",
        result.max_last_fragment_len,
    )?;
    sn_no_comment(&mut out, "average quality:", fmt_1dp(avg_quality))?;
    sn_no_comment(&mut out, "insert size average:", fmt_1dp(avg_insert))?;
    sn_no_comment(
        &mut out,
        "insert size standard deviation:",
        fmt_1dp(insert_sd),
    )?;
    sn_no_comment(&mut out, "inward oriented pairs:", result.inward_pairs)?;
    sn_no_comment(&mut out, "outward oriented pairs:", result.outward_pairs)?;
    sn_no_comment(
        &mut out,
        "pairs with other orientation:",
        result.other_orientation,
    )?;
    sn_no_comment(
        &mut out,
        "pairs on different chromosomes:",
        result.mate_diff_chr_mapq5,
    )?;

    // Percentage stats
    sn_no_comment(
        &mut out,
        "percentage of properly paired reads (%):",
        fmt_1dp(if raw_total > 0 {
            100.0 * result.properly_paired as f64 / raw_total as f64
        } else {
            0.0
        }),
    )?;

    info!("Wrote samtools stats output to {}", output_path.display());
    Ok(())
}

/// Write a single SN line with a numeric value.
fn sn<W: std::io::Write, V: std::fmt::Display>(
    out: &mut W,
    key: &str,
    value: V,
    comment: &str,
) -> Result<()> {
    writeln!(out, "SN\t{key}\t{value}\t# {comment}")?;
    Ok(())
}

/// Write a single SN line with a pre-formatted string value.
fn sn_str<W: std::io::Write>(out: &mut W, key: &str, value: &str, comment: &str) -> Result<()> {
    writeln!(out, "SN\t{key}\t{value}\t# {comment}")?;
    Ok(())
}

/// Write a single SN line without a comment (matching samtools for many lines).
fn sn_no_comment<W: std::io::Write, V: std::fmt::Display>(
    out: &mut W,
    key: &str,
    value: V,
) -> Result<()> {
    writeln!(out, "SN\t{key}\t{value}")?;
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
    fn test_stats_sn_format() {
        let mut reader =
            bam::Reader::from_path("tests/data/test.bam").expect("Failed to open test.bam");
        let mut accum = BamStatAccum::default();
        let mut record = bam::Record::new();

        while let Some(res) = reader.read(&mut record) {
            res.expect("Error reading BAM record");
            accum.process_read(&record, 30);
        }

        let result = accum.into_result();
        let tmp_path = std::env::temp_dir().join("rustqc_test_stats.txt");
        write_stats(&result, &tmp_path).expect("Failed to write stats");

        let mut contents = String::new();
        std::fs::File::open(&tmp_path)
            .unwrap()
            .read_to_string(&mut contents)
            .unwrap();
        let _ = std::fs::remove_file(&tmp_path);

        // Verify MultiQC detection header
        assert!(
            contents.contains("This file was produced by samtools stats"),
            "Missing MultiQC detection header"
        );

        // Verify SN line format
        for line in contents.lines() {
            if line.starts_with("SN\t") {
                let cols: Vec<&str> = line.splitn(4, '\t').collect();
                assert!(
                    cols.len() >= 3,
                    "SN line should have at least 3 tab-separated fields: {line}"
                );
                assert!(
                    cols[1].ends_with(':'),
                    "SN key should end with colon: {}",
                    cols[1]
                );
                // Lines with comments should start with "# "
                if cols.len() == 4 {
                    assert!(
                        cols[3].starts_with("# "),
                        "SN comment should start with '# ': {}",
                        cols[3]
                    );
                }
            }
        }
    }
}
