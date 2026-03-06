//! Samtools stats compatible output.
//!
//! Produces the Summary Numbers (SN) section and all histogram/distribution
//! sections from `samtools stats`, compatible with MultiQC and plot-bamstats.

use std::io::Write;
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

    // Average lengths — use round() to match samtools' printf("%.0f") behavior
    let avg_len = if raw_total > 0 {
        (result.total_len as f64 / raw_total as f64).round()
    } else {
        0.0
    };
    let avg_first = if result.first_fragments > 0 {
        (result.total_first_fragment_len as f64 / result.first_fragments as f64).round()
    } else {
        0.0
    };
    let avg_last = if result.last_fragments > 0 {
        (result.total_last_fragment_len as f64 / result.last_fragments as f64).round()
    } else {
        0.0
    };

    // Average quality
    let avg_quality = if result.quality_count > 0 {
        result.quality_sum / result.quality_count as f64
    } else {
        0.0
    };

    // Insert size stats — 99th percentile truncation matching samtools stats.
    // Each pair is counted once (upstream mate only), capped at 8000.
    // Samtools computes mean/SD from the "main bulk" (up to 99th percentile).
    let (avg_insert, insert_sd) = {
        // Build a flat insert-size histogram from is_hist (total column = index 0)
        let total_count: u64 = result.is_hist.values().map(|v| v[0]).sum();
        if total_count > 0 {
            // Sort insert sizes for percentile computation
            let mut sorted_isizes: Vec<(u64, u64)> =
                result.is_hist.iter().map(|(&k, v)| (k, v[0])).collect();
            sorted_isizes.sort_by_key(|&(k, _)| k);

            // Find 99th percentile cutoff (isize_main_bulk = 0.99)
            let bulk_limit = (total_count as f64 * 0.99).ceil() as u64;
            let mut cumulative = 0u64;
            let mut bulk_sum = 0.0f64;
            let mut bulk_sq_sum = 0.0f64;
            let mut bulk_count = 0u64;

            for &(isize_val, count) in &sorted_isizes {
                let remaining = bulk_limit.saturating_sub(cumulative);
                if remaining == 0 {
                    break;
                }
                let use_count = count.min(remaining);
                let v = isize_val as f64;
                bulk_sum += v * use_count as f64;
                bulk_sq_sum += v * v * use_count as f64;
                bulk_count += use_count;
                cumulative += count;
            }

            if bulk_count > 0 {
                let mean = bulk_sum / bulk_count as f64;
                let sd = if bulk_count > 1 {
                    let variance = (bulk_sq_sum / bulk_count as f64) - (mean * mean);
                    if variance > 0.0 {
                        variance.sqrt()
                    } else {
                        0.0
                    }
                } else {
                    0.0
                };
                (mean, sd)
            } else {
                (0.0, 0.0)
            }
        } else {
            (0.0, 0.0)
        }
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
    sn(
        &mut out,
        "error rate:",
        fmt_sci(error_rate),
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
    // Orientation counts: each pair is counted once (only the upstream mate
    // contributes), matching samtools stats behavior. No division needed.
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

    // =================================================================
    // Histogram/distribution sections
    // =================================================================

    // FFQ — first fragment quality per cycle
    write_ffq_lfq(&mut out, "FFQ", "first fragment", &result.ffq)?;

    // LFQ — last fragment quality per cycle
    write_ffq_lfq(&mut out, "LFQ", "last fragment", &result.lfq)?;

    // GCF — GC content of first fragments
    write_gc_content(&mut out, "GCF", "first fragments", &result.gcf)?;

    // GCL — GC content of last fragments
    write_gc_content(&mut out, "GCL", "last fragments", &result.gcl)?;

    // GCC — ACGT content per cycle (derived from FBC+LBC as percentages)
    write_gcc(
        &mut out,
        "GCC",
        "ACGT content per cycle",
        &result.fbc,
        &result.lbc,
    )?;

    // GCT — ACGT content per cycle, read-oriented
    write_gcc(
        &mut out,
        "GCT",
        "ACGT content per cycle, read-oriented (reversed for reverse reads)",
        &result.fbc_ro,
        &result.lbc_ro,
    )?;

    // FBC — ACGT raw counts per cycle, first fragments
    write_base_counts(&mut out, "FBC", "first fragments", &result.fbc)?;

    // LBC — ACGT raw counts per cycle, last fragments
    write_base_counts(&mut out, "LBC", "last fragments", &result.lbc)?;

    // FTC — total base counters, first fragments
    write_total_base_counts(&mut out, "FTC", &result.ftc)?;

    // LTC — total base counters, last fragments
    write_total_base_counts(&mut out, "LTC", &result.ltc)?;

    // IS — insert size with orientation breakdown
    write_insert_size(&mut out, &result.is_hist)?;

    // RL — read length histogram
    write_length_hist(&mut out, "RL", &result.rl_hist)?;

    // FRL — first fragment read lengths
    write_length_hist(&mut out, "FRL", &result.frl_hist)?;

    // LRL — last fragment read lengths
    write_length_hist(&mut out, "LRL", &result.lrl_hist)?;

    // MAPQ — mapping quality histogram
    write_mapq_hist(&mut out, &result.mapq_hist)?;

    // ID — indel size distribution
    write_indel_dist(&mut out, &result.id_hist)?;

    // IC — indels per cycle
    write_indel_cycle(&mut out, &result.ic)?;

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
// Histogram/distribution output helpers
// ============================================================================

/// Write FFQ or LFQ section (per-cycle quality distribution).
///
/// Each row is a cycle, columns are quality values 0..max_q.
/// `tag` is "FFQ" or "LFQ", `desc` is "first fragment" or "last fragment".
fn write_ffq_lfq<W: std::io::Write>(
    out: &mut W,
    tag: &str,
    desc: &str,
    data: &[[u64; 64]],
) -> Result<()> {
    // Determine max quality observed across all cycles
    let max_q = data
        .iter()
        .flat_map(|row| {
            row.iter()
                .enumerate()
                .rev()
                .find(|(_, &c)| c > 0)
                .map(|(i, _)| i)
        })
        .max()
        .unwrap_or(0);

    writeln!(
        out,
        "# {tag}, {desc} quality. Columns correspond to qualities and rows to cycles. First column is the cycle number."
    )?;
    writeln!(
        out,
        "# Columns correspond to qualities 0, 1, 2, ..., {max_q}"
    )?;
    writeln!(
        out,
        "# Use `plot-bamstats -p <outdir> <file>` to generate plots"
    )?;
    for (cycle, row) in data.iter().enumerate() {
        write!(out, "{tag}\t{cycle}")?;
        for &count in row.iter().take(max_q + 1) {
            write!(out, "\t{count}")?;
        }
        writeln!(out)?;
    }
    Ok(())
}

/// Write GCF or GCL section (GC content distribution, 101 buckets 0-100%).
fn write_gc_content<W: std::io::Write>(
    out: &mut W,
    tag: &str,
    desc: &str,
    data: &[u64; 101],
) -> Result<()> {
    writeln!(
        out,
        "# {tag}, GC content of {desc}. Use `plot-bamstats -p <outdir> <file>` to generate plots"
    )?;
    for (i, &count) in data.iter().enumerate() {
        writeln!(out, "{tag}\t{:.2}\t{count}", i as f64)?;
    }
    Ok(())
}

/// Write GCC or GCT section (ACGT content per cycle as percentages).
///
/// Combines first-fragment and last-fragment base counts to produce
/// per-cycle percentages of A, C, G, T, N, Other.
fn write_gcc<W: std::io::Write>(
    out: &mut W,
    tag: &str,
    desc: &str,
    fbc: &[[u64; 6]],
    lbc: &[[u64; 6]],
) -> Result<()> {
    writeln!(
        out,
        "# {tag}, {desc}. Use `plot-bamstats -p <outdir> <file>` to generate plots"
    )?;

    let max_cycles = fbc.len().max(lbc.len());
    for cycle in 0..max_cycles {
        let mut combined = [0u64; 6];
        if cycle < fbc.len() {
            for j in 0..6 {
                combined[j] += fbc[cycle][j];
            }
        }
        if cycle < lbc.len() {
            for j in 0..6 {
                combined[j] += lbc[cycle][j];
            }
        }
        let total: u64 = combined.iter().sum();
        if total > 0 {
            let pcts: Vec<f64> = combined
                .iter()
                .map(|&c| 100.0 * c as f64 / total as f64)
                .collect();
            writeln!(
                out,
                "{tag}\t{cycle}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}",
                pcts[0], pcts[1], pcts[2], pcts[3], pcts[4], pcts[5]
            )?;
        } else {
            writeln!(out, "{tag}\t{cycle}\t0.00\t0.00\t0.00\t0.00\t0.00\t0.00")?;
        }
    }
    Ok(())
}

/// Write FBC or LBC section (ACGT raw counts per cycle).
fn write_base_counts<W: std::io::Write>(
    out: &mut W,
    tag: &str,
    desc: &str,
    data: &[[u64; 6]],
) -> Result<()> {
    writeln!(
        out,
        "# {tag}, ACGT content per cycle for {desc}. Use `plot-bamstats -p <outdir> <file>` to generate plots"
    )?;
    for (cycle, row) in data.iter().enumerate() {
        writeln!(
            out,
            "{tag}\t{cycle}\t{}\t{}\t{}\t{}\t{}\t{}",
            row[0], row[1], row[2], row[3], row[4], row[5]
        )?;
    }
    Ok(())
}

/// Write FTC or LTC section (total base counters, one line).
fn write_total_base_counts<W: std::io::Write>(
    out: &mut W,
    tag: &str,
    data: &[u64; 5],
) -> Result<()> {
    writeln!(
        out,
        "# {tag}, total base counters. Use `plot-bamstats -p <outdir> <file>` to generate plots"
    )?;
    writeln!(
        out,
        "{tag}\t{}\t{}\t{}\t{}\t{}",
        data[0], data[1], data[2], data[3], data[4]
    )?;
    Ok(())
}

/// Write IS section (insert size with orientation breakdown).
fn write_insert_size<W: std::io::Write>(
    out: &mut W,
    is_hist: &std::collections::HashMap<u64, [u64; 4]>,
) -> Result<()> {
    writeln!(
        out,
        "# IS, insert size. Use `plot-bamstats -p <outdir> <file>` to generate plots"
    )?;

    if is_hist.is_empty() {
        return Ok(());
    }

    let max_isize = *is_hist.keys().max().unwrap_or(&0);
    // Cap at 8000 to match samtools stats default MAX_INSERT_SIZE
    let limit = max_isize.min(8000);
    for isize_val in 0..=limit {
        let counts = is_hist.get(&isize_val).copied().unwrap_or([0; 4]);
        writeln!(
            out,
            "IS\t{isize_val}\t{}\t{}\t{}\t{}",
            counts[0], counts[1], counts[2], counts[3]
        )?;
    }
    Ok(())
}

/// Write RL, FRL, or LRL section (read length histogram).
fn write_length_hist<W: std::io::Write>(
    out: &mut W,
    tag: &str,
    hist: &std::collections::HashMap<u64, u64>,
) -> Result<()> {
    writeln!(
        out,
        "# {tag}, read length. Use `plot-bamstats -p <outdir> <file>` to generate plots"
    )?;

    if hist.is_empty() {
        return Ok(());
    }

    let min_len = *hist.keys().min().unwrap_or(&0);
    let max_len = *hist.keys().max().unwrap_or(&0);
    for len in min_len..=max_len {
        let count = hist.get(&len).copied().unwrap_or(0);
        writeln!(out, "{tag}\t{len}\t{count}")?;
    }
    Ok(())
}

/// Write MAPQ section (mapping quality histogram).
fn write_mapq_hist<W: std::io::Write>(out: &mut W, mapq_hist: &[u64; 256]) -> Result<()> {
    writeln!(
        out,
        "# MAPQ, mapping quality. Use `plot-bamstats -p <outdir> <file>` to generate plots"
    )?;

    // Find max observed MAPQ
    let max_mapq = mapq_hist
        .iter()
        .enumerate()
        .rev()
        .find(|(_, &c)| c > 0)
        .map(|(i, _)| i)
        .unwrap_or(0);

    for (q, &count) in mapq_hist.iter().enumerate().take(max_mapq + 1) {
        writeln!(out, "MAPQ\t{q}\t{count}")?;
    }
    Ok(())
}

/// Write ID section (indel size distribution).
fn write_indel_dist<W: std::io::Write>(
    out: &mut W,
    id_hist: &std::collections::HashMap<u64, [u64; 2]>,
) -> Result<()> {
    writeln!(
        out,
        "# ID, indel size. Use `plot-bamstats -p <outdir> <file>` to generate plots"
    )?;

    if id_hist.is_empty() {
        return Ok(());
    }

    let max_len = *id_hist.keys().max().unwrap_or(&0);
    for len in 1..=max_len {
        let counts = id_hist.get(&len).copied().unwrap_or([0; 2]);
        writeln!(out, "ID\t{len}\t{}\t{}", counts[0], counts[1])?;
    }
    Ok(())
}

/// Write IC section (indels per cycle).
fn write_indel_cycle<W: std::io::Write>(out: &mut W, ic: &[[u64; 4]]) -> Result<()> {
    writeln!(
        out,
        "# IC, indels per cycle. Use `plot-bamstats -p <outdir> <file>` to generate plots"
    )?;
    for (cycle, row) in ic.iter().enumerate() {
        writeln!(
            out,
            "IC\t{cycle}\t{}\t{}\t{}\t{}",
            row[0], row[1], row[2], row[3]
        )?;
    }
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
