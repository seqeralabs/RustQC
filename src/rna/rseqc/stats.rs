//! Samtools stats compatible output.
//!
//! Produces the Summary Numbers (SN) section and all histogram/distribution
//! sections from `samtools stats`, compatible with MultiQC and plot-bamstats.

use std::collections::HashMap;
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

    // Insert size stats — 99th percentile truncation matching samtools stats.c.
    // Each pair is counted once (upstream mate only), capped at 8000.
    //
    // Upstream samtools stats.c logic (report_stats, lines 1492-1522):
    //   1. nisize = total count across ALL bins (including isize=0)
    //   2. Iterate bins: accumulate bulk count and weighted sum.
    //      Stop (including current bin) when bulk/nisize > 0.99.
    //      Then nisize = bulk (truncated count).
    //   3. mean = weighted_sum / nisize
    //   4. SD loop runs from isize=1 to ibulk-1 (skips isize=0):
    //      sd = sqrt( sum(count[i] * (i - mean)^2) / nisize )  [population SD]
    let (avg_insert, insert_sd) = {
        let total_count: u64 = result.is_hist.values().map(|v| v[0]).sum();
        if total_count > 0 {
            let mut sorted_isizes: Vec<(u64, u64)> =
                result.is_hist.iter().map(|(&k, v)| (k, v[0])).collect();
            sorted_isizes.sort_by_key(|&(k, _)| k);

            // Phase 1: Find 99th percentile cutoff and compute mean.
            // Include entire bins; stop when cumulative/total > 0.99.
            let mut bulk_count = 0u64;
            let mut weighted_sum = 0.0f64;
            let mut ibulk: u64 = 0; // exclusive upper bound for SD loop

            for &(isize_val, count) in &sorted_isizes {
                if count > 0 {
                    ibulk = isize_val + 1;
                }
                bulk_count += count;
                weighted_sum += isize_val as f64 * count as f64;

                if bulk_count as f64 / total_count as f64 > 0.99 {
                    break;
                }
            }

            if bulk_count > 0 {
                let mean = weighted_sum / bulk_count as f64;

                // Phase 2: SD — population SD, loop from isize=1 (skip 0).
                let mut sd_sum = 0.0f64;
                for &(isize_val, count) in &sorted_isizes {
                    if isize_val == 0 {
                        continue; // upstream SD loop starts at isize=1
                    }
                    if isize_val >= ibulk {
                        break;
                    }
                    let diff = isize_val as f64 - mean;
                    sd_sum += count as f64 * diff * diff;
                }
                let sd = (sd_sum / bulk_count as f64).sqrt();
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

    // COV — coverage distribution
    write_coverage_dist(&mut out, &result.cov_hist)?;

    // CHK — CRC32 checksums
    writeln!(
        out,
        "# CHK, Pair of checksum values for read names, sequences and qualities, calculated using a CRC32 algorithm."
    )?;
    writeln!(
        out,
        "CHK\t{:08x}\t{:08x}\t{:08x}",
        result.chk[0], result.chk[1], result.chk[2]
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
        write!(out, "{tag}\t{}", cycle + 1)?; // 1-based cycle numbering
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
    data: &[u64; 200],
) -> Result<()> {
    writeln!(
        out,
        "# {tag}, GC content of {desc}. Use `plot-bamstats -p <outdir> <file>` to generate plots"
    )?;
    // Output the cumulative step function with run-length encoding.
    // Matches samtools stats.c:1633-1640: skip bins with unchanged count,
    // print midpoint GC% between previous and current bin index.
    let ngc: usize = 200;
    let mut ibase_prev: usize = 0;
    for ibase in 0..ngc {
        if data[ibase] == data[ibase_prev] {
            continue;
        }
        let gc_pct = (ibase + ibase_prev) as f64 * 0.5 * 100.0 / (ngc - 1) as f64;
        writeln!(out, "{tag}\t{gc_pct:.2}\t{}", data[ibase_prev])?;
        ibase_prev = ibase;
    }
    // Print the last plateau if it has any counts
    if data[ibase_prev] > 0 {
        let gc_pct = ((ngc - 1) + ibase_prev) as f64 * 0.5 * 100.0 / (ngc - 1) as f64;
        writeln!(out, "{tag}\t{gc_pct:.2}\t{}", data[ibase_prev])?;
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
        // 1-based cycle numbering; only output A%, C%, G%, T% (4 columns, matching upstream)
        if total > 0 {
            let pcts: Vec<f64> = combined
                .iter()
                .map(|&c| 100.0 * c as f64 / total as f64)
                .collect();
            writeln!(
                out,
                "{tag}\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.2}",
                cycle + 1,
                pcts[0],
                pcts[1],
                pcts[2],
                pcts[3]
            )?;
        } else {
            writeln!(out, "{tag}\t{}\t0.00\t0.00\t0.00\t0.00", cycle + 1)?;
        }
    }
    Ok(())
}

/// Write FBC or LBC section (ACGT content per cycle as percentages).
///
/// Upstream samtools stats outputs percentages (not raw counts) with
/// 6 columns: A%, C%, G%, T%, N%, Other%.
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
        let total: u64 = row.iter().sum();
        if total > 0 {
            let pcts: Vec<f64> = row
                .iter()
                .map(|&c| 100.0 * c as f64 / total as f64)
                .collect();
            writeln!(
                out,
                "{tag}\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}",
                cycle + 1, // 1-based
                pcts[0],
                pcts[1],
                pcts[2],
                pcts[3],
                pcts[4],
                pcts[5]
            )?;
        } else {
            writeln!(
                out,
                "{tag}\t{}\t0.00\t0.00\t0.00\t0.00\t0.00\t0.00",
                cycle + 1
            )?;
        }
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

    let mut lengths: Vec<u64> = hist.keys().copied().collect();
    lengths.sort();
    for len in lengths {
        let count = hist[&len];
        if count > 0 {
            writeln!(out, "{tag}\t{len}\t{count}")?;
        }
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

    let mut sizes: Vec<u64> = id_hist.keys().copied().collect();
    sizes.sort();
    for len in sizes {
        let counts = id_hist[&len];
        if counts[0] > 0 || counts[1] > 0 {
            writeln!(out, "ID\t{len}\t{}\t{}", counts[0], counts[1])?;
        }
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
        // Only output cycles with non-zero indel counts (matching upstream)
        if row[0] > 0 || row[1] > 0 || row[2] > 0 || row[3] > 0 {
            writeln!(
                out,
                "IC\t{}\t{}\t{}\t{}\t{}",
                cycle + 1, // 1-based cycle numbering
                row[0],
                row[1],
                row[2],
                row[3]
            )?;
        }
    }
    Ok(())
}

/// Write COV section (coverage distribution).
///
/// Uses default samtools stats parameters: cov_min=1, cov_max=1000, cov_step=1.
/// Format: `COV\t[range]\t<depth>\t<count>`
fn write_coverage_dist<W: std::io::Write>(out: &mut W, cov_hist: &HashMap<u32, u64>) -> Result<()> {
    writeln!(
        out,
        "# COV, Coverage distribution. Use `plot-bamstats -p <outdir> <file>` to generate plots"
    )?;

    let cov_min: u32 = 1;
    let cov_max: u32 = 1000;
    let cov_step: u32 = 1;

    // Count positions below cov_min
    let below_min: u64 = cov_hist
        .iter()
        .filter(|(&d, _)| d < cov_min)
        .map(|(_, &c)| c)
        .sum();

    if below_min > 0 {
        writeln!(out, "COV\t[<{cov_min}]\t{}\t{below_min}", cov_min - 1)?;
    }

    // Regular bins
    let mut depth = cov_min;
    while depth <= cov_max {
        let end = depth + cov_step - 1;
        // Sum counts in this bin range
        let count: u64 = (depth..=end).filter_map(|d| cov_hist.get(&d)).sum();
        if count > 0 {
            if cov_step == 1 {
                writeln!(out, "COV\t[{depth}-{depth}]\t{depth}\t{count}")?;
            } else {
                writeln!(out, "COV\t[{depth}-{end}]\t{depth}\t{count}")?;
            }
        }
        depth += cov_step;
    }

    // Overflow bin (> cov_max)
    let above_max: u64 = cov_hist
        .iter()
        .filter(|(&d, _)| d > cov_max)
        .map(|(_, &c)| c)
        .sum();
    if above_max > 0 {
        writeln!(out, "COV\t[{cov_max}<]\t{cov_max}\t{above_max}")?;
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
