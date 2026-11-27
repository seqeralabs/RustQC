//! Samtools stats compatible output.
//!
//! Produces the Summary Numbers (SN) section and all histogram/distribution
//! sections from `samtools stats`, compatible with MultiQC and plot-bamstats.

use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use log::debug;

use super::bam_stat::{BamStatResult, GcDepthBin};

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

/// Insert size summary computed like samtools `output_stats`.
struct InsertSizeStats {
    /// Per insert size `[inward, outward, other]`, each halved (both mates
    /// were counted) and truncated, for insert sizes `0..ibulk`.
    rows: Vec<[u64; 3]>,
    /// Totals of the halved counts over all insert sizes.
    inward: u64,
    outward: u64,
    other: u64,
    avg: f64,
    sd: f64,
}

/// samtools `isize_main_bulk` default, a `float` in samtools.
const ISIZE_MAIN_BULK: f32 = 0.99;
/// samtools insert size cap (`-i`, default 8000).
const MAX_INSERT_SIZE: usize = 8000;

/// Port of the insert size part of samtools `output_stats`: halve each
/// orientation count, keep the main bulk (99%) of pairs, then compute the
/// mean and the population standard deviation over insert sizes `1..ibulk`.
fn insert_size_stats(is_hist: &HashMap<u64, [u64; 4]>) -> InsertSizeStats {
    let mut halved = vec![[0u64; 3]; MAX_INSERT_SIZE + 1];
    for (&isize, counts) in is_hist {
        if let Some(row) = halved.get_mut(isize as usize) {
            for (half, &count) in row.iter_mut().zip(&counts[1..]) {
                *half = (count as f64 * 0.5) as u64;
            }
        }
    }
    let (mut inward, mut outward, mut other) = (0u64, 0u64, 0u64);
    for row in &halved {
        inward += row[0];
        outward += row[1];
        other += row[2];
    }
    let mut nisize = inward + outward + other;

    let main_bulk = f64::from(ISIZE_MAIN_BULK);
    let mut ibulk = 0usize;
    let mut bulk = 0.0f64;
    let mut avg = 0.0f64;
    for (isize, row) in halved.iter().enumerate() {
        let num = row[0] + row[1] + row[2];
        if num > 0 {
            ibulk = isize + 1;
        }
        bulk += num as f64;
        avg += (isize as u64 * num) as f64;
        if bulk / nisize as f64 > main_bulk {
            ibulk = isize + 1;
            nisize = bulk as u64;
            break;
        }
    }
    let denom = if nisize > 0 { nisize } else { 1 } as f64;
    avg /= denom;
    let mut sd = 0.0f64;
    for (isize, row) in halved.iter().enumerate().take(ibulk).skip(1) {
        let num = (row[0] + row[1] + row[2]) as f64;
        let diff = isize as f64 - avg;
        sd += num * diff * diff / denom;
    }
    halved.truncate(ibulk);
    InsertSizeStats {
        rows: halved,
        inward,
        outward,
        other,
        avg,
        sd: sd.sqrt(),
    }
}

/// `a / b` in single precision, as samtools computes its `float` averages.
fn ratio_f32(a: u64, b: u64) -> f64 {
    if b == 0 {
        0.0
    } else {
        f64::from(a as f32 / b as f32)
    }
}

/// Write samtools stats compatible output.
///
/// The data lines match `samtools stats` (1.22.1) run without options on a
/// coordinate-sorted BAM. The output starts with a samtools-style comment
/// header (MultiQC detects "This file was produced by samtools stats").
///
/// # Arguments
/// * `result` - The computed BAM statistics
/// * `output_path` - Path to write the stats file
pub fn write_stats(result: &BamStatResult, output_path: &Path) -> Result<()> {
    let file = std::fs::File::create(output_path)
        .with_context(|| format!("Failed to create stats file: {}", output_path.display()))?;
    let mut out = std::io::BufWriter::new(file);

    // Header comment — MultiQC detects "This file was produced by samtools stats"
    writeln!(out, "# This file was produced by samtools stats and RustQC")?;
    writeln!(
        out,
        "# The command line was: rustqc rna (samtools stats compatible output)"
    )?;

    // Reads counted by samtools stats: primary reads with a sequence.
    let sequences = result.first_fragments + result.last_fragments + result.other_fragments;
    let is = insert_size_stats(&result.is_hist);

    // samtools: (float)nmismatches / nbases_mapped_cigar
    let error_rate = ratio_f32(result.mismatches, result.bases_mapped_cigar);
    // samtools: float avg_read_length = (float)total_len / sequences
    let avg_read_len = ratio_f32(result.total_len, sequences);
    let avg_quality = if result.total_len > 0 {
        result.quality_sum / result.total_len as f64
    } else {
        0.0
    };

    // CHK — CRC32 checksums (written before SN, matching upstream order)
    writeln!(
        out,
        "# CHK, Pair of checksum values for read names, sequences and qualities, calculated using a CRC32 algorithm."
    )?;
    writeln!(
        out,
        "CHK\t{:08x}\t{:08x}\t{:08x}",
        result.chk[0], result.chk[1], result.chk[2]
    )?;

    // Write SN lines
    sn(
        &mut out,
        "raw total sequences:",
        sequences,
        "excluding supplementary and secondary reads",
    )?;
    sn_no_comment(&mut out, "filtered sequences:", 0_u64)?;
    sn_no_comment(&mut out, "sequences:", sequences)?;
    sn(&mut out, "is sorted:", 1_u64, "sorted by coordinate")?;
    sn_no_comment(&mut out, "1st fragments:", result.first_fragments)?;
    sn_no_comment(&mut out, "last fragments:", result.last_fragments)?;
    sn_no_comment(&mut out, "reads mapped:", sequences - result.stats_unmapped)?;
    sn(
        &mut out,
        "reads mapped and paired:",
        result.reads_mapped_and_paired,
        "paired-end technology bit set + both mates mapped",
    )?;
    sn_no_comment(&mut out, "reads unmapped:", result.stats_unmapped)?;
    sn(
        &mut out,
        "reads properly paired:",
        result.stats_properly_paired,
        "proper-pair bit set",
    )?;
    sn(
        &mut out,
        "reads paired:",
        result.stats_paired_tech,
        "paired-end technology bit set",
    )?;
    sn(
        &mut out,
        "reads duplicated:",
        result.stats_dup_reads,
        "PCR or optical duplicate bit set",
    )?;
    sn(&mut out, "reads MQ0:", result.reads_mq0, "mapped and MQ=0")?;
    sn_no_comment(&mut out, "reads QC failed:", result.stats_qcfail)?;
    sn_no_comment(&mut out, "non-primary alignments:", result.secondary)?;
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
    sn_no_comment(&mut out, "bases duplicated:", result.stats_dup_bases)?;
    // samtools prints the unsigned counter through a (long) cast.
    sn(
        &mut out,
        "mismatches:",
        result.mismatches as i64,
        "from NM fields",
    )?;
    sn(
        &mut out,
        "error rate:",
        fmt_sci(error_rate),
        "mismatches / bases mapped (cigar)",
    )?;
    sn_no_comment(&mut out, "average length:", format!("{avg_read_len:.0}"))?;
    sn_no_comment(
        &mut out,
        "average first fragment length:",
        format!(
            "{:.0}",
            ratio_f32(result.total_first_fragment_len, result.first_fragments)
        ),
    )?;
    sn_no_comment(
        &mut out,
        "average last fragment length:",
        format!(
            "{:.0}",
            ratio_f32(result.total_last_fragment_len, result.last_fragments)
        ),
    )?;
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
    sn_no_comment(&mut out, "insert size average:", fmt_1dp(is.avg))?;
    sn_no_comment(&mut out, "insert size standard deviation:", fmt_1dp(is.sd))?;
    sn_no_comment(&mut out, "inward oriented pairs:", is.inward)?;
    sn_no_comment(&mut out, "outward oriented pairs:", is.outward)?;
    sn_no_comment(&mut out, "pairs with other orientation:", is.other)?;
    sn_no_comment(
        &mut out,
        "pairs on different chromosomes:",
        result.stats_anomalous / 2,
    )?;
    // samtools: (float)(100 * nreads_properly_paired) / sequences
    sn_no_comment(
        &mut out,
        "percentage of properly paired reads (%):",
        fmt_1dp(ratio_f32(100 * result.stats_properly_paired, sequences)),
    )?;

    // =================================================================
    // Histogram/distribution sections
    // =================================================================

    // samtools uses one maximum quality for both FFQ and LFQ.
    let max_qual = max_quality(&result.ffq).max(max_quality(&result.lfq));

    // FFQ — first fragment quality per cycle
    write_ffq_lfq(
        &mut out,
        "FFQ",
        "first fragment",
        &result.ffq,
        result.max_first_fragment_len as usize,
        max_qual,
    )?;

    // LFQ — last fragment quality per cycle
    write_ffq_lfq(
        &mut out,
        "LFQ",
        "last fragment",
        &result.lfq,
        result.max_last_fragment_len as usize,
        max_qual,
    )?;

    // GCF — GC content of first fragments
    write_gc_content(&mut out, "GCF", "first fragments", &result.gcf)?;

    // GCL — GC content of last fragments
    write_gc_content(&mut out, "GCL", "last fragments", &result.gcl)?;

    // GCC — ACGT content per cycle, read-oriented (derived from FBC_ro+LBC_ro as percentages)
    write_gcc(
        &mut out,
        "GCC",
        "ACGT content per cycle",
        &result.fbc_ro,
        &result.lbc_ro,
    )?;

    // GCT — ACGT content per cycle, read-oriented with reverse-complemented bases
    write_gct(
        &mut out,
        "GCT",
        "ACGT content per cycle, read-oriented (display complement of bases of reverse reads)",
        &result.gcc_rc,
    )?;

    // FBC — ACGT base content per cycle, first fragments (read-oriented)
    write_base_counts(&mut out, "FBC", "first fragments", &result.fbc_ro)?;

    // FTC — total base counters, first fragments
    write_total_base_counts(&mut out, "FTC", &result.ftc)?;

    // LBC — ACGT base content per cycle, last fragments (read-oriented)
    write_base_counts(&mut out, "LBC", "last fragments", &result.lbc_ro)?;

    // LTC — total base counters, last fragments
    write_total_base_counts(&mut out, "LTC", &result.ltc)?;

    // IS — insert size with orientation breakdown
    write_insert_size(&mut out, &is)?;

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

    // GCD — GC-depth
    write_gc_depth(&mut out, &result.gcd_bins, avg_read_len as f32)?;

    out.flush()?;
    debug!("Wrote samtools stats output to {}", output_path.display());
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

/// Highest quality with a non-zero count in a per-cycle quality table.
fn max_quality(data: &[[u64; 256]]) -> usize {
    data.iter()
        .flat_map(|row| row.iter().rposition(|&c| c > 0))
        .max()
        .unwrap_or(0)
}

/// Write FFQ or LFQ section (per-cycle quality distribution).
///
/// Each row is a cycle, columns are qualities `0..=max_qual + 1` (samtools
/// prints one column past the highest quality seen in either fragment).
/// samtools prints a row for every cycle up to the longest read of that
/// fragment, hard clips included; cycles past the stored ones are zero.
/// `tag` is "FFQ" or "LFQ", `desc` is "first fragment" or "last fragment".
fn write_ffq_lfq<W: std::io::Write>(
    out: &mut W,
    tag: &str,
    desc: &str,
    data: &[[u64; 256]],
    n_cycles: usize,
    max_qual: usize,
) -> Result<()> {
    // stats.c: `if (max_qual+1 < nquals) max_qual++` with nquals = 256
    let max_q = if max_qual + 1 < 256 {
        max_qual + 1
    } else {
        max_qual
    };

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
    const EMPTY: [u64; 256] = [0; 256];
    for cycle in 0..n_cycles {
        let row = data.get(cycle).unwrap_or(&EMPTY);
        write!(out, "{tag}\t{}", cycle + 1)?; // 1-based cycle numbering
        for &count in row.iter().take(max_q + 1) {
            write!(out, "\t{count}")?;
        }
        writeln!(out)?;
    }
    Ok(())
}

/// Write GCF or GCL section (GC content distribution, 200-bin cumulative step function).
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

/// Write GCC section (ACGT content per cycle as percentages).
///
/// Combines first-fragment and last-fragment base counts to produce
/// per-cycle percentages of A, C, G, T, N, Other.
///
/// Upstream samtools stats uses only A+C+G+T as the denominator (not
/// including N or Other). N and Other are then expressed as a percentage
/// of that ACGT-only total, so the six columns can sum to more than 100%.
/// Cycles where `acgt_sum == 0` are skipped entirely (matching upstream).
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
            for (count, &value) in combined.iter_mut().zip(&fbc[cycle]) {
                *count += value;
            }
        }
        if cycle < lbc.len() {
            for (count, &value) in combined.iter_mut().zip(&lbc[cycle]) {
                *count += value;
            }
        }
        // Denominator is A+C+G+T only (first 4 elements), matching upstream
        // samtools stats.c which computes acgt_sum without N or Other.
        let acgt_sum: u64 = combined[..4].iter().sum();
        if acgt_sum > 0 {
            let pcts: Vec<f64> = combined
                .iter()
                .map(|&c| 100.0 * c as f64 / acgt_sum as f64)
                .collect();
            writeln!(
                out,
                "{tag}\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}",
                cycle + 1,
                pcts[0],
                pcts[1],
                pcts[2],
                pcts[3],
                pcts[4],
                pcts[5]
            )?;
        }
        // Cycles with acgt_sum == 0 are skipped (matching upstream `if (!acgt_sum) continue`)
    }
    Ok(())
}

/// Write GCT section — ACGT content per cycle with reverse-complemented bases
/// for reverse-strand reads. Combined first+last fragments. Only 4 value columns
/// (A%, C%, G%, T%), matching upstream samtools stats format.
fn write_gct<W: std::io::Write>(
    out: &mut W,
    tag: &str,
    desc: &str,
    gcc_rc: &[[u64; 4]],
) -> std::io::Result<()> {
    writeln!(
        out,
        "# {desc}. Use `grep ^{tag} | cut -f 2-` to extract this part."
    )?;
    for (cycle, bases) in gcc_rc.iter().enumerate() {
        let total: u64 = bases.iter().sum();
        if total > 0 {
            let pcts: Vec<f64> = bases
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
        }
        // Cycles without A/C/G/T are skipped (samtools `if (!acgt_sum) continue`).
    }
    Ok(())
}

/// Write FBC or LBC section (ACGT content per cycle as percentages).
///
/// Upstream samtools stats outputs percentages (not raw counts) with
/// 6 columns: A%, C%, G%, T%, N%, Other%.
///
/// The denominator is A+C+G+T only (not including N or Other), matching
/// upstream samtools stats.c. N and Other are expressed as a percentage
/// of the ACGT-only total. Cycles where `acgt_sum == 0` are skipped.
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
        // Denominator is A+C+G+T only (first 4 elements), matching upstream
        // samtools stats.c which computes acgt_sum without N or Other.
        let acgt_sum: u64 = row[..4].iter().sum();
        if acgt_sum > 0 {
            let pcts: Vec<f64> = row
                .iter()
                .map(|&c| 100.0 * c as f64 / acgt_sum as f64)
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
        }
        // Cycles with acgt_sum == 0 are skipped (matching upstream)
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

/// Write IS section (insert size with orientation breakdown), for the insert
/// sizes in the main bulk (see [`insert_size_stats`]).
fn write_insert_size<W: std::io::Write>(out: &mut W, is: &InsertSizeStats) -> Result<()> {
    writeln!(
        out,
        "# IS, insert size. Use `plot-bamstats -p <outdir> <file>` to generate plots"
    )?;
    for (isize, &[inward, outward, other]) in is.rows.iter().enumerate() {
        writeln!(
            out,
            "IS\t{isize}\t{}\t{inward}\t{outward}\t{other}",
            inward + outward + other
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

    // Sparse format: only output non-zero entries (matching upstream samtools)
    for (q, &count) in mapq_hist.iter().enumerate().take(max_mapq + 1) {
        if count > 0 {
            writeln!(out, "MAPQ\t{q}\t{count}")?;
        }
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
// GCD — GC-depth distribution
// ============================================================================

/// samtools `gcd_percentile`: the `p`-th percentile of the sorted `depths`
/// by linear interpolation, in single precision.
///
/// The depth difference is an unsigned 32-bit subtraction in samtools; it is
/// kept as such (it can only wrap if a group mixes the raw last bin with a
/// normalised one).
fn gcd_percentile(depths: &[u32], p: u32) -> f32 {
    let n_bins = depths.len() as i32;
    let n = p as f32 * (n_bins + 1) as f32 / 100.0;
    let k = n as i32;
    if k <= 0 {
        return depths[0] as f32;
    }
    if k >= n_bins {
        return depths[(n_bins - 1) as usize] as f32;
    }
    let d = n - k as f32;
    let (lo, hi) = (depths[(k - 1) as usize], depths[k as usize]);
    lo as f32 + d * hi.wrapping_sub(lo) as f32
}

/// samtools `rint(100. * gc / depth)`: mean GC of a bin as a whole percentage,
/// computed in double precision and rounded half to even.
fn gcd_bin_gc_percent(gc: f32, depth: u32) -> f32 {
    (100.0 * f64::from(gc) / f64::from(depth)).round_ties_even() as f32
}

/// Write the GCD (GC-depth) section.
///
/// Port of samtools `output_stats` without `--ref-seq` (gc_depth_t holds a
/// `float` GC and a `uint32_t` depth; all depth arithmetic is `float`):
/// - The bins are `gcd[1..=igcd]`, with an unused zero bin `gcd[0]` in front.
/// - Every bin except the last (`gcd[igcd]`) gets its GC normalised to a
///   whole percentage; the last keeps its raw sum of per-read GC fractions.
/// - All `igcd + 1` entries are sorted by (GC, depth) and the first `igcd` are
///   printed, grouped by GC (within 0.1), with depth percentiles scaled by
///   `avg_read_len / GCD_BIN_SIZE`.
///
/// `gcd_bins` must be in BAM order (contig, then position).
fn write_gc_depth(
    out: &mut impl Write,
    gcd_bins: &[GcDepthBin],
    avg_read_len: f32,
) -> std::io::Result<()> {
    writeln!(
        out,
        "# GC-depth. Use `grep ^GCD | cut -f 2-` to extract this part. \
         The columns are: GC%, unique sequence percentiles, 10th, 25th, \
         50th, 75th and 90th depth percentile"
    )?;
    if gcd_bins.is_empty() {
        return Ok(());
    }

    let igcd = gcd_bins.len();
    let mut bins: Vec<(f32, u32)> = Vec::with_capacity(igcd + 1);
    bins.push((0.0, 0)); // gcd[0]
    for (idx, b) in gcd_bins.iter().enumerate() {
        let gc = if idx + 1 < igcd && b.depth > 0 {
            gcd_bin_gc_percent(b.gc, b.depth)
        } else {
            b.gc
        };
        bins.push((gc, b.depth));
    }
    // gcd_cmp: by GC, then depth. Equal entries are identical, so the
    // (unstable) qsort order cannot show.
    bins.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));

    const GCD_BIN_SIZE: f32 = 20_000.0;
    let scale = |depth: f32| f64::from(depth * avg_read_len / GCD_BIN_SIZE);

    let mut i = 0;
    while i < igcd {
        let gc = bins[i].0;
        let mut end = i;
        // fabs(float - float) < 0.1 (double)
        while end < igcd && f64::from((bins[end].0 - gc).abs()) < 0.1 {
            end += 1;
        }
        let depths: Vec<u32> = bins[i..end].iter().map(|b| b.1).collect();
        writeln!(
            out,
            "GCD\t{:.1}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
            f64::from(gc),
            (end + 1) as f64 * 100.0 / (igcd + 1) as f64,
            scale(gcd_percentile(&depths, 10)),
            scale(gcd_percentile(&depths, 25)),
            scale(gcd_percentile(&depths, 50)),
            scale(gcd_percentile(&depths, 75)),
            scale(gcd_percentile(&depths, 90)),
        )?;
        i = end;
    }

    Ok(())
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rna::bam::{self, Read as BamRead};
    use crate::rna::rseqc::accumulators::BamStatAccum;
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

    // --- Number formatting: Rust must print like C printf ---

    use std::ffi::{c_char, c_int, CStr, CString};

    extern "C" {
        fn snprintf(buf: *mut c_char, size: usize, fmt: *const c_char, ...) -> c_int;
    }

    /// Format one double with the C library's printf.
    fn c_fmt(fmt: &str, v: f64) -> String {
        let fmt = CString::new(fmt).unwrap();
        let mut buf = vec![0 as c_char; 128];
        // SAFETY: `buf` is writable for its full length, `fmt` is a valid
        // C string with exactly one double conversion, and `v` is a double.
        let n = unsafe { snprintf(buf.as_mut_ptr(), buf.len(), fmt.as_ptr(), v) };
        assert!(n >= 0 && (n as usize) < buf.len());
        // SAFETY: snprintf NUL-terminated the buffer.
        unsafe { CStr::from_ptr(buf.as_ptr()) }
            .to_string_lossy()
            .into_owned()
    }

    /// Deterministic pseudo-random sample of f32-derived values.
    fn sample_values() -> Vec<f64> {
        let mut v = vec![
            0.0, 0.5, 1.5, 2.5, 74.5, 98.5, 0.25, 0.35, 0.05, 37.25, 0.125, 0.375, 0.0625, 0.0225,
            1e-10, 1.0, 9.9999995, 99.995, 123456.5,
        ];
        let mut state = 0x2545F4914F6CDD1Du64;
        for _ in 0..20_000 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let a = (state >> 40) as u32;
            let b = 1 + (state & 0xFFFF) as u32;
            v.push(f64::from(a as f32 / b as f32));
            v.push(f64::from((a % 4096) as f32 * 75.0 / 20_000.0));
            v.push(f64::from(a % 1000) / 8.0);
        }
        v
    }

    #[test]
    fn fixed_point_formatting_matches_printf() {
        for v in sample_values() {
            assert_eq!(format!("{v:.0}"), c_fmt("%.0f", v), "%.0f of {v:e}");
            assert_eq!(format!("{v:.1}"), c_fmt("%.1f", v), "%.1f of {v:e}");
            assert_eq!(format!("{v:.2}"), c_fmt("%.2f", v), "%.2f of {v:e}");
            assert_eq!(format!("{v:.3}"), c_fmt("%.3f", v), "%.3f of {v:e}");
        }
    }

    #[test]
    fn scientific_formatting_matches_printf() {
        for v in sample_values() {
            assert_eq!(fmt_sci(v), c_fmt("%e", v), "%e of {v:e}");
        }
    }

    #[test]
    fn error_rate_is_single_precision() {
        // sc_pe_rev: samtools prints 8.805891e-03; the double ratio prints ...892.
        assert_eq!(fmt_sci(ratio_f32(1558, 176_927)), "8.805891e-03");
        assert_eq!(fmt_sci(1558.0 / 176_927.0), "8.805892e-03");
        assert_eq!(fmt_sci(ratio_f32(0, 0)), "0.000000e+00");
    }

    #[test]
    fn gcd_gc_percent_rounds_half_to_even_in_double() {
        // rint(): exact halves go to the even neighbour.
        assert_eq!(gcd_bin_gc_percent(0.5, 4), 12.0); // 12.5 exactly
        assert_eq!(gcd_bin_gc_percent(1.5, 4), 38.0); // 37.5 exactly
        assert_eq!(gcd_bin_gc_percent(0.0, 3), 0.0);
        // Computed in double from the float sum: 100 * 1.01f / 2 is
        // 50.49999952..., while in float 100 * 1.01f rounds to 101.0 and
        // the mean would be exactly 50.5.
        assert_eq!((100.0f32 * 1.01f32 / 2.0).round(), 51.0);
        assert_eq!(gcd_bin_gc_percent(1.01, 2), 50.0);
    }

    #[test]
    fn gcd_percentile_matches_samtools() {
        let d = [1u32, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        // n = p * (N + 1) / 100 = 1.1 → k = 1, d = 0.1
        assert_eq!(gcd_percentile(&d, 10), 1.0 + (1.1f32 - 1.0));
        assert_eq!(gcd_percentile(&d, 50), 5.5);
        assert_eq!(gcd_percentile(&d, 90), 9.0 + (9.9f32 - 9.0));
        // k <= 0 → first; k >= N → last
        assert_eq!(gcd_percentile(&[7, 9], 10), 7.0);
        assert_eq!(gcd_percentile(&[7, 9], 90), 9.0);
        assert_eq!(gcd_percentile(&[3], 50), 3.0);
        // Unsigned difference wraps like uint32_t in C.
        let d = 1.6f32 - 1.0;
        assert_eq!(
            gcd_percentile(&[5, 3, 1], 40),
            5.0 + d * 3u32.wrapping_sub(5) as f32
        );
    }

    #[test]
    fn gcd_last_bin_keeps_raw_gc_and_is_not_printed_alone() {
        // Two bins: the first is normalised to 50%, the last keeps its raw
        // GC sum (0.9). Sorted with the empty gcd[0] in front, only the first
        // igcd = 2 entries are printed, so the 50% bin is left out.
        let bins = vec![
            GcDepthBin { gc: 1.0, depth: 2 },
            GcDepthBin { gc: 0.9, depth: 3 },
        ];
        let mut out = Vec::new();
        write_gc_depth(&mut out, &bins, 100.0).unwrap();
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().filter(|l| l.starts_with("GCD")).collect();
        assert_eq!(
            lines,
            vec![
                "GCD\t0.0\t66.667\t0.000\t0.000\t0.000\t0.000\t0.000",
                "GCD\t0.9\t100.000\t0.015\t0.015\t0.015\t0.015\t0.015",
            ]
        );
    }

    #[test]
    fn insert_size_halves_each_orientation() {
        let mut h = HashMap::new();
        // isize 0 on the same contig is kept; odd counts truncate per orientation.
        h.insert(0u64, [2, 2, 0, 0]);
        h.insert(100u64, [6, 3, 3, 0]);
        h.insert(101u64, [4, 4, 0, 0]);
        let is = insert_size_stats(&h);
        assert_eq!(is.inward, 1 + 1 + 2);
        assert_eq!(is.outward, 1);
        assert_eq!(is.other, 0);
        assert_eq!(is.rows.len(), 102);
        assert_eq!(is.rows[0], [1, 0, 0]);
        assert_eq!(is.rows[100], [1, 1, 0]);
        assert_eq!(is.rows[101], [2, 0, 0]);
        // mean over all 5 halved pairs: (0 + 2 * 100 + 2 * 101) / 5
        assert!((is.avg - 80.4).abs() < 1e-12);
    }

    #[test]
    fn insert_size_main_bulk_cut() {
        let mut h = HashMap::new();
        h.insert(200u64, [198, 198, 0, 0]);
        h.insert(8000u64, [2, 2, 0, 0]);
        let is = insert_size_stats(&h);
        // 99 of 100 pairs is not > 0.99, so the cap bin is included.
        assert_eq!(is.rows.len(), 8001);
        h.insert(200u64, [398, 398, 0, 0]);
        let is = insert_size_stats(&h);
        assert_eq!(is.rows.len(), 201);
        assert_eq!(is.avg, 200.0);
        assert_eq!(is.sd, 0.0);
    }
}
