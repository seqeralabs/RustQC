//! Writers for the read QC outputs.
//!
//! Placeholder until the FastQC-compatible writer lands; the metrics are
//! reported through the seqkit-compatible table for now.

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};

use super::metrics::ReadMetrics;

/// Write the `seqkit stats -a -T` compatible table.
pub fn write_seqkit_stats(entries: &[(String, ReadMetrics)], path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .map(std::io::BufWriter::new)
        .with_context(|| format!("Failed to create read stats: {}", path.display()))?;

    writeln!(
        out,
        "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\tQ1\tQ2\tQ3\t\
         sum_gap\tN50\tN50_num\tQ20(%)\tQ30(%)\tAvgQual\tGC(%)\tsum_n"
    )?;
    for (file, m) in entries {
        let lengths = expand_lengths(m);
        let (q1, q2, q3) = tukey_quartiles(&lengths);
        let (n50, n50_num) = n50(&lengths, m.bases);
        writeln!(
            out,
            "{}\tFASTQ\tDNA\t{}\t{}\t{}\t{:.1}\t{}\t{}\t{}\t{}\t0\t{}\t{}\t{:.0}\t{:.0}\t{:.2}\t{:.2}\t{}",
            file,
            m.reads,
            m.bases,
            m.min_length(),
            m.mean_length(),
            m.max_length(),
            q1,
            q2,
            q3,
            n50,
            n50_num,
            m.q20_percent(),
            m.q30_percent(),
            m.average_quality(),
            m.gc_percent(),
            m.n_bases,
        )?;
    }
    out.flush()?;
    Ok(())
}

/// The sorted read lengths, rebuilt from the histogram.
fn expand_lengths(m: &ReadMetrics) -> Vec<u64> {
    let mut lengths = Vec::with_capacity(m.reads as usize);
    for (length, count) in &m.length_histogram {
        lengths.extend(std::iter::repeat_n(*length, *count as usize));
    }
    lengths
}

/// Median of an ascending slice, rounded half to even, as seqkit does.
fn median(sorted: &[u64]) -> u64 {
    match sorted.len() {
        0 => 0,
        n if n % 2 == 1 => sorted[n / 2],
        n => {
            let value = (sorted[n / 2 - 1] + sorted[n / 2]) as f64 / 2.0;
            let floor = value.floor();
            if (value - floor - 0.5).abs() < f64::EPSILON {
                if (floor as i64) % 2 == 0 {
                    floor as u64
                } else {
                    floor as u64 + 1
                }
            } else {
                value.round() as u64
            }
        }
    }
}

/// Tukey's quartiles, as seqkit reports them.
fn tukey_quartiles(sorted: &[u64]) -> (u64, u64, u64) {
    let n = sorted.len();
    if n == 0 {
        return (0, 0, 0);
    }
    let half = n / 2;
    (
        median(&sorted[..half]),
        median(sorted),
        median(&sorted[n - half..]),
    )
}

/// N50 and the number of *distinct lengths* needed to reach it.
///
/// The second figure is not the number of sequences, which is the obvious
/// reading and the wrong one. seqkit walks the distinct lengths from longest
/// down and counts how many it consumed: three sequences of lengths 10, 10 and
/// 3 give `N50_num` of 1, not 2, because the two tens are one length. A file
/// whose lengths are all distinct hides the difference entirely, which is how
/// this went unnoticed against a protein fixture.
fn n50(sorted_ascending: &[u64], total: u64) -> (u64, u64) {
    let mut cumulative = 0u64;
    let mut distinct = 0u64;
    let mut previous: Option<u64> = None;
    for length in sorted_ascending.iter().rev() {
        if previous != Some(*length) {
            distinct += 1;
            previous = Some(*length);
        }
        cumulative += length;
        if cumulative * 2 >= total {
            return (*length, distinct);
        }
    }
    (0, 0)
}

/// Write a FastQC-compatible `fastqc_data.txt`.
///
/// Modules RustQC computes are written with their data; the rest are omitted
/// rather than emitted empty, so a parser sees either real numbers or nothing.
///
/// Every module carries a pass, warn or fail verdict in FastQC's output. Those
/// verdicts encode FastQC's own thresholds, which are not part of the data, so
/// RustQC writes `pass` uniformly and leaves judgement to the reader. That is
/// a deliberate difference and the one place this file is not a drop-in
/// replacement.
pub fn write_fastqc_data(file: &str, m: &ReadMetrics, path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .map(std::io::BufWriter::new)
        .with_context(|| format!("Failed to create FastQC data: {}", path.display()))?;

    writeln!(out, "##FastQC\t0.12.1")?;

    writeln!(out, ">>Basic Statistics\tpass")?;
    writeln!(out, "#Measure\tValue")?;
    writeln!(out, "Filename\t{file}")?;
    writeln!(out, "File type\tConventional base calls")?;
    writeln!(out, "Encoding\tSanger / Illumina 1.9")?;
    writeln!(out, "Total Sequences\t{}", m.reads)?;
    writeln!(out, "Total Bases\t{}", format_bases(m.bases))?;
    writeln!(out, "Sequences flagged as poor quality\t0")?;
    if m.min_length() == m.max_length() {
        writeln!(out, "Sequence length\t{}", m.min_length())?;
    } else {
        writeln!(
            out,
            "Sequence length\t{}-{}",
            m.min_length(),
            m.max_length()
        )?;
    }
    writeln!(out, "%GC\t{}", m.gc_percent().round() as u64)?;
    writeln!(out, ">>END_MODULE")?;

    writeln!(out, ">>Per base sequence quality\tpass")?;
    writeln!(
        out,
        "#Base\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile"
    )?;
    for bin in m.position_bins() {
        let q = m.position_quality(bin);
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            format_bin(bin),
            q.mean,
            java_double(q.median),
            java_double(q.lower_quartile),
            java_double(q.upper_quartile),
            java_double(q.tenth),
            java_double(q.ninetieth),
        )?;
    }
    writeln!(out, ">>END_MODULE")?;

    writeln!(out, ">>Per sequence quality scores\tpass")?;
    writeln!(out, "#Quality\tCount")?;
    if let (Some(lowest), Some(highest)) = (
        m.read_quality_histogram.keys().next(),
        m.read_quality_histogram.keys().next_back(),
    ) {
        for score in *lowest..=*highest {
            let count = m.read_quality_histogram.get(&score).copied().unwrap_or(0);
            writeln!(out, "{score}\t{}", java_double(count as f64))?;
        }
    }
    writeln!(out, ">>END_MODULE")?;

    writeln!(out, ">>Per base sequence content\tpass")?;
    writeln!(out, "#Base\tG\tA\tT\tC")?;
    for bin in m.position_bins() {
        let [a, c, g, t] = m.position_composition(bin);
        writeln!(out, "{}\t{g}\t{a}\t{t}\t{c}", format_bin(bin))?;
    }
    writeln!(out, ">>END_MODULE")?;

    // FastQC spreads each read's GC contribution across neighbouring bins so
    // that a coarse discrete distribution plots smoothly, which makes its
    // counts fractional. This is the plain rounded distribution instead: a
    // different figure, reported under the same name because that is what
    // parsers look for, and called out in the documentation.
    writeln!(out, ">>Per sequence GC content\tpass")?;
    writeln!(out, "#GC Content\tCount")?;
    for percent in 0u8..=100 {
        let count = m.gc_histogram.get(&percent).copied().unwrap_or(0);
        writeln!(out, "{percent}\t{}", java_double(count as f64))?;
    }
    writeln!(out, ">>END_MODULE")?;

    writeln!(out, ">>Per base N content\tpass")?;
    writeln!(out, "#Base\tN-Count")?;
    for bin in m.position_bins() {
        writeln!(out, "{}\t{}", format_bin(bin), m.position_n_content(bin))?;
    }
    writeln!(out, ">>END_MODULE")?;

    out.flush()?;
    Ok(())
}

/// A position bin as FastQC labels it: a bare number when it spans one
/// position, a range otherwise.
fn format_bin(bin: (u64, u64)) -> String {
    if bin.0 == bin.1 {
        bin.0.to_string()
    } else {
        format!("{}-{}", bin.0, bin.1)
    }
}

/// Total bases as FastQC writes them, in whichever unit keeps it short.
fn format_bases(bases: u64) -> String {
    const UNITS: [(u64, &str); 3] = [(1_000_000_000, "Gbp"), (1_000_000, "Mbp"), (1_000, "kbp")];
    for (scale, unit) in UNITS {
        if bases >= scale {
            return format!("{:.1} {unit}", bases as f64 / scale as f64);
        }
    }
    format!("{bases} bp")
}

/// Format a double the way Java does, so an integral value keeps its `.0`.
fn java_double(value: f64) -> String {
    let text = format!("{value}");
    if text.contains('.') || text.contains('e') {
        text
    } else {
        format!("{text}.0")
    }
}
