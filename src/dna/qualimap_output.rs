//! Writers for the Qualimap `bamqc` outputs.
//!
//! The formats are reproduced from Qualimap 2.3's own output. Two details are
//! easy to miss: integers carry thousands separators, and the "Mismatches and
//! indels" section is indented by four spaces where every other section uses
//! five.
//!
//! # Figures that do not match exactly
//!
//! - `mean mapping quality` differs in the fourth decimal, 2.4179 against
//!   2.4178 on the project fixture. It is the mean of the per-window means;
//!   393 of the 397 windows match exactly and the four that do not differ by
//!   at most 0.053, which is consistent with Qualimap accumulating them
//!   differently at window boundaries.
//! - `std coverageData` differs in the fourth decimal, 154.9340 against
//!   154.9323, for the same reason.
//! - `homopolymer indels` is computed here as an indel flanked by a run of
//!   four identical bases. Qualimap's own definition was not recovered: no
//!   combination of run length from two to five, read orientation or direction
//!   reproduces its split of 7 homopolymer against 5 other indels, so this
//!   figure differs.
//! - The coverage histogram differs in 10 bins of roughly 590, always by one
//!   base and always between adjacent bins, so about five reference positions
//!   out of 40001 sit one deeper here than in Qualimap. That carries into the
//!   `coverageData >= NX` lines, which agree to within 0.003 percentage
//!   points.

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};

use super::qualimap::ContigQualimap;

/// Format an integer with thousands separators, as Qualimap does.
fn thousands(value: u64) -> String {
    let digits = value.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// Format a percentage rounded to `places` decimals, trailing zeros removed.
fn trimmed(value: f64, places: usize) -> String {
    let text = format!("{value:.places$}");
    if text.contains('.') {
        text.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        text
    }
}

/// Percentage of `part` in `whole`, guarding against an empty denominator.
fn pct(part: u64, whole: u64) -> f64 {
    if whole == 0 {
        0.0
    } else {
        // Divide before multiplying, as Qualimap does: the other order moves
        // the last two digits of the printed double.
        part as f64 / whole as f64 * 100.0
    }
}

/// Format a double the way Java's `Double.toString` does, which is what
/// Qualimap's tables carry: the shortest representation that round-trips, but
/// always with at least one digit after the point, so `0` is written `0.0`.
fn java_double(value: f64) -> String {
    let text = format!("{value}");
    if text.contains('.') || text.contains('e') || text.contains("NaN") || text.contains("inf") {
        text
    } else {
        format!("{text}.0")
    }
}

/// Write `genome_results.txt`.
pub fn write_genome_results(
    contigs: &[ContigQualimap],
    bam_path: &str,
    outfile: &Path,
) -> Result<()> {
    let mut out = std::fs::File::create(outfile)
        .map(std::io::BufWriter::new)
        .with_context(|| format!("Failed to create genome results: {}", outfile.display()))?;

    let total_length: u64 = contigs.iter().map(|c| c.length).sum();
    let mut counters = super::qualimap::QualimapCounters::default();
    for contig in contigs {
        counters.merge(&contig.counters);
    }
    let windows: usize = contigs.iter().map(|c| c.midpoints.len()).sum();

    writeln!(out, "BamQC report")?;
    writeln!(out, "-----------------------------------")?;
    writeln!(out)?;
    writeln!(out, ">>>>>>> Input")?;
    writeln!(out)?;
    writeln!(out, "     bam file = {bam_path}")?;
    writeln!(out, "     outfile = {}", outfile.display())?;
    writeln!(out)?;
    writeln!(out)?;

    writeln!(out, ">>>>>>> Reference")?;
    writeln!(out)?;
    writeln!(out, "     number of bases = {} bp", thousands(total_length))?;
    writeln!(out, "     number of contigs = {}", contigs.len())?;
    writeln!(out)?;
    writeln!(out)?;

    writeln!(out, ">>>>>>> Globals")?;
    writeln!(out)?;
    writeln!(out, "     number of windows = {windows}")?;
    writeln!(out)?;
    writeln!(out, "     number of reads = {}", thousands(counters.reads))?;
    writeln!(
        out,
        "     number of mapped reads = {} ({}%)",
        thousands(counters.mapped),
        trimmed(pct(counters.mapped, counters.reads), 2)
    )?;
    writeln!(
        out,
        "     number of secondary alignments = {}",
        thousands(counters.secondary)
    )?;
    writeln!(out)?;
    writeln!(
        out,
        "     number of mapped paired reads (first in pair) = {}",
        thousands(counters.paired_first)
    )?;
    writeln!(
        out,
        "     number of mapped paired reads (second in pair) = {}",
        thousands(counters.paired_second)
    )?;
    writeln!(
        out,
        "     number of mapped paired reads (both in pair) = {}",
        thousands(counters.paired_both)
    )?;
    writeln!(
        out,
        "     number of mapped paired reads (singletons) = {}",
        thousands(counters.singletons)
    )?;
    writeln!(out)?;
    writeln!(
        out,
        "     number of mapped bases = {} bp",
        thousands(counters.mapped_bases)
    )?;
    writeln!(
        out,
        "     number of sequenced bases = {} bp",
        thousands(counters.sequenced_bases)
    )?;
    // Qualimap reports this only when run with a reference; without one it is
    // zero, which is what RustQC always is here.
    writeln!(out, "     number of aligned bases = 0 bp")?;
    writeln!(
        out,
        "     number of duplicated reads (flagged) = {}",
        thousands(counters.duplicates)
    )?;
    writeln!(out)?;
    writeln!(out)?;

    let (insert_mean, insert_sd, insert_median) = counters.insert_size_stats();
    writeln!(out, ">>>>>>> Insert size")?;
    writeln!(out)?;
    writeln!(out, "     mean insert size = {insert_mean:.4}")?;
    writeln!(out, "     std insert size = {insert_sd:.4}")?;
    writeln!(out, "     median insert size = {insert_median}")?;
    writeln!(out)?;
    writeln!(out)?;

    let mean_mapq = if contigs.is_empty() {
        0.0
    } else {
        contigs
            .iter()
            .map(|c| c.mean_mapping_quality())
            .sum::<f64>()
            / contigs.len() as f64
    };
    writeln!(out, ">>>>>>> Mapping quality")?;
    writeln!(out)?;
    writeln!(out, "     mean mapping quality = {mean_mapq:.4}")?;
    writeln!(out)?;
    writeln!(out)?;

    let bases: u64 = counters.base_counts.iter().sum();
    writeln!(out, ">>>>>>> ACTG content")?;
    writeln!(out)?;
    for (label, index) in [("A", 0), ("C", 1), ("T", 3), ("G", 2), ("N", 4)] {
        writeln!(
            out,
            "     number of {label}'s = {} bp ({}%)",
            thousands(counters.base_counts[index]),
            trimmed(pct(counters.base_counts[index], bases), 2)
        )?;
    }
    writeln!(out)?;
    let gc = counters.base_counts[1] + counters.base_counts[2];
    writeln!(out, "     GC percentage = {}%", trimmed(pct(gc, bases), 2))?;
    writeln!(out)?;
    writeln!(out)?;

    // Note the four-space indent: this section is the odd one out.
    writeln!(out, ">>>>>>> Mismatches and indels")?;
    writeln!(out)?;
    writeln!(
        out,
        "    general error rate = {}",
        trimmed(counters.general_error_rate(), 4)
    )?;
    writeln!(
        out,
        "    number of mismatches = {}",
        thousands(counters.mismatches())
    )?;
    writeln!(
        out,
        "    number of insertions = {}",
        thousands(counters.insertions)
    )?;
    writeln!(
        out,
        "    mapped reads with insertion percentage = {}%",
        trimmed(pct(counters.reads_with_insertion, counters.mapped), 2)
    )?;
    writeln!(
        out,
        "    number of deletions = {}",
        thousands(counters.deletions)
    )?;
    writeln!(
        out,
        "    mapped reads with deletion percentage = {}%",
        trimmed(pct(counters.reads_with_deletion, counters.mapped), 2)
    )?;
    writeln!(
        out,
        "    homopolymer indels = {}%",
        trimmed(100.0 * counters.homopolymer_fraction(), 2)
    )?;
    writeln!(out)?;
    writeln!(out)?;

    let mean_coverage = if total_length == 0 {
        0.0
    } else {
        contigs
            .iter()
            .map(|c| c.coverage.iter().map(|d| f64::from(*d)).sum::<f64>())
            .sum::<f64>()
            / total_length as f64
    };
    let coverage_sd = {
        let variance = contigs
            .iter()
            .flat_map(|c| c.coverage.iter())
            .map(|d| {
                let diff = f64::from(*d) - mean_coverage;
                diff * diff
            })
            .sum::<f64>()
            / total_length.max(1) as f64;
        variance.sqrt()
    };

    writeln!(out, ">>>>>>> Coverage")?;
    writeln!(out)?;
    writeln!(out, "     mean coverageData = {mean_coverage:.4}X")?;
    writeln!(out, "     std coverageData = {coverage_sd:.4}X")?;
    writeln!(out)?;
    for (level, fraction) in genome_fraction(contigs, total_length) {
        writeln!(
            out,
            "     There is a {}% of reference with a coverageData >= {level}X",
            trimmed(fraction, 2)
        )?;
    }
    writeln!(out)?;
    writeln!(out)?;

    writeln!(out, ">>>>>>> Coverage per contig")?;
    writeln!(out)?;
    for contig in contigs {
        let covered: u64 = contig.coverage.iter().map(|d| u64::from(*d)).sum();
        writeln!(
            out,
            "\t{}\t{}\t{}\t{}\t{}",
            contig.name,
            contig.length,
            covered,
            contig.mean_coverage(),
            contig.coverage_sd()
        )?;
    }
    writeln!(out)?;
    writeln!(out)?;

    out.flush()?;
    Ok(())
}

/// Percentage of the whole reference at or above each level from 1 to 51.
fn genome_fraction(contigs: &[ContigQualimap], total_length: u64) -> Vec<(u32, f64)> {
    (1..=51)
        .map(|level| {
            let at_or_above: u64 = contigs
                .iter()
                .map(|c| c.coverage.iter().filter(|d| **d >= level).count() as u64)
                .sum();
            (level, pct(at_or_above, total_length))
        })
        .collect()
}

/// Write the twelve `raw_data_qualimapReport` tables RustQC reproduces.
///
/// Two of Qualimap's tables are not written: its GC content distribution is
/// computed over a 679-read subsample whose selection rule is not documented
/// and could not be recovered from the output, and its duplication rate
/// histogram uses a definition that does not match a read-start-position
/// count. Emitting a table under the same name with different numbers would be
/// worse than leaving it out.
pub fn write_raw_data(contigs: &[ContigQualimap], dir: &Path) -> Result<()> {
    std::fs::create_dir_all(dir)
        .with_context(|| format!("Failed to create raw data directory: {}", dir.display()))?;

    let mut counters = super::qualimap::QualimapCounters::default();
    for contig in contigs {
        counters.merge(&contig.counters);
    }

    // Per-window tables, positions given as window midpoints.
    table(
        dir,
        "coverage_across_reference.txt",
        "#Position (bp)\tCoverage\tStd",
        |out| {
            for contig in contigs {
                for i in 0..contig.midpoints.len() {
                    writeln!(
                        out,
                        "{}\t{}\t{}",
                        java_double(contig.midpoints[i]),
                        java_double(contig.window_coverage[i]),
                        java_double(contig.window_coverage_sd[i])
                    )?;
                }
            }
            Ok(())
        },
    )?;

    table(
        dir,
        "mapping_quality_across_reference.txt",
        "#Position (bp)\tmapping quality",
        |out| {
            for contig in contigs {
                for i in 0..contig.midpoints.len() {
                    writeln!(
                        out,
                        "{}\t{}",
                        java_double(contig.midpoints[i]),
                        java_double(contig.window_mapq[i])
                    )?;
                }
            }
            Ok(())
        },
    )?;

    table(
        dir,
        "insert_size_across_reference.txt",
        "#Position (bp)\tinsert size",
        |out| {
            for contig in contigs {
                for i in 0..contig.midpoints.len() {
                    writeln!(
                        out,
                        "{}\t{}",
                        java_double(contig.midpoints[i]),
                        java_double(contig.window_insert[i])
                    )?;
                }
            }
            Ok(())
        },
    )?;

    // Histograms.
    let mut coverage_histogram = std::collections::BTreeMap::new();
    let mut mapq_histogram = std::collections::BTreeMap::new();
    for contig in contigs {
        for (depth, count) in &contig.coverage_histogram {
            *coverage_histogram.entry(*depth).or_insert(0u64) += count;
        }
        for (quality, count) in &contig.mapq_histogram {
            *mapq_histogram.entry(*quality).or_insert(0u64) += count;
        }
    }

    table(
        dir,
        "coverage_histogram.txt",
        "#Coverage\tNumber of genomic locations",
        |out| {
            for (depth, count) in &coverage_histogram {
                writeln!(
                    out,
                    "{}\t{}",
                    java_double(*depth as f64),
                    java_double(*count as f64)
                )?;
            }
            Ok(())
        },
    )?;

    table(
        dir,
        "mapping_quality_histogram.txt",
        "#Mapping quality\tmapping quality",
        |out| {
            for (quality, count) in &mapq_histogram {
                writeln!(
                    out,
                    "{}\t{}",
                    java_double(*quality as f64),
                    java_double(*count as f64)
                )?;
            }
            Ok(())
        },
    )?;

    table(
        dir,
        "insert_size_histogram.txt",
        "#Insert size (bp)\tinsert size",
        |out| {
            for (size, count) in &counters.insert_sizes {
                writeln!(
                    out,
                    "{}\t{}",
                    java_double(*size as f64),
                    java_double(*count as f64)
                )?;
            }
            Ok(())
        },
    )?;

    let total_length: u64 = contigs.iter().map(|c| c.length).sum();
    table(
        dir,
        "genome_fraction_coverage.txt",
        "#Coverage (X)\tCoverage",
        |out| {
            for (level, fraction) in genome_fraction(contigs, total_length) {
                writeln!(
                    out,
                    "{}\t{}",
                    java_double(level as f64),
                    java_double(fraction)
                )?;
            }
            Ok(())
        },
    )?;

    table(
        dir,
        "mapped_reads_clipping_profile.txt",
        "#Read position (bp)\tClipping profile",
        |out| {
            for (position, count) in counters.clipping_by_position.iter().enumerate() {
                writeln!(
                    out,
                    "{}\t{}",
                    java_double(position as f64),
                    java_double(pct(*count, counters.clipped_bases))
                )?;
            }
            Ok(())
        },
    )?;

    table(
        dir,
        "mapped_reads_nucleotide_content.txt",
        "# Position (bp)\tA\tC\tG\tT\tN",
        |out| {
            for (position, counts) in counters.nucleotide_by_position.iter().enumerate() {
                let total: u64 = counts.iter().sum();
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    java_double(position as f64),
                    java_double(pct(counts[0], total)),
                    java_double(pct(counts[1], total)),
                    java_double(pct(counts[2], total)),
                    java_double(pct(counts[3], total)),
                    java_double(pct(counts[4], total)),
                )?;
            }
            Ok(())
        },
    )?;

    table(
        dir,
        "homopolymer_indels.txt",
        "#Type of indel\tNumber of indels",
        |out| {
            for (label, index) in [
                ("polyA", 0),
                ("polyC", 1),
                ("polyG", 2),
                ("polyT", 3),
                ("polyN", 4),
            ] {
                writeln!(out, "{label}\t{}", counters.homopolymer_indels[index])?;
            }
            writeln!(out, "Non-poly\t{}", counters.non_polymer_indels)?;
            Ok(())
        },
    )?;

    Ok(())
}

/// Write one raw data table with its header line.
fn table<F>(dir: &Path, name: &str, header: &str, body: F) -> Result<()>
where
    F: FnOnce(&mut dyn Write) -> Result<()>,
{
    let path = dir.join(name);
    let mut out = std::fs::File::create(&path)
        .map(std::io::BufWriter::new)
        .with_context(|| format!("Failed to create {}", path.display()))?;
    writeln!(out, "{header}")?;
    body(&mut out)?;
    out.flush()?;
    Ok(())
}

/// Write `qualimapReport.html`.
///
/// This is RustQC's own summary page rather than a copy of Qualimap's, which
/// ships a bundle of images, CSS and JavaScript. The numbers are the same ones
/// `genome_results.txt` carries; the page exists so a run has something
/// readable to open, and the raw tables remain the machine-readable source.
pub fn write_html_report(contigs: &[ContigQualimap], sample_name: &str, path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .map(std::io::BufWriter::new)
        .with_context(|| format!("Failed to create the report: {}", path.display()))?;

    let total_length: u64 = contigs.iter().map(|c| c.length).sum();
    let mut counters = super::qualimap::QualimapCounters::default();
    for contig in contigs {
        counters.merge(&contig.counters);
    }
    let mean_coverage = if total_length == 0 {
        0.0
    } else {
        contigs
            .iter()
            .map(|c| c.coverage.iter().map(|d| f64::from(*d)).sum::<f64>())
            .sum::<f64>()
            / total_length as f64
    };
    let (insert_mean, insert_sd, insert_median) = counters.insert_size_stats();

    writeln!(out, "<!doctype html>")?;
    writeln!(out, "<html lang=\"en\"><head><meta charset=\"utf-8\">")?;
    writeln!(out, "<title>BamQC report: {}</title>", escape(sample_name))?;
    writeln!(
        out,
        "<style>body{{font:14px/1.5 system-ui,sans-serif;margin:2rem;max-width:60rem}}         table{{border-collapse:collapse;margin:1rem 0}}         th,td{{border:1px solid #ccc;padding:.35rem .6rem;text-align:left}}         th{{background:#f4f4f4}}td.n{{text-align:right;font-variant-numeric:tabular-nums}}         </style></head><body>"
    )?;
    writeln!(out, "<h1>BamQC report</h1>")?;
    writeln!(
        out,
        "<p>Sample: <strong>{}</strong></p>",
        escape(sample_name)
    )?;

    let rows: Vec<(&str, String)> = vec![
        ("Reference bases", thousands(total_length)),
        ("Contigs", contigs.len().to_string()),
        ("Reads", thousands(counters.reads)),
        ("Mapped reads", thousands(counters.mapped)),
        ("Duplicated reads (flagged)", thousands(counters.duplicates)),
        ("Mapped bases", thousands(counters.mapped_bases)),
        ("Sequenced bases", thousands(counters.sequenced_bases)),
        ("Mean coverage", format!("{mean_coverage:.4}X")),
        ("Mean insert size", format!("{insert_mean:.4}")),
        ("Std insert size", format!("{insert_sd:.4}")),
        ("Median insert size", insert_median.to_string()),
        ("Mismatches", thousands(counters.mismatches())),
        ("Insertions", thousands(counters.insertions)),
        ("Deletions", thousands(counters.deletions)),
    ];
    writeln!(out, "<h2>Summary</h2><table><tbody>")?;
    for (label, value) in rows {
        writeln!(
            out,
            "<tr><th scope=\"row\">{label}</th><td class=\"n\">{value}</td></tr>"
        )?;
    }
    writeln!(out, "</tbody></table>")?;

    writeln!(out, "<h2>Coverage per contig</h2>")?;
    writeln!(
        out,
        "<table><thead><tr><th>Contig</th><th>Length</th><th>Mapped bases</th>         <th>Mean coverage</th><th>Std</th></tr></thead><tbody>"
    )?;
    for contig in contigs {
        let covered: u64 = contig.coverage.iter().map(|d| u64::from(*d)).sum();
        writeln!(
            out,
            "<tr><td>{}</td><td class=\"n\">{}</td><td class=\"n\">{}</td>             <td class=\"n\">{:.4}</td><td class=\"n\">{:.4}</td></tr>",
            escape(&contig.name),
            thousands(contig.length),
            thousands(covered),
            contig.mean_coverage(),
            contig.coverage_sd(),
        )?;
    }
    writeln!(out, "</tbody></table>")?;
    writeln!(
        out,
        "<p>Per-window and per-position tables are in \
         <code>raw_data_qualimapReport/</code>.</p>"
    )?;
    writeln!(out, "</body></html>")?;

    out.flush()?;
    Ok(())
}

/// Escape the few characters that would otherwise close a tag or attribute.
fn escape(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn doubles_are_written_the_way_java_writes_them() {
        assert_eq!(java_double(0.0), "0.0");
        assert_eq!(java_double(51.0), "51.0");
        assert_eq!(java_double(2.5), "2.5");
        assert_eq!(java_double(2.9524261893452746), "2.9524261893452746");
    }

    #[test]
    fn thousands_separators_match_qualimaps_formatting() {
        assert_eq!(thousands(0), "0");
        assert_eq!(thousands(999), "999");
        assert_eq!(thousands(1_000), "1,000");
        assert_eq!(thousands(40_001), "40,001");
        assert_eq!(thousands(670_999), "670,999");
    }

    #[test]
    fn percentages_drop_trailing_zeros() {
        assert_eq!(trimmed(2.95, 2), "2.95");
        assert_eq!(trimmed(2.50, 2), "2.5");
        assert_eq!(trimmed(2.0, 2), "2");
        assert_eq!(trimmed(15.2, 2), "15.2");
    }

    #[test]
    fn html_escaping_covers_the_characters_that_break_markup() {
        assert_eq!(escape("a<b>c&d\"e"), "a&lt;b&gt;c&amp;d&quot;e");
        assert_eq!(escape("plain"), "plain");
    }

    #[test]
    fn a_zero_denominator_gives_zero_rather_than_a_nan() {
        assert_eq!(pct(5, 0), 0.0);
        assert_eq!(pct(0, 10), 0.0);
    }
}
