//! Writers for the protein sequence QC outputs.
//!
//! Two files. The first reproduces `seqkit stats -a -T` so that anything
//! already parsing seqkit output keeps working. The second carries what seqkit
//! does not report: composition and defects.

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};

use super::defects::Defects;
use super::stats::SequenceStats;

/// Residues reported in the composition table, in the conventional order.
const REPORTED_RESIDUES: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

/// Write the `seqkit stats -a -T` compatible table.
///
/// The four sequencing-only columns, `Q20(%)`, `Q30(%)`, `AvgQual` and
/// `sum_n`, are written as seqkit writes them for a protein FASTA: zero.
/// `GC(%)` likewise, since GC content is meaningless for residues.
pub fn write_seqkit_stats(entries: &[(String, SequenceStats)], path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .map(std::io::BufWriter::new)
        .with_context(|| format!("Failed to create sequence stats: {}", path.display()))?;

    writeln!(
        out,
        "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\tQ1\tQ2\tQ3\t\
         sum_gap\tN50\tN50_num\tQ20(%)\tQ30(%)\tAvgQual\tGC(%)\tsum_n"
    )?;
    for (file, stats) in entries {
        writeln!(
            out,
            "{}\tFASTA\tProtein\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t0\t0.00\t0.00\t0",
            file,
            stats.count,
            stats.total,
            stats.min,
            format_mean(stats.mean),
            stats.max,
            stats.q1,
            stats.q2,
            stats.q3,
            stats.gaps,
            stats.n50,
            stats.n50_num,
        )?;
    }
    out.flush()?;
    Ok(())
}

/// Mean length as seqkit prints it: one decimal, and a trailing `.0` kept.
fn format_mean(value: f64) -> String {
    format!("{value:.1}")
}

/// Write the composition and defect report.
pub fn write_report(
    file: &str,
    stats: &SequenceStats,
    defects: &Defects,
    path: &Path,
) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .map(std::io::BufWriter::new)
        .with_context(|| format!("Failed to create sequence report: {}", path.display()))?;

    writeln!(out, "# RustQC protein sequence report")?;
    writeln!(out, "# file\t{file}")?;
    writeln!(out)?;

    writeln!(out, "## Composition")?;
    writeln!(out, "residue\tcount\tfraction")?;
    for residue in REPORTED_RESIDUES {
        let count = stats.composition.get(residue).copied().unwrap_or(0);
        writeln!(
            out,
            "{}\t{}\t{:.6}",
            *residue as char,
            count,
            stats.fraction(*residue)
        )?;
    }
    let other: u64 = stats
        .composition
        .iter()
        .filter(|(residue, _)| !REPORTED_RESIDUES.contains(residue))
        .map(|(_, count)| count)
        .sum();
    writeln!(out, "other\t{other}\t{:.6}", {
        if stats.total == 0 {
            0.0
        } else {
            other as f64 / stats.total as f64
        }
    })?;
    writeln!(out)?;

    writeln!(out, "## Defects")?;
    writeln!(out, "category\tcount")?;
    writeln!(out, "internal_stop\t{}", defects.internal_stops.len())?;
    writeln!(out, "non_standard_residue\t{}", defects.non_standard.len())?;
    writeln!(out, "duplicate_sequence\t{}", defects.duplicates.len())?;
    writeln!(out, "duplicate_id\t{}", defects.duplicate_ids.len())?;
    writeln!(out, "empty_sequence\t{}", defects.empty.len())?;
    writeln!(
        out,
        "missing_terminal_stop\t{}",
        defects.missing_terminal_stop.len()
    )?;
    writeln!(out, "ambiguous_residues\t{}", defects.ambiguous_residues)?;
    writeln!(out)?;

    if !defects.is_clean() {
        writeln!(out, "## Offending sequences")?;
        writeln!(out, "category\tid\tdetail")?;
        for id in &defects.internal_stops {
            writeln!(out, "internal_stop\t{id}\t")?;
        }
        for (id, residue) in &defects.non_standard {
            writeln!(out, "non_standard_residue\t{id}\t{}", *residue as char)?;
        }
        for (id, first) in &defects.duplicates {
            writeln!(out, "duplicate_sequence\t{id}\t{first}")?;
        }
        for id in &defects.duplicate_ids {
            writeln!(out, "duplicate_id\t{id}\t")?;
        }
        for id in &defects.empty {
            writeln!(out, "empty_sequence\t{id}\t")?;
        }
        for id in &defects.missing_terminal_stop {
            writeln!(out, "missing_terminal_stop\t{id}\t")?;
        }
        writeln!(out)?;
    }

    out.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::protein::sequence::Record;

    fn scratch(name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join("rustqc-protein-output-tests");
        std::fs::create_dir_all(&dir).unwrap();
        dir.join(name)
    }

    #[test]
    fn the_mean_keeps_one_decimal_even_when_whole() {
        assert_eq!(format_mean(328.4), "328.4");
        assert_eq!(format_mean(100.0), "100.0", "seqkit writes 100.0, not 100");
    }

    #[test]
    fn the_stats_table_carries_seqkits_zero_columns() {
        let records = vec![Record {
            id: "a".into(),
            residues: b"MKTAYI".to_vec(),
        }];
        let stats = SequenceStats::from_records(&records);
        let path = scratch("stats.tsv");
        write_seqkit_stats(&[("f.fa".to_string(), stats)], &path).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let row = text.lines().nth(1).unwrap();
        let fields: Vec<&str> = row.split('\t').collect();
        assert_eq!(fields[1], "FASTA");
        assert_eq!(fields[2], "Protein");
        assert_eq!(&fields[14..], &["0", "0", "0.00", "0.00", "0"]);
    }

    #[test]
    fn a_clean_report_omits_the_offenders_section() {
        let records = vec![Record {
            id: "a".into(),
            residues: b"MKTAYI".to_vec(),
        }];
        let stats = SequenceStats::from_records(&records);
        let defects = crate::protein::sequence::defects::inspect(&records, false);
        let path = scratch("clean.txt");
        write_report("f.fa", &stats, &defects, &path).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(text.contains("## Composition"));
        assert!(text.contains("## Defects"));
        assert!(
            !text.contains("## Offending sequences"),
            "nothing to list, so the section is left out"
        );
    }
}
