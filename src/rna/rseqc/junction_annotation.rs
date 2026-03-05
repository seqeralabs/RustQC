//! Splice junction annotation and classification.
//!
//! Reimplementation of RSeQC's `junction_annotation.py`: extracts splice junctions
//! from BAM files (CIGAR N-operations) and classifies them as known (annotated),
//! partial novel, or complete novel by comparing against a reference gene model
//! (BED12 or GTF).

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use indexmap::IndexMap;

// ===================================================================
// Data types
// ===================================================================

/// Classification of a splice junction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum JunctionClass {
    /// Both splice sites match known annotation positions.
    Annotated,
    /// Exactly one splice site matches a known position.
    PartialNovel,
    /// Neither splice site matches any known position.
    CompleteNovel,
}

impl JunctionClass {
    /// Label string matching RSeQC output.
    pub fn label(self) -> &'static str {
        match self {
            JunctionClass::Annotated => "annotated",
            JunctionClass::PartialNovel => "partial_novel",
            JunctionClass::CompleteNovel => "complete_novel",
        }
    }

    /// RGB color for BED12 output.
    pub fn color(self) -> &'static str {
        match self {
            JunctionClass::Annotated => "205,0,0",
            JunctionClass::PartialNovel => "0,205,0",
            JunctionClass::CompleteNovel => "0,0,205",
        }
    }
}

/// A unique splice junction identified by chromosome and intron coordinates.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Junction {
    /// Chromosome name (uppercased internally).
    pub chrom: String,
    /// 0-based intron start (first intronic base).
    pub intron_start: u64,
    /// 0-based intron end (exclusive — first base after intron).
    pub intron_end: u64,
}

/// Aggregated results from junction annotation.
#[derive(Debug)]
pub struct JunctionResults {
    /// Per-junction read counts and classification, in BAM encounter order.
    ///
    /// Uses `IndexMap` to preserve insertion order, matching the output of
    /// upstream RSeQC's Python `defaultdict` (which preserves insertion order
    /// in Python 3.7+).
    pub junctions: IndexMap<Junction, (u64, JunctionClass)>,
    /// Total splicing events (per-read N-operations, includes short introns).
    pub total_events: u64,
    /// Events classified as known.
    pub known_events: u64,
    /// Events classified as partial novel.
    pub partial_novel_events: u64,
    /// Events classified as complete novel.
    pub complete_novel_events: u64,
    /// Events filtered by min_intron threshold.
    pub filtered_events: u64,
}

/// Counts of unique junctions by classification.
#[derive(Debug)]
pub struct JunctionCounts {
    /// Number of distinct known (annotated) junctions.
    pub known: u64,
    /// Number of distinct partial-novel junctions.
    pub partial_novel: u64,
    /// Number of distinct complete-novel junctions.
    pub novel: u64,
    /// Total distinct junctions (known + partial_novel + novel).
    pub total: u64,
}

impl JunctionResults {
    /// Count distinct junctions by classification.
    pub fn junction_counts(&self) -> JunctionCounts {
        let mut known: u64 = 0;
        let mut partial_novel: u64 = 0;
        let mut novel: u64 = 0;
        for (_, class) in self.junctions.values() {
            match class {
                JunctionClass::Annotated => known += 1,
                JunctionClass::PartialNovel => partial_novel += 1,
                JunctionClass::CompleteNovel => novel += 1,
            }
        }
        JunctionCounts {
            known,
            partial_novel,
            novel,
            total: known + partial_novel + novel,
        }
    }
}

// ===================================================================
// Output formatting
// ===================================================================

/// Convert an uppercased chromosome name back to lowercase for output.
///
/// Matches RSeQC's `chrom.replace("CHR","chr")` behavior.
fn format_chrom(chrom: &str) -> String {
    chrom.replace("CHR", "chr")
}

/// Write the junction annotation XLS (TSV) file.
///
/// Format matches RSeQC's `.junction.xls` output.
pub fn write_junction_xls(results: &JunctionResults, path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(path)
        .with_context(|| format!("Failed to create file: {}", path.display()))?;

    writeln!(
        f,
        "chrom\tintron_st(0-based)\tintron_end(1-based)\tread_count\tannotation"
    )?;

    // Output junctions in BAM encounter order (IndexMap preserves insertion order,
    // matching upstream RSeQC's Python defaultdict behaviour).
    for (junction, (count, class)) in &results.junctions {
        // RSeQC output has a tab then space before the annotation label
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t {}",
            format_chrom(&junction.chrom),
            junction.intron_start,
            junction.intron_end,
            count,
            class.label()
        )?;
    }

    Ok(())
}

/// Write the junction BED12 file.
///
/// Each junction is represented with 1bp exon blocks flanking the intron,
/// matching the reference output format.
pub fn write_junction_bed(results: &JunctionResults, path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(path)
        .with_context(|| format!("Failed to create file: {}", path.display()))?;

    // Output junctions in BAM encounter order (IndexMap preserves insertion order)
    for (junction, (count, class)) in &results.junctions {
        let chrom = format_chrom(&junction.chrom);
        // BED coordinates: 1bp block before intron and 1bp block after
        let bed_start = junction.intron_start.saturating_sub(1); // 1bp before intron
        let bed_end = junction.intron_end + 1; // 1bp after intron
        let span = bed_end - bed_start;

        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\t{}\t2\t1,1\t0,{}",
            chrom,
            bed_start,
            bed_end,
            class.label(),
            count,
            bed_start,
            bed_end,
            class.color(),
            span - 1
        )?;
    }

    Ok(())
}

/// Write the R script for junction pie charts.
///
/// Generates two pie charts: splice events and splice junctions.
pub fn write_junction_plot_r(results: &JunctionResults, prefix: &str, path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(path)
        .with_context(|| format!("Failed to create file: {}", path.display()))?;

    // Event-level percentages
    let total_classified_events =
        results.known_events + results.partial_novel_events + results.complete_novel_events;
    let (e_known_pct, e_partial_pct, e_novel_pct) = if total_classified_events > 0 {
        (
            results.known_events as f64 * 100.0 / total_classified_events as f64,
            results.partial_novel_events as f64 * 100.0 / total_classified_events as f64,
            results.complete_novel_events as f64 * 100.0 / total_classified_events as f64,
        )
    } else {
        (0.0, 0.0, 0.0)
    };

    // Junction-level counts
    let jc = results.junction_counts();
    let (j_known_pct, j_partial_pct, j_novel_pct) = if jc.total > 0 {
        (
            jc.known as f64 * 100.0 / jc.total as f64,
            jc.partial_novel as f64 * 100.0 / jc.total as f64,
            jc.novel as f64 * 100.0 / jc.total as f64,
        )
    } else {
        (0.0, 0.0, 0.0)
    };

    // Events pie chart
    writeln!(f, "pdf(\"{}.splice_events.pdf\")", prefix)?;
    writeln!(
        f,
        "events=c({:.3},{:.3},{:.3})",
        e_partial_pct, e_novel_pct, e_known_pct
    )?;
    writeln!(
        f,
        "pie(events,col=c(2,3,4),init.angle=30,angle=c(60,120,150),density=c(70,70,70),main=\"splicing events\",labels=c(\"partial_novel {}%\",\"complete_novel {}%\",\"known {}%\"))",
        e_partial_pct.round() as u64,
        e_novel_pct.round() as u64,
        e_known_pct.round() as u64
    )?;
    writeln!(f, "dev.off()")?;

    // Junctions pie chart
    writeln!(f, "pdf(\"{}.splice_junction.pdf\")", prefix)?;
    writeln!(
        f,
        "junction=c({:.3},{:.3},{:.3})",
        j_partial_pct, j_novel_pct, j_known_pct
    )?;
    writeln!(
        f,
        "pie(junction,col=c(2,3,4),init.angle=30,angle=c(60,120,150),density=c(70,70,70),main=\"splicing junctions\",labels=c(\"partial_novel {}%\",\"complete_novel {}%\",\"known {}%\"))",
        j_partial_pct.round() as u64,
        j_novel_pct.round() as u64,
        j_known_pct.round() as u64
    )?;
    writeln!(f, "dev.off()")?;

    Ok(())
}

/// Write the summary statistics to a file.
///
/// Format matches RSeQC's text summary output.
pub fn write_summary(results: &JunctionResults, path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(path)
        .with_context(|| format!("Failed to create file: {}", path.display()))?;
    write_summary_to(results, &mut f)?;
    Ok(())
}

/// Print the summary to stderr (matching RSeQC behavior).
pub fn print_summary(results: &JunctionResults) {
    let mut stderr = std::io::stderr();
    // Ignore write errors to stderr
    let _ = write_summary_to(results, &mut stderr);
}

/// Write the junction annotation summary to any writer.
fn write_summary_to<W: std::io::Write>(results: &JunctionResults, f: &mut W) -> Result<()> {
    let jc = results.junction_counts();

    writeln!(
        f,
        "==================================================================="
    )?;
    writeln!(f, "Total splicing  Events:\t{}", results.total_events)?;
    writeln!(f, "Known Splicing Events:\t{}", results.known_events)?;
    writeln!(
        f,
        "Partial Novel Splicing Events:\t{}",
        results.partial_novel_events
    )?;
    writeln!(
        f,
        "Novel Splicing Events:\t{}",
        results.complete_novel_events
    )?;
    writeln!(f, "Filtered Splicing Events:\t{}", results.filtered_events)?;
    writeln!(f)?;
    writeln!(f, "Total splicing  Junctions:\t{}", jc.total)?;
    writeln!(f, "Known Splicing Junctions:\t{}", jc.known)?;
    writeln!(f, "Partial Novel Splicing Junctions:\t{}", jc.partial_novel)?;
    writeln!(f, "Novel Splicing Junctions:\t{}", jc.novel)?;
    writeln!(f)?;
    writeln!(
        f,
        "==================================================================="
    )?;

    Ok(())
}

// ===================================================================
// Unit tests
// ===================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rna::rseqc::common::ReferenceJunctions;
    use std::collections::{HashMap, HashSet};

    /// Classify a junction as annotated, partial novel, or complete novel.
    /// Kept in test module — production classification is in accumulators.rs.
    fn classify_junction(
        chrom: &str,
        intron_start: u64,
        intron_end: u64,
        reference: &ReferenceJunctions,
    ) -> JunctionClass {
        let start_known = reference
            .intron_starts
            .get(chrom)
            .is_some_and(|s| s.contains(&intron_start));
        let end_known = reference
            .intron_ends
            .get(chrom)
            .is_some_and(|s| s.contains(&intron_end));

        match (start_known, end_known) {
            (true, true) => JunctionClass::Annotated,
            (false, false) => JunctionClass::CompleteNovel,
            _ => JunctionClass::PartialNovel,
        }
    }

    #[test]
    fn test_classify_junction_annotated() {
        let mut starts = HashMap::new();
        let mut ends = HashMap::new();
        starts.insert("CHR1".to_string(), HashSet::from([100, 500]));
        ends.insert("CHR1".to_string(), HashSet::from([400, 900]));
        let reference = ReferenceJunctions {
            intron_starts: starts,
            intron_ends: ends,
        };

        assert_eq!(
            classify_junction("CHR1", 100, 400, &reference),
            JunctionClass::Annotated
        );
    }

    #[test]
    fn test_classify_junction_partial_novel() {
        let mut starts = HashMap::new();
        let mut ends = HashMap::new();
        starts.insert("CHR1".to_string(), HashSet::from([100]));
        ends.insert("CHR1".to_string(), HashSet::from([400]));
        let reference = ReferenceJunctions {
            intron_starts: starts,
            intron_ends: ends,
        };

        // Start known, end unknown
        assert_eq!(
            classify_junction("CHR1", 100, 999, &reference),
            JunctionClass::PartialNovel
        );
        // Start unknown, end known
        assert_eq!(
            classify_junction("CHR1", 999, 400, &reference),
            JunctionClass::PartialNovel
        );
    }

    #[test]
    fn test_classify_junction_complete_novel() {
        let mut starts = HashMap::new();
        let mut ends = HashMap::new();
        starts.insert("CHR1".to_string(), HashSet::from([100]));
        ends.insert("CHR1".to_string(), HashSet::from([400]));
        let reference = ReferenceJunctions {
            intron_starts: starts,
            intron_ends: ends,
        };

        assert_eq!(
            classify_junction("CHR1", 999, 888, &reference),
            JunctionClass::CompleteNovel
        );
    }

    #[test]
    fn test_classify_junction_unknown_chrom() {
        let reference = ReferenceJunctions {
            intron_starts: HashMap::new(),
            intron_ends: HashMap::new(),
        };

        assert_eq!(
            classify_junction("CHRX", 100, 400, &reference),
            JunctionClass::CompleteNovel
        );
    }

    #[test]
    fn test_junction_class_labels() {
        assert_eq!(JunctionClass::Annotated.label(), "annotated");
        assert_eq!(JunctionClass::PartialNovel.label(), "partial_novel");
        assert_eq!(JunctionClass::CompleteNovel.label(), "complete_novel");
    }

    #[test]
    fn test_junction_class_colors() {
        assert_eq!(JunctionClass::Annotated.color(), "205,0,0");
        assert_eq!(JunctionClass::PartialNovel.color(), "0,205,0");
        assert_eq!(JunctionClass::CompleteNovel.color(), "0,0,205");
    }

    #[test]
    fn test_format_chrom() {
        assert_eq!(format_chrom("CHR1"), "chr1");
        assert_eq!(format_chrom("CHRX"), "chrX");
        assert_eq!(format_chrom("6"), "6"); // numeric chroms unchanged
    }
}
