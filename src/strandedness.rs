//! Strandedness inference from BAM alignment data.
//!
//! Replaces RSeQC's `infer_experiment.py` by checking read strand orientation
//! against gene strand during the existing BAM counting pass.

use crate::counting::CountResult;
use anyhow::Result;
use log::info;
use std::io::Write;
use std::path::Path;

// ============================================================================
// Strandedness inference result
// ============================================================================

/// Inferred library strandedness protocol.
#[derive(Debug, Clone, PartialEq)]
pub enum StrandProtocol {
    /// Unstranded library (sense ≈ antisense ≈ 50%)
    Unstranded,
    /// Forward-stranded / fr-secondstrand (sense >> antisense)
    /// Corresponds to dUTP, NSR, NNSR protocols
    ForwardStranded,
    /// Reverse-stranded / fr-firststrand (antisense >> sense)
    /// Corresponds to Ligation, Standard SOLiD protocols
    ReverseStranded,
    /// Could not determine (insufficient data or ambiguous)
    Undetermined,
}

impl std::fmt::Display for StrandProtocol {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            StrandProtocol::Unstranded => write!(f, "unstranded"),
            StrandProtocol::ForwardStranded => write!(f, "forward (fr-secondstrand)"),
            StrandProtocol::ReverseStranded => write!(f, "reverse (fr-firststrand)"),
            StrandProtocol::Undetermined => write!(f, "undetermined"),
        }
    }
}

// ============================================================================
// Inference logic
// ============================================================================

/// Minimum number of informative reads required for reliable inference.
const MIN_READS_FOR_INFERENCE: u64 = 100;

/// Threshold for classifying as stranded: the dominant fraction must exceed this.
const STRANDED_THRESHOLD: f64 = 0.80;

/// Threshold for classifying as unstranded: both fractions must be within this
/// range of 0.5 (i.e., between 0.5 - tolerance and 0.5 + tolerance).
const UNSTRANDED_TOLERANCE: f64 = 0.10;

/// Infer library strandedness from sense/antisense read counts.
///
/// # Arguments
/// * `sense` - Number of reads matching the gene strand (fr-secondstrand)
/// * `antisense` - Number of reads opposing the gene strand (fr-firststrand)
/// * `undetermined` - Number of reads where strand could not be determined
///
/// # Returns
/// The inferred `StrandProtocol` and the sense/antisense/undetermined fractions.
pub fn infer_strandedness(
    sense: u64,
    antisense: u64,
    undetermined: u64,
) -> (StrandProtocol, f64, f64, f64) {
    let total = sense + antisense + undetermined;
    if total < MIN_READS_FOR_INFERENCE {
        return (StrandProtocol::Undetermined, 0.0, 0.0, 0.0);
    }

    let informative = sense + antisense;
    if informative == 0 {
        return (StrandProtocol::Undetermined, 0.0, 0.0, 1.0);
    }

    let sense_frac = sense as f64 / informative as f64;
    let antisense_frac = antisense as f64 / informative as f64;
    let undetermined_frac = undetermined as f64 / total as f64;

    let protocol = if sense_frac >= STRANDED_THRESHOLD {
        StrandProtocol::ForwardStranded
    } else if antisense_frac >= STRANDED_THRESHOLD {
        StrandProtocol::ReverseStranded
    } else if (sense_frac - 0.5).abs() <= UNSTRANDED_TOLERANCE
        && (antisense_frac - 0.5).abs() <= UNSTRANDED_TOLERANCE
    {
        StrandProtocol::Unstranded
    } else {
        StrandProtocol::Undetermined
    };

    (protocol, sense_frac, antisense_frac, undetermined_frac)
}

// ============================================================================
// Output functions
// ============================================================================

/// Write strandedness inference results in RSeQC infer_experiment.py-compatible format.
///
/// Output format matches RSeQC for compatibility with downstream tools and MultiQC:
/// ```text
/// This is PairEnd Data
///
/// Fraction of reads failed to determine: 0.0123
/// Fraction of reads explained by "1++,1--,2+-,2-+": 0.0456
/// Fraction of reads explained by "1+-,1-+,2++,2--": 0.9421
/// ```
///
/// # Arguments
/// * `path` - Output file path
/// * `result` - Count result containing strandedness counts
/// * `is_paired` - Whether the library is paired-end
pub fn write_infer_experiment(path: &Path, result: &CountResult, is_paired: bool) -> Result<()> {
    let (protocol, sense_frac, antisense_frac, undetermined_frac) = infer_strandedness(
        result.strandedness_sense,
        result.strandedness_antisense,
        result.strandedness_undetermined,
    );

    let mut writer = std::io::BufWriter::new(std::fs::File::create(path)?);

    if is_paired {
        writeln!(writer, "This is PairEnd Data")?;
    } else {
        writeln!(writer, "This is SingleEnd Data")?;
    }
    writeln!(writer)?;

    writeln!(
        writer,
        "Fraction of reads failed to determine: {:.4}",
        undetermined_frac
    )?;

    if is_paired {
        writeln!(
            writer,
            "Fraction of reads explained by \"1++,1--,2+-,2-+\": {:.4}",
            sense_frac
        )?;
        writeln!(
            writer,
            "Fraction of reads explained by \"1+-,1-+,2++,2--\": {:.4}",
            antisense_frac
        )?;
    } else {
        writeln!(
            writer,
            "Fraction of reads explained by \"++,--\": {:.4}",
            sense_frac
        )?;
        writeln!(
            writer,
            "Fraction of reads explained by \"+-,-+\": {:.4}",
            antisense_frac
        )?;
    }

    writeln!(writer)?;
    writeln!(writer, "# Inferred protocol: {}", protocol)?;

    info!(
        "Strandedness inference: {} (sense={:.1}%, antisense={:.1}%, undetermined={:.1}%)",
        protocol,
        sense_frac * 100.0,
        antisense_frac * 100.0,
        undetermined_frac * 100.0
    );

    Ok(())
}

/// Write strandedness inference results as MultiQC general stats.
///
/// Produces a MultiQC-compatible TSV with YAML header comments for
/// integration into the MultiQC general stats table.
///
/// # Arguments
/// * `path` - Output file path
/// * `result` - Count result containing strandedness counts
/// * `sample_name` - Sample name for the MultiQC table
/// * `is_paired` - Whether the library is paired-end
pub fn write_strandedness_mqc(
    path: &Path,
    result: &CountResult,
    sample_name: &str,
    is_paired: bool,
) -> Result<()> {
    let (protocol, sense_frac, antisense_frac, undetermined_frac) = infer_strandedness(
        result.strandedness_sense,
        result.strandedness_antisense,
        result.strandedness_undetermined,
    );

    let mut writer = std::io::BufWriter::new(std::fs::File::create(path)?);

    // MultiQC YAML header
    writeln!(writer, "# id: 'strandedness'")?;
    writeln!(writer, "# section_name: 'Strandedness'")?;
    writeln!(
        writer,
        "# description: 'Library strandedness inference (RSeQC infer_experiment compatible)'"
    )?;
    writeln!(writer, "# plot_type: 'generalstats'")?;
    writeln!(writer, "# pconfig:")?;
    writeln!(writer, "#     sense_pct:")?;
    writeln!(writer, "#         title: 'Sense %'")?;
    writeln!(writer, "#         min: 0")?;
    writeln!(writer, "#         max: 100")?;
    writeln!(writer, "#         suffix: '%'")?;
    writeln!(writer, "#     antisense_pct:")?;
    writeln!(writer, "#         title: 'Antisense %'")?;
    writeln!(writer, "#         min: 0")?;
    writeln!(writer, "#         max: 100")?;
    writeln!(writer, "#         suffix: '%'")?;

    // TSV data
    writeln!(
        writer,
        "Sample\tsense_pct\tantisense_pct\tundetermined_pct\tprotocol\tdata_type"
    )?;
    writeln!(
        writer,
        "{}\t{:.2}\t{:.2}\t{:.2}\t{}\t{}",
        sample_name,
        sense_frac * 100.0,
        antisense_frac * 100.0,
        undetermined_frac * 100.0,
        protocol,
        if is_paired { "PairEnd" } else { "SingleEnd" }
    )?;

    Ok(())
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_infer_forward_stranded() {
        let (protocol, sense, antisense, _) = infer_strandedness(900, 100, 0);
        assert_eq!(protocol, StrandProtocol::ForwardStranded);
        assert!((sense - 0.9).abs() < 0.001);
        assert!((antisense - 0.1).abs() < 0.001);
    }

    #[test]
    fn test_infer_reverse_stranded() {
        let (protocol, sense, antisense, _) = infer_strandedness(50, 950, 0);
        assert_eq!(protocol, StrandProtocol::ReverseStranded);
        assert!((sense - 0.05).abs() < 0.001);
        assert!((antisense - 0.95).abs() < 0.001);
    }

    #[test]
    fn test_infer_unstranded() {
        let (protocol, sense, antisense, _) = infer_strandedness(480, 520, 0);
        assert_eq!(protocol, StrandProtocol::Unstranded);
        assert!((sense - 0.48).abs() < 0.001);
        assert!((antisense - 0.52).abs() < 0.001);
    }

    #[test]
    fn test_infer_undetermined_insufficient_reads() {
        let (protocol, _, _, _) = infer_strandedness(10, 10, 5);
        assert_eq!(protocol, StrandProtocol::Undetermined);
    }

    #[test]
    fn test_infer_undetermined_ambiguous() {
        // 65/35 split: not stranded enough, not unstranded enough
        let (protocol, _, _, _) = infer_strandedness(650, 350, 0);
        assert_eq!(protocol, StrandProtocol::Undetermined);
    }

    #[test]
    fn test_infer_no_informative_reads() {
        let (protocol, _, _, undetermined) = infer_strandedness(0, 0, 500);
        assert_eq!(protocol, StrandProtocol::Undetermined);
        assert!((undetermined - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_strand_protocol_display() {
        assert_eq!(
            format!("{}", StrandProtocol::ForwardStranded),
            "forward (fr-secondstrand)"
        );
        assert_eq!(
            format!("{}", StrandProtocol::ReverseStranded),
            "reverse (fr-firststrand)"
        );
        assert_eq!(format!("{}", StrandProtocol::Unstranded), "unstranded");
        assert_eq!(format!("{}", StrandProtocol::Undetermined), "undetermined");
    }

    #[test]
    fn test_write_infer_experiment_paired() {
        let outdir = std::path::PathBuf::from("tests/output_strandedness_pe");
        let _ = std::fs::create_dir_all(&outdir);
        let path = outdir.join("test.infer_experiment.txt");

        let result = CountResult {
            gene_counts: indexmap::IndexMap::new(),
            total_reads_multi_dup: 0,
            total_reads_multi_nodup: 0,
            total_reads_unique_dup: 0,
            total_reads_unique_nodup: 0,
            stat_total_reads: 1000,
            stat_assigned: 800,
            stat_ambiguous: 50,
            stat_no_features: 150,
            stat_total_fragments: 500,
            stat_total_dup: 100,
            stat_total_multi: 50,
            strandedness_sense: 50,
            strandedness_antisense: 900,
            strandedness_undetermined: 50,
        };

        write_infer_experiment(&path, &result, true).unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("PairEnd Data"));
        assert!(content.contains("1++,1--,2+-,2-+"));
        assert!(content.contains("1+-,1-+,2++,2--"));
        assert!(content.contains("reverse (fr-firststrand)"));

        let _ = std::fs::remove_dir_all(&outdir);
    }

    #[test]
    fn test_write_infer_experiment_single() {
        let outdir = std::path::PathBuf::from("tests/output_strandedness_se");
        let _ = std::fs::create_dir_all(&outdir);
        let path = outdir.join("test.infer_experiment.txt");

        let result = CountResult {
            gene_counts: indexmap::IndexMap::new(),
            total_reads_multi_dup: 0,
            total_reads_multi_nodup: 0,
            total_reads_unique_dup: 0,
            total_reads_unique_nodup: 0,
            stat_total_reads: 1000,
            stat_assigned: 800,
            stat_ambiguous: 50,
            stat_no_features: 150,
            stat_total_fragments: 500,
            stat_total_dup: 100,
            stat_total_multi: 50,
            strandedness_sense: 850,
            strandedness_antisense: 100,
            strandedness_undetermined: 50,
        };

        write_infer_experiment(&path, &result, false).unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("SingleEnd Data"));
        assert!(content.contains("\"++,--\""));
        assert!(content.contains("\"+-,-+\""));
        assert!(content.contains("forward (fr-secondstrand)"));

        let _ = std::fs::remove_dir_all(&outdir);
    }

    #[test]
    fn test_write_strandedness_mqc() {
        let outdir = std::path::PathBuf::from("tests/output_strandedness_mqc");
        let _ = std::fs::create_dir_all(&outdir);
        let path = outdir.join("test.strandedness_mqc.txt");

        let result = CountResult {
            gene_counts: indexmap::IndexMap::new(),
            total_reads_multi_dup: 0,
            total_reads_multi_nodup: 0,
            total_reads_unique_dup: 0,
            total_reads_unique_nodup: 0,
            stat_total_reads: 1000,
            stat_assigned: 800,
            stat_ambiguous: 50,
            stat_no_features: 150,
            stat_total_fragments: 500,
            stat_total_dup: 100,
            stat_total_multi: 50,
            strandedness_sense: 100,
            strandedness_antisense: 100,
            strandedness_undetermined: 0,
        };

        write_strandedness_mqc(&path, &result, "sample1", true).unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("# id: 'strandedness'"));
        assert!(content.contains("sample1"));
        assert!(content.contains("unstranded"));
        assert!(content.contains("PairEnd"));

        let _ = std::fs::remove_dir_all(&outdir);
    }
}
