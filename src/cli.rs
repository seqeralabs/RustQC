//! Command-line interface definition for RustQC.
//!
//! Provides a subcommand-based CLI. The `rna` subcommand runs all RNA-Seq QC
//! analyses in a single pass: dupRadar duplication rate analysis, featureCounts-
//! compatible output, RSeQC-equivalent metrics (bam_stat, infer_experiment,
//! read_duplication, read_distribution, junction_annotation, junction_saturation,
//! inner_distance), TIN (Transcript Integrity Number), preseq library complexity
//! extrapolation, samtools-compatible outputs (flagstat, idxstats, stats), and
//! Qualimap RNA-seq QC. Individual tools can be disabled via the YAML config file.
//!
//! A GTF gene annotation file is required for all analyses.

use clap::{Parser, Subcommand, ValueEnum};
use serde::Deserialize;

/// Library strandedness protocol.
///
/// Determines how read strand is interpreted relative to the gene annotation
/// strand during counting. Accepted CLI values: `unstranded`, `forward`, `reverse`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, ValueEnum, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Strandedness {
    /// Count reads on either strand (library is not strand-specific).
    #[default]
    Unstranded,
    /// Forward stranded: read 1 maps to the transcript strand.
    Forward,
    /// Reverse stranded: read 2 maps to the transcript strand (e.g. dUTP).
    Reverse,
}

impl std::fmt::Display for Strandedness {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strandedness::Unstranded => write!(f, "unstranded"),
            Strandedness::Forward => write!(f, "forward"),
            Strandedness::Reverse => write!(f, "reverse"),
        }
    }
}

/// Fast quality control tools for sequencing data, written in Rust.
#[derive(Parser, Debug)]
#[command(name = "rustqc", version, about, long_about = None)]
pub struct Cli {
    /// The analysis subcommand to run.
    #[command(subcommand)]
    pub command: Commands,
}

/// Available analysis subcommands.
#[derive(Subcommand, Debug)]
pub enum Commands {
    /// RNA-Seq QC — single-pass analysis of BAM/SAM/CRAM files.
    ///
    /// Runs featureCounts, dupRadar, Qualimap, samtools stats, and RSeQC
    /// analyses in one pass. Requires a GTF annotation and duplicate-marked
    /// (not removed) alignments.
    Rna(RnaArgs),
}

/// Arguments for the `rna` subcommand.
#[derive(Parser, Debug)]
#[command(
    next_line_help = false,
    term_width = 120,
    help_template = "\
{about-with-newline}
{usage-heading} {usage}

{all-args}"
)]
pub struct RnaArgs {
    // ── Input / Output ──────────────────────────────────────────────────
    /// Duplicate-marked alignment file(s)
    #[arg(value_name = "INPUT", num_args = 1.., required = true, help_heading = "Input / Output")]
    pub input: Vec<String>,

    /// GTF gene annotation (plain or .gz)
    #[arg(short, long, value_name = "GTF", help_heading = "Input / Output")]
    pub gtf: String,

    /// Reference FASTA (required for CRAM)
    #[arg(short, long, value_name = "FASTA", help_heading = "Input / Output")]
    pub reference: Option<String>,

    /// Output directory [default: .]
    #[arg(
        short,
        long,
        default_value = ".",
        hide_default_value = true,
        help_heading = "Input / Output"
    )]
    pub outdir: String,

    /// Write outputs to a flat directory (no subdirs)
    #[arg(long, default_value_t = false, help_heading = "Input / Output")]
    pub flat_output: bool,

    /// YAML configuration file
    #[arg(short, long, value_name = "CONFIG", help_heading = "Input / Output")]
    pub config: Option<String>,

    /// JSON summary path (use "-" for stdout)
    #[arg(short = 'j', long = "json-summary", value_name = "PATH", num_args = 0..=1, default_missing_value = "", help_heading = "Input / Output")]
    pub json_summary: Option<String>,

    // ── Library ─────────────────────────────────────────────────────────
    /// Strandedness: unstranded, forward, reverse
    #[arg(short, long, value_enum, help_heading = "Library")]
    pub stranded: Option<Strandedness>,

    /// Paired-end reads
    #[arg(short, long, help_heading = "Library")]
    pub paired: bool,

    // ── General ─────────────────────────────────────────────────────────
    /// Number of threads [default: 1]
    #[arg(
        short,
        long,
        default_value_t = 1,
        hide_default_value = true,
        help_heading = "General"
    )]
    pub threads: usize,

    /// MAPQ cutoff for quality filtering [default: 30]
    #[arg(
        short = 'Q',
        long = "mapq",
        default_value_t = 30,
        hide_default_value = true,
        help_heading = "General"
    )]
    pub mapq_cut: u8,

    /// GTF attribute for biotype grouping
    #[arg(long, value_name = "ATTR", help_heading = "General")]
    pub biotype_attribute: Option<String>,

    /// Skip duplicate-marking header check
    #[arg(long, default_value_t = false, help_heading = "General")]
    pub skip_dup_check: bool,

    /// Suppress output except warnings/errors
    #[arg(
        short = 'q',
        long,
        conflicts_with = "verbose",
        help_heading = "General"
    )]
    pub quiet: bool,

    /// Show additional detail
    #[arg(short = 'v', long, conflicts_with = "quiet", help_heading = "General")]
    pub verbose: bool,

    // ── Tool parameters ─────────────────────────────────────────────────
    /// infer_experiment: sample size [default: 200000]
    #[arg(
        long = "infer-experiment-sample-size",
        value_name = "N",
        help_heading = "Tool parameters"
    )]
    pub infer_experiment_sample_size: Option<u64>,

    /// junction_annotation: min intron size [default: 50]
    #[arg(
        long = "min-intron",
        value_name = "N",
        help_heading = "Tool parameters"
    )]
    pub min_intron: Option<u64>,

    /// junction_saturation: random seed for reproducible results
    #[arg(
        long = "junction-saturation-seed",
        value_name = "N",
        help_heading = "Tool parameters"
    )]
    pub junction_saturation_seed: Option<u64>,

    /// junction_saturation: min coverage [default: 1]
    #[arg(
        long = "junction-saturation-min-coverage",
        value_name = "N",
        help_heading = "Tool parameters"
    )]
    pub junction_saturation_min_coverage: Option<u64>,

    /// junction_saturation: start % [default: 5]
    #[arg(
        long = "junction-saturation-percentile-floor",
        value_name = "N",
        help_heading = "Tool parameters"
    )]
    pub junction_saturation_percentile_floor: Option<u64>,

    /// junction_saturation: end % [default: 100]
    #[arg(
        long = "junction-saturation-percentile-ceiling",
        value_name = "N",
        help_heading = "Tool parameters"
    )]
    pub junction_saturation_percentile_ceiling: Option<u64>,

    /// junction_saturation: step % [default: 5]
    #[arg(
        long = "junction-saturation-percentile-step",
        value_name = "N",
        help_heading = "Tool parameters"
    )]
    pub junction_saturation_percentile_step: Option<u64>,

    /// inner_distance: sample size [default: 1000000]
    #[arg(
        long = "inner-distance-sample-size",
        value_name = "N",
        help_heading = "Tool parameters"
    )]
    pub inner_distance_sample_size: Option<u64>,

    /// inner_distance: lower bound [default: -250]
    #[arg(
        long = "inner-distance-lower-bound",
        value_name = "N",
        allow_hyphen_values = true,
        help_heading = "Tool parameters"
    )]
    pub inner_distance_lower_bound: Option<i64>,

    /// inner_distance: upper bound [default: 250]
    #[arg(
        long = "inner-distance-upper-bound",
        value_name = "N",
        allow_hyphen_values = true,
        help_heading = "Tool parameters"
    )]
    pub inner_distance_upper_bound: Option<i64>,

    /// inner_distance: bin width [default: 5]
    #[arg(
        long = "inner-distance-step",
        value_name = "N",
        allow_hyphen_values = true,
        help_heading = "Tool parameters"
    )]
    pub inner_distance_step: Option<i64>,

    /// TIN: random seed for reproducible results
    #[arg(long = "tin-seed", value_name = "N", help_heading = "Tool parameters")]
    pub tin_seed: Option<u64>,

    /// Skip TIN analysis
    #[arg(long, default_value_t = false, help_heading = "Tool parameters")]
    pub skip_tin: bool,

    /// Skip read duplication analysis
    #[arg(long, default_value_t = false, help_heading = "Tool parameters")]
    pub skip_read_duplication: bool,

    /// Skip preseq library complexity analysis
    #[arg(long, default_value_t = false, help_heading = "Tool parameters")]
    pub skip_preseq: bool,

    /// preseq: random seed for bootstrap CIs
    #[arg(
        long = "preseq-seed",
        value_name = "N",
        help_heading = "Tool parameters"
    )]
    pub preseq_seed: Option<u64>,

    /// preseq: max extrapolation depth
    #[arg(
        long = "preseq-max-extrap",
        value_name = "N",
        help_heading = "Tool parameters"
    )]
    pub preseq_max_extrap: Option<f64>,

    /// preseq: step size between points
    #[arg(
        long = "preseq-step-size",
        value_name = "N",
        help_heading = "Tool parameters"
    )]
    pub preseq_step_size: Option<f64>,

    /// preseq: bootstrap replicates for CIs
    #[arg(
        long = "preseq-n-bootstraps",
        value_name = "N",
        help_heading = "Tool parameters"
    )]
    pub preseq_n_bootstraps: Option<u32>,

    /// preseq: max segment length for PE merging
    #[arg(
        long = "preseq-seg-len",
        value_name = "N",
        help_heading = "Tool parameters"
    )]
    pub preseq_seg_len: Option<i64>,
}

/// Parse command-line arguments and return the Cli struct.
pub fn parse_args() -> Cli {
    Cli::parse()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rna_default_args_gtf() {
        // Test that defaults are sensible with a GTF annotation
        let cli = Cli::parse_from(["rustqc", "rna", "test.bam", "--gtf", "genes.gtf"]);
        match cli.command {
            Commands::Rna(args) => {
                assert_eq!(args.input, vec!["test.bam"]);
                assert_eq!(args.gtf, "genes.gtf");
                assert_eq!(args.stranded, None);
                assert!(!args.paired);
                assert_eq!(args.threads, 1);
                assert_eq!(args.outdir, ".");
                assert!(args.biotype_attribute.is_none());
                assert_eq!(args.mapq_cut, 30);
                assert_eq!(args.infer_experiment_sample_size, None);
                assert_eq!(args.min_intron, None);
                assert_eq!(args.inner_distance_step, None);
            }
            #[allow(unreachable_patterns)]
            _ => panic!("Expected Rna subcommand"),
        }
    }

    #[test]
    fn test_rna_multiple_bams() {
        // Test that multiple BAM files are accepted
        let cli = Cli::parse_from([
            "rustqc",
            "rna",
            "a.bam",
            "b.bam",
            "c.bam",
            "--gtf",
            "genes.gtf",
        ]);
        match cli.command {
            Commands::Rna(args) => {
                assert_eq!(args.input, vec!["a.bam", "b.bam", "c.bam"]);
                assert_eq!(args.gtf, "genes.gtf");
            }
            #[allow(unreachable_patterns)]
            _ => panic!("Expected Rna subcommand"),
        }
    }

    #[test]
    fn test_rna_gtf_all_args() {
        let cli = Cli::parse_from([
            "rustqc",
            "rna",
            "test.bam",
            "--gtf",
            "genes.gtf",
            "--stranded",
            "reverse",
            "--paired",
            "--threads",
            "4",
            "--outdir",
            "/tmp/out",
            "--reference",
            "genome.fa",
            "-Q",
            "20",
        ]);
        match cli.command {
            Commands::Rna(args) => {
                assert_eq!(args.gtf, "genes.gtf");
                assert_eq!(args.stranded, Some(Strandedness::Reverse));
                assert!(args.paired);
                assert_eq!(args.threads, 4);
                assert_eq!(args.outdir, "/tmp/out");
                assert_eq!(args.reference, Some("genome.fa".to_string()));
                assert_eq!(args.mapq_cut, 20);
            }
            #[allow(unreachable_patterns)]
            _ => panic!("Expected Rna subcommand"),
        }
    }

    #[test]
    fn test_rna_missing_gtf() {
        // --gtf is required, so omitting it should fail
        let result = Cli::try_parse_from(["rustqc", "rna", "test.bam"]);
        assert!(result.is_err(), "Expected error when --gtf is not provided");
    }

    #[test]
    fn test_rna_rseqc_params() {
        let cli = Cli::parse_from([
            "rustqc",
            "rna",
            "test.bam",
            "--gtf",
            "genes.gtf",
            "--infer-experiment-sample-size",
            "500000",
            "--min-intron",
            "100",
            "--junction-saturation-min-coverage",
            "5",
            "--inner-distance-lower-bound",
            "-500",
            "--inner-distance-upper-bound",
            "500",
            "--inner-distance-step",
            "10",
        ]);
        match cli.command {
            Commands::Rna(args) => {
                assert_eq!(args.infer_experiment_sample_size, Some(500_000));
                assert_eq!(args.min_intron, Some(100));
                assert_eq!(args.junction_saturation_min_coverage, Some(5));
                assert_eq!(args.inner_distance_lower_bound, Some(-500));
                assert_eq!(args.inner_distance_upper_bound, Some(500));
                assert_eq!(args.inner_distance_step, Some(10));
            }
            #[allow(unreachable_patterns)]
            _ => panic!("Expected Rna subcommand"),
        }
    }

    #[test]
    fn test_rna_preseq_params() {
        let cli = Cli::parse_from([
            "rustqc",
            "rna",
            "test.bam",
            "--gtf",
            "genes.gtf",
            "--preseq-max-extrap",
            "5000000000",
            "--preseq-step-size",
            "500000",
            "--preseq-n-bootstraps",
            "200",
            "--preseq-seg-len",
            "100000000",
        ]);
        match cli.command {
            Commands::Rna(args) => {
                assert!(!args.skip_preseq);
                assert_eq!(args.preseq_max_extrap, Some(5_000_000_000.0));
                assert_eq!(args.preseq_step_size, Some(500_000.0));
                assert_eq!(args.preseq_n_bootstraps, Some(200));
                assert_eq!(args.preseq_seg_len, Some(100_000_000));
            }
            #[allow(unreachable_patterns)]
            _ => panic!("Expected Rna subcommand"),
        }
    }

    #[test]
    fn test_rna_tool_seeds() {
        let cli = Cli::parse_from([
            "rustqc",
            "rna",
            "test.bam",
            "--gtf",
            "genes.gtf",
            "--preseq-seed",
            "1",
            "--tin-seed",
            "2",
            "--junction-saturation-seed",
            "3",
        ]);
        match cli.command {
            Commands::Rna(args) => {
                assert_eq!(args.preseq_seed, Some(1));
                assert_eq!(args.tin_seed, Some(2));
                assert_eq!(args.junction_saturation_seed, Some(3));
            }
            #[allow(unreachable_patterns)]
            _ => panic!("Expected Rna subcommand"),
        }
    }

    #[test]
    fn test_rna_skip_preseq() {
        let cli = Cli::parse_from([
            "rustqc",
            "rna",
            "test.bam",
            "--gtf",
            "genes.gtf",
            "--skip-preseq",
        ]);
        match cli.command {
            Commands::Rna(args) => {
                assert!(args.skip_preseq);
            }
            #[allow(unreachable_patterns)]
            _ => panic!("Expected Rna subcommand"),
        }
    }
}
