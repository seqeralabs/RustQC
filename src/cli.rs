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

use rustqc::Strandedness;
use usage::{Args, Cli, Subcommands};

/// Fast quality control tools for sequencing data, written in Rust.
#[derive(Cli, Debug)]
#[usage(
    bin = "rustqc",
    version,
    long_version = rustqc_long_version(),
    long_version_spec = "0.2.1",
    author = "Phil Ewels <phil.ewels@seqera.io>, Jon Manning <jonathan.manning@seqera.io>",
    license = "GPL-3.0-or-later",
    repository = "https://github.com/seqeralabs/RustQC",
    unknown_flags = "error",
    args_override_self = false,
    completion,
    term_width = 120
)]
pub struct Cli {
    /// The analysis subcommand to run.
    #[usage(subcommand)]
    pub command: Commands,
}

/// Available analysis subcommands.
#[derive(Subcommands, Debug)]
pub enum Commands {
    /// RNA-Seq QC — single-pass analysis of BAM/SAM/CRAM files.
    ///
    /// Runs featureCounts, dupRadar, Qualimap, samtools stats, and RSeQC
    /// analyses in one pass. Requires a GTF annotation and duplicate-marked
    /// (not removed) alignments.
    Rna(RnaArgs),
}

/// Arguments for the `rna` subcommand.
#[derive(Args, Debug)]
#[usage(unknown_flags = "error", args_override_self = false, term_width = 120)]
pub struct RnaArgs {
    // ── Input / Output ──────────────────────────────────────────────────
    /// Duplicate-marked alignment file(s)
    #[usage(arg, required, value_name = "INPUT", help_heading = "Input / Output")]
    pub input: Vec<String>,

    /// GTF gene annotation (plain or .gz)
    #[usage(
        short = 'g',
        long,
        value_name = "GTF",
        env = "RUSTQC_GTF",
        help_heading = "Input / Output"
    )]
    pub gtf: String,

    /// Reference FASTA (required for CRAM)
    #[usage(
        short = 'r',
        long,
        value_name = "FASTA",
        env = "RUSTQC_REFERENCE",
        help_heading = "Input / Output"
    )]
    pub reference: Option<String>,

    /// Output directory [default: .]
    #[usage(
        short = 'o',
        long,
        default = ".",
        hide_default_value,
        env = "RUSTQC_OUTDIR",
        help_heading = "Input / Output"
    )]
    pub outdir: String,

    /// Override sample name for output filenames (default: derived from BAM filename)
    #[usage(
        long,
        value_name = "NAME",
        env = "RUSTQC_SAMPLE_NAME",
        help_heading = "Input / Output"
    )]
    pub sample_name: Option<String>,

    /// Write outputs to a flat directory (no subdirs)
    #[usage(long, env = "RUSTQC_FLAT_OUTPUT", help_heading = "Input / Output")]
    pub flat_output: bool,

    /// YAML configuration file (see also: RUSTQC_CONFIG env var)
    #[usage(
        short = 'c',
        long,
        value_name = "CONFIG",
        help_heading = "Input / Output"
    )]
    pub config: Option<String>,

    /// JSON summary path (use "-" for stdout)
    #[usage(
        short = 'j',
        long = "json-summary",
        value_name = "PATH",
        default_missing = "",
        allow_hyphen_values,
        env = "RUSTQC_JSON_SUMMARY",
        help_heading = "Input / Output"
    )]
    pub json_summary: Option<String>,

    // ── Library ─────────────────────────────────────────────────────────
    /// Strandedness: unstranded, forward, reverse
    #[usage(
        short = 's',
        long,
        value_enum,
        env = "RUSTQC_STRANDED",
        help_heading = "Library"
    )]
    pub stranded: Option<Strandedness>,

    /// Paired-end reads
    #[usage(short = 'p', long, env = "RUSTQC_PAIRED", help_heading = "Library")]
    pub paired: bool,

    // ── General ─────────────────────────────────────────────────────────
    /// Number of threads [default: 1]
    #[usage(
        short = 't',
        long,
        default = "1",
        hide_default_value,
        env = "RUSTQC_THREADS",
        help_heading = "General"
    )]
    pub threads: usize,

    /// MAPQ cutoff for quality filtering [default: 30]
    #[usage(
        short = 'Q',
        long = "mapq",
        default = "30",
        hide_default_value,
        env = "RUSTQC_MAPQ",
        help_heading = "General"
    )]
    pub mapq_cut: u8,

    /// GTF attribute for biotype grouping
    #[usage(
        long,
        value_name = "ATTR",
        env = "RUSTQC_BIOTYPE_ATTRIBUTE",
        help_heading = "General"
    )]
    pub biotype_attribute: Option<String>,

    /// Skip duplicate-marking check
    #[usage(long, env = "RUSTQC_SKIP_DUP_CHECK", help_heading = "General")]
    pub skip_dup_check: bool,

    /// Suppress output except warnings/errors
    #[usage(
        short = 'q',
        long,
        conflicts("--verbose"),
        env = "RUSTQC_QUIET",
        help_heading = "General"
    )]
    pub quiet: bool,

    /// Show additional detail
    #[usage(
        short = 'v',
        long,
        conflicts("--quiet"),
        env = "RUSTQC_VERBOSE",
        help_heading = "General"
    )]
    pub verbose: bool,

    // ── Tool parameters ─────────────────────────────────────────────────
    /// infer_experiment: sample size [default: 200000]
    #[usage(
        long = "infer-experiment-sample-size",
        value_name = "N",
        env = "RUSTQC_INFER_EXPERIMENT_SAMPLE_SIZE",
        help_heading = "Tool parameters"
    )]
    pub infer_experiment_sample_size: Option<u64>,

    /// junction_annotation: min intron size [default: 50]
    #[usage(
        long = "min-intron",
        value_name = "N",
        env = "RUSTQC_MIN_INTRON",
        help_heading = "Tool parameters"
    )]
    pub min_intron: Option<u64>,

    /// junction_saturation: random seed for reproducible results
    #[usage(
        long = "junction-saturation-seed",
        value_name = "N",
        env = "RUSTQC_JUNCTION_SATURATION_SEED",
        help_heading = "Tool parameters"
    )]
    pub junction_saturation_seed: Option<u64>,

    /// junction_saturation: min coverage [default: 1]
    #[usage(
        long = "junction-saturation-min-coverage",
        value_name = "N",
        env = "RUSTQC_JUNCTION_SATURATION_MIN_COVERAGE",
        help_heading = "Tool parameters"
    )]
    pub junction_saturation_min_coverage: Option<u64>,

    /// junction_saturation: start % [default: 5]
    #[usage(
        long = "junction-saturation-percentile-floor",
        value_name = "N",
        env = "RUSTQC_JUNCTION_SATURATION_PERCENTILE_FLOOR",
        help_heading = "Tool parameters"
    )]
    pub junction_saturation_percentile_floor: Option<u64>,

    /// junction_saturation: end % [default: 100]
    #[usage(
        long = "junction-saturation-percentile-ceiling",
        value_name = "N",
        env = "RUSTQC_JUNCTION_SATURATION_PERCENTILE_CEILING",
        help_heading = "Tool parameters"
    )]
    pub junction_saturation_percentile_ceiling: Option<u64>,

    /// junction_saturation: step % [default: 5]
    #[usage(
        long = "junction-saturation-percentile-step",
        value_name = "N",
        env = "RUSTQC_JUNCTION_SATURATION_PERCENTILE_STEP",
        help_heading = "Tool parameters"
    )]
    pub junction_saturation_percentile_step: Option<u64>,

    /// inner_distance: sample size [default: 1000000]
    #[usage(
        long = "inner-distance-sample-size",
        value_name = "N",
        env = "RUSTQC_INNER_DISTANCE_SAMPLE_SIZE",
        help_heading = "Tool parameters"
    )]
    pub inner_distance_sample_size: Option<u64>,

    /// inner_distance: lower bound [default: -250]
    #[usage(
        long = "inner-distance-lower-bound",
        value_name = "N",
        allow_negative_numbers,
        env = "RUSTQC_INNER_DISTANCE_LOWER_BOUND",
        help_heading = "Tool parameters"
    )]
    pub inner_distance_lower_bound: Option<i64>,

    /// inner_distance: upper bound [default: 250]
    #[usage(
        long = "inner-distance-upper-bound",
        value_name = "N",
        allow_negative_numbers,
        env = "RUSTQC_INNER_DISTANCE_UPPER_BOUND",
        help_heading = "Tool parameters"
    )]
    pub inner_distance_upper_bound: Option<i64>,

    /// inner_distance: bin width [default: 5]
    #[usage(
        long = "inner-distance-step",
        value_name = "N",
        allow_negative_numbers,
        env = "RUSTQC_INNER_DISTANCE_STEP",
        help_heading = "Tool parameters"
    )]
    pub inner_distance_step: Option<i64>,

    /// TIN: random seed for reproducible results
    #[usage(
        long = "tin-seed",
        value_name = "N",
        env = "RUSTQC_TIN_SEED",
        help_heading = "Tool parameters"
    )]
    pub tin_seed: Option<u64>,

    /// Skip TIN analysis
    #[usage(long, env = "RUSTQC_SKIP_TIN", help_heading = "Tool parameters")]
    pub skip_tin: bool,

    /// Skip read duplication analysis
    #[usage(
        long,
        env = "RUSTQC_SKIP_READ_DUPLICATION",
        help_heading = "Tool parameters"
    )]
    pub skip_read_duplication: bool,

    /// Skip preseq library complexity analysis
    #[usage(long, env = "RUSTQC_SKIP_PRESEQ", help_heading = "Tool parameters")]
    pub skip_preseq: bool,

    /// preseq: random seed for bootstrap CIs
    #[usage(
        long = "preseq-seed",
        value_name = "N",
        env = "RUSTQC_PRESEQ_SEED",
        help_heading = "Tool parameters"
    )]
    pub preseq_seed: Option<u64>,

    /// preseq: max extrapolation depth
    #[usage(
        long = "preseq-max-extrap",
        value_name = "N",
        env = "RUSTQC_PRESEQ_MAX_EXTRAP",
        help_heading = "Tool parameters"
    )]
    pub preseq_max_extrap: Option<f64>,

    /// preseq: step size between points
    #[usage(
        long = "preseq-step-size",
        value_name = "N",
        env = "RUSTQC_PRESEQ_STEP_SIZE",
        help_heading = "Tool parameters"
    )]
    pub preseq_step_size: Option<f64>,

    /// preseq: bootstrap replicates for CIs
    #[usage(
        long = "preseq-n-bootstraps",
        value_name = "N",
        env = "RUSTQC_PRESEQ_N_BOOTSTRAPS",
        help_heading = "Tool parameters"
    )]
    pub preseq_n_bootstraps: Option<u32>,

    /// preseq: max segment length for PE merging
    #[usage(
        long = "preseq-seg-len",
        value_name = "N",
        env = "RUSTQC_PRESEQ_SEG_LEN",
        help_heading = "Tool parameters"
    )]
    pub preseq_seg_len: Option<i64>,
}

/// Extended `--version` text: crate version, git commit, build time, and CPU info.
///
/// Leaked so the usage derive can treat it as `&'static str` next to
/// `env!("CARGO_PKG_VERSION")` in the generated `-V` / `--version` arms.
fn rustqc_long_version() -> &'static str {
    Box::leak(
        format!(
            "{} ({}, built {})\n{}",
            env!("CARGO_PKG_VERSION"),
            env!("GIT_SHORT_HASH"),
            env!("BUILD_TIMESTAMP"),
            rustqc::cpu::cpu_info_line(),
        )
        .into_boxed_str(),
    )
}

/// Parse process argv, printing help/version and exiting the way `usage` does.
///
/// `--version` includes the git commit, build timestamp, and CPU info line.
pub fn parse_args() -> Cli {
    Cli::parse()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ffi::OsStr;

    /// Convert clap-shaped string argv (including argv0) into the `&OsStr` slice usage expects.
    fn os_args<'a>(args: &'a [&'a str]) -> Vec<&'a OsStr> {
        args.iter().copied().map(OsStr::new).collect()
    }

    fn parse_cli(args: &[&str]) -> Cli {
        Cli::try_parse_from(&os_args(args))
            .unwrap_or_else(|err| panic!("expected a successful parse for {args:?}, got {err:?}"))
    }

    fn parse_cli_err(args: &[&str]) -> bool {
        Cli::try_parse_from(&os_args(args)).is_err()
    }

    #[test]
    fn test_rna_default_args_gtf() {
        // Test that defaults are sensible with a GTF annotation
        let cli = parse_cli(&["rustqc", "rna", "test.bam", "--gtf", "genes.gtf"]);
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
        let cli = parse_cli(&[
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
        let cli = parse_cli(&[
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
        assert!(
            parse_cli_err(&["rustqc", "rna", "test.bam"]),
            "Expected error when --gtf is not provided"
        );
    }

    #[test]
    fn test_rna_rseqc_params() {
        let cli = parse_cli(&[
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
        let cli = parse_cli(&[
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
        let cli = parse_cli(&[
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
        let cli = parse_cli(&[
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

    #[test]
    fn test_rna_unknown_flag_is_rejected() {
        assert!(
            parse_cli_err(&["rustqc", "rna", "test.bam", "--gtf", "genes.gtf", "--wat"]),
            "Expected error for an unknown flag"
        );
    }

    #[test]
    fn test_rna_json_summary_bare_flag() {
        let cli = parse_cli(&[
            "rustqc",
            "rna",
            "test.bam",
            "--gtf",
            "genes.gtf",
            "--json-summary",
        ]);
        match cli.command {
            Commands::Rna(args) => {
                assert_eq!(args.json_summary.as_deref(), Some(""));
            }
            #[allow(unreachable_patterns)]
            _ => panic!("Expected Rna subcommand"),
        }
    }

    #[test]
    fn test_rna_json_summary_stdout() {
        let cli = parse_cli(&[
            "rustqc",
            "rna",
            "test.bam",
            "--gtf",
            "genes.gtf",
            "--json-summary",
            "-",
        ]);
        match cli.command {
            Commands::Rna(args) => {
                assert_eq!(args.json_summary.as_deref(), Some("-"));
            }
            #[allow(unreachable_patterns)]
            _ => panic!("Expected Rna subcommand"),
        }
    }
}
