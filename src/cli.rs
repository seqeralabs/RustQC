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

use clap::{CommandFactory, Parser, Subcommand};

use rustqc::Strandedness;

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

    /// DNA QC — single-pass analysis of BAM/SAM/CRAM files.
    ///
    /// Runs depth of coverage, samtools stats and library complexity
    /// estimation in one pass. Needs no gene annotation. Pass `--targets`
    /// to switch to targeted (exome or panel) mode.
    Dna(DnaArgs),

    /// Protein QC — sequence, coding-region and mass spectrometry analyses.
    ///
    /// Three modes taking different inputs entirely, so the mode is chosen
    /// explicitly rather than inferred from which flags were given.
    Protein(ProteinArgs),
}

/// Arguments for the `protein` subcommand.
#[derive(Parser, Debug)]
pub struct ProteinArgs {
    /// Which protein analysis to run.
    #[command(subcommand)]
    pub mode: ProteinMode,
}

/// The `protein` subcommand's modes.
#[derive(Subcommand, Debug)]
pub enum ProteinMode {
    /// Protein FASTA QC: length statistics, composition and defects.
    Sequence(ProteinSequenceArgs),
}

/// Arguments for `protein sequence`.
#[derive(Parser, Debug)]
#[command(
    next_line_help = false,
    term_width = 120,
    help_template = "\
{about-with-newline}
{usage-heading} {usage}

{all-args}"
)]
pub struct ProteinSequenceArgs {
    /// Protein FASTA file(s), plain or .gz
    #[arg(value_name = "FASTA", num_args = 1.., required = true, help_heading = "Input / Output")]
    pub input: Vec<String>,

    /// Output directory [default: .]
    #[arg(
        short,
        long,
        default_value = ".",
        hide_default_value = true,
        env = "RUSTQC_OUTDIR",
        help_heading = "Input / Output"
    )]
    pub outdir: String,

    /// Override sample name for output filenames (default: derived from filename)
    #[arg(
        long,
        value_name = "NAME",
        env = "RUSTQC_SAMPLE_NAME",
        help_heading = "Input / Output"
    )]
    pub sample_name: Option<String>,

    /// Write outputs to a flat directory (no subdirs)
    #[arg(
        long,
        default_value_t = false,
        env = "RUSTQC_FLAT_OUTPUT",
        help_heading = "Input / Output"
    )]
    pub flat_output: bool,

    /// YAML configuration file (see also: RUSTQC_CONFIG env var)
    #[arg(short, long, value_name = "CONFIG", help_heading = "Input / Output")]
    pub config: Option<String>,

    /// JSON summary path (use "-" for stdout)
    #[arg(short = 'j', long = "json-summary", value_name = "PATH", num_args = 0..=1, default_missing_value = "", env = "RUSTQC_JSON_SUMMARY", help_heading = "Input / Output")]
    pub json_summary: Option<String>,

    /// Ignore sequences shorter than this
    #[arg(
        long = "min-length",
        value_name = "N",
        default_value_t = 0,
        hide_default_value = true,
        env = "RUSTQC_MIN_LENGTH",
        help_heading = "Tool parameters"
    )]
    pub min_length: usize,

    /// Report a missing terminal stop codon as a defect
    #[arg(
        long = "expect-stop",
        default_value_t = false,
        env = "RUSTQC_EXPECT_STOP",
        help_heading = "Tool parameters"
    )]
    pub expect_stop: bool,

    /// Suppress output except warnings/errors
    #[arg(
        short = 'q',
        long,
        conflicts_with = "verbose",
        env = "RUSTQC_QUIET",
        help_heading = "General"
    )]
    pub quiet: bool,

    /// Show additional detail
    #[arg(
        short = 'v',
        long,
        conflicts_with = "quiet",
        env = "RUSTQC_VERBOSE",
        help_heading = "General"
    )]
    pub verbose: bool,
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
    #[arg(
        short,
        long,
        value_name = "GTF",
        env = "RUSTQC_GTF",
        help_heading = "Input / Output"
    )]
    pub gtf: String,

    /// Reference FASTA (required for CRAM)
    #[arg(
        short,
        long,
        value_name = "FASTA",
        env = "RUSTQC_REFERENCE",
        help_heading = "Input / Output"
    )]
    pub reference: Option<String>,

    /// Output directory [default: .]
    #[arg(
        short,
        long,
        default_value = ".",
        hide_default_value = true,
        env = "RUSTQC_OUTDIR",
        help_heading = "Input / Output"
    )]
    pub outdir: String,

    /// Override sample name for output filenames (default: derived from BAM filename)
    #[arg(
        long,
        value_name = "NAME",
        env = "RUSTQC_SAMPLE_NAME",
        help_heading = "Input / Output"
    )]
    pub sample_name: Option<String>,

    /// Write outputs to a flat directory (no subdirs)
    #[arg(
        long,
        default_value_t = false,
        env = "RUSTQC_FLAT_OUTPUT",
        help_heading = "Input / Output"
    )]
    pub flat_output: bool,

    /// YAML configuration file (see also: RUSTQC_CONFIG env var)
    #[arg(short, long, value_name = "CONFIG", help_heading = "Input / Output")]
    pub config: Option<String>,

    /// JSON summary path (use "-" for stdout)
    #[arg(short = 'j', long = "json-summary", value_name = "PATH", num_args = 0..=1, default_missing_value = "", env = "RUSTQC_JSON_SUMMARY", help_heading = "Input / Output")]
    pub json_summary: Option<String>,

    // ── Library ─────────────────────────────────────────────────────────
    /// Strandedness: unstranded, forward, reverse
    #[arg(
        short,
        long,
        value_enum,
        env = "RUSTQC_STRANDED",
        help_heading = "Library"
    )]
    pub stranded: Option<Strandedness>,

    /// Paired-end reads
    #[arg(short, long, env = "RUSTQC_PAIRED", help_heading = "Library")]
    pub paired: bool,

    // ── General ─────────────────────────────────────────────────────────
    /// Number of threads [default: 1]
    #[arg(
        short,
        long,
        default_value_t = 1,
        hide_default_value = true,
        env = "RUSTQC_THREADS",
        help_heading = "General"
    )]
    pub threads: usize,

    /// MAPQ cutoff for quality filtering [default: 30]
    #[arg(
        short = 'Q',
        long = "mapq",
        default_value_t = 30,
        hide_default_value = true,
        env = "RUSTQC_MAPQ",
        help_heading = "General"
    )]
    pub mapq_cut: u8,

    /// GTF attribute for biotype grouping
    #[arg(
        long,
        value_name = "ATTR",
        env = "RUSTQC_BIOTYPE_ATTRIBUTE",
        help_heading = "General"
    )]
    pub biotype_attribute: Option<String>,

    /// Skip duplicate-marking check
    #[arg(
        long,
        default_value_t = false,
        env = "RUSTQC_SKIP_DUP_CHECK",
        help_heading = "General"
    )]
    pub skip_dup_check: bool,

    /// Suppress output except warnings/errors
    #[arg(
        short = 'q',
        long,
        conflicts_with = "verbose",
        env = "RUSTQC_QUIET",
        help_heading = "General"
    )]
    pub quiet: bool,

    /// Show additional detail
    #[arg(
        short = 'v',
        long,
        conflicts_with = "quiet",
        env = "RUSTQC_VERBOSE",
        help_heading = "General"
    )]
    pub verbose: bool,

    // ── Tool parameters ─────────────────────────────────────────────────
    /// infer_experiment: sample size [default: 200000]
    #[arg(
        long = "infer-experiment-sample-size",
        value_name = "N",
        env = "RUSTQC_INFER_EXPERIMENT_SAMPLE_SIZE",
        help_heading = "Tool parameters"
    )]
    pub infer_experiment_sample_size: Option<u64>,

    /// junction_annotation: min intron size [default: 50]
    #[arg(
        long = "min-intron",
        value_name = "N",
        env = "RUSTQC_MIN_INTRON",
        help_heading = "Tool parameters"
    )]
    pub min_intron: Option<u64>,

    /// junction_saturation: random seed for reproducible results
    #[arg(
        long = "junction-saturation-seed",
        value_name = "N",
        env = "RUSTQC_JUNCTION_SATURATION_SEED",
        help_heading = "Tool parameters"
    )]
    pub junction_saturation_seed: Option<u64>,

    /// junction_saturation: min coverage [default: 1]
    #[arg(
        long = "junction-saturation-min-coverage",
        value_name = "N",
        env = "RUSTQC_JUNCTION_SATURATION_MIN_COVERAGE",
        help_heading = "Tool parameters"
    )]
    pub junction_saturation_min_coverage: Option<u64>,

    /// junction_saturation: start % [default: 5]
    #[arg(
        long = "junction-saturation-percentile-floor",
        value_name = "N",
        env = "RUSTQC_JUNCTION_SATURATION_PERCENTILE_FLOOR",
        help_heading = "Tool parameters"
    )]
    pub junction_saturation_percentile_floor: Option<u64>,

    /// junction_saturation: end % [default: 100]
    #[arg(
        long = "junction-saturation-percentile-ceiling",
        value_name = "N",
        env = "RUSTQC_JUNCTION_SATURATION_PERCENTILE_CEILING",
        help_heading = "Tool parameters"
    )]
    pub junction_saturation_percentile_ceiling: Option<u64>,

    /// junction_saturation: step % [default: 5]
    #[arg(
        long = "junction-saturation-percentile-step",
        value_name = "N",
        env = "RUSTQC_JUNCTION_SATURATION_PERCENTILE_STEP",
        help_heading = "Tool parameters"
    )]
    pub junction_saturation_percentile_step: Option<u64>,

    /// inner_distance: sample size [default: 1000000]
    #[arg(
        long = "inner-distance-sample-size",
        value_name = "N",
        env = "RUSTQC_INNER_DISTANCE_SAMPLE_SIZE",
        help_heading = "Tool parameters"
    )]
    pub inner_distance_sample_size: Option<u64>,

    /// inner_distance: lower bound [default: -250]
    #[arg(
        long = "inner-distance-lower-bound",
        value_name = "N",
        allow_hyphen_values = true,
        env = "RUSTQC_INNER_DISTANCE_LOWER_BOUND",
        help_heading = "Tool parameters"
    )]
    pub inner_distance_lower_bound: Option<i64>,

    /// inner_distance: upper bound [default: 250]
    #[arg(
        long = "inner-distance-upper-bound",
        value_name = "N",
        allow_hyphen_values = true,
        env = "RUSTQC_INNER_DISTANCE_UPPER_BOUND",
        help_heading = "Tool parameters"
    )]
    pub inner_distance_upper_bound: Option<i64>,

    /// inner_distance: bin width [default: 5]
    #[arg(
        long = "inner-distance-step",
        value_name = "N",
        allow_hyphen_values = true,
        env = "RUSTQC_INNER_DISTANCE_STEP",
        help_heading = "Tool parameters"
    )]
    pub inner_distance_step: Option<i64>,

    /// TIN: random seed for reproducible results
    #[arg(
        long = "tin-seed",
        value_name = "N",
        env = "RUSTQC_TIN_SEED",
        help_heading = "Tool parameters"
    )]
    pub tin_seed: Option<u64>,

    /// Skip TIN analysis
    #[arg(
        long,
        default_value_t = false,
        env = "RUSTQC_SKIP_TIN",
        help_heading = "Tool parameters"
    )]
    pub skip_tin: bool,

    /// Skip read duplication analysis
    #[arg(
        long,
        default_value_t = false,
        env = "RUSTQC_SKIP_READ_DUPLICATION",
        help_heading = "Tool parameters"
    )]
    pub skip_read_duplication: bool,

    /// Skip preseq library complexity analysis
    #[arg(
        long,
        default_value_t = false,
        env = "RUSTQC_SKIP_PRESEQ",
        help_heading = "Tool parameters"
    )]
    pub skip_preseq: bool,

    /// preseq: random seed for bootstrap CIs
    #[arg(
        long = "preseq-seed",
        value_name = "N",
        env = "RUSTQC_PRESEQ_SEED",
        help_heading = "Tool parameters"
    )]
    pub preseq_seed: Option<u64>,

    /// preseq: max extrapolation depth
    #[arg(
        long = "preseq-max-extrap",
        value_name = "N",
        env = "RUSTQC_PRESEQ_MAX_EXTRAP",
        help_heading = "Tool parameters"
    )]
    pub preseq_max_extrap: Option<f64>,

    /// preseq: step size between points
    #[arg(
        long = "preseq-step-size",
        value_name = "N",
        env = "RUSTQC_PRESEQ_STEP_SIZE",
        help_heading = "Tool parameters"
    )]
    pub preseq_step_size: Option<f64>,

    /// preseq: bootstrap replicates for CIs
    #[arg(
        long = "preseq-n-bootstraps",
        value_name = "N",
        env = "RUSTQC_PRESEQ_N_BOOTSTRAPS",
        help_heading = "Tool parameters"
    )]
    pub preseq_n_bootstraps: Option<u32>,

    /// preseq: max segment length for PE merging
    #[arg(
        long = "preseq-seg-len",
        value_name = "N",
        env = "RUSTQC_PRESEQ_SEG_LEN",
        help_heading = "Tool parameters"
    )]
    pub preseq_seg_len: Option<i64>,
}

/// Arguments for the `dna` subcommand.
///
/// Shared options keep the same long name, short flag and `RUSTQC_*`
/// environment variable as their `rna` counterparts, so wrapper scripts and
/// muscle memory carry over between the two pipelines. The differences are
/// deliberate: there is no `--gtf` and no `--stranded`, and `--mapq` defaults
/// to 0 rather than 30 because that is mosdepth's default.
#[derive(Parser, Debug)]
#[command(
    next_line_help = false,
    term_width = 120,
    help_template = "\
{about-with-newline}
{usage-heading} {usage}

{all-args}"
)]
pub struct DnaArgs {
    // ── Input / Output ──────────────────────────────────────────────────
    /// Duplicate-marked alignment file(s)
    #[arg(value_name = "INPUT", num_args = 1.., required = true, help_heading = "Input / Output")]
    pub input: Vec<String>,

    /// Reference FASTA (required for CRAM and for GC bias)
    #[arg(
        short,
        long,
        value_name = "FASTA",
        env = "RUSTQC_REFERENCE",
        help_heading = "Input / Output"
    )]
    pub reference: Option<String>,

    /// Target intervals BED; switches on targeted (exome or panel) mode
    #[arg(
        long,
        value_name = "BED",
        env = "RUSTQC_TARGETS",
        help_heading = "Input / Output"
    )]
    pub targets: Option<String>,

    /// Capture bait intervals BED [default: same as --targets]
    #[arg(
        long,
        value_name = "BED",
        env = "RUSTQC_BAITS",
        requires = "targets",
        help_heading = "Input / Output"
    )]
    pub baits: Option<String>,

    /// Output directory [default: .]
    #[arg(
        short,
        long,
        default_value = ".",
        hide_default_value = true,
        env = "RUSTQC_OUTDIR",
        help_heading = "Input / Output"
    )]
    pub outdir: String,

    /// Override sample name for output filenames (default: derived from BAM filename)
    #[arg(
        long,
        value_name = "NAME",
        env = "RUSTQC_SAMPLE_NAME",
        help_heading = "Input / Output"
    )]
    pub sample_name: Option<String>,

    /// Write outputs to a flat directory (no subdirs)
    #[arg(
        long,
        default_value_t = false,
        env = "RUSTQC_FLAT_OUTPUT",
        help_heading = "Input / Output"
    )]
    pub flat_output: bool,

    /// YAML configuration file (see also: RUSTQC_CONFIG env var)
    #[arg(short, long, value_name = "CONFIG", help_heading = "Input / Output")]
    pub config: Option<String>,

    /// JSON summary path (use "-" for stdout)
    #[arg(short = 'j', long = "json-summary", value_name = "PATH", num_args = 0..=1, default_missing_value = "", env = "RUSTQC_JSON_SUMMARY", help_heading = "Input / Output")]
    pub json_summary: Option<String>,

    // ── Library ─────────────────────────────────────────────────────────
    /// Paired-end reads
    #[arg(short, long, env = "RUSTQC_PAIRED", help_heading = "Library")]
    pub paired: bool,

    // ── General ─────────────────────────────────────────────────────────
    /// Number of threads [default: 1]
    #[arg(
        short,
        long,
        default_value_t = 1,
        hide_default_value = true,
        env = "RUSTQC_THREADS",
        help_heading = "General"
    )]
    pub threads: usize,

    /// MAPQ cutoff; reads below it are ignored [default: 0]
    #[arg(
        short = 'Q',
        long = "mapq",
        default_value_t = 0,
        hide_default_value = true,
        env = "RUSTQC_MAPQ",
        help_heading = "General"
    )]
    pub mapq_cut: u8,

    /// Skip duplicate-marking check
    #[arg(
        long,
        default_value_t = false,
        env = "RUSTQC_SKIP_DUP_CHECK",
        help_heading = "General"
    )]
    pub skip_dup_check: bool,

    /// Suppress output except warnings/errors
    #[arg(
        short = 'q',
        long,
        conflicts_with = "verbose",
        env = "RUSTQC_QUIET",
        help_heading = "General"
    )]
    pub quiet: bool,

    /// Show additional detail
    #[arg(
        short = 'v',
        long,
        conflicts_with = "quiet",
        env = "RUSTQC_VERBOSE",
        help_heading = "General"
    )]
    pub verbose: bool,

    // ── Tool parameters ─────────────────────────────────────────────────
    /// Coverage thresholds to report [default: 1,5,10,15,20,30,50]
    #[arg(
        long = "depth-thresholds",
        value_name = "N,...",
        value_delimiter = ',',
        default_values_t = vec![1u32, 5, 10, 15, 20, 30, 50],
        hide_default_value = true,
        env = "RUSTQC_DEPTH_THRESHOLDS",
        help_heading = "Tool parameters"
    )]
    pub depth_thresholds: Vec<u32>,

    /// Fixed-width window size for per-window depth
    #[arg(
        long = "window-size",
        value_name = "N",
        env = "RUSTQC_WINDOW_SIZE",
        help_heading = "Tool parameters"
    )]
    pub window_size: Option<u32>,

    /// Picard COVERAGE_CAP [default: 250]
    #[arg(
        long = "coverage-cap",
        value_name = "N",
        default_value_t = 250,
        hide_default_value = true,
        env = "RUSTQC_COVERAGE_CAP",
        help_heading = "Tool parameters"
    )]
    pub coverage_cap: u32,

    /// Picard MINIMUM_BASE_QUALITY [default: 20]
    #[arg(
        long = "min-base-quality",
        value_name = "N",
        default_value_t = 20,
        hide_default_value = true,
        env = "RUSTQC_MIN_BASE_QUALITY",
        help_heading = "Tool parameters"
    )]
    pub min_base_quality: u8,

    /// Skip the per-base depth output, by far the largest file
    #[arg(
        long,
        default_value_t = false,
        env = "RUSTQC_SKIP_PER_BASE",
        help_heading = "Tool parameters"
    )]
    pub skip_per_base: bool,

    /// Skip GC bias metrics
    #[arg(
        long,
        default_value_t = false,
        env = "RUSTQC_SKIP_GC_BIAS",
        help_heading = "Tool parameters"
    )]
    pub skip_gc_bias: bool,

    /// Cap on concurrently live per-contig depth arrays [default: derived from RAM]
    #[arg(
        long = "max-depth-workers",
        value_name = "N",
        env = "RUSTQC_MAX_DEPTH_WORKERS",
        help_heading = "Tool parameters"
    )]
    pub max_depth_workers: Option<usize>,

    /// Skip preseq library complexity analysis
    #[arg(
        long,
        default_value_t = false,
        env = "RUSTQC_SKIP_PRESEQ",
        help_heading = "Tool parameters"
    )]
    pub skip_preseq: bool,

    /// preseq: random seed for bootstrap CIs
    #[arg(
        long = "preseq-seed",
        value_name = "N",
        env = "RUSTQC_PRESEQ_SEED",
        help_heading = "Tool parameters"
    )]
    pub preseq_seed: Option<u64>,

    /// preseq: max extrapolation depth
    #[arg(
        long = "preseq-max-extrap",
        value_name = "N",
        env = "RUSTQC_PRESEQ_MAX_EXTRAP",
        help_heading = "Tool parameters"
    )]
    pub preseq_max_extrap: Option<f64>,

    /// preseq: step size between points
    #[arg(
        long = "preseq-step-size",
        value_name = "N",
        env = "RUSTQC_PRESEQ_STEP_SIZE",
        help_heading = "Tool parameters"
    )]
    pub preseq_step_size: Option<f64>,

    /// preseq: bootstrap replicates for CIs
    #[arg(
        long = "preseq-n-bootstraps",
        value_name = "N",
        env = "RUSTQC_PRESEQ_N_BOOTSTRAPS",
        help_heading = "Tool parameters"
    )]
    pub preseq_n_bootstraps: Option<u32>,

    /// preseq: max segment length for PE merging
    #[arg(
        long = "preseq-seg-len",
        value_name = "N",
        env = "RUSTQC_PRESEQ_SEG_LEN",
        help_heading = "Tool parameters"
    )]
    pub preseq_seg_len: Option<i64>,
}

/// Parse command-line arguments and return the Cli struct.
///
/// Sets a `long_version` that includes the git commit, build timestamp,
/// and CPU info line, shown when the user runs `rustqc --version`.
pub fn parse_args() -> Cli {
    use clap::FromArgMatches;

    let long_version: &'static str = Box::leak(
        format!(
            "{} ({}, built {})\n{}",
            env!("CARGO_PKG_VERSION"),
            env!("GIT_SHORT_HASH"),
            env!("BUILD_TIMESTAMP"),
            rustqc::cpu::cpu_info_line(),
        )
        .into_boxed_str(),
    );
    let cmd = Cli::command().long_version(long_version);
    let matches = cmd.get_matches();
    Cli::from_arg_matches(&matches).expect("clap arg matching failed")
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
            _ => panic!("Expected Rna subcommand"),
        }
    }

    #[test]
    fn test_dna_default_args() {
        let cli = Cli::parse_from(["rustqc", "dna", "test.bam"]);
        match cli.command {
            Commands::Dna(args) => {
                assert_eq!(args.input, vec!["test.bam"]);
                assert_eq!(args.outdir, ".");
                assert_eq!(args.threads, 1);
                assert_eq!(args.mapq_cut, 0);
                assert_eq!(args.coverage_cap, 250);
                assert_eq!(args.min_base_quality, 20);
                assert_eq!(args.depth_thresholds, vec![1, 5, 10, 15, 20, 30, 50]);
                assert_eq!(args.window_size, None);
                assert!(args.targets.is_none());
                assert!(args.baits.is_none());
                assert!(!args.skip_per_base);
                assert!(!args.skip_gc_bias);
            }
            _ => panic!("Expected Dna subcommand"),
        }
    }

    #[test]
    fn test_dna_no_gtf_required() {
        assert!(Cli::try_parse_from(["rustqc", "dna", "test.bam"]).is_ok());
    }

    #[test]
    fn test_dna_targeted_args() {
        let cli = Cli::parse_from([
            "rustqc",
            "dna",
            "a.bam",
            "b.bam",
            "--targets",
            "t.bed",
            "--baits",
            "b.bed",
            "--depth-thresholds",
            "1,10,100",
            "--window-size",
            "500",
            "--reference",
            "genome.fa",
            "-Q",
            "20",
            "--threads",
            "4",
        ]);
        match cli.command {
            Commands::Dna(args) => {
                assert_eq!(args.input, vec!["a.bam", "b.bam"]);
                assert_eq!(args.targets, Some("t.bed".to_string()));
                assert_eq!(args.baits, Some("b.bed".to_string()));
                assert_eq!(args.depth_thresholds, vec![1, 10, 100]);
                assert_eq!(args.window_size, Some(500));
                assert_eq!(args.reference, Some("genome.fa".to_string()));
                assert_eq!(args.mapq_cut, 20);
                assert_eq!(args.threads, 4);
            }
            _ => panic!("Expected Dna subcommand"),
        }
    }

    #[test]
    fn test_dna_baits_without_targets_is_rejected() {
        let result = Cli::try_parse_from(["rustqc", "dna", "test.bam", "--baits", "b.bed"]);
        assert!(
            result.is_err(),
            "--baits without --targets must be rejected"
        );
    }

    #[test]
    fn test_protein_sequence_default_args() {
        let cli = Cli::parse_from(["rustqc", "protein", "sequence", "proteome.fa"]);
        match cli.command {
            Commands::Protein(args) => match args.mode {
                ProteinMode::Sequence(args) => {
                    assert_eq!(args.input, vec!["proteome.fa"]);
                    assert_eq!(args.outdir, ".");
                    assert_eq!(args.min_length, 0);
                    assert!(!args.expect_stop);
                }
            },
            _ => panic!("Expected Protein subcommand"),
        }
    }

    #[test]
    fn test_protein_requires_a_mode() {
        let result = Cli::try_parse_from(["rustqc", "protein", "proteome.fa"]);
        assert!(
            result.is_err(),
            "the mode is explicit, so a bare file argument must be rejected"
        );
    }

    #[test]
    fn test_protein_sequence_multiple_inputs_and_flags() {
        let cli = Cli::parse_from([
            "rustqc",
            "protein",
            "sequence",
            "a.fa",
            "b.fa.gz",
            "--min-length",
            "50",
            "--expect-stop",
            "--outdir",
            "/tmp/out",
        ]);
        match cli.command {
            Commands::Protein(args) => match args.mode {
                ProteinMode::Sequence(args) => {
                    assert_eq!(args.input, vec!["a.fa", "b.fa.gz"]);
                    assert_eq!(args.min_length, 50);
                    assert!(args.expect_stop);
                    assert_eq!(args.outdir, "/tmp/out");
                }
            },
            _ => panic!("Expected Protein subcommand"),
        }
    }
}
