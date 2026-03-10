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

use clap::{Parser, Subcommand};

/// Library strandedness protocol.
///
/// Determines how read strand is interpreted relative to the gene annotation
/// strand during counting.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Strandedness {
    /// Count reads on either strand (library is not strand-specific).
    #[default]
    Unstranded,
    /// Forward stranded: read 1 maps to the transcript strand.
    Forward,
    /// Reverse stranded: read 2 maps to the transcript strand (e.g. dUTP).
    Reverse,
}

impl Strandedness {
    /// Convert from the integer convention used by featureCounts / HTSeq.
    ///
    /// Returns `None` for values outside 0..=2.
    pub fn from_u8(v: u8) -> Option<Self> {
        match v {
            0 => Some(Strandedness::Unstranded),
            1 => Some(Strandedness::Forward),
            2 => Some(Strandedness::Reverse),
            _ => None,
        }
    }
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
///
/// RustQC runs a comprehensive suite of RNA-Seq QC analyses in a single pass
/// over your BAM file(s). Use `rustqc rna` to get started.
#[derive(Parser, Debug)]
#[command(name = "rustqc", version, about, long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

/// Available analysis subcommands.
#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Run all RNA-Seq QC analyses in a single pass.
    ///
    /// Processes BAM/SAM/CRAM files through the complete RNA-Seq QC pipeline:
    ///
    /// - dupRadar: PCR duplicate rate analysis as a function of gene expression
    /// - featureCounts: Gene-level read counting with biotype summaries
    /// - bam_stat: Basic alignment statistics
    /// - infer_experiment: Library strandedness inference
    /// - read_duplication: Position- and sequence-based duplication histograms
    /// - read_distribution: Read classification into genomic features
    /// - junction_annotation: Splice junction classification
    /// - junction_saturation: Junction discovery saturation curves
    /// - inner_distance: Insert size estimation for paired-end data
    /// - TIN: Transcript Integrity Number analysis
    /// - preseq: Library complexity extrapolation (lc_extrap)
    /// - Qualimap rnaseq: Gene body coverage and read origin profiling
    /// - flagstat / idxstats / stats: samtools-compatible alignment summaries
    ///
    /// A GTF gene annotation file (--gtf) is required. Input alignment files
    /// must have duplicates marked (SAM flag 0x400) but NOT removed. Use tools
    /// like Picard MarkDuplicates or samblaster first.
    Rna(RnaArgs),
}

/// Arguments for the `rna` subcommand.
///
/// Runs all RNA-Seq QC analyses in a single pass. Tool-specific parameters
/// have sensible defaults and can also be set via the YAML config file.
/// CLI flags override config file settings.
#[derive(Parser, Debug)]
pub struct RnaArgs {
    /// Path(s) to duplicate-marked alignment file(s) (SAM/BAM/CRAM)
    #[arg(value_name = "INPUT", num_args = 1.., required = true)]
    pub input: Vec<String>,

    /// Path to a GTF gene annotation file (plain or gzip-compressed).
    ///
    /// The GTF must contain exon features with gene_id attributes. CDS features
    /// are used for read_distribution UTR/CDS classification. Gzip-compressed
    /// files (.gtf.gz) are detected and decompressed automatically.
    #[arg(short, long, value_name = "GTF")]
    pub gtf: String,

    /// Library strandedness: 0=unstranded, 1=forward, 2=reverse
    #[arg(short, long, value_parser = clap::value_parser!(u8).range(0..=2))]
    pub stranded: Option<u8>,

    /// Whether the library is paired-end
    #[arg(short, long)]
    pub paired: bool,

    /// Number of threads for parallel processing
    #[arg(short, long, default_value_t = 1)]
    pub threads: usize,

    /// Output directory for results
    #[arg(short, long, default_value = ".")]
    pub outdir: String,

    /// Path to a YAML configuration file (e.g. chromosome name mapping)
    #[arg(short, long, value_name = "CONFIG")]
    pub config: Option<String>,

    /// GTF attribute for biotype grouping (e.g. gene_biotype, gene_type).
    ///
    /// Overrides the biotype_attribute setting in the config file. If neither
    /// this flag nor the config file specifies a value, defaults to "gene_biotype".
    /// Set to empty string to disable biotype counting.
    #[arg(long, value_name = "ATTRIBUTE")]
    pub biotype_attribute: Option<String>,

    /// Path to reference FASTA file (required for CRAM input)
    #[arg(short, long, value_name = "FASTA")]
    pub reference: Option<String>,

    /// Write all output files to a flat directory (no subdirectories).
    ///
    /// By default, RustQC organises output files into subdirectories by tool:
    /// `dupradar/`, `featurecounts/`, and `rseqc/{tool}/`. When this flag is
    /// set, all output files are written directly to the output directory with
    /// no subdirectories (the legacy behaviour).
    #[arg(long, default_value_t = false)]
    pub flat_output: bool,

    /// Skip the check for duplicate-marking tools in the BAM header.
    ///
    /// By default, RustQC verifies that the BAM file has been processed by a
    /// duplicate-marking tool (e.g. Picard MarkDuplicates, samblaster) before
    /// running. Use this flag to bypass that check.
    #[arg(long, default_value_t = false)]
    pub skip_dup_check: bool,

    /// Suppress all output except warnings and errors
    #[arg(short = 'q', long, conflicts_with = "verbose")]
    pub quiet: bool,

    /// Show additional detail (per-file writes, intermediate stats)
    #[arg(short = 'v', long, conflicts_with = "quiet")]
    pub verbose: bool,

    /// Write a JSON summary to the given path (or to <outdir>/rustqc_summary.json if no path given).
    ///
    /// Use "-" to write JSON to stdout. Useful for automation and AI agents.
    #[arg(short = 'j', long = "json-summary", value_name = "PATH", num_args = 0..=1, default_missing_value = "")]
    pub json_summary: Option<String>,

    // === RSeQC tool-specific parameters ===
    /// MAPQ cutoff for read quality filtering (used by bam_stat, infer_experiment,
    /// read_duplication, junction_annotation, junction_saturation, inner_distance)
    #[arg(short = 'Q', long = "mapq", default_value_t = 30)]
    pub mapq_cut: u8,

    // --- infer_experiment ---
    /// Maximum number of reads to sample for strandedness inference [default: 200000]
    #[arg(long = "infer-experiment-sample-size")]
    pub infer_experiment_sample_size: Option<u64>,

    // --- junction_annotation ---
    /// Minimum intron size for junction annotation and saturation analysis [default: 50]
    #[arg(long = "min-intron")]
    pub min_intron: Option<u64>,

    // --- junction_saturation ---
    /// Minimum read coverage to count a known junction (junction_saturation) [default: 1]
    #[arg(long = "junction-saturation-min-coverage")]
    pub junction_saturation_min_coverage: Option<u64>,

    /// Sampling start percentage for junction saturation [default: 5]
    #[arg(long = "junction-saturation-percentile-floor")]
    pub junction_saturation_percentile_floor: Option<u64>,

    /// Sampling end percentage for junction saturation [default: 100]
    #[arg(long = "junction-saturation-percentile-ceiling")]
    pub junction_saturation_percentile_ceiling: Option<u64>,

    /// Sampling step percentage for junction saturation [default: 5]
    #[arg(long = "junction-saturation-percentile-step")]
    pub junction_saturation_percentile_step: Option<u64>,

    // --- inner_distance ---
    /// Maximum number of read pairs to sample for inner distance [default: 1000000]
    #[arg(long = "inner-distance-sample-size")]
    pub inner_distance_sample_size: Option<u64>,

    /// Lower bound of inner distance histogram [default: -250]
    #[arg(long = "inner-distance-lower-bound", allow_hyphen_values = true)]
    pub inner_distance_lower_bound: Option<i64>,

    /// Upper bound of inner distance histogram [default: 250]
    #[arg(long = "inner-distance-upper-bound", allow_hyphen_values = true)]
    pub inner_distance_upper_bound: Option<i64>,

    /// Step size (bin width) for inner distance histogram [default: 5]
    #[arg(long = "inner-distance-step", allow_hyphen_values = true)]
    pub inner_distance_step: Option<i64>,

    // --- TIN / read_duplication / preseq ---
    /// Skip the TIN (Transcript Integrity Number) analysis
    #[arg(long, default_value_t = false)]
    pub skip_tin: bool,

    /// Skip the read duplication analysis
    #[arg(long, default_value_t = false)]
    pub skip_read_duplication: bool,

    /// Skip the preseq library complexity extrapolation analysis
    #[arg(long, default_value_t = false)]
    pub skip_preseq: bool,

    /// Maximum extrapolation depth in total reads (preseq lc_extrap)
    #[arg(long = "preseq-max-extrap")]
    pub preseq_max_extrap: Option<f64>,

    /// Step size between extrapolation points (preseq lc_extrap)
    #[arg(long = "preseq-step-size")]
    pub preseq_step_size: Option<f64>,

    /// Number of bootstrap replicates for confidence intervals (preseq lc_extrap)
    #[arg(long = "preseq-n-bootstraps")]
    pub preseq_n_bootstraps: Option<u32>,
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
            "2",
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
                assert_eq!(args.stranded, Some(2)); // CLI still parses u8
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
        ]);
        match cli.command {
            Commands::Rna(args) => {
                assert!(!args.skip_preseq);
                assert_eq!(args.preseq_max_extrap, Some(5_000_000_000.0));
                assert_eq!(args.preseq_step_size, Some(500_000.0));
                assert_eq!(args.preseq_n_bootstraps, Some(200));
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
