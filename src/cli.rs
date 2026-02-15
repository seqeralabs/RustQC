//! Command-line interface definition for RustQC.
//!
//! Provides a subcommand-based CLI. The `rna` subcommand runs all RNA-Seq QC
//! analyses in a single pass: dupRadar duplication rate analysis, featureCounts-
//! compatible output, and RSeQC-equivalent metrics (bam_stat, infer_experiment,
//! read_duplication, read_distribution, junction_annotation, junction_saturation,
//! inner_distance). Individual tools can be disabled via the YAML config file.

use clap::{Parser, Subcommand};

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
    ///
    /// Provide either --gtf or --bed (mutually exclusive). With --gtf all tools
    /// run. With --bed only RSeQC tools run (dupRadar and featureCounts are
    /// skipped). Disable individual tools via the YAML config file.
    ///
    /// Input alignment files must have duplicates marked (SAM flag 0x400) but NOT
    /// removed. Use tools like Picard MarkDuplicates or samblaster first.
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
    /// When provided, all analyses are run: dupRadar, featureCounts, and all
    /// RSeQC tools. The GTF must contain exon features with gene_id attributes.
    /// CDS features are used for read_distribution UTR/CDS classification.
    /// Gzip-compressed files (.gtf.gz) are detected and decompressed automatically.
    /// Mutually exclusive with --bed.
    #[arg(
        short,
        long,
        value_name = "GTF",
        conflicts_with = "bed",
        required_unless_present = "bed"
    )]
    pub gtf: Option<String>,

    /// Path to a BED12 gene model file (plain or gzip-compressed).
    ///
    /// When provided instead of --gtf, only RSeQC tools are run (dupRadar and
    /// featureCounts are skipped because BED files lack gene-level grouping and
    /// biotype information). Gzip-compressed files (.bed.gz) are detected and
    /// decompressed automatically. Mutually exclusive with --gtf.
    #[arg(
        short,
        long,
        value_name = "BED",
        conflicts_with = "gtf",
        required_unless_present = "gtf"
    )]
    pub bed: Option<String>,

    /// Library strandedness: 0=unstranded, 1=forward, 2=reverse
    #[arg(short, long, default_value_t = 0, value_parser = clap::value_parser!(u8).range(0..=2))]
    pub stranded: u8,

    /// Whether the library is paired-end
    #[arg(short, long, default_value_t = false)]
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

    // === RSeQC tool-specific parameters ===
    /// MAPQ cutoff for read quality filtering (used by bam_stat, infer_experiment,
    /// read_duplication, junction_annotation, junction_saturation, inner_distance)
    #[arg(short = 'q', long = "mapq", default_value_t = 30)]
    pub mapq_cut: u8,

    // --- infer_experiment ---
    /// Maximum number of reads to sample for strandedness inference
    #[arg(long = "infer-experiment-sample-size", default_value_t = 200_000)]
    pub infer_experiment_sample_size: u64,

    // --- junction_annotation ---
    /// Minimum intron size for junction annotation and saturation analysis
    #[arg(long = "min-intron", default_value_t = 50)]
    pub min_intron: u64,

    // --- junction_saturation ---
    /// Minimum read coverage to count a known junction (junction_saturation)
    #[arg(long = "junction-saturation-min-coverage", default_value_t = 1)]
    pub junction_saturation_min_coverage: u64,

    /// Sampling start percentage for junction saturation
    #[arg(long = "junction-saturation-percentile-floor", default_value_t = 5)]
    pub junction_saturation_percentile_floor: u64,

    /// Sampling end percentage for junction saturation
    #[arg(long = "junction-saturation-percentile-ceiling", default_value_t = 100)]
    pub junction_saturation_percentile_ceiling: u64,

    /// Sampling step percentage for junction saturation
    #[arg(long = "junction-saturation-percentile-step", default_value_t = 5)]
    pub junction_saturation_percentile_step: u64,

    // --- inner_distance ---
    /// Maximum number of read pairs to sample for inner distance
    #[arg(long = "inner-distance-sample-size", default_value_t = 1_000_000)]
    pub inner_distance_sample_size: u64,

    /// Lower bound of inner distance histogram
    #[arg(long = "inner-distance-lower-bound", default_value_t = -250, allow_hyphen_values = true)]
    pub inner_distance_lower_bound: i64,

    /// Upper bound of inner distance histogram
    #[arg(
        long = "inner-distance-upper-bound",
        default_value_t = 250,
        allow_hyphen_values = true
    )]
    pub inner_distance_upper_bound: i64,

    /// Step size (bin width) for inner distance histogram
    #[arg(
        long = "inner-distance-step",
        default_value_t = 5,
        allow_hyphen_values = true
    )]
    pub inner_distance_step: i64,

    // --- preseq lc_extrap ---
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
                assert_eq!(args.gtf, Some("genes.gtf".to_string()));
                assert_eq!(args.stranded, 0);
                assert!(!args.paired);
                assert_eq!(args.threads, 1);
                assert_eq!(args.outdir, ".");
                assert!(args.biotype_attribute.is_none());
                assert!(args.bed.is_none());
                assert_eq!(args.mapq_cut, 30);
                assert_eq!(args.infer_experiment_sample_size, 200_000);
                assert_eq!(args.min_intron, 50);
                assert_eq!(args.inner_distance_step, 5);
            }
            #[allow(unreachable_patterns)]
            _ => panic!("Expected Rna subcommand"),
        }
    }

    #[test]
    fn test_rna_default_args_bed() {
        // Test that BED-only mode works
        let cli = Cli::parse_from(["rustqc", "rna", "test.bam", "--bed", "model.bed"]);
        match cli.command {
            Commands::Rna(args) => {
                assert_eq!(args.input, vec!["test.bam"]);
                assert!(args.gtf.is_none());
                assert_eq!(args.bed, Some("model.bed".to_string()));
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
                assert_eq!(args.gtf, Some("genes.gtf".to_string()));
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
            "-q",
            "20",
        ]);
        match cli.command {
            Commands::Rna(args) => {
                assert_eq!(args.gtf, Some("genes.gtf".to_string()));
                assert!(args.bed.is_none());
                assert_eq!(args.stranded, 2);
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
    fn test_rna_bed_all_args() {
        let cli = Cli::parse_from([
            "rustqc",
            "rna",
            "test.bam",
            "--bed",
            "model.bed",
            "--stranded",
            "2",
            "--paired",
            "--threads",
            "4",
            "--outdir",
            "/tmp/out",
            "--reference",
            "genome.fa",
            "-q",
            "20",
        ]);
        match cli.command {
            Commands::Rna(args) => {
                assert!(args.gtf.is_none());
                assert_eq!(args.bed, Some("model.bed".to_string()));
                assert_eq!(args.stranded, 2);
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
    fn test_rna_gtf_bed_mutually_exclusive() {
        // Providing both --gtf and --bed should fail
        let result = Cli::try_parse_from([
            "rustqc",
            "rna",
            "test.bam",
            "--gtf",
            "genes.gtf",
            "--bed",
            "model.bed",
        ]);
        assert!(
            result.is_err(),
            "Expected error when both --gtf and --bed are provided"
        );
    }

    #[test]
    fn test_rna_neither_gtf_nor_bed() {
        // Providing neither --gtf nor --bed should fail
        let result = Cli::try_parse_from(["rustqc", "rna", "test.bam"]);
        assert!(
            result.is_err(),
            "Expected error when neither --gtf nor --bed is provided"
        );
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
                assert_eq!(args.infer_experiment_sample_size, 500_000);
                assert_eq!(args.min_intron, 100);
                assert_eq!(args.junction_saturation_min_coverage, 5);
                assert_eq!(args.inner_distance_lower_bound, -500);
                assert_eq!(args.inner_distance_upper_bound, 500);
                assert_eq!(args.inner_distance_step, 10);
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
