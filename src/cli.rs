//! Command-line interface definition for RustQC.
//!
//! Provides a subcommand-based CLI where each QC tool is a separate subcommand.
//! The `rna` subcommand provides duplication rate analysis compatible with the
//! dupRadar R package.

use clap::{Parser, Subcommand};

/// Fast quality control tools for sequencing data, written in Rust.
///
/// RustQC provides a collection of QC analysis tools. Use a subcommand to
/// select the analysis to run.
#[derive(Parser, Debug)]
#[command(name = "rustqc", version, about, long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

/// Available analysis subcommands.
#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Analyze PCR duplicate rates in RNA-Seq data (dupRadar equivalent).
    ///
    /// Analyzes PCR duplicate rates as a function of gene expression level
    /// in RNA-Seq datasets. A Rust reimplementation of the dupRadar R package.
    ///
    /// Input alignment files (SAM/BAM/CRAM) must have duplicates marked (SAM flag 0x400)
    /// but NOT removed. Use tools like Picard MarkDuplicates or samblaster to mark
    /// duplicates first.
    Rna(RnaArgs),
}

/// Arguments for the `rna` (dupRadar) subcommand.
#[derive(Parser, Debug)]
pub struct RnaArgs {
    /// Path(s) to duplicate-marked alignment file(s) (SAM/BAM/CRAM)
    #[arg(value_name = "INPUT", num_args = 1.., required = true)]
    pub input: Vec<String>,

    /// Path to the GTF gene annotation file
    #[arg(short, long, value_name = "GTF")]
    pub gtf: String,

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

    /// Path to reference FASTA file (required for CRAM input)
    #[arg(short, long, value_name = "FASTA")]
    pub reference: Option<String>,

    /// Skip the check for duplicate-marking tools in the BAM header.
    ///
    /// By default, RustQC verifies that the BAM file has been processed by a
    /// duplicate-marking tool (e.g. Picard MarkDuplicates, samblaster) before
    /// running. Use this flag to bypass that check.
    #[arg(long, default_value_t = false)]
    pub skip_dup_check: bool,
}

/// Parse command-line arguments and return the Cli struct.
pub fn parse_args() -> Cli {
    Cli::parse()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rna_default_args() {
        // Test that defaults are sensible with a single BAM
        let cli = Cli::parse_from(["rustqc", "rna", "test.bam", "--gtf", "genes.gtf"]);
        match cli.command {
            Commands::Rna(args) => {
                assert_eq!(args.input, vec!["test.bam"]);
                assert_eq!(args.gtf, "genes.gtf");
                assert_eq!(args.stranded, 0);
                assert!(!args.paired);
                assert_eq!(args.threads, 1);
                assert_eq!(args.outdir, ".");
            }
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
        }
    }

    #[test]
    fn test_rna_all_args() {
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
        ]);
        match cli.command {
            Commands::Rna(args) => {
                assert_eq!(args.stranded, 2);
                assert!(args.paired);
                assert_eq!(args.threads, 4);
                assert_eq!(args.outdir, "/tmp/out");
                assert_eq!(args.reference, Some("genome.fa".to_string()));
            }
        }
    }
}
