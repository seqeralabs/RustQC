//! Command-line interface definition for dupRust.
//!
//! Provides argument parsing compatible with the dupRadar R package CLI wrapper,
//! accepting BAM file, GTF annotation, strandedness, and pairing information.

use clap::Parser;

/// Fast duplication rate analysis for RNA-Seq data.
///
/// dupRust analyzes PCR duplicate rates as a function of gene expression level
/// in RNA-Seq datasets. It is a Rust reimplementation of the dupRadar R package.
///
/// Input BAM files must have duplicates marked (SAM flag 0x400) but NOT removed.
/// Use tools like Picard MarkDuplicates or samblaster to mark duplicates first.
#[derive(Parser, Debug)]
#[command(name = "duprust", version, about, long_about = None)]
pub struct Args {
    /// Path to the duplicate-marked BAM file
    #[arg(value_name = "BAM")]
    pub bam: String,

    /// Path to the GTF gene annotation file
    #[arg(value_name = "GTF")]
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
}

/// Parse command-line arguments and return the Args struct.
pub fn parse_args() -> Args {
    Args::parse()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_args() {
        // Test that defaults are sensible
        let args = Args::parse_from(["duprust", "test.bam", "genes.gtf"]);
        assert_eq!(args.stranded, 0);
        assert!(!args.paired);
        assert_eq!(args.threads, 1);
        assert_eq!(args.outdir, ".");
    }

    #[test]
    fn test_all_args() {
        let args = Args::parse_from([
            "duprust",
            "test.bam",
            "genes.gtf",
            "--stranded",
            "2",
            "--paired",
            "--threads",
            "4",
            "--outdir",
            "/tmp/out",
        ]);
        assert_eq!(args.stranded, 2);
        assert!(args.paired);
        assert_eq!(args.threads, 4);
        assert_eq!(args.outdir, "/tmp/out");
    }
}
