//! RustQC — fast quality control tools for sequencing data.
//!
//! RustQC is primarily a CLI (`rustqc rna ...`) that runs a single-pass
//! RNA-Seq QC pipeline (dupRadar, featureCounts, 8 RSeQC tools, Qualimap,
//! bigWig coverage tracks, preseq, samtools-style outputs). The same analysis modules are also
//! exposed as a library so they can be embedded into other Rust programs.
//!
//! # Adding RustQC as a dependency
//!
//! ```toml
//! [dependencies]
//! rustqc = "0.2"
//! ```
//!
//! The library pulls in `rust-htslib` (linked statically), `plotters`, and
//! a small C++ component used by the preseq tool (built via `build.rs`),
//! so a working C/C++ toolchain is required at build time.
//!
//! # Modules
//!
//! - [`gtf`] — GTF gene-annotation parsing into [`gtf::Gene`] / [`gtf::Transcript`] / [`gtf::Exon`].
//! - [`io`] — shared I/O helpers (transparent gzip decompression, FNV-1a, number formatting).
//! - [`config`] — configuration types that mirror the CLI's YAML config file.
//! - [`summary`] — serializable types for the JSON run summary.
//! - [`cpu`] — CPU feature detection and binary-target identification.
//! - [`rna`] — the RNA-Seq analysis modules:
//!   - [`rna::dupradar`], [`rna::featurecounts`], [`rna::qualimap`],
//!     [`rna::bigwig`], [`rna::preseq`], [`rna::rseqc`].
//!
//! [`Strandedness`] lives at the crate root because it is used across most
//! analysis modules.
//!
//! # Stability
//!
//! The library is at `0.2.x` and the public surface is intentionally small
//! at this stage. Expect breaking changes in minor releases until `1.0`.
//! The full single-pass RNA-Seq pipeline (the `run_rna` orchestrator that
//! the binary uses) is not yet exposed as a library entry point — for now
//! library consumers drive individual analyses themselves. Pipeline-level
//! orchestration may be exposed in a future release; see issue
//! [#72](https://github.com/seqeralabs/RustQC/issues/72).
//!
//! # Examples
//!
//! Parse a GTF file and inspect the first gene:
//!
//! ```no_run
//! use rustqc::gtf;
//!
//! let genes = gtf::parse_gtf("genes.gtf", &[]).unwrap();
//! if let Some((gene_id, gene)) = genes.iter().next() {
//!     println!("{gene_id}: {} transcripts", gene.transcripts.len());
//! }
//! ```
//!
//! Use the [`Strandedness`] enum (also accepted by `serde` for YAML configs):
//!
//! ```
//! use rustqc::Strandedness;
//!
//! let s = Strandedness::Reverse;
//! assert_eq!(s.to_string(), "reverse");
//! ```

use clap::ValueEnum;
use serde::Deserialize;

pub mod config;
pub mod cpu;
pub mod gtf;
pub mod io;
pub mod rna;
pub mod summary;

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
