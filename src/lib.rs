//! RustQC — fast quality control tools for sequencing data.
//!
//! This is the library API. The companion CLI (`rustqc` binary) is built on
//! top of these same modules and provides a single-pass RNA-Seq QC pipeline
//! that runs dupRadar, featureCounts, RSeQC tools, preseq, samtools-style
//! outputs, and Qualimap analyses.
//!
//! Library consumers can drive individual analyses directly. The submodules
//! are organised by tool family:
//!
//! - [`gtf`] — GTF gene-annotation parsing.
//! - [`io`] — shared I/O helpers (transparent gzip decompression, etc.).
//! - [`config`] — configuration types (mirrors the CLI's YAML config).
//! - [`summary`] — serializable types for the JSON run summary.
//! - [`cpu`] — CPU feature detection and binary-target identification.
//! - [`rna`] — the RNA-Seq QC analysis modules (dupRadar, featureCounts,
//!   RSeQC, Qualimap, preseq, samtools-style outputs).
//!
//! The [`Strandedness`] enum lives at the crate root because it is used
//! across most analysis modules.

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
