//! RNA-Seq quality control and analysis modules.
//!
//! Contains dupRadar duplication rate analysis, featureCounts-compatible output,
//! and RSeQC tool reimplementations.

pub mod dupradar;
pub mod featurecounts;
pub mod qualimap;
pub mod rnaseq_metrics;
pub mod rseqc;

// These analyses are not RNA-specific and now live in `crate::common`.
// Re-exported here so existing `crate::rna::...` paths and the published
// 0.2.x library surface keep resolving. Drop the shims at 1.0.
pub use crate::common::{bam_flags, cpp_rng, preseq};
