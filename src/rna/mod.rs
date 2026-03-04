//! RNA-Seq quality control and analysis modules.
//!
//! Contains dupRadar duplication rate analysis, featureCounts-compatible output,
//! and RSeQC tool reimplementations.

pub mod bam_flags;
pub mod dupradar;
pub mod featurecounts;
pub mod preseq;
pub mod qualimap;
pub mod rseqc;
