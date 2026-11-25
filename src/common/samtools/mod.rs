//! samtools-compatible output writers.
//!
//! Reproduce the exact output formats of `samtools stats`, `samtools flagstat`
//! and `samtools idxstats` from the counters gathered in
//! [`crate::common::bam_stat::BamStatResult`], so that MultiQC and
//! `plot-bamstats` parse RustQC output as if samtools had produced it.

pub mod flagstat;
pub mod idxstats;
pub mod stats;
