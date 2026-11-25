//! Generic BAM/CRAM quality control — the `bamqc` equivalent of Qualimap.
//!
//! Where the `rna` command reimplements Qualimap's `rnaseq` mode (read origin,
//! transcript coverage bias, junctions — all annotation-driven), this module
//! covers the `bamqc` mode: coverage depth and breadth, GC content, insert
//! size, mapping quality and per-contig coverage. No GTF is required.

pub mod accumulator;
pub mod output;

pub use accumulator::{BamqcAccum, BamqcResult, ContigCoverage};
pub use output::write_all;
