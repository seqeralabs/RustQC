//! Picard `CollectRnaSeqMetrics` reimplementation.
//!
//! Reports where aligned bases fall relative to the annotation: coding, UTR,
//! intronic or intergenic. It belongs to the `rna` pipeline, which already
//! takes both an alignment and a gene annotation, so it costs no new
//! command-line surface.

pub mod output;
pub mod regions;

pub use regions::{Region, RegionSets};
