//! Single-pass alignment QC for DNA pipelines (`rustqc align`).
//!
//! Replaces three independent passes over a CRAM/BAM — `samtools stats`,
//! `mosdepth`, and the `bcftools mpileup` genotyping step feeding
//! NGSCheckMate — with one streaming pass that computes all three.
//!
//! The genotyping component is close to free: for a coordinate-sorted file the
//! SNP panel is a sorted array, so most reads cost a single comparison, and
//! only the tiny fraction overlapping a site needs CIGAR-resolved base
//! extraction.

pub mod depth;
pub mod output;
pub mod snp;

pub use depth::{ContigDepth, DepthAccum, Window};
pub use snp::{parse_snp_bed, AlleleCounts, SnpAccum, SnpPanel};
