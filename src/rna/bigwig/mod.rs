//! Genome-wide coverage bigWig track generation.
//!
//! Reimplements the nf-core/rnaseq bigWig pipeline (`bedtools genomecov` +
//! UCSC `bedClip` + `bedGraphToBigWig`) in a single streaming pass over the BAM.

pub mod accumulator;
pub mod output;

pub use accumulator::{GenomeCovAccum, GenomeCovResult};
pub use output::write_bigwig_tracks;
