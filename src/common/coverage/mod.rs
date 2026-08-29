//! Coverage tracks.
//!
//! Produces the per-base coverage that pipelines currently get by running
//! `bedtools genomecov` and then converting its bedGraph to bigWig. Both come
//! out of the alignment pass that is already happening.

pub mod bedgraph;

/// bigWig writing, available when built with the `bigwig` feature.
#[cfg(feature = "bigwig")]
pub mod bigwig;
