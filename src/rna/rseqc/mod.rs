//! RSeQC reimplementations.
//!
//! Standalone reimplementations of selected RSeQC Python quality control scripts
//! for RNA-Seq data analysis.

pub mod accumulators;
pub mod common;
pub mod plots;

pub mod bam_stat;
pub mod flagstat;
pub mod idxstats;
pub mod infer_experiment;
pub mod inner_distance;
pub mod junction_annotation;
pub mod junction_saturation;
pub mod read_distribution;
pub mod read_duplication;
pub mod stats;
pub mod tin;
