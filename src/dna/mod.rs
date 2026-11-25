//! DNA quality control and analysis modules.
//!
//! Contains the depth of coverage engine and the mosdepth-compatible outputs
//! built on top of it. Read-level statistics, the samtools-compatible writers
//! and preseq are shared with the RNA pipeline and live in [`crate::common`].

pub mod depth;
pub mod gc_bias;
pub mod hs_metrics;
pub mod insert_size;
pub mod intervals;
pub mod mosdepth;
pub mod wgs_metrics;
