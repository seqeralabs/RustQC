//! Analysis modules shared between the `rna` and `dna` pipelines.
//!
//! Nothing in this module is specific to a library preparation or an assay:
//! BAM flag helpers, the C++ RNG shim used for preseq bootstrap
//! reproducibility, the preseq `lc_extrap` implementation, read-level
//! alignment statistics, and the samtools-compatible output writers.

pub mod bam_flags;
pub mod bam_stat;
pub mod bam_stat_accum;
pub mod coverage;
pub mod cpp_rng;
pub mod preseq;
pub mod samtools;
