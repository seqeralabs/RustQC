//! RSeQC reimplementations.
//!
//! Standalone reimplementations of selected RSeQC Python quality control scripts
//! for RNA-Seq data analysis.

pub mod accumulators;
pub mod common;
pub mod plots;

pub mod genebody_coverage;
pub mod infer_experiment;
pub mod inner_distance;
pub mod junction_annotation;
pub mod junction_saturation;
pub mod read_distribution;
pub mod read_duplication;
pub mod read_gc;
pub mod tin;

// bam_stat and the samtools writers are read-level and assay-agnostic; they
// now live in `crate::common`. Re-exported so existing
// `crate::rna::rseqc::...` paths and the published 0.2.x library surface
// keep resolving. Drop the shims at 1.0.
pub use crate::common::bam_stat;
pub use crate::common::samtools::{flagstat, idxstats, stats};

#[cfg(test)]
mod compat_tests {
    //! Guards the re-export shims that keep the published 0.2.x paths alive.
    //! These are compile-time assertions; there is nothing to observe at runtime.

    #[test]
    fn moved_modules_are_still_reachable_from_their_old_paths() {
        let _: fn(
            &crate::rna::rseqc::bam_stat::BamStatResult,
            &std::path::Path,
        ) -> anyhow::Result<()> = crate::rna::rseqc::flagstat::write_flagstat;
        let _ = crate::rna::rseqc::accumulators::BamStatAccum::default();
        let _: u16 = crate::rna::bam_flags::BAM_FDUP;
        let _: Option<&crate::rna::preseq::PreseqAccum> = None;
    }
}
