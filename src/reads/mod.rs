//! Raw read quality control, before alignment.
//!
//! This sits apart from `rna`, `dna` and `protein` because the question is
//! different: those pipelines ask what the alignment says, this one asks
//! whether the reads coming off the instrument are usable at all. It applies
//! whatever the library was for.
//!
//! Outputs are compatible with the two tools people already parse: `seqkit
//! stats -a -T` and FastQC's `fastqc_data.txt`.

pub mod fastq;
pub mod metrics;
pub mod output;
