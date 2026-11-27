//! BAM/SAM/CRAM compatibility layer backed by [noodles](https://crates.io/crates/noodles).
//!
//! This module exposes a rust-htslib-shaped API so the rest of RustQC can read
//! alignment files without linking to htslib, while preserving samtools-identical
//! statistics output.

mod align_header;
mod align_record;
mod cigar;
mod io;
mod writer;

pub use align_header::Header;
pub use align_record::{Aux, Record, Seq};
pub use io::{FetchDefinition, IndexedReader, Reader};

/// Trait mirroring `rust_htslib::bam::Read`.
pub use io::Read;

/// Record sub-module mirroring `rust_htslib::bam::record`.
pub mod record {
    pub use super::align_record::{Aux, Record, Seq};
    pub use super::cigar::{Cigar, CigarString, CigarStringView};
}

/// Header sub-module mirroring `rust_htslib::bam::header`.
pub mod header {
    pub use super::align_header::Header;
    pub use super::writer::{HeaderRecord, HeaderView};
}

/// BAM index helpers used by integration tests.
pub mod index {
    pub use super::writer::index::*;
}

/// BAM writer helpers used by integration tests.
pub use writer::{Format, Writer};
