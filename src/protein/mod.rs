//! Protein quality control modules.
//!
//! Three modes that share an entry point and nothing else, because the
//! questions they answer take different inputs entirely:
//!
//! - [`sequence`] reads a protein FASTA and asks whether the proteome is
//!   well-formed;
//! - `coding` reads an alignment and an annotation and asks where reads fall
//!   relative to coding sequence;
//! - `spectra` reads mzML and asks whether a mass spectrometry run is healthy.

pub mod sequence;
