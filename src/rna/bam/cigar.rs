//! CIGAR types matching the rust-htslib API surface used by RustQC.

use noodles::sam::alignment::record::cigar::{op::Kind, Op};

/// A CIGAR operation, mirroring `rust_htslib::bam::record::Cigar`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Cigar {
    /// Alignment match (`M`).
    Match(u32),
    /// Insertion (`I`).
    Ins(u32),
    /// Deletion (`D`).
    Del(u32),
    /// Skipped region / intron (`N`).
    RefSkip(u32),
    /// Soft clip (`S`).
    SoftClip(u32),
    /// Hard clip (`H`).
    HardClip(u32),
    /// Padding (`P`).
    Pad(u32),
    /// Sequence match (`=`).
    Equal(u32),
    /// Sequence mismatch (`X`).
    Diff(u32),
}

impl Cigar {
    /// Length of this operation.
    pub fn len(&self) -> u32 {
        match self {
            Self::Match(n)
            | Self::Ins(n)
            | Self::Del(n)
            | Self::RefSkip(n)
            | Self::SoftClip(n)
            | Self::HardClip(n)
            | Self::Pad(n)
            | Self::Equal(n)
            | Self::Diff(n) => *n,
        }
    }

    /// Whether this operation has zero length.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Decode a noodles CIGAR operation into RustQC's CIGAR enum.
pub fn decode_op_from_op(op: Op) -> Cigar {
    match op.kind() {
        Kind::Match => Cigar::Match(op.len() as u32),
        Kind::Insertion => Cigar::Ins(op.len() as u32),
        Kind::Deletion => Cigar::Del(op.len() as u32),
        Kind::Skip => Cigar::RefSkip(op.len() as u32),
        Kind::SoftClip => Cigar::SoftClip(op.len() as u32),
        Kind::HardClip => Cigar::HardClip(op.len() as u32),
        Kind::Pad => Cigar::Pad(op.len() as u32),
        Kind::SequenceMatch => Cigar::Equal(op.len() as u32),
        Kind::SequenceMismatch => Cigar::Diff(op.len() as u32),
    }
}

/// Decode a fallible noodles CIGAR iterator item.
pub fn decode_op(op: Result<Op, std::io::Error>) -> Result<Cigar, std::io::Error> {
    op.map(decode_op_from_op)
}

/// Compute the 1-based exclusive reference end position using htslib semantics
/// (M, D, N, =, X consume reference; I, S, H, P do not).
pub fn reference_end(start: u64, ops: &[Cigar]) -> u64 {
    let mut end = start;
    for op in ops {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                end += u64::from(*len);
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                end += u64::from(*len);
            }
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    end
}

/// Owned CIGAR string (used in unit tests).
#[derive(Debug, Clone, Default)]
pub struct CigarString(pub Vec<Cigar>);

/// View over a decoded CIGAR with htslib-compatible helpers.
#[derive(Debug, Clone)]
pub struct CigarStringView {
    ops: Vec<Cigar>,
    start: u64,
}

impl CigarStringView {
    /// Build a view from decoded operations and a 0-based reference start.
    pub(crate) fn from_ops(ops: Vec<Cigar>, start: u64) -> Self {
        Self { ops, start }
    }

    /// Build a view from an owned CIGAR string and reference start (0-based).
    pub fn new(cigar: CigarString, start: u64) -> Self {
        Self {
            ops: cigar.0,
            start,
        }
    }

    /// Iterate CIGAR operations.
    pub fn iter(&self) -> impl Iterator<Item = &Cigar> {
        self.ops.iter()
    }

    /// Reference end position (htslib `end_pos` semantics).
    pub fn end_pos(&self) -> u64 {
        reference_end(self.start, &self.ops)
    }

    /// Slice view of decoded operations (rust-htslib `CigarStringView::as_ref`).
    #[allow(clippy::wrong_self_convention, clippy::should_implement_trait)]
    pub fn as_ref(&self) -> &[Cigar] {
        &self.ops
    }
}
