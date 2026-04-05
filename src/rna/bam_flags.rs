//! SAM/BAM flag bit constants and auxiliary tag helpers.
//!
//! Centralized definitions for use across all modules. Uses noodles-sam types
//! internally but exposes simple u16 constants for ergonomic bitwise operations.

use noodles_bam as bam;
use noodles_sam::alignment::record::data::field::Tag;

/// Read is paired in sequencing (0x1).
pub const BAM_FPAIRED: u16 = 0x1;
/// Read is mapped in a proper pair (0x2).
pub const BAM_FPROPER_PAIR: u16 = 0x2;
/// Read is unmapped (0x4).
pub const BAM_FUNMAP: u16 = 0x4;
/// Mate is unmapped (0x8).
pub const BAM_FMUNMAP: u16 = 0x8;
/// Read is reverse complemented (0x10).
pub const BAM_FREVERSE: u16 = 0x10;
/// Mate is reverse complemented (0x20).
pub const BAM_FMREVERSE: u16 = 0x20;
/// First in pair (0x40).
pub const BAM_FREAD1: u16 = 0x40;
/// Second in pair (0x80).
pub const BAM_FREAD2: u16 = 0x80;
/// Secondary alignment (0x100).
pub const BAM_FSECONDARY: u16 = 0x100;
/// Failed quality checks (0x200).
pub const BAM_FQCFAIL: u16 = 0x200;
/// PCR or optical duplicate (0x400).
pub const BAM_FDUP: u16 = 0x400;
/// Supplementary alignment (0x800).
pub const BAM_FSUPPLEMENTARY: u16 = 0x800;

/// Extract the flags field from a noodles BAM record as u16.
#[inline]
pub fn flags(record: &bam::Record) -> u16 {
    record.flags().bits()
}

/// Check if record is supplementary alignment.
#[inline]
pub fn is_supplementary(record: &bam::Record) -> bool {
    record.flags().is_supplementary()
}

/// Check if record is paired (segmented in noodles).
#[inline]
pub fn is_paired(record: &bam::Record) -> bool {
    record.flags().is_segmented()
}

/// Extract an integer auxiliary tag from a noodles BAM record.
#[inline]
pub fn get_aux_int(record: &bam::Record, tag_bytes: &[u8]) -> Option<i64> {
    if tag_bytes.len() != 2 {
        return None;
    }
    let tag_arr = [tag_bytes[0], tag_bytes[1]];
    let tag = Tag::try_from(tag_arr).ok()?;
    let data = record.data();
    let result = data.get(&tag)?;
    let field = result.ok()?;
    use noodles_sam::alignment::record::data::field::Value;
    match field {
        Value::Int8(n) => Some(n as i64),
        Value::Int16(n) => Some(n as i64),
        Value::Int32(n) => Some(n as i64),
        Value::UInt8(n) => Some(n as u64 as i64),
        Value::UInt16(n) => Some(n as u64 as i64),
        Value::UInt32(n) => Some(n as u64 as i64),
        _ => None,
    }
}

/// Get 0-based position from noodles record.
#[inline]
pub fn pos_0based(record: &bam::Record) -> i64 {
    record
        .alignment_start()
        .transpose()
        .ok()
        .flatten()
        .map(|p| p.get() as i64 - 1)
        .unwrap_or(-1)
}

/// Get mapping quality as u8.
#[inline]
pub fn mapping_quality(record: &bam::Record) -> u8 {
    record.mapping_quality().map(|q| q.get()).unwrap_or(0)
}

/// Get read name as Vec<u8>.
#[inline]
pub fn read_name(record: &bam::Record) -> Vec<u8> {
    record
        .name()
        .map(|n| n.to_string().into_bytes())
        .unwrap_or_default()
}

/// Get mate position (0-based).
#[inline]
pub fn mate_position_0based(record: &bam::Record) -> i64 {
    record
        .mate_alignment_start()
        .transpose()
        .ok()
        .flatten()
        .map(|p| p.get() as i64 - 1)
        .unwrap_or(-1)
}

/// Get template length (insert size).
#[inline]
pub fn template_length(record: &bam::Record) -> i32 {
    record.template_length()
}

/// Get sequence length.
#[inline]
pub fn sequence_length(record: &bam::Record) -> usize {
    record.sequence().len()
}

/// Check if mate is unmapped.
#[inline]
pub fn is_mate_unmapped(record: &bam::Record) -> bool {
    record.flags().is_mate_unmapped()
}

/// Check if record is reverse complemented.
#[inline]
pub fn is_reverse_complemented(record: &bam::Record) -> bool {
    record.flags().is_reverse_complemented()
}

/// Check if record is first in template.
#[inline]
pub fn is_first_segment(record: &bam::Record) -> bool {
    record.flags().is_first_segment()
}
