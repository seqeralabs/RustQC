//! SAM/BAM flag bit constants.
//!
//! Centralised definitions so every module uses the same constants.

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

// ============================================================
// Auxiliary tag helpers
// ============================================================

/// Extract an integer auxiliary tag from a BAM record.
///
/// Handles all integer Aux variants (U8, U16, U32, I8, I16, I32)
/// and returns the value as `i64`. Returns `None` if the tag is
/// absent or has a non-integer type.
pub fn get_aux_int(record: &rust_htslib::bam::Record, tag: &[u8]) -> Option<i64> {
    use rust_htslib::bam::record::Aux;
    match record.aux(tag) {
        Ok(Aux::U8(v)) => Some(v as i64),
        Ok(Aux::U16(v)) => Some(v as i64),
        Ok(Aux::U32(v)) => Some(v as i64),
        Ok(Aux::I8(v)) => Some(v as i64),
        Ok(Aux::I16(v)) => Some(v as i64),
        Ok(Aux::I32(v)) => Some(v as i64),
        _ => None,
    }
}
