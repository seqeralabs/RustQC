//! BAM record wrapper with a rust-htslib-compatible API backed by noodles.

use std::io;

use anyhow::{Context, Result};
use noodles::bam as noodles_bam;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::Record as AlignmentRecord;
use noodles::sam::alignment::RecordBuf;

use super::cigar::{decode_op, decode_op_from_op, Cigar, CigarStringView};

/// SAM/BAM auxiliary tag value (integer variants used by RustQC).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Aux {
    U8(u8),
    U16(u16),
    U32(u32),
    I8(i8),
    I16(i16),
    I32(i32),
}

/// Sequence view with htslib-compatible base encoding.
pub struct Seq<'a> {
    seq_packed: &'a [u8],
    seq_len: usize,
}

impl<'a> Seq<'a> {
    pub fn len(&self) -> usize {
        self.seq_len
    }

    pub fn is_empty(&self) -> bool {
        self.seq_len == 0
    }

    /// BAM 4-bit encoded sequence bytes for samtools-compatible CHK checksums.
    pub fn encoded_bytes(&self) -> &[u8] {
        let nbytes = self.seq_len.div_ceil(2);
        &self.seq_packed[..nbytes.min(self.seq_packed.len())]
    }

    /// BAM 4-bit encoded base at `i` (A=1, C=2, G=4, T=8, N=15).
    pub fn encoded_base(&self, i: usize) -> u8 {
        if i >= self.seq_len {
            return 15;
        }
        let byte = self.seq_packed[i / 2];
        if i.is_multiple_of(2) {
            byte >> 4
        } else {
            byte & 0x0f
        }
    }

    pub fn as_bytes(&self) -> Vec<u8> {
        (0..self.seq_len)
            .map(|i| b"=ACMGRSVTWYHKDBN"[self.encoded_base(i) as usize])
            .collect()
    }
}

/// A BAM alignment record with rust-htslib-compatible accessors.
#[derive(Debug, Clone)]
pub struct Record {
    inner: RecordBuf,
    seq_packed: Vec<u8>,
    qual_cache: Vec<u8>,
    cigar_cache: Vec<Cigar>,
}

impl Default for Record {
    fn default() -> Self {
        Self::new()
    }
}

impl Record {
    pub fn new() -> Self {
        Self {
            inner: RecordBuf::default(),
            seq_packed: Vec::new(),
            qual_cache: Vec::new(),
            cigar_cache: Vec::new(),
        }
    }

    pub fn set_buf(&mut self, buf: RecordBuf, seq_packed: Vec<u8>) {
        self.inner = buf;
        self.seq_packed = seq_packed;
        self.refresh_caches();
    }

    fn refresh_caches(&mut self) {
        self.qual_cache = self.inner.quality_scores().as_ref().to_vec();
        self.cigar_cache = self
            .inner
            .cigar()
            .as_ref()
            .iter()
            .map(|&op| decode_op_from_op(op))
            .collect();
    }

    pub(crate) fn from_bam(
        header: &noodles::sam::Header,
        bam: &noodles_bam::Record,
    ) -> Result<Self> {
        let mut inner = RecordBuf::default();
        inner
            .try_clone_from_alignment_record(header, bam)
            .context("failed to convert BAM record")?;
        let cigar_cache = bam
            .cigar()
            .iter()
            .map(decode_op)
            .collect::<Result<Vec<_>, io::Error>>()
            .context("failed to decode BAM CIGAR")?;
        Ok(Self {
            inner,
            seq_packed: bam.sequence().as_bytes().to_vec(),
            qual_cache: bam.quality_scores().as_ref().to_vec(),
            cigar_cache,
        })
    }

    pub(crate) fn from_buf(buf: RecordBuf, seq_packed: Option<Vec<u8>>) -> Self {
        let seq_packed = seq_packed.unwrap_or_else(|| pack_ascii_sequence(buf.sequence().as_ref()));
        let mut record = Self {
            inner: buf,
            seq_packed,
            qual_cache: Vec::new(),
            cigar_cache: Vec::new(),
        };
        record.refresh_caches();
        record
    }

    pub fn flags(&self) -> u16 {
        self.inner.flags().bits()
    }

    pub fn tid(&self) -> i32 {
        self.inner
            .reference_sequence_id()
            .map(|id| id as i32)
            .unwrap_or(-1)
    }

    pub fn pos(&self) -> i64 {
        self.inner
            .alignment_start()
            .map(|p| i64::try_from(usize::from(p) - 1).unwrap_or(-1))
            .unwrap_or(-1)
    }

    pub fn mtid(&self) -> i32 {
        self.inner
            .mate_reference_sequence_id()
            .map(|id| id as i32)
            .unwrap_or(-1)
    }

    pub fn mpos(&self) -> i64 {
        self.inner
            .mate_alignment_start()
            .map(|p| i64::try_from(usize::from(p) - 1).unwrap_or(-1))
            .unwrap_or(-1)
    }

    pub fn mapq(&self) -> u8 {
        self.inner
            .mapping_quality()
            .map(|mq| mq.get())
            .unwrap_or(255)
    }

    pub fn qname(&self) -> &[u8] {
        self.inner.name().map(|n| n.as_ref()).unwrap_or(b"*")
    }

    pub fn seq_len(&self) -> i32 {
        i32::try_from(self.inner.sequence().len()).unwrap_or(0)
    }

    pub fn insert_size(&self) -> i32 {
        self.inner.template_length()
    }

    pub fn cigar(&self) -> CigarStringView {
        CigarStringView::from_ops(self.cigar_cache.clone(), self.pos().max(0) as u64)
    }

    pub fn seq(&self) -> Seq<'_> {
        Seq {
            seq_packed: &self.seq_packed,
            seq_len: self.inner.sequence().len(),
        }
    }

    pub fn qual(&self) -> &[u8] {
        &self.qual_cache
    }

    pub fn aux(&self, tag: &[u8]) -> Result<Aux, io::Error> {
        if tag.len() != 2 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "aux tag must be 2 bytes",
            ));
        }
        let tag = Tag::new(tag[0], tag[1]);
        match self.inner.data().get(&tag) {
            Some(Value::Int8(v)) => Ok(Aux::I8(*v)),
            Some(Value::Int16(v)) => Ok(Aux::I16(*v)),
            Some(Value::Int32(v)) => Ok(Aux::I32(*v)),
            Some(Value::UInt8(v)) => Ok(Aux::U8(*v)),
            Some(Value::UInt16(v)) => Ok(Aux::U16(*v)),
            Some(Value::UInt32(v)) => Ok(Aux::U32(*v)),
            Some(_) => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "non-integer aux tag",
            )),
            None => Err(io::Error::new(io::ErrorKind::NotFound, "aux tag not found")),
        }
    }

    pub fn is_paired(&self) -> bool {
        self.inner.flags().is_segmented()
    }

    pub fn is_reverse(&self) -> bool {
        self.inner.flags().is_reverse_complemented()
    }

    pub fn is_first_in_template(&self) -> bool {
        self.inner.flags().is_first_segment()
    }

    pub fn is_unmapped(&self) -> bool {
        self.inner.flags().is_unmapped()
    }

    pub fn is_mate_unmapped(&self) -> bool {
        self.inner.flags().is_mate_unmapped()
    }

    pub fn is_secondary(&self) -> bool {
        self.inner.flags().is_secondary()
    }

    pub fn is_supplementary(&self) -> bool {
        self.inner.flags().is_supplementary()
    }

    pub fn is_quality_check_failed(&self) -> bool {
        self.inner.flags().is_qc_fail()
    }
}
fn pack_ascii_sequence(seq: &[u8]) -> Vec<u8> {
    // Full htslib `seq_nt16_table` mapping (the inverse of the `=ACMGRSVTWYHKDBN`
    // decode table), so IUPAC ambiguity codes round-trip through SAM/CRAM packing
    // and produce samtools-identical CHK checksums. Unknown bytes fall back to N.
    const TABLE: [u8; 256] = {
        let mut t = [15u8; 256];
        let symbols = b"=ACMGRSVTWYHKDBN";
        let mut code = 0usize;
        while code < symbols.len() {
            let upper = symbols[code];
            t[upper as usize] = code as u8;
            // Map the lowercase variant of each letter to the same code.
            if upper.is_ascii_uppercase() {
                t[upper.to_ascii_lowercase() as usize] = code as u8;
            }
            code += 1;
        }
        t
    };
    let mut out = Vec::with_capacity(seq.len().div_ceil(2));
    for chunk in seq.chunks(2) {
        let l = TABLE[chunk[0] as usize];
        let r = chunk.get(1).map(|b| TABLE[*b as usize]).unwrap_or(0);
        out.push((l << 4) | r);
    }
    out
}

pub(crate) fn alignment_to_record(
    header: &noodles::sam::Header,
    alignment: &dyn AlignmentRecord,
) -> Result<Record> {
    let mut inner = RecordBuf::default();
    inner
        .try_clone_from_alignment_record(header, alignment)
        .context("failed to clone alignment record")?;
    Ok(Record::from_buf(inner, None))
}

/// Export for writer tests that need BAM bytes from a SAM line.
pub(crate) fn record_buf_for_writer(record: &Record) -> &RecordBuf {
    &record.inner
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn packs_all_iupac_codes_losslessly() {
        // The full set of IUPAC symbols, in htslib 4-bit code order, must pack
        // and decode back unchanged (the inverse of the `=ACMGRSVTWYHKDBN` table).
        let seq = b"=ACMGRSVTWYHKDBN".to_vec();
        let packed = pack_ascii_sequence(&seq);
        let view = Seq {
            seq_packed: &packed,
            seq_len: seq.len(),
        };
        assert_eq!(view.as_bytes(), seq);
        for code in 0u8..16 {
            assert_eq!(view.encoded_base(code as usize), code);
        }
    }

    #[test]
    fn packs_lowercase_like_uppercase() {
        assert_eq!(
            pack_ascii_sequence(b"acgtnmrn"),
            pack_ascii_sequence(b"ACGTNMRN")
        );
    }

    #[test]
    fn packs_odd_length_trailing_nibble_is_zero() {
        // Odd-length sequences leave the final low nibble as 0 (`=`), matching
        // the BAM 4-bit packing layout.
        let packed = pack_ascii_sequence(b"A");
        assert_eq!(packed, vec![0b0001_0000]);
    }
}
