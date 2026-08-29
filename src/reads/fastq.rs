//! FASTQ reading.
//!
//! Deliberately minimal: the records stream through an accumulator and are
//! never all held at once, because a FASTQ can be tens of gigabytes.

use std::io::BufRead;
use std::path::Path;

use anyhow::{bail, Context, Result};

/// One FASTQ record, borrowed from the reader's buffers.
#[derive(Debug, Default)]
pub struct Record {
    /// The identifier, without the leading `@`.
    pub id: String,
    /// The bases, uppercased.
    pub sequence: Vec<u8>,
    /// Phred scores, already decoded from the ASCII offset.
    pub quality: Vec<u8>,
}

/// The Phred offset of Sanger and Illumina 1.8+ encoding, which is all modern
/// data uses.
pub const PHRED_OFFSET: u8 = 33;

/// Stream every record through `consume`.
///
/// Records are handed over one at a time and the buffers are reused, so peak
/// memory does not grow with the file.
pub fn for_each_record<F>(path: &Path, mut consume: F) -> Result<()>
where
    F: FnMut(&Record),
{
    let mut reader = crate::io::open_reader(path)
        .with_context(|| format!("Failed to open FASTQ: {}", path.display()))?;

    let mut header = String::new();
    let mut sequence = String::new();
    let mut plus = String::new();
    let mut quality = String::new();
    let mut record = Record::default();
    let mut number = 0u64;

    loop {
        header.clear();
        if reader.read_line(&mut header)? == 0 {
            break;
        }
        if header.trim().is_empty() {
            continue;
        }
        number += 1;

        sequence.clear();
        plus.clear();
        quality.clear();
        reader.read_line(&mut sequence)?;
        reader.read_line(&mut plus)?;
        reader.read_line(&mut quality)?;

        let header_line = header.trim_end();
        if !header_line.starts_with('@') {
            bail!(
                "{}: record {number} does not start with '@'",
                path.display()
            );
        }
        let sequence_line = sequence.trim_end();
        let quality_line = quality.trim_end();
        if sequence_line.len() != quality_line.len() {
            bail!(
                "{}: record {number} has {} bases but {} quality scores",
                path.display(),
                sequence_line.len(),
                quality_line.len()
            );
        }

        record.id.clear();
        record.id.push_str(
            header_line[1..]
                .split_whitespace()
                .next()
                .unwrap_or_default(),
        );
        record.sequence.clear();
        record
            .sequence
            .extend(sequence_line.bytes().map(|b| b.to_ascii_uppercase()));
        record.quality.clear();
        record
            .quality
            .extend(quality_line.bytes().map(|b| b.saturating_sub(PHRED_OFFSET)));

        consume(&record);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write(name: &str, contents: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join("rustqc-reads-tests");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(name);
        std::fs::write(&path, contents).unwrap();
        path
    }

    #[test]
    fn records_are_read_in_order() {
        let path = write("two.fq", "@a desc\nACGT\n+\nIIII\n@b\nTTTT\n+\n!!!!\n");
        let mut ids = Vec::new();
        let mut qualities = Vec::new();
        for_each_record(&path, |r| {
            ids.push(r.id.clone());
            qualities.push(r.quality.clone());
        })
        .unwrap();
        assert_eq!(ids, vec!["a", "b"], "the id stops at the first space");
        assert_eq!(qualities[0], vec![40, 40, 40, 40], "'I' is Phred 40");
        assert_eq!(qualities[1], vec![0, 0, 0, 0], "'!' is Phred 0");
    }

    #[test]
    fn bases_are_uppercased() {
        let path = write("case.fq", "@a\nacgt\n+\nIIII\n");
        let mut seen = Vec::new();
        for_each_record(&path, |r| seen.push(r.sequence.clone())).unwrap();
        assert_eq!(seen[0], b"ACGT");
    }

    #[test]
    fn a_length_mismatch_is_an_error_rather_than_silent_truncation() {
        let path = write("ragged.fq", "@a\nACGT\n+\nII\n");
        let result = for_each_record(&path, |_| {});
        assert!(
            result.is_err(),
            "a record with fewer quality scores than bases must be rejected"
        );
    }

    #[test]
    fn a_missing_at_sign_is_an_error() {
        let path = write("noat.fq", "a\nACGT\n+\nIIII\n");
        assert!(for_each_record(&path, |_| {}).is_err());
    }
}
