//! Protein FASTA quality control.
//!
//! Reproduces `seqkit stats -a -T` and adds what seqkit does not report:
//! amino acid composition, internal stop codons, non-standard residues and
//! duplicate sequences.

pub mod defects;
pub mod output;
pub mod stats;

use std::io::BufRead;
use std::path::Path;

use anyhow::{Context, Result};

/// One sequence read from a FASTA file.
#[derive(Debug, Clone)]
pub struct Record {
    /// Everything after `>` on the header line, up to the first whitespace.
    pub id: String,
    /// The residues, uppercased, with whitespace removed.
    pub residues: Vec<u8>,
}

impl Record {
    /// Number of residues.
    pub fn len(&self) -> usize {
        self.residues.len()
    }

    /// Whether the record holds no residues.
    pub fn is_empty(&self) -> bool {
        self.residues.is_empty()
    }
}

/// Read every record from a FASTA file, transparently handling gzip.
///
/// Residues are uppercased on the way in, so downstream code never has to
/// think about case. Gap characters are kept, because seqkit counts them.
pub fn read_fasta(path: &Path) -> Result<Vec<Record>> {
    let reader = crate::io::open_reader(path)
        .with_context(|| format!("Failed to open FASTA: {}", path.display()))?;

    let mut records: Vec<Record> = Vec::new();
    let mut current: Option<Record> = None;

    for line in reader.lines() {
        let line = line.with_context(|| format!("Failed to read FASTA: {}", path.display()))?;
        let line = line.trim_end();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            if let Some(record) = current.take() {
                records.push(record);
            }
            let id = header.split_whitespace().next().unwrap_or("").to_string();
            current = Some(Record {
                id,
                residues: Vec::new(),
            });
        } else {
            match current.as_mut() {
                Some(record) => record.residues.extend(
                    line.bytes()
                        .filter(|b| !b.is_ascii_whitespace())
                        .map(|b| b.to_ascii_uppercase()),
                ),
                None => anyhow::bail!("{}: sequence data before any header line", path.display()),
            }
        }
    }
    if let Some(record) = current {
        records.push(record);
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write(name: &str, contents: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join("rustqc-protein-tests");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(name);
        std::fs::write(&path, contents).unwrap();
        path
    }

    #[test]
    fn records_are_split_on_header_lines() {
        let path = write("two.fa", ">a desc here\nMKT\nAYI\n>b\nMSEQ\n");
        let records = read_fasta(&path).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "a", "the id stops at the first space");
        assert_eq!(records[0].residues, b"MKTAYI", "wrapped lines are joined");
        assert_eq!(records[1].residues, b"MSEQ");
    }

    #[test]
    fn residues_are_uppercased_and_stripped_of_whitespace() {
        let path = write("case.fa", ">a\nmk t\tayi\n");
        let records = read_fasta(&path).unwrap();
        assert_eq!(records[0].residues, b"MKTAYI");
    }

    #[test]
    fn blank_lines_are_ignored() {
        let path = write("blank.fa", ">a\nMKT\n\n\nAYI\n\n");
        let records = read_fasta(&path).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].residues, b"MKTAYI");
    }

    #[test]
    fn sequence_data_before_a_header_is_an_error() {
        let path = write("headless.fa", "MKTAYI\n>a\nMSEQ\n");
        assert!(
            read_fasta(&path).is_err(),
            "data before a header must be rejected rather than silently dropped"
        );
    }
}
