//! Test helpers for writing BAM files (integration tests).

use std::fs::File;
use std::io;
use std::num::NonZero;
use std::path::{Path, PathBuf};

use noodles::bam as noodles_bam;
use noodles::sam as noodles_sam;
use noodles::sam::header::record::value::map::header::{tag, Version};
use noodles::sam::header::record::value::map::{self, ReferenceSequence};
use noodles::sam::header::record::value::Map;

use noodles::sam::alignment::io::Write as AlignmentWrite;

use super::align_header::Header;
use super::align_record::{record_buf_for_writer, Record};

fn parse_sam_version(version: &str) -> Version {
    let mut parts = version.split('.');
    let major = parts.next().and_then(|s| s.parse().ok()).unwrap_or(1);
    let minor = parts.next().and_then(|s| s.parse().ok()).unwrap_or(6);
    Version::new(major, minor)
}

/// BAM output format.
#[derive(Debug, Clone, Copy)]
pub enum Format {
    /// Binary BAM.
    Bam,
}

/// SAM header builder record (rust-htslib-compatible).
#[derive(Debug, Clone)]
pub struct HeaderRecord {
    tag: [u8; 2],
    fields: Vec<(Vec<u8>, String)>,
}

impl HeaderRecord {
    /// Create a new header line builder (`@HD`, `@SQ`, etc.).
    pub fn new(tag: &[u8; 2]) -> Self {
        Self {
            tag: *tag,
            fields: Vec::new(),
        }
    }

    /// Add a `TAG:value` field and return self for chaining.
    pub fn push_tag(mut self, tag: &[u8], value: impl Into<String>) -> Self {
        self.fields.push((tag.to_vec(), value.into()));
        self
    }
}

/// Header view used when parsing SAM lines into records.
#[derive(Debug, Clone)]
pub struct HeaderView {
    header: noodles_sam::Header,
}

impl HeaderView {
    /// Build a header view from a RustQC header wrapper.
    pub fn from_header(header: &Header) -> Self {
        Self {
            header: header.noodles_header().clone(),
        }
    }
}

/// BAM writer for test fixtures.
pub struct Writer {
    inner: noodles_bam::io::Writer<noodles::bgzf::io::Writer<File>>,
    header: noodles_sam::Header,
}

impl Writer {
    /// Create a BAM writer at `path`.
    pub fn from_path(path: &Path, header: &Header, _format: Format) -> io::Result<Self> {
        let file = File::create(path)?;
        let mut inner = noodles_bam::io::Writer::new(file);
        inner.write_header(header.noodles_header())?;
        Ok(Self {
            inner,
            header: header.noodles_header().clone(),
        })
    }

    /// Write a BAM record.
    pub fn write(&mut self, record: &Record) -> io::Result<()> {
        self.inner
            .write_alignment_record(&self.header, record_buf_for_writer(record))
    }
}

impl Header {
    /// Append a `@`-prefixed header record.
    pub fn push_record(&mut self, record: HeaderRecord) {
        let mut inner = self.inner.clone();
        match &record.tag {
            b"HD" => {
                let mut version = "1.6".to_string();
                let mut sort_order = None;
                for (tag, value) in &record.fields {
                    match &tag[..] {
                        b"VN" => version = value.clone(),
                        b"SO" => sort_order = Some(value.clone()),
                        _ => {}
                    }
                }
                let version = parse_sam_version(&version);
                let mut hd_builder = Map::<map::Header>::builder().set_version(version);
                if let Some(so) = sort_order {
                    hd_builder = hd_builder.insert(tag::SORT_ORDER, so);
                }
                let hd = hd_builder.build().expect("valid SAM header map");
                inner = noodles_sam::Header::builder()
                    .set_header(hd)
                    .set_reference_sequences(inner.reference_sequences().clone())
                    .build();
            }
            b"SQ" => {
                let mut name = String::new();
                let mut len = 0u64;
                for (tag, value) in &record.fields {
                    match &tag[..] {
                        b"SN" => name = value.clone(),
                        b"LN" => len = value.parse().unwrap_or(0),
                        _ => {}
                    }
                }
                if !name.is_empty() {
                    let mut builder = noodles_sam::Header::builder();
                    if let Some(hd) = inner.header().cloned() {
                        builder = builder.set_header(hd);
                    }
                    for (existing_name, map) in inner.reference_sequences() {
                        builder = builder.add_reference_sequence(
                            String::from_utf8_lossy(existing_name.as_ref()).as_ref(),
                            map.clone(),
                        );
                    }
                    builder = builder.add_reference_sequence(
                        name.as_str(),
                        Map::<ReferenceSequence>::new(NonZero::new(len.max(1) as usize).unwrap()),
                    );
                    inner = builder.build();
                }
            }
            _ => {}
        }
        *self = Header::from_noodles(inner);
    }
}

impl Record {
    /// Parse a SAM line into a BAM record.
    pub fn from_sam(header: &HeaderView, line: &[u8]) -> io::Result<Self> {
        let sam_record = noodles_sam::Record::try_from(line)?;
        super::align_record::alignment_to_record(&header.header, &sam_record)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }
}

/// BAM index construction (test fixtures).
pub mod index {
    use super::*;

    /// Index type (BAI only in tests).
    #[derive(Debug, Clone, Copy)]
    pub enum Type {
        /// BAI index.
        Bai,
    }

    /// Build a BAI index beside a coordinate-sorted BAM file.
    pub fn build(
        bam_path: &Path,
        index_path: Option<&Path>,
        _kind: Type,
        _min_shift: u32,
    ) -> io::Result<()> {
        let index = noodles_bam::fs::index(bam_path)?;
        let index_path = index_path.map_or_else(
            || PathBuf::from(format!("{}.bai", bam_path.display())),
            |p| p.to_path_buf(),
        );
        noodles_bam::bai::fs::write(&index_path, &index)
    }
}
