//! BAM/SAM/CRAM readers with a rust-htslib-compatible API.

use std::fs::File;
use std::io::BufRead;
use std::path::{Path, PathBuf};

use super::align_header::Header;
use super::align_record::{alignment_to_record, Record};
use anyhow::{Context, Result};
use noodles::bam as noodles_bam;
use noodles::bgzf;
use noodles::core::Region;
use noodles::cram as noodles_cram;
use noodles::fasta as noodles_fasta;
use noodles::sam as noodles_sam;

fn store_record(record: &mut Record, loaded: Record) -> Result<()> {
    *record = loaded;
    Ok(())
}

/// Fetch target for indexed readers.
#[derive(Debug, Clone, Copy)]
pub enum FetchDefinition {
    /// Unmapped reads stored after mapped reads in the BAM file.
    Unmapped,
}

/// Fetch target accepted by [`IndexedReader::fetch`].
#[derive(Debug, Clone, Copy)]
pub enum FetchTarget {
    /// Reference sequence ID.
    Tid(u32),
    /// Unmapped reads bucket.
    Unmapped,
}

impl From<u32> for FetchTarget {
    fn from(value: u32) -> Self {
        Self::Tid(value)
    }
}

impl From<FetchDefinition> for FetchTarget {
    fn from(value: FetchDefinition) -> Self {
        match value {
            FetchDefinition::Unmapped => Self::Unmapped,
        }
    }
}

/// Trait mirroring `rust_htslib::bam::Read`.
pub trait Read {
    /// Read the next record into `record`.
    fn read(&mut self, record: &mut Record) -> Option<Result<()>>;
}

fn load_fasta_repository(path: &Path) -> Result<noodles_fasta::Repository> {
    let mut reader = noodles_fasta::io::reader::Builder
        .build_from_path(path)
        .with_context(|| format!("failed to open reference FASTA {}", path.display()))?;
    let mut records = Vec::new();
    for result in reader.records() {
        records.push(result.context("failed to read reference FASTA record")?);
    }
    Ok(noodles_fasta::Repository::new(records))
}

fn open_bam_reader(path: &Path) -> Result<noodles_bam::io::Reader<bgzf::io::Reader<File>>> {
    let file =
        File::open(path).with_context(|| format!("failed to open BAM file {}", path.display()))?;
    Ok(noodles_bam::io::Reader::new(file))
}

enum ReaderBackend {
    Bam {
        reader: noodles_bam::io::Reader<bgzf::io::Reader<File>>,
        header: Header,
        scratch: noodles_bam::Record,
    },
    Sam {
        reader: noodles_sam::io::Reader<Box<dyn BufRead>>,
        header: Header,
        scratch: noodles_sam::Record,
    },
    Cram {
        reader: noodles_cram::io::Reader<File>,
        header: Header,
        records: Vec<noodles_sam::alignment::RecordBuf>,
        record_idx: usize,
    },
}

/// Sequential alignment file reader (BAM/SAM/CRAM).
pub struct Reader {
    backend: ReaderBackend,
    source_path: PathBuf,
    cram_reference: Option<PathBuf>,
    decompression_threads: usize,
}

impl Reader {
    /// Open an alignment file, detecting the format from the extension.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        let path_str = path.to_string_lossy();
        let backend = if path_str.ends_with(".sam") {
            Self::open_sam(&path)?
        } else if path_str.ends_with(".cram") {
            Self::open_cram(&path, None)?
        } else {
            Self::open_bam(&path)?
        };
        Ok(Self {
            backend,
            source_path: path,
            cram_reference: None,
            decompression_threads: 0,
        })
    }

    fn open_bam(path: &Path) -> Result<ReaderBackend> {
        let mut reader = open_bam_reader(path)?;
        let header =
            Header::from_noodles(reader.read_header().context("failed to read BAM header")?);
        Ok(ReaderBackend::Bam {
            reader,
            header,
            scratch: noodles_bam::Record::default(),
        })
    }

    fn open_sam(path: &Path) -> Result<ReaderBackend> {
        let mut reader = noodles_sam::io::reader::Builder::default()
            .build_from_path(path)
            .with_context(|| format!("failed to open SAM file {}", path.display()))?;
        let header =
            Header::from_noodles(reader.read_header().context("failed to read SAM header")?);
        Ok(ReaderBackend::Sam {
            reader,
            header,
            scratch: noodles_sam::Record::default(),
        })
    }

    fn open_cram(path: &Path, reference: Option<&Path>) -> Result<ReaderBackend> {
        let repository = if let Some(fasta_path) = reference {
            load_fasta_repository(fasta_path)?
        } else {
            noodles_fasta::Repository::default()
        };
        let file = File::open(path)
            .with_context(|| format!("failed to open CRAM file {}", path.display()))?;
        let mut reader = noodles_cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_reader(file);
        let header =
            Header::from_noodles(reader.read_header().context("failed to read CRAM header")?);
        Ok(ReaderBackend::Cram {
            reader,
            header,
            records: Vec::new(),
            record_idx: 0,
        })
    }

    /// Set the reference FASTA for CRAM decoding.
    pub fn set_reference<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        let path = path.as_ref().to_path_buf();
        self.cram_reference = Some(path.clone());
        if matches!(self.backend, ReaderBackend::Cram { .. }) {
            self.backend = Self::open_cram(&self.source_path, Some(&path))?;
        }
        Ok(())
    }

    /// Enable BGZF decompression threads (BAM only).
    ///
    /// Multithreaded BGZF decompression is not yet wired through noodles; this
    /// stores the requested thread count for API compatibility.
    pub fn set_threads(&mut self, threads: usize) -> Result<()> {
        self.decompression_threads = threads;
        Ok(())
    }

    /// Return the SAM/BAM header.
    pub fn header(&self) -> &Header {
        match &self.backend {
            ReaderBackend::Bam { header, .. }
            | ReaderBackend::Sam { header, .. }
            | ReaderBackend::Cram { header, .. } => header,
        }
    }
}

impl Read for Reader {
    fn read(&mut self, record: &mut Record) -> Option<Result<()>> {
        match &mut self.backend {
            ReaderBackend::Bam {
                reader,
                header,
                scratch,
            } => match reader.read_record(scratch) {
                Ok(0) => None,
                Ok(_) => Some(
                    Record::from_bam(header.noodles_header(), scratch)
                        .and_then(|loaded| store_record(record, loaded)),
                ),
                Err(e) => Some(Err(e.into())),
            },
            ReaderBackend::Sam {
                reader,
                header,
                scratch,
            } => match reader.read_record(scratch) {
                Ok(0) => None,
                Ok(_) => Some(
                    alignment_to_record(header.noodles_header(), scratch)
                        .and_then(|loaded| store_record(record, loaded)),
                ),
                Err(e) => Some(Err(e.into())),
            },
            ReaderBackend::Cram {
                reader,
                header,
                records,
                record_idx,
            } => {
                if records.is_empty() {
                    match reader
                        .records(header.noodles_header())
                        .collect::<Result<Vec<_>, _>>()
                    {
                        Ok(v) => *records = v,
                        Err(e) => return Some(Err(e.into())),
                    }
                }
                if *record_idx >= records.len() {
                    return None;
                }
                let idx = *record_idx;
                *record_idx += 1;
                Some(store_record(
                    record,
                    Record::from_buf(records[idx].clone(), None),
                ))
            }
        }
    }
}

enum IndexedBackend {
    Bam {
        reader: noodles_bam::io::IndexedReader<bgzf::io::Reader<File>>,
        header: Header,
        query_records: Vec<noodles_bam::Record>,
        query_idx: usize,
        unmapped_records: Vec<noodles_bam::Record>,
        unmapped_idx: usize,
    },
    Cram {
        reader: noodles_cram::io::IndexedReader<File>,
        header: Header,
        query_records: Vec<noodles_sam::alignment::RecordBuf>,
        query_idx: usize,
        unmapped_records: Vec<noodles_sam::alignment::RecordBuf>,
        unmapped_idx: usize,
    },
}

/// Indexed alignment file reader (BAM/CRAM with .bai/.csi/.crai).
pub struct IndexedReader {
    backend: IndexedBackend,
    source_path: PathBuf,
    cram_reference: Option<PathBuf>,
}

impl IndexedReader {
    /// Open an indexed alignment file.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        let backend = if path.to_string_lossy().ends_with(".cram") {
            let mut reader = noodles_cram::io::indexed_reader::Builder::default()
                .build_from_path(&path)
                .with_context(|| format!("failed to open indexed CRAM {}", path.display()))?;
            let header = Header::from_noodles(
                reader
                    .read_header()
                    .context("failed to read indexed CRAM header")?,
            );
            IndexedBackend::Cram {
                reader,
                header,
                query_records: Vec::new(),
                query_idx: 0,
                unmapped_records: Vec::new(),
                unmapped_idx: 0,
            }
        } else {
            let mut reader = noodles_bam::io::indexed_reader::Builder::default()
                .build_from_path(&path)
                .with_context(|| format!("failed to open indexed BAM {}", path.display()))?;
            let header = Header::from_noodles(
                reader
                    .read_header()
                    .context("failed to read indexed BAM header")?,
            );
            IndexedBackend::Bam {
                reader,
                header,
                query_records: Vec::new(),
                query_idx: 0,
                unmapped_records: Vec::new(),
                unmapped_idx: 0,
            }
        };
        Ok(Self {
            backend,
            source_path: path,
            cram_reference: None,
        })
    }

    /// Set the reference FASTA for CRAM decoding.
    pub fn set_reference<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        let path = path.as_ref().to_path_buf();
        self.cram_reference = Some(path.clone());
        if let IndexedBackend::Cram { reader, header, .. } = &mut self.backend {
            let repository = load_fasta_repository(&path)?;
            let cram_path = self.source_path.clone();
            *reader = noodles_cram::io::indexed_reader::Builder::default()
                .set_reference_sequence_repository(repository)
                .build_from_path(&cram_path)?;
            *header = Header::from_noodles(reader.read_header()?);
        }
        Ok(())
    }

    /// Enable BGZF decompression threads (BAM only).
    ///
    /// Indexed fetching dominates runtime; this is a no-op for indexed readers.
    pub fn set_threads(&mut self, _threads: usize) -> Result<()> {
        Ok(())
    }

    /// Seek to a reference ID or the unmapped read bucket.
    pub fn fetch<T: Into<FetchTarget>>(&mut self, target: T) -> Result<()> {
        match target.into() {
            FetchTarget::Tid(tid) => self.fetch_tid(tid),
            FetchTarget::Unmapped => self.fetch_unmapped(),
        }
    }

    fn fetch_tid(&mut self, tid: u32) -> Result<()> {
        match &mut self.backend {
            IndexedBackend::Bam {
                reader,
                header,
                query_records,
                query_idx,
                unmapped_records,
                unmapped_idx,
                ..
            } => {
                *unmapped_records = Vec::new();
                *unmapped_idx = 0;
                let region = region_for_tid(header, tid);
                let mut fetched = Vec::new();
                let query = reader
                    .query(header.noodles_header(), &region)
                    .with_context(|| format!("failed to query tid {tid}"))?;
                for result in query.records() {
                    fetched.push(result.context("failed to read queried BAM record")?);
                }
                *query_records = fetched;
                *query_idx = 0;
                Ok(())
            }
            IndexedBackend::Cram {
                reader,
                header,
                query_records,
                query_idx,
                unmapped_records,
                unmapped_idx,
                ..
            } => {
                *unmapped_records = Vec::new();
                *unmapped_idx = 0;
                let region = region_for_tid(header, tid);
                let mut fetched = Vec::new();
                let query = reader.query(header.noodles_header(), &region)?;
                for result in query {
                    fetched.push(result.context("failed to read queried CRAM record")?);
                }
                *query_records = fetched;
                *query_idx = 0;
                Ok(())
            }
        }
    }

    fn fetch_unmapped(&mut self) -> Result<()> {
        match &mut self.backend {
            IndexedBackend::Bam {
                reader,
                query_records,
                query_idx,
                unmapped_records,
                unmapped_idx,
                ..
            } => {
                query_records.clear();
                *query_idx = 0;
                let mut fetched = Vec::new();
                for result in reader
                    .query_unmapped()
                    .context("failed to query unmapped BAM reads")?
                {
                    fetched.push(result.context("failed to read unmapped BAM record")?);
                }
                *unmapped_records = fetched;
                *unmapped_idx = 0;
                Ok(())
            }
            IndexedBackend::Cram {
                reader,
                header,
                query_records,
                query_idx,
                unmapped_records,
                unmapped_idx,
                ..
            } => {
                query_records.clear();
                *query_idx = 0;
                let mut fetched = Vec::new();
                for result in reader
                    .query_unmapped(header.noodles_header())
                    .context("failed to query unmapped CRAM reads")?
                {
                    fetched.push(result.context("failed to read unmapped CRAM record")?);
                }
                *unmapped_records = fetched;
                *unmapped_idx = 0;
                Ok(())
            }
        }
    }
}

impl Read for IndexedReader {
    fn read(&mut self, record: &mut Record) -> Option<Result<()>> {
        match &mut self.backend {
            IndexedBackend::Bam {
                header,
                query_records,
                query_idx,
                unmapped_records,
                unmapped_idx,
                ..
            } => {
                if *unmapped_idx < unmapped_records.len() {
                    let idx = *unmapped_idx;
                    *unmapped_idx += 1;
                    return Some(
                        Record::from_bam(header.noodles_header(), &unmapped_records[idx])
                            .and_then(|loaded| store_record(record, loaded)),
                    );
                }
                if *query_idx >= query_records.len() {
                    return None;
                }
                let idx = *query_idx;
                *query_idx += 1;
                Some(
                    Record::from_bam(header.noodles_header(), &query_records[idx])
                        .and_then(|loaded| store_record(record, loaded)),
                )
            }
            IndexedBackend::Cram {
                query_records,
                query_idx,
                unmapped_records,
                unmapped_idx,
                ..
            } => {
                if *unmapped_idx < unmapped_records.len() {
                    let idx = *unmapped_idx;
                    *unmapped_idx += 1;
                    return Some(store_record(
                        record,
                        Record::from_buf(unmapped_records[idx].clone(), None),
                    ));
                }
                if *query_idx >= query_records.len() {
                    return None;
                }
                let idx = *query_idx;
                *query_idx += 1;
                Some(store_record(
                    record,
                    Record::from_buf(query_records[idx].clone(), None),
                ))
            }
        }
    }
}

fn region_for_tid(header: &Header, tid: u32) -> Region {
    // Query the whole reference sequence by name. Building the `Region` directly
    // (rather than formatting and re-parsing `name:start-end`) avoids mangling
    // reference names that contain `:` or `-` (e.g. HLA contigs) and preserves
    // non-UTF-8 names byte-for-byte.
    Region::new(header.tid2name(tid).to_vec(), ..)
}
