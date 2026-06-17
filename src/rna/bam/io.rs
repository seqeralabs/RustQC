//! BAM/SAM/CRAM readers with a rust-htslib-compatible API.

use std::fs::File;
use std::io::BufRead;
use std::path::{Path, PathBuf};

use super::align_header::Header;
use super::align_record::{alignment_to_record, Record};
use anyhow::{Context, Result};
use noodles::bam as noodles_bam;
use noodles::bgzf;
use noodles::bgzf::io::Seek as _;
use noodles::core::region::Interval;
use noodles::core::Region;
use noodles::cram as noodles_cram;
use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
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

/// Owned iterator over decoded CRAM records. Boxed so the `query` and
/// `query_unmapped` iterators (which have distinct concrete types) can share one
/// field on the indexed cursor.
type CramRecordIter<'a> =
    Box<dyn Iterator<Item = std::io::Result<noodles_sam::alignment::RecordBuf>> + 'a>;

/// Self-referential holder for a sequential CRAM reader plus its streaming record
/// iterator.
///
/// noodles' `Records` iterator borrows the reader and header, so it cannot sit
/// next to them in a plain struct. ouroboros makes that borrow sound, letting us
/// stream one container at a time instead of collecting every record up front.
#[ouroboros::self_referencing]
struct CramReaderCursor {
    reader: noodles_cram::io::Reader<File>,
    header: noodles_sam::Header,
    #[borrows(mut reader, header)]
    #[not_covariant]
    iter: noodles_cram::io::reader::Records<'this, 'this, File>,
}

/// Self-referential holder for an indexed CRAM reader plus the iterator from its
/// most recent `fetch`. Rebuilt on each fetch via `into_heads`; before the first
/// fetch the iterator is empty.
#[ouroboros::self_referencing]
struct CramIndexedCursor {
    reader: noodles_cram::io::IndexedReader<File>,
    header: noodles_sam::Header,
    #[borrows(mut reader, header)]
    #[not_covariant]
    iter: CramRecordIter<'this>,
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
        cursor: CramReaderCursor,
        header: Header,
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
        let sam_header = reader.read_header().context("failed to read CRAM header")?;
        let header = Header::from_noodles(sam_header.clone());
        let cursor = CramReaderCursorBuilder {
            reader,
            header: sam_header,
            iter_builder: |reader, header| reader.records(header),
        }
        .build();
        Ok(ReaderBackend::Cram { cursor, header })
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
            ReaderBackend::Cram { cursor, .. } => match cursor.with_iter_mut(|iter| iter.next()) {
                Some(Ok(buf)) => Some(store_record(record, Record::from_buf(buf, None))),
                Some(Err(e)) => Some(Err(e.into())),
                None => None,
            },
        }
    }
}

/// Lazy streaming cursor for indexed BAM reads.
///
/// Mirrors the chunk state machine in `noodles_csi::io::Query` so records are
/// decoded one at a time directly from the BGZF stream, instead of buffering an
/// entire reference sequence (or all unmapped reads) in memory. This matters
/// because dupRadar fetches one chromosome per worker thread concurrently, so
/// the old buffering cost scaled as reference size times thread count.
enum BamCursor {
    /// No active fetch, or the current fetch has been exhausted.
    Done,
    /// Streaming the mapped records of a single reference via its index chunks.
    Region {
        chunks: std::vec::IntoIter<Chunk>,
        /// End virtual position of the chunk currently being read; `None` means
        /// the next chunk still needs to be seeked to.
        chunk_end: Option<bgzf::VirtualPosition>,
        reference_sequence_id: usize,
    },
    /// Streaming unmapped records from a seeked position to end of file.
    Unmapped,
}

enum IndexedBackend {
    Bam {
        reader: noodles_bam::io::IndexedReader<bgzf::io::Reader<File>>,
        header: Header,
        scratch: noodles_bam::Record,
        cursor: BamCursor,
    },
    Cram {
        /// `None` only transiently while a fetch rebuilds the cursor.
        cursor: Option<CramIndexedCursor>,
        header: Header,
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
            let sam_header = reader
                .read_header()
                .context("failed to read indexed CRAM header")?;
            let header = Header::from_noodles(sam_header.clone());
            let cursor = CramIndexedCursorBuilder {
                reader,
                header: sam_header,
                iter_builder: |_reader, _header| Box::new(std::iter::empty()) as CramRecordIter,
            }
            .build();
            IndexedBackend::Cram {
                cursor: Some(cursor),
                header,
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
                scratch: noodles_bam::Record::default(),
                cursor: BamCursor::Done,
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
        let cram_path = self.source_path.clone();
        if let IndexedBackend::Cram { cursor, header } = &mut self.backend {
            let repository = load_fasta_repository(&path)?;
            let mut reader = noodles_cram::io::indexed_reader::Builder::default()
                .set_reference_sequence_repository(repository)
                .build_from_path(&cram_path)?;
            let sam_header = reader.read_header()?;
            *header = Header::from_noodles(sam_header.clone());
            *cursor = Some(
                CramIndexedCursorBuilder {
                    reader,
                    header: sam_header,
                    iter_builder: |_reader, _header| Box::new(std::iter::empty()) as CramRecordIter,
                }
                .build(),
            );
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
            IndexedBackend::Bam { reader, cursor, .. } => {
                // Resolve the chunks covering the whole reference straight from
                // the index and stream them lazily, rather than buffering every
                // record on the reference in memory.
                let chunks = reader
                    .index()
                    .query(tid as usize, Interval::from(..))
                    .with_context(|| format!("failed to query index for tid {tid}"))?;
                *cursor = BamCursor::Region {
                    chunks: chunks.into_iter(),
                    chunk_end: None,
                    reference_sequence_id: tid as usize,
                };
                Ok(())
            }
            IndexedBackend::Cram { cursor, header } => {
                // Rebuild the cursor with a fresh region query. noodles' query
                // iterator borrows the reader, so we recover the owned reader via
                // `into_heads` and stream container-by-container from there.
                let region = region_for_tid(header, tid);
                let heads = cursor.take().expect("CRAM cursor present").into_heads();
                *cursor = Some(
                    CramIndexedCursorTryBuilder {
                        reader: heads.reader,
                        header: heads.header,
                        iter_builder: |reader, header| -> Result<CramRecordIter> {
                            let query = reader
                                .query(header, &region)
                                .with_context(|| format!("failed to query CRAM tid {tid}"))?;
                            Ok(Box::new(query))
                        },
                    }
                    .try_build()?,
                );
                Ok(())
            }
        }
    }

    fn fetch_unmapped(&mut self) -> Result<()> {
        let source_path = self.source_path.clone();
        match &mut self.backend {
            IndexedBackend::Bam { reader, cursor, .. } => {
                // Seek to the start of the unmapped read region and stream from
                // there, filtering for unmapped records as we go.
                match reader.index().last_first_record_start_position() {
                    Some(pos) => {
                        reader
                            .get_mut()
                            .seek_to_virtual_position(pos)
                            .context("failed to seek to unmapped BAM reads")?;
                    }
                    None => {
                        // The index carries no metadata pseudo-bin, so there is
                        // no recorded unmapped offset. Reopen and stream from the
                        // first record, relying on the per-record unmapped filter.
                        let mut fresh = noodles_bam::io::indexed_reader::Builder::default()
                            .build_from_path(&source_path)
                            .with_context(|| {
                                format!("failed to reopen indexed BAM {}", source_path.display())
                            })?;
                        fresh
                            .read_header()
                            .context("failed to read indexed BAM header")?;
                        *reader = fresh;
                    }
                }
                *cursor = BamCursor::Unmapped;
                Ok(())
            }
            IndexedBackend::Cram { cursor, .. } => {
                let heads = cursor.take().expect("CRAM cursor present").into_heads();
                *cursor = Some(
                    CramIndexedCursorTryBuilder {
                        reader: heads.reader,
                        header: heads.header,
                        iter_builder: |reader, header| -> Result<CramRecordIter> {
                            let query = reader
                                .query_unmapped(header)
                                .context("failed to query unmapped CRAM reads")?;
                            Ok(Box::new(query))
                        },
                    }
                    .try_build()?,
                );
                Ok(())
            }
        }
    }
}

impl Read for IndexedReader {
    fn read(&mut self, record: &mut Record) -> Option<Result<()>> {
        match &mut self.backend {
            IndexedBackend::Bam {
                reader,
                header,
                scratch,
                cursor,
            } => loop {
                // Take ownership of the cursor each iteration so `reader` can be
                // mutated (seek/read) and the new state written back without
                // fighting the borrow checker. Early returns leave the cursor as
                // `Done`, which is the correct exhausted state.
                match std::mem::replace(cursor, BamCursor::Done) {
                    BamCursor::Done => return None,
                    BamCursor::Unmapped => match reader.read_record(scratch) {
                        Ok(0) => return None,
                        Ok(_) => {
                            *cursor = BamCursor::Unmapped;
                            if scratch.flags().is_unmapped() {
                                return Some(
                                    Record::from_bam(header.noodles_header(), scratch)
                                        .and_then(|loaded| store_record(record, loaded)),
                                );
                            }
                            // Mapped record in the unmapped tail: skip and continue.
                        }
                        Err(e) => return Some(Err(e.into())),
                    },
                    BamCursor::Region {
                        mut chunks,
                        chunk_end,
                        reference_sequence_id,
                    } => match chunk_end {
                        None => match chunks.next() {
                            Some(chunk) => {
                                if let Err(e) =
                                    reader.get_mut().seek_to_virtual_position(chunk.start())
                                {
                                    return Some(Err(anyhow::Error::new(e)
                                        .context("failed to seek to BAM index chunk")));
                                }
                                *cursor = BamCursor::Region {
                                    chunks,
                                    chunk_end: Some(chunk.end()),
                                    reference_sequence_id,
                                };
                            }
                            None => return None,
                        },
                        Some(end) => {
                            if reader.get_mut().virtual_position() >= end {
                                // Reached the end of this chunk; advance to the next.
                                *cursor = BamCursor::Region {
                                    chunks,
                                    chunk_end: None,
                                    reference_sequence_id,
                                };
                                continue;
                            }
                            match reader.read_record(scratch) {
                                Ok(0) => {
                                    *cursor = BamCursor::Region {
                                        chunks,
                                        chunk_end: None,
                                        reference_sequence_id,
                                    };
                                }
                                Ok(_) => {
                                    *cursor = BamCursor::Region {
                                        chunks,
                                        chunk_end: Some(end),
                                        reference_sequence_id,
                                    };
                                    // A chunk may begin with a few records from a
                                    // neighbouring reference; keep only ours.
                                    let on_reference = match scratch.reference_sequence_id() {
                                        Some(Ok(id)) => id == reference_sequence_id,
                                        Some(Err(e)) => return Some(Err(e.into())),
                                        None => false,
                                    };
                                    if on_reference {
                                        return Some(
                                            Record::from_bam(header.noodles_header(), scratch)
                                                .and_then(|loaded| store_record(record, loaded)),
                                        );
                                    }
                                }
                                Err(e) => return Some(Err(e.into())),
                            }
                        }
                    },
                }
            },
            IndexedBackend::Cram { cursor, .. } => {
                let cursor = cursor.as_mut()?;
                match cursor.with_iter_mut(|iter| iter.next()) {
                    Some(Ok(buf)) => Some(store_record(record, Record::from_buf(buf, None))),
                    Some(Err(e)) => Some(Err(e.into())),
                    None => None,
                }
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

#[cfg(test)]
mod tests {
    use std::num::NonZero;
    use std::sync::atomic::{AtomicU64, Ordering};

    use noodles::core::Position;
    use noodles::fasta;
    use noodles::sam::alignment::io::Write as _;
    use noodles::sam::alignment::record::cigar::{op::Kind, Op};
    use noodles::sam::alignment::record::Flags;
    use noodles::sam::alignment::RecordBuf;
    use noodles::sam::header::record::value::map::ReferenceSequence;
    use noodles::sam::header::record::value::Map;
    use noodles::{cram, sam};

    use super::{Read, Reader};
    use crate::rna::bam::Record;

    static COUNTER: AtomicU64 = AtomicU64::new(0);

    fn temp_path(suffix: &str) -> std::path::PathBuf {
        let id = COUNTER.fetch_add(1, Ordering::Relaxed);
        std::env::temp_dir().join(format!(
            "rustqc_cram_test_{}_{id}_{suffix}",
            std::process::id()
        ))
    }

    // Writes a tiny reference-backed CRAM file and streams it back through the
    // ouroboros-backed sequential `Reader`, verifying the records come through in
    // order rather than being buffered or dropped.
    #[test]
    fn streams_cram_records_sequentially() {
        let reference = b"ACGTACGTACGT".to_vec();
        let cram_path = temp_path("seq.cram");
        let fasta_path = temp_path("ref.fa");

        std::fs::write(&fasta_path, b">ref0\nACGTACGTACGT\n").unwrap();

        let repository = fasta::Repository::new(vec![fasta::Record::new(
            fasta::record::Definition::new("ref0", None),
            fasta::record::Sequence::from(reference.clone()),
        )]);

        let header = sam::Header::builder()
            .add_reference_sequence(
                "ref0",
                Map::<ReferenceSequence>::new(NonZero::new(reference.len()).unwrap()),
            )
            .build();

        let expected = [(b"r1".to_vec(), 0i64), (b"r2".to_vec(), 4i64)];
        let records: Vec<RecordBuf> = expected
            .iter()
            .map(|(name, start_zero)| {
                let start = Position::new(*start_zero as usize + 1).unwrap();
                RecordBuf::builder()
                    .set_name(&name[..])
                    .set_flags(Flags::empty())
                    .set_reference_sequence_id(0)
                    .set_alignment_start(start)
                    .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect())
                    .set_sequence(b"ACGT".to_vec().into())
                    .set_quality_scores(vec![30u8; 4].into())
                    .build()
            })
            .collect();

        {
            let mut writer = cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repository)
                .build_from_path(&cram_path)
                .unwrap();
            writer.write_header(&header).unwrap();
            for record in &records {
                writer.write_alignment_record(&header, record).unwrap();
            }
            writer.try_finish(&header).unwrap();
        }

        let mut reader = Reader::from_path(&cram_path).unwrap();
        reader.set_reference(&fasta_path).unwrap();

        let mut record = Record::new();
        let mut got = Vec::new();
        while let Some(result) = reader.read(&mut record) {
            result.unwrap();
            got.push((record.qname().to_vec(), record.pos()));
        }

        let _ = std::fs::remove_file(&cram_path);
        let _ = std::fs::remove_file(&fasta_path);

        assert_eq!(got, expected.to_vec());
    }
}
