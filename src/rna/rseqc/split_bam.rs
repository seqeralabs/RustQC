//! split_bam — classify alignments against a BED12 region list (RSeQC `split_bam.py`).
//!
//! Reimplements RSeQC 5.0.4 `split_bam.py`, typically used to quantify rRNA
//! with an rRNA BED file. Counts run inside the single pass; writing the
//! `.in.bam` / `.ex.bam` / `.junk.bam` files is optional and needs one extra pass.

use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

use anyhow::{bail, ensure, Context, Result};
use log::debug;
use rust_htslib::bam::{self, Read as BamRead};

use crate::rna::bam_flags::{BAM_FMUNMAP, BAM_FQCFAIL, BAM_FUNMAP};

// ===================================================================
// BED regions
// ===================================================================

/// Half-open `[start, end)` genomic block, 0-based.
type Block = (i64, i64);

/// Exon blocks parsed from a BED file, keyed by uppercased chromosome name.
///
/// Mirrors RSeQC's `build_bitsets(BED.ParseBED(file).getExon())`: every block
/// of every BED12 record becomes an interval, chromosome names are compared
/// case-insensitively, and an alignment is "in" when a tested position falls
/// in any block. Blocks are stored sorted and merged, which gives identical
/// answers to RSeQC's interval tree for point queries.
#[derive(Debug, Default)]
pub struct SplitBamRegions {
    /// Uppercased chromosome name → sorted, non-overlapping blocks.
    chroms: HashMap<String, Vec<Block>>,
    /// Number of BED records read.
    pub num_records: usize,
    /// Number of non-empty blocks loaded (before merging).
    pub num_blocks: usize,
    /// Number of BED records with fewer than 12 columns (treated as one block).
    pub num_short_records: usize,
}

impl SplitBamRegions {
    /// Build the per-target lookup for one alignment file.
    ///
    /// # Arguments
    /// * `target_names` - Reference names from the alignment header, indexed by tid
    pub fn lookup_for_targets<S: AsRef<str>>(&self, target_names: &[S]) -> SplitBamLookup<'_> {
        let per_tid = target_names
            .iter()
            .map(|name| {
                self.chroms
                    .get(&name.as_ref().to_uppercase())
                    .map(|blocks| blocks.as_slice())
            })
            .collect();
        SplitBamLookup { per_tid }
    }

    /// Number of distinct chromosomes with at least one non-empty block.
    pub fn num_chroms(&self) -> usize {
        self.chroms.len()
    }
}

/// Per-alignment-file view of [`SplitBamRegions`], indexed by target id.
#[derive(Debug)]
pub struct SplitBamLookup<'a> {
    /// Blocks for each target id, `None` when the target has no BED blocks.
    per_tid: Vec<Option<&'a [Block]>>,
}

impl SplitBamLookup<'_> {
    /// Number of alignment targets that have at least one BED block.
    pub fn matched_targets(&self) -> usize {
        self.per_tid.iter().filter(|b| b.is_some()).count()
    }

    /// Blocks for a target id, or `None` if the target has none (or tid < 0).
    fn blocks(&self, tid: i32) -> Option<&[Block]> {
        usize::try_from(tid)
            .ok()
            .and_then(|t| self.per_tid.get(t).copied().flatten())
    }
}

/// Return true when `pos` lies in one of the sorted, disjoint `blocks`.
///
/// Equivalent to RSeQC's `len(Intersecter.find(pos, pos + 1)) >= 1`, i.e.
/// `start <= pos < end` for some block.
fn contains(blocks: &[Block], pos: i64) -> bool {
    let i = blocks.partition_point(|&(_, end)| end <= pos);
    i < blocks.len() && blocks[i].0 <= pos
}

/// Parse one integer BED field, naming the file, line and column on failure.
fn parse_field(value: &str, column: &str, path: &str, line_no: usize) -> Result<i64> {
    value.parse::<i64>().with_context(|| {
        format!("Invalid {column} '{value}' in BED file '{path}' at line {line_no}")
    })
}

/// Parse a comma-separated BED12 list (`blockSizes` / `blockStarts`).
///
/// Matches RSeQC's `[int(i) for i in field.strip(',').split(',')]`.
fn parse_block_list(value: &str, column: &str, path: &str, line_no: usize) -> Result<Vec<i64>> {
    value
        .trim_matches(',')
        .split(',')
        .map(|v| parse_field(v, column, path, line_no))
        .collect()
}

/// Parse a BED file (plain or gzip-compressed) into split_bam regions.
///
/// BED12 records contribute one block per exon (`blockStarts`/`blockSizes`),
/// exactly like RSeQC's `getExon()`; `chromEnd` is ignored for BED12, and
/// surplus block sizes or starts are dropped as Python's `zip()` does.
///
/// Inputs RSeQC would crash on are handled as follows: blank, `#`, `track`
/// and `browser` lines are skipped, and BED3–BED11 records are treated as a
/// single `[chromStart, chromEnd)` block. Negative block sizes are rejected
/// (RSeQC's interval tree asserts on them). Zero-length blocks can never
/// match a read, so they are dropped.
///
/// # Arguments
/// * `path` - Path to the BED file
///
/// # Returns
/// The parsed regions, keyed by uppercased chromosome name.
pub fn parse_bed(path: &str) -> Result<SplitBamRegions> {
    use std::io::BufRead;

    let reader = crate::io::open_reader(path)
        .with_context(|| format!("Failed to open split_bam BED file: {path}"))?;
    let mut regions = SplitBamRegions::default();

    for (idx, line) in reader.lines().enumerate() {
        let line_no = idx + 1;
        let line = line.with_context(|| format!("Failed to read BED file '{path}'"))?;
        let trimmed = line.trim();
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed.starts_with("track")
            || trimmed.starts_with("browser")
        {
            continue;
        }

        // RSeQC splits on any whitespace (`str.split()`), not only tabs.
        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        ensure!(
            fields.len() >= 3,
            "BED file '{path}' line {line_no}: expected at least 3 columns, found {}",
            fields.len()
        );
        let chrom = fields[0].to_uppercase();
        let chrom_start = parse_field(fields[1], "chromStart", path, line_no)?;

        let blocks: Vec<Block> = if fields.len() >= 12 {
            // Columns RSeQC parses as integers but does not otherwise use.
            for (col, name) in [(6, "thickStart"), (7, "thickEnd"), (9, "blockCount")] {
                parse_field(fields[col], name, path, line_no)?;
            }
            let sizes = parse_block_list(fields[10], "blockSizes", path, line_no)?;
            let starts = parse_block_list(fields[11], "blockStarts", path, line_no)?;
            starts
                .iter()
                .zip(sizes.iter())
                .map(|(&offset, &size)| (chrom_start + offset, chrom_start + offset + size))
                .collect()
        } else {
            regions.num_short_records += 1;
            let chrom_end = parse_field(fields[2], "chromEnd", path, line_no)?;
            vec![(chrom_start, chrom_end)]
        };

        for (start, end) in blocks {
            if end < start {
                bail!(
                    "BED file '{path}' line {line_no}: block end ({end}) is before its start ({start})"
                );
            }
            if end > start {
                regions
                    .chroms
                    .entry(chrom.clone())
                    .or_default()
                    .push((start, end));
                regions.num_blocks += 1;
            }
        }
        regions.num_records += 1;
    }

    ensure!(
        regions.num_records > 0,
        "BED file '{path}' contains no records"
    );

    // Sort and merge so point queries can use a binary search.
    for blocks in regions.chroms.values_mut() {
        blocks.sort_unstable();
        let mut merged: Vec<Block> = Vec::with_capacity(blocks.len());
        for &(start, end) in blocks.iter() {
            match merged.last_mut() {
                Some(last) if start <= last.1 => last.1 = last.1.max(end),
                _ => merged.push((start, end)),
            }
        }
        *blocks = merged;
    }

    debug!(
        "split_bam: {} BED records, {} blocks on {} chromosomes from {}",
        regions.num_records,
        regions.num_blocks,
        regions.chroms.len(),
        path
    );
    Ok(regions)
}

// ===================================================================
// Classification
// ===================================================================

/// Output category of a single alignment record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SplitBamClass {
    /// Read (or mate) start falls in a BED block → `.in.bam`.
    In,
    /// Mapped, QC-passing, not in any block → `.ex.bam`.
    Ex,
    /// QC-failed or unmapped → `.junk.bam`.
    Junk,
}

/// Classify one alignment record exactly as RSeQC 5.0.4 `split_bam.py` does.
///
/// - QC-failed or unmapped records are junk. No other flag is filtered, so
///   secondary, supplementary and duplicate alignments are all classified.
/// - Otherwise the record is "in" when its own 0-based start lies in a block
///   on its chromosome; or, when the mate-unmapped flag (0x8) is *not* set,
///   when the mate's start does. RSeQC tests the mate position against the
///   *read's* chromosome even if the mate maps elsewhere — reproduced here.
pub fn classify(record: &bam::Record, lookup: &SplitBamLookup) -> SplitBamClass {
    let flags = record.flags();
    if flags & (BAM_FQCFAIL | BAM_FUNMAP) != 0 {
        return SplitBamClass::Junk;
    }
    let Some(blocks) = lookup.blocks(record.tid()) else {
        return SplitBamClass::Ex;
    };
    let hit = contains(blocks, record.pos())
        || (flags & BAM_FMUNMAP == 0 && contains(blocks, record.mpos()));
    if hit {
        SplitBamClass::In
    } else {
        SplitBamClass::Ex
    }
}

// ===================================================================
// Counts / accumulator
// ===================================================================

/// split_bam record counts (the numbers RSeQC prints to stdout).
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct SplitBamCounts {
    /// All records seen (`Total records:`).
    pub total: u64,
    /// Records written to `.in.bam`.
    pub in_count: u64,
    /// Records written to `.ex.bam`.
    pub ex_count: u64,
    /// Records written to `.junk.bam` (QC-failed or unmapped).
    pub junk_count: u64,
}

impl SplitBamCounts {
    /// Count one classified record.
    fn add(&mut self, class: SplitBamClass) {
        self.total += 1;
        match class {
            SplitBamClass::In => self.in_count += 1,
            SplitBamClass::Ex => self.ex_count += 1,
            SplitBamClass::Junk => self.junk_count += 1,
        }
    }
}

/// Single-pass split_bam accumulator; one per worker, merged afterwards.
#[derive(Debug, Default)]
pub struct SplitBamAccum {
    counts: SplitBamCounts,
}

impl SplitBamAccum {
    /// Classify and count one record. Must see every record in the file.
    pub fn process_read(&mut self, record: &bam::Record, lookup: &SplitBamLookup) {
        self.counts.add(classify(record, lookup));
    }

    /// Merge counts from another worker.
    pub fn merge(&mut self, other: SplitBamAccum) {
        self.counts.total += other.counts.total;
        self.counts.in_count += other.counts.in_count;
        self.counts.ex_count += other.counts.ex_count;
        self.counts.junk_count += other.counts.junk_count;
    }

    /// Final counts.
    pub fn into_result(self) -> SplitBamCounts {
        self.counts
    }
}

// ===================================================================
// Output
// ===================================================================

/// Write the split_bam summary, byte-identical to RSeQC's stdout.
///
/// RSeQC prints four `"%-55s%d"` lines whose labels embed the `-o` prefix;
/// `prefix` plays that role here (RustQC uses the sample name).
///
/// # Arguments
/// * `counts` - Final record counts
/// * `prefix` - Output prefix shown in the `.in.bam` / `.ex.bam` / `.junk.bam` labels
/// * `path` - File to write
pub fn write_split_bam_stats(counts: &SplitBamCounts, prefix: &str, path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .with_context(|| format!("Failed to create split_bam output: {}", path.display()))?;
    let lines = [
        ("Total records:".to_string(), counts.total),
        (
            format!("{prefix}.in.bam (Alignments consumed by input gene list):"),
            counts.in_count,
        ),
        (
            format!("{prefix}.ex.bam (Alignments not consumed by input gene list):"),
            counts.ex_count,
        ),
        (
            format!("{prefix}.junk.bam (qcfailed, unmapped reads):"),
            counts.junk_count,
        ),
    ];
    for (label, n) in lines {
        writeln!(out, "{label:<55}{n}")?;
    }
    Ok(())
}

/// Write `<prefix>.in.bam`, `<prefix>.ex.bam` and `<prefix>.junk.bam`.
///
/// Runs a separate sequential pass so records keep the input order and header,
/// as `split_bam.py` produces (it copies the input header unchanged and adds no
/// `@PG` line). Returns the counts from this pass so callers can cross-check
/// them against the single-pass accumulator.
///
/// # Arguments
/// * `bam_path` - Input alignment file
/// * `reference` - Reference FASTA (CRAM input only)
/// * `regions` - Parsed BED regions
/// * `prefix` - Output path prefix (directory + sample name)
/// * `threads` - Threads for BGZF decompression/compression
pub fn write_split_bams(
    bam_path: &str,
    reference: Option<&str>,
    regions: &SplitBamRegions,
    prefix: &Path,
    threads: usize,
) -> Result<SplitBamCounts> {
    let mut reader = bam::Reader::from_path(bam_path)
        .with_context(|| format!("Failed to open alignment file: {bam_path}"))?;
    if let Some(ref_path) = reference {
        reader
            .set_reference(ref_path)
            .with_context(|| format!("Failed to set reference FASTA: {ref_path}"))?;
    }
    let pool = if threads > 1 {
        let pool = rust_htslib::tpool::ThreadPool::new(threads as u32)
            .context("Failed to create htslib thread pool")?;
        reader
            .set_thread_pool(&pool)
            .context("Failed to attach thread pool to reader")?;
        Some(pool)
    } else {
        None
    };

    let header_view = reader.header().clone();
    let target_names: Vec<String> = (0..header_view.target_count())
        .map(|tid| String::from_utf8_lossy(header_view.tid2name(tid)).into_owned())
        .collect();
    let lookup = regions.lookup_for_targets(&target_names);
    let header = bam::Header::from_template(&header_view);

    let prefix = prefix.to_string_lossy();
    let open = |suffix: &str| -> Result<bam::Writer> {
        let path = format!("{prefix}.{suffix}.bam");
        let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam)
            .with_context(|| format!("Failed to create {path}"))?;
        if let Some(ref pool) = pool {
            writer
                .set_thread_pool(pool)
                .context("Failed to attach thread pool to writer")?;
        }
        Ok(writer)
    };
    let mut in_writer = open("in")?;
    let mut ex_writer = open("ex")?;
    let mut junk_writer = open("junk")?;

    let mut counts = SplitBamCounts::default();
    let mut record = bam::Record::new();
    while let Some(result) = reader.read(&mut record) {
        result.with_context(|| format!("Error reading alignment record from {bam_path}"))?;
        let class = classify(&record, &lookup);
        let writer = match class {
            SplitBamClass::In => &mut in_writer,
            SplitBamClass::Ex => &mut ex_writer,
            SplitBamClass::Junk => &mut junk_writer,
        };
        writer
            .write(&record)
            .context("Failed to write split_bam record")?;
        counts.add(class);
    }
    Ok(counts)
}

// ===================================================================
// Tests
// ===================================================================

#[cfg(test)]
mod tests {
    use super::*;

    /// Write `content` to a unique temporary file and return its path.
    fn temp_file(name: &str, content: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("rustqc_split_bam_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(name);
        std::fs::write(&path, content).unwrap();
        path
    }

    fn regions_from(name: &str, bed: &str) -> SplitBamRegions {
        let path = temp_file(name, bed);
        parse_bed(path.to_str().unwrap()).unwrap()
    }

    #[test]
    fn test_contains_half_open() {
        let blocks = [(100, 200), (300, 400)];
        assert!(contains(&blocks, 100));
        assert!(contains(&blocks, 199));
        assert!(!contains(&blocks, 200));
        assert!(!contains(&blocks, 99));
        assert!(!contains(&blocks, 250));
        assert!(contains(&blocks, 300));
        assert!(!contains(&blocks, 400));
        assert!(!contains(&blocks, -1));
        assert!(!contains(&[], 5));
    }

    #[test]
    fn test_parse_bed12_uses_exon_blocks_only() {
        let r = regions_from(
            "blocks.bed",
            "chr1\t1000\t3000\tg\t0\t+\t1000\t3000\t0\t2\t500,500,\t0,1500,\n",
        );
        let lookup = r.lookup_for_targets(&["chr1"]);
        let b = lookup.blocks(0).unwrap();
        assert_eq!(b, &[(1000, 1500), (2500, 3000)]);
        assert!(!contains(b, 2000), "intron must not match");
    }

    #[test]
    fn test_parse_bed_case_insensitive_chrom() {
        let r = regions_from("case.bed", "Chr2\t5\t10\n");
        let lookup = r.lookup_for_targets(&["chr1", "CHR2"]);
        assert!(lookup.blocks(0).is_none());
        assert_eq!(lookup.blocks(1).unwrap(), &[(5, 10)]);
        assert_eq!(lookup.matched_targets(), 1);
        assert_eq!(r.num_short_records, 1);
    }

    #[test]
    fn test_parse_bed_zip_truncates_and_drops_zero_length() {
        // Two sizes but one start → one block (Python zip); zero-length dropped.
        let r = regions_from(
            "zip.bed",
            "chr1\t9000\t9500\ta\t0\t+\t9000\t9500\t0\t2\t100,100,\t0,\n\
             chr1\t8000\t8000\tb\t0\t+\t8000\t8000\t0\t1\t0,\t0,\n",
        );
        let lookup = r.lookup_for_targets(&["chr1"]);
        assert_eq!(lookup.blocks(0).unwrap(), &[(9000, 9100)]);
        assert_eq!(r.num_records, 2);
        assert_eq!(r.num_blocks, 1);
    }

    #[test]
    fn test_parse_bed_merges_overlaps_and_skips_headers() {
        let r = regions_from(
            "merge.bed",
            "track name=x\n# comment\nbrowser position chr1\n\n\
             chr1 10 20\nchr1 15 30\nchr1 40 50\nchr1 30 35\n",
        );
        let lookup = r.lookup_for_targets(&["chr1"]);
        assert_eq!(lookup.blocks(0).unwrap(), &[(10, 35), (40, 50)]);
    }

    #[test]
    fn test_parse_bed_rejects_bad_input() {
        let neg = temp_file("neg.bed", "chr1\t10\t20\ta\t0\t+\t10\t20\t0\t1\t-5,\t0,\n");
        assert!(parse_bed(neg.to_str().unwrap()).is_err());
        let short = temp_file("short.bed", "chr1\t10\n");
        assert!(parse_bed(short.to_str().unwrap()).is_err());
        let empty = temp_file("empty.bed", "# nothing\n");
        assert!(parse_bed(empty.to_str().unwrap()).is_err());
        let bad_int = temp_file(
            "badint.bed",
            "chr1\t10\t20\ta\t0\t+\tx\t20\t0\t1\t10,\t0,\n",
        );
        assert!(parse_bed(bad_int.to_str().unwrap()).is_err());
    }

    /// Build a record from minimal SAM fields via a one-record SAM file.
    fn records_from_sam(name: &str, body: &str) -> Vec<bam::Record> {
        let sam = format!("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:20000\n@SQ\tSN:chr2\tLN:20000\n{body}");
        let path = temp_file(name, &sam);
        let mut reader = bam::Reader::from_path(&path).unwrap();
        reader.records().map(|r| r.unwrap()).collect()
    }

    #[test]
    fn test_classify_matches_rseqc_branches() {
        let r = regions_from("cls.bed", "chr1\t1000\t1500\n");
        let lookup = r.lookup_for_targets(&["chr1", "chr2"]);
        let q = "A".repeat(10);
        let i = "I".repeat(10);
        let line = |name: &str, flag: u16, chrom: &str, pos: u32, rnext: &str, pnext: u32| {
            format!("{name}\t{flag}\t{chrom}\t{pos}\t60\t10M\t{rnext}\t{pnext}\t0\t{q}\t{i}\n")
        };
        let body = [
            line("in_self", 0, "chr1", 1101, "*", 0),
            line("ex_self", 0, "chr1", 2001, "*", 0),
            line("in_mate", 99, "chr1", 4001, "=", 1101),
            // Mate on chr2, but its position is tested against chr1 (RSeQC quirk).
            line("in_disc", 65, "chr1", 6001, "chr2", 1101),
            // Mate-unmapped flag: mate position ignored.
            line("ex_munmap", 73, "chr1", 7001, "=", 1101),
            line("ex_nochrom", 0, "chr2", 1101, "*", 0),
            line("junk_qc", 512, "chr1", 1101, "*", 0),
            line("in_secondary", 256, "chr1", 1101, "*", 0),
            "junk_unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n".to_string(),
        ]
        .concat();
        let records = records_from_sam("cls.sam", &body);
        let expected = [
            SplitBamClass::In,
            SplitBamClass::Ex,
            SplitBamClass::In,
            SplitBamClass::In,
            SplitBamClass::Ex,
            SplitBamClass::Ex,
            SplitBamClass::Junk,
            SplitBamClass::In,
            SplitBamClass::Junk,
        ];
        assert_eq!(records.len(), expected.len());
        let mut accum = SplitBamAccum::default();
        for (rec, exp) in records.iter().zip(expected) {
            let name = String::from_utf8_lossy(rec.qname()).to_string();
            assert_eq!(classify(rec, &lookup), exp, "record {name}");
            accum.process_read(rec, &lookup);
        }
        let counts = accum.into_result();
        assert_eq!(
            counts,
            SplitBamCounts {
                total: 9,
                in_count: 4,
                ex_count: 3,
                junk_count: 2,
            }
        );
    }

    #[test]
    fn test_merge_is_additive() {
        let mut a = SplitBamAccum {
            counts: SplitBamCounts {
                total: 3,
                in_count: 1,
                ex_count: 1,
                junk_count: 1,
            },
        };
        let b = SplitBamAccum {
            counts: SplitBamCounts {
                total: 2,
                in_count: 2,
                ex_count: 0,
                junk_count: 0,
            },
        };
        a.merge(b);
        assert_eq!(
            a.into_result(),
            SplitBamCounts {
                total: 5,
                in_count: 3,
                ex_count: 1,
                junk_count: 1,
            }
        );
    }

    #[test]
    fn test_stats_format_matches_rseqc() {
        let counts = SplitBamCounts {
            total: 33,
            in_count: 15,
            ex_count: 11,
            junk_count: 7,
        };
        let path = temp_file("split_bam_summary.txt", "");
        write_split_bam_stats(&counts, "split_bam_test", &path).unwrap();
        let got = std::fs::read_to_string(&path).unwrap();
        // Byte-for-byte output of RSeQC 5.0.4 split_bam.py -o split_bam_test
        let expected = "Total records:                                         33\n\
split_bam_test.in.bam (Alignments consumed by input gene list):15\n\
split_bam_test.ex.bam (Alignments not consumed by input gene list):11\n\
split_bam_test.junk.bam (qcfailed, unmapped reads):    7\n";
        assert_eq!(got, expected);
    }
}
