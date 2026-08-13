//! BED-interval read classification, equivalent to RSeQC's `split_bam.py`.
//!
//! Classifies every alignment record against a BED file of genomic intervals
//! (typically rRNA regions) into three categories — `in` (overlapping an
//! interval), `ex` (not overlapping) and `junk` (unmapped or QC-failed) — and
//! reports the counts and percentages.
//!
//! This complements the GTF/biotype route to rRNA quantification, which
//! structurally under-reports rRNA on stock GRCh38/GENCODE builds: the 45S
//! rDNA repeat is not annotated there, and rDNA multi-mappers are dropped by
//! featureCounts' default single-hit rule. Interval overlap catches both.

use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use anyhow::{Context, Result};
use coitrees::{BasicCOITree, Interval, IntervalTree};
use log::debug;

use crate::rna::bam_flags::{BAM_FQCFAIL, BAM_FUNMAP};
use rust_htslib::bam;

/// Per-chromosome interval tree over the BED regions.
type ChromIntervals = BasicCOITree<(), u32>;

/// Genomic intervals loaded from a BED file, indexed per chromosome.
///
/// `Debug` is implemented by hand because the coitrees tree type does not
/// implement it.
pub struct BedIntervals {
    /// Chromosome name → interval tree of that chromosome's BED records.
    trees: HashMap<String, ChromIntervals>,
    /// Total number of intervals loaded (for logging).
    pub num_intervals: usize,
}

impl std::fmt::Debug for BedIntervals {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BedIntervals")
            .field("chromosomes", &self.trees.len())
            .field("num_intervals", &self.num_intervals)
            .finish()
    }
}

impl BedIntervals {
    /// Return `true` when `pos` (0-based) falls inside any interval on `chrom`.
    fn contains(&self, chrom: &str, pos: i64) -> bool {
        match self.trees.get(chrom) {
            Some(tree) => {
                let p = pos as i32;
                tree.query_count(p, p) > 0
            }
            None => false,
        }
    }

    /// Return `true` when the chromosome is present in the BED file at all.
    fn has_chrom(&self, chrom: &str) -> bool {
        self.trees.contains_key(chrom)
    }
}

/// Parse a BED file (plain or gzip-compressed) into per-chromosome interval trees.
///
/// Only the first three columns (chrom, start, end) are used. BED coordinates
/// are 0-based half-open, and are stored as 0-based inclusive `[start, end-1]`
/// so that point queries against BAM positions line up.
///
/// # Arguments
/// * `path` - Path to the BED file
///
/// # Returns
/// Indexed intervals ready for overlap queries
pub fn parse_bed(path: &str) -> Result<BedIntervals> {
    let reader =
        crate::io::open_reader(path).with_context(|| format!("Failed to open BED file: {path}"))?;

    let mut per_chrom: HashMap<String, Vec<Interval<()>>> = HashMap::new();
    let mut num_intervals = 0usize;

    for (lineno, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("Failed to read line from BED file: {path}"))?;
        let trimmed = line.trim_end();
        // Skip comments, empty lines and browser/track headers
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed.starts_with("track")
            || trimmed.starts_with("browser")
        {
            continue;
        }

        let mut fields = trimmed.split('\t');
        let (Some(chrom), Some(start_s), Some(end_s)) =
            (fields.next(), fields.next(), fields.next())
        else {
            anyhow::bail!(
                "Malformed BED file '{}' at line {}: expected at least 3 tab-separated columns",
                path,
                lineno + 1
            );
        };

        let start: i64 = start_s.parse().with_context(|| {
            format!(
                "Malformed BED file '{}' at line {}: invalid start position '{}'",
                path,
                lineno + 1,
                start_s
            )
        })?;
        let end: i64 = end_s.parse().with_context(|| {
            format!(
                "Malformed BED file '{}' at line {}: invalid end position '{}'",
                path,
                lineno + 1,
                end_s
            )
        })?;
        anyhow::ensure!(
            end > start,
            "Malformed BED file '{}' at line {}: end ({}) must be greater than start ({})",
            path,
            lineno + 1,
            end,
            start
        );

        per_chrom
            .entry(chrom.to_string())
            .or_default()
            // BED is half-open; coitrees intervals are inclusive on both ends
            .push(Interval::new(start as i32, (end - 1) as i32, ()));
        num_intervals += 1;
    }

    anyhow::ensure!(
        num_intervals > 0,
        "No usable intervals found in BED file '{}'",
        path
    );

    let trees = per_chrom
        .into_iter()
        .map(|(chrom, ivs)| (chrom, ChromIntervals::new(&ivs)))
        .collect();

    debug!("Loaded {num_intervals} intervals from BED file {path}");
    Ok(BedIntervals {
        trees,
        num_intervals,
    })
}

// ===================================================================
// Accumulator
// ===================================================================

/// Per-worker accumulator classifying reads against the BED intervals.
#[derive(Debug, Default, Clone)]
pub struct SplitBamAccum {
    /// Records overlapping at least one BED interval.
    pub in_count: u64,
    /// Mapped, QC-passing records not overlapping any BED interval.
    pub ex_count: u64,
    /// Unmapped or QC-failed records.
    pub junk_count: u64,
}

impl SplitBamAccum {
    /// Classify a single alignment record.
    ///
    /// Follows RSeQC `split_bam.py`: unmapped and QC-failed records are "junk";
    /// otherwise the read's own start position, and (when the mate is mapped)
    /// the mate's start position, are tested against the intervals. A hit on
    /// either end puts the record in the "in" category.
    ///
    /// Secondary and supplementary alignments are *not* filtered out — this is
    /// deliberate, and is what lets the interval route see the rDNA
    /// multi-mappers that featureCounts' default single-hit rule discards.
    ///
    /// # Arguments
    /// * `record` - The alignment record
    /// * `chrom` - Chromosome name for this record, already mapped to BED naming
    /// * `intervals` - The BED intervals to test against
    pub fn process_read(&mut self, record: &bam::Record, chrom: &str, intervals: &BedIntervals) {
        let flags = record.flags();
        if flags & (BAM_FUNMAP | BAM_FQCFAIL) != 0 {
            self.junk_count += 1;
            return;
        }

        if !intervals.has_chrom(chrom) {
            self.ex_count += 1;
            return;
        }

        let hit = intervals.contains(chrom, record.pos())
            || (record.mtid() == record.tid()
                && record.mpos() >= 0
                && intervals.contains(chrom, record.mpos()));

        if hit {
            self.in_count += 1;
        } else {
            self.ex_count += 1;
        }
    }

    /// Merge another worker's counts into this accumulator.
    pub fn merge(&mut self, other: SplitBamAccum) {
        self.in_count += other.in_count;
        self.ex_count += other.ex_count;
        self.junk_count += other.junk_count;
    }

    /// Finalise into a result.
    pub fn into_result(self) -> SplitBamResult {
        SplitBamResult {
            in_count: self.in_count,
            ex_count: self.ex_count,
            junk_count: self.junk_count,
        }
    }
}

/// Final classification counts.
#[derive(Debug, Clone)]
pub struct SplitBamResult {
    /// Records overlapping at least one BED interval.
    pub in_count: u64,
    /// Mapped, QC-passing records not overlapping any BED interval.
    pub ex_count: u64,
    /// Unmapped or QC-failed records.
    pub junk_count: u64,
}

impl SplitBamResult {
    /// Total records classified.
    pub fn total(&self) -> u64 {
        self.in_count + self.ex_count + self.junk_count
    }

    /// Percentage of mapped, QC-passing records that fall inside the intervals.
    ///
    /// Junk records are excluded from the denominator, so the value answers
    /// "what fraction of usable alignments are rRNA?".
    pub fn percent_in(&self) -> f64 {
        let denom = self.in_count + self.ex_count;
        if denom == 0 {
            0.0
        } else {
            100.0 * self.in_count as f64 / denom as f64
        }
    }
}

// ===================================================================
// Output
// ===================================================================

/// Write the classification summary as a TSV file.
///
/// # Arguments
/// * `result` - The computed classification counts
/// * `bed_path` - BED file the intervals came from (recorded in the header)
/// * `output_path` - Path to write the summary to
pub fn write_split_bam_summary(
    result: &SplitBamResult,
    bed_path: &str,
    output_path: &Path,
) -> Result<()> {
    use std::io::Write;

    let mut out = std::fs::File::create(output_path).with_context(|| {
        format!(
            "Failed to create split_bam summary file: {}",
            output_path.display()
        )
    })?;

    let total = result.total();
    let usable = result.in_count + result.ex_count;
    let pct = |n: u64, d: u64| {
        if d == 0 {
            0.0
        } else {
            100.0 * n as f64 / d as f64
        }
    };

    writeln!(out, "# BED intervals: {bed_path}")?;
    writeln!(
        out,
        "# 'in' = alignment overlaps an interval, 'ex' = does not, \
         'junk' = unmapped or QC-failed"
    )?;
    writeln!(
        out,
        "# percent_of_usable excludes junk records from the denominator"
    )?;
    writeln!(out, "category\tcount\tpercent_of_total\tpercent_of_usable")?;
    writeln!(
        out,
        "in\t{}\t{:.4}\t{:.4}",
        result.in_count,
        pct(result.in_count, total),
        pct(result.in_count, usable)
    )?;
    writeln!(
        out,
        "ex\t{}\t{:.4}\t{:.4}",
        result.ex_count,
        pct(result.ex_count, total),
        pct(result.ex_count, usable)
    )?;
    writeln!(
        out,
        "junk\t{}\t{:.4}\tNA",
        result.junk_count,
        pct(result.junk_count, total)
    )?;
    writeln!(out, "total\t{total}\t100.0000\tNA")?;

    debug!(
        "Wrote split_bam summary to {} ({} in, {} ex, {} junk)",
        output_path.display(),
        result.in_count,
        result.ex_count,
        result.junk_count
    );
    Ok(())
}

// ===================================================================
// Tests
// ===================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn write_temp(content: &str, ext: &str) -> std::path::PathBuf {
        use std::io::Write;
        use std::sync::atomic::{AtomicU64, Ordering};
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        let id = COUNTER.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "rustqc_split_bam_test_{:?}_{}.{}",
            std::thread::current().id(),
            id,
            ext
        ));
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f.flush().unwrap();
        path
    }

    #[test]
    fn test_parse_bed_basic() {
        let path = write_temp(
            "track name=rRNA\n\
             # a comment\n\
             chr1\t100\t200\trRNA1\t0\t+\n\
             chr1\t300\t400\n\
             chr2\t50\t60\n",
            "bed",
        );
        let intervals = parse_bed(path.to_str().unwrap()).unwrap();
        let _ = std::fs::remove_file(&path);

        assert_eq!(intervals.num_intervals, 3);
        // BED is half-open: [100, 200) covers 100..=199
        assert!(intervals.contains("chr1", 100));
        assert!(intervals.contains("chr1", 199));
        assert!(!intervals.contains("chr1", 200));
        assert!(!intervals.contains("chr1", 99));
        assert!(intervals.contains("chr2", 55));
        assert!(!intervals.has_chrom("chr3"));
    }

    #[test]
    fn test_parse_bed_rejects_malformed_line() {
        let path = write_temp("chr1\t100\n", "bed");
        let err = parse_bed(path.to_str().unwrap()).unwrap_err();
        let _ = std::fs::remove_file(&path);
        assert!(
            err.to_string().contains("at least 3 tab-separated columns"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn test_parse_bed_rejects_empty_file() {
        let path = write_temp("# nothing here\n", "bed");
        let err = parse_bed(path.to_str().unwrap()).unwrap_err();
        let _ = std::fs::remove_file(&path);
        assert!(
            err.to_string().contains("No usable intervals"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn test_classification_counts_multimappers_and_mates() {
        use rust_htslib::bam::Read as BamRead;

        let bed = write_temp("chr1\t500\t600\n", "bed");
        let intervals = parse_bed(bed.to_str().unwrap()).unwrap();
        let _ = std::fs::remove_file(&bed);

        // r1: primary, starts inside the interval          -> in
        // r2: secondary (NH>1), starts inside the interval -> in (not filtered)
        // r3: primary, outside                             -> ex
        // r4: primary, outside but mate starts inside      -> in
        // r5: unmapped                                     -> junk
        let sam = "\
@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:20000\n\
r1\t99\tchr1\t501\t30\t10M\t=\t900\t0\tACGTACGTAC\tIIIIIIIIII\n\
r2\t355\tchr1\t520\t30\t10M\t=\t900\t0\tACGTACGTAC\tIIIIIIIIII\n\
r3\t99\tchr1\t900\t30\t10M\t=\t950\t0\tACGTACGTAC\tIIIIIIIIII\n\
r4\t147\tchr1\t950\t30\t10M\t=\t520\t0\tACGTACGTAC\tIIIIIIIIII\n\
r5\t77\tchr1\t960\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n";
        let sam_path = write_temp(sam, "sam");

        let mut reader = bam::Reader::from_path(&sam_path).unwrap();
        let mut accum = SplitBamAccum::default();
        let mut record = bam::Record::new();
        while let Some(res) = reader.read(&mut record) {
            res.unwrap();
            accum.process_read(&record, "chr1", &intervals);
        }
        let _ = std::fs::remove_file(&sam_path);

        let result = accum.into_result();
        assert_eq!(result.in_count, 3, "in");
        assert_eq!(result.ex_count, 1, "ex");
        assert_eq!(result.junk_count, 1, "junk");
        assert_eq!(result.total(), 5);
        assert!((result.percent_in() - 75.0).abs() < 1e-9, "percent_in");
    }
}
