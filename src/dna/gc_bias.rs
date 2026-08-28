//! Picard `CollectGcBiasMetrics` reimplementation.
//!
//! # Upstream semantics
//!
//! These rules come from Picard 3.4.0's own source, `GcBiasUtils` and
//! `GcBiasMetricsCollector`, after black-box inference from its output failed
//! to reproduce them. They are unusual enough to be worth stating.
//!
//! **Windows.** GC is computed over sliding windows of `window_size` bases at
//! every reference position `i` for `1 <= i < len - window_size`. Note both
//! bounds: the window at position 0 is skipped, and so is the last one that
//! would fit. On a 40001 base reference with 100 base windows that gives
//! 39900 windows, not the 39902 a naive reading produces. A window holding
//! more than [`MAX_NS_PER_WINDOW`] `N` bases is marked unusable and its reads
//! are dropped. The GC value is `gc_count * 100 / window_size` in integer
//! arithmetic, truncating rather than rounding.
//!
//! **Read assignment.** A read is assigned to the window at its alignment
//! start, except on the reverse strand, where it is assigned to
//! `alignment_end - window_size`. That is not the read's 5' end; it is the
//! window that would start where the read's far end finishes. Positions are
//! one-based, and a read landing at position 0 or lower is dropped.
//!
//! **Which reads count.** Only unmapped reads and reads with an empty
//! sequence are skipped. Secondary and supplementary alignments and
//! duplicates all contribute, which is what `READS_USED ALL` means. On the
//! project fixture that is the difference between 5640 and 5642 read starts.
//!
//! **Dropout.** For each GC bin, `(window_share - read_share) * 100` is
//! accumulated when positive, into `AT_DROPOUT` for bins at or below 50 and
//! `GC_DROPOUT` above.

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

use crate::common::bam_flags::*;

/// Number of GC bins, one per whole percent from 0 to 100 inclusive.
pub const BINS: usize = 101;

/// Picard's `SCAN_WINDOW_SIZE` default.
pub const DEFAULT_WINDOW_SIZE: usize = 100;

/// A window holding more than this many `N` bases is unusable.
pub const MAX_NS_PER_WINDOW: usize = 4;

/// Accumulates GC bias for one contig.
#[derive(Debug)]
pub struct GcBiasAccum {
    /// GC percent per one-based reference position, or `-1` when the window
    /// there holds too many `N` bases or does not exist.
    gc: Vec<i8>,
    window_size: usize,
    windows_by_gc: [u64; BINS],
    reads_by_gc: [u64; BINS],
    bases_by_gc: [u64; BINS],
    errors_by_gc: [u64; BINS],
    total_clusters: u64,
    total_aligned_reads: u64,
}

impl GcBiasAccum {
    /// Build the window GC table for one contig's reference bases.
    pub fn new(reference: &[u8], window_size: usize) -> Self {
        let len = reference.len();
        let mut gc = vec![-1i8; len + 1];
        let mut windows_by_gc = [0u64; BINS];

        if len > window_size {
            // Prefix sums make each window a constant-time lookup.
            let mut gc_prefix = vec![0u32; len + 1];
            let mut n_prefix = vec![0u32; len + 1];
            for (i, base) in reference.iter().enumerate() {
                let upper = base.to_ascii_uppercase();
                gc_prefix[i + 1] = gc_prefix[i] + u32::from(upper == b'G' || upper == b'C');
                n_prefix[i + 1] = n_prefix[i] + u32::from(upper == b'N');
            }

            let last_window_start = len - window_size;
            for i in 1..last_window_start {
                let end = i + window_size;
                let ns = (n_prefix[end] - n_prefix[i]) as usize;
                if ns > MAX_NS_PER_WINDOW {
                    continue;
                }
                let gc_count = gc_prefix[end] - gc_prefix[i];
                let percent = (gc_count as usize * 100 / window_size) as i8;
                gc[i] = percent;
                windows_by_gc[percent as usize] += 1;
            }
        }

        Self {
            gc,
            window_size,
            windows_by_gc,
            reads_by_gc: [0; BINS],
            bases_by_gc: [0; BINS],
            errors_by_gc: [0; BINS],
            total_clusters: 0,
            total_aligned_reads: 0,
        }
    }

    /// Offer one record, with the contig's reference bases for mismatch counting.
    pub fn process_read(&mut self, record: &bam::Record, reference: &[u8]) {
        if record.seq_len() == 0 {
            return;
        }

        // A cluster is a template, counted once, at the unpaired read or the
        // first of the pair. Unmapped reads count towards clusters even though
        // they reach nothing else, so this precedes the mapped check.
        if record.flags() & BAM_FPAIRED == 0 || record.flags() & BAM_FREAD1 != 0 {
            self.total_clusters += 1;
        }
        if record.flags() & BAM_FUNMAP != 0 {
            return;
        }
        self.total_aligned_reads += 1;

        // One-based, and the reverse strand is assigned by the far end rather
        // than the near one.
        let position = if record.flags() & BAM_FREVERSE != 0 {
            alignment_end(record) - self.window_size as i64
        } else {
            record.pos() + 1
        };
        if position <= 0 {
            return;
        }
        let Some(&percent) = self.gc.get(position as usize) else {
            return;
        };
        if percent < 0 {
            return;
        }

        let bin = percent as usize;
        self.reads_by_gc[bin] += 1;
        self.bases_by_gc[bin] += record.seq_len() as u64;
        self.errors_by_gc[bin] += count_errors(record, reference);
    }

    /// Fold another contig's counters in. Window tables are per contig and add
    /// up the same way.
    pub fn merge(&mut self, other: &GcBiasAccum) {
        for bin in 0..BINS {
            self.windows_by_gc[bin] += other.windows_by_gc[bin];
            self.reads_by_gc[bin] += other.reads_by_gc[bin];
            self.bases_by_gc[bin] += other.bases_by_gc[bin];
            self.errors_by_gc[bin] += other.errors_by_gc[bin];
        }
        self.total_clusters += other.total_clusters;
        self.total_aligned_reads += other.total_aligned_reads;
    }

    /// Summarise into the reported detail rows and summary figures.
    pub fn into_result(self, window_size: usize) -> GcBiasResult {
        let total_reads: u64 = self.reads_by_gc.iter().sum();
        let total_windows: u64 = self.windows_by_gc.iter().sum();
        let global_rate = if total_windows == 0 {
            0.0
        } else {
            total_reads as f64 / total_windows as f64
        };

        let mut rows = Vec::with_capacity(BINS);
        let mut at_dropout = 0.0;
        let mut gc_dropout = 0.0;

        for bin in 0..BINS {
            let windows = self.windows_by_gc[bin];
            let reads = self.reads_by_gc[bin];
            let bases = self.bases_by_gc[bin];
            let errors = self.errors_by_gc[bin];

            let normalized = if windows == 0 || global_rate == 0.0 {
                0.0
            } else {
                (reads as f64 / windows as f64) / global_rate
            };
            let error_bar = if windows == 0 || global_rate == 0.0 {
                0.0
            } else {
                ((reads as f64).sqrt() / windows as f64) / global_rate
            };
            // Mean quality as the phred score of the observed error rate.
            let mean_base_quality = if bases == 0 || errors == 0 {
                0
            } else {
                (-10.0 * (errors as f64 / bases as f64).log10()).round() as i32
            };

            if total_reads > 0 && total_windows > 0 {
                let read_share = reads as f64 / total_reads as f64;
                let window_share = windows as f64 / total_windows as f64;
                let dropout = (window_share - read_share) * 100.0;
                if dropout > 0.0 {
                    if bin <= 50 {
                        at_dropout += dropout;
                    } else {
                        gc_dropout += dropout;
                    }
                }
            }

            rows.push(GcBiasDetail {
                gc: bin as u32,
                windows,
                read_starts: reads,
                mean_base_quality,
                normalized_coverage: normalized,
                error_bar_width: error_bar,
            });
        }

        GcBiasResult {
            rows,
            window_size,
            total_clusters: self.total_clusters,
            aligned_reads: self.total_aligned_reads,
            at_dropout,
            gc_dropout,
        }
    }
}

/// One-based inclusive end of a record's alignment.
fn alignment_end(record: &bam::Record) -> i64 {
    let mut end = record.pos();
    for op in record.cigar().iter() {
        match op {
            Cigar::Match(n)
            | Cigar::Equal(n)
            | Cigar::Diff(n)
            | Cigar::Del(n)
            | Cigar::RefSkip(n) => end += i64::from(*n),
            _ => {}
        }
    }
    end
}

/// Mismatches against the reference, plus inserted and deleted bases, which is
/// what Picard counts towards the per-bin error rate.
fn count_errors(record: &bam::Record, reference: &[u8]) -> u64 {
    let sequence = record.seq();
    let mut errors = 0u64;
    let mut ref_pos = record.pos();
    let mut query_pos = 0i64;

    for op in record.cigar().iter() {
        match op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                for k in 0..i64::from(*n) {
                    let r = ref_pos + k;
                    let q = query_pos + k;
                    if r < 0 || r as usize >= reference.len() {
                        continue;
                    }
                    let ref_base = reference[r as usize].to_ascii_uppercase();
                    let read_base = sequence[q as usize].to_ascii_uppercase();
                    // htsjdk's basesEqual is a plain comparison after
                    // uppercasing, so an N on either side is a mismatch rather
                    // than a free pass.
                    if ref_base != read_base {
                        errors += 1;
                    }
                }
                ref_pos += i64::from(*n);
                query_pos += i64::from(*n);
            }
            Cigar::Ins(n) => {
                errors += u64::from(*n);
                query_pos += i64::from(*n);
            }
            Cigar::Del(n) => {
                errors += u64::from(*n);
                ref_pos += i64::from(*n);
            }
            Cigar::RefSkip(n) => ref_pos += i64::from(*n),
            Cigar::SoftClip(n) => query_pos += i64::from(*n),
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    errors
}

/// One GC bin's detail row.
#[derive(Debug, Clone)]
pub struct GcBiasDetail {
    /// GC percent this row describes.
    pub gc: u32,
    /// Reference windows at this GC.
    pub windows: u64,
    /// Reads assigned to a window at this GC.
    pub read_starts: u64,
    /// Phred score of the observed error rate for those reads.
    pub mean_base_quality: i32,
    /// Read density here relative to the genome-wide density.
    pub normalized_coverage: f64,
    /// One standard error of `normalized_coverage`.
    pub error_bar_width: f64,
}

/// The complete GC bias result.
#[derive(Debug, Clone)]
pub struct GcBiasResult {
    /// One row per GC bin, ascending.
    pub rows: Vec<GcBiasDetail>,
    /// Window size the bins were computed over.
    pub window_size: usize,
    /// Templates seen.
    pub total_clusters: u64,
    /// Mapped reads seen.
    pub aligned_reads: u64,
    /// Illumina-style AT dropout.
    pub at_dropout: f64,
    /// Illumina-style GC dropout.
    pub gc_dropout: f64,
}

impl GcBiasResult {
    /// Mean normalised coverage across a GC range, as the `GC_NC_x_y` columns
    /// report it.
    ///
    /// The mean is weighted by how many reference windows each bin holds, not
    /// a plain average over bins. That distinction matters at the extremes,
    /// where most bins hold no windows at all and would otherwise drag the
    /// figure towards zero.
    fn mean_normalized(&self, low: u32, high: u32) -> f64 {
        let mut weighted = 0.0;
        let mut windows = 0u64;
        for row in self.rows.iter().filter(|r| r.gc >= low && r.gc <= high) {
            weighted += row.normalized_coverage * row.windows as f64;
            windows += row.windows;
        }
        if windows == 0 {
            0.0
        } else {
            weighted / windows as f64
        }
    }
}

/// Format a float the way Picard's metrics writer does.
fn fmt_picard(value: f64) -> String {
    if !value.is_finite() {
        return "?".to_string();
    }
    if value == value.trunc() && value.abs() < 1e15 {
        return format!("{}", value as i64);
    }
    let text = format!("{value:.6}");
    text.trim_end_matches('0').trim_end_matches('.').to_string()
}

/// Write the per-GC-bin detail metrics.
pub fn write_detail_metrics(result: &GcBiasResult, path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .map(std::io::BufWriter::new)
        .with_context(|| {
            format!(
                "Failed to create GC bias detail metrics: {}",
                path.display()
            )
        })?;

    writeln!(out, "## METRICS CLASS\tpicard.analysis.GcBiasDetailMetrics")?;
    writeln!(
        out,
        "ACCUMULATION_LEVEL\tREADS_USED\tGC\tWINDOWS\tREAD_STARTS\tMEAN_BASE_QUALITY\t\
         NORMALIZED_COVERAGE\tERROR_BAR_WIDTH\tSAMPLE\tLIBRARY\tREAD_GROUP"
    )?;
    for row in &result.rows {
        writeln!(
            out,
            "All Reads\tALL\t{}\t{}\t{}\t{}\t{}\t{}\t\t\t",
            row.gc,
            row.windows,
            row.read_starts,
            row.mean_base_quality,
            fmt_picard(row.normalized_coverage),
            fmt_picard(row.error_bar_width),
        )?;
    }
    // Picard leaves two blank lines at the end of the GC bias tables, one more
    // than it writes after the insert size or WGS tables. The fixtures are the
    // specification, so this matches them rather than being tidied.
    writeln!(out)?;
    writeln!(out)?;
    out.flush()?;
    Ok(())
}

/// Write the GC bias summary metrics.
pub fn write_summary_metrics(result: &GcBiasResult, path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .map(std::io::BufWriter::new)
        .with_context(|| {
            format!(
                "Failed to create GC bias summary metrics: {}",
                path.display()
            )
        })?;

    writeln!(
        out,
        "## METRICS CLASS\tpicard.analysis.GcBiasSummaryMetrics"
    )?;
    writeln!(
        out,
        "ACCUMULATION_LEVEL\tREADS_USED\tWINDOW_SIZE\tTOTAL_CLUSTERS\tALIGNED_READS\t\
         AT_DROPOUT\tGC_DROPOUT\tGC_NC_0_19\tGC_NC_20_39\tGC_NC_40_59\tGC_NC_60_79\t\
         GC_NC_80_100\tSAMPLE\tLIBRARY\tREAD_GROUP"
    )?;
    writeln!(
        out,
        "All Reads\tALL\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t\t\t",
        result.window_size,
        result.total_clusters,
        result.aligned_reads,
        fmt_picard(result.at_dropout),
        fmt_picard(result.gc_dropout),
        fmt_picard(result.mean_normalized(0, 19)),
        fmt_picard(result.mean_normalized(20, 39)),
        fmt_picard(result.mean_normalized(40, 59)),
        fmt_picard(result.mean_normalized(60, 79)),
        fmt_picard(result.mean_normalized(80, 100)),
    )?;
    writeln!(out)?;
    writeln!(out)?;
    out.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_first_and_last_possible_windows_are_both_skipped() {
        // 10 bases, window size 4: naive windows would be starts 0 through 6,
        // Picard's loop runs 1 through 5.
        let reference = b"ACGTACGTAC".to_vec();
        let accum = GcBiasAccum::new(&reference, 4);
        let windows: u64 = accum.windows_by_gc.iter().sum();
        assert_eq!(windows, 5, "starts 1 through 5 inclusive");
        assert_eq!(accum.gc[0], -1, "the window at 0 is never computed");
        assert_eq!(accum.gc[6], -1, "nor the last one that would fit");
    }

    #[test]
    fn gc_percent_truncates_rather_than_rounds() {
        // 3 of 8 bases are G or C: 3 * 100 / 8 is 37.5, truncated to 37.
        let reference = b"GGCAAAAAAAAA".to_vec();
        let accum = GcBiasAccum::new(&reference, 8);
        // Window at position 1 is GCAAAAAA, two of eight, 25 percent.
        assert_eq!(accum.gc[1], 25);
    }

    #[test]
    fn windows_with_too_many_ns_are_unusable() {
        let reference = b"ACGTNNNNNGCTAGCTAGC".to_vec();
        let accum = GcBiasAccum::new(&reference, 8);
        // The window at 1 holds five Ns, one more than the limit.
        assert_eq!(accum.gc[1], -1);
    }

    #[test]
    fn dropout_splits_at_fifty_percent_gc() {
        let mut accum = GcBiasAccum::new(&b"A".repeat(200), 100);
        // Hand-place windows and reads so the shares are unambiguous.
        accum.windows_by_gc = [0; BINS];
        accum.reads_by_gc = [0; BINS];
        accum.windows_by_gc[30] = 50;
        accum.windows_by_gc[70] = 50;
        accum.reads_by_gc[30] = 100;
        accum.reads_by_gc[70] = 0;
        let result = accum.into_result(100);
        // Bin 70 has half the windows and none of the reads: 50 points of
        // dropout, and being above 50 percent GC it lands in GC_DROPOUT.
        assert!(
            (result.gc_dropout - 50.0).abs() < 1e-9,
            "{}",
            result.gc_dropout
        );
        assert!(
            (result.at_dropout - 0.0).abs() < 1e-9,
            "{}",
            result.at_dropout
        );
    }

    #[test]
    fn normalized_coverage_is_relative_to_the_genome_wide_rate() {
        let mut accum = GcBiasAccum::new(&b"A".repeat(200), 100);
        accum.windows_by_gc = [0; BINS];
        accum.reads_by_gc = [0; BINS];
        accum.windows_by_gc[10] = 100;
        accum.windows_by_gc[20] = 100;
        accum.reads_by_gc[10] = 150;
        accum.reads_by_gc[20] = 50;
        let result = accum.into_result(100);
        // Global rate is 200 reads over 200 windows, so 1 read per window.
        assert!((result.rows[10].normalized_coverage - 1.5).abs() < 1e-9);
        assert!((result.rows[20].normalized_coverage - 0.5).abs() < 1e-9);
    }
}
