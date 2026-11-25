//! mosdepth-equivalent per-base and per-window depth.
//!
//! Depth is computed exactly from CIGAR-aware start/end delta events. Because
//! the input is coordinate-sorted, every base before the current read's start
//! is final and is folded into the per-contig summary and the current window
//! as the file streams past, so memory stays proportional to pile-up depth
//! rather than to contig length.

use crate::rna::bam_flags::{BAM_FDUP, BAM_FQCFAIL, BAM_FSECONDARY, BAM_FUNMAP};
use rust_htslib::bam::{self, record::Cigar};
use std::collections::BTreeMap;

/// Highest depth tracked in the cumulative distribution.
///
/// Deeper positions are folded into the top bin, matching mosdepth's bounded
/// coverage array.
pub const MAX_DIST_DEPTH: u32 = 1000;

/// Records excluded from depth, matching mosdepth's default `--flag 1796`
/// (UNMAP, SECONDARY, QCFAIL, DUP).
const EXCLUDED_FLAGS: u16 = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;

/// Per-contig depth results.
#[derive(Debug, Clone)]
pub struct ContigDepth {
    /// Contig name.
    pub name: String,
    /// Contig length from the alignment header.
    pub length: u64,
    /// Σ depth over the contig (aligned bases).
    pub total_bases: u64,
    /// Minimum per-base depth (0 unless the contig is fully covered).
    pub min_depth: u32,
    /// Maximum per-base depth.
    pub max_depth: u32,
    /// depth (capped at [`MAX_DIST_DEPTH`]) → number of bases at that depth.
    pub hist: BTreeMap<u32, u64>,
}

impl ContigDepth {
    /// Mean depth across the whole contig, uncovered bases counted as zero.
    pub fn mean(&self) -> f64 {
        if self.length == 0 {
            0.0
        } else {
            self.total_bases as f64 / self.length as f64
        }
    }
}

/// One fixed-size window of the reference.
#[derive(Debug, Clone)]
pub struct Window {
    /// Contig name.
    pub chrom: String,
    /// Window start (0-based).
    pub start: u64,
    /// Window end (exclusive).
    pub end: u64,
    /// Mean depth across the window.
    pub mean: f64,
}

/// Streaming depth accumulator over a coordinate-sorted alignment file.
#[derive(Debug)]
pub struct DepthAccum {
    /// Contig names and lengths, in header order.
    contigs: Vec<(String, u64)>,
    /// Window size in bases.
    window_size: u64,
    /// Finished per-contig results.
    results: Vec<ContigDepth>,
    /// Emitted windows, in reference order.
    windows: Vec<Window>,

    // --- state for the contig being processed ---
    current_tid: i32,
    events: BTreeMap<u64, i64>,
    last_pos: u64,
    depth: i64,
    hist: BTreeMap<u32, u64>,
    total_bases: u64,
    max_depth: u32,
    /// Σ depth per window of the contig being processed.
    window_sums: Vec<u64>,
}

impl DepthAccum {
    /// Create an accumulator for a file with the given contigs.
    ///
    /// # Arguments
    /// * `contigs` - Reference names and lengths from the alignment header
    /// * `window_size` - Window size in bases (mosdepth's `--by`)
    pub fn new(contigs: Vec<(String, u64)>, window_size: u64) -> Self {
        Self {
            contigs,
            window_size: window_size.max(1),
            results: Vec::new(),
            windows: Vec::new(),
            current_tid: -1,
            events: BTreeMap::new(),
            last_pos: 0,
            depth: 0,
            hist: BTreeMap::new(),
            total_bases: 0,
            max_depth: 0,
            window_sums: Vec::new(),
        }
    }

    /// Add one alignment record.
    pub fn process_read(&mut self, record: &bam::Record) {
        if record.flags() & EXCLUDED_FLAGS != 0 {
            return;
        }
        let tid = record.tid();
        if tid < 0 {
            return;
        }
        if tid != self.current_tid {
            self.finish_contig();
            self.current_tid = tid;
        }

        let start = record.pos() as u64;
        self.advance_to(start);

        let mut ref_pos = start;
        for op in record.cigar().iter() {
            match *op {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) | Cigar::Del(len) => {
                    let len = len as u64;
                    *self.events.entry(ref_pos).or_insert(0) += 1;
                    *self.events.entry(ref_pos + len).or_insert(0) -= 1;
                    ref_pos += len;
                }
                Cigar::RefSkip(len) => ref_pos += len as u64,
                _ => {}
            }
        }
    }

    /// Fold all depth up to (but excluding) `pos` into the current contig.
    fn advance_to(&mut self, pos: u64) {
        while let Some((&event_pos, &delta)) = self.events.iter().next() {
            if event_pos >= pos {
                break;
            }
            if event_pos > self.last_pos {
                let span = event_pos - self.last_pos;
                let depth = self.depth.max(0) as u32;
                self.record(depth, span);
            }
            self.depth += delta;
            self.last_pos = event_pos;
            self.events.remove(&event_pos);
        }
        if pos > self.last_pos {
            let span = pos - self.last_pos;
            let depth = self.depth.max(0) as u32;
            self.record(depth, span);
            self.last_pos = pos;
        }
    }

    /// Record `span` consecutive bases at `depth`, splitting across windows.
    fn record(&mut self, depth: u32, span: u64) {
        if span == 0 {
            return;
        }
        if depth > 0 {
            *self.hist.entry(depth.min(MAX_DIST_DEPTH)).or_insert(0) += span;
            self.total_bases += depth as u64 * span;
            self.max_depth = self.max_depth.max(depth);
        }

        if depth == 0 {
            return;
        }

        // Split the span across window boundaries
        let mut pos = self.last_pos;
        let end = self.last_pos + span;
        while pos < end {
            let idx = (pos / self.window_size) as usize;
            let window_end = (idx as u64 + 1) * self.window_size;
            let chunk_end = window_end.min(end);
            if idx >= self.window_sums.len() {
                self.window_sums.resize(idx + 1, 0);
            }
            self.window_sums[idx] += depth as u64 * (chunk_end - pos);
            pos = chunk_end;
        }
    }

    /// Emit every window of the contig being processed, including empty ones.
    ///
    /// mosdepth tiles the whole reference when `--by <N>` is a fixed window
    /// size, so windows with no coverage are written with a mean of 0.00.
    fn flush_windows(&mut self) {
        if self.current_tid < 0 || (self.current_tid as usize) >= self.contigs.len() {
            self.window_sums.clear();
            return;
        }
        let (name, contig_len) = self.contigs[self.current_tid as usize].clone();
        let num_windows = contig_len.div_ceil(self.window_size);
        for idx in 0..num_windows {
            let start = idx * self.window_size;
            let end = ((idx + 1) * self.window_size).min(contig_len);
            let span = end - start;
            if span == 0 {
                continue;
            }
            let sum = self.window_sums.get(idx as usize).copied().unwrap_or(0);
            self.windows.push(Window {
                chrom: name.clone(),
                start,
                end,
                mean: sum as f64 / span as f64,
            });
        }
        self.window_sums.clear();
    }

    /// Finalise the contig currently being processed.
    fn finish_contig(&mut self) {
        if self.current_tid >= 0 {
            let last_event = self.events.keys().next_back().copied().unwrap_or(0);
            self.advance_to(last_event + 1);
            self.flush_windows();

            if (self.current_tid as usize) < self.contigs.len() {
                let (name, length) = self.contigs[self.current_tid as usize].clone();
                let covered: u64 = self.hist.values().sum();
                let min_depth = if covered < length {
                    0
                } else {
                    self.hist.keys().next().copied().unwrap_or(0)
                };
                self.results.push(ContigDepth {
                    name,
                    length,
                    total_bases: self.total_bases,
                    min_depth,
                    max_depth: self.max_depth,
                    hist: std::mem::take(&mut self.hist),
                });
            }
        }

        self.events.clear();
        self.hist.clear();
        self.depth = 0;
        self.last_pos = 0;
        self.total_bases = 0;
        self.max_depth = 0;
        self.window_sums.clear();
    }

    /// Finalise everything and return the per-contig results and windows.
    ///
    /// Contigs with no alignments are reported with zero coverage so the
    /// summary lists every reference in the header, as mosdepth does.
    pub fn finish(mut self) -> (Vec<ContigDepth>, Vec<Window>) {
        self.finish_contig();

        let seen: std::collections::HashSet<&str> =
            self.results.iter().map(|r| r.name.as_str()).collect();
        let mut results = Vec::with_capacity(self.contigs.len());
        for (name, length) in &self.contigs {
            if seen.contains(name.as_str()) {
                continue;
            }
            results.push(ContigDepth {
                name: name.clone(),
                length: *length,
                total_bases: 0,
                min_depth: 0,
                max_depth: 0,
                hist: BTreeMap::new(),
            });
            let num_windows = length.div_ceil(self.window_size);
            for idx in 0..num_windows {
                let start = idx * self.window_size;
                let end = ((idx + 1) * self.window_size).min(*length);
                if end > start {
                    self.windows.push(Window {
                        chrom: name.clone(),
                        start,
                        end,
                        mean: 0.0,
                    });
                }
            }
        }
        // Keep header order
        let mut all = self.results;
        all.extend(results);
        let order: std::collections::HashMap<&str, usize> = self
            .contigs
            .iter()
            .enumerate()
            .map(|(i, (name, _))| (name.as_str(), i))
            .collect();
        all.sort_by_key(|r| order.get(r.name.as_str()).copied().unwrap_or(usize::MAX));

        (all, self.windows)
    }
}

// ===================================================================
// Tests
// ===================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::Read as BamRead;

    fn accumulate(
        sam: &str,
        contigs: Vec<(String, u64)>,
        by: u64,
    ) -> (Vec<ContigDepth>, Vec<Window>) {
        use std::sync::atomic::{AtomicU64, Ordering};
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        let id = COUNTER.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "rustqc_depth_test_{:?}_{}.sam",
            std::thread::current().id(),
            id
        ));
        std::fs::write(&path, sam).unwrap();

        let mut reader = bam::Reader::from_path(&path).unwrap();
        let mut accum = DepthAccum::new(contigs, by);
        let mut record = bam::Record::new();
        while let Some(res) = reader.read(&mut record) {
            res.unwrap();
            accum.process_read(&record);
        }
        let _ = std::fs::remove_file(&path);
        accum.finish()
    }

    #[test]
    fn test_windows_and_contig_summary() {
        // Contig of 100 bases, window size 50.
        // r1 covers 1..=10, r2 covers 6..=15 (1-based) -> window 0 gets
        // 20 aligned bases, window 1 gets none.
        let sam = "\
@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:100\n\
r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n\
r2\t0\tchr1\t6\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n";
        let (contigs, windows) = accumulate(sam, vec![("chr1".to_string(), 100)], 50);

        assert_eq!(contigs.len(), 1);
        assert_eq!(contigs[0].total_bases, 20);
        assert!((contigs[0].mean() - 0.2).abs() < 1e-9);
        assert_eq!(contigs[0].max_depth, 2);
        assert_eq!(contigs[0].min_depth, 0, "contig is not fully covered");
        assert_eq!(contigs[0].hist.get(&1), Some(&10));
        assert_eq!(contigs[0].hist.get(&2), Some(&5));

        // mosdepth tiles the whole contig, so both windows are emitted
        assert_eq!(windows.len(), 2);
        assert_eq!(windows[0].start, 0);
        assert_eq!(windows[0].end, 50);
        assert!((windows[0].mean - 20.0 / 50.0).abs() < 1e-9);
        assert_eq!(windows[1].start, 50);
        assert_eq!(windows[1].mean, 0.0, "empty windows are emitted as 0");
    }

    #[test]
    fn test_window_split_across_boundary() {
        // 20-base read starting at 1-based 41 spans windows [0,50) and [50,100)
        let sam = "\
@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:100\n\
r1\t0\tchr1\t41\t60\t20M\t*\t0\t0\tACGTACGTACACGTACGTAC\tIIIIIIIIIIIIIIIIIIII\n";
        let (_, windows) = accumulate(sam, vec![("chr1".to_string(), 100)], 50);

        assert_eq!(windows.len(), 2);
        // 10 bases in each window
        assert!((windows[0].mean - 10.0 / 50.0).abs() < 1e-9);
        assert!((windows[1].mean - 10.0 / 50.0).abs() < 1e-9);
    }

    #[test]
    fn test_excluded_flags_and_empty_contigs() {
        // r1 is a duplicate, r2 secondary, r3 unmapped: none contribute.
        // chr2 has no reads at all and must still be reported.
        let sam = "\
@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:100\n\
@SQ\tSN:chr2\tLN:200\n\
r1\t1024\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n\
r2\t256\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n\
r3\t4\tchr1\t1\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n";
        let (contigs, windows) = accumulate(
            sam,
            vec![("chr1".to_string(), 100), ("chr2".to_string(), 200)],
            50,
        );

        assert_eq!(contigs.len(), 2, "every header contig is reported");
        assert_eq!(contigs[0].name, "chr1");
        assert_eq!(contigs[0].total_bases, 0);
        assert_eq!(contigs[1].name, "chr2");
        assert_eq!(contigs[1].total_bases, 0);
        // 2 windows over chr1 (100 bp) + 4 over chr2 (200 bp), all empty
        assert_eq!(windows.len(), 6);
        assert!(windows.iter().all(|w| w.mean == 0.0));
    }
}
