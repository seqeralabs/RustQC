//! Single-pass accumulation of generic BAM QC metrics.
//!
//! Collects the statistics Qualimap's `bamqc` mode reports: global read
//! counters, per-base coverage (genome-wide and per contig), insert size,
//! mapping quality, read length and base composition. No annotation is
//! required — this is genomic QC, not RNA-seq QC.

use std::collections::{BTreeMap, HashMap};

use rust_htslib::bam::{self, record::Cigar};

use crate::rna::bam_flags::*;

/// Coverage thresholds reported as "fraction of reference covered ≥ N×".
pub const COVERAGE_THRESHOLDS: [u32; 7] = [1, 5, 10, 20, 30, 40, 50];

// ===================================================================
// Depth tracking
// ===================================================================

/// Running summary of per-base depth over one contig.
#[derive(Debug, Default, Clone)]
struct DepthSummary {
    /// depth → number of reference bases at that depth (depth ≥ 1 only).
    hist: HashMap<u32, u64>,
    /// Σ depth over covered bases.
    sum_depth: u128,
    /// Σ depth² over covered bases.
    sum_sq: u128,
    /// Number of reference bases with depth ≥ 1.
    covered: u64,
}

impl DepthSummary {
    fn record(&mut self, depth: u32, span: u64) {
        if depth == 0 || span == 0 {
            return;
        }
        *self.hist.entry(depth).or_insert(0) += span;
        self.sum_depth += depth as u128 * span as u128;
        self.sum_sq += (depth as u128).pow(2) * span as u128;
        self.covered += span;
    }

    fn merge(&mut self, other: &DepthSummary) {
        for (&depth, &span) in &other.hist {
            *self.hist.entry(depth).or_insert(0) += span;
        }
        self.sum_depth += other.sum_depth;
        self.sum_sq += other.sum_sq;
        self.covered += other.covered;
    }
}

/// Exact per-base depth accumulator for a coordinate-sorted alignment file.
///
/// Reads contribute `+1` at the start and `-1` at the end of every aligned
/// block (`M`/`=`/`X`/`D`; `N` skips leave a gap). Because input is sorted by
/// start position, every base before the current read's start is final and can
/// be folded into the summary immediately, so memory stays proportional to the
/// number of overlapping reads rather than to the contig length.
#[derive(Debug, Default)]
struct DepthTracker {
    /// Pending depth deltas keyed by reference position.
    events: BTreeMap<u64, i64>,
    /// Position up to which depth has been folded into the summary.
    last_pos: u64,
    /// Current depth at `last_pos`.
    depth: i64,
    /// Summary for the contig being processed.
    current: DepthSummary,
}

impl DepthTracker {
    /// Add an aligned block covering `[start, end)` (0-based, half-open).
    fn add_block(&mut self, start: u64, end: u64) {
        if end <= start {
            return;
        }
        *self.events.entry(start).or_insert(0) += 1;
        *self.events.entry(end).or_insert(0) -= 1;
    }

    /// Fold all depth up to (but excluding) `pos` into the summary.
    fn advance_to(&mut self, pos: u64) {
        while let Some((&event_pos, &delta)) = self.events.iter().next() {
            if event_pos >= pos {
                break;
            }
            if event_pos > self.last_pos {
                let span = event_pos - self.last_pos;
                self.current.record(self.depth.max(0) as u32, span);
            }
            self.depth += delta;
            self.last_pos = event_pos;
            self.events.remove(&event_pos);
        }
        if pos > self.last_pos {
            self.current
                .record(self.depth.max(0) as u32, pos - self.last_pos);
            self.last_pos = pos;
        }
    }

    /// Fold every remaining event and return the finished contig summary.
    fn finish_contig(&mut self) -> DepthSummary {
        let last_event = self.events.keys().next_back().copied().unwrap_or(0);
        self.advance_to(last_event + 1);
        self.events.clear();
        self.depth = 0;
        self.last_pos = 0;
        std::mem::take(&mut self.current)
    }
}

// ===================================================================
// Accumulator
// ===================================================================

/// Collects every metric reported by [`BamqcResult`] in a single pass.
#[derive(Debug)]
pub struct BamqcAccum {
    // --- reference ---
    /// Contig names and lengths, in header order.
    contigs: Vec<(String, u64)>,
    /// Per-contig depth summaries, indexed by tid.
    contig_depth: Vec<DepthSummary>,
    /// Per-contig mapped read counts, indexed by tid.
    contig_reads: Vec<u64>,

    // --- global counters ---
    total_reads: u64,
    mapped_reads: u64,
    unmapped_reads: u64,
    secondary: u64,
    supplementary: u64,
    duplicates: u64,
    qc_failed: u64,
    paired_reads: u64,
    first_in_pair: u64,
    second_in_pair: u64,
    both_mates_mapped: u64,
    singletons: u64,

    // --- lengths and bases ---
    read_length_min: u64,
    read_length_max: u64,
    read_length_sum: u64,
    read_length_n: u64,
    /// Bases in the query sequences of mapped reads.
    sequenced_bases: u64,
    /// Reference bases covered by aligned blocks (M/=/X/D).
    aligned_bases: u64,
    /// Base counts: A, C, G, T, N.
    base_counts: [u64; 5],

    // --- distributions ---
    /// GC percentage (rounded to whole percent) → read count.
    gc_hist: [u64; 101],
    /// MAPQ → read count (mapped, primary).
    mapq_hist: [u64; 256],
    /// abs(TLEN) → pair count, one entry per fragment.
    insert_hist: BTreeMap<u64, u64>,

    // --- depth state ---
    depth: DepthTracker,
    /// tid of the contig currently being tracked.
    current_tid: i32,
    /// When true, duplicate-flagged records contribute to nothing but the
    /// duplicate counter (Qualimap's `--skip-duplicated`).
    skip_duplicated: bool,
}

impl BamqcAccum {
    /// Create an accumulator for a file with the given contigs.
    ///
    /// # Arguments
    /// * `contigs` - Reference names and lengths from the alignment header
    /// * `skip_duplicated` - Exclude duplicate-flagged records from all metrics
    pub fn new(contigs: Vec<(String, u64)>, skip_duplicated: bool) -> Self {
        let n = contigs.len();
        Self {
            contigs,
            contig_depth: vec![DepthSummary::default(); n],
            contig_reads: vec![0; n],
            total_reads: 0,
            mapped_reads: 0,
            unmapped_reads: 0,
            secondary: 0,
            supplementary: 0,
            duplicates: 0,
            qc_failed: 0,
            paired_reads: 0,
            first_in_pair: 0,
            second_in_pair: 0,
            both_mates_mapped: 0,
            singletons: 0,
            read_length_min: u64::MAX,
            read_length_max: 0,
            read_length_sum: 0,
            read_length_n: 0,
            sequenced_bases: 0,
            aligned_bases: 0,
            base_counts: [0; 5],
            gc_hist: [0; 101],
            mapq_hist: [0; 256],
            insert_hist: BTreeMap::new(),
            depth: DepthTracker::default(),
            current_tid: -1,
            skip_duplicated,
        }
    }

    /// Process one alignment record.
    ///
    /// Coverage counts primary mapped alignments (secondary alignments are
    /// excluded; duplicates are included, matching Qualimap's default of not
    /// skipping flagged duplicates).
    pub fn process_read(&mut self, record: &bam::Record) {
        let flags = record.flags();
        self.total_reads += 1;

        let is_unmapped = flags & BAM_FUNMAP != 0;
        let is_secondary = flags & BAM_FSECONDARY != 0;
        let is_supplementary = flags & BAM_FSUPPLEMENTARY != 0;
        let is_primary = !is_secondary && !is_supplementary;

        if is_secondary {
            self.secondary += 1;
        } else if is_supplementary {
            self.supplementary += 1;
        }
        let is_duplicate = flags & BAM_FDUP != 0;
        if is_duplicate {
            self.duplicates += 1;
            if self.skip_duplicated {
                return;
            }
        }
        if flags & BAM_FQCFAIL != 0 {
            self.qc_failed += 1;
        }

        if is_unmapped {
            self.unmapped_reads += 1;
            return;
        }
        self.mapped_reads += 1;

        if is_primary {
            if flags & BAM_FPAIRED != 0 {
                self.paired_reads += 1;
                if flags & BAM_FREAD1 != 0 {
                    self.first_in_pair += 1;
                }
                if flags & BAM_FREAD2 != 0 {
                    self.second_in_pair += 1;
                }
                if flags & BAM_FMUNMAP == 0 {
                    self.both_mates_mapped += 1;
                    // Count each fragment once, from the leftmost mate
                    let isize = record.insert_size();
                    if isize > 0 {
                        *self.insert_hist.entry(isize as u64).or_insert(0) += 1;
                    }
                } else {
                    self.singletons += 1;
                }
            }
            self.mapq_hist[record.mapq() as usize] += 1;

            let seq = record.seq();
            let len = seq.len() as u64;
            if len > 0 {
                self.read_length_min = self.read_length_min.min(len);
                self.read_length_max = self.read_length_max.max(len);
                self.read_length_sum += len;
                self.read_length_n += 1;
                self.sequenced_bases += len;

                let mut gc = 0u64;
                for i in 0..seq.len() {
                    match seq[i] {
                        b'A' | b'a' => self.base_counts[0] += 1,
                        b'C' | b'c' => {
                            self.base_counts[1] += 1;
                            gc += 1;
                        }
                        b'G' | b'g' => {
                            self.base_counts[2] += 1;
                            gc += 1;
                        }
                        b'T' | b't' => self.base_counts[3] += 1,
                        _ => self.base_counts[4] += 1,
                    }
                }
                let pct = (100.0 * gc as f64 / len as f64).round() as usize;
                self.gc_hist[pct.min(100)] += 1;
            }
        }

        if is_secondary {
            return;
        }

        // --- coverage ---
        let tid = record.tid();
        if tid < 0 {
            return;
        }
        if tid != self.current_tid {
            self.flush_contig();
            self.current_tid = tid;
        }
        if (tid as usize) < self.contig_reads.len() {
            self.contig_reads[tid as usize] += 1;
        }

        let start = record.pos() as u64;
        self.depth.advance_to(start);

        let mut ref_pos = start;
        for op in record.cigar().iter() {
            match *op {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) | Cigar::Del(len) => {
                    let len = len as u64;
                    self.depth.add_block(ref_pos, ref_pos + len);
                    self.aligned_bases += len;
                    ref_pos += len;
                }
                Cigar::RefSkip(len) => ref_pos += len as u64,
                _ => {}
            }
        }
    }

    /// Fold the current contig's pending depth into its summary.
    fn flush_contig(&mut self) {
        let summary = self.depth.finish_contig();
        if self.current_tid >= 0 && (self.current_tid as usize) < self.contig_depth.len() {
            let idx = self.current_tid as usize;
            let taken = summary;
            self.contig_depth[idx].merge(&taken);
        }
    }

    /// Finalise all pending state and produce the result.
    pub fn into_result(mut self) -> BamqcResult {
        self.flush_contig();

        let mut global = DepthSummary::default();
        for summary in &self.contig_depth {
            global.merge(summary);
        }

        let reference_bases: u64 = self.contigs.iter().map(|(_, len)| len).sum();

        let contigs = self
            .contigs
            .iter()
            .enumerate()
            .map(|(i, (name, len))| ContigCoverage {
                name: name.clone(),
                length: *len,
                mapped_reads: self.contig_reads[i],
                mapped_bases: self.contig_depth[i].sum_depth as u64,
                mean_coverage: mean_coverage(&self.contig_depth[i], *len),
                std_coverage: std_coverage(&self.contig_depth[i], *len),
            })
            .collect();

        let mapq_sum: u128 = self
            .mapq_hist
            .iter()
            .enumerate()
            .map(|(q, &n)| q as u128 * n as u128)
            .sum();
        let mapq_n: u64 = self.mapq_hist.iter().sum();

        let (insert_mean, insert_std, insert_median) = insert_size_stats(&self.insert_hist);

        BamqcResult {
            reference_bases,
            num_contigs: self.contigs.len(),
            total_reads: self.total_reads,
            mapped_reads: self.mapped_reads,
            unmapped_reads: self.unmapped_reads,
            secondary: self.secondary,
            supplementary: self.supplementary,
            duplicates: self.duplicates,
            qc_failed: self.qc_failed,
            paired_reads: self.paired_reads,
            first_in_pair: self.first_in_pair,
            second_in_pair: self.second_in_pair,
            both_mates_mapped: self.both_mates_mapped,
            singletons: self.singletons,
            read_length_min: if self.read_length_n == 0 {
                0
            } else {
                self.read_length_min
            },
            read_length_max: self.read_length_max,
            read_length_mean: if self.read_length_n == 0 {
                0.0
            } else {
                self.read_length_sum as f64 / self.read_length_n as f64
            },
            sequenced_bases: self.sequenced_bases,
            aligned_bases: self.aligned_bases,
            base_counts: self.base_counts,
            gc_hist: self.gc_hist,
            mapq_mean: if mapq_n == 0 {
                0.0
            } else {
                mapq_sum as f64 / mapq_n as f64
            },
            insert_mean,
            insert_std,
            insert_median,
            insert_hist: self.insert_hist,
            mean_coverage: mean_coverage(&global, reference_bases),
            std_coverage: std_coverage(&global, reference_bases),
            coverage_hist: global.hist.clone(),
            covered_bases: global.covered,
            contigs,
        }
    }
}

/// Mean depth over the whole span, counting uncovered bases as zero.
fn mean_coverage(summary: &DepthSummary, span: u64) -> f64 {
    if span == 0 {
        return 0.0;
    }
    summary.sum_depth as f64 / span as f64
}

/// Population standard deviation of depth, counting uncovered bases as zero.
fn std_coverage(summary: &DepthSummary, span: u64) -> f64 {
    if span == 0 {
        return 0.0;
    }
    let mean = summary.sum_depth as f64 / span as f64;
    let mean_sq = summary.sum_sq as f64 / span as f64;
    (mean_sq - mean * mean).max(0.0).sqrt()
}

/// Mean, population standard deviation and median of the insert size histogram.
fn insert_size_stats(hist: &BTreeMap<u64, u64>) -> (f64, f64, u64) {
    let total: u64 = hist.values().sum();
    if total == 0 {
        return (0.0, 0.0, 0);
    }
    let sum: u128 = hist.iter().map(|(&s, &n)| s as u128 * n as u128).sum();
    let mean = sum as f64 / total as f64;
    let var: f64 = hist
        .iter()
        .map(|(&s, &n)| {
            let d = s as f64 - mean;
            d * d * n as f64
        })
        .sum::<f64>()
        / total as f64;

    let half = total / 2;
    let mut seen = 0u64;
    let mut median = 0u64;
    for (&size, &count) in hist {
        seen += count;
        if seen > half {
            median = size;
            break;
        }
    }

    (mean, var.sqrt(), median)
}

// ===================================================================
// Result
// ===================================================================

/// Coverage summary for one contig.
#[derive(Debug, Clone)]
pub struct ContigCoverage {
    /// Contig name from the alignment header.
    pub name: String,
    /// Contig length in bases.
    pub length: u64,
    /// Primary alignments mapped to this contig.
    pub mapped_reads: u64,
    /// Sum of aligned reference bases on this contig.
    pub mapped_bases: u64,
    /// Mean depth across the contig, uncovered bases counted as zero.
    pub mean_coverage: f64,
    /// Population standard deviation of depth across the contig.
    pub std_coverage: f64,
}

/// Everything the bamqc report needs.
#[derive(Debug)]
pub struct BamqcResult {
    /// Total reference length across all contigs.
    pub reference_bases: u64,
    /// Number of contigs in the header.
    pub num_contigs: usize,
    /// Total alignment records.
    pub total_reads: u64,
    /// Records without the unmapped flag.
    pub mapped_reads: u64,
    /// Records with the unmapped flag.
    pub unmapped_reads: u64,
    /// Secondary alignments.
    pub secondary: u64,
    /// Supplementary alignments that are not also secondary.
    pub supplementary: u64,
    /// Duplicate-flagged records.
    pub duplicates: u64,
    /// QC-failed records.
    pub qc_failed: u64,
    /// Primary mapped records with the paired flag.
    pub paired_reads: u64,
    /// Primary mapped read-1 records.
    pub first_in_pair: u64,
    /// Primary mapped read-2 records.
    pub second_in_pair: u64,
    /// Primary mapped records whose mate is also mapped.
    pub both_mates_mapped: u64,
    /// Primary mapped records whose mate is unmapped.
    pub singletons: u64,
    /// Shortest query sequence among primary mapped reads.
    pub read_length_min: u64,
    /// Longest query sequence among primary mapped reads.
    pub read_length_max: u64,
    /// Mean query sequence length among primary mapped reads.
    pub read_length_mean: f64,
    /// Total query bases of primary mapped reads.
    pub sequenced_bases: u64,
    /// Total reference bases covered by aligned blocks.
    pub aligned_bases: u64,
    /// Base counts: A, C, G, T, N.
    pub base_counts: [u64; 5],
    /// GC percentage (whole percent) → read count.
    pub gc_hist: [u64; 101],
    /// Mean mapping quality of primary mapped reads.
    pub mapq_mean: f64,
    /// Mean insert size over fragments with positive TLEN.
    pub insert_mean: f64,
    /// Population standard deviation of insert size.
    pub insert_std: f64,
    /// Median insert size.
    pub insert_median: u64,
    /// abs(TLEN) → fragment count.
    pub insert_hist: BTreeMap<u64, u64>,
    /// Mean depth across the reference, uncovered bases counted as zero.
    pub mean_coverage: f64,
    /// Population standard deviation of depth across the reference.
    pub std_coverage: f64,
    /// depth → number of reference bases at that depth (depth ≥ 1).
    pub coverage_hist: HashMap<u32, u64>,
    /// Reference bases with depth ≥ 1.
    pub covered_bases: u64,
    /// Per-contig coverage summaries, in header order.
    pub contigs: Vec<ContigCoverage>,
}

impl BamqcResult {
    /// Fraction of the reference covered at depth ≥ `threshold`, as a percentage.
    pub fn genome_fraction_at(&self, threshold: u32) -> f64 {
        if self.reference_bases == 0 {
            return 0.0;
        }
        let bases: u64 = self
            .coverage_hist
            .iter()
            .filter(|(&depth, _)| depth >= threshold)
            .map(|(_, &n)| n)
            .sum();
        100.0 * bases as f64 / self.reference_bases as f64
    }

    /// Overall GC percentage of sequenced bases.
    pub fn gc_percentage(&self) -> f64 {
        let acgt: u64 = self.base_counts[..4].iter().sum();
        if acgt == 0 {
            return 0.0;
        }
        100.0 * (self.base_counts[1] + self.base_counts[2]) as f64 / acgt as f64
    }
}

// ===================================================================
// Tests
// ===================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::Read as BamRead;

    fn accumulate(sam: &str, contigs: Vec<(String, u64)>) -> BamqcResult {
        use std::sync::atomic::{AtomicU64, Ordering};
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        let id = COUNTER.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "rustqc_bamqc_test_{:?}_{}.sam",
            std::thread::current().id(),
            id
        ));
        std::fs::write(&path, sam).unwrap();

        let mut reader = bam::Reader::from_path(&path).unwrap();
        let mut accum = BamqcAccum::new(contigs, false);
        let mut record = bam::Record::new();
        while let Some(res) = reader.read(&mut record) {
            res.unwrap();
            accum.process_read(&record);
        }
        let _ = std::fs::remove_file(&path);
        accum.into_result()
    }

    #[test]
    fn test_depth_and_genome_fraction() {
        // Contig of 100 bases. Two 10-base reads at 1 and 6 (1-based) overlap
        // over 5 bases, so: 5 bases at depth 2, 10 bases at depth 1.
        let sam = "\
@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:100\n\
r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n\
r2\t0\tchr1\t6\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n";
        let result = accumulate(sam, vec![("chr1".to_string(), 100)]);

        assert_eq!(result.coverage_hist.get(&1), Some(&10), "bases at depth 1");
        assert_eq!(result.coverage_hist.get(&2), Some(&5), "bases at depth 2");
        assert_eq!(result.covered_bases, 15);
        assert_eq!(result.aligned_bases, 20);
        // mean depth = 20 aligned bases / 100 reference bases
        assert!((result.mean_coverage - 0.2).abs() < 1e-9);
        assert!((result.genome_fraction_at(1) - 15.0).abs() < 1e-9);
        assert!((result.genome_fraction_at(2) - 5.0).abs() < 1e-9);
        assert!((result.genome_fraction_at(3) - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_refskip_does_not_count_as_covered() {
        // 5M 90N 5M spans the contig but only covers 10 bases
        let sam = "\
@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:100\n\
r1\t0\tchr1\t1\t60\t5M90N5M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n";
        let result = accumulate(sam, vec![("chr1".to_string(), 100)]);
        assert_eq!(result.covered_bases, 10);
        assert_eq!(result.aligned_bases, 10);
    }

    #[test]
    fn test_per_contig_coverage_and_flags() {
        let sam = "\
@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:100\n\
@SQ\tSN:chr2\tLN:200\n\
r1\t99\tchr1\t1\t60\t10M\t=\t51\t60\tGCGCGCGCGC\tIIIIIIIIII\n\
r2\t147\tchr1\t51\t60\t10M\t=\t1\t-60\tGCGCGCGCGC\tIIIIIIIIII\n\
r3\t256\tchr2\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n\
r4\t4\tchr2\t1\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n";
        let result = accumulate(
            sam,
            vec![("chr1".to_string(), 100), ("chr2".to_string(), 200)],
        );

        assert_eq!(result.total_reads, 4);
        assert_eq!(
            result.mapped_reads, 3,
            "secondary alignment counts as mapped"
        );
        assert_eq!(result.unmapped_reads, 1);
        assert_eq!(result.secondary, 1);
        assert_eq!(result.both_mates_mapped, 2);
        assert_eq!(result.insert_median, 60, "one fragment with TLEN 60");

        // Secondary alignments are excluded from coverage
        assert_eq!(
            result.contigs[1].mapped_bases, 0,
            "chr2 has no primary coverage"
        );
        assert_eq!(result.contigs[0].mapped_bases, 20);
        assert!((result.contigs[0].mean_coverage - 0.2).abs() < 1e-9);

        // GC: two reads of GCGCGCGCGC (100%)
        assert_eq!(result.gc_hist[100], 2);
        assert!((result.gc_percentage() - 100.0).abs() < 1e-9);
    }
}
