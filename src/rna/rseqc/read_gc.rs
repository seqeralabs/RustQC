//! GC content distribution of mapped reads.
//!
//! Reimplementation of RSeQC's `read_GC.py`. For every read passing the
//! filters, the GC percentage of the query sequence is computed and binned;
//! the distribution is the standard check for GC bias introduced by library
//! preparation or sequencing.

use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use log::debug;
use rust_htslib::bam;

use crate::rna::bam_flags::{BAM_FQCFAIL, BAM_FUNMAP};

/// Number of GC-percentage buckets: one per 0.01% step over [0, 100].
///
/// RSeQC formats the percentage with `"%4.2f"`, so two decimal places is the
/// full resolution of its output keys.
const GC_BUCKETS: usize = 10_001;

// ===================================================================
// Accumulator
// ===================================================================

/// Per-worker accumulator for the read GC distribution.
///
/// Reads are filtered exactly as RSeQC's `ParseBAM.readGC` does: unmapped and
/// QC-failed records are skipped, along with records below the MAPQ cutoff.
/// Secondary alignments and duplicates are deliberately *not* filtered —
/// upstream does not filter them either.
#[derive(Debug, Clone)]
pub struct ReadGcAccum {
    /// Counts indexed by `round(gc_percent * 100)`.
    counts: Vec<u64>,
    /// MAPQ cutoff below which reads are ignored.
    mapq_cut: u8,
    /// Reads skipped because they carried no query sequence.
    pub missing_sequence: u64,
}

impl ReadGcAccum {
    /// Create an accumulator with the given MAPQ cutoff.
    pub fn new(mapq_cut: u8) -> Self {
        Self {
            counts: vec![0; GC_BUCKETS],
            mapq_cut,
            missing_sequence: 0,
        }
    }

    /// Process a single alignment record.
    pub fn process_read(&mut self, record: &bam::Record) {
        let flags = record.flags();
        if flags & (BAM_FUNMAP | BAM_FQCFAIL) != 0 {
            return;
        }
        if record.mapq() < self.mapq_cut {
            return;
        }

        let seq = record.seq();
        let len = seq.len();
        if len == 0 {
            // Records stored without SEQ ('*'). Upstream would divide by zero
            // here; we skip them and report the count instead.
            self.missing_sequence += 1;
            return;
        }

        let mut gc = 0u64;
        for i in 0..len {
            // seq() decodes to uppercase IUPAC bytes
            if matches!(seq[i], b'C' | b'G' | b'c' | b'g') {
                gc += 1;
            }
        }

        let pct = 100.0 * gc as f64 / len as f64;
        let bucket = (pct * 100.0).round() as usize;
        // pct <= 100 so bucket <= 10000, but clamp defensively
        self.counts[bucket.min(GC_BUCKETS - 1)] += 1;
    }

    /// Merge another worker's counts into this accumulator.
    pub fn merge(&mut self, other: ReadGcAccum) {
        for (a, b) in self.counts.iter_mut().zip(other.counts.iter()) {
            *a += b;
        }
        self.missing_sequence += other.missing_sequence;
    }

    /// Finalise into a sorted distribution.
    pub fn into_result(self) -> ReadGcResult {
        let distribution: Vec<(f64, u64)> = self
            .counts
            .iter()
            .enumerate()
            .filter(|(_, &count)| count > 0)
            .map(|(bucket, &count)| (bucket as f64 / 100.0, count))
            .collect();
        ReadGcResult {
            distribution,
            missing_sequence: self.missing_sequence,
        }
    }
}

/// GC distribution over all processed reads.
#[derive(Debug, Clone)]
pub struct ReadGcResult {
    /// `(gc_percent, read_count)` pairs, ascending by GC percentage.
    ///
    /// Only non-empty buckets are kept, matching upstream's sparse dictionary.
    pub distribution: Vec<(f64, u64)>,
    /// Reads skipped because they carried no query sequence.
    pub missing_sequence: u64,
}

impl ReadGcResult {
    /// Total reads contributing to the distribution.
    pub fn total_reads(&self) -> u64 {
        self.distribution.iter().map(|&(_, c)| c).sum()
    }

    /// Mean GC percentage across all reads, or `None` when empty.
    pub fn mean_gc(&self) -> Option<f64> {
        let total = self.total_reads();
        if total == 0 {
            return None;
        }
        let sum: f64 = self
            .distribution
            .iter()
            .map(|&(pct, count)| pct * count as f64)
            .sum();
        Some(sum / total as f64)
    }

    /// Distribution as a `HashMap` keyed by the formatted percentage, matching
    /// upstream's output keys.
    pub fn as_map(&self) -> HashMap<String, u64> {
        self.distribution
            .iter()
            .map(|&(pct, count)| (format!("{pct:.2}"), count))
            .collect()
    }
}

// ===================================================================
// Output
// ===================================================================

/// Write the `.GC.xls` distribution table.
///
/// Format matches RSeQC: a `GC%\tread_count` header followed by one row per
/// observed GC percentage. Rows are sorted ascending by GC percentage
/// (upstream emits Python dictionary order, which is BAM-order dependent and
/// therefore not reproducible under parallel processing).
pub fn write_gc_table(result: &ReadGcResult, output_path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(output_path)
        .with_context(|| format!("Failed to create GC table: {}", output_path.display()))?;

    writeln!(out, "GC%\tread_count")?;
    for &(pct, count) in &result.distribution {
        writeln!(out, "{pct:.2}\t{count}")?;
    }

    debug!(
        "Wrote read GC distribution to {} ({} buckets)",
        output_path.display(),
        result.distribution.len()
    );
    Ok(())
}

/// Write the R plotting script RSeQC produces alongside the table.
///
/// RustQC renders the plot itself; the script is written for drop-in
/// compatibility with pipelines that expect it.
pub fn write_gc_r_script(
    result: &ReadGcResult,
    output_prefix: &str,
    output_path: &Path,
) -> Result<()> {
    let mut out = std::fs::File::create(output_path)
        .with_context(|| format!("Failed to create GC R script: {}", output_path.display()))?;

    // Emit the distribution as rep(values, times=counts) rather than one entry
    // per read: identical to upstream's vector, but linear in the number of
    // distinct GC values instead of the number of reads.
    writeln!(out, "pdf(\"{output_prefix}.GC_plot.pdf\")")?;
    writeln!(
        out,
        "gc=rep(c({}),times=c({}))",
        result
            .distribution
            .iter()
            .map(|&(pct, _)| format!("{pct:.2}"))
            .collect::<Vec<_>>()
            .join(","),
        result
            .distribution
            .iter()
            .map(|&(_, count)| count.to_string())
            .collect::<Vec<_>>()
            .join(","),
    )?;
    writeln!(out, "hist(gc,probability=T,breaks=100,xlab=\"GC content (%)\",ylab=\"Density of Reads\",border=\"blue\",main=\"\")")?;
    writeln!(out, "lines(density(gc),col='red')")?;
    writeln!(out, "dev.off()")?;

    debug!("Wrote read GC R script to {}", output_path.display());
    Ok(())
}

// ===================================================================
// Tests
// ===================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::Read as BamRead;

    fn write_temp_sam(content: &str) -> std::path::PathBuf {
        use std::sync::atomic::{AtomicU64, Ordering};
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        let id = COUNTER.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "rustqc_read_gc_test_{:?}_{}.sam",
            std::thread::current().id(),
            id
        ));
        std::fs::write(&path, content).unwrap();
        path
    }

    fn accumulate(sam: &str, mapq_cut: u8) -> ReadGcResult {
        let path = write_temp_sam(sam);
        let mut reader = bam::Reader::from_path(&path).unwrap();
        let mut accum = ReadGcAccum::new(mapq_cut);
        let mut record = bam::Record::new();
        while let Some(res) = reader.read(&mut record) {
            res.unwrap();
            accum.process_read(&record);
        }
        let _ = std::fs::remove_file(&path);
        accum.into_result()
    }

    #[test]
    fn test_gc_percentages_and_filters() {
        // r1: GCGCGCGCGC -> 100.00
        // r2: ACGTACGTAC -> 50.00
        // r3: AAAAAAAAAA -> 0.00
        // r4: MAPQ 5, below the cutoff -> skipped
        // r5: QC-fail (0x200) -> skipped
        // r6: unmapped (0x4) -> skipped
        let sam = "\
@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:20000\n\
r1\t0\tchr1\t100\t60\t10M\t*\t0\t0\tGCGCGCGCGC\tIIIIIIIIII\n\
r2\t0\tchr1\t120\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n\
r3\t0\tchr1\t140\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n\
r4\t0\tchr1\t160\t5\t10M\t*\t0\t0\tGCGCGCGCGC\tIIIIIIIIII\n\
r5\t512\tchr1\t180\t60\t10M\t*\t0\t0\tGCGCGCGCGC\tIIIIIIIIII\n\
r6\t4\tchr1\t200\t0\t*\t*\t0\t0\tGCGCGCGCGC\tIIIIIIIIII\n";

        let result = accumulate(sam, 30);
        assert_eq!(result.total_reads(), 3, "three reads pass the filters");

        let map = result.as_map();
        assert_eq!(map.get("100.00"), Some(&1));
        assert_eq!(map.get("50.00"), Some(&1));
        assert_eq!(map.get("0.00"), Some(&1));
        assert_eq!(map.len(), 3);

        // Ascending order
        let pcts: Vec<f64> = result.distribution.iter().map(|&(p, _)| p).collect();
        assert_eq!(pcts, vec![0.0, 50.0, 100.0]);

        let mean = result.mean_gc().unwrap();
        assert!((mean - 50.0).abs() < 1e-9, "mean GC: {mean}");
    }

    #[test]
    fn test_secondary_and_duplicate_reads_are_counted() {
        // Upstream readGC does not filter secondary alignments or duplicates.
        // r1 secondary (0x100), r2 duplicate (0x400)
        let sam = "\
@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:20000\n\
r1\t256\tchr1\t100\t60\t10M\t*\t0\t0\tGCGCGCGCGC\tIIIIIIIIII\n\
r2\t1024\tchr1\t120\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n";

        let result = accumulate(sam, 30);
        assert_eq!(result.total_reads(), 2);
    }

    #[test]
    fn test_reads_without_sequence_are_reported() {
        let sam = "\
@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:20000\n\
r1\t0\tchr1\t100\t60\t10M\t*\t0\t0\t*\t*\n";

        let result = accumulate(sam, 30);
        assert_eq!(result.total_reads(), 0);
        assert_eq!(result.missing_sequence, 1);
    }

    #[test]
    fn test_merge_combines_distributions() {
        let mut a = ReadGcAccum::new(0);
        let mut b = ReadGcAccum::new(0);
        a.counts[5000] = 3;
        b.counts[5000] = 4;
        b.counts[10000] = 1;
        a.merge(b);
        let result = a.into_result();
        assert_eq!(result.distribution, vec![(50.0, 7), (100.0, 1)]);
    }
}
