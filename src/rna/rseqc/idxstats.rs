//! Samtools idxstats-compatible output.
//!
//! Produces output matching `samtools idxstats` format: a tab-delimited table
//! with one row per reference sequence plus a final `*` row for unplaced reads.

use std::path::Path;

use anyhow::{Context, Result};
use log::info;

use super::bam_stat::BamStatResult;

// ============================================================================
// Output formatting
// ============================================================================

/// Write samtools idxstats-compatible output.
///
/// Format: `ref_name\tseq_length\tmapped_count\tunmapped_count` per line,
/// with a final `*\t0\t0\t<unplaced_unmapped>` line.
///
/// The reference names and lengths come from the BAM header, so all references
/// are included even if they have zero reads.
///
/// # Arguments
/// * `result` - The computed BAM statistics
/// * `header_refs` - Reference names and lengths from the BAM header, as `(name, length)` pairs
/// * `output_path` - Path to write the idxstats file
pub fn write_idxstats(
    result: &BamStatResult,
    header_refs: &[(String, u64)],
    output_path: &Path,
) -> Result<()> {
    use std::io::Write;

    let mut out = std::fs::File::create(output_path)
        .with_context(|| format!("Failed to create idxstats file: {}", output_path.display()))?;

    // One line per reference from the BAM header
    for (tid, (name, length)) in header_refs.iter().enumerate() {
        let (mapped, unmapped) = result
            .chrom_counts
            .get(&(tid as i32))
            .copied()
            .unwrap_or((0, 0));
        writeln!(out, "{name}\t{length}\t{mapped}\t{unmapped}")?;
    }

    // Final line for unplaced unmapped reads
    writeln!(out, "*\t0\t0\t{}", result.unplaced_unmapped)?;

    info!("Wrote idxstats output to {}", output_path.display());
    Ok(())
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rna::rseqc::accumulators::BamStatAccum;
    use rust_htslib::bam::{self, Read as BamRead};
    use std::io::Read;
    #[test]
    fn test_idxstats_format() {
        let mut reader =
            bam::Reader::from_path("tests/data/test.bam").expect("Failed to open test.bam");

        // Extract header reference info
        let header = reader.header().clone();
        let header_refs: Vec<(String, u64)> = (0..header.target_count())
            .map(|tid| {
                let name = String::from_utf8_lossy(header.tid2name(tid)).to_string();
                let len = header.target_len(tid).unwrap_or(0) as u64;
                (name, len)
            })
            .collect();

        let mut accum = BamStatAccum::default();
        let mut record = bam::Record::new();

        while let Some(res) = reader.read(&mut record) {
            res.expect("Error reading BAM record");
            accum.process_read(&record, 30);
        }

        let result = accum.into_result();
        let tmp_path = std::env::temp_dir().join("rustqc_test_idxstats.txt");
        write_idxstats(&result, &header_refs, &tmp_path).expect("Failed to write idxstats");

        let mut contents = String::new();
        std::fs::File::open(&tmp_path)
            .unwrap()
            .read_to_string(&mut contents)
            .unwrap();
        let _ = std::fs::remove_file(&tmp_path);

        // Verify last line is the unplaced unmapped row
        let last_line = contents.lines().last().expect("Empty output");
        assert!(
            last_line.starts_with("*\t0\t0\t"),
            "Last line should be the unplaced unmapped row, got: {last_line}"
        );

        // Verify format: 4 tab-separated columns per line
        for line in contents.lines() {
            let cols: Vec<&str> = line.split('\t').collect();
            assert_eq!(
                cols.len(),
                4,
                "Expected 4 columns in idxstats line, got {}: {line}",
                cols.len()
            );
        }
    }
}
