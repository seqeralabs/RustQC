//! Writers for the six mosdepth-compatible output files.
//!
//! Formats are documented in the parent module. Compressed outputs are written
//! as bgzf, which is what mosdepth writes and what both `tabix` and `gunzip`
//! read. Parity against the fixtures is therefore asserted on the decompressed
//! bytes: two bgzf writers at the same level need not emit identical
//! compressed bytes, so comparing the `.gz` byte for byte would be testing the
//! compressor rather than this code.

use std::io::Write;
use std::path::Path;

use anyhow::{bail, Context, Result};
use rust_htslib::bgzf;

use super::{dist_proportions, dist_rows, merge_histograms, MosdepthResult};

/// Write `{prefix}.mosdepth.summary.txt`.
pub fn write_summary(result: &MosdepthResult, path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .map(std::io::BufWriter::new)
        .with_context(|| format!("Failed to create summary file: {}", path.display()))?;

    writeln!(out, "chrom\tlength\tbases\tmean\tmin\tmax")?;
    for contig in &result.contigs {
        writeln!(
            out,
            "{}\t{}\t{}\t{:.2}\t{}\t{}",
            contig.name,
            contig.length,
            contig.total_bases,
            contig.mean(),
            contig.min,
            contig.max
        )?;
        if result.window_size.is_some() {
            writeln!(
                out,
                "{}_region\t{}\t{}\t{:.2}\t{}\t{}",
                contig.name,
                contig.length,
                contig.total_bases,
                contig.mean(),
                contig.min,
                contig.max
            )?;
        }
    }
    writeln!(
        out,
        "total\t{}\t{}\t{:.2}\t{}\t{}",
        result.total_length(),
        result.total_bases(),
        result.mean(),
        result.min(),
        result.max()
    )?;
    if result.window_size.is_some() {
        writeln!(
            out,
            "total_region\t{}\t{}\t{:.2}\t{}\t{}",
            result.total_length(),
            result.total_bases(),
            result.mean(),
            result.min(),
            result.max()
        )?;
    }
    out.flush()?;
    Ok(())
}

/// Write `{prefix}.mosdepth.global.dist.txt`, the distribution over bases.
pub fn write_global_dist(result: &MosdepthResult, path: &Path) -> Result<()> {
    let per_contig: Vec<_> = result
        .contigs
        .iter()
        .map(|c| (c.name.as_str(), c.histogram.clone(), c.length))
        .collect();
    write_dist(&per_contig, path)
}

/// Write `{prefix}.mosdepth.region.dist.txt`, the distribution over windows
/// and their rounded mean depth.
pub fn write_region_dist(result: &MosdepthResult, path: &Path) -> Result<()> {
    let per_contig: Vec<_> = result
        .contigs
        .iter()
        .map(|c| {
            let hist = c.region_histogram();
            let total = hist.values().sum::<u64>();
            (c.name.as_str(), hist, total)
        })
        .collect();
    write_dist(&per_contig, path)
}

/// Shared body of both distribution writers.
fn write_dist(
    per_contig: &[(&str, std::collections::BTreeMap<u32, u64>, u64)],
    path: &Path,
) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .map(std::io::BufWriter::new)
        .with_context(|| format!("Failed to create distribution file: {}", path.display()))?;

    for (name, histogram, total) in per_contig {
        let rows = dist_rows(histogram);
        for (depth, proportion) in rows.iter().zip(dist_proportions(histogram, &rows, *total)) {
            writeln!(out, "{name}\t{depth}\t{proportion:.2}")?;
        }
    }

    let merged = merge_histograms(per_contig.iter().map(|(_, h, _)| h));
    let total: u64 = per_contig.iter().map(|(_, _, t)| t).sum();
    let rows = dist_rows(&merged);
    for (depth, proportion) in rows.iter().zip(dist_proportions(&merged, &rows, total)) {
        writeln!(out, "total\t{depth}\t{proportion:.2}")?;
    }

    out.flush()?;
    Ok(())
}

/// Write `{prefix}.per-base.bed.gz`, one line per run of equal depth.
pub fn write_per_base(result: &MosdepthResult, path: &Path) -> Result<()> {
    let mut lines = Vec::new();
    for contig in &result.contigs {
        for run in &contig.runs {
            lines.push(format!(
                "{}\t{}\t{}\t{}\n",
                contig.name, run.start, run.end, run.depth
            ));
        }
    }
    write_bgzf(path, &lines.concat())
}

/// Write `{prefix}.regions.bed.gz`, one line per window with its mean depth.
pub fn write_regions(result: &MosdepthResult, path: &Path) -> Result<()> {
    let mut lines = Vec::new();
    for contig in &result.contigs {
        for window in &contig.windows {
            lines.push(format!(
                "{}\t{}\t{}\t{:.2}\n",
                contig.name, window.start, window.end, window.mean
            ));
        }
    }
    write_bgzf(path, &lines.concat())
}

/// Write `{prefix}.thresholds.bed.gz`, one line per window with the number of
/// bases at or above each requested threshold.
pub fn write_thresholds(result: &MosdepthResult, path: &Path) -> Result<()> {
    let mut body = String::from("#chrom\tstart\tend\tregion");
    for threshold in &result.thresholds {
        body.push_str(&format!("\t{threshold}X"));
    }
    body.push('\n');

    for contig in &result.contigs {
        for row in &contig.thresholds {
            body.push_str(&format!(
                "{}\t{}\t{}\tunknown",
                contig.name, row.start, row.end
            ));
            for count in &row.counts {
                body.push_str(&format!("\t{count}"));
            }
            body.push('\n');
        }
    }
    write_bgzf(path, &body)
}

/// Write `contents` to `path` as bgzf, then build its `.csi` index.
fn write_bgzf(path: &Path, contents: &str) -> Result<()> {
    {
        let mut writer = bgzf::Writer::from_path(path)
            .with_context(|| format!("Failed to create bgzf file: {}", path.display()))?;
        writer
            .write_all(contents.as_bytes())
            .with_context(|| format!("Failed to write bgzf file: {}", path.display()))?;
        // The writer must be dropped, and the bgzf stream closed, before the
        // indexer reads the file back.
    }
    build_csi_index(path)
}

/// Build the `.csi` companion index for a bgzf-compressed BED file.
///
/// mosdepth writes one alongside each of its BED outputs, and `tabix` needs it
/// to seek into them. CSI rather than TBI because CSI carries no 512 Mb
/// coordinate ceiling, which matters on large contigs.
fn build_csi_index(path: &Path) -> Result<()> {
    use std::ffi::CString;

    let path_c = CString::new(path.as_os_str().as_encoded_bytes()).with_context(|| {
        format!(
            "Path is not representable as a C string: {}",
            path.display()
        )
    })?;

    // SAFETY: `path_c` is a valid NUL-terminated string that outlives the
    // call, `tbx_conf_bed` is a static provided by htslib, and the file was
    // closed above. A min_shift of 14 selects CSI, matching what mosdepth and
    // `tabix --csi` produce.
    let ret = unsafe {
        rust_htslib::htslib::tbx_index_build(
            path_c.as_ptr(),
            14,
            &raw const rust_htslib::htslib::tbx_conf_bed,
        )
    };
    if ret < 0 {
        bail!("Failed to build the CSI index for {}", path.display());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna::mosdepth::ContigDepth;
    use std::io::Read;

    fn scratch(name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join("rustqc-mosdepth-tests");
        std::fs::create_dir_all(&dir).unwrap();
        dir.join(name)
    }

    fn result_with_windows() -> MosdepthResult {
        let depths = vec![0u32, 0, 2, 2, 4, 4];
        MosdepthResult {
            contigs: vec![ContigDepth::from_depths("chr1", &depths, Some(3), &[1, 4])],
            window_size: Some(3),
            thresholds: vec![1, 4],
        }
    }

    fn read_bgzf(path: &std::path::Path) -> String {
        let mut reader = bgzf::Reader::from_path(path).unwrap();
        let mut buf = Vec::new();
        reader.read_to_end(&mut buf).unwrap();
        String::from_utf8(buf).unwrap()
    }

    #[test]
    fn summary_has_region_rows_only_when_windows_were_requested() {
        let path = scratch("summary_windows.txt");
        write_summary(&result_with_windows(), &path).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert_eq!(
            text,
            "chrom\tlength\tbases\tmean\tmin\tmax\n\
             chr1\t6\t12\t2.00\t0\t4\n\
             chr1_region\t6\t12\t2.00\t0\t4\n\
             total\t6\t12\t2.00\t0\t4\n\
             total_region\t6\t12\t2.00\t0\t4\n"
        );

        let depths = vec![0u32, 0, 2, 2, 4, 4];
        let no_windows = MosdepthResult {
            contigs: vec![ContigDepth::from_depths("chr1", &depths, None, &[])],
            window_size: None,
            thresholds: vec![],
        };
        let path = scratch("summary_nowindows.txt");
        write_summary(&no_windows, &path).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(!text.contains("_region"), "no windows means no region rows");
    }

    #[test]
    fn per_base_writes_one_line_per_run() {
        let path = scratch("per-base.bed.gz");
        write_per_base(&result_with_windows(), &path).unwrap();
        assert_eq!(
            read_bgzf(&path),
            "chr1\t0\t2\t0\nchr1\t2\t4\t2\nchr1\t4\t6\t4\n"
        );
    }

    #[test]
    fn regions_carry_two_decimal_means() {
        let path = scratch("regions.bed.gz");
        write_regions(&result_with_windows(), &path).unwrap();
        assert_eq!(read_bgzf(&path), "chr1\t0\t3\t0.67\nchr1\t3\t6\t3.33\n");
    }

    #[test]
    fn thresholds_carry_a_header_and_one_column_per_threshold() {
        let path = scratch("thresholds.bed.gz");
        write_thresholds(&result_with_windows(), &path).unwrap();
        assert_eq!(
            read_bgzf(&path),
            "#chrom\tstart\tend\tregion\t1X\t4X\n\
             chr1\t0\t3\tunknown\t1\t0\n\
             chr1\t3\t6\tunknown\t3\t2\n"
        );
    }

    #[test]
    fn global_dist_is_descending_and_ends_at_one() {
        let path = scratch("global.dist.txt");
        write_global_dist(&result_with_windows(), &path).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let chr1: Vec<&str> = text.lines().filter(|l| l.starts_with("chr1\t")).collect();
        assert_eq!(
            *chr1.first().unwrap(),
            "chr1\t4\t0.33",
            "descending from the maximum"
        );
        assert_eq!(*chr1.last().unwrap(), "chr1\t0\t1.00", "down to zero");
        assert_eq!(
            chr1.len(),
            5,
            "depths 4 down to 0, all inside the dense range"
        );
        assert!(text.contains("total\t0\t1.00"));
    }

    #[test]
    fn compressed_outputs_get_a_loadable_csi_index() {
        let path = scratch("indexed.per-base.bed.gz");
        let index = scratch("indexed.per-base.bed.gz.csi");
        let _ = std::fs::remove_file(&index);
        write_per_base(&result_with_windows(), &path).unwrap();
        assert!(index.exists(), "the .csi companion index must be written");
        // htslib refuses to open a malformed index, so opening it is the check.
        let tbx = rust_htslib::tbx::Reader::from_path(&path);
        assert!(tbx.is_ok(), "htslib could not open the indexed file");
    }
}
