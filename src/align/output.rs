//! mosdepth- and NGSCheckMate-compatible output files.

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use log::debug;

use super::depth::{ContigDepth, Window, MAX_DIST_DEPTH};
use super::snp::{SnpAccum, SnpPanel};

/// Write `{prefix}.mosdepth.summary.txt`.
///
/// Columns are mosdepth's: `chrom length bases mean min max`, one row per
/// contig plus a `total` row. Because `--by <N>` tiles the whole reference,
/// each contig also gets the `<chrom>_region` row mosdepth emits for the
/// region set, with the same values.
pub fn write_summary(contigs: &[ContigDepth], path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .with_context(|| format!("Failed to create {}", path.display()))?;
    writeln!(out, "chrom\tlength\tbases\tmean\tmin\tmax")?;

    let mut total_length = 0u64;
    let mut total_bases = 0u64;
    let mut total_max = 0u32;
    for contig in contigs {
        for name in [contig.name.clone(), format!("{}_region", contig.name)] {
            writeln!(
                out,
                "{}\t{}\t{}\t{:.2}\t{}\t{}",
                name,
                contig.length,
                contig.total_bases,
                contig.mean(),
                contig.min_depth,
                contig.max_depth
            )?;
        }
        total_length += contig.length;
        total_bases += contig.total_bases;
        total_max = total_max.max(contig.max_depth);
    }

    let total_mean = if total_length == 0 {
        0.0
    } else {
        total_bases as f64 / total_length as f64
    };
    for name in ["total", "total_region"] {
        writeln!(
            out,
            "{}\t{}\t{}\t{:.2}\t{}\t{}",
            name, total_length, total_bases, total_mean, 0, total_max
        )?;
    }

    debug!("Wrote mosdepth summary to {}", path.display());
    Ok(())
}

/// Write `{prefix}.mosdepth.global.dist.txt`.
///
/// For each contig (and `total`), emits `chrom depth proportion` from the
/// deepest observed level down to 0, where `proportion` is the fraction of
/// that contig's bases covered at **at least** that depth. mosdepth's format.
pub fn write_global_dist(contigs: &[ContigDepth], path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .with_context(|| format!("Failed to create {}", path.display()))?;

    let mut total_hist: std::collections::BTreeMap<u32, u64> = std::collections::BTreeMap::new();
    let mut total_length = 0u64;

    for contig in contigs {
        write_dist_block(&mut out, &contig.name, &contig.hist, contig.length)?;
        for (&depth, &bases) in &contig.hist {
            *total_hist.entry(depth).or_insert(0) += bases;
        }
        total_length += contig.length;
    }
    write_dist_block(&mut out, "total", &total_hist, total_length)?;

    debug!("Wrote mosdepth global distribution to {}", path.display());
    Ok(())
}

/// Emit the cumulative distribution rows for one contig.
fn write_dist_block(
    out: &mut impl Write,
    name: &str,
    hist: &std::collections::BTreeMap<u32, u64>,
    length: u64,
) -> Result<()> {
    if length == 0 {
        return Ok(());
    }
    let max_depth = hist
        .keys()
        .next_back()
        .copied()
        .unwrap_or(0)
        .min(MAX_DIST_DEPTH);
    let mut cumulative = 0u64;
    for depth in (0..=max_depth).rev() {
        cumulative += hist.get(&depth).copied().unwrap_or(0);
        let proportion = if depth == 0 {
            1.0
        } else {
            cumulative as f64 / length as f64
        };
        // mosdepth skips the sparse tail at the top of the distribution
        // (`if cum < 8e-5: continue` in its write_distribution)
        if proportion < 8e-5 {
            continue;
        }
        writeln!(out, "{name}\t{depth}\t{proportion:.2}")?;
    }
    Ok(())
}

/// Write `{prefix}.regions.bed.gz` (BGZF-compressed per-window depth).
pub fn write_regions(windows: &[Window], path: &Path) -> Result<()> {
    use rust_htslib::bgzf;

    let mut writer = bgzf::Writer::from_path(path)
        .with_context(|| format!("Failed to create {}", path.display()))?;
    for window in windows {
        let line = format!(
            "{}\t{}\t{}\t{:.2}\n",
            window.chrom, window.start, window.end, window.mean
        );
        writer
            .write_all(line.as_bytes())
            .with_context(|| format!("Failed to write {}", path.display()))?;
    }

    debug!(
        "Wrote {} depth windows to {}",
        windows.len(),
        path.display()
    );
    Ok(())
}

/// Write the NGSCheckMate VCF (BGZF-compressed).
///
/// Emits `GT:AD:DP` per site, with the genotype derived from the alternate
/// allele fraction. `ncm.py` only needs to tell hom-ref, het and hom-alt
/// apart, so no genotype likelihood model is involved.
pub fn write_ngscheckmate_vcf(
    panel: &SnpPanel,
    accum: &SnpAccum,
    sample_name: &str,
    contig_order: &[(String, u64)],
    path: &Path,
) -> Result<()> {
    use rust_htslib::bgzf;

    let mut writer = bgzf::Writer::from_path(path)
        .with_context(|| format!("Failed to create {}", path.display()))?;

    let mut header = String::new();
    header.push_str("##fileformat=VCFv4.2\n");
    header.push_str("##source=rustqc align\n");
    for (name, length) in contig_order {
        header.push_str(&format!("##contig=<ID={name},length={length}>\n"));
    }
    header.push_str("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">\n");
    header.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    header.push_str(
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles\">\n",
    );
    header.push_str("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n");
    header.push_str(&format!(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n"
    ));
    writer
        .write_all(header.as_bytes())
        .with_context(|| format!("Failed to write {}", path.display()))?;

    // Emit in header contig order, positions ascending, so the file is sorted
    for (chrom, _) in contig_order {
        let (Some(sites), Some(counts)) = (panel.by_chrom.get(chrom), accum.counts_for(chrom))
        else {
            continue;
        };
        for (site, count) in sites.iter().zip(counts.iter()) {
            let line = format!(
                "{}\t{}\t{}\t{}\t{}\t.\t.\tDP={}\tGT:AD:DP\t{}:{},{}:{}\n",
                chrom,
                site.pos + 1,
                site.id,
                site.ref_base as char,
                site.alt_base as char,
                count.depth(),
                count.genotype(),
                count.ref_count,
                count.alt_count,
                count.depth(),
            );
            writer
                .write_all(line.as_bytes())
                .with_context(|| format!("Failed to write {}", path.display()))?;
        }
    }

    debug!("Wrote NGSCheckMate VCF to {}", path.display());
    Ok(())
}

// ===================================================================
// Tests
// ===================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    #[test]
    fn test_dist_block_is_cumulative() {
        // 100-base contig: 10 bases at depth 1, 5 at depth 2
        let mut hist = BTreeMap::new();
        hist.insert(1u32, 10u64);
        hist.insert(2u32, 5u64);

        let mut out: Vec<u8> = Vec::new();
        write_dist_block(&mut out, "chr1", &hist, 100).unwrap();
        let text = String::from_utf8(out).unwrap();

        // depth 2 -> 5/100, depth 1 -> 15/100, depth 0 -> 1.00
        assert_eq!(text, "chr1\t2\t0.05\nchr1\t1\t0.15\nchr1\t0\t1.00\n");
    }

    #[test]
    fn test_summary_has_total_row() {
        let contigs = vec![
            ContigDepth {
                name: "chr1".to_string(),
                length: 100,
                total_bases: 20,
                min_depth: 0,
                max_depth: 2,
                hist: BTreeMap::new(),
            },
            ContigDepth {
                name: "chr2".to_string(),
                length: 100,
                total_bases: 0,
                min_depth: 0,
                max_depth: 0,
                hist: BTreeMap::new(),
            },
        ];
        let path = std::env::temp_dir().join(format!(
            "rustqc_summary_test_{:?}.txt",
            std::thread::current().id()
        ));
        write_summary(&contigs, &path).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let _ = std::fs::remove_file(&path);

        assert!(text.starts_with("chrom\tlength\tbases\tmean\tmin\tmax\n"));
        assert!(text.contains("chr1\t100\t20\t0.20\t0\t2\n"));
        assert!(text.contains("chr1_region\t100\t20\t0.20\t0\t2\n"));
        assert!(text.contains("total\t200\t20\t0.10\t0\t2\n"));
        assert!(text.contains("total_region\t200\t20\t0.10\t0\t2\n"));
    }
}
