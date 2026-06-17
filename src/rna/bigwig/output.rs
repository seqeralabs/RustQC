//! bigWig output generation from genome coverage accumulators.
//!
//! Converts accumulated coverage to bedGraph, clips to chromosome boundaries
//! (UCSC `bedClip`), and writes bigWig files (via the `bigtools` crate).

use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};
use bigtools::beddata::BedParserStreamingIterator;
use bigtools::{BigWigWrite, Value};

use super::accumulator::{track_to_bedgraph, GenomeCovResult, TrackAccum};

/// Write all bigWig coverage tracks for a sample.
///
/// Produces nf-core/rnaseq-compatible filenames:
/// - `{sample}.bigWig` — combined (always)
/// - `{sample}.forward.bigWig` — forward strand (stranded libraries)
/// - `{sample}.reverse.bigWig` — reverse strand (stranded libraries)
pub fn write_bigwig_tracks(
    outdir: &Path,
    sample_name: &str,
    result: &GenomeCovResult,
    chrom_sizes: &[(String, u64)],
    threads: usize,
) -> Result<Vec<(String, String)>> {
    std::fs::create_dir_all(outdir).with_context(|| {
        format!(
            "Failed to create bigwig output directory: {}",
            outdir.display()
        )
    })?;

    let mut written = Vec::new();

    let combined_path = outdir.join(format!("{sample_name}.bigWig"));
    if write_track_bigwig(&combined_path, &result.combined, chrom_sizes, threads)? {
        written.push((
            "bigWig (combined)".into(),
            combined_path.display().to_string(),
        ));
    }

    if let Some(ref forward) = result.forward {
        let path = outdir.join(format!("{sample_name}.forward.bigWig"));
        if write_track_bigwig(&path, forward, chrom_sizes, threads)? {
            written.push(("bigWig (forward)".into(), path.display().to_string()));
        }
    }

    if let Some(ref reverse) = result.reverse {
        let path = outdir.join(format!("{sample_name}.reverse.bigWig"));
        if write_track_bigwig(&path, reverse, chrom_sizes, threads)? {
            written.push(("bigWig (reverse)".into(), path.display().to_string()));
        }
    }

    Ok(written)
}

/// Write a single track as a bigWig file.
///
/// Skips writing when the track has no coverage intervals (e.g. an empty
/// per-strand track), matching upstream behaviour where `bedGraphToBigWig`
/// receives an empty bedGraph. Interval ends are already clamped to each
/// chromosome's length during accumulation, so no separate `bedClip` step is
/// needed.
fn write_track_bigwig(
    path: &Path,
    track: &TrackAccum,
    chrom_sizes: &[(String, u64)],
    threads: usize,
) -> Result<bool> {
    let bedgraph = track_to_bedgraph(track, chrom_sizes);
    if bedgraph.is_empty() {
        log::debug!(
            "Skipping empty bigWig track (no coverage intervals): {}",
            path.display()
        );
        return Ok(false);
    }
    write_bedgraph_bigwig(path, bedgraph, chrom_sizes, threads)?;
    Ok(true)
}

/// Write bedGraph entries to a bigWig file using bigtools.
fn write_bedgraph_bigwig(
    path: &Path,
    entries: Vec<(String, u32, u32, u32)>,
    chrom_sizes: &[(String, u64)],
    threads: usize,
) -> Result<()> {
    let chrom_map: HashMap<String, u32> = chrom_sizes
        .iter()
        .map(|(name, size)| (name.clone(), *size as u32))
        .collect();

    let iter = entries.into_iter().map(|(chrom, start, end, val)| {
        (
            chrom,
            Value {
                start,
                end,
                value: val as f32,
            },
        )
    });
    let data = BedParserStreamingIterator::wrap_infallible_iter(iter, true);

    let runtime = if threads <= 1 {
        tokio::runtime::Builder::new_current_thread()
            .build()
            .context("Failed to create tokio runtime for bigWig writing")?
    } else {
        tokio::runtime::Builder::new_multi_thread()
            .worker_threads(threads)
            .build()
            .context("Failed to create tokio runtime for bigWig writing")?
    };

    let outb = BigWigWrite::create_file(path.to_string_lossy().to_string(), chrom_map)
        .context("Failed to create bigWig writer")?;
    outb.write(data, runtime)
        .map_err(|e| anyhow::anyhow!("Failed to write bigWig to {}: {e}", path.display()))?;

    Ok(())
}
