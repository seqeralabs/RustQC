//! bigWig writing for coverage tracks.
//!
//! Replaces the `bedtools genomecov` into `bedGraphToBigWig` round-trip with a
//! single write from the coverage already computed in the alignment pass.
//!
//! Behind the `bigwig` cargo feature, so a build that only wants the text
//! outputs need not carry `bigtools`.

use std::collections::HashMap;
use std::path::Path;

use anyhow::{anyhow, Context, Result};
use bigtools::beddata::BedParserStreamingIterator;
use bigtools::{BigWigWrite, Value};

use super::bedgraph::Interval;

/// Write intervals as a bigWig, returning whether a file was produced.
///
/// `chrom_sizes` must name every contig the intervals refer to, and comes from
/// the alignment header. Intervals are expected in the order the contigs
/// appear there.
///
/// An empty interval list writes nothing and returns `false`. The format has
/// no representation for a track with no data, and `bigtools` rejects it
/// outright, so the alternative would be failing a run because an alignment
/// covered nothing.
///
/// Note that a gap between intervals is not zero coverage in a bigWig: it is
/// *undefined*, and reads back as `NaN`. That is the format's own semantics
/// and matches what `bedGraphToBigWig` produces from a bedGraph that omits its
/// zero-depth spans, which is what `bedtools genomecov -bg` emits.
pub fn write_bigwig(
    intervals: &[Interval],
    chrom_sizes: &[(String, u64)],
    path: &Path,
) -> Result<bool> {
    if intervals.is_empty() {
        return Ok(false);
    }
    let sizes: HashMap<String, u32> = chrom_sizes
        .iter()
        .map(|(name, length)| (name.clone(), *length as u32))
        .collect();

    let values: Vec<(String, Value)> = intervals
        .iter()
        .map(|interval| {
            (
                interval.chrom.clone(),
                Value {
                    start: interval.start,
                    end: interval.end,
                    value: interval.value,
                },
            )
        })
        .collect();

    // `false` because the intervals are already grouped by contig in header
    // order; letting the writer accept out-of-order chromosomes would hide a
    // bug in the caller rather than catching it.
    let data = BedParserStreamingIterator::wrap_infallible_iter(values.into_iter(), false);

    let writer = BigWigWrite::create_file(path, sizes)
        .with_context(|| format!("Failed to create bigWig: {}", path.display()))?;

    // One worker: the encoding is not the bottleneck next to reading the
    // alignment, and a fixed thread count keeps the output deterministic.
    let runtime = tokio::runtime::Builder::new_multi_thread()
        .worker_threads(1)
        .build()
        .context("Failed to start the bigWig writer runtime")?;

    writer
        .write(data, runtime)
        .map_err(|e| anyhow!("Failed to write bigWig {}: {e}", path.display()))?;
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scratch(name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join("rustqc-bigwig-tests");
        std::fs::create_dir_all(&dir).unwrap();
        dir.join(name)
    }

    #[test]
    fn a_written_track_reads_back_with_the_same_values() {
        use bigtools::BigWigRead;

        let intervals = vec![
            Interval {
                chrom: "chr1".into(),
                start: 10,
                end: 20,
                value: 3.0,
            },
            Interval {
                chrom: "chr1".into(),
                start: 30,
                end: 35,
                value: 7.5,
            },
        ];
        let path = scratch("roundtrip.bw");
        assert!(write_bigwig(&intervals, &[("chr1".to_string(), 100)], &path).unwrap());

        let mut reader = BigWigRead::open_file(&path).unwrap();
        let values = reader.values("chr1", 0, 100).unwrap();
        assert_eq!(values[15], 3.0, "inside the first interval");
        assert_eq!(values[32], 7.5, "inside the second");
    }

    #[test]
    fn a_gap_is_undefined_rather_than_zero() {
        use bigtools::BigWigRead;

        let intervals = vec![Interval {
            chrom: "chr1".into(),
            start: 10,
            end: 20,
            value: 3.0,
        }];
        let path = scratch("gap.bw");
        write_bigwig(&intervals, &[("chr1".to_string(), 100)], &path).unwrap();

        let mut reader = BigWigRead::open_file(&path).unwrap();
        let values = reader.values("chr1", 0, 100).unwrap();
        assert!(
            values[50].is_nan(),
            "a bigWig has no representation for zero coverage; the gap is \
             undefined and reads back as NaN, got {}",
            values[50]
        );
    }

    #[test]
    fn nothing_to_write_produces_no_file_rather_than_failing() {
        let path = scratch("empty.bw");
        let _ = std::fs::remove_file(&path);
        assert!(
            !write_bigwig(&[], &[("chr1".to_string(), 100)], &path).unwrap(),
            "an empty track reports that nothing was written"
        );
        assert!(
            !path.exists(),
            "the format cannot express an empty track, so none is created"
        );
    }
}
