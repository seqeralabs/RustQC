//! Genome coverage accumulation matching `bedtools genomecov`.
//!
//! Replicates nf-core/rnaseq bigWig generation parameters:
//! - Combined track: `-split -bg`
//! - Per-strand tracks (stranded libraries): `-split -du -strand +/- -bg`
//!
//! Coverage uses a difference-array (starts/ends) per chromosome, identical to
//! bedtools' `BedGenomeCoverage` implementation.

use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

use crate::rna::bam_flags::BAM_FREAD2;
use crate::Strandedness;

/// Strand filter for per-strand coverage tracks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrandFilter {
    /// Forward reference strand (`-strand +`).
    Plus,
    /// Reverse reference strand (`-strand -`).
    Minus,
}

/// Settings for a single coverage track.
///
/// All nf-core/rnaseq tracks use `-split` (reads are always split at splicing
/// blocks), so splitting is unconditional rather than a toggle.
#[derive(Debug, Clone, Copy)]
struct TrackSettings {
    /// Flip second-mate strand for strand-specific libraries (`-du`).
    deduplicate_strand: bool,
    /// Optional strand filter (`-strand`).
    strand: Option<StrandFilter>,
}

/// Difference-array endpoint counters for one reference position.
#[derive(Debug, Default, Clone, Copy)]
struct DiffPoint {
    starts: u32,
    ends: u32,
}

/// Per-chromosome coverage using bedtools' starts/ends difference array.
#[derive(Debug, Clone)]
pub(crate) struct ChromCoverage {
    /// Chromosome length in bases.
    size: u64,
    /// Per-base start/end counters (length == size).
    points: Vec<DiffPoint>,
}

impl ChromCoverage {
    fn new(size: u64) -> Self {
        Self {
            size,
            points: vec![DiffPoint::default(); size as usize],
        }
    }

    /// Add coverage for `[start, end]` inclusive, matching bedtools `AddCoverage`.
    fn add_interval(&mut self, start: u64, end: u64) {
        if self.size == 0 {
            return;
        }
        let chrom_size = self.size;
        if start < chrom_size {
            self.points[start as usize].starts += 1;
        }
        if end < chrom_size {
            self.points[end as usize].ends += 1;
        } else {
            self.points[(chrom_size - 1) as usize].ends += 1;
        }
    }

    fn merge(&mut self, other: &Self) {
        debug_assert_eq!(self.size, other.size);
        for (a, b) in self.points.iter_mut().zip(other.points.iter()) {
            a.starts += b.starts;
            a.ends += b.ends;
        }
    }
}

/// Accumulator for one coverage track (combined, forward, or reverse).
#[derive(Debug, Default, Clone)]
pub struct TrackAccum {
    chroms: Vec<(String, ChromCoverage)>,
    /// Maps chromosome name to its index in `chroms` for O(1) lookup in the
    /// per-read hot loop (avoids a linear name scan on every record).
    index: std::collections::HashMap<String, usize>,
}

impl TrackAccum {
    fn new(chroms: &[(String, u64)]) -> Self {
        let chroms: Vec<(String, ChromCoverage)> = chroms
            .iter()
            .map(|(name, len)| (name.clone(), ChromCoverage::new(*len)))
            .collect();
        let index = chroms
            .iter()
            .enumerate()
            .map(|(i, (name, _))| (name.clone(), i))
            .collect();
        Self { chroms, index }
    }

    fn chrom_mut(&mut self, chrom: &str) -> Option<&mut ChromCoverage> {
        let i = *self.index.get(chrom)?;
        Some(&mut self.chroms[i].1)
    }

    fn chrom(&self, chrom: &str) -> Option<&ChromCoverage> {
        let i = *self.index.get(chrom)?;
        Some(&self.chroms[i].1)
    }

    fn merge(&mut self, other: Self) {
        for (name, other_cov) in other.chroms {
            if let Some(&i) = self.index.get(&name) {
                self.chroms[i].1.merge(&other_cov);
            } else {
                self.index.insert(name.clone(), self.chroms.len());
                self.chroms.push((name, other_cov));
            }
        }
    }
}

/// Bundle of coverage accumulators for all nf-core/rnaseq bigWig tracks.
#[derive(Debug, Clone)]
pub struct GenomeCovAccum {
    /// Strand-agnostic combined coverage (`-split -bg`).
    pub combined: TrackAccum,
    /// Forward-genome-strand coverage (stranded libraries only).
    pub forward: Option<TrackAccum>,
    /// Reverse-genome-strand coverage (stranded libraries only).
    pub reverse: Option<TrackAccum>,
    /// Scratch buffer for split-read alignment blocks, reused across reads to
    /// avoid a per-record heap allocation in the hot loop.
    blocks_buf: Vec<(u64, u64)>,
}

impl GenomeCovAccum {
    /// Create accumulators for the given chromosomes.
    ///
    /// Per-strand tracks are only created for stranded libraries.
    pub fn new(stranded: Strandedness, chroms: &[(String, u64)]) -> Self {
        let (forward, reverse) = if stranded == Strandedness::Unstranded {
            (None, None)
        } else {
            (Some(TrackAccum::new(chroms)), Some(TrackAccum::new(chroms)))
        };
        Self {
            combined: TrackAccum::new(chroms),
            forward,
            reverse,
            blocks_buf: Vec::new(),
        }
    }

    /// Merge another worker's accumulators into this one.
    pub fn merge(&mut self, other: Self) {
        self.combined.merge(other.combined);
        if let (Some(ref mut f), Some(o)) = (&mut self.forward, other.forward) {
            f.merge(o);
        }
        if let (Some(ref mut r), Some(o)) = (&mut self.reverse, other.reverse) {
            r.merge(o);
        }
    }

    /// Process a mapped BAM record, updating all enabled tracks.
    pub fn process_read(&mut self, record: &bam::Record, chrom: &str, stranded: Strandedness) {
        if record.is_unmapped() {
            return;
        }

        process_track(
            &mut self.combined,
            record,
            chrom,
            TrackSettings {
                deduplicate_strand: false,
                strand: None,
            },
            &mut self.blocks_buf,
        );

        if stranded == Strandedness::Unstranded {
            return;
        }

        let (fwd_strand, rev_strand) = match stranded {
            Strandedness::Forward => (StrandFilter::Plus, StrandFilter::Minus),
            Strandedness::Reverse => (StrandFilter::Minus, StrandFilter::Plus),
            Strandedness::Unstranded => unreachable!(),
        };

        if let Some(ref mut forward) = self.forward {
            process_track(
                forward,
                record,
                chrom,
                TrackSettings {
                    deduplicate_strand: true,
                    strand: Some(fwd_strand),
                },
                &mut self.blocks_buf,
            );
        }
        if let Some(ref mut reverse) = self.reverse {
            process_track(
                reverse,
                record,
                chrom,
                TrackSettings {
                    deduplicate_strand: true,
                    strand: Some(rev_strand),
                },
                &mut self.blocks_buf,
            );
        }
    }
}

fn process_track(
    track: &mut TrackAccum,
    record: &bam::Record,
    chrom: &str,
    settings: TrackSettings,
    blocks_buf: &mut Vec<(u64, u64)>,
) {
    let Some(cov) = track.chrom_mut(chrom) else {
        return;
    };

    let mut is_reverse = record.is_reverse();
    // bedtools `-du`: flip second mate strand so both mates appear on the
    // same fragment strand (dUTP/strand-specific PE libraries).
    if settings.deduplicate_strand
        && record.is_paired()
        && !record.is_mate_unmapped()
        && record.flags() & BAM_FREAD2 != 0
    {
        is_reverse = !is_reverse;
    }

    if let Some(filter) = settings.strand {
        let want_reverse = filter == StrandFilter::Minus;
        if is_reverse != want_reverse {
            return;
        }
    }

    // All nf-core/rnaseq tracks use `-split`: accumulate per alignment block.
    genomecov_blocks(record, blocks_buf);
    for &(start, end) in blocks_buf.iter() {
        // bedtools uses inclusive end; blocks are half-open [start, end).
        cov.add_interval(start, end.saturating_sub(1));
    }
}

/// Extract alignment blocks for `-split` mode, matching bedtools `GetBamBlocks`
/// with `breakOnDeletionOps=true` (default, no `-ignoreD`).
fn genomecov_blocks(record: &bam::Record, blocks: &mut Vec<(u64, u64)>) {
    blocks.clear();
    let mut ref_pos = record.pos() as u64;

    for op in record.cigar().iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len = *len as u64;
                blocks.push((ref_pos, ref_pos + len));
                ref_pos += len;
            }
            Cigar::RefSkip(len) | Cigar::Del(len) => {
                ref_pos += *len as u64;
            }
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
}

/// Merged genome coverage ready for bedGraph / bigWig output.
#[derive(Debug)]
pub struct GenomeCovResult {
    /// Combined strand-agnostic track.
    pub combined: TrackAccum,
    /// Forward-genome-strand track (None for unstranded libraries).
    pub forward: Option<TrackAccum>,
    /// Reverse-genome-strand track (None for unstranded libraries).
    pub reverse: Option<TrackAccum>,
    /// Chromosome names and lengths in BAM header order, carried through so
    /// output does not need to re-open the BAM to recover them.
    pub chrom_sizes: Vec<(String, u64)>,
}

impl GenomeCovAccum {
    /// Finalize the accumulator into a result, attaching the BAM header-ordered
    /// chromosome sizes used for bedGraph flattening and bigWig writing.
    pub fn into_result(self, chrom_sizes: Vec<(String, u64)>) -> GenomeCovResult {
        GenomeCovResult {
            combined: self.combined,
            forward: self.forward,
            reverse: self.reverse,
            chrom_sizes,
        }
    }
}

/// Convert a chromosome's difference array to bedGraph intervals (depth > 0).
///
/// Output matches bedtools `ReportChromCoverageBedGraph` with scale=1.
pub(crate) fn chrom_to_bedgraph(chrom: &str, cov: &ChromCoverage) -> Vec<(String, u32, u32, u32)> {
    let mut out = Vec::new();
    let mut depth: i64 = 0;
    let mut last_start: i64 = -1;
    let mut last_depth: i64 = -1;

    for (pos, point) in cov.points.iter().enumerate() {
        depth += point.starts as i64;

        if depth != last_depth {
            if last_depth > 0 {
                out.push((
                    chrom.to_string(),
                    last_start as u32,
                    pos as u32,
                    last_depth as u32,
                ));
            }
            last_depth = depth;
            last_start = pos as i64;
        }

        depth -= point.ends as i64;
    }

    if last_depth > 0 {
        out.push((
            chrom.to_string(),
            last_start as u32,
            cov.size as u32,
            last_depth as u32,
        ));
    }

    out
}

/// Clip bedGraph intervals to chromosome sizes (UCSC `bedClip` behaviour).
///
/// `chrom_to_bedgraph` already bounds intervals to each chromosome's length, so
/// in the normal pipeline this is defensive: it guarantees `bedClip` parity and
/// drops any entry whose chromosome is absent from `chrom_sizes`.
pub fn clip_bedgraph(
    entries: Vec<(String, u32, u32, u32)>,
    chrom_sizes: &[(String, u64)],
) -> Vec<(String, u32, u32, u32)> {
    let sizes: std::collections::HashMap<&str, u64> =
        chrom_sizes.iter().map(|(n, s)| (n.as_str(), *s)).collect();

    entries
        .into_iter()
        .filter_map(|(chrom, start, end, val)| {
            let size = *sizes.get(chrom.as_str())? as u32;
            if start >= size {
                return None;
            }
            let clipped_end = end.min(size);
            if start >= clipped_end {
                return None;
            }
            Some((chrom, start, clipped_end, val))
        })
        .collect()
}

/// Flatten a track accumulator into bedGraph entries in chromosome order.
pub fn track_to_bedgraph(
    track: &TrackAccum,
    chrom_order: &[(String, u64)],
) -> Vec<(String, u32, u32, u32)> {
    let mut entries = Vec::new();
    for (chrom, _) in chrom_order {
        if let Some(cov) = track.chrom(chrom) {
            entries.extend(chrom_to_bedgraph(chrom, cov));
        }
    }
    entries
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chrom_to_bedgraph_simple() {
        let mut cov = ChromCoverage::new(10);
        cov.add_interval(2, 4); // positions 2,3,4 covered
        cov.add_interval(2, 4); // depth 2

        let bg = chrom_to_bedgraph("chr1", &cov);
        assert_eq!(bg, vec![("chr1".to_string(), 2, 5, 2)]);
    }

    #[test]
    fn test_clip_bedgraph() {
        let entries = vec![("chr1".to_string(), 100, 200, 5)];
        let sizes = vec![("chr1".to_string(), 150u64)];
        let clipped = clip_bedgraph(entries, &sizes);
        assert_eq!(clipped, vec![("chr1".to_string(), 100, 150, 5)]);
    }
}
