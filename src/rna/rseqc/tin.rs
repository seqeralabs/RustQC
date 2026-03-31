//! TIN (Transcript Integrity Number) analysis.
//!
//! Measures transcript integrity via Shannon entropy of read coverage
//! uniformity across sampled exonic positions. Reimplementation of
//! RSeQC's `tin.py` tool.

use rust_htslib::bam;
use std::collections::{HashMap, HashSet};
use std::hash::{BuildHasher, Hasher};
use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use indexmap::IndexMap;
use log::debug;
use rand::Rng;

use crate::gtf::Gene;

// ===================================================================
// Deterministic hash state for reproducible TIN
// ===================================================================

/// A [`BuildHasher`] that produces deterministic hashes seeded from a `u64`.
///
/// When `--tin-seed` is provided, the seed is used directly, making the
/// [`HashSet`] layout (and therefore threshold-crossing behaviour during
/// parallel accumulation) identical across runs. Without a seed, a random
/// value is drawn so the default non-deterministic behaviour is preserved.
#[derive(Clone, Debug)]
pub(crate) struct TinHashState(u64);

impl BuildHasher for TinHashState {
    type Hasher = std::collections::hash_map::DefaultHasher;

    fn build_hasher(&self) -> Self::Hasher {
        let mut h = std::collections::hash_map::DefaultHasher::new();
        h.write_u64(self.0);
        h
    }
}

impl TinHashState {
    fn new(seed: Option<u64>) -> Self {
        TinHashState(seed.unwrap_or_else(|| rand::rng().next_u64()))
    }
}

// ===================================================================
// Data structures
// ===================================================================

/// A single transcript's sampling metadata.
#[derive(Debug, Clone)]
pub struct TranscriptSampling {
    /// Gene identifier (gene_id from GTF).
    pub gene_id: String,
    /// Chromosome name (original case, for output).
    pub chrom: String,
    /// Chromosome name uppercased (for lookup matching).
    pub chrom_upper: String,
    /// Transcript start position (0-based).
    pub tx_start: u64,
    /// Transcript end position (0-based, exclusive).
    pub tx_end: u64,
    /// Sorted exon blocks as (start, end) in 0-based half-open coords.
    #[allow(dead_code)] // stored for API completeness; values used during construction
    pub exon_regions: Vec<(u64, u64)>,
    /// Total exonic bases across all exon blocks.
    #[allow(dead_code)] // stored for API completeness; values used during construction
    pub exon_length: u64,
    /// Unique genomic positions sampled within exonic regions (sorted, deduped).
    /// Used for coverage counting — each position gets one coverage slot.
    pub sampled_positions: Vec<u64>,
    /// Total position count including duplicates from upstream's `chose_bases`.
    /// Used as the denominator `l` in the TIN formula to match upstream behavior.
    /// For the small mRNA case, upstream does NOT dedup chose_bases, so
    /// `total_position_count > len(sampled_positions)` when boundary positions
    /// overlap with exonic positions.
    pub total_position_count: usize,
}

/// Index for fast lookup of sampled positions during BAM processing.
///
/// Maps chromosome → sorted list of (genomic_pos, transcript_idx, slot_idx)
/// for binary search during per-read dispatch.
#[derive(Debug, Default)]
pub struct TinIndex {
    /// All transcripts with their sampling metadata.
    pub transcripts: Vec<TranscriptSampling>,
    /// Per-chromosome sorted position entries for binary search.
    /// Each entry: (genomic_position, transcript_index, slot_index).
    pub chrom_positions: HashMap<String, Vec<(u64, u32, u32)>>,
    /// Per-chromosome transcript spans for overlap pre-filtering.
    /// Each entry: (tx_start, tx_end, transcript_index), sorted by tx_start.
    pub chrom_spans: HashMap<String, Vec<(u64, u64, u32)>>,
}

/// Result of TIN computation for a single transcript.
#[derive(Debug)]
pub struct TinResult {
    /// Gene identifier.
    pub gene_id: String,
    /// Chromosome name.
    pub chrom: String,
    /// Transcript start (0-based).
    pub tx_start: u64,
    /// Transcript end (0-based, exclusive).
    pub tx_end: u64,
    /// Computed TIN score (0-100), or 0.0 if below coverage threshold.
    pub tin: f64,
    /// Whether the transcript passed the minimum coverage threshold.
    /// Only transcripts with `passed_threshold == true` contribute to
    /// summary statistics, matching upstream RSeQC's `sample_TINs` list.
    pub passed_threshold: bool,
}

/// Collection of all TIN results for a BAM file.
#[derive(Debug)]
pub struct TinResults {
    /// Per-transcript TIN scores.
    pub transcripts: Vec<TinResult>,
}

impl TinResults {
    /// Number of transcripts with computed TIN scores.
    pub fn len(&self) -> usize {
        self.transcripts.len()
    }
}

// ===================================================================
// TIN index construction
// ===================================================================

impl TinIndex {
    /// Build TIN index from GTF gene annotations.
    ///
    /// Emits every transcript from every gene, matching RSeQC's
    /// per-transcript behaviour.
    pub fn from_genes(genes: &IndexMap<String, Gene>, sample_size: usize) -> Self {
        let mut index = TinIndex::default();

        for gene in genes.values() {
            for tx in &gene.transcripts {
                let exon_regions: Vec<(u64, u64)> = tx
                    .exons
                    .iter()
                    .map(|&(s, e)| (s.saturating_sub(1), e)) // GTF 1-based to 0-based half-open
                    .collect();
                let exon_length: u64 = exon_regions.iter().map(|(s, e)| e - s).sum();

                if exon_length == 0 {
                    continue;
                }

                let chrom = gene.chrom.clone();
                let chrom_upper = chrom.to_uppercase();
                let tx_start = exon_regions.first().map(|r| r.0).unwrap_or(0);
                let tx_end = exon_regions.last().map(|r| r.1).unwrap_or(0);

                let (sampled, total_count) =
                    sample_exonic_positions(&exon_regions, sample_size, tx_start, tx_end);
                if sampled.is_empty() {
                    continue;
                }

                index.add_transcript(TranscriptSampling {
                    gene_id: tx.transcript_id.clone(),
                    chrom,
                    chrom_upper,
                    tx_start,
                    tx_end,
                    exon_regions,
                    exon_length,
                    sampled_positions: sampled,
                    total_position_count: total_count,
                });
            }
        }

        index.build();
        index
    }

    /// Add a transcript to the index (before calling `build()`).
    fn add_transcript(&mut self, tx: TranscriptSampling) {
        let tx_idx = self.transcripts.len() as u32;
        let chrom_upper = tx.chrom_upper.clone();
        let tx_start = tx.tx_start;
        let tx_end = tx.tx_end;

        // Add sampled positions to per-chromosome position list
        for (slot_idx, &pos) in tx.sampled_positions.iter().enumerate() {
            self.chrom_positions
                .entry(chrom_upper.clone())
                .or_default()
                .push((pos, tx_idx, slot_idx as u32));
        }

        // Add transcript span
        self.chrom_spans
            .entry(chrom_upper)
            .or_default()
            .push((tx_start, tx_end, tx_idx));

        self.transcripts.push(tx);
    }

    /// Sort per-chromosome data for binary search.
    fn build(&mut self) {
        for positions in self.chrom_positions.values_mut() {
            positions.sort_unstable();
        }
        for spans in self.chrom_spans.values_mut() {
            spans.sort_unstable();
        }
        debug!(
            "TIN index: {} transcripts, {} chromosomes",
            self.transcripts.len(),
            self.chrom_positions.len()
        );
    }
}

/// Sample positions within exonic regions, matching upstream RSeQC `tin.py`'s
/// `genomic_positions()` algorithm.
///
/// Upstream uses 1-based genomic coordinates internally. All positions returned
/// here are 1-based to match pileup semantics (pysam's `pileupcolumn.pos + 1`).
///
/// When `mRNA_size <= sample_size`: returns ALL exonic bases (1-based) plus
/// `tx_start+1` and `tx_end` as boundary markers. The upstream does NOT dedup
/// these, so duplicates inflate `len(chose_bases)` which becomes the TIN
/// denominator `l`. We reproduce this by returning non-unique positions.
/// When `mRNA_size > sample_size`: returns every `step_size`-th exonic base
/// plus exon boundary positions (first and last base of each exon), deduplicated
/// and sorted.
/// Returns `(unique_positions, total_count)` where:
/// - `unique_positions`: sorted, deduplicated genomic positions for coverage counting
/// - `total_count`: length including duplicates, used as TIN denominator `l`
fn sample_exonic_positions(
    exon_regions: &[(u64, u64)],
    n: usize,
    tx_start: u64,
    tx_end: u64,
) -> (Vec<u64>, usize) {
    let mrna_size: u64 = exon_regions.iter().map(|(s, e)| e - s).sum();
    if mrna_size == 0 {
        return (Vec::new(), 0);
    }

    // Build 1-based list of all exonic positions (upstream: range(st+1, end+1))
    // exon_regions are 0-based half-open: (start, end) means bases start..end-1
    // 1-based equivalents: start+1 .. end (inclusive)
    if mrna_size <= n as u64 {
        // Small transcript: chose_bases = [tx_start+1, tx_end] + all exonic bases
        // Upstream does NOT deduplicate — duplicates inflate len(chose_bases)
        // which is used as the TIN denominator `l`.
        let mut positions: Vec<u64> = Vec::with_capacity(mrna_size as usize + 2);
        positions.push(tx_start + 1); // transcript boundary (1-based)
        positions.push(tx_end); // transcript boundary (1-based, 0-based exclusive = 1-based inclusive)
        for &(start, end) in exon_regions {
            for pos in (start + 1)..=end {
                positions.push(pos);
            }
        }
        positions.sort_unstable();
        let total_count = positions.len(); // includes duplicates

        // Dedup for coverage counting (pileup only returns unique positions)
        positions.dedup();
        (positions, total_count)
    } else {
        // Large transcript: strided sampling + exon boundaries
        // Build all 1-based exonic positions
        let mut gene_all_base: Vec<u64> = Vec::with_capacity(mrna_size as usize);
        let mut exon_bounds: Vec<u64> = Vec::new();
        for &(start, end) in exon_regions {
            for pos in (start + 1)..=end {
                gene_all_base.push(pos);
            }
            // Exon boundary positions (1-based)
            exon_bounds.push(start + 1); // first base of exon
            exon_bounds.push(end); // last base of exon
        }

        // Strided sampling: step_size = int(mRNA_size / sample_size)
        let step_size = (mrna_size as usize) / n;
        let step_size = step_size.max(1);

        let mut chose_bases: Vec<u64> = (0..gene_all_base.len())
            .step_by(step_size)
            .map(|i| gene_all_base[i])
            .collect();

        // Merge exon boundaries with sampled positions (upstream: uniqify(exon_bounds + chose_bases))
        // uniqify preserves first occurrence order and removes duplicates
        let mut all_positions = exon_bounds;
        all_positions.append(&mut chose_bases);

        let mut seen = HashSet::new();
        all_positions.retain(|x| seen.insert(*x));
        all_positions.sort_unstable();

        let total_count = all_positions.len(); // already deduped by uniqify
        (all_positions, total_count)
    }
}

// ===================================================================
// TIN accumulator
// ===================================================================

/// Per-read accumulator for TIN computation.
///
/// Tracks coverage at sampled positions and unique read start positions
/// per transcript during the BAM pass. Unique start tracking matches
/// RSeQC's `check_min_reads()`, which uses a `set()` of read start
/// positions.
///
/// To reduce memory pressure, each transcript's `unique_starts` set is
/// capped at `min_cov + 1` entries. Once the threshold is exceeded the
/// set is cleared and `exceeded_threshold` is set, since only the
/// pass/fail decision matters for the TIN score.
#[derive(Debug)]
pub struct TinAccum {
    /// Per-transcript per-slot coverage counts.
    /// Indexed as `coverage[tx_idx][slot_idx]`.
    pub coverage: Vec<Vec<u32>>,
    /// Per-transcript unique read start positions, capped at `min_cov + 1`.
    /// Once exceeded, the set is drained and `exceeded_threshold[tx_idx]`
    /// is set instead.
    pub unique_starts: Vec<HashSet<u64, TinHashState>>,
    /// Per-transcript flag: true once unique start count exceeded `min_cov`.
    /// Avoids further HashSet inserts for high-coverage transcripts.
    pub exceeded_threshold: Vec<bool>,
    /// Number of sampled slots per transcript.
    #[allow(dead_code)]
    pub n_samples: Vec<u32>,
    /// Minimum MAPQ threshold (kept for API compatibility but not used;
    /// upstream tin.py does not filter by MAPQ).
    #[allow(dead_code)]
    pub mapq_cut: u8,
    /// Minimum coverage threshold (unique read start positions).
    pub min_cov: u32,
    /// Hash state used for unique_starts sets (deterministic when seeded).
    hash_state: TinHashState,
    /// Reusable scratch buffer for aligned blocks, avoids per-read allocation.
    blocks_buf: Vec<(u64, u64)>,
    /// Last chromosome seen -- used to detect chromosome transitions and
    /// reset the cursors below.
    cursor_chrom: String,
    /// Cursor into `chrom_spans` for the current chromosome. Since the BAM
    /// is coordinate-sorted, each read's span_end is >= the previous read's,
    /// so we only need to search forward from here.
    cursor_span: usize,
    /// Cursor into `chrom_positions` for the current chromosome. The leftmost
    /// aligned-block start of each read is >= the previous read's, so we
    /// search forward from here.
    cursor_pos: usize,
}

impl TinAccum {
    /// Create a new TIN accumulator for the given index.
    ///
    /// When `seed` is `Some(n)`, the internal hash sets use a deterministic
    /// hash function seeded from `n`, making TIN results reproducible across
    /// runs. When `None`, the default random hash state is used.
    pub fn new(index: &TinIndex, mapq_cut: u8, min_cov: u32, seed: Option<u64>) -> Self {
        let n_transcripts = index.transcripts.len();
        let mut coverage = Vec::with_capacity(n_transcripts);
        let mut n_samples = Vec::with_capacity(n_transcripts);

        let hash_state = TinHashState::new(seed);

        for tx in &index.transcripts {
            let n = tx.sampled_positions.len();
            coverage.push(vec![0u32; n]);
            n_samples.push(n as u32);
        }

        let unique_starts = (0..n_transcripts)
            .map(|_| HashSet::with_hasher(hash_state.clone()))
            .collect();

        TinAccum {
            coverage,
            unique_starts,
            exceeded_threshold: vec![false; n_transcripts],
            n_samples,
            mapq_cut,
            min_cov,
            hash_state,
            blocks_buf: Vec::with_capacity(8),
            cursor_chrom: String::new(),
            cursor_span: 0,
            cursor_pos: 0,
        }
    }

    /// Process a single BAM record for TIN accumulation.
    ///
    /// Updates coverage at sampled positions that overlap the read's
    /// aligned blocks, and tracks read start positions.
    ///
    /// Matches upstream RSeQC `tin.py` filtering:
    /// - `check_min_reads()`: skips qcfail, unmapped, secondary only
    /// - `genebody_coverage()` pileup: skips is_del, qcfail, secondary, unmapped
    ///   (pysam pileup default flag_filter also excludes BAM_FDUP)
    /// - Neither function filters by MAPQ or supplementary flags
    pub fn process_read(&mut self, record: &bam::Record, chrom_upper: &str, index: &TinIndex) {
        let flags = record.flags();

        // Skip unmapped, QC-fail, secondary only (matching upstream tin.py)
        if flags & 0x4 != 0 || flags & 0x200 != 0 || flags & 0x100 != 0 {
            return;
        }

        // Fill the reusable scratch buffer with aligned blocks from CIGAR,
        // avoiding a per-read Vec allocation.
        fill_aligned_blocks(record, &mut self.blocks_buf);
        let blocks = &self.blocks_buf;
        if blocks.is_empty() {
            return;
        }

        let read_start = blocks[0].0;

        // Get chromosome positions for this chromosome
        let chrom_positions = match index.chrom_positions.get(chrom_upper) {
            Some(p) => p,
            None => return,
        };

        // Get transcript spans for this chromosome
        let chrom_spans = match index.chrom_spans.get(chrom_upper) {
            Some(s) => s,
            None => return,
        };

        // Reset cursors on chromosome transition.
        if chrom_upper != self.cursor_chrom {
            self.cursor_chrom.clear();
            self.cursor_chrom.push_str(chrom_upper);
            self.cursor_span = 0;
            self.cursor_pos = 0;
        }

        // Find overlapping transcripts via read start.
        // chrom_spans is sorted by tx_start.  We need all transcripts
        // where tx_start <= read_start < tx_end.
        //
        // Since reads are coordinate-sorted, read_start is monotonically
        // non-decreasing.  The partition_point (first span where
        // tx_start > read_start) can only move forward, so we search
        // from cursor_span instead of 0.
        let span_end = self.cursor_span
            + chrom_spans[self.cursor_span..].partition_point(|s| s.0 <= read_start);
        self.cursor_span = span_end;

        for &(_tx_start, tx_end, tx_idx) in &chrom_spans[..span_end] {
            // _tx_start <= read_start is guaranteed by the partition_point
            if read_start < tx_end {
                let idx = tx_idx as usize;
                // Skip HashSet insert if this transcript already exceeded the
                // min_cov threshold -- we only need the pass/fail decision.
                if !self.exceeded_threshold[idx] {
                    self.unique_starts[idx].insert(read_start);
                    if self.unique_starts[idx].len() > self.min_cov as usize {
                        self.exceeded_threshold[idx] = true;
                        // Free the HashSet memory; the flag is sufficient.
                        self.unique_starts[idx] = HashSet::with_hasher(self.hash_state.clone());
                    }
                }
            }
        }

        // Skip duplicate reads for coverage counting.
        // Upstream tin.py uses pysam pileup() whose default flag_filter
        // (BAM_FUNMAP|BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP = 0x704) excludes
        // duplicate-flagged reads from coverage columns. However, the upstream
        // check_min_reads() uses fetch() which does NOT filter duplicates, so
        // duplicates are allowed above for the unique_starts tracking.
        if flags & 0x400 != 0 {
            return;
        }

        // For each aligned block, find sampled positions that fall within it.
        // Sampled positions are 1-based; aligned blocks are 0-based half-open
        // [block_start, block_end). A 1-based position P is covered iff
        // P-1 >= block_start && P-1 < block_end, i.e., P > block_start && P <= block_end.
        //
        // Since reads are coordinate-sorted, block_start values are
        // monotonically non-decreasing across reads.  We advance cursor_pos
        // forward rather than binary-searching from scratch each time.
        let mut first_block = true;
        for &(block_start, block_end) in blocks.iter() {
            // Advance cursor to first position > block_start.
            // cursor_pos only moves forward across reads on the same chromosome.
            let base = self.cursor_pos;
            let pos_start = base + chrom_positions[base..].partition_point(|p| p.0 <= block_start);

            // Update cursor from the first aligned block (leftmost), since
            // the next read's first block_start >= this read's first block_start.
            if first_block {
                self.cursor_pos = pos_start;
                first_block = false;
            }

            for &(pos, tx_idx, slot_idx) in &chrom_positions[pos_start..] {
                if pos > block_end {
                    break;
                }
                // This sampled position is covered by this aligned block
                self.coverage[tx_idx as usize][slot_idx as usize] += 1;
            }
        }
    }

    /// Merge another accumulator into this one.
    pub fn merge(&mut self, other: TinAccum) {
        for (i, other_cov) in other.coverage.into_iter().enumerate() {
            for (j, count) in other_cov.into_iter().enumerate() {
                self.coverage[i][j] += count;
            }
            // If either side already exceeded the threshold, the merged
            // result also exceeds it -- no need to union the sets.
            if self.exceeded_threshold[i] || other.exceeded_threshold[i] {
                self.exceeded_threshold[i] = true;
                self.unique_starts[i] = HashSet::with_hasher(self.hash_state.clone());
            } else {
                self.unique_starts[i].extend(other.unique_starts[i].iter());
                if self.unique_starts[i].len() > self.min_cov as usize {
                    self.exceeded_threshold[i] = true;
                    self.unique_starts[i] = HashSet::with_hasher(self.hash_state.clone());
                }
            }
        }
    }

    /// Compute TIN scores for all transcripts.
    pub fn into_result(self, index: &TinIndex) -> TinResults {
        let mut transcripts = Vec::with_capacity(index.transcripts.len());

        for (i, tx) in index.transcripts.iter().enumerate() {
            // A transcript passes if the exceeded_threshold flag was set
            // (unique starts > min_cov during accumulation or merge), or
            // if the remaining set still exceeds the threshold.
            let passed =
                self.exceeded_threshold[i] || self.unique_starts[i].len() as u32 > self.min_cov;
            let coverage = &self.coverage[i];

            // Upstream: `if len(read_count) > cutoff` -- strictly greater
            if !passed {
                transcripts.push(TinResult {
                    gene_id: tx.gene_id.clone(),
                    chrom: tx.chrom.clone(),
                    tx_start: tx.tx_start,
                    tx_end: tx.tx_end,
                    tin: 0.0,
                    passed_threshold: false,
                });
                continue;
            }

            // Use total_position_count (which includes duplicates from upstream's
            // chose_bases list) as the denominator, matching upstream's l = len(pick_positions)
            let n_total = tx.total_position_count;
            let tin = compute_tin(coverage, n_total);

            transcripts.push(TinResult {
                gene_id: tx.gene_id.clone(),
                chrom: tx.chrom.clone(),
                tx_start: tx.tx_start,
                tx_end: tx.tx_end,
                tin,
                passed_threshold: true,
            });
        }

        TinResults { transcripts }
    }
}

/// Compute TIN score from coverage vector using Shannon entropy.
///
/// Matches upstream RSeQC's `tin_score()`:
///   cvg_eff = [c for c in cvg if c > 0]
///   entropy = shannon_entropy(cvg_eff)
///   tin = 100 * exp(entropy) / l
///
/// where `l` is the total number of sampled positions (including zeros),
/// NOT just the nonzero count. A transcript with many uncovered positions
/// gets a lower TIN because the denominator is larger.
fn compute_tin(coverage: &[u32], n_total_positions: usize) -> f64 {
    if n_total_positions == 0 {
        return 0.0;
    }

    let nonzero: Vec<f64> = coverage
        .iter()
        .filter(|&&c| c > 0)
        .map(|&c| c as f64)
        .collect();

    if nonzero.is_empty() {
        return 0.0;
    }

    let total: f64 = nonzero.iter().sum();
    if total == 0.0 {
        return 0.0;
    }

    // Shannon entropy over non-zero positions only (matching upstream)
    let entropy: f64 = nonzero
        .iter()
        .map(|&c| {
            let p = c / total;
            -p * p.ln()
        })
        .sum();

    // Upstream: tin = 100 * exp(entropy) / l
    // where l = total number of sampled positions (not just nonzero)
    100.0 * entropy.exp() / n_total_positions as f64
}

/// Get aligned blocks from a BAM record's CIGAR string.
///
/// Returns sorted (start, end) pairs of aligned blocks (0-based half-open).
/// Fill `buf` with aligned blocks from the CIGAR, reusing the existing
/// Vec capacity to avoid per-read heap allocation.
fn fill_aligned_blocks(record: &bam::Record, buf: &mut Vec<(u64, u64)>) {
    use rust_htslib::bam::record::Cigar;

    buf.clear();
    let mut pos = record.pos() as u64;

    for op in record.cigar().iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                buf.push((pos, pos + *len as u64));
                pos += *len as u64;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                pos += *len as u64;
            }
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
}

// ===================================================================
// Output writers
// ===================================================================

/// Write per-transcript TIN scores in RSeQC format.
///
/// Format: TSV with header geneID\tchrom\ttx_start\ttx_end\tTIN
pub fn write_tin(results: &TinResults, output_path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(output_path)
        .with_context(|| format!("Failed to create TIN output: {}", output_path.display()))?;

    writeln!(f, "geneID\tchrom\ttx_start\ttx_end\tTIN")?;

    for r in &results.transcripts {
        // Match upstream formatting:
        // - Failed check_min_reads: explicitly writes "0.0" (Python float literal)
        // - Passed but tin_score returned int 0 (empty coverage): writes "0"
        // - Non-zero TIN: writes float value
        if !r.passed_threshold {
            // Below coverage threshold → upstream writes explicit 0.0
            writeln!(
                f,
                "{}\t{}\t{}\t{}\t0.0",
                r.gene_id, r.chrom, r.tx_start, r.tx_end
            )?;
        } else if r.tin == 0.0 {
            // Passed threshold but zero TIN → upstream tin_score returns int 0
            writeln!(
                f,
                "{}\t{}\t{}\t{}\t0",
                r.gene_id, r.chrom, r.tx_start, r.tx_end
            )?;
        } else {
            writeln!(
                f,
                "{}\t{}\t{}\t{}\t{}",
                r.gene_id, r.chrom, r.tx_start, r.tx_end, r.tin
            )?;
        }
    }

    Ok(())
}

/// Write TIN summary statistics.
///
/// Format: Bam_file\tTIN(mean)\tTIN(median)\tTIN(stdev)
pub fn write_tin_summary(results: &TinResults, bam_name: &str, output_path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(output_path)
        .with_context(|| format!("Failed to create TIN summary: {}", output_path.display()))?;

    let scores: Vec<f64> = results
        .transcripts
        .iter()
        .filter(|r| r.passed_threshold)
        .map(|r| r.tin)
        .collect();

    let (mean, median, stdev) = if scores.is_empty() {
        (0.0, 0.0, 0.0)
    } else {
        let n = scores.len() as f64;
        let mean = scores.iter().sum::<f64>() / n;

        let median = crate::io::median(&scores);

        let variance = scores.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;
        let stdev = variance.sqrt();

        (mean, median, stdev)
    };

    writeln!(f, "Bam_file\tTIN(mean)\tTIN(median)\tTIN(stdev)")?;
    writeln!(f, "{}\t{}\t{}\t{}", bam_name, mean, median, stdev)?;

    Ok(())
}

// ===================================================================
// Unit tests
// ===================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sample_exonic_positions() {
        // Single exon of 10 bases: 0-based [100, 110)
        // tx_start=100, tx_end=110 (matching BED format)
        // mRNA_size = 10, sample_size = 5, so mRNA_size > sample_size
        // step_size = 10 / 5 = 2
        // gene_all_base (1-based) = [101, 102, ..., 110]
        // exon_bounds = [101, 110]
        // strided: indices 0,2,4,6,8 -> [101, 103, 105, 107, 109]
        // merged with bounds [101, 110] -> [101, 103, 105, 107, 109, 110]
        // Large mRNA case: deduped by uniqify, total_count = unique count
        let exons = vec![(100, 110)];
        let (positions, total_count) = sample_exonic_positions(&exons, 5, 100, 110);
        assert_eq!(positions.len(), 6); // 5 strided + exon end boundary 110
        assert_eq!(total_count, 6); // no duplicates in large mRNA case
                                    // Positions should be 1-based
        assert_eq!(positions[0], 101);
        assert_eq!(positions[5], 110);
    }

    #[test]
    fn test_sample_multi_exon() {
        // Two exons: 0-based [100,105) and [200,205) → 10 bases total
        // tx_start=100, tx_end=205
        // mRNA_size = 10, sample_size = 10, so mRNA_size <= sample_size
        // chose_bases = [tx_start+1, tx_end] + all exonic = [101, 205, 101..=105, 201..=205]
        // Sorted with dupes: [101, 101, 102, 103, 104, 105, 201, 202, 203, 204, 205, 205]
        // total_count = 12 (with duplicates), unique positions = 10
        let exons = vec![(100, 105), (200, 205)];
        let (positions, total_count) = sample_exonic_positions(&exons, 10, 100, 205);
        assert_eq!(positions.len(), 10); // 10 unique positions
        assert_eq!(total_count, 12); // 10 + 2 boundary duplicates
        assert_eq!(positions[0], 101);
        assert_eq!(positions[9], 205);
    }

    #[test]
    fn test_sample_tiny_exon() {
        // Exon smaller than sample size: 0-based [100, 103) → 3 bases
        // tx_start=100, tx_end=103
        // mRNA_size = 3 <= sample_size = 100
        // chose_bases = [101, 103, 101, 102, 103] -> sorted: [101, 101, 102, 103, 103]
        // total_count = 5 (with duplicates), unique positions = 3
        let exons = vec![(100, 103)];
        let (positions, total_count) = sample_exonic_positions(&exons, 100, 100, 103);
        assert_eq!(positions.len(), 3); // 3 unique exonic positions
        assert_eq!(total_count, 5); // 3 + 2 boundary duplicates
        assert_eq!(positions[0], 101);
        assert_eq!(positions[2], 103);
    }

    #[test]
    fn test_compute_tin_uniform() {
        // Uniform coverage → max TIN
        let cov = vec![10, 10, 10, 10, 10, 10, 10, 10, 10, 10];
        let tin = compute_tin(&cov, cov.len());
        assert!(
            (tin - 100.0).abs() < 0.01,
            "Uniform coverage TIN should be ~100: {tin}"
        );
    }

    #[test]
    fn test_compute_tin_uniform_with_zeros() {
        // Uniform coverage but only 5 of 10 positions covered
        // TIN = 100 * exp(H) / 10; H entropy over 5 equal values -> exp(H)=5
        // TIN = 100 * 5 / 10 = 50
        let cov = vec![10, 10, 10, 10, 10, 0, 0, 0, 0, 0];
        let tin = compute_tin(&cov, cov.len());
        assert!(
            (tin - 50.0).abs() < 0.01,
            "Half-covered uniform TIN should be ~50: {tin}"
        );
    }

    #[test]
    fn test_compute_tin_degraded() {
        // 5' degradation: high at start, low at end
        let cov = vec![100, 80, 60, 40, 20, 10, 5, 2, 1, 1];
        let tin = compute_tin(&cov, cov.len());
        assert!(
            tin > 0.0 && tin < 100.0,
            "Degraded TIN should be between 0 and 100: {tin}"
        );
    }

    #[test]
    fn test_compute_tin_all_zero() {
        let cov = vec![0, 0, 0, 0, 0];
        let tin = compute_tin(&cov, cov.len());
        assert_eq!(tin, 0.0);
    }

    #[test]
    fn test_compute_tin_single_nonzero() {
        // Only one position covered out of 5 -> entropy=0 -> exp(0)=1
        // TIN = 100 * 1 / 5 = 20
        let cov = vec![0, 0, 100, 0, 0];
        let tin = compute_tin(&cov, cov.len());
        assert!(
            (tin - 20.0).abs() < 0.01,
            "Single position coverage should give TIN=20: {tin}"
        );
    }
}
