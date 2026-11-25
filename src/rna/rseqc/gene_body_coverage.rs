//! Gene body coverage profiling (5' → 3').
//!
//! Reimplementation of RSeQC's `geneBody_coverage.py`. Each transcript's
//! exonic bases are reduced to 100 percentile positions; read coverage at
//! those positions is accumulated across transcripts into a 100-bin profile
//! running 5' → 3'. Non-uniform profiles are the standard signal for RNA
//! degradation and 3'-end bias.
//!
//! Note that this is a different measurement from the Qualimap gene body
//! coverage RustQC also produces: Qualimap bins each transcript into fixed
//! relative windows, while RSeQC samples percentile *positions* and pools the
//! raw depth across transcripts.

use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use indexmap::IndexMap;
use log::debug;
use rust_htslib::bam;

use crate::gtf::Gene;
use crate::rna::bam_flags::{BAM_FDUP, BAM_FQCFAIL, BAM_FSECONDARY, BAM_FUNMAP};

/// Number of percentile bins along the gene body.
pub const NUM_BINS: usize = 100;

/// Minimum mRNA length for a transcript to be profiled.
///
/// Matches RSeQC's `--minimum_length` default; shorter transcripts are skipped.
pub const MIN_MRNA_LENGTH: usize = 100;

// ===================================================================
// Index construction
// ===================================================================

/// Sampling metadata for one profiled transcript.
#[derive(Debug, Clone)]
struct GbcTranscript {
    /// Whether the transcript is on the reverse strand (profile is flipped).
    is_reverse: bool,
    /// Offset of this transcript's first coverage slot.
    slot_offset: usize,
    /// Number of coverage slots (unique percentile positions, ≤ [`NUM_BINS`]).
    num_slots: usize,
}

/// Percentile positions for all profiled transcripts, indexed for per-read lookup.
#[derive(Debug, Default)]
pub struct GeneBodyIndex {
    /// Per-transcript metadata, in construction order.
    transcripts: Vec<GbcTranscript>,
    /// Uppercased chromosome → sorted `(position_1based, slot)` entries.
    chrom_positions: HashMap<String, Vec<(u64, u32)>>,
    /// Total number of coverage slots across all transcripts.
    total_slots: usize,
}

impl GeneBodyIndex {
    /// Number of transcripts contributing to the profile.
    pub fn num_transcripts(&self) -> usize {
        self.transcripts.len()
    }

    /// Whether any transcript was long enough to profile.
    pub fn is_empty(&self) -> bool {
        self.transcripts.is_empty()
    }

    /// Build the index from GTF gene annotations.
    ///
    /// One representative transcript per gene is profiled — the one with the
    /// most exonic bases — matching how RustQC's TIN analysis picks a
    /// representative. RSeQC instead profiles every entry of a BED12 reference
    /// (typically a housekeeping-gene set).
    pub fn from_genes(genes: &IndexMap<String, Gene>) -> Self {
        let mut transcripts = Vec::new();
        let mut chrom_positions: HashMap<String, Vec<(u64, u32)>> = HashMap::new();
        let mut total_slots = 0usize;

        for gene in genes.values() {
            // Representative transcript: most exonic bases. Genes without
            // transcript records fall back to their merged exon list.
            let exons: Vec<(u64, u64)> = match gene
                .transcripts
                .iter()
                .max_by_key(|t| t.exons.iter().map(|(s, e)| e - s + 1).sum::<u64>())
            {
                Some(tx) => tx.exons.clone(),
                None => gene.exons.iter().map(|e| (e.start, e.end)).collect(),
            };
            if exons.is_empty() {
                continue;
            }

            // 1-based inclusive exonic positions, ascending. GTF exon
            // coordinates are already 1-based inclusive.
            let mut sorted_exons = exons;
            sorted_exons.sort_unstable();
            let mrna_length: u64 = sorted_exons.iter().map(|(s, e)| e - s + 1).sum();
            if (mrna_length as usize) < MIN_MRNA_LENGTH {
                continue;
            }

            let mut bases: Vec<u64> = Vec::with_capacity(mrna_length as usize);
            for &(start, end) in &sorted_exons {
                bases.extend(start..=end);
            }

            let mut positions = percentile_positions(&bases);
            positions.sort_unstable();
            positions.dedup();
            if positions.is_empty() {
                continue;
            }

            let slot_offset = total_slots;
            let num_slots = positions.len();
            total_slots += num_slots;

            let chrom_upper = gene.chrom.to_uppercase();
            let entries = chrom_positions.entry(chrom_upper).or_default();
            for (i, pos) in positions.into_iter().enumerate() {
                entries.push((pos, (slot_offset + i) as u32));
            }

            transcripts.push(GbcTranscript {
                is_reverse: gene.strand == '-',
                slot_offset,
                num_slots,
            });
        }

        for entries in chrom_positions.values_mut() {
            entries.sort_unstable();
        }

        debug!(
            "Built gene body coverage index: {} transcripts, {} sampled positions",
            transcripts.len(),
            total_slots
        );

        GeneBodyIndex {
            transcripts,
            chrom_positions,
            total_slots,
        }
    }
}

/// Compute the 100 percentile positions of a sorted position list.
///
/// Direct port of RSeQC's `mystat.percentile_list`: for `i` in 1..=100 it takes
/// `k = (n - 1) * i / 100` and either the exact element (when `k` is integral)
/// or the rounded linear interpolation between the neighbouring elements.
/// Lists shorter than 100 entries are returned unchanged.
fn percentile_positions(sorted: &[u64]) -> Vec<u64> {
    if sorted.is_empty() {
        return Vec::new();
    }
    if sorted.len() < NUM_BINS {
        return sorted.to_vec();
    }

    let n = sorted.len();
    let mut out = Vec::with_capacity(NUM_BINS);
    for i in 1..=NUM_BINS {
        let k = (n - 1) as f64 * i as f64 / 100.0;
        let f = k.floor();
        let c = k.ceil();
        if (f - c).abs() < f64::EPSILON {
            out.push(sorted[k as usize]);
        } else {
            let d0 = sorted[f as usize] as f64 * (c - k);
            let d1 = sorted[c as usize] as f64 * (k - f);
            out.push(round_half_to_even(d0 + d1) as u64);
        }
    }
    out
}

/// Round half to even, matching Python's built-in `round()`.
///
/// `f64::round` rounds halves away from zero, which shifts interpolated
/// percentile positions by one base relative to upstream whenever the
/// interpolation lands exactly on `.5` — enough to change a bin's coverage by
/// a read. Only non-negative values occur here.
fn round_half_to_even(x: f64) -> f64 {
    let floor = x.floor();
    if x - floor == 0.5 {
        if (floor as i64) % 2 == 0 {
            floor
        } else {
            floor + 1.0
        }
    } else {
        x.round()
    }
}

// ===================================================================
// Accumulator
// ===================================================================

/// Per-worker coverage counts at the sampled positions.
#[derive(Debug, Clone)]
pub struct GeneBodyCoverageAccum {
    /// One counter per sampled position, indexed by slot.
    counts: Vec<u32>,
}

impl GeneBodyCoverageAccum {
    /// Create an accumulator sized for the given index.
    pub fn new(index: &GeneBodyIndex) -> Self {
        Self {
            counts: vec![0; index.total_slots],
        }
    }

    /// Add a single alignment record's coverage.
    ///
    /// Filters match the ones RSeQC applies inside its pileup loop: unmapped,
    /// QC-failed, secondary and duplicate records are ignored, and positions
    /// falling in a deletion or a skipped (`N`) region are not counted.
    ///
    /// # Arguments
    /// * `record` - The alignment record
    /// * `chrom_upper` - Uppercased chromosome name for this record
    /// * `index` - The sampled-position index
    pub fn process_read(&mut self, record: &bam::Record, chrom_upper: &str, index: &GeneBodyIndex) {
        let flags = record.flags();
        if flags & (BAM_FUNMAP | BAM_FQCFAIL | BAM_FSECONDARY | BAM_FDUP) != 0 {
            return;
        }

        let Some(entries) = index.chrom_positions.get(chrom_upper) else {
            return;
        };

        // Walk the CIGAR, counting only reference-consuming aligned blocks
        // (M/=/X). Deletions and reference skips leave the read uncovered.
        let mut ref_pos = record.pos() as u64; // 0-based
        for op in record.cigar().iter() {
            use rust_htslib::bam::record::Cigar;
            match *op {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    // Convert to 1-based inclusive [start, end]
                    let start = ref_pos + 1;
                    let end = ref_pos + len as u64;
                    self.add_block(entries, start, end);
                    ref_pos += len as u64;
                }
                Cigar::Del(len) | Cigar::RefSkip(len) => ref_pos += len as u64,
                _ => {}
            }
        }
    }

    /// Increment every sampled position inside the 1-based inclusive block.
    fn add_block(&mut self, entries: &[(u64, u32)], start: u64, end: u64) {
        let lo = entries.partition_point(|&(pos, _)| pos < start);
        for &(pos, slot) in &entries[lo..] {
            if pos > end {
                break;
            }
            self.counts[slot as usize] += 1;
        }
    }

    /// Merge another worker's counts into this accumulator.
    pub fn merge(&mut self, other: GeneBodyCoverageAccum) {
        for (a, b) in self.counts.iter_mut().zip(other.counts.iter()) {
            *a += b;
        }
    }

    /// Aggregate per-position counts into the 5' → 3' profile.
    ///
    /// Reverse-strand transcripts have their coverage vector flipped before
    /// aggregation, so bin 1 is always the 5' end.
    pub fn into_result(self, index: &GeneBodyIndex) -> GeneBodyCoverageResult {
        let mut coverage = vec![0u64; NUM_BINS];
        let mut touched = vec![false; NUM_BINS];

        for tx in &index.transcripts {
            let slots = &self.counts[tx.slot_offset..tx.slot_offset + tx.num_slots];
            for (i, &count) in slots.iter().enumerate() {
                let bin = if tx.is_reverse {
                    tx.num_slots - 1 - i
                } else {
                    i
                };
                coverage[bin] += count as u64;
                touched[bin] = true;
            }
        }

        GeneBodyCoverageResult {
            coverage,
            touched,
            num_transcripts: index.num_transcripts(),
        }
    }
}

/// Aggregated 100-bin gene body coverage profile.
#[derive(Debug, Clone)]
pub struct GeneBodyCoverageResult {
    /// Total coverage per percentile bin, 5' → 3'.
    pub coverage: Vec<u64>,
    /// Whether any transcript contributed to each bin.
    touched: Vec<bool>,
    /// Number of transcripts profiled.
    pub num_transcripts: usize,
}

impl GeneBodyCoverageResult {
    /// Coverage normalised to the range [0, 1] by the maximum bin.
    ///
    /// This is what RSeQC plots, and what makes profiles comparable between
    /// samples of different depth.
    pub fn normalised(&self) -> Vec<f64> {
        let max = self.coverage.iter().copied().max().unwrap_or(0);
        if max == 0 {
            return vec![0.0; self.coverage.len()];
        }
        self.coverage
            .iter()
            .map(|&v| v as f64 / max as f64)
            .collect()
    }
}

// ===================================================================
// Output
// ===================================================================

/// Write the `.geneBodyCoverage.txt` table.
///
/// Format matches RSeQC: a `Percentile` header row listing bins 1–100,
/// followed by one row per sample.
pub fn write_gene_body_coverage(
    result: &GeneBodyCoverageResult,
    sample_name: &str,
    output_path: &Path,
) -> Result<()> {
    let mut out = std::fs::File::create(output_path).with_context(|| {
        format!(
            "Failed to create gene body coverage file: {}",
            output_path.display()
        )
    })?;

    let header: Vec<String> = (1..=NUM_BINS).map(|i| i.to_string()).collect();
    writeln!(out, "Percentile\t{}", header.join("\t"))?;

    let values: Vec<String> = result
        .coverage
        .iter()
        .zip(result.touched.iter())
        .map(|(&v, &touched)| {
            if touched {
                // Upstream accumulates Python floats here, so values print
                // with a trailing ".0"
                format!("{:.1}", v as f64)
            } else {
                "0".to_string()
            }
        })
        .collect();
    writeln!(out, "{}\t{}", sample_name, values.join("\t"))?;

    debug!(
        "Wrote gene body coverage to {} ({} transcripts)",
        output_path.display(),
        result.num_transcripts
    );
    Ok(())
}

/// Write the R plotting script RSeQC produces alongside the table.
///
/// RustQC renders the plot itself; the script is written for drop-in
/// compatibility with pipelines that expect it.
pub fn write_gene_body_r_script(
    result: &GeneBodyCoverageResult,
    sample_name: &str,
    output_prefix: &str,
    output_path: &Path,
) -> Result<()> {
    let mut out = std::fs::File::create(output_path).with_context(|| {
        format!(
            "Failed to create gene body coverage R script: {}",
            output_path.display()
        )
    })?;

    let normalised: Vec<String> = result
        .normalised()
        .iter()
        .map(|v| format!("{v:.10}"))
        .collect();

    writeln!(out, "{}<- c({})", sample_name, normalised.join(","))?;
    writeln!(out, "x <- 1:{NUM_BINS}")?;
    writeln!(out, "pdf(\"{output_prefix}.geneBodyCoverage.curves.pdf\")")?;
    writeln!(
        out,
        "plot(x,{sample_name},type='l',xlab=\"Gene body percentile (5'->3')\",\
         ylab=\"Coverage\",lwd=0.8,col=\"blue\")"
    )?;
    writeln!(out, "dev.off()")?;

    Ok(())
}

// ===================================================================
// Tests
// ===================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_percentile_positions_matches_upstream_formula() {
        // 200 contiguous positions 1..=200: k = 199 * i / 100
        let bases: Vec<u64> = (1..=200).collect();
        let per = percentile_positions(&bases);
        assert_eq!(per.len(), 100);
        // i=1  -> k=1.99  -> interpolate between bases[1]=2 and bases[2]=3
        //         d0 = 2*(2-1.99)=0.02, d1 = 3*(1.99-1)=2.97 -> round(2.99)=3
        assert_eq!(per[0], 3);
        // i=100 -> k=199 (integral) -> bases[199] = 200
        assert_eq!(per[99], 200);
        // Monotonically non-decreasing
        assert!(per.windows(2).all(|w| w[0] <= w[1]));
    }

    #[test]
    fn test_round_half_to_even_matches_python() {
        // Python's round(): 0.5 -> 0, 1.5 -> 2, 2.5 -> 2, 3.5 -> 4
        assert_eq!(round_half_to_even(0.5), 0.0);
        assert_eq!(round_half_to_even(1.5), 2.0);
        assert_eq!(round_half_to_even(2.5), 2.0);
        assert_eq!(round_half_to_even(3.5), 4.0);
        // Non-half values round normally
        assert_eq!(round_half_to_even(2.4), 2.0);
        assert_eq!(round_half_to_even(2.6), 3.0);
    }

    #[test]
    fn test_percentile_positions_short_list_returned_unchanged() {
        let bases: Vec<u64> = (1..=50).collect();
        assert_eq!(percentile_positions(&bases), bases);
    }

    #[test]
    fn test_percentile_positions_spans_exon_junctions() {
        // Two exons: 1..=100 and 1001..=1100. Interpolated points may land in
        // the intron — upstream behaves the same way.
        let mut bases: Vec<u64> = (1..=100).collect();
        bases.extend(1001..=1100);
        let per = percentile_positions(&bases);
        assert_eq!(per.len(), 100);
        assert_eq!(per[99], 1100, "last percentile is the 3'-most base");
        assert!(per[0] < 100, "first percentile is inside the first exon");
    }

    /// Build a single-gene index and count coverage from a handful of reads.
    fn single_gene_index(strand: char) -> (IndexMap<String, Gene>, GeneBodyIndex) {
        use crate::gtf::{Exon, Transcript};
        let exons = vec![Exon {
            chrom: "chr1".to_string(),
            start: 1,
            end: 1000,
            strand,
        }];
        let gene = Gene {
            gene_id: "G1".to_string(),
            chrom: "chr1".to_string(),
            start: 1,
            end: 1000,
            strand,
            exons,
            effective_length: 1000,
            attributes: Default::default(),
            transcripts: vec![Transcript {
                transcript_id: "T1".to_string(),
                chrom: "chr1".to_string(),
                start: 1,
                end: 1000,
                strand,
                exons: vec![(1, 1000)],
                cds_start: None,
                cds_end: None,
            }],
        };
        let mut genes = IndexMap::new();
        genes.insert("G1".to_string(), gene);
        let index = GeneBodyIndex::from_genes(&genes);
        (genes, index)
    }

    #[test]
    fn test_coverage_profile_forward_strand() {
        use rust_htslib::bam::Read as BamRead;
        let (_genes, index) = single_gene_index('+');
        assert_eq!(index.num_transcripts(), 1);

        // One read covering positions 1..=50 (the 5' end of the transcript)
        let sam = "\
@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:20000\n\
r1\t0\tchr1\t1\t60\t50M\t*\t0\t0\t*\t*\n";
        let path = std::env::temp_dir().join(format!(
            "rustqc_gbc_fwd_{:?}.sam",
            std::thread::current().id()
        ));
        std::fs::write(&path, sam).unwrap();

        let mut reader = bam::Reader::from_path(&path).unwrap();
        let mut accum = GeneBodyCoverageAccum::new(&index);
        let mut record = bam::Record::new();
        while let Some(res) = reader.read(&mut record) {
            res.unwrap();
            accum.process_read(&record, "CHR1", &index);
        }
        let _ = std::fs::remove_file(&path);

        let result = accum.into_result(&index);
        let covered: u64 = result.coverage.iter().sum();
        assert!(covered > 0, "read should cover sampled positions");
        // Coverage must sit at the 5' end for a + strand transcript
        let first_half: u64 = result.coverage[..50].iter().sum();
        let second_half: u64 = result.coverage[50..].iter().sum();
        assert_eq!(second_half, 0, "no coverage expected in the 3' half");
        assert_eq!(first_half, covered);
    }

    #[test]
    fn test_coverage_profile_reverse_strand_is_flipped() {
        use rust_htslib::bam::Read as BamRead;
        let (_genes, index) = single_gene_index('-');

        // Same read at the low-coordinate end; on a - strand transcript that is
        // the 3' end, so coverage must land in the high bins.
        let sam = "\
@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:20000\n\
r1\t0\tchr1\t1\t60\t50M\t*\t0\t0\t*\t*\n";
        let path = std::env::temp_dir().join(format!(
            "rustqc_gbc_rev_{:?}.sam",
            std::thread::current().id()
        ));
        std::fs::write(&path, sam).unwrap();

        let mut reader = bam::Reader::from_path(&path).unwrap();
        let mut accum = GeneBodyCoverageAccum::new(&index);
        let mut record = bam::Record::new();
        while let Some(res) = reader.read(&mut record) {
            res.unwrap();
            accum.process_read(&record, "CHR1", &index);
        }
        let _ = std::fs::remove_file(&path);

        let result = accum.into_result(&index);
        let first_half: u64 = result.coverage[..50].iter().sum();
        let second_half: u64 = result.coverage[50..].iter().sum();
        assert_eq!(first_half, 0, "5' end is the high-coordinate end here");
        assert!(second_half > 0);
    }

    #[test]
    fn test_short_transcripts_are_skipped() {
        use crate::gtf::Exon;
        let gene = Gene {
            gene_id: "G_short".to_string(),
            chrom: "chr1".to_string(),
            start: 1,
            end: 50,
            strand: '+',
            exons: vec![Exon {
                chrom: "chr1".to_string(),
                start: 1,
                end: 50,
                strand: '+',
            }],
            effective_length: 50,
            attributes: Default::default(),
            transcripts: Vec::new(),
        };
        let mut genes = IndexMap::new();
        genes.insert("G_short".to_string(), gene);
        let index = GeneBodyIndex::from_genes(&genes);
        assert!(
            index.is_empty(),
            "50 bp transcript is below the 100 bp cutoff"
        );
    }
}
