//! Targeted genotyping at known SNP positions (NGSCheckMate).
//!
//! NGSCheckMate needs allele counts at ~10K known positions, which
//! `bcftools mpileup` currently supplies by building complete pileup columns
//! across the whole genome. Only the reads overlapping those positions matter,
//! so the work folds into an existing streaming pass for almost nothing: a
//! pointer advance per read, and CIGAR-resolved base extraction for the
//! ~0.003% of reads that actually overlap a site.

use std::collections::HashMap;
use std::io::BufRead;

use anyhow::{Context, Result};
use log::debug;
use rust_htslib::bam::{self, record::Cigar};

use crate::rna::bam_flags::{BAM_FDUP, BAM_FQCFAIL, BAM_FSECONDARY, BAM_FUNMAP};

/// One SNP site to genotype.
#[derive(Debug, Clone)]
pub struct SnpSite {
    /// 0-based reference position.
    pub pos: u64,
    /// Variant identifier (BED column 4), used as the VCF ID field.
    pub id: String,
    /// Reference allele (BED column 5).
    pub ref_base: u8,
    /// Alternate allele (BED column 6).
    pub alt_base: u8,
}

/// SNP sites grouped by chromosome, each list sorted by position.
#[derive(Debug, Default)]
pub struct SnpPanel {
    /// Chromosome name → sorted sites.
    pub by_chrom: HashMap<String, Vec<SnpSite>>,
    /// Total number of sites loaded.
    pub num_sites: usize,
}

/// Parse an NGSCheckMate SNP BED file (plain or gzip-compressed).
///
/// The file must have the 6-column NGSCheckMate layout:
/// `chrom  start  end  id  ref  alt`. The reference and alternate alleles are
/// required — without them, samples could not be compared against a common
/// set of alleles.
pub fn parse_snp_bed(path: &str) -> Result<SnpPanel> {
    let reader = crate::io::open_reader(path)
        .with_context(|| format!("Failed to open SNP BED file: {path}"))?;

    let mut by_chrom: HashMap<String, Vec<SnpSite>> = HashMap::new();
    let mut num_sites = 0usize;

    for (lineno, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("Failed to read line from SNP BED: {path}"))?;
        let trimmed = line.trim_end();
        if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.starts_with("track") {
            continue;
        }

        let fields: Vec<&str> = trimmed.split('\t').collect();
        anyhow::ensure!(
            fields.len() >= 6,
            "Malformed SNP BED '{}' at line {}: expected the 6-column NGSCheckMate layout \
             (chrom, start, end, id, ref, alt), found {} column(s)",
            path,
            lineno + 1,
            fields.len()
        );

        let end: u64 = fields[2].parse().with_context(|| {
            format!(
                "Malformed SNP BED '{}' at line {}: invalid end position '{}'",
                path,
                lineno + 1,
                fields[2]
            )
        })?;
        anyhow::ensure!(
            end > 0,
            "Malformed SNP BED '{}' at line {}: end position must be greater than 0",
            path,
            lineno + 1
        );

        let ref_base = single_base(fields[4], path, lineno + 1, "reference")?;
        let alt_base = single_base(fields[5], path, lineno + 1, "alternate")?;

        by_chrom
            .entry(fields[0].to_string())
            .or_default()
            .push(SnpSite {
                // BED is half-open, so the SNP base is the last one in the interval
                pos: end - 1,
                id: fields[3].to_string(),
                ref_base,
                alt_base,
            });
        num_sites += 1;
    }

    anyhow::ensure!(num_sites > 0, "No SNP sites found in BED file '{}'", path);

    for sites in by_chrom.values_mut() {
        sites.sort_unstable_by_key(|s| s.pos);
    }

    debug!("Loaded {num_sites} SNP sites from {path}");
    Ok(SnpPanel {
        by_chrom,
        num_sites,
    })
}

/// Validate and upper-case a single-base allele field.
fn single_base(field: &str, path: &str, lineno: usize, which: &str) -> Result<u8> {
    let bytes = field.as_bytes();
    anyhow::ensure!(
        bytes.len() == 1 && matches!(bytes[0].to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'),
        "Malformed SNP BED '{}' at line {}: {} allele must be a single A/C/G/T base, found '{}'",
        path,
        lineno,
        which,
        field
    );
    Ok(bytes[0].to_ascii_uppercase())
}

// ===================================================================
// Allele counting
// ===================================================================

/// Per-site allele counts.
#[derive(Debug, Clone, Copy, Default)]
pub struct AlleleCounts {
    /// Reads supporting the reference allele.
    pub ref_count: u32,
    /// Reads supporting the alternate allele.
    pub alt_count: u32,
    /// Reads carrying any other base at the site.
    pub other_count: u32,
}

impl AlleleCounts {
    /// Total reads observed at the site.
    pub fn depth(&self) -> u32 {
        self.ref_count + self.alt_count + self.other_count
    }

    /// Alternate allele fraction over ref + alt reads.
    pub fn alt_fraction(&self) -> f64 {
        let informative = self.ref_count + self.alt_count;
        if informative == 0 {
            0.0
        } else {
            self.alt_count as f64 / informative as f64
        }
    }

    /// Genotype call from the alternate allele fraction.
    ///
    /// NGSCheckMate only needs to tell hom-ref, het and hom-alt apart, so the
    /// simple fraction thresholds used here stand in for a full genotype
    /// likelihood model.
    pub fn genotype(&self) -> &'static str {
        if self.depth() == 0 {
            return "./.";
        }
        let frac = self.alt_fraction();
        if frac < 0.15 {
            "0/0"
        } else if frac <= 0.85 {
            "0/1"
        } else {
            "1/1"
        }
    }
}

/// Accumulates allele counts at the panel's sites during a streaming pass.
#[derive(Debug)]
pub struct SnpAccum {
    /// Counts per chromosome, parallel to `SnpPanel::by_chrom` entries.
    counts: HashMap<String, Vec<AlleleCounts>>,
    /// Minimum base quality for a base to be counted.
    min_base_quality: u8,
    /// Minimum MAPQ for a read to be considered.
    min_mapq: u8,
}

impl SnpAccum {
    /// Create an accumulator for the given panel.
    pub fn new(panel: &SnpPanel, min_mapq: u8, min_base_quality: u8) -> Self {
        let counts = panel
            .by_chrom
            .iter()
            .map(|(chrom, sites)| (chrom.clone(), vec![AlleleCounts::default(); sites.len()]))
            .collect();
        Self {
            counts,
            min_base_quality,
            min_mapq,
        }
    }

    /// Add one alignment record's base calls at any overlapping SNP sites.
    ///
    /// Records that are unmapped, secondary, QC-failed, duplicate-flagged or
    /// below the MAPQ cutoff are ignored.
    pub fn process_read(&mut self, record: &bam::Record, chrom: &str, panel: &SnpPanel) {
        let flags = record.flags();
        if flags & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) != 0 {
            return;
        }
        if record.mapq() < self.min_mapq {
            return;
        }

        let (Some(sites), Some(counts)) = (panel.by_chrom.get(chrom), self.counts.get_mut(chrom))
        else {
            return;
        };
        if sites.is_empty() {
            return;
        }

        let start = record.pos() as u64;
        // Fast reject: almost every read lands past the current site pointer
        let first = sites.partition_point(|s| s.pos < start);
        if first >= sites.len() {
            return;
        }

        let seq = record.seq();
        let quals = record.qual();
        let mut ref_pos = start;
        let mut read_pos = 0usize;

        for op in record.cigar().iter() {
            match *op {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    let len = len as u64;
                    let block_end = ref_pos + len;
                    for (i, site) in sites.iter().enumerate().skip(first) {
                        if site.pos >= block_end {
                            break;
                        }
                        if site.pos < ref_pos {
                            continue;
                        }
                        let offset = (site.pos - ref_pos) as usize;
                        let idx = read_pos + offset;
                        if idx >= seq.len() {
                            continue;
                        }
                        if quals.get(idx).copied().unwrap_or(0) < self.min_base_quality {
                            continue;
                        }
                        let base = seq[idx].to_ascii_uppercase();
                        let entry = &mut counts[i];
                        if base == site.ref_base {
                            entry.ref_count += 1;
                        } else if base == site.alt_base {
                            entry.alt_count += 1;
                        } else {
                            entry.other_count += 1;
                        }
                    }
                    ref_pos = block_end;
                    read_pos += len as usize;
                }
                Cigar::Ins(len) | Cigar::SoftClip(len) => read_pos += len as usize,
                Cigar::Del(len) | Cigar::RefSkip(len) => ref_pos += len as u64,
                _ => {}
            }
        }
    }

    /// Counts for one chromosome, in panel order.
    pub fn counts_for(&self, chrom: &str) -> Option<&[AlleleCounts]> {
        self.counts.get(chrom).map(|v| v.as_slice())
    }

    /// Number of sites with at least one supporting read.
    pub fn covered_sites(&self) -> usize {
        self.counts
            .values()
            .flat_map(|v| v.iter())
            .filter(|c| c.depth() > 0)
            .count()
    }
}

// ===================================================================
// Tests
// ===================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::Read as BamRead;

    fn write_temp(content: &str, ext: &str) -> std::path::PathBuf {
        use std::sync::atomic::{AtomicU64, Ordering};
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        let id = COUNTER.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "rustqc_snp_test_{:?}_{}.{}",
            std::thread::current().id(),
            id,
            ext
        ));
        std::fs::write(&path, content).unwrap();
        path
    }

    #[test]
    fn test_parse_snp_bed() {
        let path = write_temp(
            "chr1\t99\t100\trs1\tA\tG\n\
             chr1\t199\t200\trs2\tC\tT\n\
             chr2\t9\t10\trs3\tG\tA\n",
            "bed",
        );
        let panel = parse_snp_bed(path.to_str().unwrap()).unwrap();
        let _ = std::fs::remove_file(&path);

        assert_eq!(panel.num_sites, 3);
        let chr1 = &panel.by_chrom["chr1"];
        assert_eq!(chr1.len(), 2);
        assert_eq!(chr1[0].pos, 99, "BED end - 1 is the 0-based SNP position");
        assert_eq!(chr1[0].ref_base, b'A');
        assert_eq!(chr1[0].alt_base, b'G');
        assert_eq!(chr1[0].id, "rs1");
    }

    #[test]
    fn test_parse_snp_bed_requires_alleles() {
        let path = write_temp("chr1\t99\t100\trs1\n", "bed");
        let err = parse_snp_bed(path.to_str().unwrap()).unwrap_err();
        let _ = std::fs::remove_file(&path);
        assert!(
            err.to_string().contains("6-column NGSCheckMate layout"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn test_allele_counting_and_genotype() {
        // SNP at 0-based 100 (1-based 101), ref A alt G
        let bed = write_temp("chr1\t100\t101\trs1\tA\tG\n", "bed");
        let panel = parse_snp_bed(bed.to_str().unwrap()).unwrap();
        let _ = std::fs::remove_file(&bed);

        // Reads start at 1-based 96, so the SNP is the 6th base of each read.
        // r1/r2 carry G (alt), r3 carries A (ref), r4 is a duplicate and is
        // ignored, r5 has MAPQ 5 and is ignored.
        let sam = "\
@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:20000\n\
r1\t0\tchr1\t96\t60\t10M\t*\t0\t0\tTTTTTGTTTT\tIIIIIIIIII\n\
r2\t0\tchr1\t96\t60\t10M\t*\t0\t0\tTTTTTGTTTT\tIIIIIIIIII\n\
r3\t0\tchr1\t96\t60\t10M\t*\t0\t0\tTTTTTATTTT\tIIIIIIIIII\n\
r4\t1024\tchr1\t96\t60\t10M\t*\t0\t0\tTTTTTGTTTT\tIIIIIIIIII\n\
r5\t0\tchr1\t96\t5\t10M\t*\t0\t0\tTTTTTGTTTT\tIIIIIIIIII\n";
        let sam_path = write_temp(sam, "sam");

        let mut reader = bam::Reader::from_path(&sam_path).unwrap();
        let mut accum = SnpAccum::new(&panel, 30, 13);
        let mut record = bam::Record::new();
        while let Some(res) = reader.read(&mut record) {
            res.unwrap();
            accum.process_read(&record, "chr1", &panel);
        }
        let _ = std::fs::remove_file(&sam_path);

        let counts = accum.counts_for("chr1").unwrap();
        assert_eq!(counts[0].ref_count, 1, "one ref-supporting read");
        assert_eq!(counts[0].alt_count, 2, "two alt-supporting reads");
        assert_eq!(counts[0].other_count, 0);
        assert_eq!(counts[0].depth(), 3);
        assert!((counts[0].alt_fraction() - 2.0 / 3.0).abs() < 1e-9);
        assert_eq!(counts[0].genotype(), "0/1");
        assert_eq!(accum.covered_sites(), 1);
    }

    #[test]
    fn test_genotype_thresholds() {
        let hom_ref = AlleleCounts {
            ref_count: 100,
            alt_count: 5,
            other_count: 0,
        };
        assert_eq!(hom_ref.genotype(), "0/0");

        let het = AlleleCounts {
            ref_count: 50,
            alt_count: 50,
            other_count: 0,
        };
        assert_eq!(het.genotype(), "0/1");

        let hom_alt = AlleleCounts {
            ref_count: 2,
            alt_count: 100,
            other_count: 0,
        };
        assert_eq!(hom_alt.genotype(), "1/1");

        assert_eq!(AlleleCounts::default().genotype(), "./.");
    }

    #[test]
    fn test_soft_clip_and_deletion_offsets() {
        // Read at 1-based 100 (0-based 99) with CIGAR 5S5M2D5M and sequence
        // CCCCC TTTTG AAAAA:
        //   5S  -> read bases 0..5, no reference
        //   5M  -> reference 99..103, read bases 5..10 ("TTTTG")
        //   2D  -> reference 104..105, no read bases
        //   5M  -> reference 106..110, read bases 10..15 ("AAAAA")
        // rs1 at 0-based 103 must pick up the 'G' (alt) despite the soft clip;
        // rs2 at 0-based 104 falls inside the deletion and must count nothing.
        let bed = write_temp(
            "chr1\t103\t104\trs1\tA\tG\n\
             chr1\t104\t105\trs2\tA\tG\n",
            "bed",
        );
        let panel = parse_snp_bed(bed.to_str().unwrap()).unwrap();
        let _ = std::fs::remove_file(&bed);

        // Reference positions 99..103 are the first 5M (read bases 5..9);
        // position 104 is the 5th base of that block -> read base index 9 = 'G'
        let sam = "\
@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:20000\n\
r1\t0\tchr1\t100\t60\t5S5M2D5M\t*\t0\t0\tCCCCCTTTTGAAAAA\tIIIIIIIIIIIIIII\n";
        let sam_path = write_temp(sam, "sam");

        let mut reader = bam::Reader::from_path(&sam_path).unwrap();
        let mut accum = SnpAccum::new(&panel, 0, 0);
        let mut record = bam::Record::new();
        while let Some(res) = reader.read(&mut record) {
            res.unwrap();
            accum.process_read(&record, "chr1", &panel);
        }
        let _ = std::fs::remove_file(&sam_path);

        let counts = accum.counts_for("chr1").unwrap();
        assert_eq!(
            counts[0].alt_count, 1,
            "soft clip must not shift the offset"
        );
        assert_eq!(counts[0].ref_count, 0);
        assert_eq!(counts[1].depth(), 0, "a site inside a deletion has no base");
    }
}
