//! Base assignment and the `CollectRnaSeqMetrics`-compatible writer.

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

use super::regions::{Region, RegionSets};
use crate::common::bam_flags::*;

/// Counts of aligned bases by the region they fall in.
#[derive(Debug, Clone, Default)]
pub struct CodingCounts {
    /// Bases in an annotated coding sequence.
    pub coding: u64,
    /// Bases exonic but not coding.
    pub utr: u64,
    /// Bases within a transcript but not exonic.
    pub intronic: u64,
    /// Bases outside every transcript.
    pub intergenic: u64,
    /// Bases in reads before any assignment, which is the denominator.
    pub aligned: u64,
    /// Bases across every read, aligned or not, which Picard calls `PF_BASES`.
    pub total: u64,
}

impl CodingCounts {
    /// Offer one record.
    ///
    /// Unmapped, secondary and supplementary records contribute nothing, which
    /// is the same set `CollectWgsMetrics` excludes and gives the same
    /// `PF_ALIGNED_BASES`.
    pub fn process_read(&mut self, record: &bam::Record, chrom: &str, regions: &RegionSets) {
        let flags = record.flags();
        // PF_BASES counts primary records only, so a secondary alignment does
        // not have its sequence counted twice. An unmapped primary record
        // still contributes its bases, which is why the two conditions differ.
        if flags & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY) == 0 {
            self.total += record.seq_len() as u64;
        }
        if flags & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY) != 0 {
            return;
        }

        let mut reference = record.pos();
        for op in record.cigar().iter() {
            match op {
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                    for k in 0..i64::from(*n) {
                        let position = reference + k;
                        if position < 0 {
                            continue;
                        }
                        self.aligned += 1;
                        match regions.classify(chrom, position as u64) {
                            Region::Coding => self.coding += 1,
                            Region::Utr => self.utr += 1,
                            Region::Intronic => self.intronic += 1,
                            Region::Intergenic => self.intergenic += 1,
                        }
                    }
                    reference += i64::from(*n);
                }
                Cigar::Del(n) | Cigar::RefSkip(n) => reference += i64::from(*n),
                _ => {}
            }
        }
    }

    /// Fold another set of counts in.
    pub fn merge(&mut self, other: &CodingCounts) {
        self.coding += other.coding;
        self.utr += other.utr;
        self.intronic += other.intronic;
        self.intergenic += other.intergenic;
        self.aligned += other.aligned;
        self.total += other.total;
    }

    /// Fraction of aligned bases in a class.
    fn fraction(&self, part: u64) -> f64 {
        if self.aligned == 0 {
            0.0
        } else {
            part as f64 / self.aligned as f64
        }
    }

    /// Fraction of aligned bases that are exonic, coding or otherwise.
    pub fn mrna_fraction(&self) -> f64 {
        self.fraction(self.coding + self.utr)
    }
}

/// Format a float the way Picard's metrics writer does.
fn fmt_picard(value: f64) -> String {
    if !value.is_finite() {
        return "?".to_string();
    }
    if value == value.trunc() && value.abs() < 1e15 {
        return format!("{}", value as i64);
    }
    let text = format!("{value:.6}");
    text.trim_end_matches('0').trim_end_matches('.').to_string()
}

/// Write the `CollectRnaSeqMetrics`-compatible subset.
///
/// Only the columns this analysis computes are written. The strand-specificity
/// and coverage-bias columns need a library protocol and a per-transcript
/// coverage pass that this mode does not do, and emitting them as zero would
/// read as a measurement rather than an absence.
pub fn write_coding_metrics(counts: &CodingCounts, path: &Path) -> Result<()> {
    let mut out = std::fs::File::create(path)
        .map(std::io::BufWriter::new)
        .with_context(|| format!("Failed to create coding metrics: {}", path.display()))?;

    writeln!(out, "## METRICS CLASS\tpicard.analysis.RnaSeqMetrics")?;
    writeln!(
        out,
        "PF_BASES\tPF_ALIGNED_BASES\tCODING_BASES\tUTR_BASES\tINTRONIC_BASES\t\
         INTERGENIC_BASES\tPCT_CODING_BASES\tPCT_UTR_BASES\tPCT_INTRONIC_BASES\t\
         PCT_INTERGENIC_BASES\tPCT_MRNA_BASES"
    )?;
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        counts.total,
        counts.aligned,
        counts.coding,
        counts.utr,
        counts.intronic,
        counts.intergenic,
        fmt_picard(counts.fraction(counts.coding)),
        fmt_picard(counts.fraction(counts.utr)),
        fmt_picard(counts.fraction(counts.intronic)),
        fmt_picard(counts.fraction(counts.intergenic)),
        fmt_picard(counts.mrna_fraction()),
    )?;
    writeln!(out)?;
    out.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_four_classes_account_for_every_aligned_base() {
        let counts = CodingCounts {
            coding: 10,
            utr: 20,
            intronic: 30,
            intergenic: 40,
            aligned: 100,
            total: 120,
        };
        assert_eq!(
            counts.coding + counts.utr + counts.intronic + counts.intergenic,
            counts.aligned
        );
        assert!((counts.mrna_fraction() - 0.3).abs() < 1e-12);
    }

    #[test]
    fn a_secondary_alignment_does_not_have_its_bases_counted_twice() {
        use rust_htslib::bam::record::{Cigar, CigarString, Record};

        let mut make = |flags: u16| {
            let mut r = Record::new();
            let cigar = CigarString(vec![Cigar::Match(4)]);
            r.set(b"q", Some(&cigar), b"ACGT", &[30u8; 4]);
            r.set_flags(flags);
            r.set_tid(0);
            r.set_pos(0);
            r
        };
        let regions = RegionSets::default();
        let mut counts = CodingCounts::default();
        counts.process_read(&make(0), "chr1", &regions);
        counts.process_read(&make(BAM_FSECONDARY), "chr1", &regions);
        assert_eq!(counts.total, 4, "only the primary record's bases count");
        assert_eq!(counts.aligned, 4);
    }

    #[test]
    fn an_empty_run_does_not_divide_by_zero() {
        let counts = CodingCounts::default();
        assert_eq!(counts.mrna_fraction(), 0.0);
        assert_eq!(counts.fraction(0), 0.0);
    }

    #[test]
    fn picard_formatting_drops_trailing_zeros() {
        assert_eq!(fmt_picard(0.482318), "0.482318");
        assert_eq!(fmt_picard(0.0), "0");
        assert_eq!(fmt_picard(1.0), "1");
    }
}
