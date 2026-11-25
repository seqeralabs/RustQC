//! Samtools flagstat-compatible output.
//!
//! Produces output matching `samtools flagstat` format, allowing MultiQC and
//! other downstream tools to parse the results as if they came from samtools.

use std::path::Path;

use anyhow::{Context, Result};
use log::debug;

use super::bam_stat::BamStatResult;

use crate::rna::bam_flags::*;

// ============================================================================
// Counters
// ============================================================================

/// samtools flagstat counters, each indexed `[QC-passed, QC-failed]`.
///
/// A direct port of `flagstat_loop` in samtools `bam_stat.c`: every record is
/// counted, split by its QC-fail flag (0x200).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct FlagstatCounts {
    /// All records.
    pub n_reads: [u64; 2],
    /// Primary records (neither secondary nor supplementary).
    pub n_primary: [u64; 2],
    /// Secondary records (0x100).
    pub n_secondary: [u64; 2],
    /// Supplementary records (0x800, not secondary).
    pub n_supp: [u64; 2],
    /// Duplicate records (0x400).
    pub n_dup: [u64; 2],
    /// Primary duplicate records.
    pub n_pdup: [u64; 2],
    /// Mapped records.
    pub n_mapped: [u64; 2],
    /// Primary mapped records.
    pub n_pmapped: [u64; 2],
    /// Primary records paired in sequencing (0x1).
    pub n_pair_all: [u64; 2],
    /// Primary paired records that are read1 (0x40).
    pub n_read1: [u64; 2],
    /// Primary paired records that are read2 (0x80).
    pub n_read2: [u64; 2],
    /// Primary paired, mapped records with the proper-pair flag (0x2).
    pub n_pair_good: [u64; 2],
    /// Primary paired records with itself and mate mapped.
    pub n_pair_map: [u64; 2],
    /// Primary paired records mapped with the mate unmapped.
    pub n_sgltn: [u64; 2],
    /// `n_pair_map` records whose mate is on a different reference.
    pub n_diffchr: [u64; 2],
    /// `n_diffchr` records with MAPQ >= 5.
    pub n_diffhigh: [u64; 2],
}

impl FlagstatCounts {
    /// Count one record (samtools `flagstat_loop`).
    pub fn add(&mut self, flags: u16, tid: i32, mtid: i32, mapq: u8) {
        let w = usize::from(flags & BAM_FQCFAIL != 0);
        let unmapped = flags & BAM_FUNMAP != 0;
        self.n_reads[w] += 1;
        if flags & BAM_FSECONDARY != 0 {
            self.n_secondary[w] += 1;
        } else if flags & BAM_FSUPPLEMENTARY != 0 {
            self.n_supp[w] += 1;
        } else {
            self.n_primary[w] += 1;
            if flags & BAM_FPAIRED != 0 {
                self.n_pair_all[w] += 1;
                if flags & BAM_FPROPER_PAIR != 0 && !unmapped {
                    self.n_pair_good[w] += 1;
                }
                if flags & BAM_FREAD1 != 0 {
                    self.n_read1[w] += 1;
                }
                if flags & BAM_FREAD2 != 0 {
                    self.n_read2[w] += 1;
                }
                let mate_unmapped = flags & BAM_FMUNMAP != 0;
                if mate_unmapped && !unmapped {
                    self.n_sgltn[w] += 1;
                }
                if !unmapped && !mate_unmapped {
                    self.n_pair_map[w] += 1;
                    if mtid != tid {
                        self.n_diffchr[w] += 1;
                        if mapq >= 5 {
                            self.n_diffhigh[w] += 1;
                        }
                    }
                }
            }
            if !unmapped {
                self.n_pmapped[w] += 1;
            }
            if flags & BAM_FDUP != 0 {
                self.n_pdup[w] += 1;
            }
        }
        if !unmapped {
            self.n_mapped[w] += 1;
        }
        if flags & BAM_FDUP != 0 {
            self.n_dup[w] += 1;
        }
    }

    /// Add another set of counters into this one.
    pub fn merge(&mut self, other: &FlagstatCounts) {
        let pairs: [(&mut [u64; 2], &[u64; 2]); 16] = [
            (&mut self.n_reads, &other.n_reads),
            (&mut self.n_primary, &other.n_primary),
            (&mut self.n_secondary, &other.n_secondary),
            (&mut self.n_supp, &other.n_supp),
            (&mut self.n_dup, &other.n_dup),
            (&mut self.n_pdup, &other.n_pdup),
            (&mut self.n_mapped, &other.n_mapped),
            (&mut self.n_pmapped, &other.n_pmapped),
            (&mut self.n_pair_all, &other.n_pair_all),
            (&mut self.n_read1, &other.n_read1),
            (&mut self.n_read2, &other.n_read2),
            (&mut self.n_pair_good, &other.n_pair_good),
            (&mut self.n_pair_map, &other.n_pair_map),
            (&mut self.n_sgltn, &other.n_sgltn),
            (&mut self.n_diffchr, &other.n_diffchr),
            (&mut self.n_diffhigh, &other.n_diffhigh),
        ];
        for (a, b) in pairs {
            a[0] += b[0];
            a[1] += b[1];
        }
    }
}

// ============================================================================
// Output formatting
// ============================================================================

/// samtools `percent()`: `"%.2f%%"` of `(float)n / total * 100.0`, or `N/A`.
///
/// The division is done in single precision, as in samtools.
fn percent(n: u64, total: u64) -> String {
    if total == 0 {
        "N/A".to_string()
    } else {
        let pct = f64::from(n as f32 / total as f32) * 100.0;
        format!("{pct:.2}%")
    }
}

/// Write samtools flagstat-compatible output (samtools `print_flagstat` text
/// format: `<QC-passed> + <QC-failed> <description>` per line).
///
/// # Arguments
/// * `result` - The computed BAM statistics
/// * `output_path` - Path to write the flagstat file
pub fn write_flagstat(result: &BamStatResult, output_path: &Path) -> Result<()> {
    use std::io::Write;

    let mut out = std::fs::File::create(output_path)
        .with_context(|| format!("Failed to create flagstat file: {}", output_path.display()))?;
    out.write_all(format_flagstat(&result.flagstat).as_bytes())?;

    debug!("Wrote flagstat output to {}", output_path.display());
    Ok(())
}

/// Render flagstat counters in samtools' default text format.
fn format_flagstat(s: &FlagstatCounts) -> String {
    let line = |v: &[u64; 2], desc: &str| format!("{} + {} {desc}\n", v[0], v[1]);
    let line_pct = |v: &[u64; 2], total: &[u64; 2], desc: &str| {
        format!(
            "{} + {} {desc} ({} : {})\n",
            v[0],
            v[1],
            percent(v[0], total[0]),
            percent(v[1], total[1])
        )
    };
    let mut o = String::new();
    o += &line(&s.n_reads, "in total (QC-passed reads + QC-failed reads)");
    o += &line(&s.n_primary, "primary");
    o += &line(&s.n_secondary, "secondary");
    o += &line(&s.n_supp, "supplementary");
    o += &line(&s.n_dup, "duplicates");
    o += &line(&s.n_pdup, "primary duplicates");
    o += &line_pct(&s.n_mapped, &s.n_reads, "mapped");
    o += &line_pct(&s.n_pmapped, &s.n_primary, "primary mapped");
    o += &line(&s.n_pair_all, "paired in sequencing");
    o += &line(&s.n_read1, "read1");
    o += &line(&s.n_read2, "read2");
    o += &line_pct(&s.n_pair_good, &s.n_pair_all, "properly paired");
    o += &line(&s.n_pair_map, "with itself and mate mapped");
    o += &line_pct(&s.n_sgltn, &s.n_pair_all, "singletons");
    o += &line(&s.n_diffchr, "with mate mapped to a different chr");
    o += &line(
        &s.n_diffhigh,
        "with mate mapped to a different chr (mapQ>=5)",
    );
    o
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
    fn test_flagstat_format() {
        let mut reader =
            bam::Reader::from_path("tests/data/test.bam").expect("Failed to open test.bam");
        let mut accum = BamStatAccum::default();
        let mut record = bam::Record::new();

        while let Some(res) = reader.read(&mut record) {
            res.expect("Error reading BAM record");
            accum.process_read(&record, 30);
        }

        let result = accum.into_result();
        let tmp_path = std::env::temp_dir().join("rustqc_test_flagstat.txt");
        write_flagstat(&result, &tmp_path).expect("Failed to write flagstat");

        let mut contents = String::new();
        std::fs::File::open(&tmp_path)
            .unwrap()
            .read_to_string(&mut contents)
            .unwrap();
        let _ = std::fs::remove_file(&tmp_path);

        // Verify first line format
        assert!(
            contents.contains("in total (QC-passed reads + QC-failed reads)"),
            "Missing total line in flagstat output"
        );
        // Verify mapped line has percentage
        assert!(
            contents.contains("mapped"),
            "Missing mapped line in flagstat output"
        );
    }

    #[test]
    fn percent_divides_in_single_precision() {
        // samtools: (float)n / total * 100.0. In double 1/160 prints 0.62%.
        assert_eq!(percent(1, 160), "0.63%");
        assert_eq!(percent(7, 160), "4.37%");
        assert_eq!(percent(0, 0), "N/A");
        assert_eq!(percent(5, 5), "100.00%");
    }

    #[test]
    fn counts_are_split_by_qc_fail() {
        let mut c = FlagstatCounts::default();
        let pe = BAM_FPAIRED | BAM_FREAD1;
        c.add(pe | BAM_FPROPER_PAIR, 0, 0, 60); // passed, proper pair
        c.add(pe | BAM_FQCFAIL | BAM_FPROPER_PAIR, 0, 1, 60); // failed, mate elsewhere
        c.add(pe | BAM_FUNMAP | BAM_FPROPER_PAIR, 0, 0, 0); // unmapped: not proper
        c.add(pe | BAM_FMUNMAP | BAM_FDUP, 0, 0, 3); // singleton duplicate
        c.add(BAM_FSECONDARY | BAM_FDUP, 0, 0, 0);
        c.add(BAM_FSUPPLEMENTARY | BAM_FQCFAIL, 0, 0, 0);
        let text = format_flagstat(&c);
        let expected = "\
4 + 2 in total (QC-passed reads + QC-failed reads)
3 + 1 primary
1 + 0 secondary
0 + 1 supplementary
2 + 0 duplicates
1 + 0 primary duplicates
3 + 2 mapped (75.00% : 100.00%)
2 + 1 primary mapped (66.67% : 100.00%)
3 + 1 paired in sequencing
3 + 1 read1
0 + 0 read2
1 + 1 properly paired (33.33% : 100.00%)
1 + 1 with itself and mate mapped
1 + 0 singletons (33.33% : 0.00%)
0 + 1 with mate mapped to a different chr
0 + 1 with mate mapped to a different chr (mapQ>=5)
";
        assert_eq!(text, expected);

        let mut merged = FlagstatCounts::default();
        merged.merge(&c);
        merged.merge(&c);
        assert_eq!(merged.n_reads, [8, 4]);
        assert_eq!(merged.n_diffhigh, [0, 2]);
    }
}
