//! Samtools flagstat-compatible output.
//!
//! Produces output matching `samtools flagstat` format, allowing MultiQC and
//! other downstream tools to parse the results as if they came from samtools.

use std::path::Path;

use anyhow::{Context, Result};
use log::debug;

use super::bam_stat::BamStatResult;

// ============================================================================
// Output formatting
// ============================================================================

/// Format a flagstat count line with the standard `<passed> + 0 <description>` format.
fn fmt_line(count: u64, description: &str) -> String {
    format!("{count} + 0 {description}")
}

/// Format a flagstat percentage as `(XX.XX% : N/A)`.
///
/// samtools uses `:` separated `passed : failed` percentages. Since we do not
/// separate QC-pass/fail, the failed side is always `N/A`.
fn fmt_pct(numerator: u64, denominator: u64) -> String {
    if denominator == 0 {
        "(N/A : N/A)".to_string()
    } else {
        let pct = 100.0 * numerator as f64 / denominator as f64;
        format!("({pct:.2}% : N/A)")
    }
}

/// Write samtools flagstat-compatible output.
///
/// The output format matches `samtools flagstat` exactly, with one line per
/// statistic. Each line has the form `<passed> + 0 <description>`.
///
/// # Arguments
/// * `result` - The computed BAM statistics
/// * `output_path` - Path to write the flagstat file
pub fn write_flagstat(result: &BamStatResult, output_path: &Path) -> Result<()> {
    use std::io::Write;

    let mut out = std::fs::File::create(output_path)
        .with_context(|| format!("Failed to create flagstat file: {}", output_path.display()))?;

    // Total = primary + secondary + supplementary (all QC states)
    let total = result.total_records;

    // Primary count for percentage denominators
    let primary = result.primary_count;

    // Line 1: total (QC-passed + QC-failed)
    writeln!(
        out,
        "{}",
        fmt_line(total, "in total (QC-passed reads + QC-failed reads)")
    )?;

    // Line 2: primary
    writeln!(out, "{}", fmt_line(primary, "primary"))?;

    // Line 3: secondary
    writeln!(out, "{}", fmt_line(result.secondary, "secondary"))?;

    // Line 4: supplementary
    writeln!(out, "{}", fmt_line(result.supplementary, "supplementary"))?;

    // Line 5: duplicates
    writeln!(out, "{}", fmt_line(result.duplicates, "duplicates"))?;

    // Line 6: primary duplicates
    writeln!(
        out,
        "{}",
        fmt_line(result.primary_duplicates, "primary duplicates")
    )?;

    // Line 7: mapped
    writeln!(
        out,
        "{} {}",
        fmt_line(result.mapped, "mapped"),
        fmt_pct(result.mapped, total)
    )?;

    // Line 8: primary mapped
    writeln!(
        out,
        "{} {}",
        fmt_line(result.primary_mapped, "primary mapped"),
        fmt_pct(result.primary_mapped, primary)
    )?;

    // Line 9: paired in sequencing
    writeln!(
        out,
        "{}",
        fmt_line(result.paired_flagstat, "paired in sequencing")
    )?;

    // Line 10: read1
    writeln!(out, "{}", fmt_line(result.read1_flagstat, "read1"))?;

    // Line 11: read2
    writeln!(out, "{}", fmt_line(result.read2_flagstat, "read2"))?;

    // Line 12: properly paired (percentage of paired reads)
    writeln!(
        out,
        "{} {}",
        fmt_line(result.properly_paired, "properly paired"),
        fmt_pct(result.properly_paired, result.paired_flagstat)
    )?;

    // Line 13: with itself and mate mapped
    writeln!(
        out,
        "{}",
        fmt_line(result.both_mapped, "with itself and mate mapped")
    )?;

    // Line 14: singletons (percentage of paired reads)
    writeln!(
        out,
        "{} {}",
        fmt_line(result.singletons, "singletons"),
        fmt_pct(result.singletons, result.paired_flagstat)
    )?;

    // Line 15: with mate mapped to a different chr
    writeln!(
        out,
        "{}",
        fmt_line(result.mate_diff_chr, "with mate mapped to a different chr")
    )?;

    // Line 16: with mate mapped to a different chr (mapQ>=5)
    writeln!(
        out,
        "{}",
        fmt_line(
            result.mate_diff_chr_mapq5,
            "with mate mapped to a different chr (mapQ>=5)"
        )
    )?;

    debug!("Wrote flagstat output to {}", output_path.display());
    Ok(())
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rna::bam_flags::*;
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

    // ------------------------------------------------------------------------
    // Oracle: samtools 1.22.1 `flagstat_loop` and `out_fmt_default`
    // (bam_stat.c:48-79 and 124-143), transcribed.
    // ------------------------------------------------------------------------

    /// samtools' `bam_flagstat_t`; index 0 is QC-passed, 1 is QC-failed.
    #[derive(Debug, Default)]
    struct SamtoolsFlagstat {
        n_reads: [u64; 2],
        n_mapped: [u64; 2],
        n_pair_all: [u64; 2],
        n_pair_map: [u64; 2],
        n_pair_good: [u64; 2],
        n_sgltn: [u64; 2],
        n_read1: [u64; 2],
        n_read2: [u64; 2],
        n_dup: [u64; 2],
        n_diffchr: [u64; 2],
        n_diffhigh: [u64; 2],
        n_secondary: [u64; 2],
        n_supp: [u64; 2],
        n_primary: [u64; 2],
        n_pmapped: [u64; 2],
        n_pdup: [u64; 2],
    }

    impl SamtoolsFlagstat {
        fn add(&mut self, flag: u16, tid: i32, mtid: i32, qual: u8) {
            let w = usize::from(flag & BAM_FQCFAIL != 0);
            self.n_reads[w] += 1;
            if flag & BAM_FSECONDARY != 0 {
                self.n_secondary[w] += 1;
            } else if flag & BAM_FSUPPLEMENTARY != 0 {
                self.n_supp[w] += 1;
            } else {
                self.n_primary[w] += 1;
                if flag & BAM_FPAIRED != 0 {
                    self.n_pair_all[w] += 1;
                    if flag & BAM_FPROPER_PAIR != 0 && flag & BAM_FUNMAP == 0 {
                        self.n_pair_good[w] += 1;
                    }
                    if flag & BAM_FREAD1 != 0 {
                        self.n_read1[w] += 1;
                    }
                    if flag & BAM_FREAD2 != 0 {
                        self.n_read2[w] += 1;
                    }
                    if flag & BAM_FMUNMAP != 0 && flag & BAM_FUNMAP == 0 {
                        self.n_sgltn[w] += 1;
                    }
                    if flag & BAM_FUNMAP == 0 && flag & BAM_FMUNMAP == 0 {
                        self.n_pair_map[w] += 1;
                        if mtid != tid {
                            self.n_diffchr[w] += 1;
                            if qual >= 5 {
                                self.n_diffhigh[w] += 1;
                            }
                        }
                    }
                }
                if flag & BAM_FUNMAP == 0 {
                    self.n_pmapped[w] += 1;
                }
                if flag & BAM_FDUP != 0 {
                    self.n_pdup[w] += 1;
                }
            }
            if flag & BAM_FUNMAP == 0 {
                self.n_mapped[w] += 1;
            }
            if flag & BAM_FDUP != 0 {
                self.n_dup[w] += 1;
            }
            // Not samtools: RustQC also counts a 0x900 record as supplementary.
            // Which behaviour is right is open (#125), so the oracle keeps
            // RustQC's for that one line.
            if flag & BAM_FSECONDARY != 0 && flag & BAM_FSUPPLEMENTARY != 0 {
                self.n_supp[w] += 1;
            }
        }

        fn render(&self) -> String {
            // samtools computes `(float)n / total * 100.0`.
            fn percent(n: u64, total: u64) -> String {
                if total == 0 {
                    "N/A".to_string()
                } else {
                    format!("{:.2}%", (n as f32 / total as f32) as f64 * 100.0)
                }
            }
            let line = |c: [u64; 2], d: &str| format!("{} + {} {d}\n", c[0], c[1]);
            let pct = |c: [u64; 2], t: [u64; 2], d: &str| {
                format!(
                    "{} + {} {d} ({} : {})\n",
                    c[0],
                    c[1],
                    percent(c[0], t[0]),
                    percent(c[1], t[1])
                )
            };
            [
                line(self.n_reads, "in total (QC-passed reads + QC-failed reads)"),
                line(self.n_primary, "primary"),
                line(self.n_secondary, "secondary"),
                line(self.n_supp, "supplementary"),
                line(self.n_dup, "duplicates"),
                line(self.n_pdup, "primary duplicates"),
                pct(self.n_mapped, self.n_reads, "mapped"),
                pct(self.n_pmapped, self.n_primary, "primary mapped"),
                line(self.n_pair_all, "paired in sequencing"),
                line(self.n_read1, "read1"),
                line(self.n_read2, "read2"),
                pct(self.n_pair_good, self.n_pair_all, "properly paired"),
                line(self.n_pair_map, "with itself and mate mapped"),
                pct(self.n_sgltn, self.n_pair_all, "singletons"),
                line(self.n_diffchr, "with mate mapped to a different chr"),
                line(
                    self.n_diffhigh,
                    "with mate mapped to a different chr (mapQ>=5)",
                ),
            ]
            .concat()
        }
    }

    /// (tid, mtid, mapq) placements: mate on the same contig, mate on another
    /// contig, mate on another contig with MAPQ below 5.
    const PLACEMENTS: [(i32, i32, u8); 3] = [(0, 0, 60), (0, 1, 60), (1, 0, 3)];

    fn record(flag: u16, (tid, mtid, mapq): (i32, i32, u8)) -> bam::Record {
        use rust_htslib::bam::record::{Cigar, CigarString};
        let mut rec = bam::Record::new();
        rec.set(
            b"r",
            Some(&CigarString(vec![Cigar::Match(4)])),
            b"ACGT",
            &[30; 4],
        );
        rec.set_flags(flag);
        rec.set_tid(tid);
        rec.set_pos(9);
        rec.set_mtid(mtid);
        rec.set_mpos(9);
        rec.set_mapq(mapq);
        rec
    }

    fn render_flagstat(accum: BamStatAccum, name: &str) -> String {
        let path = std::env::temp_dir().join(format!(
            "rustqc_test_flagstat_{name}_{}.txt",
            std::process::id()
        ));
        write_flagstat(&accum.into_result(), &path).expect("Failed to write flagstat");
        let contents = std::fs::read_to_string(&path).unwrap();
        let _ = std::fs::remove_file(&path);
        contents
    }

    /// Each of the 4,096 flags alone, in every placement, against the oracle.
    #[test]
    fn test_flagstat_matches_samtools_flagstat_loop_per_flag() {
        let mut mismatches = Vec::new();
        for placement in PLACEMENTS {
            for flag in 0..4096u16 {
                let mut accum = BamStatAccum::default();
                accum.process_read(&record(flag, placement), 30);
                let mut oracle = SamtoolsFlagstat::default();
                oracle.add(flag, placement.0, placement.1, placement.2);

                let got = render_flagstat(accum, "per_flag");
                let want = oracle.render();
                if got != want {
                    let diff: Vec<String> = want
                        .lines()
                        .zip(got.lines())
                        .filter(|(w, g)| w != g)
                        .map(|(w, g)| format!("want `{w}`, got `{g}`"))
                        .collect();
                    mismatches.push(format!("{flag:#06x} {placement:?}: {}", diff.join("; ")));
                }
            }
        }
        assert!(
            mismatches.is_empty(),
            "{} of {} (flag, placement) cases differ from samtools' flagstat_loop; first:\n{}",
            mismatches.len(),
            4096 * PLACEMENTS.len(),
            mismatches[..mismatches.len().min(8)].join("\n")
        );
    }

    /// All flags and placements in one file, split across two merged
    /// accumulators as the per-contig parallel path does.
    #[test]
    fn test_flagstat_matches_samtools_flagstat_loop_merged() {
        let mut parts = [BamStatAccum::default(), BamStatAccum::default()];
        let mut oracle = SamtoolsFlagstat::default();
        for placement in PLACEMENTS {
            for flag in 0..4096u16 {
                parts[placement.0 as usize].process_read(&record(flag, placement), 30);
                oracle.add(flag, placement.0, placement.1, placement.2);
            }
        }
        let [mut accum, other] = parts;
        accum.merge(other);
        assert_eq!(render_flagstat(accum, "merged"), oracle.render());
    }
}
