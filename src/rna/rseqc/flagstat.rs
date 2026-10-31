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
// Counting
// ============================================================================

/// samtools flagstat counters. Each is `[QC-passed, QC-failed]`, split on
/// 0x200 like samtools' `bam_flagstat_t`.
#[derive(Debug, Clone, Default)]
pub struct FlagstatCounts {
    /// All records.
    pub total: [u64; 2],
    /// Neither secondary nor supplementary.
    pub primary: [u64; 2],
    /// Secondary (0x100).
    pub secondary: [u64; 2],
    /// Supplementary (0x800), excluding secondary records.
    pub supplementary: [u64; 2],
    /// Duplicates (0x400).
    pub duplicates: [u64; 2],
    /// Primary duplicates.
    pub primary_duplicates: [u64; 2],
    /// Mapped (not 0x4).
    pub mapped: [u64; 2],
    /// Primary mapped.
    pub primary_mapped: [u64; 2],
    /// Primary and paired (0x1).
    pub paired: [u64; 2],
    /// Primary, paired, read1 (0x40).
    pub read1: [u64; 2],
    /// Primary, paired, read2 (0x80).
    pub read2: [u64; 2],
    /// Primary, paired, proper pair (0x2) and mapped.
    pub properly_paired: [u64; 2],
    /// Primary, paired, this and mate mapped.
    pub both_mapped: [u64; 2],
    /// Primary, paired, mapped with mate unmapped.
    pub singletons: [u64; 2],
    /// Both mapped, mate on a different reference.
    pub mate_diff_chr: [u64; 2],
    /// As `mate_diff_chr`, with MAPQ >= 5.
    pub mate_diff_chr_mapq5: [u64; 2],
}

impl FlagstatCounts {
    /// Count one record, following samtools' `flagstat_loop` (bam_stat.c).
    pub fn add(&mut self, flags: u16, tid: i32, mtid: i32, mapq: u8) {
        let w = usize::from(flags & BAM_FQCFAIL != 0);
        let is_secondary = flags & BAM_FSECONDARY != 0;
        let is_supplementary = flags & BAM_FSUPPLEMENTARY != 0;
        let is_mapped = flags & BAM_FUNMAP == 0;
        let is_dup = flags & BAM_FDUP != 0;

        self.total[w] += 1;
        if is_secondary {
            self.secondary[w] += 1;
        } else if is_supplementary {
            self.supplementary[w] += 1;
        }
        if is_mapped {
            self.mapped[w] += 1;
        }
        if is_dup {
            self.duplicates[w] += 1;
        }
        if is_secondary || is_supplementary {
            return;
        }

        self.primary[w] += 1;
        if is_mapped {
            self.primary_mapped[w] += 1;
        }
        if is_dup {
            self.primary_duplicates[w] += 1;
        }
        if flags & BAM_FPAIRED != 0 {
            self.paired[w] += 1;
            if flags & BAM_FREAD1 != 0 {
                self.read1[w] += 1;
            }
            if flags & BAM_FREAD2 != 0 {
                self.read2[w] += 1;
            }
            if is_mapped && flags & BAM_FPROPER_PAIR != 0 {
                self.properly_paired[w] += 1;
            }
            let mate_mapped = flags & BAM_FMUNMAP == 0;
            if is_mapped && mate_mapped {
                self.both_mapped[w] += 1;
                if tid != mtid {
                    self.mate_diff_chr[w] += 1;
                    if mapq >= 5 {
                        self.mate_diff_chr_mapq5[w] += 1;
                    }
                }
            }
            if is_mapped && !mate_mapped {
                self.singletons[w] += 1;
            }
        }
    }

    /// Add another set of counters into this one.
    pub fn merge(&mut self, other: &FlagstatCounts) {
        let pairs = [
            (&mut self.total, other.total),
            (&mut self.primary, other.primary),
            (&mut self.secondary, other.secondary),
            (&mut self.supplementary, other.supplementary),
            (&mut self.duplicates, other.duplicates),
            (&mut self.primary_duplicates, other.primary_duplicates),
            (&mut self.mapped, other.mapped),
            (&mut self.primary_mapped, other.primary_mapped),
            (&mut self.paired, other.paired),
            (&mut self.read1, other.read1),
            (&mut self.read2, other.read2),
            (&mut self.properly_paired, other.properly_paired),
            (&mut self.both_mapped, other.both_mapped),
            (&mut self.singletons, other.singletons),
            (&mut self.mate_diff_chr, other.mate_diff_chr),
            (&mut self.mate_diff_chr_mapq5, other.mate_diff_chr_mapq5),
        ];
        for (mine, theirs) in pairs {
            mine[0] += theirs[0];
            mine[1] += theirs[1];
        }
    }
}

// ============================================================================
// Output formatting
// ============================================================================

/// Format a flagstat count line as `<passed> + <failed> <description>`.
fn fmt_line(count: [u64; 2], description: &str) -> String {
    format!("{} + {} {description}", count[0], count[1])
}

/// Format a flagstat percentage as `(XX.XX% : YY.YY%)`, one per QC column,
/// `N/A` where the denominator is zero.
fn fmt_pct(numerator: [u64; 2], denominator: [u64; 2]) -> String {
    let pct = |n: u64, d: u64| {
        if d == 0 {
            "N/A".to_string()
        } else {
            // samtools bam_stat.c divides in float before multiplying by 100.0.
            format!("{:.2}%", (n as f32 / d as f32) as f64 * 100.0)
        }
    };
    format!(
        "({} : {})",
        pct(numerator[0], denominator[0]),
        pct(numerator[1], denominator[1])
    )
}

/// Write samtools flagstat-compatible output.
///
/// The output format matches `samtools flagstat` exactly, with one line per
/// statistic. Each line has the form `<passed> + <failed> <description>`.
///
/// # Arguments
/// * `result` - The computed BAM statistics
/// * `output_path` - Path to write the flagstat file
pub fn write_flagstat(result: &BamStatResult, output_path: &Path) -> Result<()> {
    use std::io::Write;

    let mut out = std::fs::File::create(output_path)
        .with_context(|| format!("Failed to create flagstat file: {}", output_path.display()))?;

    let fs = &result.flagstat;

    // Line 1: total (QC-passed + QC-failed)
    writeln!(
        out,
        "{}",
        fmt_line(fs.total, "in total (QC-passed reads + QC-failed reads)")
    )?;

    // Line 2: primary
    writeln!(out, "{}", fmt_line(fs.primary, "primary"))?;

    // Line 3: secondary
    writeln!(out, "{}", fmt_line(fs.secondary, "secondary"))?;

    // Line 4: supplementary
    writeln!(out, "{}", fmt_line(fs.supplementary, "supplementary"))?;

    // Line 5: duplicates
    writeln!(out, "{}", fmt_line(fs.duplicates, "duplicates"))?;

    // Line 6: primary duplicates
    writeln!(
        out,
        "{}",
        fmt_line(fs.primary_duplicates, "primary duplicates")
    )?;

    // Line 7: mapped
    writeln!(
        out,
        "{} {}",
        fmt_line(fs.mapped, "mapped"),
        fmt_pct(fs.mapped, fs.total)
    )?;

    // Line 8: primary mapped
    writeln!(
        out,
        "{} {}",
        fmt_line(fs.primary_mapped, "primary mapped"),
        fmt_pct(fs.primary_mapped, fs.primary)
    )?;

    // Line 9: paired in sequencing
    writeln!(out, "{}", fmt_line(fs.paired, "paired in sequencing"))?;

    // Line 10: read1
    writeln!(out, "{}", fmt_line(fs.read1, "read1"))?;

    // Line 11: read2
    writeln!(out, "{}", fmt_line(fs.read2, "read2"))?;

    // Line 12: properly paired (percentage of paired reads)
    writeln!(
        out,
        "{} {}",
        fmt_line(fs.properly_paired, "properly paired"),
        fmt_pct(fs.properly_paired, fs.paired)
    )?;

    // Line 13: with itself and mate mapped
    writeln!(
        out,
        "{}",
        fmt_line(fs.both_mapped, "with itself and mate mapped")
    )?;

    // Line 14: singletons (percentage of paired reads)
    writeln!(
        out,
        "{} {}",
        fmt_line(fs.singletons, "singletons"),
        fmt_pct(fs.singletons, fs.paired)
    )?;

    // Line 15: with mate mapped to a different chr
    writeln!(
        out,
        "{}",
        fmt_line(fs.mate_diff_chr, "with mate mapped to a different chr")
    )?;

    // Line 16: with mate mapped to a different chr (mapQ>=5)
    writeln!(
        out,
        "{}",
        fmt_line(
            fs.mate_diff_chr_mapq5,
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

    #[test]
    fn test_flagstat_dual_flags_with_qcfail_and_duplicate() {
        let mut accum = BamStatAccum::default();
        let mut direct = FlagstatCounts::default();
        let flags = [
            0,
            BAM_FDUP,
            BAM_FSECONDARY,
            BAM_FSUPPLEMENTARY,
            BAM_FSECONDARY | BAM_FSUPPLEMENTARY,
            BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FDUP,
        ];
        for qc in [0, BAM_FQCFAIL] {
            for flag in flags {
                let flag = flag | qc;
                accum.process_read(&record(flag, PLACEMENTS[0]), 30);
                direct.add(flag, 0, 0, 60);
            }
        }
        for fs in [&accum.flagstat, &direct] {
            assert_eq!(fs.total, [6, 6]);
            assert_eq!(fs.primary, [2, 2]);
            assert_eq!(fs.secondary, [3, 3]);
            assert_eq!(fs.supplementary, [1, 1]);
            assert_eq!(fs.duplicates, [2, 2]);
            assert_eq!(fs.primary_duplicates, [1, 1]);
            assert_eq!(fs.mapped, [6, 6]);
            assert_eq!(fs.primary_mapped, [2, 2]);
            for w in 0..2 {
                assert_eq!(
                    fs.primary[w] + fs.secondary[w] + fs.supplementary[w],
                    fs.total[w]
                );
            }
        }
        let result = accum.into_result();
        assert_eq!(result.total_records, 12);
        assert_eq!(result.primary_count, 4);
        assert_eq!(result.secondary, 6);
        assert_eq!(result.supplementary, 2);
        assert_eq!(
            result.primary_count + result.secondary + result.supplementary,
            result.total_records
        );
    }

    #[test]
    fn test_flagstat_partitioned_counts_match_serial() {
        for partitions in [1, 2, 7, 31] {
            let mut parts = vec![FlagstatCounts::default(); partitions];
            let mut oracle = SamtoolsFlagstat::default();
            for (i, (flag, placement)) in (0..4096u16)
                .flat_map(|flag| PLACEMENTS.map(|placement| (flag, placement)))
                .enumerate()
            {
                parts[i % partitions].add(flag, placement.0, placement.1, placement.2);
                oracle.add(flag, placement.0, placement.1, placement.2);
            }
            let mut merged = FlagstatCounts::default();
            for part in parts.into_iter().rev() {
                merged.merge(&part);
            }
            for w in 0..2 {
                assert_eq!(
                    merged.primary[w] + merged.secondary[w] + merged.supplementary[w],
                    merged.total[w]
                );
            }
            let mut accum = BamStatAccum::default();
            accum.flagstat.merge(&merged);
            assert_eq!(
                render_flagstat(accum, "partitioned"),
                oracle.render(),
                "{partitions} partitions"
            );
        }
    }

    #[test]
    fn test_flagstat_large_counter_percentages_match_samtools_c() {
        // samtools 1.22.1 and 1.24 bam_stat.c percent():
        // sprintf(buffer, "%.2f%%", (float)n / total * 100.0).
        // Expected strings were checked by compiling that C function unchanged.
        let cases = [
            (0, 0, "N/A"),
            (0, 1_000_000_000, "0.00%"),
            (1, 32, "3.12%"),
            (1, 6, "16.67%"),
            (16_777_217, 16_777_218, "100.00%"),
            (500_049_999, 1_000_000_000, "50.01%"),
            (999_950_000, 1_000_000_000, "99.99%"),
            (499_050_000, 1_000_000_000, "49.90%"),
            (499_150_000, 1_000_000_000, "49.92%"),
            (i64::MAX as u64 - 1, i64::MAX as u64, "100.00%"),
        ];
        for (n, d, expected) in cases {
            assert_eq!(
                fmt_pct([n, n], [d, d]),
                format!("({expected} : {expected})"),
                "{n}/{d}"
            );
        }
    }
}
