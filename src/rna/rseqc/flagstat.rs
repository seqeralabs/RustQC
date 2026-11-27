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
    use crate::rna::bam::{self, Read as BamRead};
    use crate::rna::rseqc::accumulators::BamStatAccum;
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

    // ------------------------------------------------------------------------
    // Oracle: samtools 1.22.1 `flagstat_loop` and `out_fmt_default`
    // (bam_stat.c), transcribed independently of `FlagstatCounts`.
    // Test cases adapted from seqeralabs/RustQC#175 (Nico DeVeaux).
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
                    format!("{:.2}%", f64::from(n as f32 / total as f32) * 100.0)
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
        use crate::rna::bam::header::{HeaderRecord, HeaderView};
        let mut header = bam::Header::new();
        for name in ["chr1", "chr2"] {
            header.push_record(
                HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", name)
                    .push_tag(b"LN", "1000"),
            );
        }
        let chrom = |tid: i32| ["chr1", "chr2"][tid as usize];
        let line = format!(
            "r\t{flag}\t{}\t10\t{mapq}\t4M\t{}\t10\t0\tACGT\t????",
            chrom(tid),
            chrom(mtid)
        );
        bam::Record::from_sam(&HeaderView::from_header(&header), line.as_bytes()).unwrap()
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
    fn flagstat_matches_samtools_flagstat_loop_per_flag() {
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
    fn flagstat_matches_samtools_flagstat_loop_merged() {
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

    /// Dual-flagged records (secondary + supplementary), QC-failed and
    /// duplicated: counted once, as secondary, in the right column, in both the
    /// flagstat counters and the legacy per-field counters used by `stats`.
    #[test]
    fn flagstat_dual_flags_with_qcfail_and_duplicate() {
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
            assert_eq!(fs.n_reads, [6, 6]);
            assert_eq!(fs.n_primary, [2, 2]);
            assert_eq!(fs.n_secondary, [3, 3]);
            assert_eq!(fs.n_supp, [1, 1]);
            assert_eq!(fs.n_dup, [2, 2]);
            assert_eq!(fs.n_pdup, [1, 1]);
            assert_eq!(fs.n_mapped, [6, 6]);
            assert_eq!(fs.n_pmapped, [2, 2]);
            for w in 0..2 {
                assert_eq!(
                    fs.n_primary[w] + fs.n_secondary[w] + fs.n_supp[w],
                    fs.n_reads[w]
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
    fn flagstat_partitioned_counts_match_serial() {
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
                    merged.n_primary[w] + merged.n_secondary[w] + merged.n_supp[w],
                    merged.n_reads[w]
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
    fn large_counter_percentages_match_samtools_c() {
        // samtools 1.22.1 bam_stat.c percent():
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
            assert_eq!(percent(n, d), expected, "{n}/{d}");
        }
    }
}
