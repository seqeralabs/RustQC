//! Emit raw production flagstat counters for independent differential validation.
//! Each row calls the library implementation on a fresh counter.

use rustqc::rna::rseqc::flagstat;

fn main() -> std::io::Result<()> {
    use std::io::Write;

    let mut out = std::io::BufWriter::new(std::io::stdout().lock());
    writeln!(out, "flags\ttid\tmtid\tmapq\ttotal\tprimary\tsecondary\tsupplementary\tduplicates\tprimary_duplicates\tmapped\tprimary_mapped\tpaired\tread1\tread2\tproperly_paired\tboth_mapped\tsingletons\tmate_diff_chr\tmate_diff_chr_mapq5")?;
    for flags in 0..4096u16 {
        for (tid, mtid, mapq) in [(0, 0, 60), (0, 1, 60), (1, 0, 4)] {
            let mut c = flagstat::FlagstatCounts::default();
            c.add(flags, tid, mtid, mapq);
            write!(out, "{flags}\t{tid}\t{mtid}\t{mapq}")?;
            for value in [
                c.total,
                c.primary,
                c.secondary,
                c.supplementary,
                c.duplicates,
                c.primary_duplicates,
                c.mapped,
                c.primary_mapped,
                c.paired,
                c.read1,
                c.read2,
                c.properly_paired,
                c.both_mapped,
                c.singletons,
                c.mate_diff_chr,
                c.mate_diff_chr_mapq5,
            ] {
                write!(out, "\t{},{}", value[0], value[1])?;
            }
            writeln!(out)?;
        }
    }
    out.flush()
}
