---
title: Quick Start
description: Get up and running with RustQC in minutes.
---

This guide walks you through a basic RustQC analysis from start to finish.

## RNA-seq duplicate analysis

Run all RNA-seq QC analyses (dupRadar, featureCounts, RSeQC tools including TIN, Qualimap, preseq, and samtools) in a single pass:

```bash
rustqc rna sample.markdup.bam --gtf genes.gtf -p -o results/
```

This command:

- Analyses `sample.markdup.bam` against `genes.gtf`
- Uses paired-end mode (`-p`)
- Writes all output to `results/`

### Input requirements

- **BAM file(s)** — must have PCR duplicates **marked** (SAM flag `0x400`) but **not removed**. RustQC needs both duplicate and non-duplicate reads to calculate duplication rates.
- **BAM index** (`.bai` / `.csi`) - not technically required, but highly recommended as it's needed for multi-threaded performance; without one, only a single counting thread is used. Should sit alongisde the BAM file with the same filename, plus the `.bai` suffix.
- **GTF annotation file** (`--gtf`) — provides all annotation needed for every tool. Gene models are derived internally so all analyses run. Can be plain or gzip-compressed (`.gz`) — compression is detected automatically.

Compatible duplicate-marking tools:

- [Picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)
- [samblaster](https://github.com/GregoryFaust/samblaster)
- [sambamba markdup](https://lomereiter.github.io/sambamba/)
- [biobambam bammarkduplicates](https://gitlab.com/german.tischler/biobambam2)

RustQC automatically checks the BAM header for duplicate-marking tool signatures and exits with an error if none are found. Use `--skip-dup-check` to bypass this validation if your duplicate-marking tool is not recognized.

Note that the RSeQC tools themselves do not require duplicate marking — the duplicate-marking requirement applies to the dupRadar analysis.


### Output

Output files are organized into subdirectories by tool group.
Files are generally named the same as their upstream tool equivalents.

This means that MultiQC should find them and report them as if they were created
by the original tool.

One addition is that RustQC produces SVG versions of all plots.

Use `--flat-output` to write all files directly to the output directory without subdirectories.

## Multiple BAM files

Multiple input files are accepted and processed sequentially, with each
producing its own set of output files:

```bash
rustqc rna sample1.bam sample2.bam sample3.bam \
  --gtf genes.gtf -p -t 8 -o results/
```

The annotation file is parsed once and shared across all samples. Threads
are distributed automatically among concurrent jobs.

## Next steps

- [CLI Reference](/usage/cli-reference/) for all available options
- [dupRadar Outputs](/outputs/dupradar/) for duplication analysis output details
- [featureCounts Outputs](/outputs/featurecounts/) for read counting output details
- [RSeQC Outputs](/outputs/rseqc/) for RSeQC tool output details
- [Configuration](/usage/configuration/) for YAML config options
