---
title: Quick Start
description: Get up and running with RustQC in minutes.
---

A basic RustQC analysis from install to results.

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

- BAM file(s): coordinate-sorted and duplicate-marked (not removed). See [CLI reference](/usage/cli-reference/) for supported duplicate-marking tools. The RSeQC tools do not require duplicate marking; this is only needed for the dupRadar analysis. Use `--skip-dup-check` to bypass the check.
- BAM index: a `.bai` index is not strictly required, but without one RustQC falls back to a single counting thread. For multi-threaded performance, ensure an index file is present alongside each BAM.
- GTF annotation file (`--gtf`): a gene annotation file (plain or gzip-compressed). See [CLI reference](/usage/cli-reference/) for details.

### Output

Output files are organized into subdirectories by tool group.
Files are generally named the same as their upstream tool equivalents.
This means that MultiQC should find them and report them as if they were created
by the original tool.
One addition is that RustQC produces SVG versions of all plots.

Use `--flat-output` to write all files directly to the output directory without subdirectories.

## Multiple BAM files

Multiple input files are accepted and processed in parallel, with each
producing its own set of output files:

```bash
rustqc rna sample1.bam sample2.bam sample3.bam \
  --gtf genes.gtf -p -t 8 -o results/
```

The annotation file is parsed once and shared across all samples. Threads
are distributed automatically among concurrent jobs.
