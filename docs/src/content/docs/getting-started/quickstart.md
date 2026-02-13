---
title: Quick Start
description: Get up and running with RustQC in minutes.
---

This guide walks you through a basic RustQC analysis from start to finish.

## Prerequisites

RustQC includes several tools, each with different input requirements:

- **`rna` subcommand** (dupRadar + featureCounts): Requires a **duplicate-marked** alignment file (BAM, SAM, or CRAM) and a **GTF annotation** file. Duplicates must be flagged with SAM flag 0x400 by a tool like [Picard MarkDuplicates](https://broadinstitute.github.io/picard/), [samblaster](https://github.com/GregoryFaust/samblaster), or [sambamba](https://github.com/biod/sambamba).
- **RSeQC tools** (`bam-stat`, `infer-experiment`, etc.): Require an alignment file. Most also require a **BED12 gene model** file (except `bam-stat` and `read-duplication`).

## RNA-seq duplicate analysis

Run the dupRadar + featureCounts analysis:

```bash
rustqc rna sample.markdup.bam --gtf genes.gtf -p -o results/
```

This command:
- Analyses `sample.markdup.bam` against `genes.gtf`
- Uses paired-end mode (`-p`)
- Writes all output to `results/`

### Output

After running, you will find in the output directory:

| File | Description |
|------|-------------|
| `sample.markdup_dupMatrix.txt` | Full duplication matrix (TSV) |
| `sample.markdup_intercept_slope.txt` | Fitted model parameters |
| `sample.markdup_duprateExpDens.png` | Density scatter plot |
| `sample.markdup_duprateExpBoxplot.png` | Duplication rate boxplot |
| `sample.markdup_expressionHist.png` | Expression histogram |
| `sample.markdup.featureCounts.tsv` | Gene-level read counts |
| `sample.markdup.featureCounts.tsv.summary` | Counting summary statistics |

Plus SVG versions of all plots, biotype count tables, and MultiQC-compatible report files.
See [dupRadar Outputs](/outputs/dupradar/) and [featureCounts Outputs](/outputs/featurecounts/) for full details.

## RSeQC quality control tools

RustQC reimplements seven tools from [RSeQC](https://rseqc.sourceforge.net/).
Each is a standalone subcommand:

```bash
# Basic alignment statistics
rustqc bam-stat sample.bam -o results/

# Infer library strandedness
rustqc infer-experiment sample.bam -b genes.bed -o results/

# Read duplication histograms
rustqc read-duplication sample.bam -o results/

# Read distribution across genomic features
rustqc read-distribution sample.bam -b genes.bed -o results/

# Splice junction annotation
rustqc junction-annotation sample.bam -b genes.bed -o results/

# Junction saturation analysis
rustqc junction-saturation sample.bam -b genes.bed -o results/

# Inner distance between paired-end reads
rustqc inner-distance sample.bam -b genes.bed -o results/
```

See [RSeQC Outputs](/outputs/rseqc/) for details on output files.

## Multiple BAM files

All subcommands accept multiple input files. They are processed sequentially,
with each producing its own output files:

```bash
# RNA-seq analysis with parallel processing
rustqc rna sample1.bam sample2.bam sample3.bam --gtf genes.gtf -p -t 8 -o results/

# RSeQC tools with multiple inputs
rustqc bam-stat sample1.bam sample2.bam sample3.bam -o results/
```

For the `rna` subcommand, the GTF is parsed once and shared across all samples.
Threads are distributed automatically among concurrent jobs.

## Common options

```bash
# Stranded library (forward)
rustqc rna sample.bam --gtf genes.gtf -p -s 1

# Use 8 threads
rustqc rna sample.bam --gtf genes.gtf -p -t 8

# CRAM input with reference
rustqc rna sample.cram --gtf genes.gtf -p --reference genome.fa

# Skip duplicate-marking validation
rustqc rna sample.bam --gtf genes.gtf -p --skip-dup-check

# Use a YAML config file
rustqc rna sample.bam --gtf genes.gtf -p --config config.yaml

# Custom MAPQ cutoff for RSeQC tools
rustqc bam-stat sample.bam -q 20 -o results/
```

## Next steps

- [CLI Reference](/usage/cli-reference/) for all available options and subcommands
- [dupRadar Outputs](/outputs/dupradar/) for duplication analysis output details
- [featureCounts Outputs](/outputs/featurecounts/) for read counting output details
- [RSeQC Outputs](/outputs/rseqc/) for RSeQC tool output details
- [Configuration](/usage/configuration/) for YAML config options
