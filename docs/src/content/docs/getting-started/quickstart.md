---
title: Quick Start
description: Get up and running with RustQC in minutes.
---

This guide walks you through a basic RustQC analysis from start to finish.

## Prerequisites

Before running RustQC, you need:

1. A **duplicate-marked** alignment file (BAM, SAM, or CRAM). Tools like [Picard MarkDuplicates](https://broadinstitute.github.io/picard/), [samblaster](https://github.com/GregoryFaust/samblaster), or [sambamba](https://github.com/biod/sambamba) can mark duplicates.
2. A **GTF annotation** file for your reference genome (e.g., from [Ensembl](https://www.ensembl.org/) or [GENCODE](https://www.gencodegenes.org/)).

## Basic usage

Run the RNA-seq duplicate rate analysis:

```bash
rustqc rna sample.markdup.bam --gtf genes.gtf -p -o results/
```

This command:
- Analyses `sample.markdup.bam` against `genes.gtf`
- Uses paired-end mode (`-p`)
- Writes all output to `results/`

## Output

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

Plus SVG versions of all plots and MultiQC-compatible report files.

## Multiple BAM files

Process several samples at once:

```bash
rustqc rna sample1.bam sample2.bam sample3.bam --gtf genes.gtf -p -t 8 -o results/
```

The GTF is parsed once and shared across all samples. Threads are distributed automatically.

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
```

## Next steps

- [CLI Reference](/RustQC/usage/cli-reference/) for all available options
- [Output Files](/RustQC/usage/output-files/) for detailed output descriptions
- [Interpreting Plots](/RustQC/guide/interpreting-plots/) to understand your results
- [RNA-seq Duplicate Analysis](/RustQC/guide/rna-duprate/) for the underlying methodology
