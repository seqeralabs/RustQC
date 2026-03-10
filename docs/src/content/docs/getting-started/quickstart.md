---
title: Quick Start
description: Get up and running with RustQC in minutes.
---

This guide walks you through a basic RustQC analysis from start to finish.

## Prerequisites

The `rustqc rna` command runs all analyses in a single pass. It requires:

- A **coordinate-sorted, duplicate-marked** alignment file (BAM, SAM, or CRAM). Duplicates must be flagged with SAM flag 0x400 by a tool like [Picard MarkDuplicates](https://broadinstitute.github.io/picard/), [samblaster](https://github.com/GregoryFaust/samblaster), or [sambamba](https://github.com/biod/sambamba). A BAM index (`.bai` / `.csi`) is recommended for multi-threaded performance; without one, only a single counting thread is used.
- A **GTF annotation** file (`--gtf`). Can be plain or gzip-compressed (`.gz`) — compression is detected automatically. The GTF provides all annotation needed for every analysis: dupRadar, featureCounts, all 8 RSeQC tools including TIN, Qualimap, preseq, and samtools.

## RNA-seq duplicate analysis

Run all RNA-seq QC analyses (dupRadar, featureCounts, RSeQC tools including TIN, Qualimap, preseq, and samtools) in a single pass:

```bash
rustqc rna sample.markdup.bam --gtf genes.gtf -p -o results/
```

This command:

- Analyses `sample.markdup.bam` against `genes.gtf`
- Uses paired-end mode (`-p`)
- Writes all output to `results/`

### Output

Output files are organized into subdirectories by tool group. After running, you
will find in the output directory:

| File                                                                        | Description                                                                     |
| --------------------------------------------------------------------------- | ------------------------------------------------------------------------------- |
| `dupradar/sample.markdup_dupMatrix.txt`                                     | Full duplication matrix (TSV)                                                   |
| `dupradar/sample.markdup_intercept_slope.txt`                               | Fitted model parameters                                                         |
| `dupradar/sample.markdup_duprateExpDens.png`                                | Density scatter plot                                                            |
| `dupradar/sample.markdup_duprateExpBoxplot.png`                             | Duplication rate boxplot                                                        |
| `dupradar/sample.markdup_expressionHist.png`                                | Expression histogram                                                            |
| `featurecounts/sample.markdup.featureCounts.tsv`                            | Gene-level read counts                                                          |
| `featurecounts/sample.markdup.featureCounts.tsv.summary`                    | Counting summary statistics                                                     |
| `rseqc/read_duplication/sample.markdup.DupRate_plot.png`                    | Read duplication rate plot                                                      |
| `rseqc/junction_annotation/sample.markdup.splice_events.png`                | Splice events pie chart                                                         |
| `rseqc/junction_annotation/sample.markdup.splice_junction.png`              | Splice junctions pie chart                                                      |
| `rseqc/junction_saturation/sample.markdup.junctionSaturation_plot.png`      | Junction saturation plot                                                        |
| `rseqc/inner_distance/sample.markdup.inner_distance_plot.png`               | Inner distance histogram                                                        |
| `rseqc/tin/sample.markdup.tin.xls`                                          | Transcript Integrity Numbers                                                    |
| `samtools/sample.markdup.flagstat`                                          | samtools flagstat-compatible output                                             |
| `samtools/sample.markdup.idxstats`                                          | samtools idxstats-compatible output                                             |
| `samtools/sample.markdup.stats`                                             | samtools stats compatible output (full format including all histogram sections) |
| `preseq/sample.markdup.lc_extrap.txt`                                       | Library complexity extrapolation                                                |
| `qualimap/rnaseq_qc_results.txt`                                            | Qualimap-compatible QC report (bias, read origin, SSP)                          |
| `qualimap/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt` | Gene body coverage profile                                                      |

Plus SVG versions of all plots, biotype count tables, and MultiQC-compatible
report files. Use `--flat-output` to write all files directly to the output
directory without subdirectories.

See [dupRadar Outputs](/outputs/dupradar/), [featureCounts Outputs](/outputs/featurecounts/), [RSeQC Outputs](/outputs/rseqc/), [TIN](/outputs/tin/), [Qualimap](/outputs/qualimap/), [Preseq Outputs](/outputs/preseq/), and [Samtools Outputs](/outputs/samtools/) for full details.

## RSeQC quality control tools

RustQC reimplements eight [RSeQC](https://rseqc.sourceforge.net/) tools (including TIN), all
integrated into the `rustqc rna` command. They run automatically alongside
dupRadar, featureCounts, Qualimap, preseq, and samtools analyses:

```bash
rustqc rna sample.markdup.bam --gtf genes.gtf -p -o results/
```

To disable specific RSeQC tools, use a YAML config file:

```yaml
# config.yaml -- disable inner_distance
inner_distance:
  enabled: false
```

```bash
rustqc rna sample.markdup.bam --gtf genes.gtf -p -c config.yaml -o results/
```

See [RSeQC Outputs](/outputs/rseqc/) for details on output files.

## Multiple BAM files

Multiple input files are accepted and processed sequentially, with each
producing its own set of output files:

```bash
rustqc rna sample1.bam sample2.bam sample3.bam \
  --gtf genes.gtf -p -t 8 -o results/
```

The annotation file is parsed once and shared across all samples. Threads
are distributed automatically among concurrent jobs.

## Common options

```bash
# Stranded library (forward)
rustqc rna sample.bam --gtf genes.gtf -p -s 1

# Use 8 threads
rustqc rna sample.bam --gtf genes.gtf -p -t 8

# CRAM input with reference
rustqc rna sample.cram --gtf genes.gtf -p --reference genome.fa

# Gzip-compressed annotation file (auto-detected)
rustqc rna sample.bam --gtf genes.gtf.gz -p -o results/

# Skip duplicate-marking validation
rustqc rna sample.bam --gtf genes.gtf -p --skip-dup-check

# Use a YAML config file
rustqc rna sample.bam --gtf genes.gtf -p --config config.yaml

# Custom MAPQ cutoff for RSeQC tools
rustqc rna sample.bam --gtf genes.gtf -p -q 20 -o results/
```

## Next steps

- [CLI Reference](/usage/cli-reference/) for all available options
- [dupRadar Outputs](/outputs/dupradar/) for duplication analysis output details
- [featureCounts Outputs](/outputs/featurecounts/) for read counting output details
- [RSeQC Outputs](/outputs/rseqc/) for RSeQC tool output details
- [Configuration](/usage/configuration/) for YAML config options
