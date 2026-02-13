---
title: Introduction
description: What is RustQC and why use it for RNA-seq quality control.
---

RustQC is a fast quality control toolkit for sequencing data, written in Rust. It reimplements several established bioinformatics QC tools as a single compiled binary with no runtime dependencies:

- **[dupRadar](https://bioconductor.org/packages/dupRadar/)** -- RNA-seq PCR duplicate rate analysis
- **[featureCounts](http://subread.sourceforge.net/)** -- gene-level read counting and biotype quantification
- **[RSeQC](https://rseqc.sourceforge.net/)** -- 7 RNA-seq quality control modules (bam_stat, infer_experiment, read_duplication, read_distribution, junction_annotation, junction_saturation, inner_distance)

## Why RustQC?

RNA-seq quality control typically involves running multiple tools written in R and Python, each with their own dependencies, interpreters, and runtime overhead. RustQC consolidates these into a single fast binary.

<div class="benchmark-chart">
  <img class="only-dark" src="/benchmarks/benchmark_dark.svg" alt="Benchmark: RustQC 54s vs featureCounts 16m vs dupRadar 30m" />
  <img class="only-light" src="/benchmarks/benchmark_light.svg" alt="Benchmark: RustQC 54s vs featureCounts 16m vs dupRadar 30m" />
  <p align="center"><em>Run time for dupRadar + featureCounts analysis on a 10 GB paired-end BAM</em></p>
</div>

Key advantages:

- **Speed**: Up to 33x faster than the R/Python implementations
- **Single binary**: No runtime dependencies -- no R, Python, or Bioconductor installation required
- **Single-pass architecture**: The `rna` subcommand performs read counting and duplicate analysis simultaneously, eliminating the need for a separate featureCounts run
- **Identical output**: Produces bit-for-bit identical results to the original tools (verified across 820,000+ values with zero mismatches for dupRadar)
- **Multiple input support**: Process several BAM files in a single command with automatic parallelisation
- **Modern format support**: Accepts SAM, BAM, and CRAM input files

## Available tools

### `rustqc rna` -- Duplicate rate analysis + read counting

Given a duplicate-marked BAM file and a GTF annotation:

1. **Counts reads** per gene (equivalent to Subread featureCounts)
2. **Computes duplication rates** per gene at multiple counting levels
3. **Fits a logistic regression model** relating expression to duplication rate
4. **Generates diagnostic plots**: density scatter, boxplot, and expression histogram
5. **Produces MultiQC-compatible reports** for integration into analysis pipelines
6. **Outputs featureCounts-format files** including biotype-level summaries

### RSeQC tools -- RNA-seq quality control

Seven reimplementations of [RSeQC](https://rseqc.sourceforge.net/) tools, each available as a top-level subcommand:

| Subcommand | Description |
|------------|-------------|
| `rustqc bam-stat` | Basic BAM alignment statistics (total reads, duplicates, mapping quality, etc.) |
| `rustqc infer-experiment` | Infer library strandedness from read/gene-model strand concordance |
| `rustqc read-duplication` | Position-based and sequence-based read duplication histograms |
| `rustqc read-distribution` | Classify reads across genomic features (CDS, UTR, intron, intergenic) |
| `rustqc junction-annotation` | Classify splice junctions as known, partial novel, or complete novel |
| `rustqc junction-saturation` | Assess saturation of splice junction detection at increasing read depths |
| `rustqc inner-distance` | Compute inner distance between paired-end read mates |

Most RSeQC tools require a BED12 gene model file (`-b`) instead of a GTF. All accept SAM/BAM/CRAM input and support multiple input files.

## Credits

RustQC stands on the shoulders of the original tools. If you use RustQC, please cite [dupRadar](https://bioconductor.org/packages/dupRadar/), [Subread/featureCounts](http://subread.sourceforge.net/), and [RSeQC](https://rseqc.sourceforge.net/). See [Credits & Citation](/about/credits/) for full details.
