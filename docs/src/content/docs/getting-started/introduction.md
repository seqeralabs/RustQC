---
title: Introduction
description: What is RustQC and why use it for RNA-seq quality control.
---

RustQC is a fast quality control toolkit for sequencing data, written in Rust. It currently includes an RNA-seq duplicate rate analyser that reimplements [dupRadar](https://bioconductor.org/packages/dupRadar/) (a Bioconductor R package) with built-in featureCounts-compatible read counting.

## Why RustQC?

RNA-seq experiments rely on PCR amplification, which introduces duplicate reads. Distinguishing technical duplicates (PCR artifacts) from natural duplicates (multiple reads from highly expressed genes) is critical for data quality assessment. The original [dupRadar](https://bioconductor.org/packages/dupRadar/) tool models duplication rate as a function of gene expression, providing a clear diagnostic of library quality.

RustQC reimplements this analysis with several advantages:

- **Speed**: 9--33x faster than the R implementation, depending on thread count
- **Single-pass architecture**: Performs read counting and duplicate analysis simultaneously, eliminating the need for a separate featureCounts run
- **Identical output**: Produces bit-for-bit identical results to R dupRadar (verified across 820,000+ values with zero mismatches)
- **Multiple input support**: Process several BAM files in a single command with automatic parallelisation
- **Modern format support**: Accepts SAM, BAM, and CRAM input files

## What does it do?

Given a duplicate-marked BAM file and a GTF annotation, RustQC:

1. **Counts reads** per gene (equivalent to Subread featureCounts)
2. **Computes duplication rates** per gene at multiple counting levels
3. **Fits a logistic regression model** relating expression to duplication rate
4. **Generates diagnostic plots**: density scatter, boxplot, and expression histogram
5. **Produces MultiQC-compatible reports** for integration into analysis pipelines
6. **Outputs featureCounts-format files** including biotype-level summaries

## Credits

RustQC stands on the shoulders of the original tools. If you use RustQC, please cite both [dupRadar](https://bioconductor.org/packages/dupRadar/) and [Subread/featureCounts](http://subread.sourceforge.net/). See [Credits & Citation](/RustQC/about/credits/) for full details.
