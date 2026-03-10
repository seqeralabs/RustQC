---
title: Introduction
description: What is RustQC and why use it for RNA-seq quality control.
---

RustQC is a fast quality control toolkit for sequencing data, written in Rust. It reimplements several established bioinformatics QC tools as a single compiled binary with no runtime dependencies:

- **[dupRadar](https://bioconductor.org/packages/dupRadar/)** -- RNA-seq PCR duplicate rate analysis
- **[featureCounts](http://subread.sourceforge.net/)** -- gene-level read counting and biotype quantification
- **[RSeQC](https://rseqc.sourceforge.net/)** -- 8 RNA-seq quality control modules (bam_stat, infer_experiment, read_duplication, read_distribution, junction_annotation, junction_saturation, inner_distance, TIN)
- **[preseq](https://github.com/smithlabcode/preseq)** -- library complexity extrapolation (lc_extrap)
- **[samtools](http://www.htslib.org/)** -- flagstat, idxstats, and full stats output including all histogram sections
- **[Qualimap](http://qualimap.conesalab.org/)** -- gene body coverage profiling and RNA-seq QC summary

## Why RustQC?

RNA-seq quality control typically involves running multiple tools written in R and Python, each with their own dependencies, interpreters, and runtime overhead. RustQC consolidates these into a single fast binary.

<div class="benchmark-chart">
  <img class="only-dark" src="/benchmarks/benchmark_dark.png" alt="Benchmark: RustQC ~41m vs traditional tools ~15h 34m (sequential)" />
  <img class="only-light" src="/benchmarks/benchmark_light.png" alt="Benchmark: RustQC ~41m vs traditional tools ~15h 34m (sequential)" />
  <p align="center"><em>Run time for a large paired-end RNA-seq BAM (~186M reads) on AWS. RSeQC TIN alone takes 9h 45m; RustQC completes everything in 41 minutes.</em></p>
</div>

Key advantages:

- **Speed**: Processes ~186M reads in ~41 minutes vs. ~15h 34m of sequential tool runtimes on AWS — including TIN, which takes 9h 45m alone in the traditional workflow
- **Single binary**: No runtime dependencies -- no R, Python, or Bioconductor installation required
- **Single-pass architecture**: The `rna` subcommand performs read counting and duplicate analysis simultaneously, eliminating the need for a separate featureCounts run
- **Identical output**: Produces bit-for-bit identical results to the original tools (verified across 820,000+ values with zero mismatches for dupRadar)
- **Multiple input support**: Process several BAM files in a single command with automatic parallelisation
- **Modern format support**: Accepts SAM, BAM, and CRAM input files; GTF annotation files can be plain or gzip-compressed

## Available tools

### `rustqc rna` -- RNA-Seq quality control pipeline

Given a duplicate-marked BAM file and a GTF annotation, `rustqc rna` runs all of the following in a single pass:

1. **Counts reads** per gene (equivalent to Subread featureCounts)
2. **Computes duplication rates** per gene at multiple counting levels (dupRadar)
3. **Fits a logistic regression model** relating expression to duplication rate
4. **Generates diagnostic plots**: density scatter, boxplot, and expression histogram
5. **Produces MultiQC-compatible reports** for integration into analysis pipelines
6. **Outputs featureCounts-format files** including biotype-level summaries
7. **Estimates library complexity** via extrapolation (preseq lc_extrap equivalent)
8. **Computes Transcript Integrity Number (TIN)** for RNA degradation assessment
9. **Runs Qualimap rnaseq analysis**: gene body coverage profiling, 5'/3' bias metrics, read origin classification, strand-specificity estimation, and junction analysis -- all with Qualimap-compatible MultiQC-parseable output
10. **Produces samtools-compatible outputs** (flagstat, idxstats, stats)

### RSeQC tools -- RNA-seq quality control

Eight reimplementations of [RSeQC](https://rseqc.sourceforge.net/) tools (including TIN), all integrated into the `rustqc rna` command and running automatically in the same single-pass analysis:

| Tool                | Description                                                                                |
| ------------------- | ------------------------------------------------------------------------------------------ |
| bam_stat            | Basic BAM alignment statistics (total reads, duplicates, mapping quality, etc.)            |
| infer_experiment    | Infer library strandedness from read/gene-model strand concordance                         |
| read_duplication    | Position-based and sequence-based read duplication histograms with plots                   |
| read_distribution   | Classify reads across genomic features (CDS, UTR, intron, intergenic)                      |
| junction_annotation | Classify splice junctions as known, partial novel, or complete novel, with pie chart plots |
| junction_saturation | Assess saturation of splice junction detection at increasing read depths, with plot        |
| inner_distance      | Compute inner distance between paired-end read mates, with histogram plot                  |
| tin                 | Transcript Integrity Number (TIN) for RNA degradation assessment                           |

### Additional tools

| Tool            | Equivalent                                                         | Description                                                                                                                                                                                                       |
| --------------- | ------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Qualimap rnaseq | [Qualimap](http://qualimap.conesalab.org/) rnaseq                  | Gene body coverage profiling, 5'/3' bias metrics, read origin classification (exonic/intronic/intergenic), strand-specificity estimation, and junction analysis. Qualimap-compatible output parseable by MultiQC. |
| preseq          | [preseq](http://smithlabresearch.org/software/preseq/) `lc_extrap` | Library complexity extrapolation with bootstrap confidence intervals                                                                                                                                              |
| flagstat        | `samtools flagstat`                                                | Alignment flag statistics                                                                                                                                                                                         |
| idxstats        | `samtools idxstats`                                                | Per-chromosome read counts                                                                                                                                                                                        |
| stats           | `samtools stats`                                                   | Full samtools stats output including all histogram sections                                                                                                                                                       |

A GTF file is required (`--gtf`) and all tools run automatically — transcript-level structure is extracted from the GTF. GTF files can be provided plain or gzip-compressed (`.gz`) — compression is detected automatically. Individual tools can be disabled via the YAML configuration file.

## Credits

RustQC stands on the shoulders of the original tools. If you use RustQC, please cite [dupRadar](https://bioconductor.org/packages/dupRadar/), [Subread/featureCounts](http://subread.sourceforge.net/), [RSeQC](https://rseqc.sourceforge.net/), and [preseq](http://smithlabresearch.org/software/preseq/). See [Credits & Citation](/about/credits/) for full details.
