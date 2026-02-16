---
title: Introduction
description: What is RustQC and why use it for RNA-seq quality control.
---

RustQC is a fast quality control toolkit for sequencing data, written in Rust. It reimplements several established bioinformatics QC tools as a single compiled binary with no runtime dependencies:

- **[dupRadar](https://bioconductor.org/packages/dupRadar/)** -- RNA-seq PCR duplicate rate analysis
- **[featureCounts](http://subread.sourceforge.net/)** -- gene-level read counting and biotype quantification
- **[RSeQC](https://rseqc.sourceforge.net/)** -- 7 RNA-seq quality control modules (bam_stat, infer_experiment, read_duplication, read_distribution, junction_annotation, junction_saturation, inner_distance) plus TIN (Transcript Integrity Number) analysis
- **[preseq](https://github.com/smithlabcode/preseq)** -- library complexity extrapolation (lc_extrap)
- **[samtools](http://www.htslib.org/)** -- flagstat, idxstats, and stats SN-section compatible outputs
- **[Qualimap](http://qualimap.conesalab.org/)** -- gene body coverage profiling and RNA-seq QC summary

## Why RustQC?

RNA-seq quality control typically involves running multiple tools written in R and Python, each with their own dependencies, interpreters, and runtime overhead. RustQC consolidates these into a single fast binary.

<div class="benchmark-chart">
  <img class="only-dark" src="/benchmarks/benchmark_dark.svg" alt="Benchmark: RustQC ~5m vs traditional tools ~2h 45m" />
  <img class="only-light" src="/benchmarks/benchmark_light.svg" alt="Benchmark: RustQC ~5m vs traditional tools ~2h 45m" />
  <p align="center"><em>Run time for a 10 GB paired-end BAM (dupRadar + featureCounts + RSeQC + preseq + samtools + Qualimap)</em></p>
</div>

Key advantages:

- **Speed**: Up to 21x faster than the traditional R/Python workflow
- **Single binary**: No runtime dependencies -- no R, Python, or Bioconductor installation required
- **Single-pass architecture**: The `rna` subcommand performs read counting and duplicate analysis simultaneously, eliminating the need for a separate featureCounts run
- **Identical output**: Produces bit-for-bit identical results to the original tools (verified across 820,000+ values with zero mismatches for dupRadar)
- **Multiple input support**: Process several BAM files in a single command with automatic parallelisation
- **Modern format support**: Accepts SAM, BAM, and CRAM input files; annotation files (GTF and BED) can be plain or gzip-compressed

## Available tools

### `rustqc rna` -- RNA-Seq quality control pipeline

Given a duplicate-marked BAM file and a GTF or BED12 annotation, `rustqc rna` runs all of the following in a single pass:

1. **Counts reads** per gene (equivalent to Subread featureCounts)
2. **Computes duplication rates** per gene at multiple counting levels (dupRadar)
3. **Fits a logistic regression model** relating expression to duplication rate
4. **Generates diagnostic plots**: density scatter, boxplot, and expression histogram
5. **Produces MultiQC-compatible reports** for integration into analysis pipelines
6. **Outputs featureCounts-format files** including biotype-level summaries
7. **Estimates library complexity** via extrapolation (preseq lc_extrap equivalent)
8. **Computes Transcript Integrity Number (TIN)** for RNA degradation assessment
9. **Profiles gene body coverage** with Qualimap-compatible output
10. **Produces samtools-compatible outputs** (flagstat, idxstats, stats)

### RSeQC tools -- RNA-seq quality control

Seven reimplementations of [RSeQC](https://rseqc.sourceforge.net/) tools, plus TIN analysis, all integrated into the `rustqc rna` command and running automatically in the same single-pass analysis:

| Tool | Description |
|------|-------------|
| bam_stat | Basic BAM alignment statistics (total reads, duplicates, mapping quality, etc.) |
| infer_experiment | Infer library strandedness from read/gene-model strand concordance |
| read_duplication | Position-based and sequence-based read duplication histograms with plots |
| read_distribution | Classify reads across genomic features (CDS, UTR, intron, intergenic) |
| junction_annotation | Classify splice junctions as known, partial novel, or complete novel, with pie chart plots |
| junction_saturation | Assess saturation of splice junction detection at increasing read depths, with plot |
| inner_distance | Compute inner distance between paired-end read mates, with histogram plot |
| tin | Transcript Integrity Number (TIN) for RNA degradation assessment |

### Additional tools

| Tool | Equivalent | Description |
|------|-----------|-------------|
| preseq | [preseq](http://smithlabresearch.org/software/preseq/) `lc_extrap` | Library complexity extrapolation with bootstrap confidence intervals |
| gene body coverage | [Qualimap](http://qualimap.conesalab.org/) rnaseq | Coverage profile along gene bodies with 5'/3' bias metrics |
| flagstat | `samtools flagstat` | Alignment flag statistics |
| idxstats | `samtools idxstats` | Per-chromosome read counts |
| stats | `samtools stats` | Summary number (SN) statistics |

When a GTF file is provided via `--gtf`, all tools run automatically — transcript-level structure is extracted from the GTF. Alternatively, a BED12 gene model file can be provided via `--bed` (mutually exclusive with `--gtf`), which runs the RSeQC tools, TIN, preseq, and samtools outputs, but skips dupRadar, featureCounts, and gene body coverage (they require a GTF). Both GTF and BED files can be provided plain or gzip-compressed (`.gz`) — compression is detected automatically. Individual tools can be disabled via the YAML configuration file.

## Credits

RustQC stands on the shoulders of the original tools. If you use RustQC, please cite [dupRadar](https://bioconductor.org/packages/dupRadar/), [Subread/featureCounts](http://subread.sourceforge.net/), [RSeQC](https://rseqc.sourceforge.net/), and [preseq](http://smithlabresearch.org/software/preseq/). See [Credits & Citation](/about/credits/) for full details.
