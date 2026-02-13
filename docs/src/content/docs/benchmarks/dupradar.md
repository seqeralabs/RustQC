---
title: dupRadar Benchmark
description: Performance comparison between RustQC and R dupRadar, demonstrating up to 33x speedup with zero output mismatches across all values.
---

RustQC produces output **identical** to the R dupRadar package. Every value in the duplication matrix, every count, and the fitted model parameters match exactly. The benchmarks below quantify the performance difference on real RNA-seq data.

## Small benchmark

**Input:** A test BAM file with a chr6-only GTF annotation (2,905 genes).

| Metric | dupRadar (R) | RustQC |
|--------|-------------|--------|
| **Runtime** | 2.50s | 0.25s |
| **Speedup** | -- | **10x** |

### Results comparison (small)

| Metric | dupRadar (R) | RustQC | Exact match |
|--------|------------:|-------:|:-----------:|
| **Intercept** | 0.03186 | 0.03186 | Yes |
| **Slope** | 1.60189 | 1.60189 | Yes |
| **Genes total** | 2,905 | 2,905 | Yes |
| **Genes with reads** | 636 | 636 | Yes |
| **Genes with duplicates** | 201 | 201 | Yes |
| **allCounts (unique)** | 20,449 | 20,449 | 100% |
| **filteredCounts (unique)** | 17,879 | 17,879 | 100% |
| **allCountsMulti** | 22,812 | 22,812 | 100% |
| **filteredCountsMulti** | 20,034 | 20,034 | 100% |
| **Total values compared** | 37,765 | 37,765 | -- |
| **Value mismatches** | -- | **0** | -- |

## Large benchmark

**Input:** GM12878 REP1 -- a full-size RNA-seq BAM (~10 GB) from the [nf-core/rnaseq](https://nf-co.re/rnaseq) pipeline, duplicate-marked with Picard. Paired-end, unstranded, aligned to GRCh38 (63,086 genes).

### Performance

| Tool | Runtime | Max RSS |
|------|---------|---------|
| R dupRadar | 29m 56s | N/A (Docker) |
| RustQC (1 thread) | 3m 16s | 503 MB |
| RustQC (8 threads) | 1m 03s | 893 MB |
| RustQC (10 threads) | 0m 54s | 1.3 GB |

Speedup: **~9x** with 1 thread, **~28x** with 8 threads, **~33x** with 10 threads.

### Results comparison (large)

| Metric | dupRadar (R) | RustQC | Exact match |
|--------|------------:|-------:|:-----------:|
| **Intercept** | 0.8245 | 0.8245 | Yes |
| **Slope** | 1.6774 | 1.6774 | Yes |
| **Genes total** | 63,086 | 63,086 | Yes |
| **Genes with reads (unique)** | 23,597 | 23,597 | Yes |
| **Genes with reads (multi)** | 24,719 | 24,719 | Yes |
| **allCounts (unique)** | 14,654,579 | 14,654,579 | 100% |
| **filteredCounts (unique)** | 3,599,832 | 3,599,832 | 100% |
| **allCountsMulti** | 16,089,488 | 16,089,488 | 100% |
| **filteredCountsMulti** | 4,503,920 | 4,503,920 | 100% |
| **Total values compared** | 820,118 | 820,118 | -- |
| **Value mismatches** | -- | **0** | -- |

All four count columns match exactly across all 63,086 genes for both unique and multi-mapper counts. A cell-by-cell comparison of the full duplication matrix (820,118 values) shows **zero mismatches** at a relative tolerance of 1e-6. Model fit parameters (intercept and slope) match to at least 10 significant digits.

### Side-by-side plots

The plots below compare the R dupRadar output (left) with RustQC output (right) for the large benchmark.

#### Density scatter plot

<div style="display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; margin: 1rem 0;">
  <figure style="margin: 0; text-align: center;">
    <img src="/RustQC/benchmarks/large/dupRadar/duprateExpDens.png" alt="dupRadar (R) density scatter plot" style="width: 100%;" />
    <figcaption>dupRadar (R)</figcaption>
  </figure>
  <figure style="margin: 0; text-align: center;">
    <img src="/RustQC/benchmarks/large/RustQC/GM12878_REP1.markdup.sorted_duprateExpDens.png" alt="RustQC density scatter plot" style="width: 100%;" />
    <figcaption>RustQC</figcaption>
  </figure>
</div>

#### Duplication rate boxplot

<div style="display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; margin: 1rem 0;">
  <figure style="margin: 0; text-align: center;">
    <img src="/RustQC/benchmarks/large/dupRadar/duprateExpBoxplot.png" alt="dupRadar (R) duplication rate boxplot" style="width: 100%;" />
    <figcaption>dupRadar (R)</figcaption>
  </figure>
  <figure style="margin: 0; text-align: center;">
    <img src="/RustQC/benchmarks/large/RustQC/GM12878_REP1.markdup.sorted_duprateExpBoxplot.png" alt="RustQC duplication rate boxplot" style="width: 100%;" />
    <figcaption>RustQC</figcaption>
  </figure>
</div>

#### Expression histogram

<div style="display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; margin: 1rem 0;">
  <figure style="margin: 0; text-align: center;">
    <img src="/RustQC/benchmarks/large/dupRadar/expressionHist.png" alt="dupRadar (R) expression histogram" style="width: 100%;" />
    <figcaption>dupRadar (R)</figcaption>
  </figure>
  <figure style="margin: 0; text-align: center;">
    <img src="/RustQC/benchmarks/large/RustQC/GM12878_REP1.markdup.sorted_expressionHist.png" alt="RustQC expression histogram" style="width: 100%;" />
    <figcaption>RustQC</figcaption>
  </figure>
</div>


