---
title: Performance Benchmarks
description: Benchmark results comparing RustQC to R dupRadar, showing up to 33x speedup with identical output across all values.
---

RustQC produces output identical to the R dupRadar package while running
significantly faster. The benchmarks below compare RustQC against the original
R implementation on the same input data.

In addition to dupRadar outputs, RustQC also generates featureCounts-compatible
output files (counts TSV, summary, biotype counts, and MultiQC files) in the
same pass -- the R workflow requires a separate featureCounts run.

## Small benchmark

**Input:** A test BAM file with a chr6-only GTF annotation (2,905 genes).

| Metric | dupRadar (R) | RustQC |
|--------|-------------|--------|
| **Runtime** | 2.50s | 0.25s |
| **Speedup** | -- | **10x** |
| **Intercept** | 0.03186 | 0.03186 |
| **Slope** | 1.60189 | 1.60189 |
| **Genes total** | 2,905 | 2,905 |
| **Genes with reads** | 636 | 636 |
| **Genes with duplicates** | 201 | 201 |
| **Total values compared** | 37,765 | 37,765 |
| **Value mismatches** | -- | **0** |

### Count comparison (small)

| Metric | dupRadar (R) | RustQC | Exact match |
|--------|------------:|-------:|:-----------:|
| **allCounts (unique)** | 20,449 | 20,449 | 100% |
| **filteredCounts (unique)** | 17,879 | 17,879 | 100% |
| **allCountsMulti** | 22,812 | 22,812 | 100% |
| **filteredCountsMulti** | 20,034 | 20,034 | 100% |

## Large benchmark

**Input:** GM12878 REP1 -- a full-size RNA-seq BAM (~10 GB) from the
[nf-core/rnaseq](https://nf-co.re/rnaseq) pipeline, duplicate-marked with
Picard. Paired-end, unstranded, aligned to GRCh38 (63,086 genes).

| Metric | dupRadar (R) | RustQC (1 thread) | RustQC (8 threads) | RustQC (10 threads) |
|--------|-------------|-------------------|--------------------|--------------------|
| **Runtime** | 29m 56s | 3m 16s | 1m 03s | 0m 54s |
| **Speedup** | -- | ~9x | ~28x | **~33x** |
| **Max RSS** | N/A | 503 MB | 893 MB | 1.3 GB |

| Metric | dupRadar (R) | RustQC |
|--------|-------------|--------|
| **Intercept** | 0.8245 | 0.8245 |
| **Slope** | 1.6774 | 1.6774 |
| **Genes total** | 63,086 | 63,086 |
| **Genes with reads (unique)** | 23,597 | 23,597 |
| **Genes with reads (multi)** | 24,719 | 24,719 |
| **Total values compared** | 820,118 | 820,118 |
| **Value mismatches** | -- | **0** |

### Count comparison (large)

| Metric | dupRadar (R) | RustQC | Exact match |
|--------|------------:|-------:|:-----------:|
| **allCounts (unique)** | 14,654,579 | 14,654,579 | 100% |
| **filteredCounts (unique)** | 3,599,832 | 3,599,832 | 100% |
| **allCountsMulti** | 16,089,488 | 16,089,488 | 100% |
| **filteredCountsMulti** | 4,503,920 | 4,503,920 | 100% |

All four count columns match exactly across all 63,086 genes for both unique and
multi-mapper counts. A cell-by-cell comparison of the full duplication matrix
(820,118 values) shows **zero mismatches** at a relative tolerance of 1e-6.

Model fit parameters (intercept and slope) match to at least 10 significant
digits.

## Where the speedup comes from

RustQC's performance advantage has two sources:

1. **Rust vs. R:** Compiled, zero-cost abstractions, and efficient memory
   management provide a baseline speedup for the same algorithm.
2. **Single-pass architecture:** The R dupRadar workflow requires running
   featureCounts as a separate step before dupRadar processes the counts.
   RustQC performs read counting, matrix construction, model fitting, and plot
   generation in a single pass over the BAM file.

Multi-threaded scaling depends on the number of chromosomes with mapped reads
and the evenness of their read distribution. The parallelism is applied across
chromosomes within each BAM file.

## Notes on benchmark conditions

Runtime may vary depending on hardware. The times above were measured on a
single machine (10-core Apple Silicon) for relative comparison. The R dupRadar
benchmark was run via Docker (x86 emulation on ARM Mac), so the R timings
include container overhead and emulation penalties.

To reproduce these benchmarks, see the instructions in the
[benchmark directory](https://github.com/ewels/RustQC/tree/main/benchmark) of
the repository.
