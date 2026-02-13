---
title: Combined Benchmark
description: Overall performance comparison between the traditional R workflow (featureCounts + dupRadar) and RustQC's single-pass approach.
---

The traditional RNA-seq QC workflow requires two separate tools run sequentially:
**featureCounts** for gene-level read counting and biotype quantification, then
**dupRadar** for duplication rate analysis. Biotype counting — essential for
assessing library quality (rRNA contamination, coding vs. non-coding ratios) —
is the primary reason pipelines run a standalone featureCounts step rather than
relying on dupRadar's internal counting alone. RustQC replaces both in a single
pass, producing read counts, biotype summaries, and duplication metrics together.

## Traditional workflow vs. RustQC

**Large benchmark input:** GM12878 REP1 -- a 10 GB paired-end RNA-seq BAM aligned
to GRCh38 (63,086 genes).

| Step | Traditional R workflow | RustQC |
|------|----------------------:|-------:|
| Read counting (featureCounts) | 16m 26s | -- |
| Duplication analysis (dupRadar) | 29m 56s | -- |
| Biotype summaries | Additional scripting | -- |
| **All outputs, single pass** | -- | **0m 54s** |
| **Total** | **46m 22s** | **0m 54s** |

RustQC produces all outputs -- dupRadar duplication matrix, model fit, plots,
featureCounts-compatible counts, assignment summary, biotype counts, and MultiQC
files -- in a single pass over the BAM file.

## Where the speedup comes from

1. **Single-pass architecture:** The R workflow reads the BAM file twice -- once
   for featureCounts, once for dupRadar. RustQC reads it once and produces
   everything in that single pass.
2. **Compiled Rust vs. interpreted R:** Compiled code with zero-cost abstractions
   and efficient memory management.
3. **Multi-threaded parallelism:** RustQC parallelizes across chromosomes.
   Scaling depends on the number of chromosomes with mapped reads and the
   evenness of their read distribution.

## Output equivalence

Every output file produced by RustQC matches the R tools exactly:

- **820,118 duplication matrix values** compared with zero mismatches
- **Gene-level read counts** identical across all 63,086 genes
- **Model fit parameters** (intercept and slope) match to 10+ significant digits
- **Assignment statistics** (Assigned, NoFeatures, Ambiguous) match exactly

See the [dupRadar benchmark](dupradar/) and
[featureCounts benchmark](featurecounts/) pages for detailed per-tool
comparisons, side-by-side plots, and count tables.

## Benchmark conditions

- **Hardware:** 10-core Apple Silicon Mac.
- **R benchmarks:** Run via Docker with x86 emulation on ARM Mac. Absolute R
  timings include container overhead and emulation penalties. The relative
  comparison demonstrates the magnitude of the speedup.
- **Reproducibility:** Benchmark scripts and data are in the
  [benchmark directory](https://github.com/ewels/RustQC/tree/main/benchmark).
