---
title: Combined Benchmark
description: Overall performance comparison between the traditional R workflow (featureCounts + dupRadar + RSeQC) and RustQC's single-pass approach.
---

The traditional RNA-seq QC workflow requires multiple separate tools:
**featureCounts** for gene-level read counting and biotype quantification,
**dupRadar** for duplication rate analysis, and
**RSeQC** for a suite of quality control metrics (strandedness, read distribution,
junction analysis, inner distance, and more). RustQC replaces all three in a single
pass, producing read counts, biotype summaries, duplication metrics, and RSeQC-equivalent
outputs together.

## Traditional workflow vs. RustQC

**Large benchmark input:** GM12878 REP1 -- a 10 GB paired-end RNA-seq BAM aligned
to GRCh38 (63,086 genes).

| Step | Traditional workflow | RustQC |
|------|--------------------:|-------:|
| Read counting (featureCounts) | 3m 39s | -- |
| Duplication analysis (dupRadar) | 27m 21s | -- |
| bam_stat (RSeQC) | 6m 07s | -- |
| infer_experiment (RSeQC) | 7s | -- |
| read_duplication (RSeQC) | 29m 43s | -- |
| read_distribution (RSeQC) | 6m 00s | -- |
| junction_annotation (RSeQC) | 4m 37s | -- |
| junction_saturation (RSeQC) | 6m 32s | -- |
| inner_distance (RSeQC) | 1m 09s | -- |
| Biotype summaries | Additional scripting | -- |
| **All outputs, single pass** | -- | **16m 11s** |
| **Total** | **~1h 25m** | **16m 11s** |

RustQC produces all outputs -- dupRadar duplication matrix, model fit, plots,
featureCounts-compatible counts, assignment summary, biotype counts, RSeQC-equivalent
metrics (bam_stat, infer_experiment, read_duplication, read_distribution,
junction_annotation, junction_saturation, inner_distance), and MultiQC
files -- in a single pass over the BAM file.

The traditional workflow requires 9 separate tool invocations, each reading
the BAM file independently. RustQC replaces all of them with a single command.

## Where the speedup comes from

1. **Single-pass architecture:** The R workflow reads the BAM file multiple times --
   once for featureCounts, once for dupRadar, and once per RSeQC tool. RustQC reads
   it once and produces everything in that single pass.
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
