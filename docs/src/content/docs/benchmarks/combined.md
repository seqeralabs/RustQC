---
title: Combined Benchmark
description: Overall performance comparison between the traditional workflow and RustQC's single-pass approach.
---

The traditional RNA-seq QC workflow requires multiple separate tools:
**featureCounts** for gene-level read counting and biotype quantification,
**dupRadar** for duplication rate analysis,
**RSeQC** for a suite of quality control metrics,
**samtools** for alignment statistics,
**preseq** for library complexity estimation,
**RSeQC tin.py** for transcript integrity, and
**Qualimap** for gene body coverage.
RustQC replaces all of these in a single pass, producing every output together.

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
| TIN (RSeQC tin.py) | ~20m | -- |
| samtools flagstat | ~2m | -- |
| samtools idxstats | <1s | -- |
| samtools stats | ~5m | -- |
| preseq lc_extrap | ~3m | -- |
| Qualimap rnaseq | ~10m | -- |
| Biotype summaries | Additional scripting | -- |
| **All outputs, single pass** | -- | **~5m** |
| **Total** | **~2h 5m** | **~5m** |

RustQC produces all outputs -- dupRadar duplication matrix, model fit, plots,
featureCounts-compatible counts, assignment summary, biotype counts, all 7
RSeQC-equivalent metrics, TIN scores, samtools flagstat/idxstats/stats,
preseq library complexity extrapolation, Qualimap-compatible gene body coverage,
and MultiQC files -- in a single pass over the BAM file.

The traditional workflow requires 15+ separate tool invocations, each reading
the BAM file independently. RustQC replaces all of them with a single command.

## Where the speedup comes from

1. **Single-pass architecture:** The traditional workflow reads the BAM file
   many times -- once per tool. RustQC reads it once and produces everything
   in that single pass.
2. **Compiled Rust vs. interpreted R/Python:** Compiled code with zero-cost
   abstractions and efficient memory management.
3. **Multi-threaded parallelism:** RustQC parallelizes across chromosomes.
   Scaling depends on the number of chromosomes with mapped reads and the
   evenness of their read distribution.

## Output equivalence

Every output file produced by RustQC matches the original tools:

- **820,118 duplication matrix values** compared with zero mismatches
- **Gene-level read counts** identical across all 63,086 genes
- **Model fit parameters** (intercept and slope) match to 10+ significant digits
- **Assignment statistics** (Assigned, NoFeatures, Ambiguous) match exactly
- **flagstat** all 16 metrics identical to samtools flagstat
- **idxstats** all per-chromosome counts identical to samtools idxstats
- **stats** core SN metrics (read counts, lengths, duplicates) match samtools stats
- **preseq** extrapolation curve within <0.1% of preseq v3.2.0 across entire range
- **RSeQC tools** all values match Python RSeQC (see [RSeQC benchmark](rseqc/))

See the individual benchmark pages for detailed per-tool comparisons:
[dupRadar](dupradar/), [featureCounts](featurecounts/), [RSeQC](rseqc/),
[preseq](preseq/), [Samtools](samtools/).

## Benchmark conditions

- **Hardware:** 10-core Apple Silicon Mac (M3 Max, 128 GB RAM).
- **RustQC:** 10 threads, `--gtf` mode (all tools enabled).
- **Traditional tools:** Run via Docker with x86 emulation on ARM Mac. Absolute
  timings include container overhead and emulation penalties. The relative
  comparison demonstrates the magnitude of the speedup.
- **Reproducibility:** Benchmark scripts and data are in the
  [benchmark directory](https://github.com/ewels/RustQC/tree/main/benchmark).
