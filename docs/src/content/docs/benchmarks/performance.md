---
title: Performance & Tuning
description: Per-tool cost breakdown, cumulative removal analysis, and guidance on tuning RustQC for speed by disabling unneeded tools or adjusting thread counts.
---

RustQC runs all RNA-seq QC analyses in a single pass over the BAM file.
Because every tool shares the same BAM decode and per-read dispatch,
the cost of each tool is **not** simply additive -- there is a shared baseline
cost for reading and decoding the BAM that all tools amortise together.

This page documents how much each individual tool contributes to the total
runtime, how disabling tools affects performance, and how to tune thread
counts for your hardware.

## Recommendations

- **Leave all tools enabled** unless you have a specific reason not to.
  The full run completes in ~10 minutes for a large BAM -- far faster than
  running the equivalent traditional tools separately.
- **Disabling TIN** provides the single largest speedup (roughly halves
  runtime). Disabling read_duplication alongside TIN brings it under
  5 minutes.
- **Thread count:** `-t 8` to `-t 25` is the productive range for human
  data. Beyond the number of chromosomes, extra threads are mostly idle.
  For small genomes, the ceiling is lower (see [CPU count table](#effect-of-cpu-count) below).
- **Multiple BAMs:** the useful thread count scales with the number of
  input files (e.g. 4 human BAMs can use ~100 threads).
- **Ensure a BAM index exists** (`.bai`/`.csi`). Without one, only a
  single counting thread is used regardless of `-t`.

## Per-tool cost

### Benchmark setup

The per-tool benchmarks were measured with the following setup:

| Parameter        | Value                                                         |
| ---------------- | ------------------------------------------------------------- |
| **BAM file**     | GM12878 REP1, paired-end RNA-seq, GRCh38 (~186M reads, 10 GB) |
| **Annotation**   | GENCODE v46 GTF + BED12 (63,677 genes)                        |
| **Hardware**     | Apple M1 Pro, 10 cores, 16 GB RAM                             |
| **Threads**      | 8 (`-t 8`)                                                    |
| **Build**        | `cargo build --release` (LTO, codegen-units=1, opt-level=3)   |
| **Strandedness** | Reverse (`-s 2`)                                              |
| **Mode**         | Paired-end (`-p`), all tools enabled by default               |

The **baseline** is a run with every tool disabled (BAM decode + dispatch
overhead only): **34 seconds**.

### Results

Each tool was benchmarked in isolation -- a run with only that single tool
enabled, all others disabled. The "marginal cost" is the difference between
the single-tool run and the no-tools baseline.

| Tool                  | Run time | Marginal cost |
| --------------------- | -------: | ------------: |
| read_duplication      |      95s |           61s |
| infer_experiment      |      53s |           19s |
| bam_stat              |      52s |           18s |
| samtools_stats        |      51s |           17s |
| idxstats              |      50s |           16s |
| flagstat              |      50s |           16s |
| qualimap              |      48s |           14s |
| junction_annotation   |      47s |           13s |
| inner_distance        |      45s |           11s |
| preseq                |      40s |            6s |
| tin                   |      39s |            5s |
| read_distribution     |      35s |            1s |
| junction_saturation   |      35s |            1s |
| dupradar              |      35s |            1s |
| featurecounts         |      33s |           ~0s |
| **All tools enabled** | **631s** |               |

> **Why don't the marginal costs add up to 631s?**
> Two reasons. First, the shared BAM decode baseline (34s) is paid once
> regardless of how many tools are active. Second -- and more importantly --
> **TIN** runs as a separate post-counting phase that performs random-access
> BAM lookups for every transcript. When TIN is enabled, this additional
> pass dominates the total runtime. The single-tool cost of TIN appears low
> (5s marginal) because the TIN-only run processes far fewer transcripts
> than a full run where all counting tools feed transcript-level data into
> the TIN phase.

## All tools enabled

With every tool active, the total runtime is **631 seconds** (10 minutes 31 seconds).

## Cumulative removal: stripping tools from slowest to fastest

Starting from a full run with all 15 tools, tools were removed one at a time
in order of their individual marginal cost (most expensive first). This shows
how the total runtime decreases as tools are stripped away.

| Step | Tool removed          | Tools remaining | Run time |     Delta |
| ---: | --------------------- | --------------: | -------: | --------: |
|    0 | _(none -- all tools)_ |              15 |     631s |           |
|    1 | read_duplication      |              14 |     390s |     -241s |
|    2 | infer_experiment      |              13 |     359s |      -31s |
|    3 | bam_stat              |              12 |     358s |       -1s |
|    4 | samtools_stats        |              11 |     356s |       -2s |
|    5 | idxstats              |              10 |     358s |       +2s |
|    6 | flagstat              |               9 |     340s |      -18s |
|    7 | qualimap              |               8 |     314s |      -26s |
|    8 | junction_annotation   |               7 |     300s |      -14s |
|    9 | inner_distance        |               6 |     284s |      -16s |
|   10 | preseq                |               5 |     276s |       -8s |
|   11 | tin                   |               4 |      36s | **-240s** |
|   12 | read_distribution     |               3 |      36s |        0s |
|   13 | junction_saturation   |               2 |      34s |       -2s |
|   14 | dupradar              |               1 |      34s |        0s |
|   15 | featurecounts         |               0 |      34s |        0s |

### Key observations

**TIN dominates the runtime.** The single largest drop occurs at step 11 when
TIN is removed -- the total plummets from 276s to 36s (a **240-second reduction**).
TIN runs as a separate post-counting pass that performs random-access BAM lookups
for every transcript. Even though its marginal cost when run alone appears small
(5s), in the context of a full run it becomes the dominant bottleneck because
the counting phase feeds it data for all transcripts.

**read_duplication is the most expensive counting-phase tool.** Removing it
saves 241s (step 1), cutting the runtime from 631s to 390s. This tool maintains
per-read sequence hashing for duplication detection, which is computationally
intensive.

**The "samtools" group (bam_stat, flagstat, idxstats, samtools_stats) is cheap.**
These tools add minimal overhead individually because they compute simple
per-read counters with no complex data structures.

**Post-TIN removal, the baseline overhead dominates.** Once TIN and the expensive
counting tools are stripped away, the runtime converges to the ~34s BAM decode
baseline.

## Disabling tools

Every analysis tool in RustQC can be individually disabled via the
[YAML configuration file](/usage/configuration/).

Example: disable TIN and read_duplication for a faster run:

```yaml
tin:
  enabled: false
read_duplication:
  enabled: false
```

Or disable all outputs for dupRadar (per-output granularity):

```yaml
dupradar:
  dup_matrix: false
  intercept_slope: false
  density_scatter_plot: false
  boxplot: false
  expression_histogram: false
  multiqc_intercept: false
  multiqc_curve: false
```

When a tool is disabled, its accumulators are never constructed, its per-read
processing is completely skipped, and no output files are written. This is a
true zero-cost disable, not just output suppression.

See the [Configuration File](/usage/configuration/) documentation for the
complete list of toggles.

## Effect of CPU count

RustQC parallelises by assigning chromosomes to worker threads -- each worker
opens its own BAM reader and processes its assigned chromosomes independently.
This means **the number of chromosomes is the hard ceiling** for useful threads.
Any threads beyond the chromosome count are used for BAM decompression, which
gives diminishing returns.

When multiple BAM files are provided, the thread budget is split across files,
so the ceiling scales with the number of inputs.

A BAM index (`.bai`/`.csi`) is required for parallel mode. Without one, RustQC
falls back to a single counting thread.

| Organism                  | Primary sequences | Useful threads (single BAM) |
| ------------------------- | ----------------: | --------------------------: |
| Human (GRCh38)            |                25 |                       ~25   |
| Mouse (GRCm39)            |                22 |                       ~22   |
| Zebrafish (GRCz11)        |                26 |                       ~26   |
| _Drosophila_ (BDGP6)     |                 7 |                        ~7   |
| _C. elegans_ (WBcel235)   |                 7 |                        ~7   |
| _S. cerevisiae_ (R64-1-1) |                17 |                       ~17   |

### CPU scaling benchmarks

<!-- TODO: Add CPU scaling benchmarks with actual timing data -->

_Benchmarks in progress -- this table will be updated with timing data._

| Threads | Run time | Speedup vs 1 thread |
| ------: | -------: | ------------------: |
|       1 |          |                1.0x |
|       2 |          |                     |
|       4 |          |                     |
|       8 |          |                     |
|      16 |          |                     |
|      25 |          |                     |
|      32 |          |                     |


