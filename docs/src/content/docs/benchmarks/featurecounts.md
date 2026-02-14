---
title: featureCounts Benchmark
description: Performance and accuracy comparison between RustQC's built-in read counting and Rsubread featureCounts.
---

RustQC includes built-in gene-level read counting that produces output
compatible with the Subread featureCounts format. This page compares the
counting performance and output accuracy against standalone Rsubread
featureCounts.

## Performance

**Large benchmark input:** GM12878 REP1 -- a 10 GB paired-end RNA-seq BAM
aligned to GRCh38 (63,086 genes).

| Tool | Runtime |
|------|---------|
| Rsubread featureCounts | 3m 39s |
| RustQC (all tools, single pass) | 16m 11s |

RustQC's single pass includes featureCounts-compatible counting alongside
dupRadar duplication analysis and all 7 RSeQC tools. The traditional workflow
requires running each tool separately. Standalone featureCounts timing is from
Docker with x86 emulation on ARM Mac.

## Output equivalence

RustQC's read counting uses the same algorithm as Subread featureCounts:

- **Feature type:** exon-level features grouped by `gene_id`
- **Overlap detection:** at least 1 base pair overlap
- **Strand awareness:** configurable via the `-s` / `--strandedness` flag
- **Multi-mapping:** tracked separately for unique and multi-mapper columns

### Count comparison (large benchmark)

| Metric | Rsubread featureCounts | RustQC | Exact match |
|--------|----------------------:|-------:|:-----------:|
| **allCounts (unique)** | 14,654,579 | 14,654,579 | 100% |
| **filteredCounts (unique)** | 3,599,832 | 3,599,832 | 100% |
| **allCountsMulti** | 16,089,488 | 16,089,488 | 100% |
| **filteredCountsMulti** | 4,503,920 | 4,503,920 | 100% |

Gene-level read counts are identical across all 63,086 genes. Assignment
statistics (Assigned, Unassigned_NoFeatures, Unassigned_Ambiguous) match
exactly.

The output format is directly compatible with downstream tools such as DESeq2
and MultiQC.

### Side-by-side plots

Since dupRadar visualizations are derived from the same gene-level counts, identical plots confirm that the underlying feature counting is correct. The plots below compare R dupRadar output (left) with RustQC output (right) for the large benchmark.

#### Density scatter plot

<div style="display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; margin: 1rem 0;">
  <figure style="margin: 0; text-align: center;">
    <img src="/benchmarks/large/dupRadar/duprateExpDens.png" alt="dupRadar (R) density scatter plot" style="width: 100%;" />
    <figcaption>dupRadar (R)</figcaption>
  </figure>
  <figure style="margin: 0; text-align: center;">
    <img src="/benchmarks/large/RustQC/GM12878_REP1.markdup.sorted_duprateExpDens.png" alt="RustQC density scatter plot" style="width: 100%;" />
    <figcaption>RustQC</figcaption>
  </figure>
</div>

#### Duplication rate boxplot

<div style="display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; margin: 1rem 0;">
  <figure style="margin: 0; text-align: center;">
    <img src="/benchmarks/large/dupRadar/duprateExpBoxplot.png" alt="dupRadar (R) duplication rate boxplot" style="width: 100%;" />
    <figcaption>dupRadar (R)</figcaption>
  </figure>
  <figure style="margin: 0; text-align: center;">
    <img src="/benchmarks/large/RustQC/GM12878_REP1.markdup.sorted_duprateExpBoxplot.png" alt="RustQC duplication rate boxplot" style="width: 100%;" />
    <figcaption>RustQC</figcaption>
  </figure>
</div>

#### Expression histogram

<div style="display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; margin: 1rem 0;">
  <figure style="margin: 0; text-align: center;">
    <img src="/benchmarks/large/dupRadar/expressionHist.png" alt="dupRadar (R) expression histogram" style="width: 100%;" />
    <figcaption>dupRadar (R)</figcaption>
  </figure>
  <figure style="margin: 0; text-align: center;">
    <img src="/benchmarks/large/RustQC/GM12878_REP1.markdup.sorted_expressionHist.png" alt="RustQC expression histogram" style="width: 100%;" />
    <figcaption>RustQC</figcaption>
  </figure>
</div>

## Additional outputs

Beyond the standard featureCounts counts file and summary, RustQC also
produces:

- **Biotype counts** (`.biotype_counts.tsv`) -- per-biotype read count
  summaries
- **Biotype MultiQC bargraph** (`.biotype_counts_mqc.tsv`) -- ready for
  MultiQC visualization
- **rRNA percentage** (`.biotype_counts_rrna_mqc.tsv`) -- rRNA fraction for
  MultiQC general statistics

Generating these in the traditional workflow requires additional scripting after
the featureCounts run.
