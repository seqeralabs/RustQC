---
title: Configuration File
description: YAML configuration file reference for RustQC, covering chromosome mapping, output toggles, and biotype settings.
---

The `rustqc rna` subcommand supports an optional YAML configuration file for
advanced settings that go beyond what CLI flags offer. Pass the config file with
`--config` / `-c`:

```bash
rustqc rna sample.bam --gtf genes.gtf -p -c config.yaml -o results/
# or with BED (RSeQC tools only)
rustqc rna sample.bam --bed genes.bed -p -c config.yaml -o results/
```

All sections and fields are optional. Missing fields use their default values.
Unknown fields are silently ignored, so config files remain forward-compatible.

:::note
Since all analyses run within the `rustqc rna` command, the configuration file
controls all tools: dupRadar, featureCounts, all 7 RSeQC tools, TIN, Qualimap,
preseq, and samtools. Each tool has an `enabled` toggle and tool-specific
parameter overrides where applicable.
CLI flags take precedence over config file values.
:::

## Full example

```yaml
# Library settings
stranded: 0 # 0=unstranded, 1=forward, 2=reverse
paired: true # Paired-end mode

# Chromosome name remapping
chromosome_prefix: "chr"
chromosome_mapping:
  chrM: "MT"

# Output directory layout
flat_output: false # Set to true to skip subdirectories

# dupRadar output toggles
dupradar:
  dup_matrix: true
  intercept_slope: true
  density_scatter_plot: true
  boxplot: true
  expression_histogram: true
  multiqc_intercept: true
  multiqc_curve: true

# featureCounts output toggles
featurecounts:
  counts_file: true
  summary_file: true
  biotype_counts: true
  biotype_counts_mqc: true
  biotype_rrna_mqc: true
  biotype_attribute: "gene_biotype"

# RSeQC tool toggles and settings
bam_stat:
  enabled: true
infer_experiment:
  enabled: true
  sample_size: 200000
read_duplication:
  enabled: true
read_distribution:
  enabled: true
junction_annotation:
  enabled: true
  min_intron: 50
junction_saturation:
  enabled: true
  min_intron: 50
  min_coverage: 1
  percentile_floor: 5
  percentile_ceiling: 100
  percentile_step: 5
inner_distance:
  enabled: true
  sample_size: 1000000
  lower_bound: -250
  upper_bound: 250
  step: 5

# TIN (Transcript Integrity Number)
tin:
  enabled: true
  sample_size: 100
  min_coverage: 10

# Qualimap RNA-Seq QC
qualimap:
  enabled: true

# Library complexity (preseq lc_extrap)
preseq:
  enabled: true
  max_extrap: 10000000000
  step_size: 1000000
  n_bootstraps: 100
  confidence_level: 0.95
  seed: 408
  max_terms: 100
  defects: false

# Samtools-compatible outputs
flagstat:
  enabled: true
idxstats:
  enabled: true
samtools_stats:
  enabled: true
```

## Library settings

### `stranded`

Library strandedness for strand-aware read counting:

| Value | Meaning                                                 |
| ----- | ------------------------------------------------------- |
| `0`   | Unstranded (count reads on either strand)               |
| `1`   | Forward stranded (read 1 maps to the transcript strand) |
| `2`   | Reverse stranded (read 2 maps to the transcript strand) |

**Default:** `0` (unstranded)

This can also be set via the `-s` / `--stranded` CLI flag, which takes
precedence over the config file value.

### `paired`

Enable paired-end mode. When `true`, read pairs are counted as a single
fragment.

**Default:** `false` (single-end mode)

This can also be set via the `-p` / `--paired` CLI flag. If either the CLI flag
or the config file value is `true`, paired-end mode is enabled.

```yaml
stranded: 2
paired: true
```

## Chromosome name mapping

When the chromosome names in your alignment file differ from those in the GTF,
RustQC can remap them automatically. Two mechanisms are available, and they can
be combined.

### `chromosome_prefix`

A string prefix to prepend to alignment file chromosome names before matching
against GTF names. Applied first, before explicit mapping lookups.

```yaml
# Alignment has "1", "2", "X"; GTF has "chr1", "chr2", "chrX"
chromosome_prefix: "chr"
```

### `chromosome_mapping`

An explicit mapping from GTF chromosome names (keys) to alignment file chromosome
names (values). Applied after the prefix, so explicit entries can override the
prefix for specific chromosomes.

```yaml
# After adding "chr" prefix, override the mitochondrial chromosome
chromosome_mapping:
  chrM: "MT"
```

A common use case is GENCODE GTFs (which use `chr1`, `chr2`, ...) with Ensembl
alignments (which use `1`, `2`, ...):

```yaml
chromosome_prefix: "chr"
chromosome_mapping:
  chrM: "MT"
```

## dupRadar output toggles

The `dupradar:` section controls which dupRadar output files are generated.
All outputs are **enabled by default**.

```yaml
dupradar:
  dup_matrix: true # Duplication matrix TSV
  intercept_slope: true # Intercept/slope fit results
  density_scatter_plot: true # Density scatter plot (PNG + SVG)
  boxplot: true # Duplication rate boxplot (PNG + SVG)
  expression_histogram: true # Expression histogram (PNG + SVG)
  multiqc_intercept: true # MultiQC intercept/slope file
  multiqc_curve: true # MultiQC fitted curve file
```

Set any field to `false` to skip generating that output:

```yaml
dupradar:
  boxplot: false
  expression_histogram: false
```

## featureCounts output toggles

The `featurecounts:` section controls which featureCounts-compatible output files
are generated, plus the biotype attribute setting.

```yaml
featurecounts:
  counts_file: true # featureCounts-compatible counts TSV
  summary_file: true # Assignment summary file
  biotype_counts: true # Biotype counts TSV
  biotype_counts_mqc: true # Biotype counts MultiQC bargraph file
  biotype_rrna_mqc: true # Biotype rRNA percentage MultiQC file
  biotype_attribute: "gene_biotype" # GTF attribute for biotype grouping
```

### `biotype_attribute`

The GTF attribute name used for biotype grouping. This controls how genes are
categorized in the biotype output files.

| GTF source | Typical attribute |
| ---------- | ----------------- |
| Ensembl    | `gene_biotype`    |
| GENCODE    | `gene_type`       |

**Default:** `"gene_biotype"`

This can also be set via the `--biotype-attribute` CLI flag, which takes
precedence over the config file value.

RustQC auto-detects the biotype attribute if the specified one is not found in
the GTF. If neither `gene_biotype` nor `gene_type` is present, a warning is
printed and biotype counting is skipped.

## RSeQC tool settings

Each of the 7 RSeQC tools has an `enabled` toggle (default `true`) and
tool-specific parameter overrides. Disabling a tool here prevents it from
running even when annotation is provided. CLI flags take precedence over
config file values for all parameters.

### bam_stat

```yaml
bam_stat:
  enabled: true # Set to false to skip bam_stat
```

No additional parameters. This tool does not require annotation.

### infer_experiment

```yaml
infer_experiment:
  enabled: true
  sample_size: 200000 # Number of reads to sample (default: 200000)
```

Requires annotation (`--gtf` or `--bed`). The `sample_size` can also be set via
`--infer-experiment-sample-size`.

### read_duplication

```yaml
read_duplication:
  enabled: true # Set to false to skip read_duplication
```

No additional parameters. This tool does not require annotation.

### read_distribution

```yaml
read_distribution:
  enabled: true # Set to false to skip read_distribution
```

No additional parameters. Requires annotation (`--gtf` or `--bed`).

### junction_annotation

```yaml
junction_annotation:
  enabled: true
  min_intron: 50 # Minimum intron length in bases (default: 50)
```

Requires annotation (`--gtf` or `--bed`). The `min_intron` can also be set via `--min-intron`.

### junction_saturation

```yaml
junction_saturation:
  enabled: true
  min_intron: 50 # Minimum intron length in bases (default: 50)
  min_coverage: 1 # Minimum read count to consider a junction (default: 1)
  percentile_floor: 5 # Sampling start percentage (default: 5)
  percentile_ceiling: 100 # Sampling end percentage (default: 100)
  percentile_step: 5 # Sampling step size (default: 5)
```

Requires annotation (`--gtf` or `--bed`). These parameters can also be set via CLI flags:
`--min-intron`, `--junction-saturation-min-coverage`,
`--junction-saturation-percentile-floor`, `--junction-saturation-percentile-ceiling`,
`--junction-saturation-percentile-step`.

### inner_distance

```yaml
inner_distance:
  enabled: true
  sample_size: 1000000 # Number of reads to sample (default: 1000000)
  lower_bound: -250 # Histogram lower bound (default: -250)
  upper_bound: 250 # Histogram upper bound (default: 250)
  step: 5 # Histogram bin width (default: 5)
```

Requires annotation (`--gtf` or `--bed`). These parameters can also be set via CLI flags:
`--inner-distance-sample-size`, `--inner-distance-lower-bound`,
`--inner-distance-upper-bound`, `--inner-distance-step`.

### tin

```yaml
tin:
  enabled: true
  sample_size: 100 # Equally-spaced positions to sample per transcript (default: 100)
  min_coverage: 10 # Minimum read-start count to compute TIN (default: 10)
```

Requires annotation (`--gtf` or `--bed`). The TIN (Transcript Integrity Number)
measures transcript integrity via Shannon entropy of read coverage uniformity.

## qualimap

```yaml
qualimap:
  enabled: true # Set to false to skip Qualimap RNA-Seq QC analysis
```

Requires annotation (`--gtf` only). Runs the Qualimap RNA-Seq QC analysis:
gene body coverage profiling (100 percentile bins, 5'→3'), 5'/3' bias metrics,
read origin classification (exonic/intronic/intergenic), strand-specificity
estimation, and splice junction motif counting. Produces Qualimap-compatible
output files parseable by MultiQC.

## preseq

```yaml
preseq:
  enabled: true
  max_extrap: 10000000000 # Maximum extrapolation depth (default: 1e10)
  step_size: 1000000 # Step size between extrapolation points (default: 1e6)
  n_bootstraps: 100 # Bootstrap replicates for confidence intervals (default: 100)
  confidence_level: 0.95 # CI confidence level (default: 0.95)
  seed: 408 # Random seed for reproducibility (default: 408, matching upstream preseq)
  max_terms: 100 # Maximum terms in power series (default: 100)
  defects: false # Use defects model for problematic histograms (default: false)
```

Runs in both GTF and BED modes (only needs BAM fragment info). The `max_extrap`,
`step_size`, and `n_bootstraps` can also be set via CLI flags (`--preseq-max-extrap`,
`--preseq-step-size`, `--preseq-n-bootstraps`). Use `--skip-preseq` to disable entirely.

## samtools

### flagstat / idxstats / samtools_stats

```yaml
flagstat:
  enabled: true # samtools flagstat-compatible output
idxstats:
  enabled: true # samtools idxstats-compatible output
samtools_stats:
  enabled: true # samtools stats compatible output (full format including all histogram sections)
```

These produce samtools-compatible output files in the `samtools/` subdirectory.
They share the same BAM statistics accumulator as `bam_stat` -- enabling any of
them ensures the statistics are collected.
