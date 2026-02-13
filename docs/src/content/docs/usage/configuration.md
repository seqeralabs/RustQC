---
title: Configuration File
description: YAML configuration file reference for RustQC, covering chromosome mapping, output toggles, and biotype settings.
---

The `rustqc rna` subcommand supports an optional YAML configuration file for
advanced settings that go beyond what CLI flags offer. Pass the config file with
`--config` / `-c`:

```bash
rustqc rna sample.bam --gtf genes.gtf -p -c config.yaml -o results/
```

All sections and fields are optional. Missing fields use their default values.
Unknown fields are silently ignored, so config files remain forward-compatible.

:::note
The configuration file currently applies only to the `rna` subcommand (dupRadar
and featureCounts analysis). The RSeQC subcommands are configured entirely via
their CLI flags.
:::

## Full example

```yaml
# Chromosome name remapping
chromosome_prefix: "chr"
chromosome_mapping:
  chrM: "MT"

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
  dup_matrix: true              # Duplication matrix TSV
  intercept_slope: true         # Intercept/slope fit results
  density_scatter_plot: true    # Density scatter plot (PNG + SVG)
  boxplot: true                 # Duplication rate boxplot (PNG + SVG)
  expression_histogram: true    # Expression histogram (PNG + SVG)
  multiqc_intercept: true       # MultiQC intercept/slope file
  multiqc_curve: true           # MultiQC fitted curve file
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
  counts_file: true           # featureCounts-compatible counts TSV
  summary_file: true          # Assignment summary file
  biotype_counts: true        # Biotype counts TSV
  biotype_counts_mqc: true    # Biotype counts MultiQC bargraph file
  biotype_rrna_mqc: true      # Biotype rRNA percentage MultiQC file
  biotype_attribute: "gene_biotype"  # GTF attribute for biotype grouping
```

### `biotype_attribute`

The GTF attribute name used for biotype grouping. This controls how genes are
categorized in the biotype output files.

| GTF source | Typical attribute |
|------------|-------------------|
| Ensembl | `gene_biotype` |
| GENCODE | `gene_type` |

**Default:** `"gene_biotype"`

This can also be set via the `--biotype-attribute` CLI flag, which takes
precedence over the config file value.

RustQC auto-detects the biotype attribute if the specified one is not found in
the GTF. If neither `gene_biotype` nor `gene_type` is present, a warning is
printed and biotype counting is skipped.
