---
title: featureCounts Output
description: How RustQC produces featureCounts-compatible read count files and biotype summaries.
---

RustQC includes built-in gene-level read counting that produces output compatible with the [Subread featureCounts](http://subread.sourceforge.net/) format. This happens in the same pass as the duplicate rate analysis, eliminating the need for a separate featureCounts run.

## How it works

During the BAM processing pass, RustQC simultaneously:

1. Counts reads assigned to each gene (using the same algorithm as featureCounts)
2. Tracks duplication rates for the dupRadar analysis
3. Aggregates counts by gene biotype

This single-pass approach is one of the key reasons RustQC is faster than running dupRadar and featureCounts separately.

## Output files

### Counts file (`.featureCounts.tsv`)

A tab-separated file with the standard featureCounts format:

```
# Program:RustQC v0.1.0; Command: rustqc rna sample.bam --gtf genes.gtf -p
Geneid  Chr     Start   End     Strand  Length  sample.bam
ENSG00000000003 chrX    100627108;100629986;...  100636806;100637104;... -;-;... 3768    521
ENSG00000000005 chrX    100584936;100585053;...  100585091;100599885;... +;+;... 1339    0
```

This file is directly compatible with tools that accept featureCounts output, such as [DESeq2](https://bioconductor.org/packages/DESeq2/) and [MultiQC](https://multiqc.info/).

### Summary file (`.featureCounts.tsv.summary`)

Assignment statistics in featureCounts summary format:

```
Status  sample.bam
Assigned        22812
Unassigned_Unmapped     0
Unassigned_NoFeatures   1227
Unassigned_Ambiguous    2395
```

### Biotype counts (`.biotype_counts.tsv`)

Gene counts aggregated by biotype:

```
biotype count
protein_coding  18234
lincRNA 2103
processed_pseudogene    891
```

### Biotype MultiQC files

- `.biotype_counts_mqc.tsv` - Biotype distribution for MultiQC bargraph
- `.biotype_rrna_mqc.tsv` - rRNA percentage for MultiQC general statistics table

## Biotype attribute detection

RustQC automatically detects the biotype attribute in your GTF file:

- **Ensembl GTFs** use `gene_biotype`
- **GENCODE GTFs** use `gene_type`

If neither is found, biotype outputs are skipped with a warning. You can override the auto-detection with the `--biotype-attribute` flag or the `featurecounts.biotype_attribute` config option.

## Configuring outputs

Use a YAML configuration file to control which featureCounts outputs are generated:

```yaml
featurecounts:
  counts_file: true          # .featureCounts.tsv
  summary_file: true         # .featureCounts.tsv.summary
  biotype_counts: true       # .biotype_counts.tsv
  biotype_counts_mqc: true   # .biotype_counts_mqc.tsv
  biotype_rrna_mqc: true     # .biotype_rrna_mqc.tsv
```

All outputs default to `true`. Set any to `false` to skip generation.

## Differences from Subread featureCounts

RustQC's read counting follows the same algorithm as featureCounts with these default settings:

- Feature type: `exon`
- Attribute: `gene_id`
- Overlap detection: at least 1 base overlap
- Multi-mapping reads: counted (for multi columns), excluded (for unique columns)
- Strand-aware counting based on the `-s` flag

The key difference is that RustQC performs counting and duplicate analysis in a single pass, while the traditional workflow requires running featureCounts separately before dupRadar.
