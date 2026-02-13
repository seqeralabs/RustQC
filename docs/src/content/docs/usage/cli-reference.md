---
title: CLI Reference
description: Complete command-line interface reference for RustQC, including all options, flags, and usage examples.
---

RustQC uses a subcommand-based CLI. The `rna` subcommand provides RNA-seq
duplicate rate analysis (dupRadar equivalent) with integrated featureCounts-compatible
output.

## Synopsis

```
rustqc rna <INPUT>... --gtf <GTF> [OPTIONS]
```

## Positional arguments

### `<INPUT>...`

One or more paths to duplicate-marked alignment files. Accepted formats are
**SAM**, **BAM**, and **CRAM**.

Multiple files can be passed and will be processed sequentially, with each file
producing its own set of output files. If threads are specified with `-t`, reads
within each file are counted in parallel across chromosomes.

```bash
# Single file
rustqc rna sample.bam --gtf genes.gtf

# Multiple files
rustqc rna sample1.bam sample2.bam sample3.bam --gtf genes.gtf
```

## Required options

### `--gtf <GTF>` / `-g <GTF>`

Path to the GTF gene annotation file. Used to define gene models for read
counting and expression-level calculation. The GTF must contain `exon` features
with a `gene_id` attribute.

## Optional flags

### `-o, --outdir <DIR>`

Output directory for all result files. Created if it does not exist.

**Default:** `.` (current working directory)

### `-s, --stranded <0|1|2>`

Library strandedness for strand-aware read counting:

| Value | Meaning |
|-------|---------|
| `0` | Unstranded (count reads on either strand) |
| `1` | Forward stranded (read 1 maps to the transcript strand) |
| `2` | Reverse stranded (read 2 maps to the transcript strand) |

**Default:** `0` (unstranded)

### `-p, --paired`

Enable paired-end mode. When set, read pairs are counted as a single fragment.
Both mates must overlap the gene for the pair to be assigned.

**Default:** single-end mode

### `-t, --threads <N>`

Number of threads for parallel processing. Parallelism is applied across
chromosomes within each BAM file, so the effective speedup depends on the number
of chromosomes with mapped reads.

**Default:** `1`

### `-r, --reference <PATH>`

Path to a reference FASTA file. Required when using **CRAM** input files, as CRAM
files store sequences relative to a reference.

### `--skip-dup-check`

Skip the pre-flight validation that checks for duplicate-marking tool signatures
in the BAM header. By default, RustQC inspects `@PG` header lines for known
duplicate-marking tools (Picard MarkDuplicates, samblaster, sambamba, biobambam,
etc.) and exits with an error if none are found.

Use this flag if your BAM was marked by an unrecognized tool, or if you want to
run on a file without duplicate marking for testing purposes.

### `--biotype-attribute <NAME>`

GTF attribute name to use for biotype grouping in the featureCounts biotype
output files.

- Ensembl GTFs typically use `gene_biotype`
- GENCODE GTFs typically use `gene_type`

If not specified, RustQC defaults to `gene_biotype` and auto-detects the
attribute. If the specified attribute is not found in the GTF, a warning is
printed and biotype counting is skipped.

Set to an empty string to disable biotype counting entirely.

### `-c, --config <PATH>`

Path to a YAML configuration file for advanced settings such as chromosome name
mapping and per-output-file toggles. See the
[Configuration](/RustQC/usage/configuration/) page for the full reference.

## Examples

Basic paired-end analysis:

```bash
rustqc rna sample.bam --gtf genes.gtf -p -o results/
```

Reverse-stranded library with 8 threads:

```bash
rustqc rna sample.bam --gtf genes.gtf -p -s 2 -t 8 -o results/
```

CRAM input with reference:

```bash
rustqc rna sample.cram --gtf genes.gtf -p -r genome.fa -o results/
```

GENCODE GTF with explicit biotype attribute:

```bash
rustqc rna sample.bam --gtf gencode.v46.gtf -p \
  --biotype-attribute gene_type -o results/
```

Using a YAML config for chromosome name mapping:

```bash
rustqc rna sample.bam --gtf genes.gtf -p -c config.yaml -o results/
```

Multiple BAM files with parallel processing:

```bash
rustqc rna *.bam --gtf genes.gtf -p -t 8 -o results/
```

## Exit codes

| Code | Meaning |
|------|---------|
| `0` | Success |
| `1` | Error (missing input, invalid arguments, BAM processing failure, etc.) |

Error messages are printed to stderr with context about the failure.
