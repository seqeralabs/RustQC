---
title: CLI Reference
description: Complete command-line interface reference for RustQC, including all options, flags, and usage examples.
---

RustQC provides a single `rna` subcommand that runs all RNA-Seq QC analyses in
one pass over the BAM file.

## `rna`

RNA-seq quality control: duplicate rate analysis (dupRadar equivalent),
featureCounts-compatible read counting with biotype summaries, and 7
RSeQC-equivalent tools (bam_stat, infer_experiment, read_duplication,
read_distribution, junction_annotation, junction_saturation, inner_distance).

### Synopsis

```
rustqc rna <INPUT>... (--gtf <GTF> | --bed <BED>) [OPTIONS]
```

### Positional arguments

#### `<INPUT>...`

One or more paths to duplicate-marked alignment files. Accepted formats are
**SAM**, **BAM**, and **CRAM**.

Multiple files can be passed and will be processed in parallel, with each file
producing its own set of output files. Threads are divided evenly among
concurrent jobs.

```bash
# Single file with GTF (all analyses)
rustqc rna sample.bam --gtf genes.gtf

# Single file with BED (RSeQC tools only; dupRadar/featureCounts skipped)
rustqc rna sample.bam --bed genes.bed

# Multiple files
rustqc rna sample1.bam sample2.bam sample3.bam --gtf genes.gtf
```

### Annotation options

Exactly one of `--gtf` or `--bed` must be provided. They are **mutually exclusive**.

#### `--gtf <GTF>` / `-g <GTF>`

Path to a GTF gene annotation file (plain or gzip-compressed). The GTF must
contain `exon` features with a `gene_id` attribute. When a GTF is provided,
**all analyses run**: dupRadar, featureCounts, and all 7 RSeQC tools.
Transcript-level structure (exon blocks, CDS features) is extracted automatically
and used by the RSeQC tools that previously required a separate BED file.

Gzip compression is detected automatically by inspecting the file header (magic
bytes), so the `.gz` extension is not required.

#### `-b, --bed <BED>`

Path to a BED12-format gene model file (plain or gzip-compressed). When a BED
file is provided **without** a GTF, only the 7 RSeQC tools run (plus bam_stat
and read_duplication). dupRadar and featureCounts are skipped because BED files
lack gene-level grouping and biotype information.

Gzip compression is detected automatically by inspecting the file header (magic
bytes), so the `.gz` extension is not required.

Individual tools can be disabled via the [configuration file](/usage/configuration/).

#### `-o, --outdir <DIR>`

Output directory for all result files. Created if it does not exist.

**Default:** `.` (current working directory)

#### `-s, --stranded <0|1|2>`

Library strandedness for strand-aware read counting:

| Value | Meaning |
|-------|---------|
| `0` | Unstranded (count reads on either strand) |
| `1` | Forward stranded (read 1 maps to the transcript strand) |
| `2` | Reverse stranded (read 2 maps to the transcript strand) |

**Default:** `0` (unstranded)

#### `-p, --paired`

Enable paired-end mode. When set, read pairs are counted as a single fragment.
Both mates must overlap the gene for the pair to be assigned.

**Default:** single-end mode

#### `-t, --threads <N>`

Number of threads for parallel processing. Parallelism is applied across
chromosomes within each BAM file, so the effective speedup depends on the number
of chromosomes with mapped reads.

**Default:** `1`

#### `-q, --mapq <N>`

Minimum mapping quality (MAPQ) threshold used by all RSeQC tools. Reads below
this threshold are excluded from the "uniquely mapped" counts.

**Default:** `30`

#### `-r, --reference <PATH>`

Path to a reference FASTA file. Required when using **CRAM** input files, as CRAM
files store sequences relative to a reference.

#### `--skip-dup-check`

Skip the pre-flight validation that checks for duplicate-marking tool signatures
in the BAM header. By default, RustQC inspects `@PG` header lines for known
duplicate-marking tools (Picard MarkDuplicates, samblaster, sambamba, biobambam,
etc.) and exits with an error if none are found.

Use this flag if your BAM was marked by an unrecognized tool, or if you want to
run on a file without duplicate marking for testing purposes.

#### `--biotype-attribute <NAME>`

GTF attribute name to use for biotype grouping in the featureCounts biotype
output files.

- Ensembl GTFs typically use `gene_biotype`
- GENCODE GTFs typically use `gene_type`

If not specified, RustQC defaults to `gene_biotype` and auto-detects the
attribute. If the specified attribute is not found in the GTF, a warning is
printed and biotype counting is skipped.

#### `--flat-output`

Write all output files directly into the output directory instead of organizing
them into subdirectories. By default, RustQC creates `dupradar/`,
`featurecounts/`, and `rseqc/<tool>/` subdirectories under the output directory.
With `--flat-output`, all files are written to the top-level output directory.

This can also be set in the [configuration file](/usage/configuration/) as
`flat_output: true`.

**Default:** `false` (nested subdirectories)

#### `-c, --config <PATH>`

Path to a YAML configuration file for advanced settings such as chromosome name
mapping, per-output-file toggles, and enabling/disabling individual tools. See
the [Configuration](/usage/configuration/) page for the full reference.

### RSeQC tool options

These flags control parameters for specific RSeQC-equivalent analyses. Each
tool runs by default as part of `rustqc rna` and can be disabled via the
[configuration file](/usage/configuration/).

#### infer_experiment

| Option | Default | Description |
|--------|---------|-------------|
| `--infer-experiment-sample-size <N>` | `200000` | Maximum number of reads to sample |

#### junction_annotation / junction_saturation

| Option | Default | Description |
|--------|---------|-------------|
| `--min-intron <N>` | `50` | Minimum intron size to consider (shared by both tools) |

#### junction_saturation

| Option | Default | Description |
|--------|---------|-------------|
| `--junction-saturation-min-coverage <N>` | `1` | Minimum reads for a junction to count as detected |
| `--junction-saturation-percentile-floor <N>` | `5` | Sampling start percentage |
| `--junction-saturation-percentile-ceiling <N>` | `100` | Sampling end percentage |
| `--junction-saturation-percentile-step <N>` | `5` | Sampling step percentage |

#### inner_distance

| Option | Default | Description |
|--------|---------|-------------|
| `--inner-distance-sample-size <N>` | `1000000` | Maximum read pairs to sample |
| `--inner-distance-lower-bound <N>` | `-250` | Histogram lower bound |
| `--inner-distance-upper-bound <N>` | `250` | Histogram upper bound |
| `--inner-distance-step <N>` | `5` | Histogram bin width |

### Preseq options

These flags control parameters for the preseq library complexity extrapolation.
Preseq runs by default and can be skipped entirely with `--skip-preseq`.

| Option | Default | Description |
|--------|---------|-------------|
| `--skip-preseq` | `false` | Skip the preseq library complexity analysis entirely |
| `--preseq-max-extrap <N>` | `1e10` | Maximum extrapolation depth in total reads |
| `--preseq-step-size <N>` | `1e6` | Step size between extrapolation points |
| `--preseq-n-bootstraps <N>` | `100` | Number of bootstrap replicates for confidence intervals |

### Examples

```bash
# Basic paired-end analysis with GTF (all tools)
rustqc rna sample.bam --gtf genes.gtf -p -o results/

# BED-only mode (RSeQC tools only; dupRadar/featureCounts skipped)
rustqc rna sample.bam --bed genes.bed -p -o results/

# Reverse-stranded library with 8 threads
rustqc rna sample.bam --gtf genes.gtf -p -s 2 -t 8 -o results/

# CRAM input with reference
rustqc rna sample.cram --gtf genes.gtf -p -r genome.fa -o results/

# GENCODE GTF with explicit biotype attribute
rustqc rna sample.bam --gtf gencode.v46.gtf -p \
  --biotype-attribute gene_type -o results/

# Custom junction saturation sampling range
rustqc rna sample.bam --gtf genes.gtf -p \
  --junction-saturation-percentile-floor 10 \
  --junction-saturation-percentile-ceiling 100 \
  --junction-saturation-percentile-step 10 -o results/

# Custom inner distance histogram bounds
rustqc rna sample.bam --gtf genes.gtf -p \
  --inner-distance-lower-bound -500 --inner-distance-upper-bound 500 \
  --inner-distance-step 10 -o results/

# Multiple BAM files with parallel processing
rustqc rna *.bam --gtf genes.gtf -p -t 8 -o results/

# Gzip-compressed annotation files (auto-detected)
rustqc rna sample.bam --gtf genes.gtf.gz -p -o results/
rustqc rna sample.bam --bed genes.bed.gz -p -o results/
```

---

## Exit codes

| Code | Meaning |
|------|---------|
| `0` | Success |
| `1` | Error (missing input, invalid arguments, BAM processing failure, etc.) |

Error messages are printed to stderr with context about the failure.
