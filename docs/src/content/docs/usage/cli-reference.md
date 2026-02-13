---
title: CLI Reference
description: Complete command-line interface reference for RustQC, including all subcommands, options, flags, and usage examples.
---

RustQC uses a subcommand-based CLI. Each subcommand provides a specific QC
analysis tool, all operating on SAM/BAM/CRAM alignment files.

## Available subcommands

| Subcommand | Description | Equivalent to |
|------------|-------------|---------------|
| [`rna`](#rna) | RNA-seq duplicate rate analysis with integrated read counting | dupRadar + featureCounts |
| [`bam-stat`](#bam-stat) | Basic BAM alignment statistics | RSeQC `bam_stat.py` |
| [`infer-experiment`](#infer-experiment) | Infer library strandedness | RSeQC `infer_experiment.py` |
| [`read-duplication`](#read-duplication) | Position- and sequence-based duplication histograms | RSeQC `read_duplication.py` |
| [`read-distribution`](#read-distribution) | Read distribution across genomic features | RSeQC `read_distribution.py` |
| [`junction-annotation`](#junction-annotation) | Classify splice junctions as known/novel | RSeQC `junction_annotation.py` |
| [`junction-saturation`](#junction-saturation) | Saturation analysis for splice junctions | RSeQC `junction_saturation.py` |
| [`inner-distance`](#inner-distance) | Inner distance between paired-end reads | RSeQC `inner_distance.py` |

---

## `rna`

RNA-seq duplicate rate analysis (dupRadar equivalent) with integrated
featureCounts-compatible output and biotype counting.

### Synopsis

```
rustqc rna <INPUT>... --gtf <GTF> [OPTIONS]
```

### Positional arguments

#### `<INPUT>...`

One or more paths to duplicate-marked alignment files. Accepted formats are
**SAM**, **BAM**, and **CRAM**.

Multiple files can be passed and will be processed in parallel, with each file
producing its own set of output files. Threads are divided evenly among
concurrent jobs.

```bash
# Single file
rustqc rna sample.bam --gtf genes.gtf

# Multiple files
rustqc rna sample1.bam sample2.bam sample3.bam --gtf genes.gtf
```

### Required options

#### `--gtf <GTF>` / `-g <GTF>`

Path to the GTF gene annotation file. Used to define gene models for read
counting and expression-level calculation. The GTF must contain `exon` features
with a `gene_id` attribute.

### Optional flags

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

#### `-c, --config <PATH>`

Path to a YAML configuration file for advanced settings such as chromosome name
mapping and per-output-file toggles. See the
[Configuration](/usage/configuration/) page for the full reference.

### Examples

```bash
# Basic paired-end analysis
rustqc rna sample.bam --gtf genes.gtf -p -o results/

# Reverse-stranded library with 8 threads
rustqc rna sample.bam --gtf genes.gtf -p -s 2 -t 8 -o results/

# CRAM input with reference
rustqc rna sample.cram --gtf genes.gtf -p -r genome.fa -o results/

# GENCODE GTF with explicit biotype attribute
rustqc rna sample.bam --gtf gencode.v46.gtf -p \
  --biotype-attribute gene_type -o results/

# Using a YAML config for chromosome name mapping
rustqc rna sample.bam --gtf genes.gtf -p -c config.yaml -o results/

# Multiple BAM files with parallel processing
rustqc rna *.bam --gtf genes.gtf -p -t 8 -o results/
```

---

## `bam-stat`

Basic BAM alignment statistics. Reimplements [RSeQC](https://rseqc.sourceforge.net/)
`bam_stat.py`: a single-pass BAM reader that collects fundamental alignment
metrics (total records, QC-failed, duplicates, mapping quality distribution,
splice reads, proper pairs, etc.).

### Synopsis

```
rustqc bam-stat <INPUT>... [OPTIONS]
```

### Positional arguments

#### `<INPUT>...`

One or more SAM/BAM/CRAM files.

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `-q, --mapq <N>` | `30` | Minimum MAPQ for the "uniquely mapped" count |
| `-o, --outdir <DIR>` | `.` | Output directory |
| `-r, --reference <PATH>` | -- | Reference FASTA (required for CRAM) |

### Example

```bash
rustqc bam-stat sample.bam -o results/
```

### Output

`{sample}.bam_stat.txt` -- A formatted text report matching the RSeQC
`bam_stat.py` output format. See [RSeQC Outputs](/outputs/rseqc/) for
details.

---

## `infer-experiment`

Infer library strandedness from RNA-seq alignments. Reimplements RSeQC
`infer_experiment.py`. Samples reads that overlap gene models and determines the
fraction consistent with each strand protocol.

### Synopsis

```
rustqc infer-experiment <INPUT>... -b <BED> [OPTIONS]
```

### Positional arguments

#### `<INPUT>...`

One or more SAM/BAM/CRAM files.

### Required options

#### `-b, --bed <BED>`

BED12-format gene model file. Each line defines a transcript with exon block
information. Used to determine gene strand for strandedness inference.

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `-q, --mapq <N>` | `30` | Minimum MAPQ for reads to sample |
| `-s, --sample-size <N>` | `200000` | Maximum number of reads to sample |
| `-o, --outdir <DIR>` | `.` | Output directory |
| `-r, --reference <PATH>` | -- | Reference FASTA (required for CRAM) |

### Example

```bash
rustqc infer-experiment sample.bam -b genes.bed -o results/
```

### Output

`{sample}.infer_experiment.txt` -- Strandedness fractions (fraction failing to
determine, fraction consistent with forward/reverse protocols). See
[RSeQC Outputs](/outputs/rseqc/) for interpretation guidance.

---

## `read-duplication`

Position-based and sequence-based duplication histograms. Reimplements RSeQC
`read_duplication.py`. Computes how many reads map to the same position or share
identical sequences, producing occurrence histograms.

### Synopsis

```
rustqc read-duplication <INPUT>... [OPTIONS]
```

### Positional arguments

#### `<INPUT>...`

One or more SAM/BAM/CRAM files.

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `-q, --mapq <N>` | `30` | Minimum MAPQ |
| `-o, --outdir <DIR>` | `.` | Output directory |
| `-r, --reference <PATH>` | -- | Reference FASTA (required for CRAM) |

### Example

```bash
rustqc read-duplication sample.bam -o results/
```

### Output

- `{sample}.pos.DupRate.xls` -- Position-based duplication histogram
- `{sample}.seq.DupRate.xls` -- Sequence-based duplication histogram

See [RSeQC Outputs](/outputs/rseqc/) for details.

---

## `read-distribution`

Read distribution across genomic features. Reimplements RSeQC
`read_distribution.py`. Classifies BAM read tags into CDS exons, UTRs, introns,
and intergenic regions using a BED12 gene model.

### Synopsis

```
rustqc read-distribution <INPUT>... -b <BED> [OPTIONS]
```

### Positional arguments

#### `<INPUT>...`

One or more SAM/BAM/CRAM files.

### Required options

#### `-b, --bed <BED>`

BED12-format gene model file.

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --outdir <DIR>` | `.` | Output directory |
| `-r, --reference <PATH>` | -- | Reference FASTA (required for CRAM) |

### Example

```bash
rustqc read-distribution sample.bam -b genes.bed -o results/
```

### Output

`{sample}.read_distribution.txt` -- Tabular report with total reads, total
tags, and per-region breakdown (CDS Exons, 5'UTR, 3'UTR, Introns, TSS/TES
flanking regions). See [RSeQC Outputs](/outputs/rseqc/) for details.

---

## `junction-annotation`

Splice junction annotation and classification. Reimplements RSeQC
`junction_annotation.py`. Extracts splice junctions from CIGAR N-operations and
classifies them as known (annotated), partial novel, or complete novel by
comparing against a BED12 gene model.

### Synopsis

```
rustqc junction-annotation <INPUT>... -b <BED> [OPTIONS]
```

### Positional arguments

#### `<INPUT>...`

One or more SAM/BAM/CRAM files.

### Required options

#### `-b, --bed <BED>`

BED12-format gene model file providing reference splice site annotations.

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `-m, --min-intron <N>` | `50` | Minimum intron size to consider |
| `-q, --mapq <N>` | `30` | Minimum MAPQ |
| `-o, --outdir <DIR>` | `.` | Output directory |
| `-r, --reference <PATH>` | -- | Reference FASTA (required for CRAM) |

### Example

```bash
rustqc junction-annotation sample.bam -b genes.bed -o results/
```

### Output

- `{sample}.junction.xls` -- TSV of all junctions with annotation status
- `{sample}.junction.bed` -- BED12 with colour-coded junctions (red=known, green=partial novel, blue=novel)
- `{sample}.junction_plot.r` -- R script for pie chart visualisation
- `{sample}.junction_annotation.txt` -- Summary counts

See [RSeQC Outputs](/outputs/rseqc/) for details.

---

## `junction-saturation`

Saturation analysis for splice junctions. Reimplements RSeQC
`junction_saturation.py`. Subsamples splice junctions at increasing percentages
of total reads and reports how many known/novel/total unique junctions are
detected at each level.

### Synopsis

```
rustqc junction-saturation <INPUT>... -b <BED> [OPTIONS]
```

### Positional arguments

#### `<INPUT>...`

One or more SAM/BAM/CRAM files.

### Required options

#### `-b, --bed <BED>`

BED12-format gene model file providing reference splice site annotations.

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `-m, --min-intron <N>` | `50` | Minimum intron size |
| `-v, --min-coverage <N>` | `1` | Minimum reads for a junction to count as known |
| `-l, --percentile-floor <N>` | `5` | Sampling start percentage |
| `-u, --percentile-ceiling <N>` | `100` | Sampling end percentage |
| `-s, --percentile-step <N>` | `5` | Sampling step percentage |
| `-q, --mapq <N>` | `30` | Minimum MAPQ |
| `-o, --outdir <DIR>` | `.` | Output directory |
| `-r, --reference <PATH>` | -- | Reference FASTA (required for CRAM) |

### Example

```bash
rustqc junction-saturation sample.bam -b genes.bed -o results/

# Custom sampling range
rustqc junction-saturation sample.bam -b genes.bed -l 10 -u 100 -s 10
```

### Output

- `{sample}.junctionSaturation_plot.r` -- R script for saturation curve
- `{sample}.junctionSaturation_summary.txt` -- Tabular data (percentage vs junction counts)

See [RSeQC Outputs](/outputs/rseqc/) for details.

---

## `inner-distance`

Inner distance between paired-end reads. Reimplements RSeQC
`inner_distance.py`. Computes the inner distance between read pairs, classifying
each pair by gene model overlap and producing a histogram.

### Synopsis

```
rustqc inner-distance <INPUT>... -b <BED> [OPTIONS]
```

### Positional arguments

#### `<INPUT>...`

One or more SAM/BAM/CRAM files. Must be paired-end data.

### Required options

#### `-b, --bed <BED>`

BED12-format gene model file.

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `-k, --sample-size <N>` | `1000000` | Maximum read pairs to sample |
| `-l, --lower-bound <N>` | `-250` | Histogram lower bound |
| `-u, --upper-bound <N>` | `250` | Histogram upper bound |
| `-s, --step <N>` | `5` | Histogram bin width |
| `-q, --mapq <N>` | `30` | Minimum MAPQ |
| `-o, --outdir <DIR>` | `.` | Output directory |
| `-r, --reference <PATH>` | -- | Reference FASTA (required for CRAM) |

### Example

```bash
rustqc inner-distance sample.bam -b genes.bed -o results/

# Wider histogram range with larger bins
rustqc inner-distance sample.bam -b genes.bed -l -500 -u 500 -s 10
```

### Output

- `{sample}.inner_distance.txt` -- Per-pair detail (read name, distance, classification)
- `{sample}.inner_distance_freq.txt` -- Histogram (bin start, bin end, count)
- `{sample}.inner_distance_plot.r` -- R script for histogram and density plot
- `{sample}.inner_distance_summary.txt` -- Summary counts by classification

See [RSeQC Outputs](/outputs/rseqc/) for details.

---

## Common options

Several options appear across most subcommands:

| Option | Description | Subcommands |
|--------|-------------|-------------|
| `-o, --outdir` | Output directory | All |
| `-r, --reference` | Reference FASTA for CRAM | All |
| `-q, --mapq` | MAPQ cutoff | All except `rna`, `read-distribution` |
| `-b, --bed` | BED12 gene model | All RSeQC tools except `bam-stat`, `read-duplication` |

## Exit codes

| Code | Meaning |
|------|---------|
| `0` | Success |
| `1` | Error (missing input, invalid arguments, BAM processing failure, etc.) |

Error messages are printed to stderr with context about the failure.
