<h1 align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/public/RustQC-logo-darkbg.svg">
  <source media="(prefers-color-scheme: light)" srcset="docs/public/RustQC-logo.svg">
  <img width="500" src="docs/public/RustQC-logo.svg" alt="RustQC">
</picture>
</h1>

<h4 align="center">Fast genomics quality control tools for sequencing data, written in Rust.</h4>

---

**RustQC** is a suite of fast QC tools for sequencing data. Currently it includes the `rna` subcommand — a reimplementation of [dupRadar](https://github.com/ssayols/dupRadar) for assessing PCR duplicate rates in RNA-Seq datasets, with [featureCounts](http://subread.sourceforge.net/)-compatible output and biotype counting.

It analyzes duplicate-marked alignment files (SAM/BAM/CRAM) to compute per-gene duplication rates as a function of expression level. In a single pass, it produces the same outputs as the original **dupRadar** R/Bioconductor package, featureCounts-format count files, and biotype-level QC summaries — all significantly faster and compiled to a single static binary with no runtime dependencies.

<p align="center">
<picture>
   <source media="(prefers-color-scheme: dark)" srcset="docs/public/benchmarks/benchmark_dark.svg">
   <img src="docs/public/benchmarks/benchmark_light.svg" alt="Benchmark: RustQC 54s vs featureCounts 16m vs dupRadar 30m" width="600">
</picture>
</p>

<p align="center"><em>Run time for a 10 GB paired-end BAM</em></p>

## Comparison with dupRadar

| Feature | dupRadar (R) | RustQC |
|---------|-------------|--------|
| Language | R | Rust |
| Dependencies | R, Bioconductor, Rsubread | None (static binary) |
| Read counting | 4 separate featureCounts calls | Single-pass alignment reading |
| Speed | ~30 min for 10 GB BAM | <1 min for 10 GB BAM |
| Memory | High (R overhead) | Low |
| Output format | Identical | Identical |

Benchmark results:

| Metric | dupRadar (R) | RustQC |
| --- | --- | --- |
| **Runtime** | 29m 56s | 0m 54s (~33x faster) |
| **Intercept** | 0.8245 | 0.8245 |
| **Slope** | 1.6774 | 1.6774 |

All gene counts match **exactly** across all 63,086 genes. Cell-by-cell comparison of the full duplication matrix (820,118 values) shows **zero mismatches**.

See the [benchmark documentation](https://ewels.github.io/RustQC/benchmarks/dupradar/) for detailed results and the [benchmark directory](benchmark/) for replication instructions.

### Density scatter plots

<table>
<tr><th>dupRadar (R)</th><th>RustQC</th></tr>
<tr>
<td><img src="benchmark/large/dupRadar/duprateExpDens.png" width="400"></td>
<td><img src="benchmark/large/RustQC/GM12878_REP1.markdup.sorted_duprateExpDens.png" width="400"></td>
</tr>
</table>

### Boxplots

<table>
<tr><th>dupRadar (R)</th><th>RustQC</th></tr>
<tr>
<td><img src="benchmark/large/dupRadar/duprateExpBoxplot.png" width="400"></td>
<td><img src="benchmark/large/RustQC/GM12878_REP1.markdup.sorted_duprateExpBoxplot.png" width="400"></td>
</tr>
</table>

### Expression histograms

<table>
<tr><th>dupRadar (R)</th><th>RustQC</th></tr>
<tr>
<td><img src="benchmark/large/dupRadar/expressionHist.png" width="400"></td>
<td><img src="benchmark/large/RustQC/GM12878_REP1.markdup.sorted_expressionHist.png" width="400"></td>
</tr>
</table>

## Installation

### Pre-built binaries

Download a pre-built binary for your platform from the [Releases](https://github.com/ewels/RustQC/releases) page:

```bash
# Linux (x86_64)
curl -fsSL https://github.com/ewels/RustQC/releases/latest/download/rustqc-linux-x86_64.tar.gz | tar xz
chmod +x ./rustqc
sudo mv rustqc /usr/local/bin/

# Linux (aarch64)
curl -fsSL https://github.com/ewels/RustQC/releases/latest/download/rustqc-linux-aarch64.tar.gz | tar xz
chmod +x ./rustqc
sudo mv rustqc /usr/local/bin/

# macOS (Apple Silicon)
curl -fsSL https://github.com/ewels/RustQC/releases/latest/download/rustqc-macos-aarch64.tar.gz | tar xz
chmod +x ./rustqc
sudo mv rustqc /usr/local/bin/

# macOS (Intel)
curl -fsSL https://github.com/ewels/RustQC/releases/latest/download/rustqc-macos-x86_64.tar.gz | tar xz
chmod +x ./rustqc
sudo mv rustqc /usr/local/bin/
```

### Docker

```bash
docker run --rm -v "$PWD":/data ghcr.io/ewels/rustqc:latest \
  rna /data/sample.markdup.bam --gtf /data/genes.gtf --outdir /data/results
```

Available tags: `latest`, or a specific version (e.g., `0.1.0`).

### From source

Requires Rust toolchain and C build dependencies (see [CONTRIBUTING.md](CONTRIBUTING.md#prerequisites)).

```bash
cargo build --release
```

The binary will be at `target/release/rustqc`.

## Usage

```bash
rustqc rna <INPUT>... --gtf <GTF> [OPTIONS]
```

### Required arguments

| Argument | Description |
|----------|-------------|
| `<INPUT>...` | One or more duplicate-marked alignment files (SAM/BAM/CRAM). Duplicates must be flagged (SAM flag 0x400), not removed. BAM/CRAM files should be sorted and indexed for parallel processing. |
| `--gtf <GTF>` | Path to a GTF gene annotation file (e.g., from Ensembl or UCSC). |

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--stranded <0\|1\|2>` | `0` | Library strandedness: `0` = unstranded, `1` = stranded (forward), `2` = reverse-stranded |
| `--paired` | `false` | Set if the library is paired-end |
| `--threads <N>` | `1` | Number of threads for parallel alignment processing |
| `--outdir <DIR>` | `.` | Output directory |
| `--reference <FASTA>` / `-r` | none | Reference FASTA file (required for CRAM input) |
| `--config <FILE>` | none | Path to a YAML configuration file (see [Configuration](#configuration)) |
| `--biotype-attribute <NAME>` | none | Override the GTF attribute used for biotype counting (default from config: `gene_biotype`). Use `gene_type` for GENCODE GTFs. |
| `--skip-dup-check` | `false` | Skip verification that duplicates have been marked in the BAM file (see [Duplicate marking](#duplicate-marking)) |

### Duplicate marking

RustQC requires that the input BAM file has been processed by a duplicate-marking tool such as [Picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates), [samblaster](https://github.com/GregoryFaust/samblaster), or [sambamba markdup](https://lomereiter.github.io/sambamba/). These tools set the SAM flag `0x400` on PCR/optical duplicate reads, which RustQC uses to compute duplication rates.

Before processing, RustQC checks the BAM `@PG` header lines for known duplicate-marking programs. If none are found, it exits with an error explaining how to mark duplicates. As a secondary safeguard, if processing completes but zero duplicate-flagged reads are found among mapped reads, RustQC also exits with an error.

If you are confident that your BAM file has duplicates correctly flagged despite the header check failing, you can bypass the verification with `--skip-dup-check`.

When multiple BAM files are provided, they are processed in parallel using the available threads. The GTF is parsed once and shared across all BAM files. Threads are divided evenly among the parallel BAM jobs.

### Example

```bash
# Single BAM, single-end, unstranded
rustqc rna sample.markdup.bam --gtf genes.gtf --outdir results/

# Single BAM, paired-end, reverse-stranded
rustqc rna sample.markdup.bam --gtf genes.gtf --paired --stranded 2 --outdir results/

# CRAM input with reference FASTA
rustqc rna sample.markdup.cram --gtf genes.gtf --reference genome.fa --outdir results/

# SAM input (single-threaded only, no index)
rustqc rna sample.markdup.sam --gtf genes.gtf --outdir results/

# Multiple alignment files in parallel
rustqc rna sample1.bam sample2.bam sample3.bam --gtf genes.gtf --paired --threads 12 --outdir results/

# With chromosome name mapping (e.g. Ensembl alignment + UCSC GTF)
rustqc rna sample.markdup.bam --gtf genes.gtf --paired --config config.yaml --outdir results/

# GENCODE GTF with gene_type attribute for biotypes
rustqc rna sample.markdup.bam --gtf gencode.gtf --paired --biotype-attribute gene_type --outdir results/
```

## Configuration

An optional YAML configuration file can be provided with `--config` to control runtime behaviour. The file is designed to be extensible — unknown fields are silently ignored, so config files remain forward-compatible.

### Chromosome name mapping

When the alignment and GTF files use different chromosome naming conventions (e.g. Ensembl `1, 2, X` vs. UCSC `chr1, chr2, chrX`), RustQC will detect the mismatch and exit with a helpful error. You can resolve this with either a prefix or explicit mapping.

**Prefix** — prepend a string to every alignment chromosome name before matching:

```yaml
chromosome_prefix: "chr"
```

**Explicit mapping** — for fine-grained control, map individual GTF names to alignment file names:

```yaml
chromosome_mapping:
  chr1: "1"
  chr2: "2"
  chrX: "X"
  chrM: "MT"
```

Both options can be combined. The prefix is applied first, then explicit mappings override specific names.

### Output control

Output files are grouped by tool and can be individually enabled or disabled. All outputs are enabled by default.

```yaml
dupradar:
  dup_matrix: true
  intercept_slope: true
  density_scatter_plot: true
  boxplot: true
  expression_histogram: true
  multiqc_intercept: true
  multiqc_curve: true

featurecounts:
  counts_file: true            # featureCounts-format gene counts
  summary_file: true           # featureCounts-format assignment summary
  biotype_counts: true         # per-biotype count table
  biotype_counts_mqc: true     # MultiQC biotype bargraph
  biotype_rrna_mqc: true       # MultiQC rRNA percentage
  biotype_attribute: "gene_biotype"  # GTF attribute for biotype grouping
```

The `biotype_attribute` controls which GTF attribute is used to group genes into biotypes. Ensembl GTFs use `gene_biotype` (the default), while GENCODE GTFs use `gene_type`. This can also be overridden on the command line with `--biotype-attribute`.

If the configured biotype attribute is not found in the GTF, biotype outputs are skipped with a warning.

### Example config file

```yaml
# Prepend "chr" to all alignment chromosome names
chromosome_prefix: "chr"

# Override the mitochondrial chromosome mapping
# (prefix would produce "chrMT", but GTF uses "chrM")
chromosome_mapping:
  chrM: "MT"

# Disable dupRadar plots, keep only the matrix
dupradar:
  density_scatter_plot: false
  boxplot: false
  expression_histogram: false

# Use GENCODE biotype attribute
featurecounts:
  biotype_attribute: "gene_type"
```

## Output files

For an input file named `sample.bam` (or `sample.cram`, etc.), the following files are generated. All outputs can be individually enabled or disabled via the [configuration file](#output-control).

### dupRadar outputs

| File | Description |
|------|-------------|
| `sample_dupMatrix.txt` | Tab-separated duplication matrix (14 columns, one row per gene) |
| `sample_duprateExpDens.{png,svg}` | Density scatter plot of duplication rate vs. expression |
| `sample_duprateExpBoxplot.{png,svg}` | Boxplot of duplication rate by expression quantile bins |
| `sample_expressionHist.{png,svg}` | Histogram of gene expression levels (log10 RPK) |
| `sample_intercept_slope.txt` | Logistic regression fit parameters (intercept and slope) |
| `sample_dup_intercept_mqc.txt` | MultiQC general stats format with intercept value |
| `sample_duprateExpDensCurve_mqc.txt` | MultiQC line graph data for the fitted curve |

### featureCounts outputs

| File | Description |
|------|-------------|
| `sample.featureCounts.tsv` | Gene-level read counts in featureCounts format (7 columns: Geneid, Chr, Start, End, Strand, Length, counts) |
| `sample.featureCounts.tsv.summary` | Assignment summary statistics (Assigned, Unassigned_NoFeatures, Unassigned_Ambiguity, etc.) |
| `sample.biotype_counts.tsv` | Read counts aggregated by biotype (e.g., protein_coding, lncRNA, rRNA) |
| `sample.biotype_counts_mqc.tsv` | MultiQC bargraph format for biotype counts |
| `sample.biotype_counts_rrna_mqc.tsv` | MultiQC general stats with rRNA percentage |

The featureCounts output files use the same format as [Subread featureCounts](http://subread.sourceforge.net/), making them compatible with downstream tools that expect that format. The biotype outputs replicate the [nf-core/rnaseq](https://nf-co.re/rnaseq) biotype QC functionality.

### Duplication matrix columns

| Column | Description |
|--------|-------------|
| `ID` | Gene identifier |
| `geneLength` | Effective gene length (non-overlapping exon bases) |
| `allCountsMulti` | Total read count (including multimappers and duplicates) |
| `filteredCountsMulti` | Read count excluding duplicates (including multimappers) |
| `dupRateMulti` | Duplication rate with multimappers |
| `dupsPerIdMulti` | Number of duplicate reads with multimappers |
| `RPKMulti` | Reads per kilobase with multimappers |
| `RPKMMulti` | RPKM with multimappers |
| `allCounts` | Total read count (unique mappers only) |
| `filteredCounts` | Read count excluding duplicates (unique mappers only) |
| `dupRate` | Duplication rate (unique mappers only) |
| `dupsPerId` | Number of duplicate reads (unique mappers only) |
| `RPK` | Reads per kilobase (unique mappers only) |
| `RPKM` | RPKM (unique mappers only) |

## Performance tuning

RustQC uses multi-threaded alignment processing when `--threads` is set above 1. Within a single file, chromosomes are distributed across threads and processed in parallel, typically achieving near-linear speedup. When multiple alignment files are provided, they are also processed in parallel — the available threads are divided evenly among concurrent jobs. For a single sample, `--threads 4` is a good starting point. For multiple samples, use enough threads to keep all jobs busy (e.g., `--threads 12` for 3 files gives each 4 threads). Multi-threading requires an indexed file (`.bai`/`.csi` for BAM, `.crai` for CRAM). SAM files are always processed single-threaded.

For maximum performance when building from source, you can enable CPU-specific optimizations:

```bash
# Build with native CPU instruction set (AVX2, etc.)
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

For an additional 5-20% speedup on frequently-used machines, Profile-Guided Optimization (PGO) can be used:

```bash
# Step 1: Build with profiling instrumentation
RUSTFLAGS="-Cprofile-generate=/tmp/pgo-data" cargo build --release

# Step 2: Run on representative data to collect profiles
target/release/rustqc rna sample.bam --gtf genes.gtf --paired --threads 4 -o /tmp/pgo-run

# Step 3: Merge profile data
llvm-profdata merge -o /tmp/pgo-data/merged.profdata /tmp/pgo-data

# Step 4: Rebuild with profile data
RUSTFLAGS="-Cprofile-use=/tmp/pgo-data/merged.profdata -Cllvm-args=-pgo-warn-missing-function" \
  cargo build --release
```

Note: PGO profiles are machine-specific and `target-cpu=native` produces non-portable binaries. Pre-built release binaries use generic optimizations that work on all machines.

## How it works

1. **GTF parsing**: Reads gene annotations, computes effective gene lengths from non-overlapping exon bases, and extracts additional attributes (e.g., biotype) for downstream grouping.
2. **Read counting**: Reads the alignment file (SAM/BAM/CRAM) once, assigning each read to a gene based on exon overlap. Four count modes are tracked simultaneously:
   - With/without multimappers x with/without duplicates
   - Assignment statistics (assigned, ambiguous, no features) are tracked for the featureCounts summary.
3. **featureCounts output**: Writes gene-level counts and summary statistics in the standard featureCounts format.
4. **Biotype counting**: Aggregates assigned read counts by a configurable GTF attribute (e.g., `gene_biotype`), producing biotype count tables and MultiQC-compatible rRNA QC metrics.
5. **Duplication matrix**: Computes RPK, RPKM, and duplication rates for each gene in all four modes.
6. **Logistic regression**: Fits a binomial GLM (`dupRate ~ log10(RPK)`) using iteratively reweighted least squares (IRLS) to model the relationship between expression and duplication.
7. **Plots**: Generates density scatter, boxplot, and histogram visualizations.
8. **MultiQC integration**: Outputs files compatible with [MultiQC](https://multiqc.info/) for pipeline reporting (dupRadar intercept, fit curve, biotype bargraph, rRNA percentage).

## Interpreting the results

- **Intercept** (exp(beta0)): Indicates duplication rate at low expression. Low values = good quality. High values = PCR artifact problems.
- **Slope** (exp(beta1)): Rate at which duplication increases with expression. Single-end libraries typically have higher slope than paired-end.
- **Density plot**: Good samples show low duplication (bottom of y-axis) at low expression (left), with duplication rising naturally only at very high expression.
- **1 read/bp threshold** (red dashed line): At RPK=1000, a 1kb gene has 1000 reads, meaning roughly 1 read per base pair -- near the theoretical maximum for unique reads.

## References

- Sayols S, Scherzinger D, Klein H (2016). dupRadar: a Bioconductor package for the assessment of PCR artifacts in RNA-Seq data. *BMC Bioinformatics*, 17, 428. doi:10.1186/s12859-016-1276-2
- Original R package: https://github.com/ssayols/dupRadar

## License

MIT License. See [LICENSE](LICENSE) for details.
