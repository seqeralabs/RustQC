<h1 align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/public/RustQC-logo-darkbg.svg">
  <source media="(prefers-color-scheme: light)" srcset="docs/public/RustQC-logo.svg">
  <img width="500" src="docs/public/RustQC-logo.svg" alt="RustQC">
</picture>
</h1>

<h4 align="center">Fast genomics quality control tools for sequencing data, written in Rust.</h4>

---

**RustQC** is a growing suite of fast QC tools for sequencing data, compiled to a single static binary with no runtime dependencies. It currently includes:

- **`rustqc rna`** -- A single-command RNA-Seq QC pipeline that runs all analyses in one pass: [dupRadar](https://github.com/ssayols/dupRadar) duplicate rate analysis, [featureCounts](http://subread.sourceforge.net/)-compatible read counting with biotype summaries, 7 [RSeQC](https://rseqc.sourceforge.net/) quality control tools (bam_stat, infer_experiment, read_duplication, read_distribution, junction_annotation, junction_saturation, inner_distance), [preseq](http://smithlabresearch.org/software/preseq/) library complexity estimation, [TIN](https://rseqc.sourceforge.net/#tin-py) transcript integrity analysis, gene body coverage profiling, and [samtools](http://www.htslib.org/)-compatible flagstat/idxstats/stats output.

All tools accept SAM/BAM/CRAM input and support processing multiple files in a single command. Annotation files (GTF and BED) can be provided plain or gzip-compressed (`.gz`).

<p align="center">
<picture>
   <source media="(prefers-color-scheme: dark)" srcset="docs/public/benchmarks/benchmark_dark.svg">
   <img src="docs/public/benchmarks/benchmark_light.svg" alt="Benchmark: RustQC ~5m vs traditional tools ~2h 45m (dupRadar + featureCounts + 9 RSeQC tools + preseq + samtools + Qualimap)" width="600">
</picture>
</p>

<p align="center"><em>Run time for a 10 GB paired-end BAM (dupRadar + featureCounts + RSeQC + preseq + samtools + Qualimap)</em></p>

## Available tools

### `rustqc rna` -- RNA-Seq quality control

Performs dupRadar-equivalent duplicate rate analysis, featureCounts-compatible read counting, 7 RSeQC-equivalent QC analyses, preseq library complexity estimation, TIN transcript integrity analysis, gene body coverage profiling, and samtools-compatible summary statistics -- all in a single pass. Given a duplicate-marked BAM and either a GTF annotation or a BED12 gene model, it computes per-gene duplication rates, fits a logistic regression model, generates diagnostic plots, produces gene-level count files with biotype summaries, extrapolates library complexity, and runs comprehensive QC metrics.

| Feature | dupRadar (R) | RustQC |
|---------|-------------|--------|
| Language | R | Rust |
| Dependencies | R, Bioconductor, Rsubread | None (static binary) |
| Read counting | 4 separate featureCounts calls | Single-pass alignment reading |
| Speed | ~30 min for 10 GB BAM | ~4 min for 10 GB BAM (all tools) |
| Memory | High (R overhead) | Low |
| Output format | Identical | Identical |

All gene counts match **exactly** across all 63,086 genes. Cell-by-cell comparison of the full duplication matrix (820,118 values) shows **zero mismatches**. See the [benchmark documentation](https://ewels.github.io/RustQC/benchmarks/dupradar/) for detailed results.

### RSeQC tools (integrated)

The `rustqc rna` command also includes reimplementations of 7 popular [RSeQC](https://rseqc.sourceforge.net/) Python QC tools, all running in the same pipeline:

| Tool | RSeQC equivalent | Description |
|------|-----------------|-------------|
| bam_stat | `bam_stat.py` | Basic alignment statistics (total reads, duplicates, mapping quality, splice reads, etc.) |
| infer_experiment | `infer_experiment.py` | Infer library strandedness from read/gene-model overlap |
| read_duplication | `read_duplication.py` | Position-based and sequence-based duplication histograms |
| read_distribution | `read_distribution.py` | Read distribution across genomic features (CDS, UTR, intron, intergenic) |
| junction_annotation | `junction_annotation.py` | Classify splice junctions as known, partial novel, or complete novel |
| junction_saturation | `junction_saturation.py` | Saturation analysis of detected splice junctions |
| inner_distance | `inner_distance.py` | Inner distance distribution for paired-end reads |

All RSeQC tools run by default when annotation is provided via `--gtf` or `--bed`. With a GTF file, all tools run alongside dupRadar and featureCounts. With a BED file only, the RSeQC tools, TIN, and preseq run but dupRadar, featureCounts, and gene body coverage are skipped (they require a GTF). Individual tools can be disabled via the [configuration file](#configuration).

### Additional tools (integrated)

| Tool | Upstream equivalent | Description |
|------|-------------------|-------------|
| preseq | [`preseq lc_extrap`](http://smithlabresearch.org/software/preseq/) | Library complexity extrapolation -- estimates distinct molecules at increasing sequencing depths |
| TIN | [`tin.py`](https://rseqc.sourceforge.net/#tin-py) | Transcript Integrity Number -- measures transcript degradation via coverage uniformity |
| Gene body coverage | [Qualimap](http://qualimap.conesalab.org/) rnaseq | Coverage profile along gene bodies (5'→3'), Qualimap-compatible output for MultiQC |
| flagstat | [`samtools flagstat`](http://www.htslib.org/) | Alignment flag summary statistics, identical format to samtools |
| idxstats | [`samtools idxstats`](http://www.htslib.org/) | Per-chromosome read counts, identical format to samtools |
| stats | [`samtools stats`](http://www.htslib.org/) | Summary Numbers (SN) section, MultiQC-compatible format |

## Density scatter plots

<table>
<tr><th>dupRadar (R)</th><th>RustQC</th></tr>
<tr>
<td><img src="benchmark/dupRadar/large/duprateExpDens.png" width="400"></td>
<td><img src="benchmark/RustQC/large/dupradar/GM12878_REP1.markdup.sorted_duprateExpDens.png" width="400"></td>
</tr>
</table>

### Boxplots

<table>
<tr><th>dupRadar (R)</th><th>RustQC</th></tr>
<tr>
<td><img src="benchmark/dupRadar/large/duprateExpBoxplot.png" width="400"></td>
<td><img src="benchmark/RustQC/large/dupradar/GM12878_REP1.markdup.sorted_duprateExpBoxplot.png" width="400"></td>
</tr>
</table>

### Expression histograms

<table>
<tr><th>dupRadar (R)</th><th>RustQC</th></tr>
<tr>
<td><img src="benchmark/dupRadar/large/expressionHist.png" width="400"></td>
<td><img src="benchmark/RustQC/large/dupradar/GM12878_REP1.markdup.sorted_expressionHist.png" width="400"></td>
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

### RNA duplicate rate analysis

```bash
rustqc rna <INPUT>... (--gtf <GTF> | --bed <BED>) [OPTIONS]
```

#### Required arguments

| Argument | Description |
|----------|-------------|
| `<INPUT>...` | One or more duplicate-marked alignment files (SAM/BAM/CRAM). Duplicates must be flagged (SAM flag 0x400), not removed. BAM/CRAM files should be sorted and indexed for parallel processing. |
| `--gtf <GTF>` | Path to a GTF gene annotation file (plain or gzip-compressed). Runs all analyses (dupRadar + featureCounts + all 7 RSeQC tools). Mutually exclusive with `--bed`. |
| `--bed <BED>` | Path to a BED12 gene model file (plain or gzip-compressed). Runs RSeQC tools only (dupRadar and featureCounts are skipped). Mutually exclusive with `--gtf`. |

#### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--stranded <0\|1\|2>` | `0` | Library strandedness: `0` = unstranded, `1` = stranded (forward), `2` = reverse-stranded |
| `--paired` | `false` | Set if the library is paired-end |
| `--threads <N>` | `1` | Number of threads for parallel alignment processing |
| `--outdir <DIR>` | `.` | Output directory |
| `--reference <FASTA>` / `-r` | none | Reference FASTA file (required for CRAM input) |
| `--config <FILE>` | none | Path to a YAML configuration file (see [Configuration](#configuration)) |
| `--biotype-attribute <NAME>` | none | Override the GTF attribute used for biotype counting (default from config: `gene_biotype`). Use `gene_type` for GENCODE GTFs. |
| `--flat-output` | `false` | Write all output files directly to `--outdir` without subdirectories (legacy layout) |
| `--skip-dup-check` | `false` | Skip verification that duplicates have been marked in the BAM file (see [Duplicate marking](#duplicate-marking)) |

#### Examples

```bash
# Single BAM with GTF (all analyses: dupRadar + featureCounts + all 7 RSeQC tools)
rustqc rna sample.markdup.bam --gtf genes.gtf --outdir results/

# Single BAM with BED (RSeQC tools only, no dupRadar/featureCounts)
rustqc rna sample.markdup.bam --bed genes.bed --outdir results/

# Paired-end, reverse-stranded
rustqc rna sample.markdup.bam --gtf genes.gtf --paired --stranded 2 --outdir results/

# CRAM input with reference FASTA
rustqc rna sample.markdup.cram --gtf genes.gtf --reference genome.fa --outdir results/

# Multiple alignment files in parallel
rustqc rna sample1.bam sample2.bam sample3.bam --gtf genes.gtf --paired --threads 12 --outdir results/

# GENCODE GTF with gene_type attribute for biotypes
rustqc rna sample.markdup.bam --gtf gencode.gtf --paired --biotype-attribute gene_type --outdir results/

# Gzip-compressed annotation files (auto-detected)
rustqc rna sample.markdup.bam --gtf genes.gtf.gz --paired --outdir results/
rustqc rna sample.markdup.bam --bed genes.bed.gz --paired --outdir results/
```

### RSeQC tools

The 7 RSeQC tools are integrated into `rustqc rna` and run automatically with either `--gtf` or `--bed`. With a GTF file, all analyses run (dupRadar + featureCounts + all 7 RSeQC tools). With a BED file, only the 7 RSeQC tools run.

RSeQC-specific options:

| Option | Default | Description |
|--------|---------|-------------|
| `-q` / `--mapq <N>` | `30` | Minimum mapping quality for RSeQC tools |
| `--infer-experiment-sample-size <N>` | `200000` | Max reads to sample for strandedness inference |
| `--min-intron <N>` | `50` | Minimum intron size for junction tools |
| `--inner-distance-lower-bound <N>` | `-250` | Inner distance histogram lower bound |
| `--inner-distance-upper-bound <N>` | `250` | Inner distance histogram upper bound |

#### Duplicate marking

RustQC's `rna` subcommand requires that the input BAM file has been processed by a duplicate-marking tool such as [Picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates), [samblaster](https://github.com/GregoryFaust/samblaster), or [sambamba markdup](https://lomereiter.github.io/sambamba/). These tools set the SAM flag `0x400` on PCR/optical duplicate reads, which RustQC uses to compute duplication rates.

Before processing, RustQC checks the BAM `@PG` header lines for known duplicate-marking programs. If none are found, it exits with an error explaining how to mark duplicates. As a secondary safeguard, if processing completes but zero duplicate-flagged reads are found among mapped reads, RustQC also exits with an error.

If you are confident that your BAM file has duplicates correctly flagged despite the header check failing, you can bypass the verification with `--skip-dup-check`.

The RSeQC analyses within `rustqc rna` do not require duplicate marking -- only the dupRadar component does.

## Configuration

An optional YAML configuration file can be provided to `rustqc rna` with `--config` to control runtime behaviour. The file is designed to be extensible -- unknown fields are silently ignored, so config files remain forward-compatible.

### Chromosome name mapping

When the alignment and GTF files use different chromosome naming conventions (e.g. Ensembl `1, 2, X` vs. UCSC `chr1, chr2, chrX`), RustQC will detect the mismatch and exit with a helpful error. You can resolve this with either a prefix or explicit mapping.

**Prefix** -- prepend a string to every alignment chromosome name before matching:

```yaml
chromosome_prefix: "chr"
```

**Explicit mapping** -- for fine-grained control, map individual GTF names to alignment file names:

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

By default, outputs are organized into subdirectories by tool group. Use `--flat-output` to write all files directly to `--outdir` without subdirectories (the legacy layout).

### dupRadar outputs (`dupradar/`)

For an input file named `sample.bam`, the following files are generated in the `dupradar/` subdirectory. All outputs can be individually enabled or disabled via the [configuration file](#output-control).

| File | Description |
|------|-------------|
| `dupradar/sample_dupMatrix.txt` | Tab-separated duplication matrix (14 columns, one row per gene) |
| `dupradar/sample_duprateExpDens.{png,svg}` | Density scatter plot of duplication rate vs. expression |
| `dupradar/sample_duprateExpBoxplot.{png,svg}` | Boxplot of duplication rate by expression quantile bins |
| `dupradar/sample_expressionHist.{png,svg}` | Histogram of gene expression levels (log10 RPK) |
| `dupradar/sample_intercept_slope.txt` | Logistic regression fit parameters (intercept and slope) |
| `dupradar/sample_dup_intercept_mqc.txt` | MultiQC general stats format with intercept value |
| `dupradar/sample_duprateExpDensCurve_mqc.txt` | MultiQC line graph data for the fitted curve |

### featureCounts outputs (`featurecounts/`)

| File | Description |
|------|-------------|
| `featurecounts/sample.featureCounts.tsv` | Gene-level read counts in featureCounts format (7 columns: Geneid, Chr, Start, End, Strand, Length, counts) |
| `featurecounts/sample.featureCounts.tsv.summary` | Assignment summary statistics (Assigned, Unassigned_NoFeatures, Unassigned_Ambiguity, etc.) |
| `featurecounts/sample.biotype_counts.tsv` | Read counts aggregated by biotype (e.g., protein_coding, lncRNA, rRNA) |
| `featurecounts/sample.biotype_counts_mqc.tsv` | MultiQC bargraph format for biotype counts |
| `featurecounts/sample.biotype_counts_rrna_mqc.tsv` | MultiQC general stats with rRNA percentage |

### RSeQC outputs (`rseqc/`)

RSeQC outputs are organized into per-tool subdirectories under `rseqc/`.

| Tool | Output files | Description |
|------|-------------|-------------|
| bam_stat | `rseqc/bam_stat/sample.bam_stat.txt` | Text report with alignment statistics |
| infer_experiment | `rseqc/infer_experiment/sample.infer_experiment.txt` | Strandedness fractions |
| read_duplication | `rseqc/read_duplication/sample.pos.DupRate.xls`, `...seq.DupRate.xls`, `...DupRate_plot.{png,svg}` | Position-based and sequence-based duplication histograms with plot |
| read_distribution | `rseqc/read_distribution/sample.read_distribution.txt` | Per-region read distribution table |
| junction_annotation | `rseqc/junction_annotation/sample.junction.xls`, `...junction.bed`, `...splice_events.{png,svg}`, `...splice_junction.{png,svg}`, `...junction_plot.r`, `...junction_annotation.txt` | Junction classifications, pie charts, and summary |
| junction_saturation | `rseqc/junction_saturation/sample.junctionSaturation_plot.{png,svg}`, `...junctionSaturation_plot.r`, `...junctionSaturation_summary.txt` | Saturation curve plot and data |
| inner_distance | `rseqc/inner_distance/sample.inner_distance.txt`, `...inner_distance_freq.txt`, `...inner_distance_plot.{png,svg}`, `...inner_distance_plot.r`, `...inner_distance_summary.txt` | Per-pair distances, histogram plot, and summary |
| TIN | `rseqc/tin/sample.tin.xls`, `...summary.txt` | Transcript Integrity Number per gene and summary statistics |

### preseq outputs (`preseq/`)

| File | Description |
|------|-------------|
| `preseq/sample.lc_extrap.txt` | Library complexity extrapolation curve (TSV: TOTAL_READS, EXPECTED_DISTINCT, LOWER_CI, UPPER_CI) |

### samtools outputs (`samtools/`)

| File | Description |
|------|-------------|
| `samtools/sample.flagstat` | samtools flagstat-compatible alignment flag summary |
| `samtools/sample.idxstats` | samtools idxstats-compatible per-chromosome read counts |
| `samtools/sample.stats` | samtools stats SN-section compatible summary numbers (MultiQC-parseable) |

### Gene body coverage outputs (`qualimap/`)

| File | Description |
|------|-------------|
| `qualimap/coverage_profile_along_genes_(total).txt` | Gene body coverage profile (100 percentile bins, 5'→3') |
| `qualimap/rnaseq_qc_results.txt` | Qualimap rnaseq-compatible QC results (MultiQC-parseable) |

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

RustQC uses multi-threaded alignment processing when `--threads` is set above 1 (currently supported by `rustqc rna`). Within a single file, chromosomes are distributed across threads and processed in parallel, typically achieving near-linear speedup. When multiple alignment files are provided, they are also processed in parallel -- the available threads are divided evenly among concurrent jobs. For a single sample, `--threads 4` is a good starting point. For multiple samples, use enough threads to keep all jobs busy (e.g., `--threads 12` for 3 files gives each 4 threads). Multi-threading requires an indexed file (`.bai`/`.csi` for BAM, `.crai` for CRAM). SAM files are always processed single-threaded.

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

### `rustqc rna`

1. **GTF parsing**: Reads gene annotations (plain or gzip-compressed), computes effective gene lengths from non-overlapping exon bases, and extracts additional attributes (e.g., biotype) for downstream grouping.
2. **Read counting**: Reads the alignment file (SAM/BAM/CRAM) once, assigning each read to a gene based on exon overlap. Four count modes are tracked simultaneously:
   - With/without multimappers x with/without duplicates
   - Assignment statistics (assigned, ambiguous, no features) are tracked for the featureCounts summary.
3. **featureCounts output**: Writes gene-level counts and summary statistics in the standard featureCounts format.
4. **Biotype counting**: Aggregates assigned read counts by a configurable GTF attribute (e.g., `gene_biotype`), producing biotype count tables and MultiQC-compatible rRNA QC metrics.
5. **Duplication matrix**: Computes RPK, RPKM, and duplication rates for each gene in all four modes.
6. **Logistic regression**: Fits a binomial GLM (`dupRate ~ log10(RPK)`) using iteratively reweighted least squares (IRLS) to model the relationship between expression and duplication.
7. **Plots**: Generates density scatter, boxplot, and histogram visualizations.
8. **MultiQC integration**: Outputs files compatible with [MultiQC](https://multiqc.info/) for pipeline reporting (dupRadar intercept, fit curve, biotype bargraph, rRNA percentage).

### RSeQC tools

Each RSeQC tool is a single-pass BAM reader that processes the alignment file and produces output compatible with the original Python tools. They use a BED12 gene model for genomic feature annotation. Key implementation details:

- **bam_stat**: Flag-based read classification in a single pass (duplicates, mapping quality, splice events, proper pairs).
- **infer_experiment**: Samples reads overlapping gene intervals and determines strand protocol fractions.
- **read_duplication**: Builds position-based and sequence-based occurrence histograms with duplication rate plot.
- **read_distribution**: Classifies read tags by midpoint overlap with CDS, UTR, intron, and intergenic regions.
- **junction_annotation**: Extracts CIGAR N-operations and classifies against known splice sites from the gene model. Generates splice event and junction pie charts.
- **junction_saturation**: Subsamples junctions at increasing fractions to build a saturation curve with plot.
- **inner_distance**: Computes mRNA-level or genomic inner distance between paired-end reads with transcript-aware classification and histogram plot.

## Interpreting the results

- **Intercept** (exp(beta0)): Indicates duplication rate at low expression. Low values = good quality. High values = PCR artifact problems.
- **Slope** (exp(beta1)): Rate at which duplication increases with expression. Single-end libraries typically have higher slope than paired-end.
- **Density plot**: Good samples show low duplication (bottom of y-axis) at low expression (left), with duplication rising naturally only at very high expression.
- **1 read/bp threshold** (red dashed line): At RPK=1000, a 1kb gene has 1000 reads, meaning roughly 1 read per base pair -- near the theoretical maximum for unique reads.

## References

- Sayols S, Scherzinger D, Klein H (2016). dupRadar: a Bioconductor package for the assessment of PCR artifacts in RNA-Seq data. *BMC Bioinformatics*, 17, 428. doi:10.1186/s12859-016-1276-2
- Original dupRadar R package: https://github.com/ssayols/dupRadar
- Wang L, Wang S, Li W (2012). RSeQC: quality control of RNA-seq experiments. *Bioinformatics*, 28(16), 2184-2185. doi:10.1093/bioinformatics/bts356
- RSeQC: https://rseqc.sourceforge.net/

## License

MIT License. See [LICENSE](LICENSE) for details.
