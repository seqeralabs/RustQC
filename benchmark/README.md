# RustQC Benchmarks 🧬 🦀

Comparison of [dupRadar](https://bioconductor.org/packages/dupRadar/) (R/Bioconductor) and RustQC on the same input data.

## Small benchmark

A small test BAM file (`test.bam`) with a chr6-only GTF annotation, included in this repository.

### Results

| Metric | dupRadar (R) | RustQC |
| --- | --- | --- |
| **Intercept** | 0.03 | 0.03 |
| **Slope** | 1.60 | 1.60 |
| **Genes total** | 2,905 | 2,905 |
| **Genes with reads** | 636 | 636 |

#### Count comparison

| Metric | dupRadar (R) | RustQC | Exact match |
| --- | ---: | ---: | ---: |
| **allCounts (unique)** | 20,449 | 20,449 | 100% |
| **filteredCounts (unique)** | 17,879 | 17,879 | 100% |
| **allCountsMulti** | 22,812 | 22,812 | **100%** |
| **filteredCountsMulti** | 20,034 | 20,034 | **100%** |

### Replication

```bash
# dupRadar (R)
Rscript benchmark/small/run_dupRadar_R.R

# RustQC
cargo run --release -- rna benchmark/small/test.bam --gtf benchmark/small/chr6.gtf -p -o benchmark/small/RustQC
```

## Large benchmark

**GM12878 REP1** — a full-size RNA-seq BAM from the [nf-core/rnaseq](https://nf-co.re/rnaseq) pipeline, duplicate-marked with picard.
Paired-end, unstranded, aligned to GRCh38 (Ensembl chromosome names).

### Input data

| File | Size | URL |
| ---- | ---- | --- |
| BAM  | ~10 GB | <https://nf-core-awsmegatests.s3-eu-west-1.amazonaws.com/rnaseq/results-3816d48abd9fab2eee41775b60b4eb8745e1fcaa/aligner_star_salmon/star_salmon/GM12878_REP1.markdup.sorted.bam> |
| GTF  | ~1.5 GB | <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz> |

### Results

| Metric | dupRadar (R) | RustQC (1 thread) | RustQC (8 threads) | RustQC (10 threads) |
| --- | --- | --- | --- | --- |
| **Runtime** | 23m 48s | 3m 20s (~7x) | 1m 04s (~22x) | 0m 53s (~27x) |
| **Speedup** | — | **~7x** | **~22x** | **~27x** |
| **Intercept** | 0.8245 | 0.8245 | 0.8245 | 0.8245 |
| **Slope** | 1.6774 | 1.6774 | 1.6774 | 1.6774 |
| **Genes total** | 63,086 | 63,086 | 63,086 | 63,086 |
| **Genes with reads (unique)** | 23,597 | 23,597 | 23,597 | 23,597 |
| **Genes with reads (multi)** | 24,719 | 24,719 | 24,719 | 24,719 |

### Count comparison

| Metric | dupRadar (R) | RustQC | Exact match |
| --- | ---: | ---: | ---: |
| **allCounts (unique)** | 14,654,579 | 14,654,579 | **100%** |
| **filteredCounts (unique)** | 3,599,832 | 3,599,832 | **100%** |
| **allCountsMulti** | 16,089,488 | 16,089,488 | **100%** |
| **filteredCountsMulti** | 4,503,920 | 4,503,920 | **100%** |

All four count columns match exactly across all 63,086 genes — both unique-mapper and multi-mapper counts.

Model fit parameters (**intercept** and **slope**) match to the displayed precision.

### Replication

#### 1. Download input files

```bash
mkdir -p benchmark/large

# Download BAM (~10 GB)
curl -L -o benchmark/large/GM12878_REP1.markdup.sorted.bam \
  "https://nf-core-awsmegatests.s3-eu-west-1.amazonaws.com/rnaseq/results-3816d48abd9fab2eee41775b60b4eb8745e1fcaa/aligner_star_salmon/star_salmon/GM12878_REP1.markdup.sorted.bam"

# Download GTF
curl -L -o benchmark/large/genes.gtf.gz \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz"
gunzip benchmark/large/genes.gtf.gz

# Index BAM
samtools index benchmark/large/GM12878_REP1.markdup.sorted.bam
```

#### 2. Run dupRadar (R)

Requires R with `dupRadar` and `Rsubread` installed.

```bash
Rscript benchmark/large/run_dupRadar_R.R
```

#### 3. Run RustQC

The alignment file uses Ensembl chromosome names (`1`, `2`, ...) but the GENCODE GTF uses UCSC names (`chr1`, `chr2`, ...).
A config file is used to add the `chr` prefix to alignment chromosome names:

```yaml
# benchmark/large/config.yaml
chromosome_prefix: "chr"
```

```bash
cargo build --release

# Single-threaded
./target/release/rustqc rna \
  benchmark/large/GM12878_REP1.markdup.sorted.bam \
  --gtf benchmark/large/genes.gtf \
  -p \
  -o benchmark/large/RustQC \
  -c benchmark/large/config.yaml

# Multi-threaded (8 threads)
./target/release/rustqc rna \
  benchmark/large/GM12878_REP1.markdup.sorted.bam \
  --gtf benchmark/large/genes.gtf \
  -p \
  -t 8 \
  -o benchmark/large/RustQC \
  -c benchmark/large/config.yaml
```

### Known differences

Both benchmarks achieve **100% exact match** across all four count columns (unique and multi-mapper) and all genes. Model fit parameters match to at least 10 significant digits.

Runtime may vary depending on hardware — the times above were measured on a single machine (10-core Apple Silicon) for relative comparison. Multi-threaded scaling depends on the number of chromosomes with mapped reads and the evenness of their read distribution.
