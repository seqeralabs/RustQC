# 🧬 dupRust Benchmarks 🦀

Comparison of [dupRadar](https://bioconductor.org/packages/dupRadar/) (R/Bioconductor) and dupRust on the same input data.

## Small benchmark

A small test BAM file (`test.bam`) with a chr6-only GTF annotation, included in this repository.

### Results

| Metric | dupRadar (R) | dupRust |
| --- | --- | --- |
| **Intercept** | 0.03 | 0.03 |
| **Slope** | 1.6 | 1.5 |
| **Genes total** | 2,905 | 2,905 |
| **Genes with reads** | 621 | 621 |

### Replication

```bash
# dupRadar (R)
Rscript benchmark/run_dupRadar_R.R

# dupRust
cargo run --release -- benchmark/test.bam benchmark/chr6.gtf -o benchmark/rust_output
```

---

## Large benchmark

**GM12878 REP1** — a full-size RNA-seq BAM from the [nf-core/rnaseq](https://nf-co.re/rnaseq) pipeline, duplicate-marked with picard.
Paired-end, unstranded, aligned to GRCh38 (Ensembl chromosome names).

### Input data

| File | Size | URL |
| ---- | ---- | --- |
| BAM  | ~10 GB | <https://nf-core-awsmegatests.s3-eu-west-1.amazonaws.com/rnaseq/results-3816d48abd9fab2eee41775b60b4eb8745e1fcaa/aligner_star_salmon/star_salmon/GM12878_REP1.markdup.sorted.bam> |
| GTF  | ~1.5 GB | <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz> |

### Results

| Metric | dupRadar (R) | dupRust |
| --- | --- | --- |
| **Runtime** | 1,428 s (23.8 min) | 178 s (3.0 min) |
| **Speedup** | — | **~8x** |
| **Intercept** | 0.8245 | 0.8213 |
| **Slope** | 1.6774 | 1.6859 |
| **Genes total** | 63,086 | 63,086 |
| **Genes with reads** | 24,719 | 23,159 |

### Count comparison

| Metric | dupRadar (R) | dupRust | Ratio |
| --- | ---: | ---: | ---: |
| **Total allCountsMulti** | 16,097,284 | 14,237,148 | 0.88 |
| **Total filteredCountsMulti** | 4,671,832 | 4,020,122 | 0.86 |
| **Total allCounts (unique)** | 14,710,916 | 13,694,224 | 0.93 |
| **Total filteredCounts (unique)** | 4,305,768 | 3,830,538 | 0.89 |

#### Count correlation (excluding mitochondrial genes)

| Metric | Pearson r |
| --- | --- |
| allCountsMulti | 0.973 |
| allCounts (unique) | 0.985 |

#### Per-gene exact match rates

| Metric | Exact matches | Percentage |
| --- | ---: | ---: |
| allCountsMulti | 51,670 / 63,086 | 81.9% |
| allCounts (unique) | 56,436 / 63,086 | 89.5% |

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

#### 3. Run dupRust

The BAM uses Ensembl chromosome names (`1`, `2`, ...) but the GENCODE GTF uses UCSC names (`chr1`, `chr2`, ...).
A config file is used to add the `chr` prefix to BAM chromosome names:

```yaml
# benchmark/large/config.yaml
chromosome_prefix: "chr"
```

```bash
cargo build --release
./target/release/duprust \
  benchmark/large/GM12878_REP1.markdup.sorted.bam \
  benchmark/large/genes.gtf \
  -p \
  -o benchmark/large/rust_output \
  -c benchmark/large/config.yaml
```

### Known differences

- **1,221 genes** have reads in R but not Rust — these are on alternative contigs / patches not present in the BAM reference.
- Count differences stem from Rsubread's `featureCounts` (used by dupRadar) vs dupRust's CIGAR-aware exon overlap counting. featureCounts runs 4 separate counting passes with different duplicate/multimap filters; dupRust does a single pass tracking all flags simultaneously.
- Model fit parameters (intercept, slope) agree within ~0.5%.
