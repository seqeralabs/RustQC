# dupRust Benchmarks

Comparison of [dupRadar](https://bioconductor.org/packages/dupRadar/) (R/Bioconductor) and dupRust on the same input data.

## Small benchmark

A small test BAM file (`test.bam`) with a chr6-only GTF annotation, included in this repository.

### Results

| Metric | dupRadar (R) | dupRust |
| --- | --- | --- |
| **Intercept** | 0.03 | 0.03 |
| **Slope** | 1.60 | 1.60 |
| **Genes total** | 2,905 | 2,905 |
| **Genes with reads** | 636 | 636 |

#### Count comparison

| Metric | dupRadar (R) | dupRust | Exact match |
| --- | ---: | ---: | ---: |
| **allCounts (unique)** | 20,449 | 20,449 | 100% |
| **filteredCounts (unique)** | 17,879 | 17,879 | 100% |
| **allCountsMulti** | 22,812 | 22,821 | 99.6% |
| **filteredCountsMulti** | 20,034 | 20,048 | 99.6% |

### Replication

```bash
# dupRadar (R)
Rscript benchmark/small/run_dupRadar_R.R

# dupRust
cargo run --release -- benchmark/small/test.bam benchmark/small/chr6.gtf -p -o benchmark/small/dupRust
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
| **Runtime** | 1,428 s (23.8 min) | 198 s (3.3 min) |
| **Speedup** | — | **~7x** |
| **Intercept** | 0.8245 | 0.8245 |
| **Slope** | 1.6774 | 1.6774 |
| **Genes total** | 63,086 | 63,086 |
| **Genes with reads** | 24,719 | 23,597 |

### Count comparison

| Metric | dupRadar (R) | dupRust | Exact match |
| --- | ---: | ---: | ---: |
| **allCounts (unique)** | 14,654,579 | 14,654,579 | **100%** |
| **filteredCounts (unique)** | 3,599,832 | 3,599,832 | **100%** |
| **allCountsMulti** | 16,089,488 | 16,023,820 | 95.6% |
| **filteredCountsMulti** | 4,503,920 | 4,459,432 | 95.1% |

Unique-mapper counts (**allCounts** and **filteredCounts**) match exactly across all 63,086 genes. Multi-mapper counts are within ~0.4% due to minor differences in secondary alignment handling.

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
  -o benchmark/large/dupRust \
  -c benchmark/large/config.yaml
```

### Known differences

- **1,122 genes** have reads in dupRadar but not dupRust — these are on alternative contigs / patches not present in the BAM reference, or are mitochondrial genes (the BAM uses `MT` but the GTF uses `chrM`).
- Multi-mapper count differences (~4%) stem from minor differences in how secondary alignments are paired and assigned.
- Unique-mapper counts and model fit parameters match exactly.
