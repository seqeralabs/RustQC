# RustQC Benchmarks

Benchmark data and scripts for comparing RustQC against
[dupRadar](https://bioconductor.org/packages/dupRadar/) (R/Bioconductor) and
[Subread featureCounts](http://subread.sourceforge.net/).

For detailed results, tables, and side-by-side plot comparisons, see the
documentation:

- [dupRadar benchmarks](https://ewels.github.io/RustQC/benchmarks/dupradar/)
- [featureCounts benchmarks](https://ewels.github.io/RustQC/benchmarks/featurecounts/)

## Benchmark data

### Small benchmark

Included in this repository. A test BAM file with a chr6-only GTF annotation
(2,905 genes).

- `small/test.bam` + `small/chr6.gtf`
- `small/dupRadar/` — R dupRadar reference output
- `small/RustQC/` — RustQC output

### Large benchmark

GM12878 REP1 — a full-size RNA-seq BAM (~10 GB) from the
[nf-core/rnaseq](https://nf-co.re/rnaseq) pipeline, duplicate-marked with
Picard. Paired-end, unstranded, aligned to GRCh38.

| File | Size | URL |
| ---- | ---- | --- |
| BAM  | ~10 GB | <https://nf-core-awsmegatests.s3-eu-west-1.amazonaws.com/rnaseq/results-3816d48abd9fab2eee41775b60b4eb8745e1fcaa/aligner_star_salmon/star_salmon/GM12878_REP1.markdup.sorted.bam> |
| GTF  | ~1.5 GB | <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz> |

## Reproducing benchmarks

### 1. Download large benchmark input files

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

### 2. Run dupRadar (R)

Requires R with `dupRadar` and `Rsubread` installed.

```bash
# Small
Rscript benchmark/small/run_dupRadar_R.R

# Large
Rscript benchmark/large/run_dupRadar_R.R
```

### 3. Run RustQC

```bash
cargo build --release

# Small
./target/release/rustqc rna benchmark/small/test.bam \
  --gtf benchmark/small/chr6.gtf -p --skip-dup-check \
  -o benchmark/small/RustQC

# Large (the GENCODE GTF uses UCSC chrom names while the BAM uses Ensembl names)
./target/release/rustqc rna benchmark/large/GM12878_REP1.markdup.sorted.bam \
  --gtf benchmark/large/genes.gtf -p -t 10 \
  -o benchmark/large/RustQC \
  -c benchmark/large/config.yaml \
  --biotype-attribute gene_type
```

> **Note:** The config file adds a `chr` prefix to alignment chromosome names
> to match the GENCODE GTF. The `--biotype-attribute gene_type` flag is needed
> because GENCODE uses `gene_type` instead of the Ensembl default `gene_biotype`.

### 4. Compare results

```bash
# Compare duplication matrices cell-by-cell
python3 -c "
import csv
with open('benchmark/large/dupRadar/dupMatrix.txt') as rf, \
     open('benchmark/large/RustQC/GM12878_REP1.markdup.sorted_dupMatrix.txt') as rustf:
    r = list(csv.reader(rf, delimiter='\t'))
    rust = list(csv.reader(rustf, delimiter='\t'))
    mismatches = sum(1 for i in range(1, len(r)) for j in range(1, len(r[i]))
                     if r[i][j] != 'NA' and rust[i][j] != 'NA'
                     and abs(float(r[i][j]) - float(rust[i][j])) / max(abs(float(r[i][j])), 1e-15) > 1e-6)
    total = sum(1 for i in range(1, len(r)) for j in range(1, len(r[i])))
    print(f'{total} values compared, {mismatches} mismatches')
"
```
