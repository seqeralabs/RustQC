# RustQC Benchmarks

Benchmark data and scripts for comparing RustQC against
[dupRadar](https://bioconductor.org/packages/dupRadar/) (R/Bioconductor),
[Subread featureCounts](http://subread.sourceforge.net/), and
[RSeQC](https://rseqc.sourceforge.net/).

For detailed results, tables, and side-by-side plot comparisons, see the
documentation:

- [Combined benchmarks](https://ewels.github.io/RustQC/benchmarks/combined/)
- [dupRadar benchmarks](https://ewels.github.io/RustQC/benchmarks/dupradar/)
- [featureCounts benchmarks](https://ewels.github.io/RustQC/benchmarks/featurecounts/)

## Latest results (large dataset)

Run via `python3 benchmark/run_benchmarks.py --all` on a 10-core Apple Silicon
Mac. Upstream tools run via Docker with x86 emulation; RustQC runs natively.

| Tool | Wall time |
|------|----------:|
| dupRadar (R) | 27m 21s |
| featureCounts (Subread) | 3m 39s |
| bam_stat (RSeQC) | 6m 07s |
| infer_experiment (RSeQC) | 7s |
| read_duplication (RSeQC) | 29m 43s |
| read_distribution (RSeQC) | 6m 00s |
| junction_annotation (RSeQC) | 4m 37s |
| junction_saturation (RSeQC) | 6m 32s |
| inner_distance (RSeQC) | 1m 09s |
| **Traditional total** | **1h 25m** |
| **RustQC (10 threads)** | **16m 11s** |

> **Note:** Docker x86 emulation on ARM inflates the upstream tool timings.
> The key takeaway is that RustQC replaces 9 separate tool invocations
> (each requiring a full BAM pass) with a single command and single BAM pass.

## Directory structure

```
benchmark/
  input/            — Shared input files (BAM, GTF, BED, config)
    large/          — Full-size GM12878 dataset (~10 GB BAM)
    small/          — Small chr6 test dataset
  dupRadar/         — R dupRadar reference output
    large/          — R script + output for large dataset
    small/          — R script + output for small dataset
  RustQC/           — RustQC output (dupRadar + featureCounts)
    large/          — Output for large dataset
    small/          — Output for small dataset
  RSeQC/            — Python RSeQC reference output
    large/          — Reference output for large dataset
    small/          — Reference output for small dataset
```

## Benchmark data

### Small benchmark

Included in this repository. A test BAM file with chr6 reads, a chr6-only GTF
annotation (2,905 genes), and a BED12 gene model for RSeQC tools.

- `input/small/test.bam` + `input/small/chr6.gtf` + `input/small/chr6.bed`
- `dupRadar/small/` — R dupRadar reference output
- `RustQC/small/` — RustQC output

### Large benchmark

GM12878 REP1 — a full-size RNA-seq BAM (~10 GB) from the
[nf-core/rnaseq](https://nf-co.re/rnaseq) pipeline, duplicate-marked with
Picard. Paired-end, unstranded, aligned to GRCh37 (Ensembl chromosome names).

| File | Size | URL |
| ---- | ---- | --- |
| BAM  | ~10 GB | <https://nf-core-awsmegatests.s3-eu-west-1.amazonaws.com/rnaseq/results-3816d48abd9fab2eee41775b60b4eb8745e1fcaa/aligner_star_salmon/star_salmon/GM12878_REP1.markdup.sorted.bam> |
| GTF  | ~1.5 GB | <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz> |

The dupRadar/featureCounts pipeline uses a GENCODE v46 GTF (`input/large/genes.gtf`). The
RSeQC tools use a BED12 gene model (`input/large/genes.bed`) converted from a
genome-matched GTF with matching Ensembl chromosome names.

### RSeQC reference outputs

Reference outputs for validating RSeQC reimplementations are stored in
`RSeQC/small/` and `RSeQC/large/`. These were generated with
[RSeQC 5.0.4](https://rseqc.sourceforge.net/) run via Docker
(`--platform linux/amd64`).

**Small BAM** — All 7 tools have reference output:

| Tool | Key metrics |
| ---- | ----------- |
| `bam_stat` | 43,476 total reads, 6,032 duplicates |
| `infer_experiment` | Undetermined (50/50 strand split) |
| `read_duplication` | Position-based and sequence-based histograms |
| `read_distribution` | 68,660 tags, 66,693 assigned (49,589 CDS) |
| `junction_annotation` | 3,261 junctions (2,982 known, 88 partial novel, 191 novel) |
| `junction_saturation` | Saturation curve, 20 sampling points |
| `inner_distance` | 20,861 pairs (mean -38.85, median -72.5) |

**Large BAM** — 6 of 7 tools have reference output (`read_duplication`
pending):

| Tool | Key metrics |
| ---- | ----------- |
| `bam_stat` | 185,718,543 total records, 133,912,519 duplicates, 39,827,099 unique (MAPQ>=30) |
| `infer_experiment` | Reverse stranded (92.2% / 1.2%) |
| `read_distribution` | 55,374,023 tags, 52,400,513 assigned (33.3M CDS, 8.5M 3'UTR) |
| `junction_annotation` | 256,466 junctions (178,797 known, 50,936 partial novel, 26,733 novel) |
| `junction_saturation` | 256,466 junctions at 100% (163,710 known, 92,756 novel) |
| `inner_distance` | 1,000,000 pairs sampled (mean 29.43, median 27.5, SD 32.80) |

## Running benchmarks

The simplest way to run benchmarks is with the included script:

```bash
# Run everything (requires Docker for upstream tools)
python3 benchmark/run_benchmarks.py --all

# Run only specific tools
python3 benchmark/run_benchmarks.py --rustqc --dupradar

# Run only the small dataset
python3 benchmark/run_benchmarks.py --all --small-only
```

The script profiles CPU/memory usage with `psrecord`, saves timing results to
`benchmark/profiling/results.json`, per-tool profile logs and plots under
`benchmark/profiling/<tool>/<dataset>/`, and regenerates the SVG bar charts
in `docs/public/benchmarks/`.

Dependencies: `pip install psrecord matplotlib psutil`

## Manual reproduction

### 1. Download large benchmark input files

```bash
mkdir -p benchmark/input/large

# Download BAM (~10 GB)
curl -L -o benchmark/input/large/GM12878_REP1.markdup.sorted.bam \
  "https://nf-core-awsmegatests.s3-eu-west-1.amazonaws.com/rnaseq/results-3816d48abd9fab2eee41775b60b4eb8745e1fcaa/aligner_star_salmon/star_salmon/GM12878_REP1.markdup.sorted.bam"

# Download GTF
curl -L -o benchmark/input/large/genes.gtf.gz \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz"
gunzip benchmark/input/large/genes.gtf.gz

# Index BAM
samtools index benchmark/input/large/GM12878_REP1.markdup.sorted.bam
```

### 2. Generate BED12 for RSeQC tools

Convert the GTF to BED12 format using UCSC tools. This can be done via Docker
with the [wave CLI](https://github.com/seqeralabs/wave-cli):

```bash
wave --conda ucsc-gtftogenepred ucsc-genepredtobed

# Then run the conversion (substitute the wave image name):
docker run --rm -v $(pwd)/benchmark/input/large:/data <wave-image> \
  bash -c 'gtfToGenePred /data/genes_matched.gtf /tmp/genes.genePred && \
           genePredToBed /tmp/genes.genePred /data/genes.bed'
```

### 3. Run dupRadar (R)

Requires R with `dupRadar` and `Rsubread` installed.

```bash
# Small
Rscript benchmark/dupRadar/small/run_dupRadar_R.R

# Large
Rscript benchmark/dupRadar/large/run_dupRadar_R.R
```

### 4. Run RustQC

```bash
cargo build --release

# Small
./target/release/rustqc rna benchmark/input/small/test.bam \
  --gtf benchmark/input/small/chr6.gtf -p --skip-dup-check \
  -o benchmark/RustQC/small

# Large (the GENCODE GTF uses UCSC chrom names while the BAM uses Ensembl names)
./target/release/rustqc rna benchmark/input/large/GM12878_REP1.markdup.sorted.bam \
  --gtf benchmark/input/large/genes.gtf -p -t 10 \
  -o benchmark/RustQC/large \
  -c benchmark/input/large/config.yaml \
  --biotype-attribute gene_type
```

> **Note:** The config file adds a `chr` prefix to alignment chromosome names
> to match the GENCODE GTF. The `--biotype-attribute gene_type` flag is needed
> because GENCODE uses `gene_type` instead of the Ensembl default `gene_biotype`.

### 5. RSeQC tools (RustQC reimplementations)

All 7 RSeQC tools are integrated into `rustqc rna` and run automatically when
`--gtf` is provided (the commands in step 4 already include all RSeQC analyses).
No separate `--bed` file is needed — all required data is derived from the GTF.

### 6. Generate RSeQC Python reference outputs

Requires RSeQC 5.0.4 via Docker:

```bash
RSEQC_IMG="wave.seqera.io/wt/ea3e9f972b6e/wave/build:rseqc-5.0.4--14c99cde3bff8d57"

# bam_stat
docker run --rm --platform linux/amd64 -v $(pwd)/benchmark/input/small:/data \
  $RSEQC_IMG bam_stat.py -i /data/test.bam > benchmark/RSeQC/small/bam_stat.txt 2>&1

# infer_experiment
docker run --rm --platform linux/amd64 -v $(pwd)/benchmark/input/small:/data \
  $RSEQC_IMG infer_experiment.py -i /data/test.bam -r /data/chr6.bed \
  > benchmark/RSeQC/small/infer_experiment.txt 2>&1

# read_duplication
docker run --rm --platform linux/amd64 -v $(pwd)/benchmark/input/small:/data -v $(pwd)/benchmark/RSeQC/small:/out \
  $RSEQC_IMG read_duplication.py -i /data/test.bam -o /out/read_duplication

# read_distribution
docker run --rm --platform linux/amd64 -v $(pwd)/benchmark/input/small:/data \
  $RSEQC_IMG read_distribution.py -i /data/test.bam -r /data/chr6.bed \
  > benchmark/RSeQC/small/read_distribution.txt 2>&1

# junction_annotation
docker run --rm --platform linux/amd64 -v $(pwd)/benchmark/input/small:/data -v $(pwd)/benchmark/RSeQC/small:/out \
  $RSEQC_IMG junction_annotation.py -i /data/test.bam -r /data/chr6.bed -o /out/junction_annotation

# junction_saturation
docker run --rm --platform linux/amd64 -v $(pwd)/benchmark/input/small:/data -v $(pwd)/benchmark/RSeQC/small:/out \
  $RSEQC_IMG junction_saturation.py -i /data/test.bam -r /data/chr6.bed -o /out/junction_saturation

# inner_distance
docker run --rm --platform linux/amd64 -v $(pwd)/benchmark/input/small:/data -v $(pwd)/benchmark/RSeQC/small:/out \
  $RSEQC_IMG inner_distance.py -i /data/test.bam -r /data/chr6.bed -o /out/inner_distance
```

### 7. Compare results

```bash
# Compare duplication matrices cell-by-cell
python3 -c "
import csv
with open('benchmark/dupRadar/large/dupMatrix.txt') as rf, \
     open('benchmark/RustQC/large/GM12878_REP1.markdup.sorted_dupMatrix.txt') as rustf:
    r = list(csv.reader(rf, delimiter='\t'))
    rust = list(csv.reader(rustf, delimiter='\t'))
    mismatches = sum(1 for i in range(1, len(r)) for j in range(1, len(r[i]))
                     if r[i][j] != 'NA' and rust[i][j] != 'NA'
                     and abs(float(r[i][j]) - float(rust[i][j])) / max(abs(float(r[i][j])), 1e-15) > 1e-6)
    total = sum(1 for i in range(1, len(r)) for j in range(1, len(r[i])))
    print(f'{total} values compared, {mismatches} mismatches')
"
```
