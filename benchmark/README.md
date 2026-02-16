# RustQC Benchmarks

Benchmark data and scripts for comparing RustQC against
[dupRadar](https://bioconductor.org/packages/dupRadar/) (R/Bioconductor),
[Subread featureCounts](http://subread.sourceforge.net/),
[RSeQC](https://rseqc.sourceforge.net/),
[preseq](https://github.com/smithlabcode/preseq),
[samtools](http://www.htslib.org/), and
[Qualimap](http://qualimap.conesalab.org/).

For detailed results, tables, and side-by-side plot comparisons, see the
documentation:

- [Combined benchmarks](https://ewels.github.io/RustQC/benchmarks/combined/)
- [dupRadar benchmarks](https://ewels.github.io/RustQC/benchmarks/dupradar/)
- [featureCounts benchmarks](https://ewels.github.io/RustQC/benchmarks/featurecounts/)
- [Preseq benchmarks](https://ewels.github.io/RustQC/benchmarks/preseq/)
- [Samtools benchmarks](https://ewels.github.io/RustQC/benchmarks/samtools/)

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
| preseq lc_extrap | ~4m |
| samtools flagstat | ~1m |
| samtools idxstats | ~1s |
| samtools stats | ~3m |
| tin.py (RSeQC) | ~30m |
| **Traditional total** | **~2h 45m** |
| **RustQC (10 threads)** | **~5m** |

> **Note:** Docker x86 emulation on ARM inflates the upstream tool timings.
> The key takeaway is that RustQC replaces 9 separate tool invocations
> (each requiring a full BAM pass) with a single command and single BAM pass.

## Directory structure

```
benchmark/
  input/            — Shared input files (BAM, GTF, BED, config)
    large/          — Full-size GM12878 dataset (~10 GB BAM)
    small/          — Small chr6 test dataset
  RustQC/           — RustQC output (all tools, single-pass)
    large/          — Output for large dataset
    small/          — Output for small dataset
  dupRadar/         — R dupRadar reference output
    large/          — R script + output for large dataset
    small/          — R script + output for small dataset
  rseqc/            — Python RSeQC reference output
    large/          — Reference output for large dataset
      bam_stat/
      infer_experiment/
      read_distribution/
      read_duplication/
      junction_annotation/
      junction_saturation/
      inner_distance/
      tin/
    small/          — Reference output for small dataset
      (same tool subdirectories)
  samtools/         — samtools reference output
    large/          — flagstat, idxstats, stats
    small/
  preseq/           — preseq lc_extrap reference output
    large/          — PE and SE reference extrapolations
    small/
  qualimap/         — Qualimap rnaseq reference output
    large/          — Gene body coverage, QC results
    small/
```

## Benchmark data

### Small benchmark

Included in this repository. A test BAM file with chr6 reads, a chr6-only GTF
annotation (2,905 genes), and a BED12 gene model for RSeQC tools.

- `input/small/test.bam` + `input/small/chr6.gtf` + `input/small/chr6.bed`
- `dupRadar/small/` — R dupRadar reference output
- `rseqc/small/` — Python RSeQC reference output (per-tool subdirectories)
- `samtools/small/` — samtools reference output
- `preseq/small/` — preseq reference output
- `qualimap/small/` — Qualimap reference output
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
`rseqc/small/` and `rseqc/large/` in per-tool subdirectories. These were generated with
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
  $RSEQC_IMG bam_stat.py -i /data/test.bam > benchmark/rseqc/small/bam_stat/bam_stat.txt 2>&1

# infer_experiment
docker run --rm --platform linux/amd64 -v $(pwd)/benchmark/input/small:/data \
  $RSEQC_IMG infer_experiment.py -i /data/test.bam -r /data/chr6.bed \
  > benchmark/rseqc/small/infer_experiment/infer_experiment.txt 2>&1

# read_duplication
docker run --rm --platform linux/amd64 -v $(pwd)/benchmark/input/small:/data -v $(pwd)/benchmark/rseqc/small/read_duplication:/out \
  $RSEQC_IMG read_duplication.py -i /data/test.bam -o /out/read_duplication

# read_distribution
docker run --rm --platform linux/amd64 -v $(pwd)/benchmark/input/small:/data \
  $RSEQC_IMG read_distribution.py -i /data/test.bam -r /data/chr6.bed \
  > benchmark/rseqc/small/read_distribution/read_distribution.txt 2>&1

# junction_annotation
docker run --rm --platform linux/amd64 -v $(pwd)/benchmark/input/small:/data -v $(pwd)/benchmark/rseqc/small/junction_annotation:/out \
  $RSEQC_IMG junction_annotation.py -i /data/test.bam -r /data/chr6.bed -o /out/junction_annotation

# junction_saturation
docker run --rm --platform linux/amd64 -v $(pwd)/benchmark/input/small:/data -v $(pwd)/benchmark/rseqc/small/junction_saturation:/out \
  $RSEQC_IMG junction_saturation.py -i /data/test.bam -r /data/chr6.bed -o /out/junction_saturation

# inner_distance
docker run --rm --platform linux/amd64 -v $(pwd)/benchmark/input/small:/data -v $(pwd)/benchmark/rseqc/small/inner_distance:/out \
  $RSEQC_IMG inner_distance.py -i /data/test.bam -r /data/chr6.bed -o /out/inner_distance
```

### 7. Generate samtools reference outputs

Reference outputs for validating the samtools reimplementations are stored in
`samtools/small/` and `samtools/large/`.

```bash
SAMTOOLS_IMG="quay.io/biocontainers/samtools:1.21--h50ea8bc_0"

# Small
docker run --rm -v $(pwd)/benchmark/input/small:/data $SAMTOOLS_IMG \
  samtools flagstat /data/test.bam > benchmark/samtools/small/flagstat.txt
docker run --rm -v $(pwd)/benchmark/input/small:/data $SAMTOOLS_IMG \
  samtools idxstats /data/test.bam > benchmark/samtools/small/idxstats.txt
docker run --rm -v $(pwd)/benchmark/input/small:/data $SAMTOOLS_IMG \
  samtools stats /data/test.bam > benchmark/samtools/small/stats.txt

# Large
docker run --rm -v $(pwd)/benchmark/input/large:/data $SAMTOOLS_IMG \
  samtools flagstat /data/GM12878_REP1.markdup.sorted.bam > benchmark/samtools/large/flagstat.txt
docker run --rm -v $(pwd)/benchmark/input/large:/data $SAMTOOLS_IMG \
  samtools idxstats /data/GM12878_REP1.markdup.sorted.bam > benchmark/samtools/large/idxstats.txt
docker run --rm -v $(pwd)/benchmark/input/large:/data $SAMTOOLS_IMG \
  samtools stats /data/GM12878_REP1.markdup.sorted.bam > benchmark/samtools/large/stats.txt
```

### 8. Generate preseq reference outputs

Reference outputs for validating the preseq `lc_extrap` reimplementation are stored
in `preseq/small/` and `preseq/large/`. These were generated with
[preseq 3.2.0](https://github.com/smithlabcode/preseq) run via Docker.

Two reference files are provided per dataset:

| File | Mode | Command |
| ---- | ---- | ------- |
| `lc_extrap_se.txt` | SE (default) | `preseq lc_extrap -bam -seed 1` |
| `lc_extrap_pe.txt` | PE (nf-core style) | `preseq lc_extrap -bam -pe -seed 1 -seg_len 100000000` |

The **PE reference** matches the invocation used by
[nf-core/rnaseq](https://nf-co.re/rnaseq), which auto-adds `-pe` for
paired-end data with a large `-seg_len` window. Both references use the
coordinate-sorted BAM directly (no name-sorting required).

```bash
PRESEQ_IMG="quay.io/biocontainers/preseq:3.2.0--hdcf5f25_6"

# Small
docker run --rm --platform linux/amd64 \
  -v $(pwd)/benchmark:/data \
  $PRESEQ_IMG preseq lc_extrap -bam -pe -seed 1 -seg_len 100000000 \
  -output /data/preseq/small/lc_extrap_pe.txt /data/input/small/test.bam

# Large — SE reference
docker run --rm --platform linux/amd64 \
  -v $(pwd)/benchmark:/data \
  $PRESEQ_IMG preseq lc_extrap -bam -seed 1 \
  -output /data/preseq/large/lc_extrap_se.txt \
  /data/input/large/GM12878_REP1.markdup.sorted.bam

# Large — PE reference (nf-core/rnaseq style)
docker run --rm --platform linux/amd64 \
  -v $(pwd)/benchmark:/data \
  $PRESEQ_IMG preseq lc_extrap -bam -pe -seed 1 -seg_len 100000000 \
  -output /data/preseq/large/lc_extrap_pe.txt \
  /data/input/large/GM12878_REP1.markdup.sorted.bam
```

> **Note:** preseq's `-pe` mode on coordinate-sorted BAMs uses the `-seg_len`
> parameter to find mates within a window. Without `-pe`, preseq treats every
> mapped read independently using only `{tid, pos}` as the unique identifier,
> even for paired-end BAMs. The nf-core/rnaseq pipeline always passes `-pe`
> for paired-end data.

### 9. Generate Qualimap reference outputs

Reference outputs for validating the gene body coverage reimplementation are
stored in `qualimap/small/` and `qualimap/large/`. These were generated with
[Qualimap 2.3](http://qualimap.conesalab.org/) run via Docker.

```bash
QUALIMAP_IMG="quay.io/biocontainers/qualimap:2.3--hdfd78af_0"

# Small (requires uncompressed GTF)
gunzip -k benchmark/input/small/chr6.gtf.gz
docker run --rm -v $(pwd)/benchmark:/data \
  $QUALIMAP_IMG qualimap rnaseq \
  -bam /data/input/small/test.bam -gtf /data/input/small/chr6.gtf \
  -outdir /data/qualimap_tmp_small --java-mem-size=4G -pe
cp benchmark/qualimap_tmp_small/rnaseq_qc_results.txt benchmark/qualimap/small/
cp "benchmark/qualimap_tmp_small/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt" benchmark/qualimap/small/

# Large
gunzip -k benchmark/input/large/genes.gtf.gz
docker run --rm -v $(pwd)/benchmark:/data \
  $QUALIMAP_IMG qualimap rnaseq \
  -bam /data/input/large/GM12878_REP1.markdup.sorted.bam -gtf /data/input/large/genes.gtf \
  -outdir /data/qualimap_tmp_large --java-mem-size=8G -pe
cp benchmark/qualimap_tmp_large/rnaseq_qc_results.txt benchmark/qualimap/large/
cp "benchmark/qualimap_tmp_large/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt" benchmark/qualimap/large/
```

### 10. Compare results

```bash
# Compare duplication matrices cell-by-cell
python3 -c "
import csv
with open('benchmark/dupRadar/large/dupMatrix.txt') as rf, \
     open('benchmark/RustQC/large/dupradar/GM12878_REP1.markdup.sorted_dupMatrix.txt') as rustf:
    r = list(csv.reader(rf, delimiter='\t'))
    rust = list(csv.reader(rustf, delimiter='\t'))
    mismatches = sum(1 for i in range(1, len(r)) for j in range(1, len(r[i]))
                     if r[i][j] != 'NA' and rust[i][j] != 'NA'
                     and abs(float(r[i][j]) - float(rust[i][j])) / max(abs(float(r[i][j])), 1e-15) > 1e-6)
    total = sum(1 for i in range(1, len(r)) for j in range(1, len(r[i])))
    print(f'{total} values compared, {mismatches} mismatches')
"
```
