# RustQC Changelog 🧬 🦀

## [Unreleased]

### RSeQC reimplementations

Seven new subcommands reimplementing the most widely-used
[RSeQC](https://rseqc.sourceforge.net/) Python tools as fast, single-pass Rust
equivalents. Output formats are compatible with the Python originals.

- **`rustqc bam-stat`** -- Basic alignment statistics (total records, duplicates,
  mapping quality distribution, splice reads, proper pairs)
- **`rustqc infer-experiment`** -- Infer library strandedness from read/gene-model
  overlap (samples reads against a BED12 gene model)
- **`rustqc read-duplication`** -- Position-based and sequence-based read duplication
  histograms
- **`rustqc read-distribution`** -- Read distribution across genomic features (CDS,
  5'/3' UTR, introns, intergenic, TSS/TES flanking regions)
- **`rustqc junction-annotation`** -- Classify splice junctions as known, partial
  novel, or complete novel against a reference gene model
- **`rustqc junction-saturation`** -- Saturation analysis of detected splice junctions
  at increasing subsampling percentages
- **`rustqc inner-distance`** -- Inner distance distribution for paired-end reads with
  transcript-aware classification

All RSeQC tools support multiple input files, MAPQ filtering (`-q`), and CRAM
input with `--reference`.

### Features

- **featureCounts-compatible output**: Generate counts TSV and summary files in the same
  format as Subread featureCounts, produced in the same single-pass as dupRadar analysis
- **Biotype counting**: Aggregate read counts by gene biotype (e.g., `protein_coding`,
  `lncRNA`, `rRNA`) with MultiQC-compatible output files for bargraph and rRNA QC
- **Configurable biotype attribute**: `--biotype-attribute` CLI flag (or YAML config) to
  specify which GTF attribute to use for biotype grouping (default: `gene_biotype`,
  auto-detects `gene_type` for GENCODE GTFs)
- **YAML output control**: Nested `dupradar:` and `featurecounts:` config sections to
  individually enable/disable each output file (all enabled by default)
- **GTF attribute extraction**: Parser now extracts configurable extra attributes
  (e.g., `gene_biotype`, `gene_type`) stored in `Gene.attributes`
- **Plot improvements**: Plots are now square with proper aspect ratios and improved
  anti-aliasing

### Internal

- **Module reorganisation**: Flat module structure replaced with nested
  `src/rna/{dupradar, featurecounts, rseqc}/` layout
- **New dependencies**: `coitrees` (interval trees), `rayon` (parallelism),
  `rand`/`rand_chacha` (reproducible RNG for junction saturation)

### New Output Files

| File | Description |
|------|-------------|
| `{sample}.featureCounts.tsv` | featureCounts-format counts (Geneid, Chr, Start, End, Strand, Length, Count) |
| `{sample}.featureCounts.tsv.summary` | Assignment summary with all 14 featureCounts status categories |
| `{sample}.biotype_counts.tsv` | Per-biotype read counts |
| `{sample}.biotype_counts_mqc.tsv` | MultiQC bargraph of biotype counts |
| `{sample}.biotype_counts_rrna_mqc.tsv` | MultiQC general stats with rRNA percentage |
| `{sample}.bam_stat.txt` | Alignment statistics (bam-stat) |
| `{sample}.infer_experiment.txt` | Strandedness fractions (infer-experiment) |
| `{sample}.pos.DupRate.xls` | Position-based duplication histogram (read-duplication) |
| `{sample}.seq.DupRate.xls` | Sequence-based duplication histogram (read-duplication) |
| `{sample}.read_distribution.txt` | Per-region read distribution (read-distribution) |
| `{sample}.junction.xls` | Junction classifications (junction-annotation) |
| `{sample}.junction.bed` | Color-coded junction BED (junction-annotation) |
| `{sample}.junction_annotation.txt` | Junction summary (junction-annotation) |
| `{sample}.junctionSaturation_plot.r` | Saturation curve R script (junction-saturation) |
| `{sample}.inner_distance.txt` | Per-pair inner distances (inner-distance) |
| `{sample}.inner_distance_freq.txt` | Distance histogram (inner-distance) |

## [Version 0.1.0](https://github.com/ewels/RustQC/releases/tag/v0.1.0) - 2026-02-13

Initial release of RustQC — fast quality control tools for sequencing data, written in Rust. The first module (`rna`) is a reimplementation of [dupRadar](https://github.com/ssayols/dupRadar) for assessing PCR duplicate rates in RNA-Seq datasets.

### Breaking Changes

- **GTF is now a named flag** — use `--gtf <GTF>` (or `-g <GTF>`) instead of passing the GTF as a positional argument. This was required to support variadic BAM inputs.

### Features

- Drop-in replacement for R dupRadar with identical numerical output
- **Multiple BAM file support** — pass multiple BAM files as positional arguments and process them in parallel; duplicate BAM stems are detected and rejected
- Single-end and paired-end library support
- Strand-aware counting (unstranded, forward, reverse-stranded)
- SAM, BAM, and CRAM input support (auto-detected; `--reference` for CRAM)
- Multi-threaded alignment processing across chromosomes (`--threads`)
- 14-column duplication matrix, density scatter plot, boxplot, and expression histogram
- MultiQC-compatible output files
- YAML configuration for chromosome name mapping
- Cross-platform builds (Linux x86_64/aarch64, macOS x86_64/aarch64)
- Docker container at `ghcr.io/ewels/rustqc`

### Performance

- ~7x faster than R dupRadar single-threaded, ~27x faster with 10 threads
- Parallel alignment processing using rayon thread pool
- Multi-BAM parallelism — GTF is parsed once and shared, threads are distributed across BAM files
- Cache-oblivious interval trees (coitrees) for fast overlap queries
- Interned gene IDs and reusable buffers to minimise allocations
