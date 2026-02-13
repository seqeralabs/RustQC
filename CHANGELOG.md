# RustQC Changelog 🧬 🦀

## [Unreleased]

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

### New Output Files

| File | Description |
|------|-------------|
| `{sample}.featureCounts.tsv` | featureCounts-format counts (Geneid, Chr, Start, End, Strand, Length, Count) |
| `{sample}.featureCounts.tsv.summary` | Assignment summary with all 14 featureCounts status categories |
| `{sample}.biotype_counts.tsv` | Per-biotype read counts |
| `{sample}.biotype_counts_mqc.tsv` | MultiQC bargraph of biotype counts |
| `{sample}.biotype_counts_rrna_mqc.tsv` | MultiQC general stats with rRNA percentage |

## [Version 0.1.0](https://github.com/ewels/RustQC/releases/tag/v0.1.0) - 2026-02-13

Initial release of RustQC — fast quality control tools for sequencing data, written in Rust. The first module (`rna`) is a reimplementation of [dupRadar](https://github.com/ssayols/dupRadar) for assessing PCR duplicate rates in RNA-Seq datasets.

### Features

- Drop-in replacement for R dupRadar with identical numerical output
- Single-end and paired-end library support
- Strand-aware counting (unstranded, forward, reverse-stranded)
- Multi-threaded BAM processing across chromosomes (`--threads`)
- 14-column duplication matrix, density scatter plot, boxplot, and expression histogram
- MultiQC-compatible output files
- YAML configuration for chromosome name mapping
- Cross-platform builds (Linux x86_64/aarch64, macOS x86_64/aarch64)
- Docker container at `ghcr.io/ewels/rustqc`

### Performance

- ~7x faster than R dupRadar single-threaded, ~27x faster with 10 threads
- Parallel BAM processing using rayon thread pool
- Cache-oblivious interval trees (coitrees) for fast overlap queries
- Interned gene IDs and reusable buffers to minimise allocations
