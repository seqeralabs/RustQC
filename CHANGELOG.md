# RustQC Changelog 🧬 🦀

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
