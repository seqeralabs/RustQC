# RustQC Changelog

## [Version 0.1.0](https://github.com/ewels/RustQC/releases/tag/v0.1.0) - 2026-02-13

Initial release of RustQC -- fast quality control tools for sequencing data, written in Rust.

### dupRadar

Reimplementation of [dupRadar](https://github.com/ssayols/dupRadar) for assessing
PCR duplicate rates in RNA-Seq datasets, with identical numerical output to the R original.

- 14-column duplication matrix with density scatter plot, boxplot, and expression histogram
- Single-end and paired-end library support, strand-aware counting
- Multi-threaded alignment processing across chromosomes (`--threads`)

### featureCounts

- featureCounts-compatible gene counts output (TSV + summary) generated in the same single-pass analysis as dupRadar
- Per-biotype read counting and duplication rate analysis
- MultiQC-compatible output files for biotype counts

### RSeQC

Seven subcommands reimplementing the most widely-used
[RSeQC](https://rseqc.sourceforge.net/) Python tools as fast, single-pass Rust equivalents.
Output formats are compatible with the Python originals.

- `bam-stat` -- alignment statistics
- `infer-experiment` -- library strandedness inference
- `read-duplication` -- position-based and sequence-based duplication histograms
- `read-distribution` -- read distribution across genomic features
- `junction-annotation` -- splice junction classification
- `junction-saturation` -- junction saturation curves
- `inner-distance` -- insert size estimation for paired-end reads

### General

- SAM, BAM, and CRAM input support (auto-detected)
- Multiple BAM file support -- pass multiple files as positional arguments
- YAML configuration for output control and chromosome name mapping
- Cross-platform builds (Linux x86_64/aarch64, macOS x86_64/aarch64)
- Docker container at `ghcr.io/ewels/rustqc`
