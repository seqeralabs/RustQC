# RustQC Changelog

## [Version 0.1.1](https://github.com/seqeralabs/RustQC/releases/tag/v0.1.1) - 2026-04-02

### Bug fixes

- Fix featureCounts summary to use gene-level stats; add biotype summary (#66)
- Fix inner_distance histogram to include overflow bucket in bulk cutoff loop (#67)

### Other changes

- Add crates.io publishing to release workflow (#62)
- Documentation fixes (#70)

## [Version 0.1.0](https://github.com/seqeralabs/RustQC/releases/tag/v0.1.0) - 2026-04-01

Initial release of RustQC -- fast quality control tools for sequencing data, written in Rust.

A single `rustqc rna` command runs 15 QC analyses in one pass over the BAM file, producing output that is format- and numerically identical to the upstream tools and fully compatible with [MultiQC](https://multiqc.info/).

### Tools

- **dupRadar** -- PCR duplicate rate vs. expression analysis with density scatter plots, boxplots, and expression histograms. 14-column duplication matrix with logistic-regression model fitting matching the [R dupRadar](https://github.com/ssayols/dupRadar) package.
- **featureCounts** -- Gene-level read counting with assignment summary, compatible with [Subread featureCounts](http://subread.sourceforge.net/). Includes per-biotype read counting and MultiQC integration.
- **RSeQC** (8 tools) -- [RSeQC](https://rseqc.sourceforge.net/)-compatible implementations of bam_stat, infer_experiment, read_duplication, read_distribution, junction_annotation, junction_saturation, inner_distance, and TIN (Transcript Integrity Number). Includes native plot generation (PNG + SVG) with no R dependency.
- **preseq** -- Library complexity extrapolation (`lc_extrap`) matching [preseq](http://smithlabresearch.org/software/preseq/) v3.2.0, including C++ RNG FFI for reproducible bootstrap sampling.
- **Qualimap rnaseq** -- Gene body coverage profiling, read genomic origin, junction analysis, and transcript coverage bias matching [Qualimap](http://qualimap.conesalab.org/).
- **samtools** -- flagstat, idxstats, and full stats output (SN section + all histogram sections) matching [samtools](http://www.htslib.org/).

### Features

- Single static binary with no runtime dependencies
- SAM, BAM, and CRAM input support (auto-detected)
- Multi-threaded alignment processing with htslib thread pool
- GTF annotation support (gzip-compressed files accepted)
- YAML configuration for output control, chromosome name mapping, and tool parameters
- Multiple BAM file support via positional arguments
- `--sample-name` flag to override BAM-derived sample name in output filenames
- Per-tool seed flags (`--tin-seed`, `--junction-saturation-seed`, `--preseq-seed`) for reproducible results
- Docker container at `ghcr.io/seqeralabs/rustqc`
- Cross-platform builds (Linux x86_64/aarch64, macOS x86_64/aarch64)
