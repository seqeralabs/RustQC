# RustQC Changelog

## Unreleased

### Features

- New `rustqc dna` subcommand for DNA (WGS) quality control: depth of coverage
  with mosdepth-compatible outputs, samtools-compatible stats, flagstat and
  idxstats, and preseq library complexity, all in a single pass over the
  alignment with one worker per contig, plus Picard-compatible
  CollectWgsMetrics, CollectInsertSizeMetrics and CollectGcBiasMetrics.
  Passing `--targets` switches on targeted mode and Picard-compatible
  CollectHsMetrics. Qualimap-compatible `bamqc` output rounds it out, with
  `genome_results.txt`, the raw data tables and an HTML summary. Validated
  against mosdepth 0.3.14, samtools 1.24, Picard 3.4.0 and Qualimap 2.3.

- New `rustqc protein sequence` mode: protein FASTA quality control, with
  length statistics reproducing `seqkit stats -a` exactly, plus amino acid
  composition and defect detection (internal stop codons, non-standard
  residues, duplicate sequences and identifiers) that seqkit does not report.
- New `rustqc protein spectra` mode: mass spectrometry run QC from mzML, with
  per-level spectrum and peak counts, total ion current, retention time range
  and precursor charge distribution. Behind the `proteomics` cargo feature,
  which is on by default; `--no-default-features` drops it along with mzdata.

### Changed

- Internal: assay-agnostic analyses (BAM flag helpers, read-level statistics,
  the samtools stats/flagstat/idxstats writers, preseq) moved from `rna` to a
  new `common` module. The old `rustqc::rna::...` paths still resolve through
  re-exports, so this is not a breaking change for library users.

## [Version 0.2.1](https://github.com/seqeralabs/RustQC/releases/tag/v0.2.1) - 2026-04-09

### Bug fixes

- Fix SIMD builds by scoping rustflags to target triple with explicit `--target` (#90, #91)

### Other changes

- Trigger releases via workflow dispatch instead of tag push (#92)

## [Version 0.2.0](https://github.com/seqeralabs/RustQC/releases/tag/v0.2.0) - 2026-04-09

### Features

- Ship SIMD-optimized binaries with CPU detection and upgrade hints (#81)
- Write `CITATIONS.md` with upstream tool versions (#87)
- Add XDG config discovery and deep-merge support (#88)

### Bug fixes

- Replace header-based duplicate check with flag-based detection (#84)
- Use `.log` extension for junction_annotation output (#80)

### Other changes

- Bump docker/login-action from 4.0.0 to 4.1.0 (#78)
- Bump vite from 7.3.1 to 7.3.2 in docs (#77)
- Bump defu from 6.1.4 to 6.1.6 in docs (#74)

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
