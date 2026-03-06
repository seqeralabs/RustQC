# RustQC Changelog

## [Unreleased]

### Added

- **TIN (Transcript Integrity Number) calculation**: RSeQC `tin.py`-compatible transcript integrity analysis
- **preseq library complexity estimation**: `lc_extrap` reimplementation for extrapolating library complexity curves
- **samtools flagstat compatible output**: Alignment flag statistics matching samtools flagstat format
- **samtools idxstats compatible output**: Per-reference alignment statistics matching samtools idxstats format
- **samtools stats compatible output**: Full SN-section output with ALL histogram sections (FFQ, LFQ, GCF, GCL, GCC, GCT, FBC, LBC, FTC, LTC, IS, RL, FRL, LRL, MAPQ, ID, IC)
- **Qualimap rnaseq compatible output**: Gene body coverage profiling with coverage, bias, and counts output matching Qualimap format
- **Simultaneous `--gtf` and `--bed` support**: GTF used for dupRadar/featureCounts, BED used for read_distribution when both are provided
- **`--skip-preseq` flag**: Disable preseq library complexity estimation
- **`--flat-output` flag**: Flat output directory structure (all files in one directory)
- **`--skip-dup-check` flag**: Skip duplicate-marking verification checks
- **Gzip-compressed annotation support**: GTF and BED annotation files can now be provided as `.gz` files. Compression is detected automatically by inspecting the file header (magic bytes), so both plain and gzip-compressed files work transparently with `--gtf` and `--bed`
- **Preseq CLI parameters**: `--seed`, `--n-bootstraps`, `--confidence-level`, `--max-extrap`, `--step-size`
- **htslib threading support**: Multi-threaded BAM I/O via htslib thread pool
- **coitrees interval tree**: Fast overlap queries using cache-oblivious interval trees

### Changed

- **preseq**: Rewritten PE fragment counting to match preseq v3.2.0 — mates are now merged by read name into genomic intervals `(chrom, min_start, max_end)` instead of counting individual reads with `(tid, pos, mtid, mpos)` keys. Fragment counts now match upstream exactly.
- **preseq**: Bootstrap RNG now uses C++ FFI shim (`std::mt19937` + `std::binomial_distribution`) instead of Rust `rand_mt`/`rand_distr`, producing byte-identical bootstrap samples when compiled with the same C++ standard library as upstream preseq
- **preseq**: Added `max_segment_length` config option (default 100M, matching nf-core's `-seg_len` setting) to control maximum merged PE fragment length

### Fixed

- **preseq**: Use `tid>=0` filter matching upstream (not flag-based filtering)
- **preseq**: Default seed changed to 1 matching upstream default
- **preseq**: Discard bootstrap on failed degree selection instead of using max_terms
- **preseq**: Remove incorrect `(j+1)/N` scaling in defects mode PS coefficients
- **preseq**: First output row uses integers (`0\t0\t0\t0`) matching upstream
- **preseq**: `search_max_val` always 100.0 matching upstream
- **preseq**: `interpolate()` uses dense histogram with `i32` accumulator matching upstream's `std::accumulate(..., 0)` integer truncation
- **preseq**: `evaluate_cf()` returns `f64` unconditionally (no validity checks), matching upstream
- **preseq**: `extrapolate()` returns `f64` unconditionally; removed `target<=n` guard and finite/non-negative checks
- **preseq**: `compute_curve()` no longer breaks early on extrapolation failure; no `.max(0.0)` clamping
- **preseq**: Added `check_yield_estimates_stability()` matching upstream's stability validation (non-negative, finite, monotonically increasing, concave)
- **preseq**: Bootstrap loop uses `check_yield_estimates_stability` on full yield vector instead of simple `all(finite && >=0)` check
- **qualimap**: Aligned pairs count and 5'-3' bias now match upstream Qualimap exactly
- **dupradar**: RPKM calculation now matches upstream R dupRadar exactly
- **read_duplication**: Added `DupRate_plot.r` R script output matching RSeQC format
- **samtools stats**: Full output with all histogram sections (FFQ, LFQ, GCF, GCL, GCC, GCT, FBC, LBC, FTC, LTC, IS, RL, FRL, LRL, MAPQ, ID, IC, COV, CHK) matching upstream format
- **samtools stats**: Fixed cycle numbering to be 1-based matching samtools
- **samtools stats**: Fixed FBC/LBC sections to output percentages (not raw counts)
- **samtools stats**: Fixed GCC columns to match samtools column order
- **samtools stats**: Fixed GCF/GCL sections to use non-uniform binning matching samtools
- **samtools stats**: Fixed COV section depth distribution
- **samtools stats**: Added CHK section with CRC32 checksums
- **samtools stats**: CIGAR operation caching for performance
- **samtools stats**: MQ0 counts primary mapped reads only
- **samtools stats**: Insert size uses paired+both-mapped filter (not proper-pair), 99th percentile truncation, capped at 8000
- **samtools stats**: Orientation uses upstream `is_fst*pos_fst` algorithm, counts both mates, divides by 2
- **samtools stats**: CIGAR mapped bases includes insertions (Ins)
- **samtools stats**: Average length uses `round()` matching `printf %.0f`
- **samtools stats**: Average quality computes per-base average (not per-read average of averages)
- **samtools stats**: Filtered sequences always 0 (no CLI filtering)
- **infer_experiment**: Use transcript-level intervals instead of gene-level spans
- **infer_experiment**: Include supplementary alignments matching upstream
- **infer_experiment**: Add soft clips to query length matching pysam qlen
- **infer_experiment**: Fix parallel over-sampling with proportional scale-down
- **infer_experiment**: Output 2 blank lines before header matching upstream format
- **qualimap**: Include NaN in 5'-3' bias ratio array matching Java behavior
- **qualimap**: Threshold-based transcript selection matching `pickTranscripts()` tie handling
- **dupradar**: RPKM denominator excludes secondary alignments for unique-mapper mode
- **dupradar**: Count NaN dupRate genes in `nRegionsDuplication` matching R behavior
- **dupradar**: Classify unmatched mates as `Unassigned_Singleton`
- **dupradar**: Logistic fit uses deviance-based convergence and per-observation starting values matching R's `glm`

## [Version 0.1.0](https://github.com/seqeralabs/RustQC/releases/tag/v0.1.0) - 2026-02-13

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

Seven [RSeQC](https://rseqc.sourceforge.net/) tools integrated into the
`rustqc rna` command, running automatically in the same single-pass analysis.
Output formats are compatible with the Python originals. Individual tools can
be disabled via the YAML config file.

- bam_stat -- alignment statistics
- infer_experiment -- library strandedness inference
- read_duplication -- position-based and sequence-based duplication histograms with duplication rate plot (PNG + SVG)
- read_distribution -- read distribution across genomic features
- junction_annotation -- splice junction classification with splice event and junction pie charts (PNG + SVG)
- junction_saturation -- junction saturation curves with saturation plot (PNG + SVG)
- inner_distance -- insert size estimation for paired-end reads with histogram and density plot (PNG + SVG)
- Native plot generation for all applicable RSeQC tools -- no R required

### General

- SAM, BAM, and CRAM input support (auto-detected)
- Multiple BAM file support -- pass multiple files as positional arguments
- YAML configuration for output control and chromosome name mapping
- Cross-platform builds (Linux x86_64/aarch64, macOS x86_64/aarch64)
- Docker container at `ghcr.io/seqeralabs/rustqc`
