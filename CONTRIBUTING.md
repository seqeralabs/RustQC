# Contributing to RustQC

Thanks for your interest in contributing to RustQC! This document covers how to get started.

## Getting started

### Prerequisites

You need a working Rust toolchain (stable) and the following system libraries for building `rust-htslib`:

- cmake, clang
- zlib, libbz2, liblzma
- libcurl, libssl
- libfontconfig (for plot rendering)

On Ubuntu/Debian:

```bash
sudo apt-get install cmake clang zlib1g-dev libbz2-dev liblzma-dev \
  libcurl4-openssl-dev libssl-dev libfontconfig1-dev
```

On macOS (with Homebrew):

```bash
brew install cmake xz
```

### Building

```bash
cargo build              # debug build
cargo build --release    # optimized release build
```

### Running tests

```bash
cargo test               # all tests (unit + integration)
cargo test --lib         # unit tests only
cargo test --test integration_test  # integration tests only
cargo test test_name     # single test by name substring
```

### Linting and formatting

Both are enforced in CI and must pass before merging:

```bash
cargo fmt                # auto-format code
cargo fmt --check        # check formatting without changes
cargo clippy -- -D warnings  # lint with warnings as errors
```

### Pre-commit hooks

The repository includes a [pre-commit](https://pre-commit.com/) configuration
that runs formatting and linting checks before each commit.
We recommend [prek](https://github.com/j178/prek) (a fast Rust rewrite of pre-commit):

```bash
brew install prek   # or: cargo install prek
prek install
```

This runs `cargo fmt --check` and `cargo clippy -- -D warnings` automatically
on every commit, plus basic file hygiene checks. If a hook fails, fix the issue
and re-stage your changes.

To run all hooks manually:

```bash
prek run --all-files
```

## Project structure

RustQC is a **binary crate** with a nested module structure. All top-level modules are declared in `main.rs` (no `lib.rs`).

```
src/
  main.rs                        Entry point, dispatches subcommands
  cli.rs                         CLI argument parsing (clap derive)
  config.rs                      YAML configuration loading (serde)
  io.rs                          Shared I/O utilities (gzip-transparent file reading)
  gtf.rs                         GTF annotation file parser
  rna/
    mod.rs                       Re-exports all submodules (bam_flags, cpp_rng, dupradar,
                                   featurecounts, preseq, qualimap, rseqc)
    bam_flags.rs                 BAM flag constants
    cpp_rng.rs                   C++ RNG FFI shim for preseq bootstrap reproducibility
    dupradar/
      mod.rs                     Re-exports counting, dupmatrix, fitting, plots
      counting.rs                BAM read counting engine (4-mode counting)
      dupmatrix.rs               Duplication matrix construction and TSV output
      fitting.rs                 Logistic regression via IRLS
      plots.rs                   Plot generation (density scatter, boxplot, histogram)
    featurecounts/
      mod.rs                     Re-exports output
      output.rs                  featureCounts-format output and biotype counting
    preseq.rs                    preseq lc_extrap library complexity extrapolation
    qualimap/
      mod.rs                     Re-exports accumulator, coverage, index, output, plots, report
      accumulator.rs             Gene body coverage accumulator
      coverage.rs                Coverage computation logic
      index.rs                   Transcript interval index for coverage lookups
      output.rs                  Qualimap-compatible output file generation
      plots.rs                   Gene body coverage plot generation
      report.rs                  HTML/text report generation
    rseqc/
      mod.rs                     Re-exports all RSeQC modules + common helpers
      accumulators.rs            Shared RSeQC accumulator infrastructure (read dispatch)
      common.rs                  Shared junction/intron extraction, from_genes builders
      plots.rs                   RSeQC native plot generation (PNG/SVG)
      bam_stat.rs                Basic BAM alignment statistics
      flagstat.rs                samtools flagstat-compatible output
      idxstats.rs                samtools idxstats-compatible output
      infer_experiment.rs        Library strandedness inference
      inner_distance.rs          Inner distance for paired-end reads
      junction_annotation.rs     Splice junction classification
      junction_saturation.rs     Junction saturation analysis
      read_distribution.rs       Read distribution across genomic features
      read_duplication.rs        Position- and sequence-based duplication histograms
      stats.rs                   samtools stats full output (SN + all histogram sections)
      tin.rs                     TIN (Transcript Integrity Number) analysis
tests/
  integration_test.rs            Integration tests vs R reference output
  data/                          Test BAM/GTF input files
  expected/                      R-generated reference outputs
  create_test_data.R             R script to regenerate test data
```

Inter-module access uses `crate::` paths (e.g., `use crate::gtf::Gene;`).

## Making changes

1. Fork the repository and create a branch from `main`.
2. Make your changes, following the code style below.
3. Add or update tests as appropriate.
4. Run `cargo fmt`, `cargo clippy -- -D warnings`, and `cargo test` locally (or use `prek run --all-files` for fmt + clippy).
5. Open a pull request against `main`.

## Code style

The project uses **default rustfmt** (no config file) and **default clippy** with `-D warnings`. See [AGENTS.md](AGENTS.md) for detailed style guidelines, but the key points are:

- **Error handling**: Use `anyhow::Result` with `.context()` / `bail!()` / `ensure!()`. No custom error types. No `unwrap()` or `expect()` in production code.
- **Documentation**: Every source file needs `//!` module docs. All public items need `///` doc comments.
- **Imports**: Three groups -- crate-internal, third-party, std -- each as single `use` statements.
- **Naming**: Standard Rust conventions -- `CamelCase` types, `snake_case` functions/variables, `SCREAMING_SNAKE_CASE` constants.
- **Derives**: `#[derive(Debug)]` on all structs. Add `Clone`, `Default`, `Deserialize` only as needed.
- **Collections**: `IndexMap` when insertion order matters (gene ordering), `HashMap` otherwise.

## Test data

Test reference data in `tests/expected/` is generated by `tests/create_test_data.R` using R's dupRadar package. Do not edit these files by hand. If test inputs change, re-run the R script to regenerate references.

Benchmark reference outputs and reproduction scripts are in the [RustQC-benchmarks](https://github.com/seqeralabs/RustQC-benchmarks) Nextflow pipeline.

## Numerical accuracy

RustQC aims to match the output of the original tools exactly:

- **dupRadar**: Float formatting follows R conventions (15 significant digits, `NA` for NaN, trailing-zero trimming). Integration tests verify this against R reference output with tight tolerances.
- **RSeQC**: Output formats replicate Python RSeQC conventions. Validation is done by comparing against reference outputs from RSeQC 5.0.4.

## CI

All pull requests must pass three checks:

1. `cargo test --release` on Ubuntu and macOS
2. `cargo fmt --check`
3. `cargo clippy -- -D warnings`

## Reporting issues

Please open a [GitHub issue](https://github.com/seqeralabs/RustQC/issues) with:

- What you expected to happen
- What actually happened
- Steps to reproduce (SAM/BAM/CRAM/GTF files or minimal examples if possible)
- Your OS and Rust version (`rustc --version`)

## License

By contributing, you agree that your contributions will be licensed under the [GNU General Public License v3.0 or later](LICENSE).
