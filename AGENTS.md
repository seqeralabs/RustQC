# AGENTS.md — RustQC

> Fast quality control tools for sequencing data, written in Rust. The
> **`rustqc rna`** command runs all RNA-Seq QC analyses in a single pass:
> - **dupRadar** duplicate rate analysis (reimplementation of
>   [dupRadar](https://github.com/ssayols/dupRadar))
> - **featureCounts**-compatible gene-level read counting with biotype summaries
> - **7 RSeQC tools** integrated as built-in analyses: bam_stat, infer_experiment,
>   read_duplication, read_distribution, junction_annotation, junction_saturation,
>   inner_distance (reimplementations of [RSeQC](https://rseqc.sourceforge.net/))
> - **TIN** (Transcript Integrity Number) analysis (reimplementation of RSeQC `tin.py`)
> - **preseq** library complexity extrapolation (reimplementation of
>   [preseq](https://github.com/smithlabcode/preseq) `lc_extrap`)
> - **samtools-compatible** flagstat, idxstats, and stats SN outputs
> - **Gene body coverage** profiling with Qualimap-compatible output
>
> Binary crate (`rustqc`), Rust edition 2021.

## Build / Lint / Test Commands

```bash
# Build
cargo build              # debug build
cargo build --release    # optimized release build (LTO, strip, opt-level 3)

# Format (enforced in CI — default rustfmt, no config file)
cargo fmt                # auto-format
cargo fmt --check        # check only (CI uses this)

# Lint (enforced in CI — default clippy, warnings are errors)
cargo clippy -- -D warnings

# Test — all unit + integration tests
cargo test               # debug mode
cargo test --release     # release mode (CI uses this)

# Run a single test by name (substring match)
cargo test test_dup_rate_calculation
cargo test test_dup_matrix_exact_match -- --nocapture   # with stdout

# Run only unit tests (skip integration tests)
cargo test --lib

# Run only integration tests
cargo test --test integration_test

# Run a specific integration test
cargo test --test integration_test test_intercept_slope_match
```

## Project Structure

```
src/
  main.rs           — Entry point, dispatches subcommands
  cli.rs            — CLI argument parsing (clap derive, subcommand structure)
  config.rs         — YAML configuration loading (serde), nested tool configs
  io.rs             — Shared I/O utilities (gzip-transparent file reading)
  gtf.rs            — GTF annotation file parser (with configurable attribute extraction)
  rna/
    mod.rs          — Re-exports dupradar, featurecounts, rseqc sub-modules
    dupradar/
      mod.rs        — Re-exports counting, dupmatrix, fitting, plots
      counting.rs   — BAM read counting engine (largest module)
      dupmatrix.rs  — Duplication matrix construction & TSV output
      fitting.rs    — Logistic regression via IRLS
      plots.rs      — Plot generation: density scatter, boxplot, histogram
    featurecounts/
      mod.rs        — Re-exports output
      output.rs     — featureCounts-format output & biotype counting
    genebody.rs     — Gene body coverage profiling, Qualimap-compatible output
    preseq.rs       — preseq lc_extrap library complexity extrapolation
    rseqc/
      mod.rs                — Re-exports all RSeQC modules + common helpers
      accumulators.rs       — Shared RSeQC accumulator infrastructure (read dispatch)
      common.rs             — Shared junction/intron extraction, from_genes/from_bed builders
      bam_stat.rs           — bam_stat.py reimplementation
      flagstat.rs           — samtools flagstat-compatible output
      idxstats.rs           — samtools idxstats-compatible output
      infer_experiment.rs   — infer_experiment.py reimplementation
      inner_distance.rs     — inner_distance.py reimplementation
      junction_annotation.rs — junction_annotation.py reimplementation
      junction_saturation.rs — junction_saturation.py reimplementation
      plots.rs              — RSeQC plot generation (duplication, junctions, etc.)
      read_distribution.rs  — read_distribution.py reimplementation
      read_duplication.rs   — read_duplication.py reimplementation
      stats.rs              — samtools stats SN-section compatible output
      tin.rs                — TIN (Transcript Integrity Number) analysis
tests/
  integration_test.rs  — 13 integration tests vs R dupRadar reference output
  data/                — Test BAM/GTF input files
  expected/            — R-generated reference outputs
  create_test_data.R   — R script to regenerate test data + references
```

Nested module structure — top-level modules (`cli`, `config`, `io`, `gtf`, `rna`) declared
in `main.rs`, no `lib.rs`. The `rna` module contains sub-modules for each tool group.
Inter-module access uses `crate::` paths (e.g., `use crate::rna::dupradar::counting::GeneCounts;`).

The CLI uses a single subcommand:
- `rustqc rna <BAM>... (--gtf <GTF> | --bed <BED>) [OPTIONS]`

The `--gtf` and `--bed` flags are **mutually exclusive** — provide one or the other:

- **`--gtf`** (recommended): Runs all analyses — dupRadar duplicate rate analysis,
  featureCounts-compatible gene counting, all 7 RSeQC-equivalent tools
  (bam_stat, infer_experiment, read_duplication, read_distribution,
  junction_annotation, junction_saturation, inner_distance), TIN analysis,
  preseq library complexity extrapolation, samtools-compatible outputs
  (flagstat, idxstats, stats), and gene body coverage profiling with
  Qualimap-compatible output. The GTF parser extracts transcript-level
  structure (exons + CDS features) to derive all data needed by every tool.
- **`--bed`**: Runs the 7 RSeQC tools + bam_stat + read_duplication + TIN +
  preseq + samtools outputs (flagstat, idxstats, stats). DupRadar,
  featureCounts, and gene body coverage are skipped (they require gene-level
  annotation with biotype/exon attributes that BED12 cannot provide).

Individual tools can be disabled via the YAML config file (each has an `enabled`
toggle). Tool-specific parameters (e.g., `--min-intron`, `--inner-distance-lower-bound`)
are available as CLI flags and in the config file.

Multiple BAM files can be passed as positional arguments and are processed sequentially.
Parallel processing is supported via `--threads`.

## Code Style

### Formatting

- **Default `rustfmt`** — no `rustfmt.toml` exists. Do not create one.
- 4-space indentation, ~100 char line width.
- Trailing commas on all multi-line constructs.
- Chained method calls break to new line with indent.

### Imports

Three groups (crate-internal, third-party, std), though blank-line separation
between groups is not strictly enforced. Each `use` is a single statement:

```rust
use crate::gtf::Gene;
use anyhow::{Context, Result};
use indexmap::IndexMap;
use log::{debug, info};
use std::collections::HashMap;
```

Localized `use` inside function bodies is acceptable for narrow imports
(e.g., `use std::io::Write;`).

### Naming

| Kind              | Convention             | Examples                                       |
|-------------------|------------------------|-------------------------------------------------|
| Types / Structs   | `CamelCase`            | `GeneCounts`, `DupMatrix`, `FitResult`          |
| Functions/Methods | `snake_case`           | `count_reads`, `build_index`, `format_float`    |
| Constants         | `SCREAMING_SNAKE_CASE` | `BAM_FDUP`, `DENSITY_COLORS`, `SCALE`           |
| Modules           | `snake_case`           | `dupmatrix`, `counting`, `bam_stat`, `inner_distance` |
| Variables/Fields  | `snake_case`           | `gene_counts`, `dup_rate_multi`, `is_dup`       |
| Type aliases      | `CamelCase`            | `MateBufferKey`                                 |

### Error Handling

- **`anyhow::Result<T>`** for all fallible functions. No custom error types.
- Propagate with `?` operator.
- Add context with `.context("msg")` or `.with_context(|| format!(...))`.
- Use `anyhow::bail!()` for early error returns.
- Use `anyhow::ensure!()` for precondition checks.
- **`unwrap()` / `expect()`** are restricted to test code only. In production code,
  `unwrap()` is acceptable only when a prior guard makes it provably safe (add a comment).

### Documentation

- **Every source file** starts with `//!` module doc comment (2-4 lines).
- **All public items** (structs, fields, functions, methods) get `///` doc comments.
- Complex functions include `# Arguments` and `# Returns` sections.
- Inline `//` comments explain complex logic, domain-specific behavior, and
  references to R equivalents.
- Long files use section dividers: `// ===` for major sections, `// ---` for sub-sections.

### Types and Derives

- `#[derive(Debug)]` on all structs.
- Add `Clone`, `Default`, `Deserialize` as needed — keep derives minimal.
- Public structs expose `pub` fields. Private helper structs keep fields private.
- Numeric conventions: `u64` for counts/positions, `f64` for metrics, `u8` for flags/strandedness.
- `IndexMap` when insertion order matters (gene ordering); `HashMap` for unordered lookups.

### Clippy

- Default clippy settings with `-D warnings` (deny all warnings).
- Targeted `#[allow(...)]` annotations are acceptable with justification:
  - `#[allow(dead_code)]` for fields kept for API completeness.
  - `#[allow(clippy::too_many_arguments)]` when refactoring would reduce clarity.

### Tests

- **Unit tests** co-located in each source file inside `#[cfg(test)] mod tests { use super::*; ... }`.
- **Integration tests** in `tests/integration_test.rs` — run the binary as a subprocess
  and compare output against R reference files in `tests/expected/`.
- Test naming: `test_<description_in_snake_case>`.
- Use `assert_eq!` with descriptive messages for exact comparisons.
- Use `assert!((val - expected).abs() < tolerance)` for float comparisons.
- No dev-dependencies — tests use only std + crate dependencies.

## CI Pipeline

GitHub Actions (`.github/workflows/ci.yml`) runs on push to `main` and all PRs:

1. **Test** — `cargo test --release` on Ubuntu and macOS
2. **Format** — `cargo fmt --check`
3. **Clippy** — `cargo clippy -- -D warnings`

All three must pass. Uses `dtolnay/rust-toolchain@stable` and `Swatinem/rust-cache@v2`.

## Key Dependencies

| Crate          | Purpose                          |
|----------------|----------------------------------|
| `clap` v4      | CLI argument parsing (derive)    |
| `rust-htslib`  | BAM file I/O (statically linked) |
| `plotters`     | Chart generation (PNG + SVG)     |
| `serde`        | YAML config deserialization      |
| `anyhow`       | Error handling                   |
| `log`          | Logging facade                   |
| `env_logger`   | Log output backend               |
| `indexmap`     | Insertion-order-preserving maps  |
| `coitrees`     | Cache-oblivious interval trees   |
| `rayon`        | Data parallelism                 |
| `rand` / `rand_chacha` | Reproducible random sampling |
| `flate2`       | Gzip decompression (annotation files) |

## Duplicate Marking Validation

RustQC verifies that BAM files have been processed by a duplicate-marking tool before
running the RNA duplication analysis. This is implemented in two layers:

1. **Pre-flight `@PG` header check** (`counting::verify_duplicates_marked`): Parses the
   BAM header text for `@PG` lines and checks the `ID:` and `PN:` fields against
   `KNOWN_DUP_MARKERS` (Picard MarkDuplicates, samblaster, sambamba, biobambam, etc.).
   Exits with `anyhow::bail!()` if no known tool is found.
2. **Post-hoc zero-duplicates check**: After counting, if `total_dup == 0` with
   `total_mapped > 0`, bails with an error — catches cases where the header is present
   but duplicates weren't actually flagged.

Both checks are skipped when `--skip-dup-check` is passed (stored as `RnaArgs.skip_dup_check`,
forwarded to `count_reads()` as the `skip_dup_check: bool` parameter).

## Notes for Agents

- The codebase is a pure binary crate with no library target.
- Release builds use aggressive optimization (`lto = true`, `codegen-units = 1`, `strip = true`).
- Test data is generated by `tests/create_test_data.R` — do not modify `tests/expected/` by hand.
- Float output formatting must match R's behavior (15 significant digits, "NA" for NaN, trailing-zero trimming).
- The pipeline processes BAM files which can be very large — performance matters.
- System dependencies needed for building: cmake, zlib, bz2, lzma, curl, ssl, clang (for `rust-htslib`).
- **IMPORTANT:** When benchmarks are re-run, verify that all results referenced in both `benchmark/README.md`
  and the top-level `README.md` are updated to reflect the new numbers (timings, percentages, etc.).
  These documents must always accurately reflect the latest benchmark data.
- The YAML config nests output toggles under `dupradar:` and `featurecounts:` keys.
  Each output file can be individually enabled/disabled (all default to `true`).
- The YAML config also has sections for each of the 7 RSeQC tools (`bam_stat:`,
  `infer_experiment:`, `read_duplication:`, `read_distribution:`, `junction_annotation:`,
  `junction_saturation:`, `inner_distance:`). Each has an `enabled: bool` toggle
  and tool-specific parameter overrides. CLI flags take precedence over config values.
- The YAML config also has sections for `tin:`, `preseq:`, `genebody_coverage:`,
  `flagstat:`, `idxstats:`, and `samtools_stats:`. Each has an `enabled: bool` toggle.
  Preseq has additional parameters: `max_extrap`, `step_size`, `n_bootstraps`,
  `confidence_level`, `seed`, `max_terms`, `defects`. TIN has `sample_size` and
  `min_coverage`.
- The `featurecounts.rs` module produces featureCounts-compatible output files (counts TSV,
  summary, biotype counts, and MultiQC files). These are generated in the same pass as
  the dupRadar analysis — no separate featureCounts run is needed.
- The genebody module (`genebody.rs`) produces Qualimap-compatible output files
  (`coverage_profile_along_genes_(total).txt` and `rnaseq_qc_results.txt`). These
  files are written to a `qualimap/` subdirectory and are parseable by MultiQC as
  Qualimap rnaseq reports. Gene body coverage only runs in GTF mode.
- The GTF parser (`gtf.rs`) accepts extra attribute names to extract (e.g., `gene_biotype`,
  `gene_type`) and stores them in `Gene.attributes`. The biotype attribute name is
  configurable via `--biotype-attribute` CLI flag or `featurecounts.biotype_attribute`
  in the YAML config.
- GENCODE GTFs use `gene_type` while Ensembl GTFs use `gene_biotype`. The tool
  auto-detects which is present and falls back gracefully with a warning if neither is found.
- Both GTF and BED annotation files can be provided plain or gzip-compressed (`.gz`).
  Detection is based on file magic bytes, not the file extension. The shared `io` module
  (`src/io.rs`) provides `open_reader()` and `read_to_string()` helpers that transparently
  handle both formats. All annotation file readers use these helpers.
