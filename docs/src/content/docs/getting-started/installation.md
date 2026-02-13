---
title: Installation
description: How to install RustQC and its system dependencies, build from source, and verify your installation.
---

RustQC is built from source using the Rust toolchain. It statically links against
[htslib](https://github.com/samtools/htslib) for SAM/BAM/CRAM I/O, so a few system
libraries are required at compile time.

## Prerequisites

### Rust toolchain

Install Rust via [rustup](https://rustup.rs/) if you do not already have it:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

### System dependencies

RustQC depends on [rust-htslib](https://github.com/rust-bio/rust-htslib), which
requires the following system libraries for building:

**macOS** (Homebrew):

```bash
brew install cmake zlib bzip2 xz curl openssl
```

**Ubuntu / Debian**:

```bash
sudo apt install cmake zlib1g-dev libbz2-dev liblzma-dev \
  libcurl4-openssl-dev libssl-dev clang
```

**Fedora / RHEL**:

```bash
sudo dnf install cmake zlib-devel bzip2-devel xz-devel \
  libcurl-devel openssl-devel clang
```

In summary, you need: **cmake**, **zlib**, **bz2**, **lzma**, **curl**, **ssl**, and **clang**.

## Build from source

Clone the repository and build in release mode:

```bash
git clone https://github.com/ewels/RustQC.git
cd RustQC
cargo build --release
```

The compiled binary is at `target/release/rustqc`. You can copy it to a directory
on your `PATH`, or run it directly:

```bash
./target/release/rustqc --version
```

Release builds use link-time optimization (LTO) and symbol stripping for maximum
performance and a small binary size.

## Quick start

Run RustQC on a single BAM file with a GTF annotation:

```bash
rustqc rna sample.bam --gtf genes.gtf -p -o output/
```

Process multiple BAM files in parallel:

```bash
rustqc rna sample1.bam sample2.bam sample3.bam \
  --gtf genes.gtf -p -t 8 -o output/
```

Use CRAM input with a reference FASTA:

```bash
rustqc rna sample.cram --gtf genes.gtf -p -r reference.fa -o output/
```

RustQC accepts **SAM**, **BAM**, and **CRAM** alignment files. Multiple input files
can be passed as positional arguments and are processed in parallel when threads
are available.

## Input requirements

Your alignment files **must** have PCR duplicates **marked** (SAM flag `0x400`) but
**not removed**. Use one of these tools before running RustQC:

- [Picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)
- [samblaster](https://github.com/GregoryFaust/samblaster)
- [sambamba markdup](https://lomereiter.github.io/sambamba/)
- [biobambam bammarkduplicates](https://gitlab.com/german.tischler/biobambam2)

RustQC automatically checks the BAM header for duplicate-marking tool signatures
and exits with an error if none are found. Use `--skip-dup-check` to bypass this
validation if needed.

## Next steps

- [CLI Reference](/RustQC/usage/cli-reference/) -- full list of options and flags
- [Output Files](/RustQC/usage/output-files/) -- description of every output file
- [Interpreting Plots](/RustQC/guide/interpreting-plots/) -- how to read the generated plots
- [Configuration](/RustQC/usage/configuration/) -- YAML config file options
