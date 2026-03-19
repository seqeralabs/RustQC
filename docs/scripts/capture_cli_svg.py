#!/usr/bin/env python3
"""Capture RustQC CLI output as SVG files using Rich.

Generates two SVGs for documentation:
1. Full CLI output from a normal run (no warnings)
2. Strandedness mismatch warning (trimmed to just the warning box)

Uses short symlinks in a temp directory so the file paths printed by
RustQC are concise and look good in the SVG screenshots.

Usage:
    uv run --with rich docs/scripts/capture_cli_svg.py
"""

import os
import subprocess
import sys
import tempfile

from pathlib import Path
from rich.console import Console
from rich.text import Text

# Paths — adjust these if your test data is elsewhere
RUSTQC = Path("./target/release/rustqc").resolve()
GTF = Path("../RustQC-benchmarks/test-data/rna/small/chr6.gtf.gz").resolve()
BAM = Path("../RustQC-benchmarks/test-data/rna/small/test.bam").resolve()

TITLE = "rustqc rna"


def run_rustqc(stranded: str) -> str:
    """Run RustQC and return combined stdout+stderr as a string.

    Creates a temporary directory with short symlinks so that the paths
    printed by RustQC in its output are concise.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create short symlinks
        gtf_link = os.path.join(tmpdir, "genes.gtf.gz")
        bam_link = os.path.join(tmpdir, "sample.bam")
        bai_src = str(BAM) + ".bai"
        bai_link = os.path.join(tmpdir, "sample.bam.bai")
        outdir = os.path.join(tmpdir, "results")

        os.symlink(str(GTF), gtf_link)
        os.symlink(str(BAM), bam_link)
        if os.path.exists(bai_src):
            os.symlink(bai_src, bai_link)

        cmd = [
            str(RUSTQC),
            "rna",
            "--stranded",
            stranded,
            "-p",
            "--gtf",
            gtf_link,
            bam_link,
            "--outdir",
            outdir,
        ]
        env = {
            "CLICOLOR_FORCE": "1",
            "PATH": "/usr/bin:/bin:/usr/local/bin",
        }
        result = subprocess.run(
            cmd,
            capture_output=True,
            env=env,
            timeout=120,
        )
        output = (result.stdout + result.stderr).decode("utf-8", errors="replace")

    # Replace the temp dir path with a clean short path for display
    output = output.replace(tmpdir, "/data")
    return output


def write_svg(output: str, path: str, title: str = TITLE):
    """Render ANSI output as SVG via Rich and write to file."""
    console = Console(record=True, width=100)
    text = Text.from_ansi(output)
    console.print(text)
    svg = console.export_svg(title=title)
    with open(path, "w") as f:
        f.write(svg)
    print(f"SVG written to {path}")


def main():
    # 1. Normal run — full output, correct strandedness (reverse for this data)
    print("Capturing normal run...")
    output_normal = run_rustqc("reverse")
    write_svg(
        output_normal,
        "docs/public/examples/cli-output.svg",
    )

    # 2. Mismatch run — wrong strandedness to trigger warning
    print("Capturing strandedness mismatch warning...")
    output_mismatch = run_rustqc("forward")

    # Trim to only show the warning box and lines after it
    lines = output_mismatch.splitlines(keepends=True)
    start = 0
    for i, line in enumerate(lines):
        if "Warning" in line and "\u256d" in line:
            # Include one blank line before the box for spacing
            start = max(0, i - 1)
            break
    output_mismatch = "".join(lines[start:]).rstrip() + "\n"
    write_svg(
        output_mismatch,
        "docs/public/examples/cli-strandedness-warning.svg",
    )


if __name__ == "__main__":
    main()
