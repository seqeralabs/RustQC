#!/usr/bin/env python3
"""
Benchmark runner for RustQC vs. upstream tools.

Runs benchmarks for dupRadar (R), featureCounts (Subread), RSeQC (Python),
preseq, samtools, Qualimap, and RustQC. Profiles CPU and memory usage with
psrecord.

Usage:
    python3 benchmark/run_benchmarks.py --all
    python3 benchmark/run_benchmarks.py --rustqc --dupradar
    python3 benchmark/run_benchmarks.py --all --large-only

Requirements: Python 3.8+, psrecord, matplotlib, Docker.
Run from the repository root.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# Docker images (Wave-built)
# ---------------------------------------------------------------------------
DUPRADAR_IMG = (
    "wave.seqera.io/wt/66d7208af012/wave/build:bioconductor-dupradar--57fe475bbf67ac88"
)
SUBREAD_IMG = "wave.seqera.io/wt/30bc38d2ad74/wave/build:subread--d7eae49d8c0e1b66"
RSEQC_IMG = "wave.seqera.io/wt/ea3e9f972b6e/wave/build:rseqc-5.0.4--14c99cde3bff8d57"
RBASE_IMG = "wave.seqera.io/wt/d63a19ed5f77/wave/build:r-base--6fa87e5d876579f1"
PRESEQ_IMG = "quay.io/biocontainers/preseq:3.2.0--hdcf5f25_6"
SAMTOOLS_IMG = "quay.io/biocontainers/samtools:1.21--h50ea8bc_0"
QUALIMAP_IMG = "quay.io/biocontainers/qualimap:2.3--hdfd78af_0"
DOCKER_PLATFORM = "linux/amd64"

# ---------------------------------------------------------------------------
# Datasets
# ---------------------------------------------------------------------------
DATASETS = {
    "small": dict(
        bam="benchmark/input/small/test.bam",
        gtf="benchmark/input/small/chr6.gtf",
        bed="benchmark/input/small/chr6.bed",
        paired=True,
        stranded=0,
        skip_dup_check=True,
        config=None,
        biotype_attribute=None,
    ),
    "large": dict(
        bam="benchmark/input/large/GM12878_REP1.markdup.sorted.bam",
        gtf="benchmark/input/large/genes.gtf",
        bed="benchmark/input/large/genes.bed",
        paired=True,
        stranded=0,
        skip_dup_check=False,
        config="benchmark/input/large/config.yaml",
        biotype_attribute="gene_type",
    ),
}

RSEQC_TOOLS = [
    "bam_stat",
    "infer_experiment",
    "read_duplication",
    "read_distribution",
    "junction_annotation",
    "junction_saturation",
    "inner_distance",
]

# RSeQC tools that produce R plot scripts, and the script filenames they generate.
# Paths are relative to the output directory (prefix is the -o value).
RSEQC_PLOT_SCRIPTS = {
    "read_duplication": ["read_duplication.DupRate_plot.r"],
    "junction_annotation": ["junction_annotation.junction_plot.r"],
    "junction_saturation": ["junction_saturation.junctionSaturation_plot.r"],
    "inner_distance": ["inner_distance.inner_distance_plot.r"],
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def fmt(seconds: float) -> str:
    """Format seconds as e.g. '2m 30s' or '54s'."""
    if seconds < 60:
        return f"{seconds:.0f}s"
    m, s = int(seconds // 60), int(seconds % 60)
    return f"{m}m {s:02d}s" if s else f"{m}m"


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def docker_pull(image: str) -> None:
    """Pull image if not already present."""
    r = subprocess.run(["docker", "image", "inspect", image], capture_output=True)
    if r.returncode != 0:
        print(f"  Pulling {image}")
        subprocess.run(
            ["docker", "pull", "--platform", DOCKER_PLATFORM, image], check=True
        )


def run_profiled(cmd: list[str], profile_dir: Path, label: str) -> dict:
    """
    Run a command with psrecord profiling. Returns a dict with
    wall_seconds, peak_rss_mb, and paths to the profiling outputs.
    """
    ensure_dir(profile_dir)
    log_path = profile_dir / "psrecord.log"
    plot_path = profile_dir / "profile.png"

    # Write command to a temp script so psrecord can exec it
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".sh", delete=False, prefix="bench_"
    ) as f:
        f.write("#!/bin/bash\nset -euo pipefail\n")
        # Quote args that need it
        parts = []
        for a in cmd:
            if any(c in a for c in " '\"\t\n;$(){}"):
                parts.append(f"'{a}'")
            else:
                parts.append(a)
        f.write(" ".join(parts) + "\n")
        script = f.name
    os.chmod(script, 0o755)

    short_cmd = " ".join(cmd[:5]) + ("..." if len(cmd) > 5 else "")
    print(f"  Running: {short_cmd}")

    start = time.time()
    try:
        subprocess.run(
            [
                "psrecord",
                f"bash {script}",
                "--log",
                str(log_path),
                "--plot",
                str(plot_path),
                "--include-children",
                "--interval",
                "1",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start
        print(f"  FAILED ({fmt(elapsed)})")
        if e.stderr:
            print(f"  stderr: {e.stderr[:500]}")
        raise
    finally:
        os.unlink(script)

    elapsed = time.time() - start

    # Parse peak RSS from psrecord log
    peak_rss = 0.0
    if log_path.exists():
        for line in log_path.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 3:
                try:
                    peak_rss = max(peak_rss, float(parts[2]))
                except ValueError:
                    pass

    result = dict(
        tool=label,
        wall_seconds=round(elapsed, 2),
        wall_formatted=fmt(elapsed),
        peak_rss_mb=round(peak_rss, 1) if peak_rss > 0 else None,
    )
    rss_str = f"{peak_rss:.0f} MB" if peak_rss > 0 else "N/A"
    print(f"  Done: {fmt(elapsed)}  (peak RSS: {rss_str})")
    return result


# ---------------------------------------------------------------------------
# Benchmark runners
# ---------------------------------------------------------------------------


def run_dupradar(ds_name: str, ds: dict, root: Path) -> dict:
    """Run R dupRadar via Docker."""
    print(f"\n{'=' * 60}")
    print(f"dupRadar (R) — {ds_name}")
    print(f"{'=' * 60}")
    docker_pull(DUPRADAR_IMG)

    outdir = ensure_dir(root / "benchmark" / "dupRadar" / ds_name)
    bam = (root / ds["bam"]).resolve()
    gtf = (root / ds["gtf"]).resolve()

    r_script = f"""\
library(dupRadar)
dm <- analyzeDuprates("/data/input.bam", "/data/input.gtf",
                      stranded={ds["stranded"]}, paired={str(ds["paired"]).upper()}, threads=1)
write.table(dm, file="/data/output/dupMatrix.txt", sep="\\t", row.names=FALSE, quote=FALSE)
fit <- duprateExpFit(DupMatrix=dm)
write.table(data.frame(intercept=fit$intercept, slope=fit$slope),
            file="/data/output/intercept_slope.txt", sep="\\t", row.names=FALSE, quote=FALSE)
png("/data/output/duprateExpDens.png", width=800, height=800); duprateExpDensPlot(DupMatrix=dm); dev.off()
pdf("/data/output/duprateExpDens.pdf"); duprateExpDensPlot(DupMatrix=dm); dev.off()
png("/data/output/duprateExpBoxplot.png", width=800, height=800); duprateExpBoxplot(DupMatrix=dm); dev.off()
pdf("/data/output/duprateExpBoxplot.pdf"); duprateExpBoxplot(DupMatrix=dm); dev.off()
png("/data/output/expressionHist.png", width=800, height=800); expressionHist(DupMatrix=dm); dev.off()
pdf("/data/output/expressionHist.pdf"); expressionHist(DupMatrix=dm); dev.off()
"""

    cmd = [
        "docker",
        "run",
        "--rm",
        "--platform",
        DOCKER_PLATFORM,
        "-v",
        f"{bam}:/data/input.bam:ro",
        "-v",
        f"{bam}.bai:/data/input.bam.bai:ro",
        "-v",
        f"{gtf}:/data/input.gtf:ro",
        "-v",
        f"{outdir}:/data/output",
        DUPRADAR_IMG,
        "Rscript",
        "-e",
        r_script,
    ]
    r = run_profiled(
        cmd, root / "benchmark" / "profiling" / "dupRadar" / ds_name, "dupRadar"
    )
    r["dataset"] = ds_name
    return r


def run_featurecounts(ds_name: str, ds: dict, root: Path) -> dict:
    """Run Subread featureCounts via Docker."""
    print(f"\n{'=' * 60}")
    print(f"featureCounts (Subread) — {ds_name}")
    print(f"{'=' * 60}")
    docker_pull(SUBREAD_IMG)

    outdir = ensure_dir(root / "benchmark" / "featureCounts" / ds_name)
    bam = (root / ds["bam"]).resolve()
    gtf = (root / ds["gtf"]).resolve()

    fc_args = [
        "featureCounts",
        "-a",
        "/data/input.gtf",
        "-o",
        "/data/output/featureCounts.tsv",
        "-T",
        "1",
        "-s",
        str(ds["stranded"]),
    ]
    if ds["paired"]:
        fc_args.append("-p")
    fc_args.append("/data/input.bam")

    cmd = [
        "docker",
        "run",
        "--rm",
        "--platform",
        DOCKER_PLATFORM,
        "-v",
        f"{bam}:/data/input.bam:ro",
        "-v",
        f"{bam}.bai:/data/input.bam.bai:ro",
        "-v",
        f"{gtf}:/data/input.gtf:ro",
        "-v",
        f"{outdir}:/data/output",
        SUBREAD_IMG,
    ] + fc_args
    r = run_profiled(
        cmd,
        root / "benchmark" / "profiling" / "featureCounts" / ds_name,
        "featureCounts",
    )
    r["dataset"] = ds_name
    return r


def run_rseqc_plots(tool: str, outdir: Path) -> None:
    """Run the R plot scripts produced by an RSeQC tool."""
    scripts = RSEQC_PLOT_SCRIPTS.get(tool, [])
    if not scripts:
        return

    docker_pull(RBASE_IMG)
    for script_name in scripts:
        script = outdir / script_name
        if not script.exists():
            print(f"  Warning: expected R script not found: {script}")
            continue
        print(f"  Running R plot: {script_name}")
        cmd = [
            "docker",
            "run",
            "--rm",
            "--platform",
            DOCKER_PLATFORM,
            "-v",
            f"{outdir}:/data/output",
            RBASE_IMG,
            "Rscript",
            f"/data/output/{script_name}",
        ]
        subprocess.run(cmd, check=True, capture_output=True)


def run_rseqc(ds_name: str, ds: dict, root: Path) -> list[dict]:
    """Run all 7 RSeQC tools via Docker, one at a time.

    Each tool's output is written to its own subdirectory under
    ``benchmark/rseqc/{ds_name}/{tool}/``.
    """
    print(f"\n{'=' * 60}")
    print(f"RSeQC (Python) — {ds_name}")
    print(f"{'=' * 60}")
    docker_pull(RSEQC_IMG)

    rseqc_root = root / "benchmark" / "rseqc" / ds_name
    bam = (root / ds["bam"]).resolve()
    bai = (root / (ds["bam"] + ".bai")).resolve()
    bed = (root / ds["bed"]).resolve()

    # Tool commands. Some write to stdout, some to stderr, some to files.
    # Wrap in bash to capture everything.  Output prefix is just the tool
    # name since each tool gets its own mounted directory.
    tool_cmds = {
        "bam_stat": "bam_stat.py -i /data/input.bam > /data/output/bam_stat.txt 2>&1",
        "infer_experiment": "infer_experiment.py -i /data/input.bam -r /data/input.bed > /data/output/infer_experiment.txt 2>&1",
        "read_duplication": "read_duplication.py -i /data/input.bam -o /data/output/read_duplication > /data/output/read_duplication.txt 2>&1",
        "read_distribution": "read_distribution.py -i /data/input.bam -r /data/input.bed > /data/output/read_distribution.txt 2>&1",
        "junction_annotation": "junction_annotation.py -i /data/input.bam -r /data/input.bed -o /data/output/junction_annotation 2> /data/output/junction_annotation.txt",
        "junction_saturation": "junction_saturation.py -i /data/input.bam -r /data/input.bed -o /data/output/junction_saturation 2> /data/output/junction_saturation.txt",
        "inner_distance": "inner_distance.py -i /data/input.bam -r /data/input.bed -o /data/output/inner_distance 2> /data/output/inner_distance.txt",
    }

    results = []
    for tool in RSEQC_TOOLS:
        print(f"\n--- {tool} ---")
        tool_outdir = ensure_dir(rseqc_root / tool)
        docker_cmd = [
            "docker",
            "run",
            "--rm",
            "--platform",
            DOCKER_PLATFORM,
            "-v",
            f"{bam}:/data/input.bam:ro",
            "-v",
            f"{bai}:/data/input.bam.bai:ro",
            "-v",
            f"{bed}:/data/input.bed:ro",
            "-v",
            f"{tool_outdir}:/data/output",
            RSEQC_IMG,
            "bash",
            "-c",
            tool_cmds[tool],
        ]
        r = run_profiled(
            docker_cmd,
            root / "benchmark" / "profiling" / "rseqc" / tool / ds_name,
            tool,
        )
        r["dataset"] = ds_name
        results.append(r)

        # Run any R plot scripts this tool generated
        if tool in RSEQC_PLOT_SCRIPTS:
            run_rseqc_plots(tool, tool_outdir)

    return results


def run_preseq(ds_name: str, ds: dict, root: Path) -> dict:
    """Run preseq lc_extrap via Docker."""
    print(f"\n{'=' * 60}")
    print(f"preseq lc_extrap — {ds_name}")
    print(f"{'=' * 60}")
    docker_pull(PRESEQ_IMG)

    outdir = ensure_dir(root / "benchmark" / "preseq" / ds_name)
    bam = (root / ds["bam"]).resolve()

    preseq_args = ["preseq", "lc_extrap", "-bam", "-seed", "1"]
    if ds["paired"]:
        preseq_args.append("-pe")
    preseq_args.extend(["-o", "/data/output/lc_extrap.txt", "/data/input.bam"])

    cmd = [
        "docker",
        "run",
        "--rm",
        "--platform",
        DOCKER_PLATFORM,
        "-v",
        f"{bam}:/data/input.bam:ro",
        "-v",
        f"{bam}.bai:/data/input.bam.bai:ro",
        "-v",
        f"{outdir}:/data/output",
        PRESEQ_IMG,
    ] + preseq_args
    r = run_profiled(
        cmd, root / "benchmark" / "profiling" / "preseq" / ds_name, "preseq"
    )
    r["dataset"] = ds_name
    return r


def run_samtools(ds_name: str, ds: dict, root: Path) -> list[dict]:
    """Run samtools flagstat, idxstats, and stats via Docker."""
    print(f"\n{'=' * 60}")
    print(f"samtools — {ds_name}")
    print(f"{'=' * 60}")
    docker_pull(SAMTOOLS_IMG)

    outdir = ensure_dir(root / "benchmark" / "samtools" / ds_name)
    bam = (root / ds["bam"]).resolve()

    docker_base = [
        "docker",
        "run",
        "--rm",
        "--platform",
        DOCKER_PLATFORM,
        "-v",
        f"{bam}:/data/input.bam:ro",
        "-v",
        f"{bam}.bai:/data/input.bam.bai:ro",
        "-v",
        f"{outdir}:/data/output",
        SAMTOOLS_IMG,
    ]

    tool_cmds = {
        "samtools_flagstat": "samtools flagstat /data/input.bam > /data/output/flagstat.txt",
        "samtools_idxstats": "samtools idxstats /data/input.bam > /data/output/idxstats.txt",
        "samtools_stats": "samtools stats /data/input.bam > /data/output/stats.txt",
    }

    results = []
    for tool, tool_cmd in tool_cmds.items():
        print(f"\n--- {tool} ---")
        cmd = docker_base + ["bash", "-c", tool_cmd]
        r = run_profiled(
            cmd,
            root / "benchmark" / "profiling" / "samtools" / tool / ds_name,
            tool,
        )
        r["dataset"] = ds_name
        results.append(r)

    return results


def run_tin(ds_name: str, ds: dict, root: Path) -> dict:
    """Run RSeQC tin.py via Docker."""
    print(f"\n{'=' * 60}")
    print(f"RSeQC: tin.py — {ds_name}")
    print(f"{'=' * 60}")
    docker_pull(RSEQC_IMG)

    outdir = ensure_dir(root / "benchmark" / "rseqc" / ds_name / "tin")
    bam = (root / ds["bam"]).resolve()
    bai = (root / (ds["bam"] + ".bai")).resolve()
    bed = (root / ds["bed"]).resolve()

    cmd = [
        "docker",
        "run",
        "--rm",
        "--platform",
        DOCKER_PLATFORM,
        "-v",
        f"{bam}:/data/input.bam:ro",
        "-v",
        f"{bai}:/data/input.bam.bai:ro",
        "-v",
        f"{bed}:/data/input.bed:ro",
        "-v",
        f"{outdir}:/data/output",
        "-w",
        "/data/output",
        RSEQC_IMG,
        "tin.py",
        "-i",
        "/data/input.bam",
        "-r",
        "/data/input.bed",
    ]
    r = run_profiled(cmd, root / "benchmark" / "profiling" / "tin" / ds_name, "tin")
    r["dataset"] = ds_name
    return r


def run_qualimap(ds_name: str, ds: dict, root: Path) -> dict:
    """Run Qualimap rnaseq via Docker for gene body coverage."""
    print(f"\n{'=' * 60}")
    print(f"Qualimap rnaseq — {ds_name}")
    print(f"{'=' * 60}")
    docker_pull(QUALIMAP_IMG)

    outdir = ensure_dir(root / "benchmark" / "qualimap" / ds_name)
    bam = (root / ds["bam"]).resolve()
    gtf = (root / ds["gtf"]).resolve()

    cmd = [
        "docker",
        "run",
        "--rm",
        "--platform",
        DOCKER_PLATFORM,
        "-v",
        f"{bam}:/data/input.bam:ro",
        "-v",
        f"{bam}.bai:/data/input.bam.bai:ro",
        "-v",
        f"{gtf}:/data/input.gtf:ro",
        "-v",
        f"{outdir}:/data/output",
        QUALIMAP_IMG,
        "qualimap",
        "rnaseq",
        "-bam",
        "/data/input.bam",
        "-gtf",
        "/data/input.gtf",
        "-outdir",
        "/data/output",
        "--java-mem-size=8G",
    ]
    r = run_profiled(
        cmd, root / "benchmark" / "profiling" / "qualimap" / ds_name, "qualimap"
    )
    r["dataset"] = ds_name
    return r


def run_rustqc(ds_name: str, ds: dict, root: Path) -> dict:
    """Run RustQC natively."""
    print(f"\n{'=' * 60}")
    print(f"RustQC — {ds_name}")
    print(f"{'=' * 60}")

    outdir = ensure_dir(root / "benchmark" / "RustQC" / ds_name)
    bam = root / ds["bam"]
    gtf = root / ds["gtf"]

    cmd = [
        str(root / "target" / "release" / "rustqc"),
        "rna",
        str(bam),
        "--gtf",
        str(gtf),
        "-p",
        "-t",
        "10",
        "-o",
        str(outdir),
    ]
    if ds["skip_dup_check"]:
        cmd.append("--skip-dup-check")
    if ds["config"]:
        cmd.extend(["-c", str(root / ds["config"])])
    if ds["biotype_attribute"]:
        cmd.extend(["--biotype-attribute", ds["biotype_attribute"]])

    r = run_profiled(
        cmd, root / "benchmark" / "profiling" / "RustQC" / ds_name, "RustQC"
    )
    r["dataset"] = ds_name
    return r


# ---------------------------------------------------------------------------
# SVG bar chart generation
# ---------------------------------------------------------------------------


def generate_svgs(results: list[dict], root: Path) -> None:
    """Generate dark/light SVG bar charts from large-dataset results."""
    large = [r for r in results if r.get("dataset") == "large"]
    if not large:
        print("No large-dataset results — skipping SVG generation.")
        return

    # Order: RustQC first, then upstream tools sorted by category
    order = (
        ["RustQC", "featureCounts", "dupRadar"]
        + RSEQC_TOOLS
        + [
            "tin",
            "preseq",
            "samtools_flagstat",
            "samtools_idxstats",
            "samtools_stats",
            "qualimap",
        ]
    )
    # Pretty display names for SVG labels
    display_names = {
        "RustQC": "RustQC",
        "featureCounts": "featureCounts",
        "dupRadar": "dupRadar",
        "bam_stat": "RSeQC: bam_stat",
        "infer_experiment": "RSeQC: infer_experiment",
        "read_duplication": "RSeQC: read_duplication",
        "read_distribution": "RSeQC: read_distribution",
        "junction_annotation": "RSeQC: junction_annotation",
        "junction_saturation": "RSeQC: junction_saturation",
        "inner_distance": "RSeQC: inner_distance",
        "tin": "RSeQC: tin.py",
        "preseq": "preseq lc_extrap",
        "samtools_flagstat": "samtools flagstat",
        "samtools_idxstats": "samtools idxstats",
        "samtools_stats": "samtools stats",
        "qualimap": "Qualimap rnaseq",
    }
    bars: list[tuple[str, float]] = []
    for name in order:
        match = [r for r in large if r["tool"] == name]
        if match:
            bars.append((display_names.get(name, name), match[0]["wall_seconds"]))
    if not bars:
        return

    for theme in ("dark", "light"):
        _write_svg(bars, root, theme)


def _write_svg(bars: list[tuple[str, float]], root: Path, theme: str) -> None:
    # Catppuccin Mocha (dark) / Latte (light)
    if theme == "dark":
        grid, label_c, text_c, accent, bar_c = (
            "#313244",
            "#a6adc8",
            "#cdd6f4",
            "#cba6f7",
            "#7f849c",
        )
    else:
        grid, label_c, text_c, accent, bar_c = (
            "#ccd0da",
            "#6c6f85",
            "#4c4f69",
            "#8839ef",
            "#7c7f93",
        )

    label_w = 180
    chart_l = label_w + 10
    chart_r = 720
    cw = chart_r - chart_l
    bh, gap, top = 24, 8, 10
    bot = 28
    n = len(bars)
    h = top + n * (bh + gap) - gap + bot
    chart_bot = top + n * (bh + gap) - gap

    max_s = max(s for _, s in bars)
    # Choose tick interval
    for limit, iv in [(60, 15), (300, 60), (600, 120), (1800, 300), (3600, 600)]:
        if max_s <= limit:
            tick_iv = iv
            break
    else:
        tick_iv = 1800
    max_t = ((int(max_s * 1.15) // tick_iv) + 1) * tick_iv

    def x(s):
        return chart_l + (s / max_t) * cw

    def tick_label(s):
        if s == 0:
            return "0s"
        return f"{s // 60}m" if s >= 60 else f"{s}s"

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 720 {h}" '
        f'width="720" height="{h}" '
        f"style=\"font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;\">"
    ]

    # Gridlines + tick labels
    t = 0
    while t <= max_t:
        xv = x(t)
        lines.append(
            f'  <line x1="{xv:.1f}" y1="{top}" x2="{xv:.1f}" y2="{chart_bot}" stroke="{grid}" stroke-width="1"/>'
        )
        lines.append(
            f'  <text x="{xv:.1f}" y="{chart_bot + 16}" text-anchor="middle" fill="{label_c}" font-size="12">{tick_label(t)}</text>'
        )
        t += tick_iv

    # Bars
    for i, (name, secs) in enumerate(bars):
        y = top + i * (bh + gap)
        w = max((secs / max_t) * cw, 2)
        is_rqc = name == "RustQC"
        fill = accent if is_rqc else bar_c
        wt = "bold" if is_rqc else "normal"
        nc = text_c if is_rqc else label_c
        ty = y + bh / 2 + 5

        lines.append(
            f'  <text x="{chart_l - 12}" y="{ty:.1f}" text-anchor="end" fill="{nc}" font-size="14" font-weight="{wt}">{name}</text>'
        )
        lines.append(
            f'  <rect x="{chart_l}" y="{y:.1f}" width="{w:.1f}" height="{bh}" fill="{fill}"/>'
        )
        lines.append(
            f'  <text x="{chart_l + w + 8:.1f}" y="{ty:.1f}" text-anchor="start" fill="{text_c}" font-size="14" font-weight="{wt}">{fmt(secs)}</text>'
        )

    lines.append("</svg>")

    out = ensure_dir(root / "docs" / "public" / "benchmarks") / f"benchmark_{theme}.svg"
    out.write_text("\n".join(lines) + "\n")
    print(f"  Wrote {out}")


# ---------------------------------------------------------------------------
# Results summary
# ---------------------------------------------------------------------------


def print_and_save_results(results: list[dict], root: Path) -> None:
    summary_dir = ensure_dir(root / "benchmark" / "profiling")

    # JSON
    json_path = summary_dir / "results.json"
    json_path.write_text(json.dumps(results, indent=2) + "\n")

    # Print table
    print(f"\n{'=' * 70}")
    print("RESULTS")
    print(f"{'=' * 70}")
    print(f"{'Tool':<25} {'Dataset':<8} {'Wall Time':>12} {'Peak RSS':>12}")
    print(f"{'-' * 25} {'-' * 8} {'-' * 12} {'-' * 12}")
    for r in results:
        rss = f"{r['peak_rss_mb']:.0f} MB" if r.get("peak_rss_mb") else "N/A"
        print(f"{r['tool']:<25} {r['dataset']:<8} {r['wall_formatted']:>12} {rss:>12}")
    print(f"{'=' * 70}")
    print(f"\nSaved to {json_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description="Run RustQC benchmarks.")
    g = parser.add_argument_group("Tools (opt-in)")
    g.add_argument("--all", action="store_true", help="Run all benchmarks")
    g.add_argument("--rustqc", action="store_true")
    g.add_argument("--dupradar", action="store_true")
    g.add_argument("--featurecounts", action="store_true")
    g.add_argument("--rseqc", action="store_true")
    g.add_argument("--preseq", action="store_true")
    g.add_argument("--samtools", action="store_true")
    g.add_argument("--tin", action="store_true")
    g.add_argument("--qualimap", action="store_true")
    d = parser.add_argument_group("Dataset")
    d.add_argument("--small-only", action="store_true")
    d.add_argument("--large-only", action="store_true")
    parser.add_argument("--skip-svg", action="store_true", help="Skip SVG generation")
    args = parser.parse_args()

    tool_selected = (
        args.all
        or args.rustqc
        or args.dupradar
        or args.featurecounts
        or args.rseqc
        or args.preseq
        or args.samtools
        or args.tin
        or args.qualimap
    )
    if not tool_selected:
        parser.error(
            "Select tools: --all, --rustqc, --dupradar, --featurecounts, --rseqc, "
            "--preseq, --samtools, --tin, --qualimap"
        )

    # Must run from repo root
    root = Path.cwd()
    assert (root / "Cargo.toml").exists(), (
        "Run from the repository root (where Cargo.toml is)"
    )

    # Check prerequisites
    if args.rustqc or args.all:
        assert (root / "target/release/rustqc").exists(), (
            "Build first: cargo build --release"
        )
    needs_docker = (
        args.dupradar
        or args.featurecounts
        or args.rseqc
        or args.preseq
        or args.samtools
        or args.tin
        or args.qualimap
        or args.all
    )
    if needs_docker:
        r = subprocess.run(["docker", "info"], capture_output=True)
        assert r.returncode == 0, "Docker is not running"

    # Which datasets?
    ds_names = ["small", "large"]
    if args.small_only:
        ds_names = ["small"]
    elif args.large_only:
        ds_names = ["large"]

    # Verify inputs exist
    for name in ds_names:
        ds = DATASETS[name]
        for key in ("bam", "gtf", "bed"):
            p = root / ds[key]
            assert p.exists(), f"Missing input: {p}\nSee benchmark/README.md for setup."

    results: list[dict] = []

    for name in ds_names:
        ds = DATASETS[name]
        if args.dupradar or args.all:
            results.append(run_dupradar(name, ds, root))
        if args.featurecounts or args.all:
            results.append(run_featurecounts(name, ds, root))
        if args.rseqc or args.all:
            results.extend(run_rseqc(name, ds, root))
        if args.preseq or args.all:
            results.append(run_preseq(name, ds, root))
        if args.samtools or args.all:
            results.extend(run_samtools(name, ds, root))
        if args.tin or args.all:
            results.append(run_tin(name, ds, root))
        if args.qualimap or args.all:
            results.append(run_qualimap(name, ds, root))
        if args.rustqc or args.all:
            results.append(run_rustqc(name, ds, root))

    if results:
        print_and_save_results(results, root)
        if not args.skip_svg:
            print("\nGenerating SVG bar charts...")
            generate_svgs(results, root)

    print("\nDone!")


if __name__ == "__main__":
    main()
