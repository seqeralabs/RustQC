#!/usr/bin/env bash
set -euo pipefail

# CPU scaling benchmark for RustQC on local hardware (Apple M1 Pro, 10 cores)
# Matches the per-tool benchmark setup: -p -s 2 --biotype-attribute gene_type

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

RUSTQC="$REPO_DIR/target/release/rustqc"
BAM="$REPO_DIR/benchmark/input/large/GM12878_REP1.markdup.sorted.bam"
GTF="$REPO_DIR/benchmark/input/large/genes.gtf"
CFG="$REPO_DIR/benchmark/input/large/config.yaml"
OUTDIR="$REPO_DIR/benchmark/cpu_scaling_results"
RESULTS="$OUTDIR/results.csv"

MAX_THREADS=10

# Preflight checks
for f in "$RUSTQC" "$BAM" "$BAM.bai" "$GTF" "$CFG"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Missing file: $f" >&2
        exit 1
    fi
done

mkdir -p "$OUTDIR"
echo "threads,wall_seconds" > "$RESULTS"

for t in $(seq 1 "$MAX_THREADS"); do
    DIR="$OUTDIR/t${t}"

    # Skip if already completed
    if [ -d "$DIR" ] && [ "$(ls -A "$DIR" 2>/dev/null)" ]; then
        echo "Skipping t=$t (output already exists in $DIR)"
        # Try to grab the time from the log if it exists
        if [ -f "$DIR/time.txt" ]; then
            ELAPSED=$(cat "$DIR/time.txt")
            echo "$t,$ELAPSED" >> "$RESULTS"
        fi
        continue
    fi

    rm -rf "$DIR"
    mkdir -p "$DIR"

    echo "Running t=$t at $(date +%H:%M:%S)..."
    START=$(python3 -c 'import time; print(time.time())')

    "$RUSTQC" rna "$BAM" \
        --gtf "$GTF" \
        -p -s 2 \
        -t "$t" \
        -o "$DIR" \
        -c "$CFG" \
        --biotype-attribute gene_type \
        2>/dev/null

    END=$(python3 -c 'import time; print(time.time())')
    ELAPSED=$(python3 -c "print(f'{$END - $START:.1f}')")

    echo "$ELAPSED" > "$DIR/time.txt"
    echo "$t,$ELAPSED" >> "$RESULTS"
    echo "  t=$t done: ${ELAPSED}s"
done

echo ""
echo "=== RESULTS ==="
cat "$RESULTS"
