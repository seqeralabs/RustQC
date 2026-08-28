#!/usr/bin/env bash
# Regenerate the DNA test inputs and the reference outputs they are compared against.
#
# Inputs come from nf-core/test-datasets (a real human chr22 slice, 40 kb).
# The upstream BAM is not duplicate-marked, so this script marks duplicates
# with samtools; RustQC requires duplicate-marked input.
#
# The reference outputs are produced by the pinned tool versions recorded in
# tests/expected/dna/VERSIONS.txt. Regenerating with a different version will
# make the parity tests fail, which is the intended behaviour: fixtures and
# tool versions travel together.
set -euo pipefail

MOSDEPTH_VERSION="0.3.14"
SAMTOOLS_VERSION="1.24"

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
data="$here/data/dna"
expected="$here/expected/dna"
base="https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens"

have() { command -v "$1" >/dev/null || { echo "missing tool: $1" >&2; exit 1; }; }
have samtools; have mosdepth; have curl

check_version() {
  local tool="$1" want="$2" got
  got="$($tool --version 2>&1 | head -1 | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1)"
  if [[ "$got" != "$want" ]]; then
    echo "$tool version $got does not match the pinned $want" >&2
    echo "Install the pinned version, or update VERSIONS.txt and the fixtures together." >&2
    exit 1
  fi
}
check_version samtools "$SAMTOOLS_VERSION"
check_version mosdepth "$MOSDEPTH_VERSION"

mkdir -p "$data" "$expected"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

curl -sSfL -o "$tmp/upstream.bam"      "$base/illumina/bam/test.paired_end.sorted.bam"
curl -sSfL -o "$data/genome.fasta"     "$base/genome/genome.fasta"
curl -sSfL -o "$data/genome.fasta.fai" "$base/genome/genome.fasta.fai"

# Mark duplicates: name-sort, add mate tags, coordinate-sort, then markdup.
samtools sort -n -o "$tmp/ns.bam" "$tmp/upstream.bam"
samtools fixmate -m "$tmp/ns.bam" "$tmp/fm.bam"
samtools sort -o "$tmp/cs.bam" "$tmp/fm.bam"
samtools markdup -S "$tmp/cs.bam" "$data/test.dna.bam"
samtools index "$data/test.dna.bam"

mosdepth --by 500 --thresholds 1,5,10,15,20,30,50 "$expected/test" "$data/test.dna.bam"
samtools stats    "$data/test.dna.bam" > "$expected/test.stats.txt"
samtools flagstat "$data/test.dna.bam" > "$expected/test.flagstat.txt"
samtools idxstats "$data/test.dna.bam" > "$expected/test.idxstats.txt"

printf 'mosdepth\t%s\nsamtools\t%s\n' "$MOSDEPTH_VERSION" "$SAMTOOLS_VERSION" > "$expected/VERSIONS.txt"

echo "Regenerated $(find "$data" "$expected" -type f | wc -l | tr -d ' ') files."
