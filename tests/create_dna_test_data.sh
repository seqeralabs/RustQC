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
PICARD_VERSION="3.4.0"
SAMTOOLS_VERSION="1.24"

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
data="$here/data/dna"
expected="$here/expected/dna"
base="https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens"

have() { command -v "$1" >/dev/null || { echo "missing tool: $1" >&2; exit 1; }; }
have samtools; have mosdepth; have curl; have java

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
curl -sSfL -o "$data/targets.bed"      "$base/genome/genome.multi_intervals.bed"

# Mark duplicates: name-sort, add mate tags, coordinate-sort, then markdup.
samtools sort -n -o "$tmp/ns.bam" "$tmp/upstream.bam"
samtools fixmate -m "$tmp/ns.bam" "$tmp/fm.bam"
samtools sort -o "$tmp/cs.bam" "$tmp/fm.bam"
samtools markdup -S "$tmp/cs.bam" "$data/test.dna.bam"
samtools index "$data/test.dna.bam"

mosdepth --by 500 --thresholds 1,5,10,15,20,30,50 "$expected/test" "$data/test.dna.bam"

# Picard is a jar rather than a command, so it is fetched by version instead of
# version-checked. The JVM locale is pinned: a French default locale writes
# "3,531312" where an English one writes "3.531312", which would make the
# fixtures depend on the machine that produced them.
picard_jar="$tmp/picard-$PICARD_VERSION.jar"
curl -sSfL -o "$picard_jar" \
  "https://github.com/broadinstitute/picard/releases/download/$PICARD_VERSION/picard.jar"
picard() {
  java -Duser.language=en -Duser.country=US -jar "$picard_jar" "$@" 2>/dev/null
}

picard CollectWgsMetrics \
  -I "$data/test.dna.bam" \
  -O "$expected/test.wgs_metrics.txt" \
  -R "$data/genome.fasta"

picard CollectInsertSizeMetrics \
  -I "$data/test.dna.bam" \
  -O "$expected/test.insert_size_metrics.txt" \
  -H "$tmp/insert_size_histogram.pdf"

# The chart output needs R, so it goes to the scratch directory and is not
# compared against; only the two metrics tables are fixtures.
# Picard consumes interval lists rather than BED, so the targets are converted
# with Picard's own tool. The two conventions differ: BED is zero-based
# half-open, an interval list one-based inclusive.
picard CreateSequenceDictionary -R "$data/genome.fasta" -O "$tmp/genome.dict"
picard BedToIntervalList \
  -I "$data/targets.bed" \
  -O "$tmp/targets.interval_list" \
  -SD "$tmp/genome.dict"

picard CollectHsMetrics \
  -I "$data/test.dna.bam" \
  -O "$expected/test.hs_metrics.txt" \
  -R "$data/genome.fasta" \
  -BI "$tmp/targets.interval_list" \
  -TI "$tmp/targets.interval_list"

picard CollectGcBiasMetrics \
  -I "$data/test.dna.bam" \
  -O "$expected/test.gc_bias.detail_metrics.txt" \
  -S "$expected/test.gc_bias.summary_metrics.txt" \
  -CHART "$tmp/gc_bias.pdf" \
  -R "$data/genome.fasta"

# Picard stamps a start time and the full command line, absolute paths and all,
# into the first four lines of every metrics file. Those are dropped: they would
# change on every regeneration and say nothing about the numbers. The
# "## METRICS CLASS" and "## HISTOGRAM" markers further down are part of the
# format and are kept.
for f in "$expected/test.wgs_metrics.txt" "$expected/test.insert_size_metrics.txt" \
         "$expected/test.gc_bias.detail_metrics.txt" "$expected/test.gc_bias.summary_metrics.txt" \
         "$expected/test.hs_metrics.txt"; do
  sed -e '/^## htsjdk\.samtools\.metrics\.StringHeader$/d' -e '/^# /d' "$f" \
    | sed -e '/./,$!d' > "$f.tmp" && mv "$f.tmp" "$f"
done
samtools stats    "$data/test.dna.bam" > "$expected/test.stats.txt"
samtools flagstat "$data/test.dna.bam" > "$expected/test.flagstat.txt"
samtools idxstats "$data/test.dna.bam" > "$expected/test.idxstats.txt"

printf 'mosdepth\t%s\nsamtools\t%s\npicard\t%s\n' \
  "$MOSDEPTH_VERSION" "$SAMTOOLS_VERSION" "$PICARD_VERSION" > "$expected/VERSIONS.txt"

echo "Regenerated $(find "$data" "$expected" -type f | wc -l | tr -d ' ') files."
