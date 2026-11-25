#!/usr/bin/env bash
# Regenerate the protein test inputs and the reference outputs they are
# compared against.
#
# Inputs come from nf-core/test-datasets: two small real protein FASTA files.
# The reference statistics come from seqkit, pinned to the version recorded in
# tests/expected/protein/VERSIONS.txt. Regenerating with a different version
# will make the parity tests fail, which is the intended behaviour.
set -euo pipefail

SEQKIT_VERSION="2.13.0"

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
data="$here/data/protein"
expected="$here/expected/protein"
base="https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/proteomics/database"

have() { command -v "$1" >/dev/null || { echo "missing tool: $1" >&2; exit 1; }; }
have seqkit; have curl

got="$(seqkit version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"
if [[ "$got" != "$SEQKIT_VERSION" ]]; then
  echo "seqkit version $got does not match the pinned $SEQKIT_VERSION" >&2
  echo "Install the pinned version, or update VERSIONS.txt and the fixtures together." >&2
  exit 1
fi

mkdir -p "$data" "$expected"
curl -sSfL -o "$data/yeast_UPS_mini.fasta"          "$base/yeast_UPS_mini.fasta"
curl -sSfL -o "$data/protein_mini_with_cazymes.faa" "$base/protein_mini_with_cazymes.faa"

# -T gives tab-separated output, which is what RustQC reproduces. The file
# column holds the path it was run with, so it is rewritten to the bare
# basename to keep the fixture independent of where it was generated.
for f in yeast_UPS_mini.fasta protein_mini_with_cazymes.faa; do
  seqkit stats -a -T "$data/$f" \
    | awk -v n="$f" 'BEGIN{FS=OFS="\t"} NR==1{print; next} {$1=n; print}' \
    > "$expected/${f%.*}.seqkit.tsv"
done

printf 'seqkit\t%s\n' "$SEQKIT_VERSION" > "$expected/VERSIONS.txt"

echo "Regenerated $(find "$data" "$expected" -type f | wc -l | tr -d ' ') files."
