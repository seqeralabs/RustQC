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
PYTEOMICS_VERSION="5.0.1"

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
data="$here/data/protein"
expected="$here/expected/protein"
base="https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/proteomics/database"

have() { command -v "$1" >/dev/null || { echo "missing tool: $1" >&2; exit 1; }; }
have seqkit; have curl; have python3

got="$(seqkit version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"
if [[ "$got" != "$SEQKIT_VERSION" ]]; then
  echo "seqkit version $got does not match the pinned $SEQKIT_VERSION" >&2
  echo "Install the pinned version, or update VERSIONS.txt and the fixtures together." >&2
  exit 1
fi

mkdir -p "$data" "$expected"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
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

# The mzML fixture comes from mzdata's own test data: 48 spectra, 14 MS1 and
# 34 MS2, which is enough to exercise every metric. nf-core's proteomics
# fixtures are either a single profile scan or far past the size budget.
curl -sSfL -o "$data/small.mzML" \
  "https://raw.githubusercontent.com/mobiusklein/mzdata/main/test/data/small.mzML"

# The reference figures come from pyteomics, installed into a throwaway
# environment so the host's Python is left alone. Total ion current is summed
# from the peak intensities rather than read from the header, matching what
# RustQC reports: the two differ when a profile spectrum has been centroided.
python3 -m venv "$tmp/venv"
"$tmp/venv/bin/pip" install --quiet "pyteomics==$PYTEOMICS_VERSION" numpy psims lxml
"$tmp/venv/bin/python" - "$data/small.mzML" "$expected/small.pyteomics.tsv" <<'PYEOF'
import collections, sys
from pyteomics import mzml

source, destination = sys.argv[1], sys.argv[2]
levels = collections.Counter(); peaks = collections.Counter()
tic = collections.Counter(); lo = {}; hi = {}
rts = []; charges = collections.Counter(); nocharge = 0; nprecursors = 0; mzs = []

with mzml.read(source) as reader:
    for spectrum in reader:
        level = int(spectrum["ms level"])
        intensity = spectrum.get("intensity array")
        n = len(intensity) if intensity is not None else 0
        levels[level] += 1
        peaks[level] += n
        tic[level] += float(intensity.sum()) if n else 0.0
        lo[level] = n if level not in lo else min(lo[level], n)
        hi[level] = n if level not in hi else max(hi[level], n)
        start = spectrum.get("scanList", {}).get("scan", [{}])[0].get("scan start time")
        if start is not None:
            rts.append(float(start))
        for precursor in spectrum.get("precursorList", {}).get("precursor", []):
            for ion in precursor.get("selectedIonList", {}).get("selectedIon", []):
                nprecursors += 1
                charge = ion.get("charge state")
                if charge is None:
                    nocharge += 1
                else:
                    charges[int(charge)] += 1
                mz = ion.get("selected ion m/z")
                if mz is not None:
                    mzs.append(float(mz))

rows = ["metric\tvalue"]
rows.append(f"spectra\t{sum(levels.values())}")
rows.append(f"peaks\t{sum(peaks.values())}")
rows.append(f"rt_min\t{min(rts):.6f}")
rows.append(f"rt_max\t{max(rts):.6f}")
rows.append(f"precursors\t{nprecursors}")
rows.append(f"precursors_without_charge\t{nocharge}")
if mzs:
    rows.append(f"precursor_mz_min\t{min(mzs):.4f}")
    rows.append(f"precursor_mz_max\t{max(mzs):.4f}")
for level in sorted(levels):
    rows.append(f"ms{level}_spectra\t{levels[level]}")
    rows.append(f"ms{level}_peaks\t{peaks[level]}")
    rows.append(f"ms{level}_min_peaks\t{lo[level]}")
    rows.append(f"ms{level}_max_peaks\t{hi[level]}")
    rows.append(f"ms{level}_total_ion_current\t{tic[level]:.4f}")
open(destination, "w").write("\n".join(rows) + "\n")
PYEOF

printf 'seqkit\t%s\npyteomics\t%s\n' "$SEQKIT_VERSION" "$PYTEOMICS_VERSION" > "$expected/VERSIONS.txt"

echo "Regenerated $(find "$data" "$expected" -type f | wc -l | tr -d ' ') files."
