#!/usr/bin/env bash
# Regenerate the split_bam parity fixture and its RSeQC reference outputs.
#
# Inputs (hand-designed, checked in):
#   data/split_bam_test.sam  - 33 records covering every branch of split_bam.py
#                              (each record's XC:Z tag holds its intended class)
#   data/split_bam_test.bed  - BED12 with multi-block, zero-length, overlapping,
#                              zip-truncated, mixed-case and absent-chromosome entries
#
# Outputs:
#   data/split_bam_test.bam(.bai)          - coordinate-sorted, indexed fixture
#   expected/split_bam/split_bam_test.*    - RSeQC split_bam.py stdout + 3 BAMs
#
# Requires samtools and RSeQC 5.0.4 (split_bam.py) on PATH.
set -euo pipefail
cd "$(dirname "$0")"

samtools sort --no-PG -o data/split_bam_test.bam data/split_bam_test.sam
samtools index data/split_bam_test.bam

mkdir -p expected/split_bam
cd expected/split_bam
split_bam.py \
    -i ../../data/split_bam_test.bam \
    -r ../../data/split_bam_test.bed \
    -o split_bam_test \
    > split_bam_test.split_bam.txt
