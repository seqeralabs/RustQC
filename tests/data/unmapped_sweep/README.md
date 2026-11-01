# Unmapped-sweep CRAM fixtures

`sweep.sam` has two mapped reads on chr1, a mapped read on chr2 with its unmapped mate placed at the same position (flag 0x85, tid 1), and two unplaced reads (no reference). The CRAMs are that file written by samtools 1.22.1 against `sweep.fa`:

```bash
samtools view -C -T sweep.fa -o single_ref_slices.cram sweep.sam
samtools view -C -T sweep.fa --output-fmt-option multi_seq_per_slice=1 -o multi_ref_slice.cram sweep.sam
samtools index single_ref_slices.cram
samtools index multi_ref_slice.cram
```

`single_ref_slices.cram` has one slice per reference plus one for the unplaced reads; `multi_ref_slice.cram` holds all six records in one multi-reference slice, so its index points the unplaced-read entry at a slice that also holds placed records.

`placed_only.sam` drops the two unplaced reads from `sweep.sam`. It pins the empty-unplaced-sweep case; regenerate with samtools 1.22.1:

```bash
samtools view --no-PG -C -T sweep.fa -o placed_only.cram placed_only.sam
samtools index placed_only.cram
```
