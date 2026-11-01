# Unmapped-sweep CRAM fixtures

`sweep.sam` has two mapped reads on chr1, a mapped read on chr2 with its unmapped mate placed at the same position (flag 0x85, tid 1), and two unplaced reads (no reference). The CRAMs are that file written by samtools 1.22.1 against `sweep.fa`:

```bash
samtools view -C -T sweep.fa -o single_ref_slices.cram sweep.sam
samtools view -C -T sweep.fa --output-fmt-option multi_seq_per_slice=1 -o multi_ref_slice.cram sweep.sam
samtools index single_ref_slices.cram
samtools index multi_ref_slice.cram
```

`single_ref_slices.cram` has one slice per reference plus one for the unplaced reads; `multi_ref_slice.cram` holds all six records in one multi-reference slice, so its index points the unplaced-read entry at a slice that also holds placed records.

`same_reference_multi_slice.cram` contains one chr1 container with two slices plus the unplaced-read container. Both chr1 slices hold r0 and r1, deliberately retaining repeated names as distinct physical records. Its CRAI has two chr1 entries at the same container offset. Regenerate with:

```bash
python3 make_multi_slice.py
samtools index same_reference_multi_slice.cram
samtools quickcheck -v same_reference_multi_slice.cram
samtools view -T sweep.fa same_reference_multi_slice.cram
```

The generator reuses the original fixture's encoded slices, updates the second slice's record counter and the container's counts, landmarks, size and CRC. It omits the chr2 container, keeping the existing unplaced container's record counter correct. Samtools 1.22.1 decodes six records: r0, r1, r0, r1, u1, u2. Record-level or read-name deduplication would incorrectly discard two legitimate records.

`placed_only.sam` drops the two unplaced reads from `sweep.sam`. It pins the empty-unplaced-sweep case; regenerate with samtools 1.22.1:

```bash
samtools view --no-PG -C -T sweep.fa -o placed_only.cram placed_only.sam
samtools index placed_only.cram
```
