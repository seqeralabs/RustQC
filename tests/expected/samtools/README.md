# samtools reference outputs

`all_flags.flagstat` and `all_flags.stats.sn` are `samtools flagstat` and the `SN` lines of
`samtools stats` from samtools 1.22.1 (bioconda build `h96c455f_0`) for the records
`all_flags_sam_lines()` in `tests/integration_test.rs` writes: all 4096 flag values in three
mate placements. To regenerate:

```bash
{ printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n'
  awk 'BEGIN{for(f=0;f<4096;f++){printf "f%da\t%d\tchr1\t10\t60\t4M\t=\t10\t0\tACGT\tIIII\n",f,f; printf "f%db\t%d\tchr1\t10\t60\t4M\tchr2\t10\t0\tACGT\tIIII\n",f,f}; for(f=0;f<4096;f++) printf "f%dc\t%d\tchr2\t10\t3\t4M\tchr1\t10\t0\tACGT\tIIII\n",f,f}'
} > all_flags.sam
samtools flagstat all_flags.sam > tests/expected/samtools/all_flags.flagstat
samtools stats all_flags.sam | grep '^SN' > tests/expected/samtools/all_flags.stats.sn
```

Do not edit these files by hand.
