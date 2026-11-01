# samtools reference outputs

`all_flags.flagstat` is `samtools flagstat` 1.22.1 output (1.24 prints the same) for the records `all_flags_sam_lines()` in `tests/integration_test.rs` writes. To regenerate:

```bash
{ printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n'
  awk 'BEGIN{for(f=0;f<4096;f++){printf "f%da\t%d\tchr1\t10\t60\t4M\t=\t10\t0\tACGT\tIIII\n",f,f; printf "f%db\t%d\tchr1\t10\t60\t4M\tchr2\t10\t0\tACGT\tIIII\n",f,f}; for(f=0;f<4096;f++) printf "f%dc\t%d\tchr2\t10\t3\t4M\tchr1\t10\t0\tACGT\tIIII\n",f,f}'
} > all_flags.sam
samtools flagstat all_flags.sam > tests/expected/samtools/all_flags.flagstat
```
