# TODO: Generate large Qualimap reference output

## Problem

The large Qualimap reference output (`benchmark/qualimap/large/`) is missing.
Qualimap's RNA-seq mode name-sorts the entire BAM file, which requires ~50GB of
temporary disk space for the 10GB GM12878 BAM (185M reads). This filled the disk
during the first attempt.

## Prerequisites

- ~60GB free disk space
- Docker installed
- The Qualimap container: `quay.io/biocontainers/qualimap:2.3--hdfd78af_0`

## Steps

1. **Decompress the GTF** (Qualimap can't read gzipped GTF):

   ```bash
   gunzip -k benchmark/input/large/genes.gtf.gz
   ```

2. **Run Qualimap**:

   ```bash
   mkdir -p benchmark/qualimap/large
   docker run --rm \
     -v "$(pwd)/benchmark:/data" \
     quay.io/biocontainers/qualimap:2.3--hdfd78af_0 \
     qualimap rnaseq \
       -bam /data/input/large/GM12878_REP1.markdup.sorted.bam \
       -gtf /data/input/large/genes.gtf \
       -outdir /data/qualimap_tmp_large \
       --java-mem-size=8G \
       -pe
   ```

   This will take a long time (1-2 hours). The name-sort phase is the bottleneck.

3. **Copy the key output files**:

   ```bash
   cp benchmark/qualimap_tmp_large/rnaseq_qc_results.txt benchmark/qualimap/large/
   cp "benchmark/qualimap_tmp_large/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt" benchmark/qualimap/large/
   ```

4. **Clean up**:

   ```bash
   rm -rf benchmark/qualimap_tmp_large
   rm -f benchmark/input/large/genes.gtf
   ```

5. **Verify** — compare with RustQC output:

   ```bash
   diff benchmark/qualimap/large/rnaseq_qc_results.txt \
        benchmark/RustQC/large/qualimap/rnaseq_qc_results.txt
   diff benchmark/qualimap/large/coverage_profile_along_genes_\(total\).txt \
        benchmark/RustQC/large/qualimap/coverage_profile_along_genes_\(total\).txt
   ```

6. **Commit**:

   ```bash
   git add benchmark/qualimap/large/
   git commit -m "Add large Qualimap reference output"
   git push
   ```

7. **Delete this file** once done.

## Small dataset reference

Already generated and committed at `benchmark/qualimap/small/`.
