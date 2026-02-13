---
title: RSeQC Outputs
description: Output files produced by RustQC's RSeQC-compatible quality control tools.
---

RustQC reimplements seven tools from the [RSeQC](https://rseqc.sourceforge.net/)
package. Each tool produces output files that match the format and content of the
original Python implementation.

All RSeQC tools use the input filename stem as a prefix for output files. For
example, `rustqc bam-stat sample.bam -o results/` produces
`results/sample.bam_stat.txt`.

## bam-stat

Basic alignment statistics from a single-pass BAM scan.

| File | Description |
|------|-------------|
| `{stem}.bam_stat.txt` | Formatted text report with total records, QC failures, duplicates, mapping quality distribution, splice reads, proper pairs, and more |

The output format matches `bam_stat.py` exactly, including the same section
headings and number formatting. Key metrics include:

- Total records, QC-failed, duplicates, non-primary, unmapped
- Unique and multi-mapped read counts at the configured MAPQ cutoff
- Read 1 / Read 2 counts, forward/reverse strand counts
- Splice / non-splice read counts
- Proper pair and paired-on-different-chromosome counts
- MAPQ distribution histogram

## infer-experiment

Library strandedness inference by sampling reads overlapping gene models.

| File | Description |
|------|-------------|
| `{stem}.infer_experiment.txt` | Strandedness fractions: failed-to-determine, and the two strand protocols |

The output reports the fraction of reads consistent with each strand protocol:

- **Fraction failed to determine** -- reads that could not be assigned to either protocol
- **Fraction "1++,1--,2+-,2-+"** (forward stranded) -- reads consistent with the same-strand protocol
- **Fraction "1+-,1-+,2++,2--"** (reverse stranded) -- reads consistent with the opposite-strand protocol

For paired-end data, the labels are `PairEnd` with `1++,1--,2+-,2-+` and
`1+-,1-+,2++,2--`. For single-end: `SingleEnd` with `++,--` and `+-,-+`.

**Interpreting results:**
- Both fractions near 50% = unstranded library
- First fraction near 100% = forward stranded (e.g., Ligation protocol)
- Second fraction near 100% = reverse stranded (e.g., dUTP protocol, most common)

## read-duplication

Position-based and sequence-based duplication rate histograms.

| File | Description |
|------|-------------|
| `{stem}.pos.DupRate.xls` | Position-based duplication histogram (TSV: Occurrence, UniqReadNumber, ReadNumber) |
| `{stem}.seq.DupRate.xls` | Sequence-based duplication histogram (TSV: Occurrence, UniqReadNumber, ReadNumber) |

Each file is a tab-separated table where each row represents a duplication
level (number of times a read was seen). The columns are:

- **Occurrence** -- the duplication count (1 = unique, 2 = seen twice, etc.)
- **UniqReadNumber** -- number of unique read groups at this duplication level
- **ReadNumber** -- total reads consumed (Occurrence x UniqReadNumber)

Position-based deduplication groups reads by alignment position (chromosome,
start, CIGAR-derived exon blocks). Sequence-based deduplication groups reads by
the actual read sequence.

## read-distribution

Classification of reads across genomic feature types.

| File | Description |
|------|-------------|
| `{stem}.read_distribution.txt` | Tabular report with total reads, total tags, and per-region breakdown |

The output includes:

- **Total Reads** and **Total Tags** (CIGAR M-block midpoints)
- A table of genomic regions with columns: Group, Total_bases, Tag_count, Tags/Kb
- Regions: CDS_Exons, 5'UTR_Exons, 3'UTR_Exons, Introns, TSS_up_1kb, TSS_up_5kb, TSS_up_10kb, TES_down_1kb, TES_down_5kb, TES_down_10kb
- Tags assigned to each region (with priority: CDS > UTR > Intron > Intergenic)

## junction-annotation

Splice junction classification against a reference gene model.

| File | Description |
|------|-------------|
| `{stem}.junction.xls` | TSV with all observed junctions: chrom, intron_start(0-based), intron_end(1-based), read_count, annotation_status |
| `{stem}.junction.bed` | BED12 file with color-coded junctions (red = known, green = partial novel, blue = complete novel) |
| `{stem}.junction_plot.r` | R script for generating splice event and junction pie charts |
| `{stem}.junction_annotation.txt` | Summary: total/known/partial novel/complete novel event and junction counts |

Junctions are classified by comparing splice sites (CIGAR N-operations)
against the reference BED12 gene model:

- **Known (Annotated)** -- both donor and acceptor sites match known introns
- **Partial novel** -- one splice site matches, the other is novel
- **Complete novel** -- neither splice site matches any known intron

## junction-saturation

Splice junction discovery rate at increasing sequencing depths.

| File | Description |
|------|-------------|
| `{stem}.junctionSaturation_plot.r` | R script for saturation curve plots |
| `{stem}.junctionSaturation_summary.txt` | TSV: percent_of_reads, known_junctions, novel_junctions, all_junctions |

The tool subsamples reads at configurable percentages (default: 5%, 10%, ...,
100%) and counts how many unique known and novel junctions are detected at each
level. This reveals whether sequencing depth is sufficient for comprehensive
junction detection. A saturated library will show a plateau in the curve; an
unsaturated library will show continuing growth.

## inner-distance

Fragment inner distance for paired-end RNA-seq libraries.

| File | Description |
|------|-------------|
| `{stem}.inner_distance.txt` | Per-pair detail: readpair_name, inner_distance, classification |
| `{stem}.inner_distance_freq.txt` | Histogram: lower_bound, upper_bound, count |
| `{stem}.inner_distance_plot.r` | R script for histogram and density plot |
| `{stem}.inner_distance_summary.txt` | Summary counts by pair classification |

The inner distance is defined as the gap between the end of read 1 and the
start of read 2 on the mRNA transcript. Negative values indicate read overlap.
Pairs are classified as:

- **sameTranscript=Yes, sameExon=Yes** -- both reads on the same exon
- **sameTranscript=Yes, sameExon=No** -- reads on the same transcript, different exons (distance calculated on mRNA)
- **sameTranscript=No** -- reads on different transcripts or no transcript overlap
- **readPairOverlap** -- reads overlap each other
- **nonExonic** -- one or both reads not on exons
- **sameChrom=No** -- reads on different chromosomes
- **unknownChromosome** -- chromosome not in the gene model

The histogram bins are configurable via `--lower-bound`, `--upper-bound`, and
`--step` (defaults: -250 to 250, step 5).

## Compatibility with RSeQC

All output files are designed to be drop-in replacements for the corresponding
RSeQC Python tool output. The generated R scripts produce identical plots when
executed. File formats, column names, and numeric precision match the Python
originals to facilitate migration from RSeQC to RustQC without downstream
pipeline changes.
