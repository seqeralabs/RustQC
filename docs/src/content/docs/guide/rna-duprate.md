---
title: RNA-seq Duplicate Rate Analysis
description: Understanding PCR duplicate rates in RNA-seq data, the dupRadar model, and how RustQC implements expression-dependent duplication analysis.
---

RustQC's `rna` subcommand performs duplicate rate analysis as a function of gene
expression level. This is a Rust reimplementation of the
[dupRadar](https://bioconductor.org/packages/dupRadar/) Bioconductor R package,
producing identical output.

## What is dupRadar?

In DNA sequencing, PCR duplicates are straightforward to identify and remove:
reads with identical alignment coordinates are likely amplification artifacts.
In **RNA-seq**, the situation is more nuanced.

Highly expressed genes produce many RNA fragments. When sequencing depth is
sufficient, some of these independent fragments will share the same start and
end positions purely by chance. These are **natural duplicates** -- biologically
real, not PCR artifacts. Removing them would distort expression estimates.

**dupRadar** addresses this by analyzing duplication rates in the context of
expression level. Rather than applying a blanket duplicate removal, it models
the expected relationship between how highly a gene is expressed and what
fraction of its reads are marked as duplicates.

The key insight: at low expression levels, duplicates are almost certainly PCR
artifacts. At high expression levels, a baseline level of duplication is
expected and biologically meaningful.

## Why it matters

Duplicate rate analysis helps answer critical quality questions about an
RNA-seq library:

- **Was the library over-amplified?** If low-expression genes show high
  duplication rates, too many PCR cycles were used.
- **Is the library complex enough?** Low-complexity libraries produce high
  duplication across all expression levels.
- **Can I trust the expression estimates?** If duplication follows the expected
  expression-dependent pattern, the data is likely reliable.

This information is essential for experiment QC, troubleshooting library
preparation protocols, and deciding whether to include a sample in downstream
analysis.

## The logistic regression model

dupRadar fits a logistic regression to the relationship between expression level
and duplication rate. RustQC implements this same model using iteratively
reweighted least squares (IRLS).

The model:

```
duplication_rate = 1 / (1 + exp(-(intercept + slope * x)))
```

where `x` is `log10(reads per kilobase)`.

### Interpreting the parameters

**Intercept:** Controls the baseline duplication rate at zero expression.

- Values close to **0** (or negative) indicate good library quality -- genes
  with few reads have few duplicates.
- Values above **0.5** indicate significant PCR over-amplification -- even
  lowly expressed genes have high duplication rates.

**Slope:** Controls how steeply duplication rises with expression.

- Positive values mean duplication increases with expression, which is the
  **expected biological pattern**.
- Very high slope values with a low intercept indicate a well-prepared library
  where duplication is driven by expression level rather than PCR artifacts.
- A low slope combined with a high intercept suggests uniform over-amplification.

### Practical interpretation

| Intercept | Slope | Interpretation |
|-----------|-------|----------------|
| Low (~0) | Moderate-high | Good library. Duplicates are expression-dependent. |
| High (>0.5) | Low | Over-amplified. Duplicates everywhere. |
| High (>0.5) | High | Over-amplified, but expression signal still visible. |
| Low (~0) | Very low | Good complexity, very little duplication at any level. |

## How RustQC implements this

RustQC performs the full dupRadar analysis in a single pass over each BAM file:

1. **Read counting:** Reads are assigned to genes based on the GTF annotation,
   tracking both total counts and non-duplicate counts for unique and
   multi-mapper modes.
2. **Matrix construction:** The 14-column duplication matrix is built, computing
   duplication rates, RPK, and RPKM for each gene.
3. **Model fitting:** A logistic regression is fitted via IRLS to the
   expression-vs-duplication data.
4. **Visualization:** Three diagnostic plots are generated (density scatter,
   boxplot, histogram).
5. **featureCounts output:** In the same pass, featureCounts-compatible files
   and biotype counts are produced.

This single-pass architecture is a key performance advantage over the R workflow,
which requires running featureCounts as a separate step before dupRadar.

## References

- **dupRadar:** Sayols S, Scherzinger D, Klein H. dupRadar: a Bioconductor
  package for the assessment of PCR artifacts in RNA-Seq data. *BMC
  Bioinformatics*. 2016;17(1):428.
  [Bioconductor page](https://bioconductor.org/packages/dupRadar/)
- **featureCounts:** Liao Y, Smyth GK, Shi W. featureCounts: an efficient
  general purpose program for assigning sequence reads to genomic features.
  *Bioinformatics*. 2014;30(7):923-930.
  [Subread/featureCounts](http://subread.sourceforge.net/)
