<h1 align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/public/RustQC-logo-darkbg.svg">
  <source media="(prefers-color-scheme: light)" srcset="docs/public/RustQC-logo.svg">
  <img width="500" src="docs/public/RustQC-logo.svg" alt="RustQC">
</picture>
</h1>

<h4 align="center">Fast genomics quality control tools for sequencing data, written in Rust.</h4>

<p align="center">
  <a href="https://github.com/seqeralabs/RustQC/actions/workflows/ci.yml"><img src="https://github.com/seqeralabs/RustQC/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
  <a href="https://github.com/seqeralabs/RustQC/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-MIT-green" alt="License"></a>
</p>

<p align="center">
  <a href="https://rustqc.netlify.app/">Documentation</a> &bull;
  <a href="https://rustqc.netlify.app/getting-started/quickstart/">Quickstart</a> &bull;
  <a href="https://rustqc.netlify.app/benchmarks/combined/">Benchmarks</a> &bull;
  <a href="https://github.com/seqeralabs/RustQC/releases">Releases</a>
</p>

---

**RustQC** is a suite of fast QC tools for sequencing data, compiled to a single static binary with no runtime dependencies.

<p align="center">
<picture>
   <source media="(prefers-color-scheme: dark)" srcset="docs/public/benchmarks/benchmark_dark.png">
   <img src="docs/public/benchmarks/benchmark_light.png" alt="Benchmark: RustQC ~14m 54s vs traditional tools ~15h 34m sequential (dupRadar + featureCounts + 8 RSeQC tools incl. TIN + preseq + samtools + Qualimap)" width="600">
</picture>
</p>

<p align="center"><em>Run time for a large paired-end RNA-seq BAM (~186M reads) on AWS.</em></p>

It currently includes:

- `rustqc rna` is a single-command RNA-Seq QC tool that runs all QC analyses in one pass. Designed to slot into the [nf-core/rnaseq pipeline](https://nf-co.re/rnaseq/), but works anywhere:

| Tool                | Reimplements                                                                            | Description                                                           |
| ------------------- | --------------------------------------------------------------------------------------- | --------------------------------------------------------------------- |
| dupRadar            | [dupRadar](https://github.com/ssayols/dupRadar)                                         | PCR duplicate rate vs. expression analysis with density scatter plots |
| featureCounts       | [featureCounts](http://subread.sourceforge.net/)                                        | Gene-level read counting with biotype summaries                       |
| bam_stat            | [RSeQC](https://rseqc.sourceforge.net/#bam-stat-py) `bam_stat.py`                       | Basic alignment statistics                                            |
| infer_experiment    | [RSeQC](https://rseqc.sourceforge.net/#infer-experiment-py) `infer_experiment.py`       | Library strandedness inference                                        |
| read_duplication    | [RSeQC](https://rseqc.sourceforge.net/#read-duplication-py) `read_duplication.py`       | Position- and sequence-based duplication histograms                   |
| read_distribution   | [RSeQC](https://rseqc.sourceforge.net/#read-distribution-py) `read_distribution.py`     | Read distribution across genomic features                             |
| junction_annotation | [RSeQC](https://rseqc.sourceforge.net/#junction-annotation-py) `junction_annotation.py` | Splice junction classification                                        |
| junction_saturation | [RSeQC](https://rseqc.sourceforge.net/#junction-saturation-py) `junction_saturation.py` | Splice junction saturation analysis                                   |
| inner_distance      | [RSeQC](https://rseqc.sourceforge.net/#inner-distance-py) `inner_distance.py`           | Paired-end inner distance distribution                                |
| TIN                 | [RSeQC](https://rseqc.sourceforge.net/#tin-py) `tin.py`                                 | Transcript Integrity Number                                           |
| preseq              | [preseq](http://smithlabresearch.org/software/preseq/) `lc_extrap`                      | Library complexity extrapolation                                      |
| Qualimap rnaseq     | [Qualimap](http://qualimap.conesalab.org/) `rnaseq`                                     | Gene body coverage, read origin, strand specificity                   |
| flagstat            | [samtools](http://www.htslib.org/) `flagstat`                                           | Alignment flag summary                                                |
| idxstats            | [samtools](http://www.htslib.org/) `idxstats`                                           | Per-chromosome read counts                                            |
| stats               | [samtools](http://www.htslib.org/) `stats`                                              | Full samtools stats output including all histogram sections           |

All outputs are format- and numerically identical to the upstream tools, and compatible with [MultiQC](https://multiqc.info/) for reporting.

## Quick start

```bash
# Install (Linux x86_64 example -- see docs for all platforms)
curl -fsSL https://github.com/seqeralabs/RustQC/releases/latest/download/rustqc-linux-x86_64.tar.gz | tar xz --strip-components=1
sudo mv ./rustqc /usr/local/bin/

# Run RNA-Seq QC
rustqc rna sample.markdup.bam --gtf genes.gtf --paired --outdir results/
```

```bash
# Or use Docker
docker run --rm -v "$PWD":/data ghcr.io/seqeralabs/rustqc:latest \
  rustqc rna /data/sample.markdup.bam --gtf /data/genes.gtf --outdir /data/results
```

See the [documentation](https://rustqc.netlify.app/) for full usage details, configuration options, output file descriptions, and benchmark results.

## License

MIT License. See [LICENSE](LICENSE) for details.
