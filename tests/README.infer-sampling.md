Run the native reproduction with:

```sh
cargo test --release --test infer_sampling
```

On unchanged v0.2.1/main, the explicit cap is ignored. The 300,000-record SE
witness has 156,000 sense and 44,000 antisense records in its first contig,
followed by 100,000 sense records on a larger second contig. Cap 1 must report
1.0000/0.0000 and cap 200,000 must report 0.7800/0.2200; the broken program uses
all records and reports 0.8533/0.1467. That population change can change a
downstream strandedness/QC decision.

The tests invoke Cargo's compile-time binary path. They generate indexed BAMs
and a reference-backed indexed CRAM with rust-htslib, including the FASTA index;
no noodles-only fixture, external executable or test dependency is required.
Expected positive-cap fractions and the exact zero frame were checked against
native RSeQC 5.0.4. Paired cap 3 stops midpair after three records. Other controls
cover filtered and ambiguous reads, supplementary eligibility, CLI/YAML
precedence, prefix/mapping resolution, zero, maximum cap, current unset/all-read
behavior and threads 1/2/4/8. A two-input control repeats both argument orders
to check that each file independently gets its cap. New commands clear ambient
RUSTQC_INFER_EXPERIMENT_SAMPLE_SIZE. Full flagstat equality protects other analyses.

The default is intentionally unchanged. A production pipeline benefits only
after using a fixed artifact and requesting the intended finite cap; this PR
does not change pipeline artifacts, flags or thresholds. Unknown-strand GTF
eligibility and whole-pipeline equivalence are separate concerns.
