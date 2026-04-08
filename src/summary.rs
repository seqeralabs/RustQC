//! Machine-readable JSON summary output.
//!
//! Provides serializable structs that capture all key metrics from a RustQC
//! run. Written to a JSON file via `--json-summary` for automation and AI agents.

use serde::Serialize;

/// Top-level summary of a complete RustQC run.
#[derive(Debug, Serialize)]
pub struct RunSummary {
    /// RustQC version string.
    pub version: String,
    /// Short git commit hash.
    pub commit: String,
    /// Binary microarchitecture target (e.g. "x86-64-v3", "aarch64-neoverse-v1").
    pub binary_target: String,
    /// CPU features detected at runtime (e.g. ["AVX2", "SSE4.2", "POPCNT"]).
    pub cpu_features: Vec<String>,
    /// UTC timestamp when the run started.
    pub timestamp_start: String,
    /// UTC timestamp when the run finished.
    pub timestamp_end: String,
    /// Total wall-clock runtime in seconds.
    pub runtime_seconds: f64,
    /// Per-input-file summaries.
    pub inputs: Vec<InputSummary>,
}

/// Summary of a single BAM file's processing.
#[derive(Debug, Serialize)]
pub struct InputSummary {
    /// Path to the BAM/SAM/CRAM file.
    pub bam_file: String,
    /// Processing status: "success" or "failed".
    pub status: String,
    /// Error message if processing failed.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    /// Wall-clock time for this file in seconds.
    pub runtime_seconds: f64,
    /// Counting statistics (if successful).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub counting: Option<CountingSummary>,
    /// dupRadar summary (if successful and enabled).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub dupradar: Option<DupradarSummary>,
    /// List of output files written.
    pub outputs: Vec<OutputFile>,
}

/// Counting statistics for a single BAM file.
#[derive(Debug, Serialize)]
pub struct CountingSummary {
    /// Total reads/records in the BAM.
    pub total_reads: u64,
    /// Mapped reads.
    pub mapped_reads: u64,
    /// Total fragments (single-end reads or paired-end pairs).
    pub fragments: u64,
    /// Fragments assigned to exactly one gene.
    pub assigned: u64,
    /// Fragments overlapping no annotated gene.
    pub no_features: u64,
    /// Fragments overlapping multiple genes.
    pub ambiguous: u64,
    /// Duplicate-flagged reads.
    pub duplicates: u64,
    /// Multimapping reads (NH > 1).
    pub multimappers: u64,
    /// Percentage of fragments assigned.
    pub assigned_pct: f64,
    /// Duplicate rate as percentage of mapped reads.
    pub duplicate_pct: f64,
}

/// dupRadar summary for a single BAM file.
#[derive(Debug, Serialize)]
pub struct DupradarSummary {
    /// Total genes in annotation.
    pub total_genes: u64,
    /// Genes with at least one read.
    pub genes_with_reads: u64,
    /// Genes with at least one duplicate.
    pub genes_with_duplication: u64,
    /// Logistic regression intercept (if fit succeeded).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub intercept: Option<f64>,
    /// Logistic regression slope (if fit succeeded).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub slope: Option<f64>,
}

/// A single output file written during processing.
#[derive(Debug, Serialize)]
pub struct OutputFile {
    /// Tool name that produced this file.
    pub tool: String,
    /// Path to the output file.
    pub path: String,
}
