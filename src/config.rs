//! Configuration file support for RustQC.
//!
//! Supports an optional YAML configuration file that can provide settings
//! like chromosome name mappings between alignment file and GTF references,
//! per-tool output configuration, and tool enable/disable toggles.

use crate::cli::Strandedness;
use anyhow::{Context, Result};
use serde::Deserialize;
use std::collections::HashMap;
use std::path::Path;

/// Top-level configuration structure.
///
/// Mirrors the CLI hierarchy: each subcommand (e.g. `rna`) has its own
/// configuration section. This allows the config file to grow naturally
/// as new subcommands are added.
///
/// Example:
/// ```yaml
/// rna:
///   chromosome_prefix: "chr"
///   bam_stat:
///     enabled: true
/// ```
#[derive(Debug, Deserialize, Default)]
#[serde(default)]
pub struct Config {
    /// RNA-Seq QC configuration (matches the `rna` subcommand).
    #[serde(default)]
    pub rna: RnaConfig,
}

/// RNA-Seq QC configuration.
///
/// Contains all settings for the `rustqc rna` subcommand. Tool-specific
/// settings are nested under their tool name (e.g. `dupradar:`,
/// `featurecounts:`, `bam_stat:`).
///
/// Designed to be extensible — new sections can be added as optional fields
/// without breaking existing config files.
#[derive(Debug, Deserialize, Default)]
#[serde(default)]
pub struct RnaConfig {
    /// Prefix to prepend to alignment file chromosome names before matching to GTF names.
    ///
    /// Applied before explicit chromosome_mapping lookups. For example, if the
    /// alignment file has "1", "2", "X" and the GTF has "chr1", "chr2", "chrX", set:
    /// ```yaml
    /// chromosome_prefix: "chr"
    /// ```
    #[serde(default)]
    pub chromosome_prefix: Option<String>,

    /// Chromosome name mapping from GTF names to alignment file names.
    ///
    /// Keys are chromosome names as they appear in the GTF file,
    /// values are the corresponding names in the alignment file (SAM/BAM/CRAM).
    /// Applied after chromosome_prefix (so explicit mappings can override
    /// the prefix for specific chromosomes like chrM -> MT).
    ///
    /// Example:
    /// ```yaml
    /// chromosome_mapping:
    ///   chr1: "1"
    ///   chr2: "2"
    ///   chrX: "X"
    /// ```
    pub chromosome_mapping: HashMap<String, String>,

    /// Library strandedness for strand-aware read counting.
    ///
    /// - `unstranded` = count reads on either strand
    /// - `forward` = forward stranded (read 1 maps to the transcript strand)
    /// - `reverse` = reverse stranded (read 2 maps to the transcript strand)
    ///
    /// The CLI `-s` / `--stranded` flag takes precedence over this setting.
    /// **Default:** `unstranded`.
    #[serde(default)]
    pub stranded: Option<Strandedness>,

    /// Enable paired-end mode.
    ///
    /// When `true`, read pairs are counted as a single fragment.
    /// The CLI `-p` / `--paired` flag takes precedence over this setting.
    /// **Default:** `false` (single-end mode).
    #[serde(default)]
    pub paired: Option<bool>,

    /// Override sample name for output filenames.
    ///
    /// By default, the sample name is derived from the BAM file stem. When set,
    /// this value is used instead for all output filenames.
    /// The CLI `--sample-name` flag takes precedence over this setting.
    #[serde(default)]
    pub sample_name: Option<String>,

    /// Write all output files to a flat directory (no subdirectories).
    ///
    /// By default (`false`), output files are organised into subdirectories by
    /// tool: `dupradar/`, `featurecounts/`, and `rseqc/{tool}/`. When `true`,
    /// all files are written directly to the output directory (legacy behaviour).
    /// The CLI `--flat-output` flag enables flat output regardless of this
    /// setting (either source being `true` produces flat output).
    #[serde(default)]
    pub flat_output: bool,

    /// dupRadar-specific output configuration.
    #[serde(default)]
    pub dupradar: DupradarConfig,

    /// featureCounts-compatible output configuration.
    #[serde(default)]
    pub featurecounts: FeatureCountsConfig,

    /// bam_stat tool configuration.
    #[serde(default)]
    pub bam_stat: BamStatConfig,

    /// infer_experiment tool configuration.
    #[serde(default)]
    pub infer_experiment: InferExperimentConfig,

    /// read_duplication tool configuration.
    #[serde(default)]
    pub read_duplication: ReadDuplicationConfig,

    /// read_distribution tool configuration.
    #[serde(default)]
    pub read_distribution: ReadDistributionConfig,

    /// junction_annotation tool configuration.
    #[serde(default)]
    pub junction_annotation: JunctionAnnotationConfig,

    /// junction_saturation tool configuration.
    #[serde(default)]
    pub junction_saturation: JunctionSaturationConfig,

    /// inner_distance tool configuration.
    #[serde(default)]
    pub inner_distance: InnerDistanceConfig,

    /// samtools flagstat-compatible output configuration.
    #[serde(default)]
    pub flagstat: FlagstatConfig,

    /// samtools idxstats-compatible output configuration.
    #[serde(default)]
    pub idxstats: IdxstatsConfig,

    /// TIN (Transcript Integrity Number) tool configuration.
    #[serde(default)]
    pub tin: TinConfig,

    /// samtools stats-compatible output configuration (SN section).
    #[serde(default)]
    pub samtools_stats: SamtoolsStatsConfig,

    /// preseq lc_extrap library complexity extrapolation configuration.
    #[serde(default)]
    pub preseq: PreseqConfig,

    /// Qualimap RNA-Seq QC configuration.
    #[serde(default)]
    pub qualimap: QualimapConfig,
}

// ============================================================================
// dupRadar configuration
// ============================================================================

/// Configuration for dupRadar outputs.
///
/// Controls which dupRadar output files are generated.
/// All outputs are enabled by default.
///
/// Example:
/// ```yaml
/// dupradar:
///   dup_matrix: true
///   intercept_slope: true
///   density_scatter_plot: true
///   boxplot: true
///   expression_histogram: true
///   multiqc_intercept: true
///   multiqc_curve: true
/// ```
#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct DupradarConfig {
    /// Write the duplication matrix TSV file.
    pub dup_matrix: bool,

    /// Write the intercept/slope fit results file.
    pub intercept_slope: bool,

    /// Generate the density scatter plot (PNG + SVG).
    pub density_scatter_plot: bool,

    /// Generate the duplication rate boxplot (PNG + SVG).
    pub boxplot: bool,

    /// Generate the expression histogram (PNG + SVG).
    pub expression_histogram: bool,

    /// Write the MultiQC intercept file.
    pub multiqc_intercept: bool,

    /// Write the MultiQC curve file.
    pub multiqc_curve: bool,
}

impl Default for DupradarConfig {
    fn default() -> Self {
        Self {
            dup_matrix: true,
            intercept_slope: true,
            density_scatter_plot: true,
            boxplot: true,
            expression_histogram: true,
            multiqc_intercept: true,
            multiqc_curve: true,
        }
    }
}

// ============================================================================
// featureCounts configuration
// ============================================================================

/// Configuration for featureCounts-compatible outputs.
///
/// Controls which featureCounts output files are generated and the
/// biotype counting behaviour.
///
/// Example:
/// ```yaml
/// featurecounts:
///   counts_file: true
///   summary_file: true
///   biotype_counts: true
///   biotype_counts_mqc: true
///   biotype_rrna_mqc: true
///   biotype_attribute: "gene_biotype"
/// ```
#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct FeatureCountsConfig {
    /// Write the featureCounts-compatible counts file.
    pub counts_file: bool,

    /// Write the featureCounts summary file (.summary).
    pub summary_file: bool,

    /// Write the biotype counts TSV file.
    pub biotype_counts: bool,

    /// Write the biotype counts MultiQC bargraph file.
    pub biotype_counts_mqc: bool,

    /// Write the biotype rRNA percentage MultiQC file.
    pub biotype_rrna_mqc: bool,

    /// GTF attribute name to use for biotype grouping.
    ///
    /// Defaults to `"gene_biotype"` (Ensembl convention).
    /// Use `"gene_type"` for GENCODE GTF files.
    pub biotype_attribute: String,
}

impl Default for FeatureCountsConfig {
    fn default() -> Self {
        Self {
            counts_file: true,
            summary_file: true,
            biotype_counts: true,
            biotype_counts_mqc: true,
            biotype_rrna_mqc: true,
            biotype_attribute: "gene_biotype".to_string(),
        }
    }
}

// ============================================================================
// RSeQC tool configurations
// ============================================================================

/// Configuration for bam_stat output.
///
/// Example:
/// ```yaml
/// bam_stat:
///   enabled: true
/// ```
#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct BamStatConfig {
    /// Whether to run bam_stat analysis. Defaults to true.
    pub enabled: bool,
}

impl Default for BamStatConfig {
    fn default() -> Self {
        Self { enabled: true }
    }
}

/// Configuration for infer_experiment output.
///
/// Example:
/// ```yaml
/// infer_experiment:
///   enabled: true
///   sample_size: 200000
/// ```
#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct InferExperimentConfig {
    /// Whether to run infer_experiment analysis. Defaults to true.
    pub enabled: bool,

    /// Maximum number of reads to sample for strandedness inference.
    /// Can be overridden by `--infer-experiment-sample-size` CLI flag.
    pub sample_size: Option<u64>,
}

impl Default for InferExperimentConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            sample_size: None,
        }
    }
}

/// Configuration for read_duplication output.
///
/// Example:
/// ```yaml
/// read_duplication:
///   enabled: true
/// ```
#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct ReadDuplicationConfig {
    /// Whether to run read_duplication analysis. Defaults to true.
    pub enabled: bool,
}

impl Default for ReadDuplicationConfig {
    fn default() -> Self {
        Self { enabled: true }
    }
}

/// Configuration for read_distribution output.
///
/// Example:
/// ```yaml
/// read_distribution:
///   enabled: true
/// ```
#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct ReadDistributionConfig {
    /// Whether to run read_distribution analysis. Defaults to true.
    pub enabled: bool,
}

impl Default for ReadDistributionConfig {
    fn default() -> Self {
        Self { enabled: true }
    }
}

/// Configuration for junction_annotation output.
///
/// Example:
/// ```yaml
/// junction_annotation:
///   enabled: true
///   min_intron: 50
/// ```
#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct JunctionAnnotationConfig {
    /// Whether to run junction_annotation analysis. Defaults to true.
    pub enabled: bool,

    /// Minimum intron size for junction filtering.
    /// Can be overridden by `--min-intron` CLI flag.
    pub min_intron: Option<u64>,
}

impl Default for JunctionAnnotationConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            min_intron: None,
        }
    }
}

/// Configuration for junction_saturation output.
///
/// Example:
/// ```yaml
/// junction_saturation:
///   enabled: true
///   min_coverage: 1
///   percentile_floor: 5
///   percentile_ceiling: 100
///   percentile_step: 5
/// ```
#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct JunctionSaturationConfig {
    /// Whether to run junction_saturation analysis. Defaults to true.
    pub enabled: bool,

    /// Minimum read coverage to count a known junction.
    /// Can be overridden by `--junction-saturation-min-coverage` CLI flag.
    pub min_coverage: Option<u64>,

    /// Sampling start percentage.
    /// Can be overridden by `--junction-saturation-percentile-floor` CLI flag.
    pub percentile_floor: Option<u64>,

    /// Sampling end percentage.
    /// Can be overridden by `--junction-saturation-percentile-ceiling` CLI flag.
    pub percentile_ceiling: Option<u64>,

    /// Sampling step percentage.
    /// Can be overridden by `--junction-saturation-percentile-step` CLI flag.
    pub percentile_step: Option<u64>,

    /// Random seed for the observation shuffle in saturation sampling.
    /// When set via `--junction-saturation-seed`, replaces the default hard-coded seed (42).
    pub seed: Option<u64>,
}

impl Default for JunctionSaturationConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            min_coverage: None,
            percentile_floor: None,
            percentile_ceiling: None,
            percentile_step: None,
            seed: None,
        }
    }
}

/// Configuration for inner_distance output.
///
/// Example:
/// ```yaml
/// inner_distance:
///   enabled: true
///   sample_size: 1000000
///   lower_bound: -250
///   upper_bound: 250
///   step: 5
/// ```
#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct InnerDistanceConfig {
    /// Whether to run inner_distance analysis. Defaults to true.
    pub enabled: bool,

    /// Maximum number of read pairs to sample.
    /// Can be overridden by `--inner-distance-sample-size` CLI flag.
    pub sample_size: Option<u64>,

    /// Lower bound of the inner distance histogram.
    /// Can be overridden by `--inner-distance-lower-bound` CLI flag.
    pub lower_bound: Option<i64>,

    /// Upper bound of the inner distance histogram.
    /// Can be overridden by `--inner-distance-upper-bound` CLI flag.
    pub upper_bound: Option<i64>,

    /// Bin width for the inner distance histogram.
    /// Can be overridden by `--inner-distance-step` CLI flag.
    pub step: Option<i64>,
}

impl Default for InnerDistanceConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            sample_size: None,
            lower_bound: None,
            upper_bound: None,
            step: None,
        }
    }
}

// ============================================================================
// samtools-compatible output configurations
// ============================================================================

/// Configuration for samtools flagstat-compatible output.
///
/// When enabled, produces a file matching `samtools flagstat` output format,
/// which is parseable by MultiQC.
///
/// Example:
/// ```yaml
/// flagstat:
///   enabled: true
/// ```
#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct FlagstatConfig {
    /// Whether to generate flagstat output. Defaults to true.
    pub enabled: bool,
}

impl Default for FlagstatConfig {
    fn default() -> Self {
        Self { enabled: true }
    }
}

/// Configuration for TIN (Transcript Integrity Number) analysis.
///
/// TIN measures per-transcript coverage uniformity using Shannon entropy.
/// Values range from 0 (completely degraded) to 100 (perfectly uniform).
///
/// Example:
/// ```yaml
/// tin:
///   enabled: true
///   sample_size: 100
///   min_coverage: 10
/// ```
#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct TinConfig {
    /// Whether to run TIN analysis. Defaults to true.
    pub enabled: bool,
    /// Number of equally-spaced positions to sample per transcript.
    pub sample_size: Option<u32>,
    /// Minimum number of reads covering a transcript to compute TIN.
    pub min_coverage: Option<u32>,
    /// Random seed for reproducible TIN results. When set, the internal
    /// hash state used for read-start tracking is seeded deterministically,
    /// ensuring identical output across runs. Without a seed the default
    /// random hash state may produce slightly different results for PE
    /// samples due to non-deterministic hash ordering.
    pub seed: Option<u64>,
}

impl Default for TinConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            sample_size: None,
            min_coverage: None,
            seed: None,
        }
    }
}

/// Configuration for Qualimap RNA-Seq QC.
///
/// When enabled, produces Qualimap-compatible output files including
/// `rnaseq_qc_results.txt`, coverage profiles (total/high/low), plots,
/// and an HTML report. Uses Qualimap-compatible counting logic:
/// enclosure-based gene assignment with M-only CIGAR parsing.
///
/// Example:
/// ```yaml
/// qualimap:
///   enabled: true
/// ```
#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct QualimapConfig {
    /// Whether to produce Qualimap RNA-Seq QC output. Defaults to true.
    pub enabled: bool,
}

impl Default for QualimapConfig {
    fn default() -> Self {
        Self { enabled: true }
    }
}

/// Configuration for samtools idxstats-compatible output.
///
/// When enabled, produces a file matching `samtools idxstats` output format,
/// which is parseable by MultiQC.
///
/// Example:
/// ```yaml
/// idxstats:
///   enabled: true
/// ```
#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct IdxstatsConfig {
    /// Whether to generate idxstats output. Defaults to true.
    pub enabled: bool,
}

impl Default for IdxstatsConfig {
    fn default() -> Self {
        Self { enabled: true }
    }
}

/// Configuration for samtools stats-compatible output (SN summary numbers section).
///
/// When enabled, produces a file matching the `SN` (Summary Numbers) section
/// of `samtools stats` output, which is parseable by MultiQC.
///
/// Example:
/// ```yaml
/// samtools_stats:
///   enabled: true
/// ```
#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct SamtoolsStatsConfig {
    /// Whether to generate samtools stats SN output. Defaults to true.
    pub enabled: bool,
}

impl Default for SamtoolsStatsConfig {
    fn default() -> Self {
        Self { enabled: true }
    }
}

// ============================================================================
// preseq lc_extrap configuration
// ============================================================================

/// Configuration for preseq lc_extrap library complexity extrapolation.
///
/// Estimates the expected number of distinct molecules as a function of
/// sequencing depth using Good-Toulmin rational function extrapolation
/// with bootstrap confidence intervals.
///
/// Example:
/// ```yaml
/// preseq:
///   enabled: true
///   max_extrap: 10000000000
///   step_size: 1000000
///   n_bootstraps: 100
///   confidence_level: 0.95
///   seed: 1
///   max_terms: 100
///   max_segment_length: 100000000
///   defects: false
/// ```
#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct PreseqConfig {
    /// Whether to run preseq lc_extrap analysis. Defaults to true.
    pub enabled: bool,

    /// Maximum extrapolation depth in total reads. Defaults to 1e10.
    pub max_extrap: f64,

    /// Step size between extrapolation points (in reads). Defaults to 1e6.
    pub step_size: f64,

    /// Number of bootstrap replicates for confidence intervals. Defaults to 100.
    pub n_bootstraps: u32,

    /// Confidence level for bootstrap intervals (e.g. 0.95 for 95%). Defaults to 0.95.
    pub confidence_level: f64,

    /// Random seed for bootstrap reproducibility. Defaults to 408
    /// (matching upstream preseq v3.2.0).
    pub seed: u64,

    /// Maximum number of terms in the power series / continued fraction. Defaults to 100.
    pub max_terms: usize,

    /// Maximum merged PE fragment length (bp). Defaults to 100000000 (effectively
    /// unlimited). Merged PE fragments longer than this are split back into
    /// individual reads. Corresponds to preseq's `-seg_len` option.
    pub max_segment_length: i64,

    /// Use the defects model for extrapolation. Defaults to false.
    ///
    /// When true, uses a modified rational function approximation that can
    /// handle certain problematic histograms where the standard method fails.
    pub defects: bool,
}

impl Default for PreseqConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            max_extrap: 1e10,
            step_size: 1e6,
            n_bootstraps: 100,
            confidence_level: 0.95,
            seed: 408,
            max_terms: 100,
            max_segment_length: 100_000_000,
            defects: false,
        }
    }
}

// ============================================================================
// Config implementation
// ============================================================================

impl Config {
    /// Load configuration from a YAML file.
    pub fn from_file(path: &Path) -> Result<Self> {
        let contents = std::fs::read_to_string(path)
            .with_context(|| format!("Failed to read config file: {}", path.display()))?;
        let config: Config = serde_yaml_ng::from_str(&contents)
            .with_context(|| format!("Failed to parse config file: {}", path.display()))?;
        Ok(config)
    }
}

impl RnaConfig {
    /// Build a reverse mapping: alignment chromosome name -> GTF chromosome name.
    ///
    /// The config file maps GTF -> alignment, but at lookup time we need
    /// alignment -> GTF (since the index is keyed by GTF names).
    /// Explicit chromosome_mapping entries take priority over the prefix.
    pub fn alignment_to_gtf_mapping(&self) -> HashMap<String, String> {
        self.chromosome_mapping
            .iter()
            .map(|(gtf_name, aln_name)| (aln_name.clone(), gtf_name.clone()))
            .collect()
    }

    /// Returns the chromosome prefix if configured.
    pub fn chromosome_prefix(&self) -> Option<&str> {
        self.chromosome_prefix.as_deref()
    }

    /// Returns true if there is any chromosome name remapping configured
    /// (either an explicit mapping or a prefix).
    pub fn has_chromosome_mapping(&self) -> bool {
        !self.chromosome_mapping.is_empty() || self.chromosome_prefix.is_some()
    }

    /// Returns true if any featureCounts output is enabled.
    pub fn any_featurecounts_output(&self) -> bool {
        let fc = &self.featurecounts;
        fc.counts_file || fc.summary_file
    }

    /// Returns true if any biotype output is enabled.
    pub fn any_biotype_output(&self) -> bool {
        let fc = &self.featurecounts;
        fc.biotype_counts || fc.biotype_counts_mqc || fc.biotype_rrna_mqc
    }

    /// Returns true if any dupRadar output is enabled.
    pub fn any_dupradar_output(&self) -> bool {
        let dr = &self.dupradar;
        dr.dup_matrix
            || dr.intercept_slope
            || dr.density_scatter_plot
            || dr.boxplot
            || dr.expression_histogram
            || dr.multiqc_intercept
            || dr.multiqc_curve
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- Tests for the top-level Config (rna: wrapper) ---

    #[test]
    fn test_top_level_rna_wrapper() {
        let yaml = r#"
rna:
  chromosome_prefix: "chr"
  stranded: reverse
  paired: true
  bam_stat:
    enabled: false
"#;
        let config: Config = serde_yaml_ng::from_str(yaml).unwrap();
        assert_eq!(config.rna.chromosome_prefix, Some("chr".to_string()));
        assert_eq!(config.rna.stranded, Some(Strandedness::Reverse));
        assert_eq!(config.rna.paired, Some(true));
        assert!(!config.rna.bam_stat.enabled);
    }

    #[test]
    fn test_empty_top_level_config() {
        let config: Config = serde_yaml_ng::from_str("").unwrap();
        // rna section defaults to RnaConfig::default()
        assert!(config.rna.chromosome_mapping.is_empty());
        assert!(config.rna.bam_stat.enabled);
        assert!(config.rna.preseq.enabled);
    }

    #[test]
    fn test_unknown_top_level_fields_ignored() {
        let yaml = r#"
rna:
  chromosome_prefix: "chr"
future_subcommand:
  key: value
"#;
        let config: Config = serde_yaml_ng::from_str(yaml).unwrap();
        assert_eq!(config.rna.chromosome_prefix, Some("chr".to_string()));
    }

    // --- Tests for RnaConfig (inner struct) ---

    #[test]
    fn test_empty_rna_config() {
        let config: RnaConfig = serde_yaml_ng::from_str("").unwrap();
        assert!(config.chromosome_mapping.is_empty());
        assert!(!config.has_chromosome_mapping());
        // stranded/paired default to None (defer to CLI)
        assert_eq!(config.stranded, None);
        assert_eq!(config.paired, None);
        // flat_output defaults to false (nested subdirectories)
        assert!(!config.flat_output);
        // Defaults: all outputs enabled
        assert!(config.dupradar.dup_matrix);
        assert!(config.featurecounts.counts_file);
        assert_eq!(config.featurecounts.biotype_attribute, "gene_biotype");
        // RSeQC tools all enabled by default
        assert!(config.bam_stat.enabled);
        assert!(config.infer_experiment.enabled);
        assert!(config.read_duplication.enabled);
        assert!(config.read_distribution.enabled);
        assert!(config.junction_annotation.enabled);
        assert!(config.junction_saturation.enabled);
        assert!(config.inner_distance.enabled);
        // samtools-compatible outputs all enabled by default
        assert!(config.flagstat.enabled);
        assert!(config.idxstats.enabled);
        assert!(config.samtools_stats.enabled);
        // preseq enabled by default with standard defaults
        assert!(config.preseq.enabled);
        assert!((config.preseq.max_extrap - 1e10).abs() < 1.0);
        assert!((config.preseq.step_size - 1e6).abs() < 1.0);
        assert_eq!(config.preseq.n_bootstraps, 100);
        assert!((config.preseq.confidence_level - 0.95).abs() < 1e-10);
        assert_eq!(config.preseq.seed, 408);
        assert_eq!(config.preseq.max_terms, 100);
        assert!(!config.preseq.defects);
    }

    #[test]
    fn test_stranded_paired_config() {
        // Defaults: None (defer to CLI)
        let config: RnaConfig = serde_yaml_ng::from_str("").unwrap();
        assert_eq!(config.stranded, None);
        assert_eq!(config.paired, None);

        // Explicit values
        let yaml = "stranded: reverse\npaired: true\n";
        let config: RnaConfig = serde_yaml_ng::from_str(yaml).unwrap();
        assert_eq!(config.stranded, Some(Strandedness::Reverse));
        assert_eq!(config.paired, Some(true));

        // Unstranded is a valid explicit value
        let yaml = "stranded: unstranded\npaired: false\n";
        let config: RnaConfig = serde_yaml_ng::from_str(yaml).unwrap();
        assert_eq!(config.stranded, Some(Strandedness::Unstranded));
        assert_eq!(config.paired, Some(false));

        // Forward stranded
        let yaml = "stranded: forward\n";
        let config: RnaConfig = serde_yaml_ng::from_str(yaml).unwrap();
        assert_eq!(config.stranded, Some(Strandedness::Forward));
    }

    #[test]
    fn test_flat_output_config() {
        let yaml = "flat_output: true\n";
        let config: RnaConfig = serde_yaml_ng::from_str(yaml).unwrap();
        assert!(config.flat_output);

        let yaml = "flat_output: false\n";
        let config: RnaConfig = serde_yaml_ng::from_str(yaml).unwrap();
        assert!(!config.flat_output);
    }

    #[test]
    fn test_chromosome_mapping() {
        let yaml = r#"
chromosome_mapping:
  chr1: "1"
  chr2: "2"
  chrX: "X"
  chrM: "MT"
"#;
        let config: RnaConfig = serde_yaml_ng::from_str(yaml).unwrap();
        assert_eq!(config.chromosome_mapping.len(), 4);
        assert_eq!(config.chromosome_mapping.get("chr1").unwrap(), "1");
        assert_eq!(config.chromosome_mapping.get("chrM").unwrap(), "MT");

        let reverse = config.alignment_to_gtf_mapping();
        assert_eq!(reverse.get("1").unwrap(), "chr1");
        assert_eq!(reverse.get("MT").unwrap(), "chrM");
    }

    #[test]
    fn test_unknown_rna_fields_ignored() {
        let yaml = r#"
chromosome_mapping:
  chr1: "1"
future_setting: true
another_section:
  key: value
"#;
        let config: RnaConfig = serde_yaml_ng::from_str(yaml).unwrap();
        assert_eq!(config.chromosome_mapping.len(), 1);
    }

    #[test]
    fn test_nested_tool_config() {
        let yaml = r#"
dupradar:
  dup_matrix: true
  boxplot: false
featurecounts:
  counts_file: true
  summary_file: false
  biotype_attribute: "gene_type"
"#;
        let config: RnaConfig = serde_yaml_ng::from_str(yaml).unwrap();
        assert!(config.dupradar.dup_matrix);
        assert!(!config.dupradar.boxplot);
        assert!(config.featurecounts.counts_file);
        assert!(!config.featurecounts.summary_file);
        assert_eq!(config.featurecounts.biotype_attribute, "gene_type");
    }

    #[test]
    fn test_disable_all_dupradar() {
        let yaml = r#"
dupradar:
  dup_matrix: false
  intercept_slope: false
  density_scatter_plot: false
  boxplot: false
  expression_histogram: false
  multiqc_intercept: false
  multiqc_curve: false
"#;
        let config: RnaConfig = serde_yaml_ng::from_str(yaml).unwrap();
        assert!(!config.any_dupradar_output());
    }

    #[test]
    fn test_disable_rseqc_tools() {
        let yaml = r#"
bam_stat:
  enabled: false
infer_experiment:
  enabled: false
read_duplication:
  enabled: false
read_distribution:
  enabled: false
junction_annotation:
  enabled: false
junction_saturation:
  enabled: false
inner_distance:
  enabled: false
"#;
        let config: RnaConfig = serde_yaml_ng::from_str(yaml).unwrap();
        assert!(!config.bam_stat.enabled);
        assert!(!config.infer_experiment.enabled);
        assert!(!config.read_duplication.enabled);
        assert!(!config.read_distribution.enabled);
        assert!(!config.junction_annotation.enabled);
        assert!(!config.junction_saturation.enabled);
        assert!(!config.inner_distance.enabled);
    }

    #[test]
    fn test_rseqc_tool_params() {
        let yaml = r#"
infer_experiment:
  enabled: true
  sample_size: 500000
junction_saturation:
  enabled: true
  min_coverage: 5
  percentile_floor: 10
  percentile_ceiling: 95
  percentile_step: 10
inner_distance:
  enabled: true
  sample_size: 2000000
  lower_bound: -500
  upper_bound: 500
  step: 10
"#;
        let config: RnaConfig = serde_yaml_ng::from_str(yaml).unwrap();
        assert_eq!(config.infer_experiment.sample_size, Some(500_000));
        assert_eq!(config.junction_saturation.min_coverage, Some(5));
        assert_eq!(config.junction_saturation.percentile_floor, Some(10));
        assert_eq!(config.junction_saturation.percentile_ceiling, Some(95));
        assert_eq!(config.junction_saturation.percentile_step, Some(10));
        assert_eq!(config.inner_distance.sample_size, Some(2_000_000));
        assert_eq!(config.inner_distance.lower_bound, Some(-500));
        assert_eq!(config.inner_distance.upper_bound, Some(500));
        assert_eq!(config.inner_distance.step, Some(10));
    }

    #[test]
    fn test_preseq_config() {
        let yaml = r#"
preseq:
  enabled: true
  seed: 1
  max_segment_length: 500000
  max_extrap: 5000000000
  step_size: 500000
  n_bootstraps: 50
  confidence_level: 0.99
  max_terms: 50
  defects: true
"#;
        let config: RnaConfig = serde_yaml_ng::from_str(yaml).unwrap();
        assert!(config.preseq.enabled);
        assert_eq!(config.preseq.seed, 1);
        assert_eq!(config.preseq.max_segment_length, 500_000);
        assert_eq!(config.preseq.max_extrap, 5_000_000_000.0);
        assert_eq!(config.preseq.step_size, 500_000.0);
        assert_eq!(config.preseq.n_bootstraps, 50);
        assert_eq!(config.preseq.confidence_level, 0.99);
        assert_eq!(config.preseq.max_terms, 50);
        assert!(config.preseq.defects);
    }
}
