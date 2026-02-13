//! Configuration file support for RustQC.
//!
//! Supports an optional YAML configuration file that can provide settings
//! like chromosome name mappings between alignment file and GTF references,
//! per-tool output configuration, and tool enable/disable toggles.

use anyhow::{Context, Result};
use serde::Deserialize;
use std::collections::HashMap;
use std::path::Path;

/// Top-level configuration structure.
///
/// Designed to be extensible — new sections can be added as optional fields
/// without breaking existing config files. Tool-specific settings are nested
/// under their tool name (e.g. `dupradar:`, `featurecounts:`, `bam_stat:`).
#[derive(Debug, Deserialize, Default)]
#[serde(default)]
pub struct Config {
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
/// Requires a BED file to be provided via `--bed`.
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
///   plot: true
/// ```
#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct ReadDuplicationConfig {
    /// Whether to run read_duplication analysis. Defaults to true.
    pub enabled: bool,

    /// Generate the read duplication plot (PNG + SVG).
    pub plot: bool,
}

impl Default for ReadDuplicationConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            plot: true,
        }
    }
}

/// Configuration for read_distribution output.
///
/// Requires a BED file to be provided via `--bed`.
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
/// Requires a BED file to be provided via `--bed`.
///
/// Example:
/// ```yaml
/// junction_annotation:
///   enabled: true
///   min_intron: 50
///   splice_events_plot: true
///   splice_junction_plot: true
/// ```
#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct JunctionAnnotationConfig {
    /// Whether to run junction_annotation analysis. Defaults to true.
    pub enabled: bool,

    /// Minimum intron size for junction filtering.
    /// Can be overridden by `--min-intron` CLI flag.
    pub min_intron: Option<u64>,

    /// Generate the splice events pie chart (PNG + SVG).
    pub splice_events_plot: bool,

    /// Generate the splice junctions pie chart (PNG + SVG).
    pub splice_junction_plot: bool,
}

impl Default for JunctionAnnotationConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            min_intron: None,
            splice_events_plot: true,
            splice_junction_plot: true,
        }
    }
}

/// Configuration for junction_saturation output.
///
/// Requires a BED file to be provided via `--bed`.
///
/// Example:
/// ```yaml
/// junction_saturation:
///   enabled: true
///   min_coverage: 1
///   percentile_floor: 5
///   percentile_ceiling: 100
///   percentile_step: 5
///   plot: true
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

    /// Generate the junction saturation plot (PNG + SVG).
    pub plot: bool,
}

impl Default for JunctionSaturationConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            min_coverage: None,
            percentile_floor: None,
            percentile_ceiling: None,
            percentile_step: None,
            plot: true,
        }
    }
}

/// Configuration for inner_distance output.
///
/// Requires a BED file to be provided via `--bed`.
///
/// Example:
/// ```yaml
/// inner_distance:
///   enabled: true
///   sample_size: 1000000
///   lower_bound: -250
///   upper_bound: 250
///   step: 5
///   plot: true
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

    /// Generate the inner distance plot (PNG + SVG).
    pub plot: bool,
}

impl Default for InnerDistanceConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            sample_size: None,
            lower_bound: None,
            upper_bound: None,
            step: None,
            plot: true,
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
        let config: Config = serde_yml::from_str(&contents)
            .with_context(|| format!("Failed to parse config file: {}", path.display()))?;
        Ok(config)
    }

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

    #[test]
    fn test_empty_config() {
        let config: Config = serde_yml::from_str("").unwrap();
        assert!(config.chromosome_mapping.is_empty());
        assert!(!config.has_chromosome_mapping());
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
        // RSeQC plot toggles default to true
        assert!(config.read_duplication.plot);
        assert!(config.junction_annotation.splice_events_plot);
        assert!(config.junction_annotation.splice_junction_plot);
        assert!(config.junction_saturation.plot);
        assert!(config.inner_distance.plot);
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
        let config: Config = serde_yml::from_str(yaml).unwrap();
        assert_eq!(config.chromosome_mapping.len(), 4);
        assert_eq!(config.chromosome_mapping.get("chr1").unwrap(), "1");
        assert_eq!(config.chromosome_mapping.get("chrM").unwrap(), "MT");

        let reverse = config.alignment_to_gtf_mapping();
        assert_eq!(reverse.get("1").unwrap(), "chr1");
        assert_eq!(reverse.get("MT").unwrap(), "chrM");
    }

    #[test]
    fn test_unknown_fields_ignored() {
        let yaml = r#"
chromosome_mapping:
  chr1: "1"
future_setting: true
another_section:
  key: value
"#;
        let config: Config = serde_yml::from_str(yaml).unwrap();
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
        let config: Config = serde_yml::from_str(yaml).unwrap();
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
        let config: Config = serde_yml::from_str(yaml).unwrap();
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
        let config: Config = serde_yml::from_str(yaml).unwrap();
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
        let config: Config = serde_yml::from_str(yaml).unwrap();
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
}
