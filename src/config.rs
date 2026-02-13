//! Configuration file support for RustQC.
//!
//! Supports an optional YAML configuration file that can provide settings
//! like chromosome name mappings between alignment file and GTF references
//! and per-tool output configuration.

use anyhow::{Context, Result};
use serde::Deserialize;
use std::collections::HashMap;
use std::path::Path;

/// Top-level configuration structure.
///
/// Designed to be extensible — new sections can be added as optional fields
/// without breaking existing config files. Tool-specific settings are nested
/// under their tool name (e.g. `dupradar:`, `featurecounts:`).
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

    /// Strandedness inference output configuration.
    #[serde(default)]
    pub strandedness: StrandednessConfig,
}

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

/// Configuration for strandedness inference outputs.
///
/// Controls whether strandedness inference results (RSeQC infer_experiment.py
/// compatible) are generated.
///
/// Example:
/// ```yaml
/// strandedness:
///   infer_experiment: true
///   multiqc_strandedness: true
/// ```
#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct StrandednessConfig {
    /// Write the RSeQC-compatible infer_experiment.txt file.
    pub infer_experiment: bool,

    /// Write the MultiQC strandedness general stats file.
    pub multiqc_strandedness: bool,
}

impl Default for StrandednessConfig {
    fn default() -> Self {
        Self {
            infer_experiment: true,
            multiqc_strandedness: true,
        }
    }
}

impl Config {
    /// Load configuration from a YAML file.
    pub fn from_file(path: &Path) -> Result<Self> {
        let contents = std::fs::read_to_string(path)
            .with_context(|| format!("Failed to read config file: {}", path.display()))?;
        let config: Config = serde_yaml::from_str(&contents)
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

    /// Returns true if any strandedness inference output is enabled.
    pub fn any_strandedness_output(&self) -> bool {
        let s = &self.strandedness;
        s.infer_experiment || s.multiqc_strandedness
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
        let config: Config = serde_yaml::from_str("").unwrap();
        assert!(config.chromosome_mapping.is_empty());
        assert!(!config.has_chromosome_mapping());
        // Defaults: all outputs enabled
        assert!(config.dupradar.dup_matrix);
        assert!(config.featurecounts.counts_file);
        assert_eq!(config.featurecounts.biotype_attribute, "gene_biotype");
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
        let config: Config = serde_yaml::from_str(yaml).unwrap();
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
        let config: Config = serde_yaml::from_str(yaml).unwrap();
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
        let config: Config = serde_yaml::from_str(yaml).unwrap();
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
        let config: Config = serde_yaml::from_str(yaml).unwrap();
        assert!(!config.any_dupradar_output());
    }
}
