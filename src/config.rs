//! Configuration file support for dupRust.
//!
//! Supports an optional YAML configuration file that can provide settings
//! like chromosome name mappings between BAM and GTF references.

use anyhow::{Context, Result};
use serde::Deserialize;
use std::collections::HashMap;
use std::path::Path;

/// Top-level configuration structure.
///
/// Designed to be extensible — new sections can be added as optional fields
/// without breaking existing config files.
#[derive(Debug, Deserialize, Default)]
#[serde(default)]
pub struct Config {
    /// Prefix to prepend to BAM chromosome names before matching to GTF names.
    ///
    /// Applied before explicit chromosome_mapping lookups. For example, if the
    /// BAM has "1", "2", "X" and the GTF has "chr1", "chr2", "chrX", set:
    /// ```yaml
    /// chromosome_prefix: "chr"
    /// ```
    #[serde(default)]
    pub chromosome_prefix: Option<String>,

    /// Chromosome name mapping from GTF names to BAM names.
    ///
    /// Keys are chromosome names as they appear in the GTF file,
    /// values are the corresponding names in the BAM file.
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

    /// Build a reverse mapping: BAM chromosome name -> GTF chromosome name.
    ///
    /// The config file maps GTF -> BAM, but at lookup time we need BAM -> GTF
    /// (since the index is keyed by GTF names).
    /// Explicit chromosome_mapping entries take priority over the prefix.
    pub fn bam_to_gtf_mapping(&self) -> HashMap<String, String> {
        self.chromosome_mapping
            .iter()
            .map(|(gtf_name, bam_name)| (bam_name.clone(), gtf_name.clone()))
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_config() {
        let config: Config = serde_yaml::from_str("").unwrap();
        assert!(config.chromosome_mapping.is_empty());
        assert!(!config.has_chromosome_mapping());
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

        let reverse = config.bam_to_gtf_mapping();
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
}
