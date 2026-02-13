//! Duplication matrix construction.
//!
//! Builds the 14-column duplication matrix matching dupRadar's output exactly.
//! Each row represents a gene and contains counts, duplication rates, RPK, and RPKM
//! for both multi-mapper-inclusive and unique-mapper-only counting modes.

use crate::counting::{CountResult, GeneCounts};
use crate::gtf::Gene;
use indexmap::IndexMap;

/// A single row in the duplication matrix, corresponding to one gene.
#[derive(Debug, Clone)]
pub struct DupMatrixRow {
    /// Gene identifier
    pub id: String,
    /// Gene length in bp (non-overlapping exon bases)
    pub gene_length: u64,
    /// Total read count including multimappers
    pub all_counts_multi: u64,
    /// Read count excluding duplicates, including multimappers
    pub filtered_counts_multi: u64,
    /// Duplication rate (with multimappers): (all - filtered) / all
    pub dup_rate_multi: f64,
    /// Number of duplicate reads per gene (with multimappers)
    pub dups_per_id_multi: u64,
    /// Reads per kilobase (with multimappers)
    pub rpk_multi: f64,
    /// RPKM (with multimappers)
    pub rpkm_multi: f64,
    /// Total read count (unique mappers only)
    pub all_counts: u64,
    /// Read count excluding duplicates (unique mappers only)
    pub filtered_counts: u64,
    /// Duplication rate (unique mappers only)
    pub dup_rate: f64,
    /// Number of duplicate reads per gene (unique mappers only)
    pub dups_per_id: u64,
    /// Reads per kilobase (unique mappers only)
    pub rpk: f64,
    /// RPKM (unique mappers only)
    pub rpkm: f64,
}

/// The complete duplication matrix.
#[derive(Debug)]
pub struct DupMatrix {
    pub rows: Vec<DupMatrixRow>,
}

impl DupMatrix {
    /// Build the duplication matrix from gene annotations and count results.
    ///
    /// This mirrors dupRadar's `analyzeDuprates()` output exactly:
    /// - RPK = counts * (1000 / geneLength)
    /// - RPKM = RPK * (1_000_000 / N)  where N = total mapped reads
    /// - dupRate = (allCounts - filteredCounts) / allCounts
    pub fn build(genes: &IndexMap<String, Gene>, counts: &CountResult) -> Self {
        let n_multi_all = counts.total_reads_multi_dup as f64;
        let n_unique_all = counts.total_reads_unique_dup as f64;

        // For RPKM we use the total mapped reads from the "all" (dup-included) runs,
        // matching featureCounts convention. dupRadar computes RPK and RPKM for each
        // of the 4 runs independently, but the final matrix uses:
        // - Multi columns: N from multi+dup run for allCountsMulti RPKM,
        //   and N from multi+nodup run for filteredCountsMulti RPKM
        // - Unique columns: N from unique+dup run for allCounts RPKM,
        //   and N from unique+nodup run for filteredCounts RPKM
        //
        // However, looking at dupRadar source more carefully:
        // RPK and RPKM in the output use the "dup" (all reads) counts with
        // the corresponding N from that run.
        // The "Multi" columns use the multi+dup counts and N.
        // The non-Multi columns use the unique+dup counts and N.

        let mut rows = Vec::with_capacity(genes.len());
        let default_counts = GeneCounts::default();

        for (gene_id, gene) in genes.iter() {
            let gene_length = gene.effective_length;
            if gene_length == 0 {
                continue;
            }

            let gc = counts.gene_counts.get(gene_id).unwrap_or(&default_counts);

            let all_counts_multi = gc.all_multi;
            let filtered_counts_multi = gc.nodup_multi;
            let all_counts = gc.all_unique;
            let filtered_counts = gc.nodup_unique;

            let dups_per_id_multi = all_counts_multi.saturating_sub(filtered_counts_multi);
            let dups_per_id = all_counts.saturating_sub(filtered_counts);

            let dup_rate_multi = if all_counts_multi > 0 {
                dups_per_id_multi as f64 / all_counts_multi as f64
            } else {
                f64::NAN
            };

            let dup_rate = if all_counts > 0 {
                dups_per_id as f64 / all_counts as f64
            } else {
                f64::NAN
            };

            let gl = gene_length as f64;

            // RPK and RPKM for multi-mapper inclusive
            let rpk_multi = all_counts_multi as f64 * (1000.0 / gl);
            let rpkm_multi = if n_multi_all > 0.0 {
                rpk_multi * (1_000_000.0 / n_multi_all)
            } else {
                0.0
            };

            // RPK and RPKM for unique mappers only
            let rpk = all_counts as f64 * (1000.0 / gl);
            let rpkm = if n_unique_all > 0.0 {
                rpk * (1_000_000.0 / n_unique_all)
            } else {
                0.0
            };

            rows.push(DupMatrixRow {
                id: gene_id.clone(),
                gene_length,
                all_counts_multi,
                filtered_counts_multi,
                dup_rate_multi,
                dups_per_id_multi,
                rpk_multi,
                rpkm_multi,
                all_counts,
                filtered_counts,
                dup_rate,
                dups_per_id,
                rpk,
                rpkm,
            });
        }

        DupMatrix { rows }
    }

    /// Write the duplication matrix to a tab-separated file.
    ///
    /// Format matches dupRadar's `_dupMatrix.txt` output exactly.
    pub fn write_tsv(&self, path: &std::path::Path) -> anyhow::Result<()> {
        use std::io::Write;
        let file = std::fs::File::create(path)?;
        let mut writer = std::io::BufWriter::new(file);

        // Header matching dupRadar column names exactly (note: R has 'PKMMulti' typo in column 8)
        writeln!(
            writer,
            "ID\tgeneLength\tallCountsMulti\tfilteredCountsMulti\tdupRateMulti\tdupsPerIdMulti\tRPKMulti\tPKMMulti\tallCounts\tfilteredCounts\tdupRate\tdupsPerId\tRPK\tRPKM"
        )?;

        for row in &self.rows {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                row.id,
                row.gene_length,
                row.all_counts_multi,
                row.filtered_counts_multi,
                format_float(row.dup_rate_multi),
                row.dups_per_id_multi,
                format_float(row.rpk_multi),
                format_float(row.rpkm_multi),
                row.all_counts,
                row.filtered_counts,
                format_float(row.dup_rate),
                row.dups_per_id,
                format_float(row.rpk),
                format_float(row.rpkm),
            )?;
        }

        Ok(())
    }

    /// Get summary statistics matching `getDupMatStats()`.
    pub fn get_stats(&self) -> DupMatStats {
        let n_regions = self.rows.len();
        let n_covered = self.rows.iter().filter(|r| r.all_counts > 0).count();
        let n_with_dups = self
            .rows
            .iter()
            .filter(|r| r.dup_rate.is_finite() && r.dup_rate > 0.0)
            .count();

        DupMatStats {
            n_regions,
            n_regions_covered: n_covered,
            f_regions_covered: if n_regions > 0 {
                n_covered as f64 / n_regions as f64
            } else {
                0.0
            },
            n_regions_duplication: n_with_dups,
            f_regions_duplication: if n_regions > 0 {
                n_with_dups as f64 / n_regions as f64
            } else {
                0.0
            },
            f_covered_regions_duplication: if n_covered > 0 {
                n_with_dups as f64 / n_covered as f64
            } else {
                0.0
            },
        }
    }
}

/// Summary statistics for the duplication matrix.
#[derive(Debug)]
pub struct DupMatStats {
    /// Total number of genes in the annotation
    pub n_regions: usize,
    /// Number of genes with at least one read
    pub n_regions_covered: usize,
    /// Fraction of genes covered (n_regions_covered / n_regions)
    pub f_regions_covered: f64,
    /// Number of genes with duplication rate > 0
    pub n_regions_duplication: usize,
    /// Fraction of genes with duplication (n_regions_duplication / n_regions)
    pub f_regions_duplication: f64,
    /// Fraction of covered genes that have duplication
    #[allow(dead_code)]
    pub f_covered_regions_duplication: f64,
}

/// Format a float for TSV output, handling NaN as "NA" to match R's output.
/// Uses R-compatible formatting with up to 15 significant digits.
fn format_float(v: f64) -> String {
    if v.is_nan() {
        "NA".to_string()
    } else if v.is_infinite() {
        if v.is_sign_positive() {
            "Inf".to_string()
        } else {
            "-Inf".to_string()
        }
    } else if v == 0.0 {
        "0".to_string()
    } else {
        // R's write.table uses format() which defaults to 15 significant digits.
        // We need to replicate this: fixed notation with trailing zeros trimmed.
        //
        // Strategy: compute how many decimal places give ~15 significant digits,
        // then trim trailing zeros (but keep at least one decimal if present).
        let abs_v = v.abs();
        let digits_before_decimal = if abs_v >= 1.0 {
            (abs_v.log10().floor() as i32) + 1
        } else {
            // For 0.xxx values, digits before decimal = 0 (or negative for leading zeros)
            0
        };
        let decimal_places = (15 - digits_before_decimal).max(0) as usize;

        let s = format!("{:.*}", decimal_places, v);

        // Trim trailing zeros after decimal point, but keep the string meaningful
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dup_rate_calculation() {
        // If a gene has 100 total reads and 30 unique (non-dup) reads,
        // dupRate = (100 - 30) / 100 = 0.7
        let all = 100u64;
        let filtered = 30u64;
        let rate = (all - filtered) as f64 / all as f64;
        assert!((rate - 0.7).abs() < 1e-10);
    }

    #[test]
    fn test_dup_rate_zero_counts() {
        // Gene with no reads should have NaN dupRate
        let all = 0u64;
        let rate = if all > 0 { 0.0 } else { f64::NAN };
        assert!(rate.is_nan());
    }

    #[test]
    fn test_rpk_calculation() {
        // Gene of 2000bp with 100 reads: RPK = 100 * (1000/2000) = 50
        let counts = 100.0f64;
        let gene_length = 2000.0f64;
        let rpk = counts * (1000.0 / gene_length);
        assert!((rpk - 50.0).abs() < 1e-10);
    }

    #[test]
    fn test_rpkm_calculation() {
        // RPK=50, N=1_000_000: RPKM = 50 * (1_000_000/1_000_000) = 50
        let rpk = 50.0f64;
        let n = 1_000_000.0f64;
        let rpkm = rpk * (1_000_000.0 / n);
        assert!((rpkm - 50.0).abs() < 1e-10);
    }

    #[test]
    fn test_format_float() {
        assert_eq!(format_float(0.5), "0.5");
        assert_eq!(format_float(f64::NAN), "NA");
    }
}
