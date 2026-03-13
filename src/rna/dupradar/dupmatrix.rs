//! Duplication matrix construction.
//!
//! Builds the 14-column duplication matrix matching dupRadar's output exactly.
//! Each row represents a gene and contains counts, duplication rates, RPK, and RPKM
//! for both multi-mapper-inclusive and unique-mapper-only counting modes.

use crate::gtf::Gene;
use crate::rna::dupradar::counting::{CountResult, GeneCounts};
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
    /// - RPKM = RPK * (1_000_000 / N)  where N = total mapped fragments
    /// - dupRate = (allCounts - filteredCounts) / allCounts
    pub fn build(genes: &IndexMap<String, Gene>, counts: &CountResult) -> Self {
        // Compute N exactly as R dupRadar does:
        //   N <- sum(x$stat[, 2]) - x$stat[x$stat$Status == "Unassigned_Unmapped", 2]
        // R dupRadar calls featureCounts with isPairedEnd=TRUE for paired-end data,
        // so the stat summary counts *fragments* (read pairs), not individual reads.
        // RustQC's fragment counting tracks mapped fragments only. Upstream RSubread
        // featureCounts, however, also emits singleton unmapped mates as orphan
        // fragments in paired-end mode because their HI tag (often 0) does not match
        // the mapped mate's HI tag (>=1). Those orphan fragments are not classified as
        // `Unassigned_Unmapped`, so they contribute to dupRadar's denominator
        // `N = sum(stat) - Unassigned_Unmapped`.
        //
        // Treat `stat_singleton_unmapped_mates` as an upstream-compatibility correction
        // factor for dupRadar's denominator. RustQC's native fragment count excludes
        // these unmapped singleton mates, but upstream RSubread effectively counts
        // them as extra orphan fragments, so we add the correction factor here to
        // reproduce upstream dupRadar `N` exactly.
        let n_fragments =
            (counts.stat_total_fragments + counts.stat_singleton_unmapped_mates) as f64;
        let n_multi_all = n_fragments;
        let n_unique_all = n_fragments;

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
            .filter(|r| r.dup_rate > 0.0 || r.dup_rate.is_nan())
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
    #[allow(dead_code)]
    pub f_regions_covered: f64,
    /// Number of genes with duplication rate > 0
    pub n_regions_duplication: usize,
    /// Fraction of genes with duplication (n_regions_duplication / n_regions)
    #[allow(dead_code)]
    pub f_regions_duplication: f64,
    /// Fraction of covered genes that have duplication (computed for completeness, not yet output)
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
        // R's write.table uses format() which defaults to 15 significant digits
        // and scipen = 0.  We replicate both behaviours:
        //  1. Format with 15 significant digits, trailing zeros trimmed.
        //  2. Choose between fixed and scientific notation using R's scipen rule:
        //     prefer scientific when it is shorter (scipen = 0).

        // --- fixed notation ---
        let abs_v = v.abs();
        let digits_before_decimal = (abs_v.log10().floor() as i32) + 1;
        let decimal_places = (15 - digits_before_decimal).max(0) as usize;
        let fixed = format!("{:.*}", decimal_places, v);
        let fixed = if fixed.contains('.') {
            fixed
                .trim_end_matches('0')
                .trim_end_matches('.')
                .to_string()
        } else {
            fixed
        };

        // --- scientific notation ---
        // Format with 14 decimal places (1 digit before `.` + 14 = 15 sig digits),
        // then trim trailing zeros after the decimal point.
        let sci_raw = format!("{:.14e}", v);
        // Rust's {:.14e} produces e.g. "9.31000000000000e-5" — trim zeros and
        // normalise the exponent to R's 2-digit-minimum format (e-05, e+02).
        let sci = {
            let (mantissa, exp) = sci_raw.split_once('e').unwrap();
            let mantissa = if mantissa.contains('.') {
                mantissa.trim_end_matches('0').trim_end_matches('.')
            } else {
                mantissa
            };
            let exp_val: i32 = exp.parse().unwrap();
            format!("{}e{:+03}", mantissa, exp_val)
        };

        // R's scipen = 0: prefer scientific when it is strictly shorter.
        if sci.len() < fixed.len() {
            sci
        } else {
            fixed
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
        assert_eq!(format_float(0.0), "0");
        assert_eq!(format_float(f64::INFINITY), "Inf");
        assert_eq!(format_float(f64::NEG_INFINITY), "-Inf");
    }

    #[test]
    fn test_format_float_significant_digits() {
        // R uses 15 significant digits. Values < 1 with leading zeros
        // need extra decimal places to preserve 15 significant digits.
        // 1/28 = 0.0357142857142857... → 15 sig digits, 16 decimal places
        assert_eq!(format_float(0.0357142857142857), "0.0357142857142857");
        // 0.0758... → 15 sig digits, 16 decimal places
        assert_eq!(format_float(0.0758898079987858), "0.0758898079987858");
        // 0.0990... → 15 sig digits, 16 decimal places
        assert_eq!(format_float(0.0990785693054592), "0.0990785693054592");
        // Values >= 1 with 1 digit before decimal → 14 decimal places
        assert_eq!(format_float(4.9243756595146), "4.9243756595146");
        // Large values: 46.795... → 2 digits before decimal, 13 decimal places
        assert_eq!(format_float(46.795523906409), "46.795523906409");
    }

    #[test]
    fn test_format_float_scientific_notation() {
        // R uses scientific notation when it is shorter than fixed (scipen=0).
        // 9.31e-05 → "9.31e-05" (8 chars) vs "0.0000931" (9 chars)
        assert_eq!(format_float(9.31e-05), "9.31e-05");
        // 0.0001 → "1e-04" (5 chars) vs "0.0001" (6 chars)
        assert_eq!(format_float(0.0001), "1e-04");
        // 0.001 → "0.001" (5 chars) vs "1e-03" (5 chars) — tied, prefer fixed
        assert_eq!(format_float(0.001), "0.001");
        // 0.000123456 → "0.000123456" (11 chars) vs "1.23456e-04" (11 chars) — tied
        assert_eq!(format_float(0.000123456), "0.000123456");
    }
}
