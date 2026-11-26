//! Number formatting that reproduces Qualimap's Java output.
//!
//! Qualimap 2.3 formats every number in `rnaseq_qc_results.txt` and in the
//! summary tables of `qualimapReport.html` through `StringUtilsSwing`
//! (`org.bioinfo.ngs.qc.qualimap.gui.utils`), which calls
//! `org.bioinfo.commons.utils.StringUtils.decimalFormat(value, pattern)`, i.e.
//! `new java.text.DecimalFormat(pattern).format(value)`:
//!
//! | method               | pattern                             | example      |
//! |----------------------|-------------------------------------|--------------|
//! | `formatLong`         | `###,###,####,###,…,###.##`         | `1,234,567`  |
//! | `formatDecimal`      | `###,###,###,###,###,###,###.##`    | `1,234.5`    |
//! | `formatPercentage`   | `###.##` followed by `%`            | `60%`        |
//!
//! `DecimalFormat` with these patterns prints at most two fraction digits,
//! drops trailing zeros (so whole numbers have no decimal point), always prints
//! at least one integer digit, groups the integer part by three (grouping size
//! is the last group in the pattern; `###.##` has no grouping) and rounds
//! with `RoundingMode.HALF_EVEN`.

/// Doubles at or above 2^53 are whole numbers whose exact decimal expansion
/// can have more digits than Java prints.
const EXACT_INTEGER_LIMIT: f64 = 9_007_199_254_740_992.0;

/// Java `DecimalFormat` output for NaN (CLDR locale data, Java 9+).
///
/// The nf-core Qualimap image (`biocontainers/qualimap:2.3--hdfd78af_0`)
/// ships OpenJDK 21.0.1 with `LANG=C.UTF-8`, so this and [`JAVA_INFINITY`]
/// are what it prints; Java 18+ also writes files as UTF-8 by default.
const JAVA_NAN: &str = "NaN";

/// Java `DecimalFormat` output for infinity.
const JAVA_INFINITY: &str = "\u{221e}";

/// Qualimap `StringUtilsSwing.formatLong`: an integer with `,` every three digits.
///
/// `formatLong` passes the `long` to `decimalFormat(double, String)`, so counts
/// above 2^53 lose precision exactly as a double does.
pub(crate) fn format_long(n: u64) -> String {
    decimal_format(n as f64, true)
}

/// Qualimap `StringUtilsSwing.formatDecimal`: up to two decimals, grouped by three.
pub(crate) fn format_decimal(v: f64) -> String {
    decimal_format(v, true)
}

/// Qualimap `StringUtilsSwing.formatPercentage`: up to two decimals and a `%` sign.
pub(crate) fn format_percentage(v: f64) -> String {
    let mut s = decimal_format(v, false);
    s.push('%');
    s
}

/// A percentage computed the way Qualimap does it: `(100. * part) / total`.
///
/// The operand order matters for the last bit of the result, and therefore for
/// how a value close to a rounding boundary is printed.
pub(crate) fn percentage_of(part: u64, total: u64) -> f64 {
    (100.0 * part as f64) / total as f64
}

/// Java's `Double.toString(v)` (Java 19+), which Qualimap uses when it writes
/// chart data such as the coverage profiles in `raw_data_qualimapReport/`.
///
/// The digits are the shortest that round-trip, as Rust prints them. Values
/// in `[1e-3, 1e7)` are plain decimals with at least one fraction digit
/// (`0.0`, `12.0`, `0.0041`); others use computerized scientific notation
/// (`1.0E-4`, `1.2345678E7`).
pub(crate) fn java_double_to_string(v: f64) -> String {
    if v.is_nan() {
        return "NaN".to_string();
    }
    let sign = if v.is_sign_negative() { "-" } else { "" };
    if v.is_infinite() {
        return format!("{sign}Infinity");
    }
    let abs = v.abs();
    let (mantissa, exponent) = if abs == 0.0 || (1e-3..1e7).contains(&abs) {
        (abs.to_string(), None)
    } else {
        let sci = format!("{abs:e}");
        let (m, e) = sci.split_once('e').expect("exponent in {:e} output");
        (m.to_string(), Some(e.to_string()))
    };
    let point = if mantissa.contains('.') { "" } else { ".0" };
    match exponent {
        None => format!("{sign}{mantissa}{point}"),
        Some(e) => format!("{sign}{mantissa}{point}E{e}"),
    }
}

/// `new DecimalFormat(pattern).format(v)` for Qualimap's two decimal patterns:
/// at most two fraction digits, at least one integer digit, HALF_EVEN rounding,
/// optional grouping by three.
///
/// `DecimalFormat` rounds the double's `FloatingDecimal` digits (17 significant
/// digits at most) and corrects ties for whether those digits are above or
/// below the exact value (`java.text.DigitList.shouldRoundUp`). While a
/// hundredth is far larger than the spacing between doubles, as for every
/// percentage, ratio and bias Qualimap prints, that equals rounding the exact
/// binary value half to even, which is what Rust's `{:.2}` does. Whole numbers
/// from 2^53 up are printed from their shortest digits, as Java does.
fn decimal_format(v: f64, grouping: bool) -> String {
    if v.is_nan() {
        return JAVA_NAN.to_string();
    }
    // Java prints a minus sign for any negative value, including -0.0 and
    // values that round to zero.
    let sign = if v.is_sign_negative() { "-" } else { "" };
    if v.is_infinite() {
        return format!("{sign}{JAVA_INFINITY}");
    }
    let abs = v.abs();
    let text = if abs >= EXACT_INTEGER_LIMIT {
        format!("{abs}")
    } else {
        format!("{abs:.2}")
    };
    let (int_part, frac_part) = text.split_once('.').unwrap_or((&text, ""));
    let frac_part = frac_part.trim_end_matches('0');
    let int_part = if grouping {
        group_thousands(int_part)
    } else {
        int_part.to_string()
    };
    if frac_part.is_empty() {
        format!("{sign}{int_part}")
    } else {
        format!("{sign}{int_part}.{frac_part}")
    }
}

/// Insert `,` between groups of three integer digits.
fn group_thousands(digits: &str) -> String {
    let len = digits.len();
    let mut out = String::with_capacity(len + len / 3);
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (len - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Expected strings below were produced by Qualimap 2.3's own
    /// `StringUtilsSwing` (bioinfo-commons 0.10.1) on Java 21.
    #[test]
    fn format_long_groups_by_three() {
        assert_eq!(format_long(0), "0");
        assert_eq!(format_long(123), "123");
        assert_eq!(format_long(1234), "1,234");
        assert_eq!(format_long(1_234_567), "1,234,567");
        assert_eq!(format_long(1_000_000_000), "1,000,000,000");
        // Above 2^53 Java formats the count as a double.
        assert_eq!(format_long(9_007_199_254_740_993), "9,007,199,254,740,992");
        assert_eq!(format_long(u64::MAX >> 1), "9,223,372,036,854,776,000");
    }

    #[test]
    fn percentage_drops_trailing_zeros() {
        assert_eq!(format_percentage(60.0), "60%");
        assert_eq!(format_percentage(100.0), "100%");
        assert_eq!(format_percentage(0.0), "0%");
        assert_eq!(format_percentage(37.5), "37.5%");
        assert_eq!(format_percentage(86.77), "86.77%");
        assert_eq!(format_percentage(0.05), "0.05%");
        assert_eq!(format_percentage(1.0 / 3.0 * 100.0), "33.33%");
    }

    #[test]
    fn percentage_is_not_grouped() {
        assert_eq!(format_percentage(1234.5), "1234.5%");
    }

    #[test]
    fn decimal_drops_trailing_zeros_and_groups() {
        assert_eq!(format_decimal(0.5), "0.5");
        assert_eq!(format_decimal(1.0), "1");
        assert_eq!(format_decimal(1.3), "1.3");
        assert_eq!(format_decimal(0.44), "0.44");
        assert_eq!(format_decimal(2.0 / 3.0), "0.67");
        assert_eq!(format_decimal(1234.5), "1,234.5");
        assert_eq!(format_decimal(1e6 + 0.125), "1,000,000.12");
        assert_eq!(format_decimal(1e15), "1,000,000,000,000,000");
        assert_eq!(format_decimal(1e20), "100,000,000,000,000,000,000");
        assert_eq!(
            format_decimal(123_456_789_012_345_678.0),
            "123,456,789,012,345,680"
        );
    }

    #[test]
    fn rounding_is_half_even_on_the_exact_binary_value() {
        // Exact binary ties go to the even neighbour.
        assert_eq!(format_decimal(0.125), "0.12");
        assert_eq!(format_decimal(0.375), "0.38");
        assert_eq!(format_decimal(0.625), "0.62");
        assert_eq!(format_decimal(0.875), "0.88");
        assert_eq!(format_percentage(3.125), "3.12%");
        // Not ties: the double is just above or below the decimal literal.
        assert_eq!(format_decimal(0.015), "0.01"); // 0.01499999...
        assert_eq!(format_decimal(0.035), "0.04"); // 0.03500000...04
        assert_eq!(format_decimal(1.005), "1"); // 1.00499999...
        assert_eq!(format_decimal(2.675), "2.67"); // 2.67499999...
        assert_eq!(format_decimal(0.995), "0.99"); // 0.99499999...
        assert_eq!(format_decimal(9.995), "9.99"); // 9.99499999...
        assert_eq!(format_percentage(99.995), "100%"); // 99.99500000...03
                                                       // Carry into the integer part and a new group.
        assert_eq!(format_decimal(999.999), "1,000");
    }

    #[test]
    fn special_values_match_java() {
        assert_eq!(format_decimal(f64::NAN), "NaN");
        assert_eq!(format_percentage(f64::NAN), "NaN%");
        assert_eq!(format_decimal(f64::INFINITY), "\u{221e}");
        assert_eq!(format_decimal(f64::NEG_INFINITY), "-\u{221e}");
        assert_eq!(format_decimal(-0.0), "-0");
        assert_eq!(format_decimal(-0.004), "-0");
        assert_eq!(format_decimal(-1.5), "-1.5");
        assert_eq!(format_decimal(-1_234_567.891), "-1,234,567.89");
        assert_eq!(format_decimal(f64::from_bits(1)), "0");
    }

    #[test]
    fn strand_estimation_is_formatted_as_a_float() {
        // Qualimap divides in single precision: fwd / (float)(fwd + rev).
        let ssp = |fwd: u64, rev: u64| f64::from(fwd as f32 / (fwd + rev) as f32);
        assert_eq!(format_decimal(ssp(1, 1)), "0.5");
        assert_eq!(format_decimal(ssp(1, 7)), "0.12"); // 0.125f is a tie
        assert_eq!(format_decimal(ssp(3, 0)), "1");
    }

    #[test]
    fn double_to_string_matches_java() {
        assert_eq!(java_double_to_string(0.0), "0.0");
        assert_eq!(java_double_to_string(-0.0), "-0.0");
        assert_eq!(java_double_to_string(1.0), "1.0");
        assert_eq!(java_double_to_string(99.0), "99.0");
        assert_eq!(java_double_to_string(0.1), "0.1");
        assert_eq!(
            java_double_to_string(0.004131012097964001),
            "0.004131012097964001"
        );
        assert_eq!(java_double_to_string(1234.5), "1234.5");
        assert_eq!(java_double_to_string(0.001), "0.001");
        assert_eq!(java_double_to_string(0.000_999), "9.99E-4");
        assert_eq!(java_double_to_string(1e-4), "1.0E-4");
        assert_eq!(java_double_to_string(9_999_999.0), "9999999.0");
        assert_eq!(java_double_to_string(1e7), "1.0E7");
        assert_eq!(java_double_to_string(12_345_678.9), "1.23456789E7");
        assert_eq!(java_double_to_string(-2.5e-5), "-2.5E-5");
        assert_eq!(java_double_to_string(f64::NAN), "NaN");
        assert_eq!(java_double_to_string(f64::INFINITY), "Infinity");
        assert_eq!(java_double_to_string(f64::NEG_INFINITY), "-Infinity");
    }

    #[test]
    fn percentage_of_uses_qualimaps_operand_order() {
        assert_eq!(percentage_of(1, 3), (100.0 * 1.0) / 3.0);
        assert_eq!(format_percentage(percentage_of(6, 10)), "60%");
        assert_eq!(format_percentage(percentage_of(1639, 1889)), "86.77%");
    }
}
