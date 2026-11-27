//! Plain number formatting for Qualimap's text and HTML output.
//!
//! Numbers have at most two decimals, trailing zeros (and a trailing `.`) are
//! dropped, and counts are grouped with `,` every three digits. Percentages are
//! not grouped and end in `%`.

/// An integer with `,` every three digits.
pub(crate) fn format_long(n: u64) -> String {
    group_thousands(&n.to_string())
}

/// At most two decimals, with `,` every three integer digits.
pub(crate) fn format_decimal(v: f64) -> String {
    format_number(v, true)
}

/// At most two decimals and a `%` sign, without grouping.
pub(crate) fn format_percentage(v: f64) -> String {
    format!("{}%", format_number(v, false))
}

fn format_number(v: f64, grouping: bool) -> String {
    if !v.is_finite() {
        return v.to_string();
    }
    let text = format!("{:.2}", v.abs());
    let (int_part, frac_part) = text.split_once('.').unwrap_or((&text, ""));
    let frac_part = frac_part.trim_end_matches('0');
    let int_part = if grouping {
        group_thousands(int_part)
    } else {
        int_part.to_string()
    };
    let sign = if v.is_sign_negative() { "-" } else { "" };
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

    #[test]
    fn format_long_groups_by_three() {
        assert_eq!(format_long(0), "0");
        assert_eq!(format_long(123), "123");
        assert_eq!(format_long(1234), "1,234");
        assert_eq!(format_long(1_234_567), "1,234,567");
        assert_eq!(format_long(1_000_000_000), "1,000,000,000");
        assert_eq!(format_long(u64::MAX), "18,446,744,073,709,551,615");
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
        assert_eq!(format_percentage(99.995), "100%");
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
        assert_eq!(format_decimal(2.0 / 3.0), "0.67");
        assert_eq!(format_decimal(1234.5), "1,234.5");
        assert_eq!(format_decimal(1e15), "1,000,000,000,000,000");
        assert_eq!(format_decimal(999.999), "1,000");
        assert_eq!(format_decimal(-1_234_567.891), "-1,234,567.89");
    }

    #[test]
    fn rounds_to_two_decimals() {
        assert_eq!(format_decimal(0.015), "0.01");
        assert_eq!(format_decimal(0.035), "0.04");
        assert_eq!(format_decimal(1.005), "1");
        assert_eq!(format_decimal(2.675), "2.67");
        assert_eq!(format_decimal(-0.004), "-0");
    }

    #[test]
    fn non_finite_values_are_parseable() {
        assert_eq!(format_decimal(f64::NAN), "NaN");
        assert_eq!(format_percentage(f64::NAN), "NaN%");
        assert_eq!(format_decimal(f64::INFINITY), "inf");
    }
}
