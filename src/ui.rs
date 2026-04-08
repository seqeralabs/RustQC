//! Terminal UI: styled output, progress bars, and verbosity control.
//!
//! Centralises all user-facing terminal output for RustQC. Provides styled
//! headers, progress bars, summary boxes, and consistent formatting across
//! Normal, Verbose, and Quiet modes.

use console::Style;
use indicatif::{ProgressBar, ProgressStyle};
use std::time::Duration;

// ============================================================================
// Verbosity
// ============================================================================

/// Controls how much output is printed to the terminal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verbosity {
    /// Suppress all output except warnings and errors.
    Quiet,
    /// Default: header, config, progress, summary, output checklist.
    Normal,
    /// Show additional detail (per-file writes, intermediate stats, build steps).
    Verbose,
}

// ============================================================================
// Ui
// ============================================================================

/// Terminal UI handle. All user-facing output flows through this struct.
///
/// Methods are no-ops in `Quiet` mode (except `warn` and `error`, which always
/// print). `Verbose` mode shows additional detail beyond `Normal`.
///
/// **Note:** when processing multiple BAMs in parallel, output lines from
/// different files may interleave on stderr.
#[derive(Debug)]
pub struct Ui {
    /// Current verbosity level.
    verbosity: Verbosity,
    // Reusable styles
    /// Bold white for section headers.
    style_header: Style,
    /// Cyan for labels / keys.
    style_label: Style,
    /// Dim for paths and secondary info.
    style_dim: Style,
    /// Green for success markers.
    style_green: Style,
    /// Yellow for warnings.
    style_yellow: Style,
    /// Red for errors.
    style_red: Style,
    /// Bold for emphasis.
    style_bold: Style,
    /// Dim for box borders.
    style_box: Style,
}

impl Ui {
    /// Create a new UI handle with the given verbosity.
    pub fn new(verbosity: Verbosity) -> Self {
        Self {
            verbosity,
            style_header: Style::new().bold(),
            style_label: Style::new().cyan(),
            style_dim: Style::new().dim(),
            style_green: Style::new().green(),
            style_yellow: Style::new().yellow(),
            style_red: Style::new().red().bold(),
            style_bold: Style::new().bold(),
            style_box: Style::new().dim(),
        }
    }

    /// Whether we are in quiet mode.
    fn is_quiet(&self) -> bool {
        self.verbosity == Verbosity::Quiet
    }

    /// Whether we are in verbose mode.
    pub fn is_verbose(&self) -> bool {
        self.verbosity == Verbosity::Verbose
    }

    // ========================================================================
    // Header & config
    // ========================================================================

    /// Print the startup banner, optionally including a CPU info line.
    pub fn header(&self, version: &str, commit: &str, build: &str, cpu_info: Option<&str>) {
        if self.is_quiet() {
            return;
        }
        eprintln!();
        eprintln!(
            "  {} {}",
            self.style_header.apply_to(format!("RustQC v{version}")),
            self.style_dim
                .apply_to(format!("({commit}, built {build})")),
        );
        if let Some(info) = cpu_info {
            eprintln!("  {}", self.style_dim.apply_to(info));
        }
        eprintln!();
    }

    /// Print a key-value config line (e.g. "  Input:      sample.bam").
    pub fn config(&self, key: &str, value: &str) {
        if self.is_quiet() {
            return;
        }
        eprintln!(
            "  {:<13}{}",
            self.style_label.apply_to(format!("{key}:")),
            value,
        );
    }

    /// Print a blank line separator.
    pub fn blank(&self) {
        if self.is_quiet() {
            return;
        }
        eprintln!();
    }

    // ========================================================================
    // Sections & steps
    // ========================================================================

    /// Print a section header (e.g. "  Processing sample.bam").
    pub fn section(&self, text: &str) {
        if self.is_quiet() {
            return;
        }
        eprintln!("  {}", self.style_header.apply_to(text));
    }

    /// Print an indented step line (visible in Normal + Verbose).
    pub fn step(&self, text: &str) {
        if self.is_quiet() {
            return;
        }
        eprintln!("    {text}");
    }

    /// Print an indented detail line (only visible in Verbose).
    pub fn detail(&self, text: &str) {
        if self.verbosity != Verbosity::Verbose {
            return;
        }
        eprintln!("    {}", self.style_dim.apply_to(text));
    }

    // ========================================================================
    // Progress bar
    // ========================================================================

    /// Create a progress bar for BAM read counting.
    ///
    /// Returns a spinner-style bar (total is unknown upfront) that shows the
    /// read count and elapsed time. In Quiet mode, returns a hidden bar.
    pub fn progress_bar(&self) -> ProgressBar {
        if self.is_quiet() {
            return ProgressBar::hidden();
        }
        let pb = ProgressBar::new_spinner();
        pb.set_style(
            ProgressStyle::with_template("    {spinner:.cyan} {msg}  {elapsed:.dim}")
                .expect("valid template")
                .tick_strings(&[
                    "\u{28fe}", "\u{28fd}", "\u{28fb}", "\u{28f7}", "\u{28ef}", "\u{28df}",
                    "\u{28bf}", "\u{287f}", "\u{28fe}",
                ]),
        );
        pb.enable_steady_tick(Duration::from_millis(100));
        pb
    }

    /// Finish the progress bar with a final styled line.
    pub fn finish_progress(&self, pb: &ProgressBar, reads: u64, duration: Duration) {
        if self.is_quiet() {
            pb.finish_and_clear();
            return;
        }
        pb.finish_and_clear();
        eprintln!(
            "    {} {} reads processed  {}",
            self.style_green.apply_to("\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}\u{2501}"),
            self.style_bold.apply_to(format_count(reads)),
            self.style_dim.apply_to(format_duration(duration)),
        );
    }

    // ========================================================================
    // Summary box
    // ========================================================================

    /// Print the per-BAM summary box with key metrics.
    ///
    /// Takes a title and a list of (label, bold_value, annotation) rows.
    #[allow(clippy::too_many_arguments)]
    pub fn summary_box(&self, title: &str, rows: &[(&str, String, String)]) {
        if self.is_quiet() {
            return;
        }
        // Adapt box width to the terminal, clamped to [40, 80].
        let inner_width = console::Term::stderr()
            .size()
            .1
            .saturating_sub(6) // 2 indent + 2 border chars + 2 padding
            as usize;
        let inner_width = inner_width.clamp(40, 80);
        let border_h = "\u{2500}".repeat(inner_width + 2);

        eprintln!();
        eprintln!(
            "  {}",
            self.style_box
                .apply_to(format!("\u{250c}{border_h}\u{2510}"))
        );
        // Title row
        let title_padded = format!(" {title:<width$} ", width = inner_width);
        eprintln!(
            "  {}{}{} ",
            self.style_box.apply_to("\u{2502}"),
            self.style_bold.apply_to(title_padded),
            self.style_box.apply_to("\u{2502}"),
        );
        // Blank separator
        eprintln!(
            "  {}",
            self.style_box.apply_to(format!(
                "\u{2502}{:width$}\u{2502}",
                "",
                width = inner_width + 2
            ))
        );
        // Data rows
        for (label, value, annotation) in rows {
            eprintln!(
                "  {}{}{}",
                self.style_box.apply_to("\u{2502}"),
                format_summary_row(
                    &self.style_label,
                    &self.style_bold,
                    label,
                    value,
                    annotation,
                    inner_width,
                ),
                self.style_box.apply_to("\u{2502}"),
            );
        }
        // Bottom border
        eprintln!(
            "  {}",
            self.style_box
                .apply_to(format!("\u{2514}{border_h}\u{2518}"))
        );
        eprintln!();
    }

    // ========================================================================
    // Output checklist
    // ========================================================================

    /// Print a tool output line: "  ✓ tool_name    path".
    pub fn output_item(&self, tool: &str, path: &str) {
        if self.is_quiet() {
            return;
        }
        eprintln!(
            "    {} {:<22}{}",
            self.style_green.apply_to("\u{2713}"),
            tool,
            self.style_dim.apply_to(path),
        );
    }

    /// Print a tool output sub-detail (verbose only).
    pub fn output_detail(&self, text: &str) {
        if self.verbosity != Verbosity::Verbose {
            return;
        }
        eprintln!("      {}", self.style_dim.apply_to(text));
    }

    /// Print a group header in the output checklist (e.g. "samtools", "RSeQC").
    pub fn output_group(&self, name: &str) {
        if self.is_quiet() {
            return;
        }
        eprintln!("    {}", self.style_label.apply_to(format!("{name}:")));
    }

    // ========================================================================
    // Multi-BAM summary
    // ========================================================================

    /// Print a per-BAM result in the multi-BAM summary.
    pub fn bam_result_ok(&self, name: &str, duration: Duration) {
        if self.is_quiet() {
            return;
        }
        eprintln!(
            "    {} {:<40}{}",
            self.style_green.apply_to("\u{2713}"),
            name,
            self.style_dim.apply_to(format_duration(duration)),
        );
    }

    /// Print a failed BAM result in the multi-BAM summary.
    pub fn bam_result_err(&self, name: &str, error: &str) {
        if self.is_quiet() {
            return;
        }
        eprintln!(
            "    {} {:<40}{}",
            self.style_red.apply_to("\u{2717}"),
            name,
            self.style_red.apply_to(format!("failed: {error}")),
        );
    }

    // ========================================================================
    // Finish
    // ========================================================================

    /// Print a completion line with a label (e.g. BAM stem or "RustQC run finished").
    pub fn finish(&self, label: &str, duration: Duration) {
        if self.is_quiet() {
            return;
        }
        eprintln!(
            "  {} {} {}",
            self.style_green.apply_to("\u{2713}"),
            self.style_header.apply_to(label),
            self.style_dim
                .apply_to(format!("finished in {}", format_duration(duration))),
        );
        eprintln!();
    }

    // ========================================================================
    // Warnings & errors
    // ========================================================================

    /// Print a styled warning (always visible, even in quiet mode).
    pub fn warn(&self, msg: &str) {
        eprintln!(
            "  {} {}",
            self.style_yellow.apply_to("\u{26a0}"),
            self.style_yellow.apply_to(msg),
        );
    }

    /// Print a multi-line warning inside a yellow box with Unicode border.
    ///
    /// The box is drawn with rounded corners and a "Warning" label in the top
    /// border. The first line is rendered **bold** as a heading. Remaining lines
    /// are padded to the box width. Always visible, even in quiet mode.
    pub fn warn_box(&self, lines: &[&str]) {
        // Find the widest line to size the box (minimum 40 chars)
        let content_width = lines.iter().map(|l| l.len()).max().unwrap_or(0).max(40);
        let y = &self.style_yellow;
        let yb = Style::new().yellow().bold();

        // Top border:  ╭─ Warning ─────╮
        let top_label = " Warning ";
        let remaining = content_width + 1 - top_label.len(); // +1 for inner padding
        eprintln!(
            "  {}{}{}",
            y.apply_to("╭─"),
            y.apply_to(top_label),
            y.apply_to(format!("{}╮", "─".repeat(remaining))),
        );
        // Content lines: │ text │
        // First line is bold, rest are normal yellow
        for (i, line) in lines.iter().enumerate() {
            let styled = if i == 0 {
                format!("{}", yb.apply_to(line))
            } else {
                format!("{}", y.apply_to(line))
            };
            // Pad with spaces to fill the box width (use raw len for padding calc)
            let pad = content_width.saturating_sub(line.len());
            eprintln!(
                "  {} {}{} {}",
                y.apply_to("│"),
                styled,
                " ".repeat(pad),
                y.apply_to("│"),
            );
        }
        // Bottom border: ╰──────────╯
        eprintln!(
            "  {}",
            y.apply_to(format!("╰{}╯", "─".repeat(content_width + 2))),
        );
    }

    /// Print a styled error (always visible, even in quiet mode).
    pub fn error(&self, msg: &str) {
        eprintln!(
            "  {} {}",
            self.style_red.apply_to("\u{2717} error:"),
            self.style_red.apply_to(msg),
        );
    }
}

// ============================================================================
// Formatting helpers
// ============================================================================

/// Format a summary box row with styled label, value, and annotation.
///
/// Returns the full inner content string (already padded to `inner_width + 2`).
fn format_summary_row(
    style_label: &Style,
    style_bold: &Style,
    label: &str,
    value: &str,
    annotation: &str,
    inner_width: usize,
) -> String {
    let styled_label = style_label.apply_to(format!("{label:<16}"));
    let styled_value = style_bold.apply_to(format!("{value:>6}"));
    let content = if annotation.is_empty() {
        format!("   {styled_label}{styled_value}")
    } else {
        format!("   {styled_label}{styled_value}  {annotation}")
    };
    // We need to pad based on visible length, not byte length (ANSI codes add bytes).
    // Use console::measure_text_width for accurate visible width.
    let visible_width = console::measure_text_width(&content);
    let target = inner_width + 2;
    if visible_width < target {
        format!("{content}{:width$}", "", width = target - visible_width)
    } else {
        content
    }
}

/// Format a count with SI prefix (e.g. 48200000 → "48.2M").
///
/// Values below 1000 are shown as-is. Values above use K/M/G/T suffixes
/// with one decimal place.
pub fn format_count(n: u64) -> String {
    use number_prefix::NumberPrefix;
    match NumberPrefix::decimal(n as f64) {
        NumberPrefix::Standalone(n) => format!("{n}"),
        NumberPrefix::Prefixed(prefix, n) => {
            // Map SI prefixes to short single-char suffixes
            let suffix = match prefix {
                number_prefix::Prefix::Kilo => "K",
                number_prefix::Prefix::Mega => "M",
                number_prefix::Prefix::Giga => "G",
                number_prefix::Prefix::Tera => "T",
                _ => return format!("{:.1}{prefix:?}", n),
            };
            format!("{n:.1}{suffix}")
        }
    }
}

/// Format a percentage string (e.g. "83.3%").
pub fn format_pct(n: u64, total: u64) -> String {
    if total == 0 {
        return "(0.0%)".to_string();
    }
    format!("({:.1}%)", n as f64 / total as f64 * 100.0)
}

/// Format a duration as human-friendly mm:ss or h:mm:ss.
///
/// - Under 60s: "45.2s"
/// - Under 1h: "1:23"
/// - Over 1h: "1:02:34"
pub fn format_duration(d: Duration) -> String {
    let total_secs = d.as_secs_f64();
    if total_secs < 60.0 {
        return format!("{total_secs:.1}s");
    }
    let total_secs = d.as_secs();
    let hours = total_secs / 3600;
    let minutes = (total_secs % 3600) / 60;
    let seconds = total_secs % 60;
    if hours > 0 {
        format!("{hours}:{minutes:02}:{seconds:02}")
    } else {
        format!("{minutes}:{seconds:02}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_count_small() {
        assert_eq!(format_count(0), "0");
        assert_eq!(format_count(42), "42");
        assert_eq!(format_count(999), "999");
    }

    #[test]
    fn test_format_count_thousands() {
        assert_eq!(format_count(1000), "1.0K");
        assert_eq!(format_count(1500), "1.5K");
        assert_eq!(format_count(50000), "50.0K");
    }

    #[test]
    fn test_format_count_millions() {
        assert_eq!(format_count(1_000_000), "1.0M");
        assert_eq!(format_count(48_200_000), "48.2M");
        assert_eq!(format_count(50_000_000), "50.0M");
    }

    #[test]
    fn test_format_count_billions() {
        assert_eq!(format_count(1_000_000_000), "1.0G");
        assert_eq!(format_count(5_000_000_000), "5.0G");
    }

    #[test]
    fn test_format_pct() {
        assert_eq!(format_pct(833, 1000), "(83.3%)");
        assert_eq!(format_pct(0, 0), "(0.0%)");
        assert_eq!(format_pct(1000, 1000), "(100.0%)");
    }

    #[test]
    fn test_format_duration_seconds() {
        assert_eq!(format_duration(Duration::from_secs_f64(0.5)), "0.5s");
        assert_eq!(format_duration(Duration::from_secs_f64(45.2)), "45.2s");
        assert_eq!(format_duration(Duration::from_secs_f64(59.9)), "59.9s");
    }

    #[test]
    fn test_format_duration_minutes() {
        assert_eq!(format_duration(Duration::from_secs(60)), "1:00");
        assert_eq!(format_duration(Duration::from_secs(83)), "1:23");
        assert_eq!(format_duration(Duration::from_secs(3599)), "59:59");
    }

    #[test]
    fn test_format_duration_hours() {
        assert_eq!(format_duration(Duration::from_secs(3600)), "1:00:00");
        assert_eq!(format_duration(Duration::from_secs(3754)), "1:02:34");
    }
}
