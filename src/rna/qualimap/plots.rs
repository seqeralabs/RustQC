//! Qualimap RNA-Seq QC plot generation.
//!
//! Produces JFreeChart-style plots matching Qualimap's visual output:
//! coverage profile line charts, coverage histogram, and pie charts.

use anyhow::{Context, Result};
use log::info;
use plotters::prelude::*;
use plotters_svg::SVGBackend;
use std::path::Path;

// ============================= Constants =======================================

/// Scale factor for PNG output.
/// Qualimap produces 1024×768 images — we match that exactly.
const SCALE: u32 = 1;

/// Helper to scale a dimension by SCALE.
const fn s(v: u32) -> u32 {
    v * SCALE
}

/// Base chart width in pixels (matches Qualimap's GRAPHIC_TO_SAVE_WIDTH).
const WIDTH: u32 = 1024;

/// Base chart height in pixels (matches Qualimap's GRAPHIC_TO_SAVE_HEIGHT).
const HEIGHT: u32 = 768;

// ============================= JFreeChart Colors ===============================

/// Chart background color: rgb(230, 230, 230) — JFreeChart default.
const CHART_BG: RGBColor = RGBColor(230, 230, 230);

/// Plot area background: white.
const PLOT_BG: RGBColor = RGBColor(255, 255, 255);

/// Gridline color: warm tan rgb(192, 168, 144) — matches Qualimap reference output.
const GRIDLINE_COLOR: RGBColor = RGBColor(192, 168, 144);

/// Title/subtitle text color: dark gray rgb(64, 64, 64).
const TITLE_COLOR: RGBColor = RGBColor(64, 64, 64);

/// Coverage profile line color: pure red.
const COVERAGE_LINE_COLOR: RGBColor = RGBColor(255, 0, 0);

/// Histogram bar border: dark gray rgb(50, 50, 50).
const HISTOGRAM_BAR_BORDER: RGBColor = RGBColor(50, 50, 50);

/// Qualimap pie palette — sampled from actual Qualimap output PNGs.
const PIE_COLOR_0: RGBColor = RGBColor(0xEC, 0x62, 0x5C); // #EC625C — red (Exonic / Known)
const PIE_COLOR_1: RGBColor = RGBColor(0x55, 0x55, 0xF6); // #5555F6 — blue (Intronic / Partly known)
const PIE_COLOR_2: RGBColor = RGBColor(0x8A, 0xFC, 0x6E); // #8AFC6E — green (Intergenic / Novel)

/// Pie chart shadow color: gray rgb(128, 128, 128).
const PIE_SHADOW_COLOR: RGBColor = RGBColor(0x80, 0x80, 0x80);

/// Pie label background: light yellow.
const PIE_LABEL_BG: RGBColor = RGBColor(255, 255, 204);

/// Pie label border: dark gray.
const PIE_LABEL_BORDER: RGBColor = RGBColor(153, 153, 153);

// ============================= Coverage Profile ================================

/// Generate the "Coverage Profile Along Genes" line chart.
///
/// Produces both PNG and SVG files. The chart shows normalized coverage
/// across 100 transcript position bins (5' to 3').
///
/// # Arguments
/// * `profile` - 100-element array of coverage values (from `compute_coverage_profile`)
/// * `title` - Chart title (e.g., "Coverage Profile Along Genes (Total)")
/// * `sample_name` - BAM file name shown as subtitle
/// * `output_path` - Path for the PNG file (SVG uses same stem)
pub fn coverage_profile_plot(
    profile: &[f64; 100],
    title: &str,
    sample_name: &str,
    output_path: &Path,
) -> Result<()> {
    // PNG output
    {
        let root = BitMapBackend::new(output_path, (s(WIDTH), s(HEIGHT))).into_drawing_area();
        render_coverage_profile(&root, profile, title, sample_name, SCALE as f64)?;
        root.present()
            .context("Failed to write coverage profile PNG")?;
    }

    // SVG output
    let svg_path = output_path.with_extension("svg");
    {
        let root = SVGBackend::new(&svg_path, (WIDTH, HEIGHT)).into_drawing_area();
        render_coverage_profile(&root, profile, title, sample_name, 1.0)?;
        root.present()
            .context("Failed to write coverage profile SVG")?;
    }

    info!("Wrote coverage profile: {}", output_path.display());
    Ok(())
}

/// Render the coverage profile line chart on any drawing backend.
fn render_coverage_profile<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    profile: &[f64; 100],
    title: &str,
    sample_name: &str,
    pxs: f64,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    let ps = |v: f64| (v * pxs) as u32;

    // Fill chart background
    root.fill(&CHART_BG).map_err(|e| anyhow::anyhow!("{}", e))?;

    // Compute Y-axis range
    let y_max = profile.iter().cloned().fold(0.0f64, f64::max);
    // Round up to a nice number for the Y-axis ceiling
    let y_ceil = nice_ceil(y_max);

    // Draw title and subtitle — sizes scaled up to match Qualimap visual appearance.
    // Plotters' sans-serif renders smaller than Java's SansSerif at the same pt size.
    let title_font = ("sans-serif", ps(24.0) as f64)
        .into_font()
        .style(FontStyle::Bold)
        .color(&TITLE_COLOR);
    let subtitle_font = ("sans-serif", ps(16.0) as f64)
        .into_font()
        .color(&TITLE_COLOR);

    let title_y = ps(10.0) as i32;
    let subtitle_y = title_y + ps(28.0) as i32;

    // Center title
    let (w, _h) = root.dim_in_pixel();
    root.draw_text(
        title,
        &title_font.into_text_style(root),
        (
            (w as i32 - estimate_text_width(title, ps(24.0) as f64) as i32) / 2,
            title_y,
        ),
    )
    .map_err(|e| anyhow::anyhow!("{}", e))?;

    // Center subtitle
    root.draw_text(
        sample_name,
        &subtitle_font.into_text_style(root),
        (
            (w as i32 - estimate_text_width(sample_name, ps(16.0) as f64) as i32) / 2,
            subtitle_y,
        ),
    )
    .map_err(|e| anyhow::anyhow!("{}", e))?;

    // Create chart area below title region
    let chart_top = ps(55.0);
    let chart_area = root.clone().margin(chart_top, ps(20.0), ps(20.0), ps(20.0));

    let mut chart = ChartBuilder::on(&chart_area)
        .x_label_area_size(ps(40.0))
        .y_label_area_size(ps(70.0))
        .build_cartesian_2d(0.0..100.0f64, 0.0..y_ceil)
        .map_err(|e| anyhow::anyhow!("{}", e))?;

    // White plot background
    chart
        .plotting_area()
        .fill(&PLOT_BG)
        .map_err(|e| anyhow::anyhow!("{}", e))?;

    // Configure mesh to match Qualimap gridline style.
    // Qualimap has gridlines every 5 units on X (0, 5, 10, ..., 100).
    chart
        .configure_mesh()
        .x_desc("Transcript position")
        .y_desc("Counts")
        .x_labels(21) // 0, 5, 10, ..., 100 = 21 labels
        .label_style(("sans-serif", ps(15.0) as f64).into_font().color(&BLACK))
        .axis_desc_style(("sans-serif", ps(16.0) as f64).into_font().color(&BLACK))
        .x_label_formatter(&|v| format!("{}", *v as i32))
        .y_label_formatter(&|v| format_y_label(*v))
        .axis_style(BLACK.stroke_width(1))
        .light_line_style(GRIDLINE_COLOR.mix(0.0)) // hide intermediate gridlines
        .bold_line_style(GRIDLINE_COLOR.stroke_width(1))
        .set_all_tick_mark_size(ps(4.0))
        .draw()
        .map_err(|e| anyhow::anyhow!("{}", e))?;

    // Draw the coverage line
    let data: Vec<(f64, f64)> = profile
        .iter()
        .enumerate()
        .map(|(i, &v)| (i as f64, v))
        .collect();

    chart
        .draw_series(LineSeries::new(
            data,
            COVERAGE_LINE_COLOR.stroke_width(ps(2.0)),
        ))
        .map_err(|e| anyhow::anyhow!("{}", e))?;

    Ok(())
}

// ============================= Coverage Histogram ===============================

/// Generate the "Coverage Histogram (0-50X)" bar chart.
///
/// # Arguments
/// * `histogram` - 51-element array where `histogram[i]` = number of transcripts
///   with mean coverage level `i` (bin 50 includes all transcripts with mean >= 50)
/// * `sample_name` - BAM file name shown as subtitle
/// * `output_path` - Path for the PNG file
pub fn coverage_histogram_plot(
    histogram: &[u64; 51],
    sample_name: &str,
    output_path: &Path,
) -> Result<()> {
    // PNG output
    {
        let root = BitMapBackend::new(output_path, (s(WIDTH), s(HEIGHT))).into_drawing_area();
        render_coverage_histogram(&root, histogram, sample_name, SCALE as f64)?;
        root.present()
            .context("Failed to write coverage histogram PNG")?;
    }

    // SVG output
    let svg_path = output_path.with_extension("svg");
    {
        let root = SVGBackend::new(&svg_path, (WIDTH, HEIGHT)).into_drawing_area();
        render_coverage_histogram(&root, histogram, sample_name, 1.0)?;
        root.present()
            .context("Failed to write coverage histogram SVG")?;
    }

    info!("Wrote coverage histogram: {}", output_path.display());
    Ok(())
}

/// Render the coverage histogram bar chart.
fn render_coverage_histogram<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    histogram: &[u64; 51],
    sample_name: &str,
    pxs: f64,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    let ps = |v: f64| (v * pxs) as u32;

    root.fill(&CHART_BG).map_err(|e| anyhow::anyhow!("{}", e))?;

    let y_max = *histogram.iter().max().unwrap_or(&1) as f64;
    let y_ceil = nice_ceil(y_max);

    let title = "Coverage Histogram (0-50X)";
    let title_font = ("sans-serif", ps(24.0) as f64)
        .into_font()
        .style(FontStyle::Bold)
        .color(&TITLE_COLOR);
    let subtitle_font = ("sans-serif", ps(16.0) as f64)
        .into_font()
        .color(&TITLE_COLOR);

    let (w, _h) = root.dim_in_pixel();
    let title_y = ps(10.0) as i32;
    let subtitle_y = title_y + ps(28.0) as i32;

    root.draw_text(
        title,
        &title_font.into_text_style(root),
        (
            (w as i32 - estimate_text_width(title, ps(24.0) as f64) as i32) / 2,
            title_y,
        ),
    )
    .map_err(|e| anyhow::anyhow!("{}", e))?;

    root.draw_text(
        sample_name,
        &subtitle_font.into_text_style(root),
        (
            (w as i32 - estimate_text_width(sample_name, ps(16.0) as f64) as i32) / 2,
            subtitle_y,
        ),
    )
    .map_err(|e| anyhow::anyhow!("{}", e))?;

    let chart_area = root.clone().margin(ps(55.0), ps(20.0), ps(20.0), ps(20.0));

    let mut chart = ChartBuilder::on(&chart_area)
        .x_label_area_size(ps(40.0))
        .y_label_area_size(ps(70.0))
        .build_cartesian_2d((-0.5f64)..50.5f64, 0.0..y_ceil)
        .map_err(|e| anyhow::anyhow!("{}", e))?;

    // White plot area + gridlines
    chart
        .plotting_area()
        .fill(&PLOT_BG)
        .map_err(|e| anyhow::anyhow!("{}", e))?;

    // Qualimap shows gridlines every 5 units on X-axis (0, 5, 10, ..., 50)
    chart
        .configure_mesh()
        .x_desc("Coverage (X)")
        .y_desc("Number of transcripts")
        .x_labels(11) // 0, 5, 10, ..., 50
        .label_style(("sans-serif", ps(15.0) as f64).into_font().color(&BLACK))
        .axis_desc_style(("sans-serif", ps(16.0) as f64).into_font().color(&BLACK))
        .x_label_formatter(&|v| {
            let iv = v.round() as i32;
            if (0..=50).contains(&iv) && (v - iv as f64).abs() < 0.01 {
                format!("{iv}")
            } else {
                String::new()
            }
        })
        .y_label_formatter(&|v| format_y_label(*v))
        .axis_style(BLACK.stroke_width(1))
        .light_line_style(GRIDLINE_COLOR.mix(0.0)) // hide intermediate gridlines
        .bold_line_style(GRIDLINE_COLOR.stroke_width(1))
        .set_all_tick_mark_size(ps(4.0))
        .draw()
        .map_err(|e| anyhow::anyhow!("{}", e))?;

    // Draw bars
    let bar_width = 0.45; // Half-width of each bar
    chart
        .draw_series(histogram.iter().enumerate().map(|(i, &count)| {
            let x0 = i as f64 - bar_width;
            let x1 = i as f64 + bar_width;
            Rectangle::new(
                [(x0, 0.0), (x1, count as f64)],
                ShapeStyle {
                    color: RGBColor(100, 100, 250).mix(0.59),
                    filled: true,
                    stroke_width: ps(1.0),
                },
            )
        }))
        .map_err(|e| anyhow::anyhow!("{}", e))?;

    // Draw bar borders on top
    chart
        .draw_series(histogram.iter().enumerate().map(|(i, &count)| {
            let x0 = i as f64 - bar_width;
            let x1 = i as f64 + bar_width;
            Rectangle::new(
                [(x0, 0.0), (x1, count as f64)],
                ShapeStyle {
                    color: HISTOGRAM_BAR_BORDER.to_rgba(),
                    filled: false,
                    stroke_width: ps(1.0),
                },
            )
        }))
        .map_err(|e| anyhow::anyhow!("{}", e))?;

    Ok(())
}

// ============================= Pie Charts ======================================

/// Data for a single pie chart slice.
pub struct PieSlice<'a> {
    /// Label text (e.g., "Exonic").
    pub label: &'a str,
    /// Value (will be converted to percentage internally).
    pub value: f64,
    /// Slice color.
    pub color: RGBColor,
}

/// Generate the "Reads Genomic Origin" pie chart.
///
/// # Arguments
/// * `exonic` - Number of exonic reads
/// * `intronic` - Number of intronic reads
/// * `intergenic` - Number of intergenic reads
/// * `sample_name` - BAM file name shown as subtitle
/// * `output_path` - Path for the PNG file
pub fn reads_genomic_origin_plot(
    exonic: u64,
    intronic: u64,
    intergenic: u64,
    sample_name: &str,
    output_path: &Path,
) -> Result<()> {
    let total = exonic + intronic + intergenic;
    if total == 0 {
        return Ok(());
    }

    let slices = vec![
        PieSlice {
            label: "Exonic",
            value: exonic as f64,
            color: PIE_COLOR_0,
        },
        PieSlice {
            label: "Intronic",
            value: intronic as f64,
            color: PIE_COLOR_1,
        },
        PieSlice {
            label: "Intergenic",
            value: intergenic as f64,
            color: PIE_COLOR_2,
        },
    ];

    pie_chart("Reads Genomic Origin", sample_name, &slices, output_path)
}

/// Generate the "Junction Analysis" pie chart.
///
/// # Arguments
/// * `known` - Number of known junction events
/// * `partly_known` - Number of partly known junction events
/// * `novel` - Number of novel junction events
/// * `sample_name` - BAM file name shown as subtitle
/// * `output_path` - Path for the PNG file
pub fn junction_analysis_plot(
    known: u64,
    partly_known: u64,
    novel: u64,
    sample_name: &str,
    output_path: &Path,
) -> Result<()> {
    let total = known + partly_known + novel;
    if total == 0 {
        return Ok(());
    }

    let slices = vec![
        PieSlice {
            label: "Known",
            value: known as f64,
            color: PIE_COLOR_0,
        },
        PieSlice {
            label: "Partly known",
            value: partly_known as f64,
            color: PIE_COLOR_1,
        },
        PieSlice {
            label: "Novel",
            value: novel as f64,
            color: PIE_COLOR_2,
        },
    ];

    pie_chart("Junction Analysis", sample_name, &slices, output_path)
}

/// Generate a Qualimap-style pie chart (dual PNG + SVG output).
fn pie_chart(
    title: &str,
    sample_name: &str,
    slices: &[PieSlice<'_>],
    output_path: &Path,
) -> Result<()> {
    // PNG output
    {
        let root = BitMapBackend::new(output_path, (s(WIDTH), s(HEIGHT))).into_drawing_area();
        render_pie_chart(&root, title, sample_name, slices, SCALE as f64)?;
        root.present().context("Failed to write pie chart PNG")?;
    }

    // SVG output
    let svg_path = output_path.with_extension("svg");
    {
        let root = SVGBackend::new(&svg_path, (WIDTH, HEIGHT)).into_drawing_area();
        render_pie_chart(&root, title, sample_name, slices, 1.0)?;
        root.present().context("Failed to write pie chart SVG")?;
    }

    info!("Wrote pie chart: {}", output_path.display());
    Ok(())
}

/// Render a Qualimap-style pie chart with JFreeChart visual styling.
///
/// Features: drop shadow, label boxes with leader lines, JFreeChart default palette.
fn render_pie_chart<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    title: &str,
    sample_name: &str,
    slices: &[PieSlice<'_>],
    pxs: f64,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    let ps = |v: f64| (v * pxs) as i32;
    let psu = |v: f64| (v * pxs) as u32;

    root.fill(&CHART_BG).map_err(|e| anyhow::anyhow!("{}", e))?;

    let (w, h) = root.dim_in_pixel();
    let w = w as i32;
    let h = h as i32;

    // Draw title and subtitle — font sizes matched to Qualimap reference
    let title_font = ("sans-serif", ps(24.0) as f64)
        .into_font()
        .style(FontStyle::Bold)
        .color(&TITLE_COLOR);
    let subtitle_font = ("sans-serif", ps(16.0) as f64)
        .into_font()
        .color(&TITLE_COLOR);

    let title_y = ps(12.0);
    let subtitle_y = title_y + ps(28.0);

    root.draw_text(
        title,
        &title_font.into_text_style(root),
        (
            (w - estimate_text_width(title, ps(24.0) as f64) as i32) / 2,
            title_y,
        ),
    )
    .map_err(|e| anyhow::anyhow!("{}", e))?;

    root.draw_text(
        sample_name,
        &subtitle_font.into_text_style(root),
        (
            (w - estimate_text_width(sample_name, ps(16.0) as f64) as i32) / 2,
            subtitle_y,
        ),
    )
    .map_err(|e| anyhow::anyhow!("{}", e))?;

    // Qualimap draws a plot frame: E6E6E6 filled rectangle with 1px black inner border,
    // then the white pie area sits inside it.
    let frame_margin = ps(8.0);
    let frame_top = subtitle_y + ps(12.0);
    let frame_bottom = h - frame_margin;
    let frame_left = frame_margin;
    let frame_right = w - frame_margin;

    // E6E6E6 plot frame background
    root.draw(&Rectangle::new(
        [(frame_left, frame_top), (frame_right, frame_bottom)],
        ShapeStyle {
            color: CHART_BG.to_rgba(),
            filled: true,
            stroke_width: 0,
        },
    ))
    .map_err(|e| anyhow::anyhow!("{}", e))?;

    // White interior
    let inner_pad = ps(1.0);
    root.draw(&Rectangle::new(
        [
            (frame_left + inner_pad, frame_top + inner_pad),
            (frame_right - inner_pad, frame_bottom - inner_pad),
        ],
        ShapeStyle {
            color: PLOT_BG.to_rgba(),
            filled: true,
            stroke_width: 0,
        },
    ))
    .map_err(|e| anyhow::anyhow!("{}", e))?;

    // 1px black border around frame
    root.draw(&Rectangle::new(
        [(frame_left, frame_top), (frame_right, frame_bottom)],
        ShapeStyle {
            color: BLACK.to_rgba(),
            filled: false,
            stroke_width: 1,
        },
    ))
    .map_err(|e| anyhow::anyhow!("{}", e))?;

    // Compute pie geometry within the white interior
    let plot_margin = ps(30.0);
    let plot_top = frame_top + ps(10.0);
    let plot_w = frame_right - frame_left - 2 * plot_margin;
    let plot_h = frame_bottom - plot_top - plot_margin;
    let cx = frame_left + (frame_right - frame_left) / 2;
    let cy = plot_top + plot_h / 2;
    let radius = (plot_w.min(plot_h) as f64 * 0.38) as i32;

    // Compute total and percentages
    let total: f64 = slices.iter().map(|s| s.value).sum();
    if total <= 0.0 {
        return Ok(());
    }

    // Draw drop shadow (offset filled circle, subtle gray)
    let shadow_offset = ps(4.0);
    let n_shadow_pts = 360;
    let shadow_points: Vec<(i32, i32)> = (0..=n_shadow_pts)
        .map(|i| {
            let angle = 2.0 * std::f64::consts::PI * i as f64 / n_shadow_pts as f64;
            (
                cx + shadow_offset + (radius as f64 * angle.cos()) as i32,
                cy + shadow_offset + (radius as f64 * angle.sin()) as i32,
            )
        })
        .collect();
    root.draw(&Polygon::new(
        shadow_points,
        ShapeStyle {
            color: PIE_SHADOW_COLOR.mix(0.5),
            filled: true,
            stroke_width: 0,
        },
    ))
    .map_err(|e| anyhow::anyhow!("{}", e))?;

    // Draw pie slices
    // JFreeChart starts at 12 o'clock (90 degrees) and goes clockwise
    let mut start_angle = std::f64::consts::FRAC_PI_2;

    // Collect label info for drawing after all slices
    let mut label_info: Vec<(f64, String, RGBColor)> = Vec::new();

    for slice in slices {
        let pct = slice.value / total * 100.0;
        if pct <= 0.0 {
            continue;
        }
        let sweep = slice.value / total * 2.0 * std::f64::consts::PI;

        // Build polygon: center -> arc points -> center
        let n_segments = ((sweep / 0.02) as usize).max(4);
        let mut points = Vec::with_capacity(n_segments + 3);
        points.push((cx, cy));
        for j in 0..=n_segments {
            let angle = start_angle - j as f64 * sweep / n_segments as f64;
            points.push((
                cx + (radius as f64 * angle.cos()) as i32,
                cy - (radius as f64 * angle.sin()) as i32,
            ));
        }
        points.push((cx, cy));

        // Filled slice — fully opaque
        root.draw(&Polygon::new(points, slice.color.filled()))
            .map_err(|e| anyhow::anyhow!("{}", e))?;

        // Store label info
        let mid_angle = start_angle - sweep / 2.0;
        let label_text = format!(
            "{} - {}%",
            slice.label,
            trim_trailing_zeros(&format!("{:.2}", pct))
        );
        label_info.push((mid_angle, label_text, slice.color));

        start_angle -= sweep;
    }

    // Draw labels with leader lines and background boxes
    for (mid_angle, label_text, _color) in &label_info {
        // Leader line start (on pie edge)
        let edge_x = cx + (radius as f64 * mid_angle.cos()) as i32;
        let edge_y = cy - (radius as f64 * mid_angle.sin()) as i32;

        // Leader line end (outside pie)
        let label_radius = radius as f64 * 1.15;
        let label_x = cx + (label_radius * mid_angle.cos()) as i32;
        let label_y = cy - (label_radius * mid_angle.sin()) as i32;

        // Draw leader line
        root.draw(&PathElement::new(
            vec![(edge_x, edge_y), (label_x, label_y)],
            ShapeStyle {
                color: RGBColor(80, 80, 80).to_rgba(),
                filled: false,
                stroke_width: psu(1.0),
            },
        ))
        .map_err(|e| anyhow::anyhow!("{}", e))?;

        // Label box dimensions — tighter padding to match JFreeChart
        let font_size = ps(14.0) as f64;
        let text_w = estimate_text_width(label_text, font_size) as i32;
        let text_h = ps(14.0);
        let pad_x = ps(3.0);
        let pad_y = ps(2.0);

        // Position label box: extend outward from pie center
        let (box_x, box_y) = if mid_angle.cos() >= 0.0 {
            // Right side: label extends to the right
            (label_x + ps(3.0), label_y - text_h / 2 - pad_y)
        } else {
            // Left side: label extends to the left
            (
                label_x - text_w - 2 * pad_x - ps(3.0),
                label_y - text_h / 2 - pad_y,
            )
        };

        // Draw label background box
        root.draw(&Rectangle::new(
            [
                (box_x, box_y),
                (box_x + text_w + 2 * pad_x, box_y + text_h + 2 * pad_y),
            ],
            ShapeStyle {
                color: PIE_LABEL_BG.to_rgba(),
                filled: true,
                stroke_width: psu(1.0),
            },
        ))
        .map_err(|e| anyhow::anyhow!("{}", e))?;

        // Draw label border
        root.draw(&Rectangle::new(
            [
                (box_x, box_y),
                (box_x + text_w + 2 * pad_x, box_y + text_h + 2 * pad_y),
            ],
            ShapeStyle {
                color: PIE_LABEL_BORDER.to_rgba(),
                filled: false,
                stroke_width: psu(1.0),
            },
        ))
        .map_err(|e| anyhow::anyhow!("{}", e))?;

        // Draw label text
        let text_style = ("sans-serif", font_size)
            .into_font()
            .color(&BLACK)
            .into_text_style(root);

        root.draw_text(label_text, &text_style, (box_x + pad_x, box_y + pad_y))
            .map_err(|e| anyhow::anyhow!("{}", e))?;
    }

    Ok(())
}

// ============================= Helpers =========================================

/// Format Y-axis labels with comma separators (e.g., 675000 → "675,000").
fn format_y_label(val: f64) -> String {
    if val == 0.0 {
        return "0".to_string();
    }

    // For small values (< 1), show decimal
    if val.abs() < 1.0 {
        let s = format!("{:.2}", val);
        return s.trim_end_matches('0').trim_end_matches('.').to_string();
    }

    let n = val.round() as i64;
    let s = n.to_string();
    let bytes = s.as_bytes();
    let len = bytes.len();
    let negative = bytes[0] == b'-';
    let start = if negative { 1 } else { 0 };
    let digit_len = len - start;

    if digit_len <= 3 {
        return s;
    }

    let mut result = String::with_capacity(len + (digit_len - 1) / 3);
    for (i, &b) in bytes.iter().enumerate() {
        if i > start && (len - i).is_multiple_of(3) {
            result.push(',');
        }
        result.push(b as char);
    }
    result
}

/// Compute a "nice" ceiling value for Y-axis range.
///
/// Rounds up to the nearest nice number (multiple of 2500, 5000, 10000, etc.)
/// to match JFreeChart's auto-range behavior.
fn nice_ceil(val: f64) -> f64 {
    if val <= 0.0 {
        return 1.0;
    }

    // For very small values (< 1), use decimal rounding
    if val < 1.0 {
        let magnitude = 10.0f64.powf(val.log10().floor());
        let normalized = val / magnitude;
        let nice = if normalized <= 1.0 {
            1.0
        } else if normalized <= 2.0 {
            2.0
        } else if normalized <= 5.0 {
            5.0
        } else {
            10.0
        };
        return nice * magnitude * 1.05;
    }

    // For larger values, find a nice round number
    let magnitude = 10.0f64.powf(val.log10().floor());
    let normalized = val / magnitude;

    let nice = if normalized <= 1.0 {
        1.0
    } else if normalized <= 1.5 {
        1.5
    } else if normalized <= 2.0 {
        2.0
    } else if normalized <= 2.5 {
        2.5
    } else if normalized <= 3.0 {
        3.0
    } else if normalized <= 4.0 {
        4.0
    } else if normalized <= 5.0 {
        5.0
    } else if normalized <= 7.5 {
        7.5
    } else {
        10.0
    };

    nice * magnitude
}

/// Estimate text width in pixels for a given font size.
///
/// Uses approximate character width as 0.6 * font_size (sans-serif average).
fn estimate_text_width(text: &str, font_size: f64) -> f64 {
    text.len() as f64 * font_size * 0.6
}

/// Trim trailing zeros from a decimal string (but keep at least one decimal digit).
fn trim_trailing_zeros(s: &str) -> String {
    if !s.contains('.') {
        return s.to_string();
    }
    let trimmed = s.trim_end_matches('0');
    if trimmed.ends_with('.') {
        // Keep at least one decimal digit
        format!("{trimmed}0")
    } else {
        trimmed.to_string()
    }
}

// ============================= Tests ===========================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_y_label() {
        assert_eq!(format_y_label(0.0), "0");
        assert_eq!(format_y_label(500.0), "500");
        assert_eq!(format_y_label(1000.0), "1,000");
        assert_eq!(format_y_label(675000.0), "675,000");
        assert_eq!(format_y_label(1234567.0), "1,234,567");
    }

    #[test]
    fn test_nice_ceil() {
        assert_eq!(nice_ceil(675000.0), 750000.0);
        assert_eq!(nice_ceil(87500.0), 100000.0);
        assert_eq!(nice_ceil(210000.0), 250000.0);
        assert_eq!(nice_ceil(0.2), 0.21000000000000002); // small value
    }

    #[test]
    fn test_trim_trailing_zeros() {
        assert_eq!(trim_trailing_zeros("76.00"), "76.0");
        assert_eq!(trim_trailing_zeros("18.38"), "18.38");
        assert_eq!(trim_trailing_zeros("5.20"), "5.2");
        assert_eq!(trim_trailing_zeros("100.00"), "100.0");
    }
}
