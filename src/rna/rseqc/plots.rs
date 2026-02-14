//! Native plot generation for RSeQC analyses.
//!
//! Generates PNG + SVG plots directly using the `plotters` crate, replacing the
//! R script intermediate step. Visual output matches the original RSeQC Python
//! tool plots as closely as possible.

use anyhow::{Context, Result};
use log::info;
use plotters::prelude::*;
use plotters_svg::SVGBackend;
use std::path::Path;

use super::inner_distance::InnerDistanceResult;
use super::junction_annotation::JunctionResults;
use super::junction_saturation::SaturationResult;
use super::read_duplication::ReadDuplicationResult;

// ============================================================================
// Constants
// ============================================================================

/// Scale factor for high-resolution PNG output.
const SCALE: u32 = 4;

/// Scale a pixel dimension for PNG output.
const fn s(v: u32) -> u32 {
    v * SCALE
}

/// Default plot dimensions (SVG). PNG is SCALE× larger.
const WIDTH: u32 = 480;
const HEIGHT: u32 = 480;

// ============================================================================
// Inner distance histogram + density plot
// ============================================================================

/// Generate the inner distance histogram with density overlay (PNG + SVG).
///
/// Matches the plot produced by RSeQC's inner_distance.py R script:
/// - Blue-bordered histogram bars (probability density scale)
/// - Red density curve overlay
/// - Title shows mean and standard deviation
///
/// # Arguments
/// * `result` — inner distance analysis results
/// * `step` — bin width used for histogram construction
/// * `output_path` — path for the PNG output; SVG is written alongside
pub fn inner_distance_plot(
    result: &InnerDistanceResult,
    step: i64,
    output_path: &Path,
) -> Result<()> {
    // PNG (high-res)
    {
        let root = BitMapBackend::new(output_path, (s(WIDTH), s(HEIGHT))).into_drawing_area();
        render_inner_distance(&root, result, step, SCALE as f64)?;
        root.present()
            .context("Failed to write inner distance PNG")?;
    }

    // SVG
    let svg_path = output_path.with_extension("svg");
    {
        let root = SVGBackend::new(&svg_path, (WIDTH, HEIGHT)).into_drawing_area();
        render_inner_distance(&root, result, step, 1.0)?;
        root.present()
            .context("Failed to write inner distance SVG")?;
    }

    info!("Wrote inner distance plot: {}", output_path.display());
    Ok(())
}

/// Render the inner distance histogram with density overlay.
fn render_inner_distance<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    result: &InnerDistanceResult,
    step: i64,
    pxs: f64,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    let ps = |v: f64| (v * pxs) as u32;

    root.fill(&WHITE)?;

    // --- Compute histogram data ---
    let bins = &result.histogram;
    if bins.is_empty() {
        return Ok(());
    }

    let total_count: u64 = bins.iter().map(|&(_, _, c)| c).sum();
    if total_count == 0 {
        return Ok(());
    }

    // Convert to probability density (matches R's hist(probability=T))
    let bin_width = step as f64;
    let densities: Vec<f64> = bins
        .iter()
        .map(|&(_, _, c)| c as f64 / (total_count as f64 * bin_width))
        .collect();

    let x_min = bins.first().map(|&(s, _, _)| s).unwrap_or(0) as f64;
    let x_max = bins.last().map(|&(_, e, _)| e).unwrap_or(0) as f64;
    let y_max = densities.iter().copied().fold(0.0_f64, f64::max) * 1.1; // 10% headroom

    // --- Compute density curve (kernel density estimation) ---
    // Replicate R's density(fragsize, bw = 2 * step)
    let bw = 2.0 * step as f64;
    let density_curve = compute_density_curve(bins, total_count, bw, x_min, x_max);

    // Update y_max to accommodate density curve
    let density_y_max = density_curve
        .iter()
        .map(|&(_, y)| y)
        .fold(0.0_f64, f64::max);
    let y_max = y_max.max(density_y_max * 1.1);

    // --- Compute mean and SD for title ---
    let (mean, sd) = compute_weighted_mean_sd(bins, step);

    let title = format!("Mean= {:.4} ; SD= {:.4}", mean, sd);

    // --- Build chart ---
    let mut chart = ChartBuilder::on(root)
        .margin(ps(10.0))
        .x_label_area_size(ps(35.0))
        .y_label_area_size(ps(50.0))
        .caption(&title, ("sans-serif", ps(14.0)))
        .build_cartesian_2d(x_min..x_max, 0.0..y_max)?;

    chart
        .configure_mesh()
        .x_desc("mRNA insert size (bp)")
        .y_desc("Density")
        .label_style(("sans-serif", ps(11.0)))
        .axis_desc_style(("sans-serif", ps(12.0)))
        .x_label_formatter(&|v| format!("{}", *v as i64))
        .draw()?;

    // --- Draw histogram bars ---
    for (i, &(start, end, _count)) in bins.iter().enumerate() {
        let density = densities[i];
        if density > 0.0 {
            // Gray filled bar
            chart.draw_series(std::iter::once(Rectangle::new(
                [(start as f64, 0.0), (end as f64, density)],
                ShapeStyle {
                    color: RGBAColor(211, 211, 211, 1.0),
                    filled: true,
                    stroke_width: 0,
                },
            )))?;
            // Blue border (matches R's border="blue")
            chart.draw_series(std::iter::once(Rectangle::new(
                [(start as f64, 0.0), (end as f64, density)],
                ShapeStyle {
                    color: BLUE.mix(1.0),
                    filled: false,
                    stroke_width: ps(1.0),
                },
            )))?;
        }
    }

    // --- Draw density curve (red) ---
    chart.draw_series(LineSeries::new(
        density_curve.iter().copied(),
        ShapeStyle {
            color: RED.mix(1.0),
            filled: false,
            stroke_width: ps(1.5),
        },
    ))?;

    Ok(())
}

/// Compute a Gaussian kernel density estimate.
///
/// Replicates R's `density(fragsize, bw=...)` where fragsize is the expanded
/// set of bin centers repeated by their counts.
fn compute_density_curve(
    bins: &[(i64, i64, u64)],
    total_count: u64,
    bandwidth: f64,
    x_min: f64,
    x_max: f64,
) -> Vec<(f64, f64)> {
    let n_points = 512; // R default
    let margin = 3.0 * bandwidth;
    let eval_min = x_min - margin;
    let eval_max = x_max + margin;
    let dx = (eval_max - eval_min) / (n_points - 1) as f64;

    let mut curve = Vec::with_capacity(n_points);

    for i in 0..n_points {
        let x = eval_min + i as f64 * dx;
        let mut density = 0.0;

        for &(start, end, count) in bins {
            if count == 0 {
                continue;
            }
            let center = (start as f64 + end as f64) / 2.0;
            let z = (x - center) / bandwidth;
            // Gaussian kernel
            let k = (-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt();
            density += count as f64 * k;
        }

        density /= total_count as f64 * bandwidth;
        curve.push((x, density));
    }

    curve
}

/// Compute weighted mean and standard deviation from histogram bins.
fn compute_weighted_mean_sd(bins: &[(i64, i64, u64)], step: i64) -> (f64, f64) {
    let mut sum = 0.0_f64;
    let mut sum_sq = 0.0_f64;
    let mut n = 0_u64;

    for &(start, _end, count) in bins {
        let center = start as f64 + step as f64 / 2.0;
        sum += center * count as f64;
        sum_sq += center * center * count as f64;
        n += count;
    }

    if n == 0 {
        return (0.0, 0.0);
    }

    let mean = sum / n as f64;
    let variance = sum_sq / n as f64 - mean * mean;
    let sd = if variance > 0.0 { variance.sqrt() } else { 0.0 };

    (mean, sd)
}

// ============================================================================
// Junction annotation pie charts
// ============================================================================

/// Generate junction annotation pie charts (PNG + SVG).
///
/// Matches the plots produced by RSeQC's junction_annotation.py R script:
/// - Two pie charts: splicing events and splicing junctions
/// - Three categories: known (blue), partial_novel (red), complete_novel (green)
/// - Labels show rounded integer percentages
///
/// # Arguments
/// * `results` — junction annotation analysis results
/// * `prefix` — output path prefix (e.g. `outdir/sample`)
pub fn junction_annotation_plot(results: &JunctionResults, prefix: &str) -> Result<()> {
    // --- Splice events pie chart ---
    {
        let total_classified_events =
            results.known_events + results.partial_novel_events + results.complete_novel_events;
        let (e_known_pct, e_partial_pct, e_novel_pct) = if total_classified_events > 0 {
            (
                results.known_events as f64 * 100.0 / total_classified_events as f64,
                results.partial_novel_events as f64 * 100.0 / total_classified_events as f64,
                results.complete_novel_events as f64 * 100.0 / total_classified_events as f64,
            )
        } else {
            (0.0, 0.0, 0.0)
        };

        let events_slices = [
            ("known", e_known_pct, RGBColor(0, 0, 205)),
            ("partial_novel", e_partial_pct, RGBColor(205, 0, 0)),
            ("complete_novel", e_novel_pct, RGBColor(0, 205, 0)),
        ];

        let events_png = format!("{}.splice_events.png", prefix);
        let events_path = Path::new(&events_png);

        // PNG
        {
            let root = BitMapBackend::new(events_path, (s(WIDTH), s(HEIGHT))).into_drawing_area();
            render_pie_chart(&root, "Splicing Events", &events_slices, SCALE as f64)?;
            root.present()
                .context("Failed to write splice events PNG")?;
        }
        // SVG
        {
            let svg_path = events_path.with_extension("svg");
            let root = SVGBackend::new(&svg_path, (WIDTH, HEIGHT)).into_drawing_area();
            render_pie_chart(&root, "Splicing Events", &events_slices, 1.0)?;
            root.present()
                .context("Failed to write splice events SVG")?;
        }

        info!("Wrote splice events plot: {}", events_png);
    }

    // --- Splice junctions pie chart ---
    {
        let jc = results.junction_counts();
        let (j_known_pct, j_partial_pct, j_novel_pct) = if jc.total > 0 {
            (
                jc.known as f64 * 100.0 / jc.total as f64,
                jc.partial_novel as f64 * 100.0 / jc.total as f64,
                jc.novel as f64 * 100.0 / jc.total as f64,
            )
        } else {
            (0.0, 0.0, 0.0)
        };

        let junctions_slices = [
            ("known", j_known_pct, RGBColor(0, 0, 205)),
            ("partial_novel", j_partial_pct, RGBColor(205, 0, 0)),
            ("complete_novel", j_novel_pct, RGBColor(0, 205, 0)),
        ];

        let junctions_png = format!("{}.splice_junction.png", prefix);
        let junctions_path = Path::new(&junctions_png);

        // PNG
        {
            let root =
                BitMapBackend::new(junctions_path, (s(WIDTH), s(HEIGHT))).into_drawing_area();
            render_pie_chart(&root, "Splicing Junctions", &junctions_slices, SCALE as f64)?;
            root.present()
                .context("Failed to write splice junctions PNG")?;
        }
        // SVG
        {
            let svg_path = junctions_path.with_extension("svg");
            let root = SVGBackend::new(&svg_path, (WIDTH, HEIGHT)).into_drawing_area();
            render_pie_chart(&root, "Splicing Junctions", &junctions_slices, 1.0)?;
            root.present()
                .context("Failed to write splice junctions SVG")?;
        }

        info!("Wrote splice junctions plot: {}", junctions_png);
    }

    Ok(())
}

/// Render a pie chart with labeled slices.
///
/// The `plotters` crate does not have built-in pie chart support, so this
/// is drawn manually using filled arc segments and text labels.
fn render_pie_chart<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    title: &str,
    slices: &[(&str, f64, RGBColor)],
    pxs: f64,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    let ps = |v: f64| (v * pxs) as i32;

    root.fill(&WHITE)?;

    let (w, h) = root.dim_in_pixel();
    let w = w as i32;
    let h = h as i32;

    // Title
    let title_style = TextStyle::from(("sans-serif", ps(16.0) as f64).into_font()).color(&BLACK);
    root.draw_text(title, &title_style, (w / 2 - ps(60.0), ps(10.0)))?;

    let cx = w / 2;
    let cy = h / 2 + ps(10.0);
    let radius = (w.min(h) as f64 * 0.32) as i32;

    // Filter out zero slices
    let total: f64 = slices.iter().map(|(_, pct, _)| pct).sum();
    if total <= 0.0 {
        return Ok(());
    }

    // Draw pie segments using filled polygons
    // RSeQC uses init.angle=30 degrees, so start at 30° from top
    let init_angle: f64 = 30.0_f64.to_radians();
    let mut current_angle = init_angle;

    for &(label, pct, color) in slices {
        if pct <= 0.0 {
            continue;
        }

        let sweep = pct / 100.0 * 2.0 * std::f64::consts::PI;
        let end_angle = current_angle + sweep;

        // Draw filled arc as polygon
        let mut points: Vec<(i32, i32)> = vec![(cx, cy)];
        let n_segments = ((sweep.abs() / 0.02) as usize).max(4);
        for i in 0..=n_segments {
            let angle = current_angle + sweep * (i as f64 / n_segments as f64);
            let px = cx + (radius as f64 * angle.cos()) as i32;
            let py = cy - (radius as f64 * angle.sin()) as i32;
            points.push((px, py));
        }

        // Draw the filled polygon
        root.draw(&Polygon::new(
            points,
            ShapeStyle {
                color: color.mix(0.85),
                filled: true,
                stroke_width: ps(1.0) as u32,
            },
        ))?;

        // Draw label at midpoint of arc, outside the pie
        let mid_angle = current_angle + sweep / 2.0;
        let label_r = radius as f64 * 1.25;
        let lx = cx as f64 + label_r * mid_angle.cos();
        let ly = cy as f64 - label_r * mid_angle.sin();

        let label_text = format!("{} {}%", label, pct.round() as u64);
        let label_style =
            TextStyle::from(("sans-serif", ps(10.0) as f64).into_font()).color(&BLACK);

        // Adjust text position based on which side of the pie
        let text_x = if mid_angle.cos() >= 0.0 {
            lx as i32
        } else {
            lx as i32 - ps(80.0)
        };
        root.draw_text(&label_text, &label_style, (text_x, ly as i32 - ps(5.0)))?;

        current_angle = end_angle;
    }

    Ok(())
}

// ============================================================================
// Junction saturation line plot
// ============================================================================

/// Generate the junction saturation line plot (PNG + SVG).
///
/// Matches the plot produced by RSeQC's junction_saturation.py R script:
/// - Three lines with connected points (type='o' in R)
/// - Blue: all junctions, Red: known junctions, Green: novel junctions
/// - Y-axis values divided by 1000
///
/// # Arguments
/// * `result` — junction saturation analysis results
/// * `output_path` — path for the PNG output; SVG is written alongside
pub fn junction_saturation_plot(result: &SaturationResult, output_path: &Path) -> Result<()> {
    // PNG (high-res)
    {
        let root = BitMapBackend::new(output_path, (s(WIDTH), s(HEIGHT))).into_drawing_area();
        render_junction_saturation(&root, result, SCALE as f64)?;
        root.present()
            .context("Failed to write junction saturation PNG")?;
    }

    // SVG
    let svg_path = output_path.with_extension("svg");
    {
        let root = SVGBackend::new(&svg_path, (WIDTH, HEIGHT)).into_drawing_area();
        render_junction_saturation(&root, result, 1.0)?;
        root.present()
            .context("Failed to write junction saturation SVG")?;
    }

    info!("Wrote junction saturation plot: {}", output_path.display());
    Ok(())
}

/// Render the junction saturation line plot.
fn render_junction_saturation<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    result: &SaturationResult,
    pxs: f64,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    let ps = |v: f64| (v * pxs) as u32;

    root.fill(&WHITE)?;

    if result.percentages.is_empty() {
        return Ok(());
    }

    // Build data series (y values divided by 1000, matching R script)
    let x_vals: Vec<f64> = result.percentages.iter().map(|&p| p as f64).collect();
    let all_y: Vec<f64> = result
        .all_counts
        .iter()
        .map(|&c| c as f64 / 1000.0)
        .collect();
    let known_y: Vec<f64> = result
        .known_counts
        .iter()
        .map(|&c| c as f64 / 1000.0)
        .collect();
    let novel_y: Vec<f64> = result
        .novel_counts
        .iter()
        .map(|&c| c as f64 / 1000.0)
        .collect();

    // Compute axis ranges
    let x_min = x_vals.first().copied().unwrap_or(0.0);
    let x_max = x_vals.last().copied().unwrap_or(100.0);

    // Use true min/max across all three series (not first/last) because
    // individual series may not be monotonic at higher min_coverage values.
    let y_min = all_y
        .iter()
        .chain(known_y.iter())
        .chain(novel_y.iter())
        .copied()
        .fold(f64::INFINITY, f64::min)
        .min(0.0); // ensure 0 is included
    let y_min = if y_min.is_infinite() { 0.0 } else { y_min };
    let y_max = all_y
        .iter()
        .chain(known_y.iter())
        .chain(novel_y.iter())
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let y_max = if y_max.is_infinite() { 1.0 } else { y_max };

    // Add some padding to y range
    let y_range_pad = (y_max - y_min).max(1.0) * 0.05;
    let y_min = (y_min - y_range_pad).max(0.0);
    let y_max = y_max + y_range_pad;

    // --- Build chart ---
    let mut chart = ChartBuilder::on(root)
        .margin(ps(10.0))
        .x_label_area_size(ps(35.0))
        .y_label_area_size(ps(50.0))
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

    chart
        .configure_mesh()
        .x_desc("Percent of total reads")
        .y_desc("Number of splicing junctions (x1000)")
        .label_style(("sans-serif", ps(11.0)))
        .axis_desc_style(("sans-serif", ps(12.0)))
        .draw()?;

    // --- Draw lines with points (matching R's type='o') ---
    let point_size = ps(3.0);

    // Blue: all junctions (drawn first as in R)
    chart
        .draw_series(LineSeries::new(
            x_vals.iter().zip(all_y.iter()).map(|(&x, &y)| (x, y)),
            ShapeStyle {
                color: BLUE.mix(1.0),
                filled: false,
                stroke_width: ps(1.0),
            },
        ))?
        .label("All junctions")
        .legend(move |(x, y)| {
            PathElement::new(
                vec![(x, y), (x + 20, y)],
                ShapeStyle {
                    color: BLUE.mix(1.0),
                    filled: false,
                    stroke_width: 2,
                },
            )
        });

    chart.draw_series(
        x_vals
            .iter()
            .zip(all_y.iter())
            .map(|(&x, &y)| Circle::new((x, y), point_size, BLUE.filled())),
    )?;

    // Red: known junctions
    chart
        .draw_series(LineSeries::new(
            x_vals.iter().zip(known_y.iter()).map(|(&x, &y)| (x, y)),
            ShapeStyle {
                color: RED.mix(1.0),
                filled: false,
                stroke_width: ps(1.0),
            },
        ))?
        .label("Known junctions")
        .legend(move |(x, y)| {
            PathElement::new(
                vec![(x, y), (x + 20, y)],
                ShapeStyle {
                    color: RED.mix(1.0),
                    filled: false,
                    stroke_width: 2,
                },
            )
        });

    chart.draw_series(
        x_vals
            .iter()
            .zip(known_y.iter())
            .map(|(&x, &y)| Circle::new((x, y), point_size, RED.filled())),
    )?;

    // Green: novel junctions
    chart
        .draw_series(LineSeries::new(
            x_vals.iter().zip(novel_y.iter()).map(|(&x, &y)| (x, y)),
            ShapeStyle {
                color: GREEN.mix(1.0),
                filled: false,
                stroke_width: ps(1.0),
            },
        ))?
        .label("Novel junctions")
        .legend(move |(x, y)| {
            PathElement::new(
                vec![(x, y), (x + 20, y)],
                ShapeStyle {
                    color: GREEN.mix(1.0),
                    filled: false,
                    stroke_width: 2,
                },
            )
        });

    chart.draw_series(
        x_vals
            .iter()
            .zip(novel_y.iter())
            .map(|(&x, &y)| Circle::new((x, y), point_size, GREEN.filled())),
    )?;

    // --- Draw legend ---
    chart
        .configure_series_labels()
        .position(SeriesLabelPosition::UpperLeft)
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK.mix(0.5))
        .label_font(("sans-serif", ps(10.0)))
        .draw()?;

    Ok(())
}

// ============================================================================
// Read duplication plot
// ============================================================================

/// Generate the read duplication plot (PNG + SVG).
///
/// Matches the plot produced by RSeQC's read_duplication.py:
/// - Two curves: position-based (blue) and sequence-based (red)
/// - X-axis: occurrence count (1 to 500)
/// - Y-axis: number of uniquely-mapped reads
/// - Connected points (type='o' in R)
///
/// # Arguments
/// * `result` — read duplication analysis results
/// * `output_path` — path for the PNG output; SVG is written alongside
pub fn read_duplication_plot(result: &ReadDuplicationResult, output_path: &Path) -> Result<()> {
    // PNG (high-res)
    {
        let root = BitMapBackend::new(output_path, (s(WIDTH), s(HEIGHT))).into_drawing_area();
        render_read_duplication(&root, result, SCALE as f64)?;
        root.present()
            .context("Failed to write read duplication PNG")?;
    }

    // SVG
    let svg_path = output_path.with_extension("svg");
    {
        let root = SVGBackend::new(&svg_path, (WIDTH, HEIGHT)).into_drawing_area();
        render_read_duplication(&root, result, 1.0)?;
        root.present()
            .context("Failed to write read duplication SVG")?;
    }

    info!("Wrote read duplication plot: {}", output_path.display());
    Ok(())
}

/// Render the read duplication line plot.
fn render_read_duplication<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    result: &ReadDuplicationResult,
    pxs: f64,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    let ps = |v: f64| (v * pxs) as u32;

    root.fill(&WHITE)?;

    // Cap at 500 occurrences (matching RSeQC)
    let max_occ: u64 = 500;

    // Build data points for position-based and sequence-based.
    // DupHistogram is a BTreeMap, so iteration yields keys in sorted order.
    let pos_data: Vec<(f64, f64)> = result
        .pos_histogram
        .iter()
        .filter(|(&occ, _)| occ >= 1 && occ <= max_occ)
        .map(|(&occ, &count)| (occ as f64, count as f64))
        .collect();

    let seq_data: Vec<(f64, f64)> = result
        .seq_histogram
        .iter()
        .filter(|(&occ, _)| occ >= 1 && occ <= max_occ)
        .map(|(&occ, &count)| (occ as f64, count as f64))
        .collect();

    if pos_data.is_empty() && seq_data.is_empty() {
        return Ok(());
    }

    // Compute axis ranges
    let x_max = max_occ as f64;
    let y_max = pos_data
        .iter()
        .chain(seq_data.iter())
        .map(|&(_, y)| y)
        .fold(0.0_f64, f64::max)
        * 1.05;

    let y_max = if y_max <= 0.0 { 1.0 } else { y_max };

    // --- Build chart ---
    let mut chart = ChartBuilder::on(root)
        .margin(ps(10.0))
        .x_label_area_size(ps(35.0))
        .y_label_area_size(ps(55.0))
        .build_cartesian_2d(1.0_f64..x_max, 0.0..y_max)?;

    chart
        .configure_mesh()
        .x_desc("Occurrence")
        .y_desc("Number of uniquely mapped reads")
        .label_style(("sans-serif", ps(11.0)))
        .axis_desc_style(("sans-serif", ps(12.0)))
        .draw()?;

    let point_size = ps(2.0);

    // Blue: sequence-based duplication (matches RSeQC coloring)
    if !seq_data.is_empty() {
        chart
            .draw_series(LineSeries::new(
                seq_data.iter().copied(),
                ShapeStyle {
                    color: BLUE.mix(1.0),
                    filled: false,
                    stroke_width: ps(1.0),
                },
            ))?
            .label("Sequence-based")
            .legend(move |(x, y)| {
                PathElement::new(
                    vec![(x, y), (x + 20, y)],
                    ShapeStyle {
                        color: BLUE.mix(1.0),
                        filled: false,
                        stroke_width: 2,
                    },
                )
            });

        chart.draw_series(
            seq_data
                .iter()
                .map(|&(x, y)| Circle::new((x, y), point_size, BLUE.filled())),
        )?;
    }

    // Red: mapping/position-based duplication (matches RSeQC coloring)
    if !pos_data.is_empty() {
        chart
            .draw_series(LineSeries::new(
                pos_data.iter().copied(),
                ShapeStyle {
                    color: RED.mix(1.0),
                    filled: false,
                    stroke_width: ps(1.0),
                },
            ))?
            .label("Mapping-based")
            .legend(move |(x, y)| {
                PathElement::new(
                    vec![(x, y), (x + 20, y)],
                    ShapeStyle {
                        color: RED.mix(1.0),
                        filled: false,
                        stroke_width: 2,
                    },
                )
            });

        chart.draw_series(
            pos_data
                .iter()
                .map(|&(x, y)| Circle::new((x, y), point_size, RED.filled())),
        )?;
    }

    // --- Draw legend ---
    chart
        .configure_series_labels()
        .position(SeriesLabelPosition::UpperRight)
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK.mix(0.5))
        .label_font(("sans-serif", ps(10.0)))
        .draw()?;

    Ok(())
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_weighted_mean_sd() {
        // Symmetric distribution centered at 0
        let bins = vec![(-10, -5, 10), (-5, 0, 20), (0, 5, 20), (5, 10, 10)];
        let (mean, sd) = compute_weighted_mean_sd(&bins, 5);
        assert!(
            mean.abs() < 1.0,
            "Mean should be near 0 for symmetric data, got {}",
            mean
        );
        assert!(sd > 0.0, "SD should be positive, got {}", sd);
    }

    #[test]
    fn test_compute_weighted_mean_sd_empty() {
        let bins: Vec<(i64, i64, u64)> = vec![];
        let (mean, sd) = compute_weighted_mean_sd(&bins, 5);
        assert_eq!(mean, 0.0);
        assert_eq!(sd, 0.0);
    }

    #[test]
    fn test_compute_density_curve_length() {
        let bins = vec![(-10, -5, 10), (-5, 0, 20), (0, 5, 20), (5, 10, 10)];
        let curve = compute_density_curve(&bins, 60, 10.0, -10.0, 10.0);
        assert_eq!(curve.len(), 512, "Density curve should have 512 points");
        // All density values should be non-negative
        for &(_, y) in &curve {
            assert!(y >= 0.0, "Density should be non-negative, got {}", y);
        }
    }

    #[test]
    fn test_compute_density_curve_integrates_roughly_to_one() {
        let bins = vec![(0, 5, 50), (5, 10, 30), (10, 15, 20)];
        let total: u64 = 100;
        let curve = compute_density_curve(&bins, total, 5.0, 0.0, 15.0);

        // Approximate integral using trapezoidal rule
        let dx = (curve.last().unwrap().0 - curve.first().unwrap().0) / (curve.len() - 1) as f64;
        let integral: f64 = curve.windows(2).map(|w| (w[0].1 + w[1].1) / 2.0 * dx).sum();

        // Should be approximately 1.0 (within tolerance due to edge effects)
        assert!(
            (integral - 1.0).abs() < 0.15,
            "Density should integrate to ~1.0, got {}",
            integral
        );
    }
}
