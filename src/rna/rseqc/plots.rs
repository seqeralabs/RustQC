//! Native plot generation for RSeQC analyses.
//!
//! Generates PNG + SVG plots directly using the `plotters` crate, replacing the
//! R script intermediate step. Visual output matches the original RSeQC Python
//! tool plots as closely as possible.

use anyhow::{Context, Result};
use log::debug;
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
/// - Mean ± SD as main title (matching upstream R's `main=paste(...)`)
/// - No gridlines, black border around plot area
/// - X-axis range auto-detected from actual data (first/last non-zero bins),
///   matching upstream R's `hist()` auto-scaling behaviour
///
/// # Arguments
/// * `result` — inner distance analysis results
/// * `step` — bin width used for histogram construction
/// * `lower_bound` — inner distance lower bound (fallback x-axis minimum)
/// * `upper_bound` — inner distance upper bound (fallback x-axis maximum)
/// * `sample_name` — sample identifier shown as subtitle below the main title
/// * `output_path` — path for the PNG output; SVG is written alongside
pub fn inner_distance_plot(
    result: &InnerDistanceResult,
    step: i64,
    lower_bound: i64,
    upper_bound: i64,
    sample_name: &str,
    output_path: &Path,
) -> Result<()> {
    // PNG (high-res)
    {
        let root = BitMapBackend::new(output_path, (s(WIDTH), s(HEIGHT))).into_drawing_area();
        render_inner_distance(
            &root,
            result,
            step,
            lower_bound,
            upper_bound,
            sample_name,
            SCALE as f64,
        )?;
        root.present()
            .context("Failed to write inner distance PNG")?;
    }

    // SVG
    let svg_path = output_path.with_extension("svg");
    {
        let root = SVGBackend::new(&svg_path, (WIDTH, HEIGHT)).into_drawing_area();
        render_inner_distance(
            &root,
            result,
            step,
            lower_bound,
            upper_bound,
            sample_name,
            1.0,
        )?;
        root.present()
            .context("Failed to write inner distance SVG")?;
    }

    debug!("Wrote inner distance plot: {}", output_path.display());
    Ok(())
}

/// Render the inner distance histogram with density overlay.
///
/// Replicates the upstream RSeQC R script:
/// ```r
/// hist(fragsize, probability=T, breaks=<num_bins>,
///      xlab="mRNA insert size (bp)",
///      main=paste(c("Mean=", frag_mean, ";", "SD=", frag_sd), collapse=""),
///      border="blue")
/// lines(density(fragsize, bw=<2*step>), col='red')
/// ```
///
/// R's default `hist()` draws no gridlines and `plot()` draws a black box
/// around the plot area. The x-axis range is auto-detected from the actual
/// data (first/last non-zero bins), matching R's `hist()` auto-scaling.
#[allow(clippy::too_many_arguments)]
fn render_inner_distance<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    result: &InnerDistanceResult,
    step: i64,
    lower_bound: i64,
    upper_bound: i64,
    sample_name: &str,
    pxs: f64,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    use plotters::style::text_anchor::{HPos, Pos, VPos};

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

    // X-axis range: auto-detect from actual data, matching upstream RSeQC's
    // behaviour. The upstream R script calls hist(fragsize, ...) without an
    // explicit xlim, so R auto-scales the x-axis to the range of the data.
    // We find the first and last non-zero bins and use their edges.
    let first_nonzero = bins.iter().position(|&(_, _, c)| c > 0);
    let last_nonzero = bins.iter().rposition(|&(_, _, c)| c > 0);
    let (x_min, x_max) = match (first_nonzero, last_nonzero) {
        (Some(first), Some(last)) => (bins[first].0 as f64, bins[last].1 as f64),
        _ => (lower_bound as f64, upper_bound as f64),
    };
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

    // --- Title + subtitle ---
    // Bold descriptive title, then sample name in smaller font underneath.
    // Upstream R uses: main=paste(c("Mean=",frag_mean,";","SD=",frag_sd),collapse="")
    let title_height = ps(50.0);
    let (top_area, plot_area) = root.split_vertically(title_height);
    let cx = (top_area.dim_in_pixel().0 as i32) / 2;
    let title_text = format!("Mean={:.4};SD={:.4}", mean, sd);
    top_area.draw(&Text::new(
        title_text,
        (cx, (pxs * 4.0) as i32),
        ("sans-serif", ps(14.0))
            .into_font()
            .style(FontStyle::Bold)
            .color(&BLACK)
            .pos(Pos::new(HPos::Center, VPos::Top)),
    ))?;
    top_area.draw(&Text::new(
        sample_name.to_string(),
        (cx, (pxs * 22.0) as i32),
        ("sans-serif", ps(11.0))
            .into_font()
            .color(&BLACK)
            .pos(Pos::new(HPos::Center, VPos::Top)),
    ))?;

    // --- Build chart ---
    // No secondary x-axis at the top (upstream R only has a bottom x-axis).
    let mut chart = ChartBuilder::on(&plot_area)
        .margin(ps(10.0))
        .x_label_area_size(ps(35.0))
        .y_label_area_size(ps(50.0))
        .build_cartesian_2d(x_min..x_max, 0.0..y_max)?;

    // Configure mesh: no gridlines, integer axis labels (matching upstream R defaults)
    chart
        .configure_mesh()
        .x_desc("mRNA insert size (bp)")
        .y_desc("Density")
        .label_style(("sans-serif", ps(11.0)))
        .axis_desc_style(("sans-serif", ps(12.0)))
        .x_label_formatter(&|v| format!("{}", *v as i64))
        .disable_mesh() // no gridlines (matching upstream R default)
        .draw()?;

    // Draw a black border around the plot area (matching R's default box())
    chart.draw_series(std::iter::once(Rectangle::new(
        [(x_min, 0.0), (x_max, y_max)],
        ShapeStyle {
            color: BLACK.mix(1.0),
            filled: false,
            stroke_width: ps(1.0),
        },
    )))?;

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
/// - Three categories in R slice order: partial_novel (red), complete_novel (green), known (blue)
/// - Colours match R's default `palette()`: `#DF536B`, `#61D04F`, `#2297E6`
/// - Labels show rounded integer percentages
///
/// # Arguments
/// * `results` — junction annotation analysis results
/// * `prefix` — output path prefix (e.g. `outdir/sample`)
/// * `sample_name` — sample identifier shown as subtitle below each pie chart title
pub fn junction_annotation_plot(
    results: &JunctionResults,
    prefix: &str,
    sample_name: &str,
) -> Result<()> {
    // R default palette colours:  palette()[2] = #DF536B, [3] = #61D04F, [4] = #2297E6
    let col_partial_novel = RGBColor(223, 83, 107); // R palette()[2] — red / pink
    let col_complete_novel = RGBColor(97, 208, 79); // R palette()[3] — green
    let col_known = RGBColor(34, 151, 230); // R palette()[4] — blue

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

        // Slice order matches R: c(partial_novel, complete_novel, known)
        let events_slices = [
            ("partial_novel", e_partial_pct, col_partial_novel),
            ("complete_novel", e_novel_pct, col_complete_novel),
            ("known", e_known_pct, col_known),
        ];

        let events_png = format!("{}.splice_events.png", prefix);
        let events_path = Path::new(&events_png);

        // PNG
        {
            let root = BitMapBackend::new(events_path, (s(WIDTH), s(HEIGHT))).into_drawing_area();
            render_pie_chart(
                &root,
                "Splicing Events",
                sample_name,
                &events_slices,
                SCALE as f64,
            )?;
            root.present()
                .context("Failed to write splice events PNG")?;
        }
        // SVG
        {
            let svg_path = events_path.with_extension("svg");
            let root = SVGBackend::new(&svg_path, (WIDTH, HEIGHT)).into_drawing_area();
            render_pie_chart(&root, "Splicing Events", sample_name, &events_slices, 1.0)?;
            root.present()
                .context("Failed to write splice events SVG")?;
        }

        debug!("Wrote splice events plot: {}", events_png);
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

        // Slice order matches R: c(partial_novel, complete_novel, known)
        let junctions_slices = [
            ("partial_novel", j_partial_pct, col_partial_novel),
            ("complete_novel", j_novel_pct, col_complete_novel),
            ("known", j_known_pct, col_known),
        ];

        let junctions_png = format!("{}.splice_junction.png", prefix);
        let junctions_path = Path::new(&junctions_png);

        // PNG
        {
            let root =
                BitMapBackend::new(junctions_path, (s(WIDTH), s(HEIGHT))).into_drawing_area();
            render_pie_chart(
                &root,
                "Splicing Junctions",
                sample_name,
                &junctions_slices,
                SCALE as f64,
            )?;
            root.present()
                .context("Failed to write splice junctions PNG")?;
        }
        // SVG
        {
            let svg_path = junctions_path.with_extension("svg");
            let root = SVGBackend::new(&svg_path, (WIDTH, HEIGHT)).into_drawing_area();
            render_pie_chart(
                &root,
                "Splicing Junctions",
                sample_name,
                &junctions_slices,
                1.0,
            )?;
            root.present()
                .context("Failed to write splice junctions SVG")?;
        }

        debug!("Wrote splice junctions plot: {}", junctions_png);
    }

    Ok(())
}

/// Render a pie chart with labeled slices.
///
/// The `plotters` crate does not have built-in pie chart support, so this
/// is drawn manually using filled arc segments and text labels.
///
/// Matches R's `pie()` behaviour:
/// - `init.angle = 30` (30° counter-clockwise from 3-o'clock / east)
/// - Leader lines from 1.0× to 1.05× radius
/// - Labels at 1.1× radius, right-aligned on the left side (`adj=1`)
/// - Slices drawn counter-clockwise
fn render_pie_chart<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    title: &str,
    subtitle: &str,
    slices: &[(&str, f64, RGBColor)],
    pxs: f64,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    use plotters::style::text_anchor::{HPos, Pos, VPos};

    let ps = |v: f64| (v * pxs) as i32;

    root.fill(&WHITE)?;

    let (w, h) = root.dim_in_pixel();
    let w = w as i32;
    let h = h as i32;

    // Title — bold, centred horizontally
    let title_style = TextStyle::from(("sans-serif", ps(16.0) as f64).into_font())
        .color(&BLACK)
        .pos(Pos::new(HPos::Center, VPos::Top));
    root.draw_text(title, &title_style, (w / 2, ps(6.0)))?;

    // Subtitle — sample name, smaller font, centred below title
    let subtitle_style = TextStyle::from(("sans-serif", ps(12.0) as f64).into_font())
        .color(&BLACK)
        .pos(Pos::new(HPos::Center, VPos::Top));
    root.draw_text(subtitle, &subtitle_style, (w / 2, ps(24.0)))?;

    let cx = w / 2;
    let cy = h / 2 + ps(14.0);
    let radius = (w.min(h) as f64 * 0.32) as i32;

    // Filter out zero slices
    let total: f64 = slices.iter().map(|(_, pct, _)| pct).sum();
    if total <= 0.0 {
        return Ok(());
    }

    // Draw pie segments using filled polygons.
    // R's pie() with init.angle=30: start 30° counter-clockwise from east.
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

        // Draw the filled polygon (full opacity, matching R)
        root.draw(&Polygon::new(
            points,
            ShapeStyle {
                color: color.to_rgba(),
                filled: true,
                stroke_width: ps(1.0) as u32,
            },
        ))?;

        // Label positioning matches R's pie():
        //   leader line from 1.0× to 1.05× radius
        //   text at 1.1× radius
        //   adj = ifelse(P$x < 0, 1, 0) — right-align on left side
        let mid_angle = current_angle + sweep / 2.0;
        let cos_mid = mid_angle.cos();
        let sin_mid = mid_angle.sin();

        // Leader line: from pie edge (1.0r) to 1.05r
        let line_start_x = cx as f64 + radius as f64 * cos_mid;
        let line_start_y = cy as f64 - radius as f64 * sin_mid;
        let line_end_x = cx as f64 + radius as f64 * 1.05 * cos_mid;
        let line_end_y = cy as f64 - radius as f64 * 1.05 * sin_mid;
        root.draw(&PathElement::new(
            vec![
                (line_start_x as i32, line_start_y as i32),
                (line_end_x as i32, line_end_y as i32),
            ],
            BLACK.stroke_width(ps(1.0) as u32),
        ))?;

        // Text at 1.1× radius, anchored like R's adj parameter
        let label_r = radius as f64 * 1.1;
        let lx = cx as f64 + label_r * cos_mid;
        let ly = cy as f64 - label_r * sin_mid;

        let label_text = format!("{} {}%", label, pct as u64);

        // R: adj = ifelse(P$x < 0, 1, 0)
        //   P$x < 0 ⟹ left side ⟹ adj=1 ⟹ right-aligned
        //   P$x >= 0 ⟹ right side ⟹ adj=0 ⟹ left-aligned
        let h_anchor = if cos_mid < 0.0 {
            HPos::Right
        } else {
            HPos::Left
        };
        let label_style = TextStyle::from(("sans-serif", ps(10.0) as f64).into_font())
            .color(&BLACK)
            .pos(Pos::new(h_anchor, VPos::Center));

        root.draw_text(&label_text, &label_style, (lx as i32, ly as i32))?;

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
/// - Three lines with connected open-circle points (type='o', pch=1 in R)
/// - Blue: all junctions, Red: known junctions, Green: novel junctions
/// - Y-axis values divided by 1000
/// - No gridlines, black border around plot area
/// - X-axis labels every 20% (matching R's default tick selection for 0–100)
/// - Integer axis labels (no decimals)
///
/// # Arguments
/// * `result` — junction saturation analysis results
/// * `sample_name` — sample identifier shown in the plot title
/// * `output_path` — path for the PNG output; SVG is written alongside
pub fn junction_saturation_plot(
    result: &SaturationResult,
    sample_name: &str,
    output_path: &Path,
) -> Result<()> {
    // PNG (high-res)
    {
        let root = BitMapBackend::new(output_path, (s(WIDTH), s(HEIGHT))).into_drawing_area();
        render_junction_saturation(&root, result, sample_name, SCALE as f64)?;
        root.present()
            .context("Failed to write junction saturation PNG")?;
    }

    // SVG
    let svg_path = output_path.with_extension("svg");
    {
        let root = SVGBackend::new(&svg_path, (WIDTH, HEIGHT)).into_drawing_area();
        render_junction_saturation(&root, result, sample_name, 1.0)?;
        root.present()
            .context("Failed to write junction saturation SVG")?;
    }

    debug!("Wrote junction saturation plot: {}", output_path.display());
    Ok(())
}

/// Render the junction saturation line plot.
///
/// Replicates the upstream RSeQC R script:
/// ```r
/// plot(x, z/1000, xlab='percent of total reads',
///      ylab='Number of splicing junctions (x1000)',
///      type='o', col='blue', ylim=c(n,m))
/// points(x, y/1000, type='o', col='red')
/// points(x, w/1000, type='o', col='green')
/// legend(5, m, legend=c("All junctions","known junctions","novel junctions"),
///        col=c("blue","red","green"), lwd=1, pch=1)
/// ```
///
/// R's `type='o'` draws open circles (pch=1) connected by lines.
/// R's default `plot()` draws a black box around the plot area with no gridlines.
fn render_junction_saturation<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    result: &SaturationResult,
    sample_name: &str,
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
        .x_label_area_size(ps(40.0))
        .y_label_area_size(ps(55.0))
        .caption(sample_name, ("sans-serif", ps(14.0)))
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

    // Configure mesh: no gridlines, integer axis labels, x-axis every 20%
    // (matching R's default tick selection for 0–100 range).
    chart
        .configure_mesh()
        .x_desc("percent of total reads")
        .y_desc("Number of splicing junctions (x1000)")
        .label_style(("sans-serif", ps(13.0)))
        .axis_desc_style(("sans-serif", ps(14.0)))
        .x_labels(6) // 0, 20, 40, 60, 80, 100
        .x_label_formatter(&|v| format!("{}", *v as i64))
        .y_label_formatter(&|v| format!("{}", *v as i64))
        .disable_mesh() // no gridlines (matching upstream R default)
        .draw()?;

    // Draw a black border around the plot area (matching R's default box())
    chart.draw_series(std::iter::once(Rectangle::new(
        [(x_min, y_min), (x_max, y_max)],
        ShapeStyle {
            color: BLACK.mix(1.0),
            filled: false,
            stroke_width: ps(1.0),
        },
    )))?;

    // --- Draw lines with open-circle points (matching R's type='o', pch=1) ---
    // R's pch=1 is an open (unfilled) circle.
    let point_size = ps(3.0);
    let stroke = ps(1.0);

    // Blue: all junctions (drawn first as in R)
    chart
        .draw_series(LineSeries::new(
            x_vals.iter().zip(all_y.iter()).map(|(&x, &y)| (x, y)),
            ShapeStyle {
                color: BLUE.mix(1.0),
                filled: false,
                stroke_width: stroke,
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
            .map(|(&x, &y)| Circle::new((x, y), point_size, BLUE.stroke_width(stroke))),
    )?;

    // Red: known junctions
    chart
        .draw_series(LineSeries::new(
            x_vals.iter().zip(known_y.iter()).map(|(&x, &y)| (x, y)),
            ShapeStyle {
                color: RED.mix(1.0),
                filled: false,
                stroke_width: stroke,
            },
        ))?
        .label("known junctions")
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
            .map(|(&x, &y)| Circle::new((x, y), point_size, RED.stroke_width(stroke))),
    )?;

    // Green: novel junctions
    chart
        .draw_series(LineSeries::new(
            x_vals.iter().zip(novel_y.iter()).map(|(&x, &y)| (x, y)),
            ShapeStyle {
                color: GREEN.mix(1.0),
                filled: false,
                stroke_width: stroke,
            },
        ))?
        .label("novel junctions")
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
            .map(|(&x, &y)| Circle::new((x, y), point_size, GREEN.stroke_width(stroke))),
    )?;

    // --- Draw legend ---
    chart
        .configure_series_labels()
        .position(SeriesLabelPosition::UpperLeft)
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK.mix(0.5))
        .label_font(("sans-serif", ps(11.0)))
        .draw()?;

    Ok(())
}

// ============================================================================
// Read duplication plot
// ============================================================================

/// Generate the read duplication plot (PNG + SVG).
///
/// Matches the plot produced by RSeQC's read_duplication.py R script:
/// - Position-based data plotted in blue (× markers), sequence-based in red (● markers)
/// - X-axis: "Occurrence of read" (integer labels at 0, 100, 200, 300, 400, 500)
/// - Y-axis: "Number of Reads (log10)" with log10 exponent labels (0, 1, 2, …)
/// - Legend inside the plot area (matching upstream's `legend(300, ym, ...)`)
/// - No gridlines, black border around plot area
///
/// Note: the upstream RSeQC has a label swap bug in its legend — it labels the
/// blue (position-based) series as "Sequence-based" and the red (sequence-based)
/// series as "Mapping-based". We intentionally use CORRECT labels: blue is
/// "Mapping-based" (position/mapping data) and red is "Sequence-based".
///
/// # Arguments
/// * `result` — read duplication analysis results
/// * `sample_name` — sample identifier shown in the plot title
/// * `output_path` — path for the PNG output; SVG is written alongside
pub fn read_duplication_plot(
    result: &ReadDuplicationResult,
    sample_name: &str,
    output_path: &Path,
) -> Result<()> {
    // PNG (high-res)
    {
        let root = BitMapBackend::new(output_path, (s(WIDTH), s(HEIGHT))).into_drawing_area();
        render_read_duplication(&root, result, sample_name, SCALE as f64)?;
        root.present()
            .context("Failed to write read duplication PNG")?;
    }

    // SVG
    let svg_path = output_path.with_extension("svg");
    {
        let root = SVGBackend::new(&svg_path, (WIDTH, HEIGHT)).into_drawing_area();
        render_read_duplication(&root, result, sample_name, 1.0)?;
        root.present()
            .context("Failed to write read duplication SVG")?;
    }

    debug!("Wrote read duplication plot: {}", output_path.display());
    Ok(())
}

/// Render the read duplication line plot.
///
/// Based on the upstream RSeQC R script:
/// ```r
/// plot(pos_occ, log10(pos_uniqRead), ylab='Number of Reads (log10)',
///      xlab='Occurrence of read', pch=4, cex=0.8, col='blue',
///      xlim=c(1,500), yaxt='n')
/// points(seq_occ, log10(seq_uniqRead), pch=20, cex=0.8, col='red')
/// ym = floor(max(log10(pos_uniqRead)))
/// legend(300, ym, legend=c('Mapping-based','Sequence-based'),
///        col=c('blue','red'), pch=c(4,20))
/// axis(side=2, at=0:ym, labels=0:ym)
/// axis(side=4, at=c(log10(pos_uniqRead[1..4])),
///      labels=c(round(pos_uniqRead[1..4]*100/sum(pos_uniqRead*pos_occ))))
/// mtext(4, text="Reads %", line=2)
/// ```
///
/// Note: the upstream RSeQC legend labels are swapped relative to the plotted data
/// (blue is position-based data labelled "Sequence-based", red is sequence-based
/// data labelled "Mapping-based"). We intentionally use CORRECT labels:
/// blue = "Mapping-based" (position data), red = "Sequence-based".
fn render_read_duplication<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    result: &ReadDuplicationResult,
    sample_name: &str,
    pxs: f64,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    use plotters::style::text_anchor::{HPos, Pos, VPos};

    let ps = |v: f64| (v * pxs) as u32;

    root.fill(&WHITE)?;

    // Cap at 500 occurrences (matching RSeQC's default up_bound)
    let max_occ: u64 = 500;

    // Build data points for position-based and sequence-based.
    // DupHistogram is a BTreeMap, so iteration yields keys in sorted order.
    // Upstream R script applies log10() to the y values for plotting.
    let pos_data: Vec<(f64, f64)> = result
        .pos_histogram
        .iter()
        .filter(|(&occ, _)| occ >= 1 && occ <= max_occ)
        .map(|(&occ, &count)| (occ as f64, (count as f64).log10()))
        .collect();

    let seq_data: Vec<(f64, f64)> = result
        .seq_histogram
        .iter()
        .filter(|(&occ, _)| occ >= 1 && occ <= max_occ)
        .map(|(&occ, &count)| (occ as f64, (count as f64).log10()))
        .collect();

    if pos_data.is_empty() && seq_data.is_empty() {
        return Ok(());
    }

    // Compute axis ranges.
    // x-axis: 0 to 500 (matching upstream's xlim=c(1,500) but R shows 0 as first tick)
    let x_max = max_occ as f64;

    // y-axis: 0 to ym where ym = floor(max(log10(pos_uniqRead)))
    let y_max_raw = pos_data
        .iter()
        .chain(seq_data.iter())
        .map(|&(_, y)| y)
        .fold(f64::NEG_INFINITY, f64::max);
    let ym = if y_max_raw.is_infinite() || y_max_raw < 1.0 {
        1.0
    } else {
        y_max_raw.floor()
    };
    // Add a small margin above ym so the top data points are visible
    let y_plot_max = ym + 0.5;

    // --- Build chart ---
    // Reserve right margin for the secondary "Reads %" axis.
    // Upstream R uses par(mar=c(5,4,4,5)) to allocate right margin space.
    let mut chart = ChartBuilder::on(root)
        .margin(ps(10.0))
        .margin_right(ps(55.0))
        .x_label_area_size(ps(35.0))
        .y_label_area_size(ps(50.0))
        .caption(sample_name, ("sans-serif", ps(14.0)))
        .build_cartesian_2d(0.0_f64..x_max, 0.0_f64..y_plot_max)?;

    // Configure mesh: no gridlines, integer x-axis labels at 0/100/200/300/400/500,
    // y-axis labels as log10 exponents (0, 1, 2, ..., ym).
    chart
        .configure_mesh()
        .x_desc("Occurrence of read")
        .y_desc("Number of Reads (log10)")
        .label_style(("sans-serif", ps(11.0)))
        .axis_desc_style(("sans-serif", ps(12.0)))
        .x_labels(6) // 0, 100, 200, 300, 400, 500
        .x_label_formatter(&|v| format!("{}", *v as i64))
        .y_label_formatter(&|v| format!("{}", *v as i64))
        .disable_mesh() // no gridlines (matching upstream R default)
        .draw()?;

    // Draw a black border around the plot area (matching R's default box())
    chart.draw_series(std::iter::once(Rectangle::new(
        [(0.0_f64, 0.0_f64), (x_max, y_plot_max)],
        ShapeStyle {
            color: BLACK.mix(1.0),
            filled: false,
            stroke_width: ps(1.0),
        },
    )))?;

    let point_size = ps(2.0);

    // Blue: position-based data plotted first (matching upstream R script which
    // calls plot(pos_occ, ..., col='blue') first).
    // Upstream RSeQC mislabels this as "Sequence-based"; we use the correct label.
    if !pos_data.is_empty() {
        chart
            .draw_series(LineSeries::new(
                pos_data.iter().copied(),
                ShapeStyle {
                    color: BLUE.mix(1.0),
                    filled: false,
                    stroke_width: ps(1.0),
                },
            ))?
            .label("Mapping-based")
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
            pos_data
                .iter()
                .map(|&(x, y)| Cross::new((x, y), point_size, BLUE.mix(1.0))),
        )?;
    }

    // Red: sequence-based data plotted second (matching upstream R script which
    // calls points(seq_occ, ..., col='red') second).
    // Upstream RSeQC mislabels this as "Mapping-based"; we use the correct label.
    if !seq_data.is_empty() {
        chart
            .draw_series(LineSeries::new(
                seq_data.iter().copied(),
                ShapeStyle {
                    color: RED.mix(1.0),
                    filled: false,
                    stroke_width: ps(1.0),
                },
            ))?
            .label("Sequence-based")
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
            seq_data
                .iter()
                .map(|&(x, y)| Circle::new((x, y), point_size, RED.filled())),
        )?;
    }

    // --- Draw legend inside the plot area (matching upstream's legend(300, ym, ...)) ---
    chart
        .configure_series_labels()
        .position(SeriesLabelPosition::UpperRight)
        .margin(ps(10.0) as i32)
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK.mix(0.5))
        .label_font(("sans-serif", ps(10.0)))
        .draw()?;

    // -----------------------------------------------------------------------
    // Secondary y-axis: "Reads %" on the right side
    //
    // Replicates the upstream RSeQC R code:
    //   axis(side=4, at=c(log10(pos_uniqRead[1..4])),
    //        labels=c(round(pos_uniqRead[i]*100 / sum(pos_uniqRead*pos_occ))))
    //   mtext(4, text="Reads %", line=2)
    //
    // The tick positions are log10(uniqReadCount) for the first 4 histogram
    // entries, and labels show the percentage of total reads at each level.
    // Total reads = sum(occ * uniqReadCount) across the position histogram.
    // -----------------------------------------------------------------------

    // Compute total reads from position-based histogram:
    // sum(pos_uniqRead * pos_occ) = sum(count * occurrence)
    let total_reads: u64 = result
        .pos_histogram
        .iter()
        .map(|(&occ, &count)| occ * count)
        .sum();

    if total_reads > 0 {
        // Collect the first 4 entries from the position histogram (sorted by occ).
        // These correspond to R's pos_uniqRead[1], pos_uniqRead[2], etc.
        let tick_entries: Vec<(u64, u64)> = result
            .pos_histogram
            .iter()
            .take(4)
            .map(|(&occ, &count)| (occ, count))
            .collect();

        if !tick_entries.is_empty() {
            // Build (log10_y, pct_label) pairs for each tick
            let ticks: Vec<(f64, String)> = tick_entries
                .iter()
                .filter(|&&(_, count)| count > 0)
                .map(|&(_occ, count)| {
                    let log10_y = (count as f64).log10();
                    let pct = (count as f64 * 100.0 / total_reads as f64).round();
                    (log10_y, format!("{}", pct as i64))
                })
                .filter(|(log10_y, _)| *log10_y >= 0.0 && *log10_y <= y_plot_max)
                .collect();

            if !ticks.is_empty() {
                let pa = chart.plotting_area();
                let tick_len = ps(4.0) as i32;
                let label_gap = ps(3.0) as i32;

                // Get the right edge pixel x by mapping the max x-coordinate
                let right_edge_x = pa.map_coordinate(&(x_max, 0.0)).0;

                // Draw the right-side vertical axis line (from y=0 to y=y_plot_max)
                let top_py = pa.map_coordinate(&(x_max, y_plot_max)).1;
                let bot_py = pa.map_coordinate(&(x_max, 0.0)).1;
                root.draw(&PathElement::new(
                    vec![(right_edge_x, top_py), (right_edge_x, bot_py)],
                    BLACK.stroke_width(ps(1.0)),
                ))?;

                // Draw tick marks and percentage labels on the right axis
                let label_style = TextStyle::from(("sans-serif", ps(11.0) as f64).into_font())
                    .color(&BLACK)
                    .pos(Pos::new(HPos::Left, VPos::Center));

                for (log10_y, label) in &ticks {
                    let py = pa.map_coordinate(&(x_max, *log10_y)).1;

                    // Tick mark extending rightward from the plot border
                    root.draw(&PathElement::new(
                        vec![(right_edge_x, py), (right_edge_x + tick_len, py)],
                        BLACK.stroke_width(ps(1.0)),
                    ))?;

                    // Percentage label to the right of the tick
                    root.draw_text(
                        label,
                        &label_style,
                        (right_edge_x + tick_len + label_gap, py),
                    )?;
                }

                // Draw the "Reads %" axis label, rotated 90° on the right side.
                // Positioned to the right of the tick labels, matching R's
                // mtext(4, text="Reads %", line=2).
                let mid_py = (top_py + bot_py) / 2;
                let axis_label_x = right_edge_x + ps(40.0) as i32;
                let axis_label_style = TextStyle::from(("sans-serif", ps(12.0) as f64).into_font())
                    .color(&BLACK)
                    .transform(FontTransform::Rotate90)
                    .pos(Pos::new(HPos::Center, VPos::Center));
                root.draw_text("Reads %", &axis_label_style, (axis_label_x, mid_py))?;
            }
        }
    }

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
