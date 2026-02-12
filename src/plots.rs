//! Plot generation for dupRust.
//!
//! Generates three types of plots matching dupRadar's R output style:
//! 1. Density scatter plot of duplication rate vs expression (primary plot)
//! 2. Boxplot of duplication rate by expression quantile bins
//! 3. Histogram of expression (RPK) distribution
//!
//! All plots are output as both high-resolution PNG and SVG files,
//! styled to visually match the original R dupRadar package output.

use crate::dupmatrix::DupMatrix;
use crate::fitting::FitResult;
use anyhow::Result;
use plotters::prelude::*;
use plotters_svg::SVGBackend;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Resolution scale factor for PNG output. Base dimensions (640×480) are
/// multiplied by this value to produce high-resolution images.
const SCALE: u32 = 2;

/// Convenience: scale a base pixel value by [`SCALE`].
const fn s(v: u32) -> u32 {
    v * SCALE
}

/// Density color palette matching R's dupRadar `duprateExpDensPlot()` which uses
/// `colorRampPalette(c("cyan", "blue", "green", "yellow", "red"))`.
const DENSITY_COLORS: [(u8, u8, u8); 5] = [
    (0, 255, 255), // cyan       (lowest density)
    (0, 0, 255),   // blue
    (0, 255, 0),   // green
    (255, 255, 0), // yellow
    (255, 0, 0),   // red        (highest density)
];

// ---------------------------------------------------------------------------
// Low-level helpers
// ---------------------------------------------------------------------------

/// Interpolate the density color palette at position `t ∈ [0, 1]`.
fn density_color(t: f64) -> RGBColor {
    let t = t.clamp(0.0, 1.0);
    let n = DENSITY_COLORS.len() - 1;
    let idx = t * n as f64;
    let i = (idx.floor() as usize).min(n - 1);
    let frac = if i == n - 1 && idx >= n as f64 {
        1.0
    } else {
        idx - i as f64
    };
    let r = DENSITY_COLORS[i].0 as f64 * (1.0 - frac) + DENSITY_COLORS[i + 1].0 as f64 * frac;
    let g = DENSITY_COLORS[i].1 as f64 * (1.0 - frac) + DENSITY_COLORS[i + 1].1 as f64 * frac;
    let b = DENSITY_COLORS[i].2 as f64 * (1.0 - frac) + DENSITY_COLORS[i + 1].2 as f64 * frac;
    RGBColor(r as u8, g as u8, b as u8)
}
/// 2-D kernel-density estimation on a grid, matching R's `densCols(nbin=500)`.
///
/// Matches R's `grDevices::densCols()` which internally uses
/// `.smoothScatterCalcDensity()` with bandwidth = IQR_90/25 and
/// `KernSmooth::bkde2D` with `gridsize`, `tau=3.4` truncation.
fn estimate_density(x: &[f64], y: &[f64], nbins: usize) -> Vec<f64> {
    if x.is_empty() {
        return vec![];
    }
    let x_min = x.iter().cloned().fold(f64::INFINITY, f64::min);
    let x_max = x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let y_min = y.iter().cloned().fold(f64::INFINITY, f64::min);
    let y_max = y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let x_range = x_max - x_min;
    let y_range = y_max - y_min;
    if x_range == 0.0 || y_range == 0.0 {
        return vec![1.0; x.len()];
    }

    // R's densCols bandwidth: diff(quantile(x, c(0.05, 0.95))) / 25
    // This is the 90% interquantile range divided by 25.
    let mut x_sorted = x.to_vec();
    x_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut y_sorted = y.to_vec();
    y_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let bw_x = (quantile(&x_sorted, 0.95) - quantile(&x_sorted, 0.05)) / 25.0;
    let bw_y = (quantile(&y_sorted, 0.95) - quantile(&y_sorted, 0.05)) / 25.0;

    // Guard against zero bandwidth
    let bw_x = if bw_x > 0.0 { bw_x } else { x_range / 25.0 };
    let bw_y = if bw_y > 0.0 { bw_y } else { y_range / 25.0 };
    // R's bkde2D extends the grid range by 1.5*bandwidth beyond the data range
    let grid_x_min = x_min - 1.5 * bw_x;
    let grid_x_max = x_max + 1.5 * bw_x;
    let grid_y_min = y_min - 1.5 * bw_y;
    let grid_y_max = y_max + 1.5 * bw_y;
    let grid_x_range = grid_x_max - grid_x_min;
    let grid_y_range = grid_y_max - grid_y_min;

    let x_step = grid_x_range / (nbins - 1) as f64;
    let y_step = grid_y_range / (nbins - 1) as f64;

    // Convert bandwidths to bin units for Gaussian smoothing
    let bw_x_bins = (bw_x / x_step).max(1.0);
    let bw_y_bins = (bw_y / y_step).max(1.0);
    // R's bkde2D uses tau=3.4 as the kernel truncation parameter
    let tau = 3.4;
    let radius_x = (bw_x_bins * tau).ceil() as i32;
    let radius_y = (bw_y_bins * tau).ceil() as i32;
    let sigma2_x = 2.0 * bw_x_bins * bw_x_bins;
    let sigma2_y = 2.0 * bw_y_bins * bw_y_bins;

    // 2-D histogram on grid (nbins x nbins, matching R's bkde2D gridsize)
    let grid_size = nbins;
    let mut grid = vec![vec![0.0f64; grid_size]; grid_size];
    // Use 2D linear binning (linbin2D), matching R's bkde2D:
    // each point contributes proportionally to the 4 surrounding grid cells
    for i in 0..x.len() {
        let fx = (x[i] - grid_x_min) / x_step;
        let fy = (y[i] - grid_y_min) / y_step;
        let ix = fx.floor() as i32;
        let iy = fy.floor() as i32;
        let sx = fx - ix as f64; // fractional x
        let sy = fy - iy as f64; // fractional y
                                 // Distribute weight to 4 corners
        for (dx, wx) in [(0i32, 1.0 - sx), (1, sx)] {
            for (dy, wy) in [(0i32, 1.0 - sy), (1, sy)] {
                let gx = ix + dx;
                let gy = iy + dy;
                if gx >= 0 && (gx as usize) < grid_size && gy >= 0 && (gy as usize) < grid_size {
                    grid[gx as usize][gy as usize] += wx * wy;
                }
            }
        }
    }

    // Anisotropic Gaussian smoothing (matching bkde2D with tau=3.4)
    let mut smoothed = vec![vec![0.0f64; grid_size]; grid_size];
    #[allow(clippy::needless_range_loop)] // bx/by used as integer coordinates for offset arithmetic
    for bx in 0..grid_size {
        for by in 0..grid_size {
            if grid[bx][by] == 0.0 {
                continue;
            }
            let c = grid[bx][by];
            for dx in -radius_x..=radius_x {
                for dy in -radius_y..=radius_y {
                    let nx = bx as i32 + dx;
                    let ny = by as i32 + dy;
                    if nx >= 0 && (nx as usize) < grid_size && ny >= 0 && (ny as usize) < grid_size
                    {
                        let w = (-(dx * dx) as f64 / sigma2_x - (dy * dy) as f64 / sigma2_y).exp();
                        smoothed[nx as usize][ny as usize] += c * w;
                    }
                }
            }
        }
    }

    // Per-point density lookup using nearest grid point (matching R's
    // mkBreaks + cut approach in densCols).
    let mut dens: Vec<f64> = (0..x.len())
        .map(|i| {
            let bx = ((x[i] - grid_x_min) / x_step).round() as usize;
            let by = ((y[i] - grid_y_min) / y_step).round() as usize;
            let bx = bx.min(grid_size - 1);
            let by = by.min(grid_size - 1);
            smoothed[bx][by]
        })
        .collect();
    let mn = dens.iter().cloned().fold(f64::INFINITY, f64::min);
    let mx = dens.iter().cloned().fold(0.0f64, f64::max);
    let range = mx - mn;

    if range > 0.0 {
        for d in dens.iter_mut() {
            *d = (*d - mn) / range; // linear normalization matching R's densCols
        }
    }
    dens
}

/// Format a log₁₀(RPK) tick value as a human-readable label matching R's
/// axis style: `"0.01"`, `"0.1"`, `"1"`, `"10"`, `"100"`, `"1000"`.
fn format_rpk_tick(log_val: f64) -> String {
    let val = 10.0_f64.powf(log_val);
    if val >= 10000.0 {
        // Match R's format: 1e+04, 1e+05
        let exp = log_val.round() as i32;
        format!("1e+{:02}", exp)
    } else if val >= 1.0 {
        format!("{}", val as u64)
    } else if val >= 0.01 {
        let s = format!("{:.2}", val);
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        format!("{:.0e}", val)
    }
}

/// Compute a quantile from a **sorted** slice using linear interpolation.
fn quantile(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    if p <= 0.0 {
        return sorted[0];
    }
    if p >= 1.0 {
        return *sorted.last().unwrap();
    }
    let n = sorted.len();
    let idx = p * (n - 1) as f64;
    let lo = idx.floor() as usize;
    let hi = idx.ceil() as usize;
    let frac = idx.fract();
    if lo == hi {
        sorted[lo]
    } else {
        sorted[lo] * (1.0 - frac) + sorted[hi] * frac
    }
}

/// Draw a dotted vertical line on a pixel-coordinate drawing area.
///
/// plotters has no native dashed-line support, so we draw small segments
/// separated by gaps.
#[allow(clippy::too_many_arguments)]
fn draw_dotted_vline<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    x: i32,
    y_top: i32,
    y_bot: i32,
    color: &RGBColor,
    sw: u32,
    dash: i32,
    gap: i32,
) {
    let mut y = y_top;
    while y < y_bot {
        let ye = (y + dash).min(y_bot);
        root.draw(&PathElement::new(
            vec![(x, y), (x, ye)],
            color.stroke_width(sw),
        ))
        .ok();
        y = ye + gap;
    }
}

// ============================================================================
// Plot 1 – Density scatter
// ============================================================================

/// Internal: render the density scatter plot onto an arbitrary backend.
fn render_density_scatter<DB: DrawingBackend>(
    root: DrawingArea<DB, plotters::coord::Shift>,
    dm: &DupMatrix,
    fit: &FitResult,
    rpkm_threshold: Option<f64>,
    pxs: f64, // pixel scale (SCALE for PNG, 1 for SVG)
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    let ps = |v: f64| (v * pxs) as u32;

    root.fill(&WHITE)?;

    // ── data ───────────────────────────────────────────────────────────
    let mut xd: Vec<f64> = Vec::new();
    let mut yd: Vec<f64> = Vec::new();
    for r in &dm.rows {
        if r.rpk > 0.0 && r.dup_rate.is_finite() {
            xd.push(r.rpk.log10());
            yd.push(r.dup_rate * 100.0);
        }
    }
    if xd.is_empty() {
        anyhow::bail!("No valid data points for density scatter plot");
    }

    let densities = estimate_density(&xd, &yd, 500);

    // X-axis range: R does round(range(log10(RPK))) → e.g. -1..3
    let raw_min = xd.iter().cloned().fold(f64::INFINITY, f64::min);
    let raw_max = xd.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    // R's plot() auto-pads the data range by ~4%.  Use rounded ticks but data-based range.
    let tick_min = raw_min.round() as i32; // for tick positions
    let tick_max = raw_max.round() as i32;
    let pad = (raw_max - raw_min) * 0.04;
    let x_min = (raw_min - pad).max(tick_min as f64); // don't go below first tick
    let x_max = (raw_max + pad).min(tick_max as f64 + 1.0); // a little past last tick

    // ── chart frame ────────────────────────────────────────────────────
    let mut chart = ChartBuilder::on(&root)
        .margin_top(ps(15.0))
        .margin_right(ps(15.0))
        .margin_bottom(ps(10.0))
        .margin_left(ps(10.0))
        .x_label_area_size(ps(40.0))
        .y_label_area_size(ps(55.0))
        .build_cartesian_2d(x_min..x_max, 0.0f64..100.0)?;

    chart
        .configure_mesh()
        .disable_mesh()
        .x_desc("expression (reads/kbp)")
        .y_desc("% duplicate reads")
        .x_label_formatter(&|v| {
            let r = (*v * 10.0).round() / 10.0;
            if (r - r.round()).abs() < 0.01 {
                format_rpk_tick(r)
            } else {
                String::new()
            }
        })
        .y_labels(21) // step=5 on 0-100, then filter to multiples of 25
        .y_label_formatter(&|v| {
            let iv = v.round() as i32;
            if (0..=100).contains(&iv) && iv % 25 == 0 && (*v - iv as f64).abs() < 0.1 {
                format!("{}", iv)
            } else {
                String::new()
            }
        })
        .axis_desc_style(("sans-serif", ps(14.0)))
        .label_style(("sans-serif", ps(12.0)))
        .draw()?;

    // ── points (data order, matching R's plot() behavior) ─────────────
    // R draws points in data order with pch=20 cex=0.25 → tiny filled dots.
    // R's pch=20 with cex=0.25 draws tiny filled circles.
    // Circle radius 1 at our scale gives the closest match.
    for i in 0..xd.len() {
        let c = density_color(densities[i]);
        chart.draw_series(std::iter::once(Circle::new((xd[i], yd[i]), 1, c.filled())))?;
    }

    // ── fit curve: R uses col='black', lwd=2, lty=3 (dotted) ──────────
    // Draw as a series of very short dashes to approximate R's dotted style.
    // R lty=3 draws round dots. We approximate with 2-point-on, 2-point-off.
    let n_curve = (1200.0 * pxs) as usize; // more points at higher resolution
    let curve_pts: Vec<(f64, f64)> = (0..=n_curve)
        .map(|i| {
            let x = x_min + (x_max - x_min) * i as f64 / n_curve as f64;
            let y = (fit.predict(x) * 100.0).clamp(0.0, 100.0);
            (x, y)
        })
        .collect();
    // Draw dotted: short on/off segments
    // R lty=3 (dotted): ~2pt on, ~4pt off. Each curve point spans ~0.3px,
    // so 2pt ≈ 6 points on, 4pt ≈ 12 points off at 2x scale.
    let seg_on_len = (2.0 * pxs).round().max(2.0) as usize;
    let seg_off_len = (4.0 * pxs).round().max(2.0) as usize;
    let mut seg_idx = 0;
    let mut seg_on = true;
    let mut seg_pts: Vec<(f64, f64)> = Vec::new();
    let curve_sw = (pxs * 2.0).round().max(1.0) as u32; // match R lwd=2
    for pt in &curve_pts {
        if seg_on {
            seg_pts.push(*pt);
        }
        seg_idx += 1;
        let limit = if seg_on { seg_on_len } else { seg_off_len };
        if seg_idx >= limit {
            if seg_on && seg_pts.len() >= 2 {
                chart.draw_series(LineSeries::new(
                    std::mem::take(&mut seg_pts),
                    BLACK.stroke_width(curve_sw),
                ))?;
            }
            seg_pts.clear();
            seg_idx = 0;
            seg_on = !seg_on;
        }
    }
    if seg_on && seg_pts.len() >= 2 {
        chart.draw_series(LineSeries::new(seg_pts, BLACK.stroke_width(curve_sw)))?;
    }

    // ── threshold lines (pixel-space dashed lines) ─────────────────────
    let pa = chart.plotting_area();
    let pa_dim = pa.dim_in_pixel();
    let pa_pos = pa.get_base_pixel();
    let pa_x = (pa_pos.0, pa_pos.0 + pa_dim.0 as i32);
    let pa_y = (pa_pos.1, pa_pos.1 + pa_dim.1 as i32);

    let x_to_px = |dx: f64| -> i32 {
        let frac = (dx - x_min) / (x_max - x_min);
        pa_x.0 + (frac * (pa_x.1 - pa_x.0) as f64) as i32
    };

    // R uses lty=2 (dashed): short dash, standard line width
    let dash = (pxs * 4.0) as i32;
    let gap = (pxs * 4.0) as i32;
    let sw = (pxs.round() as u32).max(1); // ~1pt matching R's lty=2 default lwd=1

    // "1 read/bp" → RPK = 1000 → log₁₀ = 3. R: col='red', lty=2 (dashed)
    let rpk_1k = 3.0f64;
    if rpk_1k >= x_min && rpk_1k <= x_max {
        draw_dotted_vline(&root, x_to_px(rpk_1k), pa_y.0, pa_y.1, &RED, sw, dash, gap);
    }

    // RPKM threshold. R: col='green', lty=2 (dashed)
    if let Some(rpk_pos) = rpkm_threshold {
        if rpk_pos > 0.0 {
            let lr = rpk_pos.log10();
            if lr.is_finite() && lr >= x_min && lr <= x_max {
                draw_dotted_vline(&root, x_to_px(lr), pa_y.0, pa_y.1, &GREEN, sw, dash, gap);
            }
        }
    }

    // ── legend: R has Int/Sl at top-left, thresholds at bottom-right ───
    let fsz = ps(11.0);
    let line_h = (pxs * 13.0) as i32;
    let legend_sw = sw; // same thin stroke as threshold lines

    // Top-left legend box: Int/Sl
    {
        let lw = (pxs * 70.0) as i32;
        let lh = (pxs * 32.0) as i32;
        let lx = pa_x.0 + (pxs * 8.0) as i32;
        let ly = pa_y.0 + (pxs * 8.0) as i32;
        root.draw(&Rectangle::new(
            [(lx, ly), (lx + lw, ly + lh)],
            ShapeStyle {
                color: WHITE.to_rgba(),
                filled: true,
                stroke_width: 0,
            },
        ))?;
        root.draw(&Rectangle::new(
            [(lx, ly), (lx + lw, ly + lh)],
            BLACK.stroke_width(legend_sw),
        ))?;
        let tx = lx + (pxs * 6.0) as i32;
        let mut cy = ly + (pxs * 5.0) as i32;
        root.draw(&Text::new(
            format!("Int: {:.2}", fit.intercept),
            (tx, cy),
            ("sans-serif", fsz).into_font().color(&BLACK),
        ))?;
        cy += line_h;
        root.draw(&Text::new(
            format!("Sl: {:.2}", fit.slope),
            (tx, cy),
            ("sans-serif", fsz).into_font().color(&BLACK),
        ))?;
    }

    // Bottom-right legend box: threshold lines
    {
        let lw = (pxs * 90.0) as i32;
        let lh = (pxs * 32.0) as i32;
        let lx = pa_x.1 - lw - (pxs * 8.0) as i32;
        let ly = pa_y.1 - lh - (pxs * 8.0) as i32;
        let txt_x = lx + (pxs * 22.0) as i32;
        let samp_x = lx + (pxs * 4.0) as i32;
        root.draw(&Rectangle::new(
            [(lx, ly), (lx + lw, ly + lh)],
            ShapeStyle {
                color: WHITE.to_rgba(),
                filled: true,
                stroke_width: 0,
            },
        ))?;
        root.draw(&Rectangle::new(
            [(lx, ly), (lx + lw, ly + lh)],
            BLACK.stroke_width(legend_sw),
        ))?;
        let mut cy = ly + (pxs * 5.0) as i32;
        // "1 read/bp" with red dashed sample line
        draw_dotted_vline(
            &root,
            samp_x + (pxs * 7.0) as i32,
            cy,
            cy + line_h - 2,
            &RED,
            sw,
            (pxs * 3.0) as i32,
            (pxs * 2.0) as i32,
        );
        root.draw(&Text::new(
            "1 read/bp",
            (txt_x, cy),
            ("sans-serif", fsz).into_font().color(&BLACK),
        ))?;
        cy += line_h;
        // "0.5 RPKM" with green dashed sample line
        draw_dotted_vline(
            &root,
            samp_x + (pxs * 7.0) as i32,
            cy,
            cy + line_h - 2,
            &GREEN,
            sw,
            (pxs * 3.0) as i32,
            (pxs * 2.0) as i32,
        );
        root.draw(&Text::new(
            "0.5 RPKM",
            (txt_x, cy),
            ("sans-serif", fsz).into_font().color(&BLACK),
        ))?;
    }

    root.present()?;
    Ok(())
}

/// Generate the density scatter plot of duplication rate vs expression.
///
/// Produces both `<output_path>.png` (high-res) and `<output_path>.svg`.
pub fn density_scatter_plot(
    dm: &DupMatrix,
    fit: &FitResult,
    rpkm_threshold: Option<f64>,
    output_path: &std::path::Path,
) -> Result<()> {
    // PNG (high-res)
    let png_path = output_path.with_extension("png");
    let root = BitMapBackend::new(&png_path, (s(480), s(480))).into_drawing_area();
    render_density_scatter(root, dm, fit, rpkm_threshold, SCALE as f64)?;

    // SVG
    let svg_path = output_path.with_extension("svg");
    let root = SVGBackend::new(&svg_path, (480, 480)).into_drawing_area();
    render_density_scatter(root, dm, fit, rpkm_threshold, 1.0)?;

    Ok(())
}

// ============================================================================
// Plot 2 – Boxplot
// ============================================================================

/// Internal: render the boxplot onto an arbitrary backend.
fn render_boxplot<DB: DrawingBackend>(
    root: DrawingArea<DB, plotters::coord::Shift>,
    dm: &DupMatrix,
    pxs: f64,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    let ps = |v: f64| (v * pxs) as u32;

    root.fill(&WHITE)?;

    let step = 0.05;
    let n_bins = (1.0 / step) as usize;

    // Collect RPK for ALL genes (including zeros) — matches R's behaviour
    // where quantile breaks are computed on the full dataset.
    let mut all_rpk: Vec<f64> = dm.rows.iter().map(|r| r.rpk).collect();
    all_rpk.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    if all_rpk.is_empty() {
        anyhow::bail!("No valid data for boxplot");
    }

    // Build (RPK, dupRate) pairs for all genes (dupRate may be NaN for zero-count genes)
    let all_gd: Vec<(f64, f64)> = dm.rows.iter().map(|r| (r.rpk, r.dup_rate)).collect();

    let mut labels: Vec<String> = Vec::new();
    let mut bins: Vec<Vec<f64>> = Vec::new();

    for bi in 0..n_bins {
        let pl = bi as f64 * step;
        let ph = (bi + 1) as f64 * step;
        let ql = quantile(&all_rpk, pl);
        let qh = quantile(&all_rpk, ph);

        // R uses: RPK >= lower_quantile AND RPK < upper_quantile
        // First bin includes the minimum value (lower bound is inclusive)
        let bin_genes: Vec<(f64, f64)> = all_gd
            .iter()
            .filter(|(r, _)| {
                if bi == n_bins - 1 {
                    *r >= ql && *r <= qh // last bin: include upper bound
                } else {
                    *r >= ql && *r < qh
                }
            })
            .copied()
            .collect();

        // Only include finite dupRate values for the boxplot (skip NaN from zero-count genes)
        let vals: Vec<f64> = bin_genes
            .iter()
            .filter(|(_, d)| d.is_finite())
            .map(|(_, d)| *d)
            .collect();

        let mean_rpk: f64 = if bin_genes.is_empty() {
            f64::NAN
        } else {
            let rpk_vals: Vec<f64> = bin_genes.iter().map(|(r, _)| *r).collect();
            rpk_vals.iter().sum::<f64>() / rpk_vals.len() as f64
        };

        let rpk_s = if mean_rpk.is_nan() || mean_rpk == 0.0 {
            "NaN".into()
        } else {
            format!("{:.1}", mean_rpk)
        };
        labels.push(format!(
            "{} - {} % / {}",
            (pl * 100.0) as u32,
            (ph * 100.0) as u32,
            rpk_s
        ));
        bins.push(vals);
    }

    let mut chart = ChartBuilder::on(&root)
        .margin_top(ps(15.0))
        .margin_right(ps(15.0))
        .margin_bottom(ps(10.0))
        .margin_left(ps(10.0))
        .x_label_area_size(ps(150.0))
        .y_label_area_size(ps(50.0))
        .build_cartesian_2d(0.0f64..n_bins as f64, 0.0f64..1.0)?;

    chart
        .configure_mesh()
        .disable_mesh()
        .y_desc("duplication (%)")
        .x_desc("mean expression (reads/kbp)")
        .x_label_formatter(&|_v| String::new()) // Disable auto x labels - we draw manually
        .x_labels(n_bins)
        .y_labels(6)
        .y_label_formatter(&|v| format!("{:.1}", v))
        .axis_desc_style(("sans-serif", ps(13.0)))
        .label_style(("sans-serif", ps(11.0)))
        .draw()?;

    // Draw rotated x-axis labels manually
    {
        let plot_area = chart.plotting_area();
        let font_size = ps(5.0) as f64;
        for (i, label) in labels.iter().enumerate() {
            // Map x position to pixel coordinate (center of bin)
            let x_coord = i as f64 + 0.5;
            let (px, py) = plot_area.map_coordinate(&(x_coord, 0.0f64));
            // Draw text below the plot area, rotated 90 degrees
            let text_style = ("sans-serif", font_size)
                .into_font()
                .color(&BLACK)
                .transform(FontTransform::Rotate270);
            root.draw_text(label, &text_style, (px - 3, py + 5))?;
        }
    }

    let gray_fill = RGBAColor(190, 190, 190, 1.0);

    for (idx, vals) in bins.iter().enumerate() {
        if vals.is_empty() {
            continue;
        }
        let mut sv = vals.clone();
        sv.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let q1 = quantile(&sv, 0.25);
        let med = quantile(&sv, 0.5);
        let q3 = quantile(&sv, 0.75);
        let iqr = q3 - q1;
        let wl = sv
            .iter()
            .find(|&&v| v >= q1 - 1.5 * iqr)
            .copied()
            .unwrap_or(q1);
        let wh = sv
            .iter()
            .rev()
            .find(|&&v| v <= q3 + 1.5 * iqr)
            .copied()
            .unwrap_or(q3);

        let bl = idx as f64 + 0.2;
        let br = idx as f64 + 0.8;
        let cx = idx as f64 + 0.5;
        let cap = 0.1;

        // box fill
        chart.draw_series(std::iter::once(Rectangle::new(
            [(bl, q1), (br, q3)],
            ShapeStyle {
                color: gray_fill,
                filled: true,
                stroke_width: ps(1.0),
            },
        )))?;
        // box border
        chart.draw_series(std::iter::once(Rectangle::new(
            [(bl, q1), (br, q3)],
            BLACK.stroke_width(ps(1.0)),
        )))?;
        // median
        chart.draw_series(LineSeries::new(
            vec![(bl, med), (br, med)],
            BLACK.stroke_width(ps(2.0)),
        ))?;
        // whiskers + caps
        for &(from, to, cap_y) in &[(q1, wl, wl), (q3, wh, wh)] {
            chart.draw_series(LineSeries::new(
                vec![(cx, from), (cx, to)],
                BLACK.stroke_width(ps(1.0)),
            ))?;
            chart.draw_series(LineSeries::new(
                vec![(bl + cap, cap_y), (br - cap, cap_y)],
                BLACK.stroke_width(ps(1.0)),
            ))?;
        }
        // outliers (open circles)
        for &v in &sv {
            if v < wl || v > wh {
                chart.draw_series(std::iter::once(Circle::new(
                    (cx, v),
                    ps(3.0),
                    BLACK.stroke_width(ps(1.0)),
                )))?;
            }
        }
    }

    root.present()?;
    Ok(())
}

/// Generate the duplication rate boxplot by expression quantile bin.
///
/// Produces both `<output_path>.png` (high-res) and `<output_path>.svg`.
pub fn duprate_boxplot(dm: &DupMatrix, output_path: &std::path::Path) -> Result<()> {
    // PNG
    let root = BitMapBackend::new(output_path, (s(480), s(480))).into_drawing_area();
    render_boxplot(root, dm, SCALE as f64)?;
    // SVG
    let svg = output_path.with_extension("svg");
    let root = SVGBackend::new(&svg, (480, 480)).into_drawing_area();
    render_boxplot(root, dm, 1.0)?;
    Ok(())
}

// ============================================================================
// Plot 3 – Expression histogram
// ============================================================================

/// Internal: render the expression histogram onto an arbitrary backend.
fn render_histogram<DB: DrawingBackend>(
    root: DrawingArea<DB, plotters::coord::Shift>,
    dm: &DupMatrix,
    pxs: f64,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    let ps = |v: f64| (v * pxs) as u32;

    root.fill(&WHITE)?;

    let log_rpk: Vec<f64> = dm
        .rows
        .iter()
        .filter(|r| r.rpk > 0.0)
        .map(|r| r.rpk.log10())
        .collect();
    if log_rpk.is_empty() {
        anyhow::bail!("No valid data for expression histogram");
    }

    // Match R's hist(breaks=100): R uses ~100 breaks over the data range.
    // Compute bin width dynamically to produce approximately 100 bins.
    let raw_min = log_rpk.iter().cloned().fold(f64::INFINITY, f64::min);
    let raw_max = log_rpk.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let target_bins = 100;
    let bw = ((raw_max - raw_min) / target_bins as f64 * 20.0).round() / 20.0; // round to nearest 0.05
    let bw = bw.max(0.05); // minimum bin width
    let x_min = (raw_min / bw).floor() * bw;
    let x_max = (raw_max / bw).ceil() * bw;
    let n_bins = ((x_max - x_min) / bw).round() as usize;

    let mut hist = vec![0u32; n_bins];
    for &v in &log_rpk {
        let b = ((v - x_min) / bw).floor() as usize;
        hist[b.min(n_bins - 1)] += 1;
    }
    let y_max = *hist.iter().max().unwrap_or(&1) as f64;

    let mut chart = ChartBuilder::on(&root)
        .margin_top(ps(15.0))
        .margin_right(ps(15.0))
        .margin_bottom(ps(10.0))
        .margin_left(ps(10.0))
        .x_label_area_size(ps(40.0))
        .y_label_area_size(ps(45.0))
        .build_cartesian_2d(x_min..x_max, 0.0..y_max * 1.08)?;

    chart
        .configure_mesh()
        .disable_mesh()
        .x_desc("reads per kilobase (RPK)")
        .y_desc("Frequency")
        .x_label_formatter(&|v| {
            let r = (*v * 10.0).round() / 10.0;
            if (r - r.round()).abs() < 0.01 {
                format_rpk_tick(r)
            } else {
                String::new()
            }
        })
        .y_labels(6)
        .y_label_formatter(&|v| {
            if *v == v.floor() && *v >= 0.0 {
                format!("{}", *v as i32)
            } else {
                String::new()
            }
        })
        .axis_desc_style(("sans-serif", ps(14.0)))
        .label_style(("sans-serif", ps(12.0)))
        .draw()?;

    let gray = RGBAColor(190, 190, 190, 1.0);
    for (i, &c) in hist.iter().enumerate() {
        if c == 0 {
            continue;
        }
        let x0 = x_min + i as f64 * bw;
        let x1 = x0 + bw;
        chart.draw_series(std::iter::once(Rectangle::new(
            [(x0, 0.0), (x1, c as f64)],
            ShapeStyle {
                color: gray,
                filled: true,
                stroke_width: 0,
            },
        )))?;
        chart.draw_series(std::iter::once(Rectangle::new(
            [(x0, 0.0), (x1, c as f64)],
            BLACK.stroke_width(1),
        )))?;
    }

    // Threshold at RPK = 1000 (log₁₀ = 3). R uses col='red'.
    let thr = 3.0f64;
    if thr >= x_min && thr <= x_max {
        chart.draw_series(LineSeries::new(
            vec![(thr, 0.0), (thr, y_max * 1.08)],
            RED.stroke_width(ps(1.0)),
        ))?;
    }

    root.present()?;
    Ok(())
}

/// Generate expression (RPK) histogram.
///
/// Produces both `<output_path>.png` (high-res) and `<output_path>.svg`.
pub fn expression_histogram(dm: &DupMatrix, output_path: &std::path::Path) -> Result<()> {
    // PNG
    let root = BitMapBackend::new(output_path, (s(480), s(480))).into_drawing_area();
    render_histogram(root, dm, SCALE as f64)?;
    // SVG
    let svg = output_path.with_extension("svg");
    let root = SVGBackend::new(&svg, (480, 480)).into_drawing_area();
    render_histogram(root, dm, 1.0)?;
    Ok(())
}

// ============================================================================
// File output helpers
// ============================================================================

/// Write the intercept and slope values (matching R's `label\tvalue` format).
pub fn write_intercept_slope(fit: &FitResult, path: &std::path::Path) -> Result<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "intercept\t{}", fit.intercept)?;
    writeln!(f, "slope\t{}", fit.slope)?;
    Ok(())
}

/// Write a MultiQC-compatible general-stats intercept file.
pub fn write_mqc_intercept(
    fit: &FitResult,
    sample_name: &str,
    path: &std::path::Path,
) -> Result<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "# id: dupRadar")?;
    writeln!(f, "# plot_type: 'generalstats'")?;
    writeln!(f, "# pconfig:")?;
    writeln!(f, "#     dupRadar_intercept:")?;
    writeln!(f, "#         title: 'dupRadar int'")?;
    writeln!(f, "#         namespace: 'dupRadar'")?;
    writeln!(
        f,
        "#         description: 'dupRadar duplication rate at low read counts'"
    )?;
    writeln!(f, "#         max: 100")?;
    writeln!(f, "#         min: 0")?;
    writeln!(f, "#         format: '{{:.2f}}'")?;
    writeln!(f, "Sample\tdupRadar_intercept")?;
    writeln!(f, "{}\t{}", sample_name, fit.intercept)?;
    Ok(())
}

/// Write a MultiQC-compatible line-graph curve file.
pub fn write_mqc_curve(fit: &FitResult, dm: &DupMatrix, path: &std::path::Path) -> Result<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "# id: 'dupradar'")?;
    writeln!(f, "# section_name: 'dupRadar'")?;
    writeln!(f, "# description: 'Duplication rate vs expression'")?;
    writeln!(f, "# plot_type: 'linegraph'")?;
    writeln!(f, "# pconfig:")?;
    writeln!(f, "#     title: 'dupRadar General Linear Model'")?;
    writeln!(f, "#     xlab: 'Expression (reads/kbp)'")?;
    writeln!(f, "#     ylab: 'Duplication rate (%)'")?;
    writeln!(f, "#     ymin: 0")?;
    writeln!(f, "#     ymax: 100")?;
    writeln!(f, "#     xlog: True")?;

    let rpks: Vec<f64> = dm
        .rows
        .iter()
        .filter(|r| r.rpk > 0.0)
        .map(|r| r.rpk)
        .collect();
    if rpks.is_empty() {
        return Ok(());
    }

    let mn = rpks.iter().copied().fold(f64::INFINITY, f64::min);
    let mx = rpks.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let lmin = mn.log10();
    let lmax = mx.log10();
    let n = 100;

    writeln!(f, "RPK\tDuplication Rate (%)")?;
    for i in 0..=n {
        let lr = lmin + (lmax - lmin) * (i as f64) / (n as f64);
        let rpk = 10.0_f64.powf(lr);
        let pct = fit.predict_rpk(rpk) * 100.0;
        writeln!(f, "{}\t{}", rpk, pct)?;
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
    fn test_density_color_endpoints() {
        let c0 = density_color(0.0);
        assert_eq!((c0.0, c0.1, c0.2), (0, 255, 255)); // cyan
        let c1 = density_color(1.0);
        assert_eq!((c1.0, c1.1, c1.2), (255, 0, 0)); // red
                                                     // Mid-point should be green
        let c_mid = density_color(0.5);
        assert_eq!((c_mid.0, c_mid.1, c_mid.2), (0, 255, 0)); // green
    }

    #[test]
    fn test_quantile() {
        let d = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert!((quantile(&d, 0.0) - 1.0).abs() < 1e-10);
        assert!((quantile(&d, 0.5) - 3.0).abs() < 1e-10);
        assert!((quantile(&d, 1.0) - 5.0).abs() < 1e-10);
        assert!((quantile(&d, 0.25) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_estimate_density() {
        let x = vec![0.0, 0.1, 0.0, 0.1, 5.0];
        let y = vec![0.0, 0.1, 0.1, 0.0, 5.0];
        let d = estimate_density(&x, &y, 100);
        assert_eq!(d.len(), 5);
        assert!(d[0] > d[4]);
    }

    #[test]
    fn test_format_rpk_tick() {
        assert_eq!(format_rpk_tick(-2.0), "0.01");
        assert_eq!(format_rpk_tick(-1.0), "0.1");
        assert_eq!(format_rpk_tick(0.0), "1");
        assert_eq!(format_rpk_tick(1.0), "10");
        assert_eq!(format_rpk_tick(2.0), "100");
        assert_eq!(format_rpk_tick(3.0), "1000");
    }
}
