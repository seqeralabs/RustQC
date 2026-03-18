# Plot Tweaks Report

Progress report for each item in TODO-plot-tweaks.md.

## dupRadar

### Data differences across all 3 plots

**Investigation**: The data differences between RustQC and upstream R dupRadar are primarily due to:

1. **Expression histogram binning**: RustQC used a custom binning approach (rounding bin width
   to nearest 0.05) which produced ~2x more bins than R's `pretty()` algorithm. This made the
   histogram appear jagged/bimodal with peak frequencies around 500 instead of R's ~1000.

2. **Density scatter intercept/slope**: RustQC shows Int: 0.72, Sl: 1.56 vs R's Int: 0.82,
   Sl: 1.54. This is a data-level difference from counting, not a plotting issue. The scatter
   plot visual appearance (density coloring, gradient, overall shape) is very similar.

3. **Boxplot**: Minor differences in last-bin boundary handling (Rust includes max RPK in last
   bin, R excludes it) and NA handling. Visually very similar.

**Status**: The histogram binning was the main visual issue — fixed below.

### Expression histogram bar widths — FIXED

**Problem**: Histogram bins didn't match R's `hist(breaks=100)` output because R treats `breaks`
as a suggestion and uses its `pretty()` algorithm internally to compute "nice" breakpoints.
RustQC was rounding bin width to nearest 0.05, producing different bin counts/widths.

**Fix**: Implemented R's `pretty()` algorithm as `r_pretty_breaks()` function in
`src/rna/dupradar/plots.rs`. The function:
- Computes raw step size as `(hi - lo) / n`
- Rounds to nearest "nice" number from {1, 2, 5, 10} × 10^k
- Floors lo and ceils hi to multiples of the nice step
- Returns evenly-spaced breakpoints

**Result**: New histogram peaks at ~950 (vs R's ~970), matching closely. The old histogram
peaked at ~500 with a jagged/bimodal appearance due to having twice as many bins. The shape,
bar widths, and overall distribution now match the upstream R output.

**Verification**: Built, clippy clean, all 217 tests pass. Generated full-size plot on
GM12878_REP1 data (201.6M reads) and visually confirmed match against R output.

## RSeQC

### read_duplication — FIXED (7 items)

**Problems & Fixes**:

1. **Missing filename title**: Added `sample_name` parameter to `read_duplication_plot()`.
   Title now shows clean sample name (e.g., "GM12878_REP1") matching upstream.

2. **Remove gridlines**: Added `.disable_mesh()` to chart mesh configuration.

3. **Black border around plot area**: Draw a `Rectangle` from (0,0) to (x_max, y_max)
   with black unfilled stroke.

4. **Y-axis log10 labels**: Changed from `.log_scale()` (showing 10, 100, 1000) to
   pre-computing log10 values with a linear axis showing integer exponents (0, 1, 2, ...7).
   Y-axis label changed to "Number of Reads (log10)" matching upstream.

5. **X-axis clean integers**: Changed to 6 labels (0, 100, 200, 300, 400, 500) with
   integer formatting. X range adjusted to 0..500 matching upstream.

6. **Legend position**: Moved inside plot area with margin, matching upstream's overlay style.

7. **Data/label verification**: The upstream RSeQC has a known bug where the legend labels are
   swapped relative to the plotted data (position-based data labeled "Sequence-based" and
   vice versa). Changed RustQC to match upstream's visual output for compatibility:
   - Blue × (Cross) markers → "Sequence-based" (actually position-based data)
   - Red ● (Circle) markers → "Mapping-based" (actually sequence-based data)

**Note**: Upstream also has a right-side "Reads %" axis which RustQC doesn't implement yet.
This wasn't listed in the TODO items.

**Verification**: Built, clippy clean, all 217 tests pass. Generated full-size plot and
confirmed visual match: title, axes, labels, legend placement, marker styles, and data
curves all match upstream closely.

### junction_annotation — FIXED (3 items)

**Problems & Fixes**:

1. **Colours wrong**: Updated to exact R `palette()` colours:
   - `partial_novel`: `#DF536B` / RGB(223, 83, 107) — pink/red
   - `complete_novel`: `#61D04F` / RGB(97, 208, 79) — green
   - `known`: `#2297E6` / RGB(34, 151, 230) — blue
   Also removed 0.85 opacity mixing — now full opacity to match R.

2. **Labels falling off plot area**: Changed label radius from 1.25× to 1.1× (matching R's
   `text(1.1 * P$x, 1.1 * P$y, ...)`). Added proper text anchoring: `HPos::Right` for
   left-side labels, `HPos::Left` for right-side labels (matching R's `adj` parameter).
   Added leader lines from 1.0× to 1.05× radius like R.

3. **Pie chart rotation**: Changed slice order from `known, partial_novel, complete_novel`
   to `partial_novel, complete_novel, known` matching the R script's `c()` order.

**Note**: R uses hatched/density fill patterns; RustQC uses solid fill. This is a minor
cosmetic difference not listed in the TODO. Title casing also differs slightly
("Splicing Events" vs "splicing events").

**Verification**: Built, clippy clean, all 217 tests pass. Visual comparison confirms
colours, label placement, and rotation now match upstream.

### junction_saturation — FIXED (7 items)

**Problems & Fixes**:

1. **Remove grid lines**: Added `.disable_mesh()`.
2. **Black box around plot area**: Draw `Rectangle` border matching R's default `box()`.
3. **Empty markers**: Changed from filled to outline-only circle markers (`stroke_width`
   instead of `filled()`) matching R's `pch=1`.
4. **No decimals on axis labels**: Added integer-format label formatters for both axes.
5. **Title with BAM filename**: Added `sample_name` parameter, shows clean sample name.
6. **X-axis every 20**: Changed to 6 labels (0, 20, 40, 60, 80, 100).
7. **Larger font sizes**: Increased label, axis description, caption, and legend fonts.

Also matched upstream legend labels exactly: "All junctions", "known junctions",
"novel junctions" (lowercase for known/novel).

**Verification**: Built, clippy clean, all 217 tests pass. Visual comparison confirms
all 7 items are fixed: no gridlines, black border, open markers, integer labels, title,
x-axis at every 20, and readable font sizes.

### inner_distance — FIXED (3 items + bonus fixes)

**Problems & Fixes**:

1. **Remove gridlines**: Added `.disable_mesh()`.
2. **X-axis range**: Changed from data-derived range to using `lower_bound`/`upper_bound`
   config parameters (default: -250 to 250), matching the configured analysis bounds.
3. **Plot title**: Changed to "Mean=...;SD=..." format matching upstream R's
   `main=paste(c("Mean=", frag_mean, ";", "SD=", frag_sd), collapse="")`.
   No separate sample name — upstream doesn't show it either.

**Bonus fixes**:
- Added black box border around plot area (matching R defaults)
- Removed spurious secondary x-axis at top (was causing text overlap)
- Cleaned up subtitle layout — Mean/SD is now the main caption

**Verification**: Built, clippy clean, all 217 tests pass. Visual comparison shows
clean layout matching upstream: Mean/SD title, no gridlines, black border, correct
x-axis range, blue histogram with red density curve.
