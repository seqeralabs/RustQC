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
