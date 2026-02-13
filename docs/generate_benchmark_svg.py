#!/usr/bin/env python3
"""Generate benchmark bar chart SVGs for the RustQC documentation.

Creates a horizontal bar chart comparing RustQC performance against
R featureCounts and R dupRadar (standalone and combined), in both
light and dark color schemes.
Output SVGs are saved to docs/public/benchmarks/.

Usage:
    python docs/generate_benchmark_svg.py
"""

import os

# --- Benchmark data ---
# (label, seconds, is_accent)
# Update these values when re-running benchmarks.

BENCHMARK_DATA = [
    ("RustQC", 54, True),
    ("featureCounts", 986, False),
    ("dupRadar", 1796, False),
]

# --- Color schemes ---

THEMES = {
    "dark": {
        "bg": "#1e1e2e",  # Catppuccin Mocha base
        "text": "#cdd6f4",  # Catppuccin Mocha text
        "text_muted": "#a6adc8",  # Catppuccin Mocha subtext0
        "bar_accent": "#cba6f7",  # Catppuccin Mocha mauve
        "bar_other": "#7f849c",  # Catppuccin Mocha overlay1
        "grid": "#313244",  # Catppuccin Mocha surface0
        "label_bold": "#cdd6f4",
    },
    "light": {
        "bg": "#eff1f5",  # Catppuccin Latte base
        "text": "#4c4f69",  # Catppuccin Latte text
        "text_muted": "#6c6f85",  # Catppuccin Latte subtext0
        "bar_accent": "#8839ef",  # Catppuccin Latte mauve
        "bar_other": "#7c7f93",  # Catppuccin Latte overlay1
        "grid": "#ccd0da",  # Catppuccin Latte surface0
        "label_bold": "#4c4f69",
    },
}


def format_time(seconds):
    """Format seconds into a human-readable string."""
    if seconds < 60:
        return f"{seconds}s"
    minutes = seconds // 60
    secs = seconds % 60
    if secs == 0:
        return f"{minutes}m"
    return f"{minutes}m {secs:02d}s"


def generate_svg(data, theme_name, width=720, bar_height=24, bar_gap=8):
    """Generate an SVG horizontal bar chart."""
    theme = THEMES[theme_name]
    n = len(data)

    # Layout
    label_width = 140
    value_label_space = 80
    chart_left = label_width
    chart_right = width - value_label_space
    chart_width = chart_right - chart_left
    top_pad = 10
    bottom_pad = 28
    total_bar_height = bar_height + bar_gap
    chart_height = n * total_bar_height + top_pad + bottom_pad

    max_val = max(d[1] for d in data)

    # Grid line positions (nice round numbers)
    if max_val <= 120:
        step = 30
    elif max_val <= 300:
        step = 60
    elif max_val <= 600:
        step = 120
    elif max_val <= 1200:
        step = 300
    else:
        step = 600
    grid_vals = list(range(0, int(max_val * 1.15) + step, step))

    def x_pos(val):
        scale_max = grid_vals[-1] if grid_vals else max_val
        return chart_left + (val / scale_max) * chart_width

    lines = []
    lines.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {chart_height}" '
        f'width="{width}" height="{chart_height}" '
        f"style=\"font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, "
        f'sans-serif;">'
    )

    # Grid lines
    for gv in grid_vals:
        gx = x_pos(gv)
        lines.append(
            f'  <line x1="{gx:.1f}" y1="{top_pad}" '
            f'x2="{gx:.1f}" y2="{chart_height - bottom_pad}" '
            f'stroke="{theme["grid"]}" stroke-width="1"/>'
        )

    # Grid labels
    for gv in grid_vals:
        gx = x_pos(gv)
        label_y = chart_height - bottom_pad + 16
        lines.append(
            f'  <text x="{gx:.1f}" y="{label_y}" '
            f'text-anchor="middle" fill="{theme["text_muted"]}" '
            f'font-size="12">{format_time(gv)}</text>'
        )

    # Bars and labels
    for i, (label, value, is_accent) in enumerate(data):
        y = top_pad + i * total_bar_height + bar_gap / 2
        bar_w = max(
            (value / (grid_vals[-1] if grid_vals else max_val)) * chart_width, 4
        )
        bar_color = theme["bar_accent"] if is_accent else theme["bar_other"]
        label_weight = "bold" if is_accent and i == 0 else "normal"
        label_color = (
            theme["label_bold"] if is_accent and i == 0 else theme["text_muted"]
        )

        # Label (right-aligned)
        lines.append(
            f'  <text x="{chart_left - 12}" y="{y + bar_height / 2 + 5}" '
            f'text-anchor="end" fill="{label_color}" '
            f'font-size="14" font-weight="{label_weight}">{label}</text>'
        )

        # Bar
        lines.append(
            f'  <rect x="{chart_left}" y="{y:.1f}" '
            f'width="{bar_w:.1f}" height="{bar_height}" '
            f'fill="{bar_color}"/>'
        )

        # Value label (after bar)
        val_x = chart_left + bar_w + 8
        val_color = theme["label_bold"] if is_accent and i == 0 else theme["text"]
        val_weight = "bold" if is_accent and i == 0 else "normal"
        lines.append(
            f'  <text x="{val_x:.1f}" y="{y + bar_height / 2 + 5}" '
            f'text-anchor="start" fill="{val_color}" '
            f'font-size="14" font-weight="{val_weight}">{format_time(value)}</text>'
        )

    lines.append("</svg>")
    return "\n".join(lines)


def main():
    out_dir = os.path.join(os.path.dirname(__file__), "public", "benchmarks")
    os.makedirs(out_dir, exist_ok=True)

    for theme_name in ("dark", "light"):
        svg = generate_svg(BENCHMARK_DATA, theme_name)
        filename = f"benchmark_{theme_name}.svg"
        path = os.path.join(out_dir, filename)
        with open(path, "w") as f:
            f.write(svg)
        print(f"  wrote {path}")

    print("Done.")


if __name__ == "__main__":
    main()
