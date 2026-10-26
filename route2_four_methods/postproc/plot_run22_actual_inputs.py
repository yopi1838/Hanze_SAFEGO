"""Plot the saved 3DEC input velocity tables for Run 22, Methods 1–4."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[2]
METHODS = [
    ("1_fixed_wall", "Method 1", "fixed target · wall-period carrier", "#28779B"),
    ("2_fixed_record", "Method 2", "fixed target · record-period carrier", "#E47B2B"),
    ("3_adaptive_wall", "Method 3", "adaptive target · wall-period carrier", "#3D966B"),
    ("4_adaptive_record", "Method 4", "adaptive target · record-period carrier", "#A64F94"),
]


def load_log_row(method_dir: Path) -> dict[str, str]:
    with (method_dir / "route2_method_log.csv").open(newline="", encoding="utf-8-sig") as f:
        rows = csv.DictReader(f)
        return next(row for row in rows if int(row["run"]) == 22)


def load_table(method_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    matches = list(method_dir.glob("Run22_FR76_s1p00_*_vel.txt"))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected one saved Run 22 velocity table in {method_dir}, found {len(matches)}")
    values = np.loadtxt(matches[0], skiprows=2)
    if values.ndim != 2 or values.shape[1] < 2:
        raise ValueError(f"Unexpected velocity table format: {matches[0]}")
    return values[:, 0], values[:, 1]


def main() -> None:
    series: list[tuple[np.ndarray, np.ndarray]] = []
    info: list[dict[str, str]] = []
    for method_key, _name, _description, _color in METHODS:
        method_dir = ROOT / "route2_four_methods" / method_key
        row = load_log_row(method_dir)
        t, v = load_table(method_dir)
        series.append((t, v))
        info.append(row)

    xmax = max(float(t[-1]) for t, _ in series)
    vmax = max(float(np.max(np.abs(v))) for _, v in series)
    output = ROOT / "route2_four_methods" / "postproc" / "fig_run22_actual_input_velocity_four_methods.png"
    # Draw at 2x resolution for a crisp, portable PNG without relying on an
    # unavailable plotting package in the bundled Python runtime.
    scale = 2
    width, height = 1400, 920
    image = Image.new("RGB", (width * scale, height * scale), "white")
    draw = ImageDraw.Draw(image)
    font_path = Path("C:/Windows/Fonts/segoeui.ttf")
    bold_path = Path("C:/Windows/Fonts/segoeuib.ttf")

    def font(size: int, bold: bool = False):
        path = bold_path if bold else font_path
        try:
            return ImageFont.truetype(str(path), size * scale)
        except OSError:
            return ImageFont.load_default()

    def text(x: float, y: float, label: str, size: int, fill: str = "#222222", bold: bool = False):
        draw.text((round(x * scale), round(y * scale)), label, font=font(size, bold), fill=fill)

    text(72, 35, "Run 22 — actual 3DEC input velocity", 27, bold=True)
    text(74, 82, "Applied table velocity (m/s) · FR76 × 1.00 · native 0.005 s samples · common Wave A carrier = 0.110 s", 13, "#555555")

    ylimit = vmax * 1.14
    panel_positions = [(72, 145), (718, 145), (72, 512), (718, 512)]
    panel_w, panel_h = 610, 340
    tick_x = np.arange(0.0, xmax + 1e-9, 0.1)
    if tick_x[-1] < xmax - 0.03:
        tick_x = np.append(tick_x, xmax)
    tick_y = np.linspace(-ylimit, ylimit, 5)

    for (method_key, name, description, color), (t, v), row, (px, py) in zip(METHODS, series, info, panel_positions):
        carrier_b = float(row["carrier_B_grid_s"])
        duration = float(row["table_duration_s"])
        text(px, py, f"{name}  |  {description}", 15, color, True)
        text(px, py + 27, f"Wave B carrier: {carrier_b:.3f} s     Input: {duration:.3f} s", 11, "#555555")

        left, right = px + 66, px + panel_w - 16
        top, bottom = py + 70, py + panel_h - 47
        plot_w, plot_h = right - left, bottom - top

        def map_x(value: float) -> float:
            return left + (value / xmax) * plot_w

        def map_y(value: float) -> float:
            return top + ((ylimit - value) / (2 * ylimit)) * plot_h

        for val in tick_y:
            yy = map_y(float(val))
            draw.line((round(left * scale), round(yy * scale), round(right * scale), round(yy * scale)), fill="#E1E4E6", width=scale)
            text(px + 4, yy - 8, f"{val:.2f}", 9, "#666666")
        for val in tick_x:
            xx = map_x(float(val))
            draw.line((round(xx * scale), round(top * scale), round(xx * scale), round(bottom * scale)), fill="#ECEEEF", width=scale)
            label = f"{val:.1f}" if abs(float(val) * 10 - round(float(val) * 10)) < 1e-8 else f"{val:.2f}"
            text(xx - 10, bottom + 7, label, 9, "#666666")

        draw.line((round(left * scale), round(top * scale), round(left * scale), round(bottom * scale)), fill="#777777", width=2 * scale)
        draw.line((round(left * scale), round(map_y(0.0) * scale), round(right * scale), round(map_y(0.0) * scale)), fill="#888888", width=scale)
        points = [(round(map_x(float(tx)) * scale), round(map_y(float(vx)) * scale)) for tx, vx in zip(t, v)]
        draw.line(points, fill=color, width=2 * scale, joint="curve")
        text(left + plot_w / 2 - 28, bottom + 27, "Time (s)", 9, "#555555")
        text(right - 140, top + 4, f"peak |v| = {np.max(np.abs(v)):.3f} m/s", 9, "#444444")

    text(74, 875, "Wave A and Wave B are combined in each table. The method sets the carrier periods; the fit procedure adjusts the two component amplitudes.", 11, "#555555")
    image.save(output, format="PNG", dpi=(220, 220), optimize=True)
    print(output)


if __name__ == "__main__":
    main()
