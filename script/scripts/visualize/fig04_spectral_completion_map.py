"""Generate Fig. 4: scoped k>0 spectral completion map (four-topology overview).

Four-panel schematic:
  (a) T^3   — commutative Fourier sector, |k|^2 dispersion
  (b) S^3   — l-symbolic tensor harmonics, (l-1)(l+3)
  (c) Nil^3 — toroidal union Heisenberg Landau tower
  (d) Sol^3 — boost Schrodinger V(w); compact mainline vs noncompact sideline

Output: LaTeX/figures/fig04_spectral_completion_map.svg
"""

from __future__ import annotations

import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER_ROOT = HERE.parents[2]
FIGURE_DIR = PAPER_ROOT / "LaTeX" / "figures"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT = FIGURE_DIR / "fig04_spectral_completion_map.svg"

# Canvas
W, H = 920, 640
MARGIN = 28
GAP = 18
COL_W = (W - 2 * MARGIN - GAP) // 2
ROW_H = (H - 2 * MARGIN - GAP - 52) // 2  # reserve header for H1-H3 bar

# Colors (paper07 palette, aligned with mermaid figures)
C_MAIN = "#2d4a6a"
C_OK = "#2d6a2d"
C_WARN = "#6a5a2d"
C_SIDE = "#6a2d2d"
C_FILL_OK = "#e8f4e8"
C_FILL_MAIN = "#e8f0f8"
C_FILL_WARN = "#f8f4e8"
C_FILL_SIDE = "#f8e8e8"
C_GRID = "#cccccc"
C_TEXT = "#222222"
C_MUTED = "#555555"


def esc(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def rect(x: float, y: float, w: float, h: float, fill: str, stroke: str, sw: float = 1.2) -> str:
    return (
        f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{sw:.1f}" rx="6" ry="6"/>'
    )


def text(
    x: float,
    y: float,
    content: str,
    *,
    size: float = 11,
    anchor: str = "start",
    fill: str = C_TEXT,
    weight: str = "normal",
) -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="DejaVu Sans, Arial, sans-serif" '
        f'font-size="{size:.1f}" font-weight="{weight}" fill="{fill}" '
        f'text-anchor="{anchor}">{esc(content)}</text>'
    )


def line(x1: float, y1: float, x2: float, y2: float, stroke: str, sw: float = 1.2, dash: str | None = None) -> str:
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    return (
        f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
        f'stroke="{stroke}" stroke-width="{sw:.1f}"{dash_attr}/>'
    )


def polyline_points(xs: list[float], ys: list[float]) -> str:
    return " ".join(f"{x:.1f},{y:.1f}" for x, y in zip(xs, ys))


def panel_origin(col: int, row: int) -> tuple[float, float]:
    x0 = MARGIN + col * (COL_W + GAP)
    y0 = MARGIN + 52 + row * (ROW_H + GAP)
    return x0, y0


def draw_axes(px: float, py: float, pw: float, ph: float, xlabel: str, ylabel: str) -> str:
    left = px + 42
    bottom = py + ph - 28
    right = px + pw - 14
    top = py + 36
    parts = [
        line(left, bottom, right, bottom, C_GRID, 1.0),
        line(left, bottom, left, top, C_GRID, 1.0),
        text((left + right) / 2, py + ph - 8, xlabel, size=10, anchor="middle", fill=C_MUTED),
        text(px + 10, (top + bottom) / 2, ylabel, size=10, anchor="middle", fill=C_MUTED),
    ]
    return "\n  ".join(parts)


def panel_frame(col: int, row: int, title: str, fill: str, stroke: str) -> tuple[str, float, float, float, float]:
    x0, y0 = panel_origin(col, row)
    frame = rect(x0, y0, COL_W, ROW_H, fill, stroke)
    title_el = text(x0 + 12, y0 + 18, title, size=11.5, weight="bold")
    return frame + "\n  " + title_el, x0, y0, COL_W, ROW_H


def panel_a() -> str:
    frame, x0, y0, pw, ph = panel_frame(0, 0, "(a) T^3 — commutative Fourier sector", C_FILL_OK, C_OK)
    left = x0 + 42
    bottom = y0 + ph - 28
    right = x0 + pw - 14
    top = y0 + 36
    k_max = 3.0
    lam_max = 9.0

    def map_x(k: float) -> float:
        return left + (k / k_max) * (right - left)

    def map_y(lam: float) -> float:
        return bottom - (lam / lam_max) * (bottom - top)

    xs = [map_x(k) for k in [0.4, 1.0, 1.7, 2.4]]
    ys = [map_y(k * k) for k in [0.4, 1.0, 1.7, 2.4]]
    curve = (
        f'<polyline fill="none" stroke="{C_OK}" stroke-width="2.0" '
        f'points="{polyline_points(xs, ys)}"/>'
    )
    notes = [
        text(x0 + 52, y0 + 52, "lambda ~ |k|^2", size=10, fill=C_OK),
        text(x0 + 52, y0 + 68, "TT=2, positive kinetic", size=9.5, fill=C_MUTED),
        text(x0 + 52, y0 + 84, "no extra zeros for k != 0", size=9.5, fill=C_MUTED),
    ]
    return "\n  ".join([frame, draw_axes(x0, y0, pw, ph, "k", "eigenvalue"), curve, *notes])


def panel_b() -> str:
    frame, x0, y0, pw, ph = panel_frame(1, 0, "(b) S^3 — l-symbolic tensor harmonics", C_FILL_MAIN, C_MAIN)
    left = x0 + 42
    bottom = y0 + ph - 28
    right = x0 + pw - 14
    top = y0 + 36
    l_vals = [2, 3, 4, 5, 6]
    eigs = [(l - 1) * (l + 3) for l in l_vals]
    e_max = max(eigs) * 1.1

    def map_x(i: int) -> float:
        span = right - left
        return left + (i + 0.5) * span / len(l_vals)

    def map_y(e: float) -> float:
        return bottom - (e / e_max) * (bottom - top)

    bars = []
    for i, (l, e) in enumerate(zip(l_vals, eigs)):
        bx = map_x(i) - 14
        by = map_y(e)
        bh = bottom - by
        bars.append(
            rect(bx, by, 28, bh, C_MAIN, C_MAIN, sw=0.8)
            + "\n  "
            + text(map_x(i), bottom + 14, f"l={l}", size=9, anchor="middle", fill=C_MUTED)
        )
    notes = [
        text(x0 + 52, y0 + 52, "eigenvalue = (l-1)(l+3), l >= 2", size=10, fill=C_MAIN),
        text(x0 + 52, y0 + 68, "closed-form all l; C=epsilon frame", size=9.5, fill=C_MUTED),
    ]
    return "\n  ".join([frame, draw_axes(x0, y0, pw, ph, "l", "TT eigenvalue"), *bars, *notes])


def panel_c() -> str:
    frame, x0, y0, pw, ph = panel_frame(0, 1, "(c) Nil^3 — toroidal union Heisenberg tower", C_FILL_WARN, C_WARN)
    left = x0 + 42
    bottom = y0 + ph - 28
    right = x0 + pw - 14
    top = y0 + 36
    k_abs = 2.5

    def map_x(k: float) -> float:
        return left + (k / k_abs) * (right - left)

    def map_y(lam: float) -> float:
        lam_max = 2 * math.pi * k_abs * 5
        return bottom - (lam / lam_max) * (bottom - top)

    # Toroidal branch (k=0): horizontal offset marker
    tor_y = map_y(4 * math.pi ** 2)
    tor = (
        line(left + 8, tor_y, left + 55, tor_y, C_WARN, 1.6)
        + "\n  "
        + text(left + 58, tor_y + 4, "toroidal: 4pi^2(m^2+n^2)", size=9, fill=C_WARN)
    )

    # Heisenberg tower for one k
    tower_x = map_x(1.6)
    levels = [2 * math.pi * 1.6 * (2 * n + 1) for n in range(4)]
    tower = []
    for n, lam in enumerate(levels):
        y = map_y(lam)
        tower.append(line(tower_x - 10, y, tower_x + 10, y, C_WARN, 1.4))
        tower.append(text(tower_x + 14, y + 4, f"N={n}", size=8.5, fill=C_MUTED))
    notes = [
        text(x0 + 52, y0 + 52, "lambda_N(k) = 2pi|k|(2N+1)", size=10, fill=C_WARN),
        text(x0 + 52, y0 + 68, "positive gap at N=0; no contamination", size=9.5, fill=C_MUTED),
    ]
    return "\n  ".join([frame, draw_axes(x0, y0, pw, ph, "|k|", "eigenvalue"), tor, *tower, *notes])


def panel_d() -> str:
    frame, x0, y0, pw, ph = panel_frame(1, 1, "(d) Sol^3 — boost Schrodinger edge", C_FILL_MAIN, C_MAIN)
    left = x0 + 42
    bottom = y0 + ph - 28
    right = x0 + pw - 14
    top = y0 + 36
    w_min, w_max = -2.0, 2.0
    p, q, R = 1.0, 1.0, 1.0

    def V(w: float) -> float:
        return q ** 2 * math.exp(2 * w / R) + p ** 2 * math.exp(-2 * w / R)

    v_min = 2 * abs(p * q)
    v_max = V(w_min) * 1.05
    ws = [w_min + (w_max - w_min) * i / 120 for i in range(121)]

    def map_x(w: float) -> float:
        return left + (w - w_min) / (w_max - w_min) * (right - left)

    def map_y(v: float) -> float:
        return bottom - (v - v_min * 0.85) / (v_max - v_min * 0.85) * (bottom - top)

    xs = [map_x(w) for w in ws]
    ys = [map_y(V(w)) for w in ws]
    curve = (
        f'<polyline fill="none" stroke="{C_MAIN}" stroke-width="2.0" '
        f'points="{polyline_points(xs, ys)}"/>'
    )
    gap_y = map_y(v_min)
    gap = (
        line(left, gap_y, right, gap_y, C_OK, 1.2, dash="5 4")
        + "\n  "
        + text(right - 4, gap_y - 6, "inf V = 2|pq| > 0 (bulk)", size=9, anchor="end", fill=C_OK)
    )

    # Compact mainline box
    main_box = (
        rect(x0 + 52, y0 + 48, pw - 64, 34, C_FILL_OK, C_OK, sw=1.0)
        + "\n  "
        + text(x0 + 60, y0 + 62, "Compact mainline: modified Mathieu, discrete lambda_1 > 0", size=9.2, fill=C_OK)
        + "\n  "
        + text(x0 + 60, y0 + 76, "Anosov gluing; representative numerical edge", size=8.8, fill=C_MUTED)
    )

    # Noncompact sideline (dashed, separate)
    side_box = (
        rect(x0 + 52, y0 + 92, pw - 64, 34, C_FILL_SIDE, C_SIDE, sw=1.0)
        + "\n  "
        + text(x0 + 60, y0 + 106, "Noncompact sideline: pq=0, gapless continuum", size=9.2, fill=C_SIDE)
        + "\n  "
        + text(x0 + 60, y0 + 120, "amenability signature; not mixed into mainline", size=8.8, fill=C_MUTED)
    )

    return "\n  ".join(
        [frame, draw_axes(x0, y0, pw, ph, "w", "V(w)"), curve, gap, main_box, side_box]
    )


def header_bar() -> str:
    x0 = MARGIN
    y0 = MARGIN
    w = W - 2 * MARGIN
    h = 40
    parts = [
        rect(x0, y0, w, h, "#f4f6fa", C_MAIN),
        text(x0 + 14, y0 + 16, "Scoped k>0 spectral completion — verification gates", size=12, weight="bold"),
        text(x0 + 14, y0 + 32, "H1 propagation  |  H2 no-ghost  |  H3 non-contamination of d_zero", size=10, fill=C_MUTED),
    ]
    return "\n  ".join(parts)


def build_svg() -> str:
    body = "\n  ".join([header_bar(), panel_a(), panel_b(), panel_c(), panel_d()])
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">
  <rect width="100%" height="100%" fill="#ffffff"/>
  {body}
</svg>
"""


def main() -> None:
    svg = build_svg()
    OUTPUT.write_text(svg, encoding="utf-8")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()