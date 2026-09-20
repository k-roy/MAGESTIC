#!/usr/bin/env python3
"""Shared style + drawing helpers for MAGESTIC documentation schematics.

Figures are hand-authored SVG (crisp at any zoom, diff-able, theme-swappable),
rasterized to PNG at 2x via cairosvg — the same approach RECTIFY's docs use. Each
``generate_<name>.py`` imports this module, builds an SVG string by composing the
helpers below against a palette (light or dark), and calls :func:`render`, which
writes ``<name>.svg`` / ``<name>.png`` (light) and ``<name>_dark.svg`` /
``<name>_dark.png``.

Every figure carries a small grey provenance footer (script filename + render
date + package version) per the lab figure convention, so each PNG is
self-documenting.

Palette is brand-aligned (Indigo-Cyan, matching docs/branding/magestic_mark.svg)
over a restrained slate neutral set for legibility.
"""
from __future__ import annotations

import datetime
import html
import subprocess
from dataclasses import dataclass
from pathlib import Path

FONT = "Inter, Helvetica Neue, Helvetica, Arial, sans-serif"
MONO = "SF Mono, Menlo, Consolas, monospace"


# ───────────────────────────────────────────────────────────── palettes
# Keys are shared between light/dark so a figure body is theme-agnostic.
PAL_LIGHT = dict(
    bg="#ffffff",
    title="#14233b", heading="#33445f", label="#5b6b86", muted="#93a0b5",
    border="#c7d2e0", divider="#e6ecf3",
    # brand families (fill_l = pale fill for chips)
    cyan="#0a8fb0",  cyan_l="#dcf3fb",
    blue="#1f6fe0",  blue_l="#dde9fb",
    indigo="#3b3fd6", indigo_l="#e3e4fb",
    teal="#0fae96",  teal_l="#d3f5ef",   # homology arm / mint
    edit="#e23b78",  edit_l="#fbd9e7",   # the variant edit
    amber="#cf8a1f", amber_l="#fbedcf",
    slate="#33486a", slate_l="#eef2f7",
    # nucleotide-ish neutrals for DNA strands
    dna="#7d8ca6", dna_l="#eaeff5",
    ink_on_fill="#14233b",
)

PAL_DARK = dict(
    bg="#0f172a",
    title="#f1f5f9", heading="#cbd5e1", label="#94a3b8", muted="#64748b",
    border="#33425c", divider="#1e293b",
    cyan="#34c6e6",  cyan_l="#06323f",
    blue="#60a5fa",  blue_l="#11244a",
    indigo="#8b8ef6", indigo_l="#1e1f4a",
    teal="#34d6bd",  teal_l="#053a32",
    edit="#ff6aa0",  edit_l="#43122a",
    amber="#fbbf24", amber_l="#3a2706",
    slate="#cbd5e1", slate_l="#1e293b",
    dna="#9aa8c0", dna_l="#1a2436",
    ink_on_fill="#f1f5f9",
)


@dataclass
class Canvas:
    """Accumulates SVG fragments for one figure at a fixed size + palette."""
    w: int
    h: int
    pal: dict
    parts: list = None

    def __post_init__(self):
        self.parts = []

    def add(self, frag: str):
        self.parts.append(frag)
        return self

    def svg(self, footer: str | None = None) -> str:
        p = self.pal
        head = (
            f'<?xml version="1.0" encoding="utf-8"?>\n'
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.w}" '
            f'height="{self.h}" viewBox="0 0 {self.w} {self.h}">\n'
            f"<defs>"
            f'<marker id="arr" markerWidth="7" markerHeight="7" refX="6" refY="3.5" '
            f'orient="auto" markerUnits="userSpaceOnUse">'
            f'<path d="M0,0 L0,7 L7,3.5 z" fill="{p["label"]}"/></marker>'
            f'<marker id="arrb" markerWidth="7" markerHeight="7" refX="6" refY="3.5" '
            f'orient="auto" markerUnits="userSpaceOnUse">'
            f'<path d="M0,0 L0,7 L7,3.5 z" fill="{p["edit"]}"/></marker>'
            f"</defs>\n"
            f'<rect width="{self.w}" height="{self.h}" fill="{p["bg"]}"/>\n'
            f'<g font-family="{FONT}">\n'
        )
        body = "\n".join(self.parts)
        foot = footer or ""
        return head + body + "\n" + foot + "\n</g>\n</svg>\n"


# ───────────────────────────────────────────────────────────── primitives
def esc(s: str) -> str:
    return html.escape(str(s), quote=True)


def rrect(x, y, w, h, r, fill, stroke=None, sw=1.5, opacity=None, dash=None):
    s = (f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
         f'rx="{r}" ry="{r}" fill="{fill}"')
    if stroke:
        s += f' stroke="{stroke}" stroke-width="{sw}"'
    if dash:
        s += f' stroke-dasharray="{dash}"'
    if opacity is not None:
        s += f' opacity="{opacity}"'
    return s + "/>"


def line(x1, y1, x2, y2, stroke, sw=1.5, dash=None, cap="round"):
    s = (f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
         f'stroke="{stroke}" stroke-width="{sw}" stroke-linecap="{cap}"')
    if dash:
        s += f' stroke-dasharray="{dash}"'
    return s + "/>"


def arrow(x1, y1, x2, y2, stroke, sw=2, dash=None, marker="arr"):
    s = (f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
         f'stroke="{stroke}" stroke-width="{sw}" stroke-linecap="round" '
         f'marker-end="url(#{marker})"')
    if dash:
        s += f' stroke-dasharray="{dash}"'
    return s + "/>"


def text(x, y, s, size=12, color="#000", weight="400", anchor="start",
         family=None, italic=False, spacing=None):
    extra = ""
    if family:
        extra += f' font-family="{family}"'
    if italic:
        extra += ' font-style="italic"'
    if spacing is not None:
        extra += f' letter-spacing="{spacing}"'
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{color}" '
            f'font-weight="{weight}" text-anchor="{anchor}"{extra}>{esc(s)}</text>')


def chip(cx, cy, w, h, label, color, fill, pal, size=12, weight="600",
         sub=None, r=8):
    """Centered rounded chip with a label (and optional sub-label)."""
    x, y = cx - w / 2, cy - h / 2
    out = [rrect(x, y, w, h, r, fill, stroke=color, sw=1.6)]
    if sub:
        out.append(text(cx, cy - 2, label, size=size, color=color,
                        weight=weight, anchor="middle"))
        out.append(text(cx, cy + 13, sub, size=size - 3, color=pal["label"],
                        weight="400", anchor="middle"))
    else:
        out.append(text(cx, cy + size * 0.34, label, size=size, color=color,
                        weight=weight, anchor="middle"))
    return "\n".join(out)


def cli(x, y, cmd, pal, anchor="start"):
    """Monospace CLI command tag."""
    return text(x, y, cmd, size=10.5, color=pal["muted"], family=MONO,
                anchor=anchor)


def protein(cx, cy, w, h, label, color, fill, pal, size=11, weight="700",
            text_color=None):
    """A rounded protein blob with a centered label (LexA, FHA, MCP, RT, …)."""
    r = min(w, h) * 0.38
    out = [rrect(cx - w / 2, cy - h / 2, w, h, r, fill, stroke=color, sw=2)]
    if label:
        out.append(text(cx, cy + size * 0.34, label, size=size,
                        color=text_color or color, weight=weight, anchor="middle"))
    return "\n".join(out)


def arc(cx, cy, r, a0_deg, a1_deg, color, sw, cap="butt"):
    """Stroked ring-segment from a0 to a1 degrees (0=east, CW since y is down)."""
    import math
    a0, a1 = math.radians(a0_deg), math.radians(a1_deg)
    x0, y0 = cx + r * math.cos(a0), cy + r * math.sin(a0)
    x1, y1 = cx + r * math.cos(a1), cy + r * math.sin(a1)
    large = 1 if (a1_deg - a0_deg) % 360 > 180 else 0
    return (f'<path d="M{x0:.2f} {y0:.2f} A{r} {r} 0 {large} 1 {x1:.2f} {y1:.2f}" '
            f'stroke="{color}" stroke-width="{sw}" fill="none" '
            f'stroke-linecap="{cap}"/>')


def stem_loop(x, y, h, color, sw=2.4, w=7):
    """A small RNA hairpin (stem + loop) rooted at (x,y), opening upward."""
    return (f'<path d="M{x:.1f} {y:.1f} L{x:.1f} {y-h:.1f} '
            f'A{w/2:.1f} {w/2:.1f} 0 1 1 {x+w:.1f} {y-h:.1f} L{x+w:.1f} {y:.1f}" '
            f'stroke="{color}" stroke-width="{sw}" fill="none" '
            f'stroke-linejoin="round"/>')


def barcode(x, y, h, pal, widths=(3, 2, 4, 2, 3, 5, 2), colors=None, gap=2):
    """A little barcode glyph (the MAGESTIC motif). Returns (svg, total_width)."""
    if colors is None:
        colors = [pal["teal"], pal["dna"], pal["edit"], pal["dna"],
                  pal["blue"], pal["edit"], pal["teal"]]
    out, cx = [], x
    for w, c in zip(widths, colors):
        out.append(rrect(cx, y, w, h, 1, c))
        cx += w + gap
    return "\n".join(out), cx - x - gap


# ───────────────────────────────────────────────────────────── render
def provenance(w, h, pal, script_name, version="MAGESTIC v3.0.0"):
    today = datetime.date.today().isoformat()
    txt = f"{script_name} · {version} · schematic · rendered {today}"
    return text(w - 8, h - 7, txt, size=6, color=pal["muted"], anchor="end")


def render(build_fn, name: str, outdir: Path, scale: int = 2):
    """Build light + dark SVGs via build_fn(pal) and rasterize each to PNG.

    build_fn must accept a palette dict and return a Canvas. The provenance
    footer + width/height come from the returned Canvas.
    """
    outdir = Path(outdir)
    script = Path(name).name
    for suffix, pal in (("", PAL_LIGHT), ("_dark", PAL_DARK)):
        canvas: Canvas = build_fn(pal)
        foot = provenance(canvas.w, canvas.h, pal, f"{name}.py")
        svg_str = canvas.svg(footer=foot)
        svg_path = outdir / f"{name}{suffix}.svg"
        png_path = outdir / f"{name}{suffix}.png"
        svg_path.write_text(svg_str)
        subprocess.run(
            ["cairosvg", str(svg_path), "-o", str(png_path),
             "--output-width", str(canvas.w * scale)],
            check=True,
        )
        print("wrote", svg_path.name, "+", png_path.name)
