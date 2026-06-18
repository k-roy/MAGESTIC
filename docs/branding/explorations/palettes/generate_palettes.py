#!/usr/bin/env python3
"""Explore color schemes on the MAGESTIC donor-bridge mark.

The two finalists Kevin liked — Concept A (flat, simple) and Tidewave
(vibrant, gradient-rich) — share the SAME geometry: a genomic duplex with a
double-strand break, a donor (homology arms + central edit) recruited down to
bridge the cut, and an integrated barcode below. The only real variable left is
COLOR. This renders that one geometry in two styles x several palettes:

    style FLAT     — solid fills, no shadow/glow (Concept A look; crispest at
                     favicon/avatar size)
    style VIBRANT  — gradients + glow + drop shadow (Tidewave look; pops as a
                     GitHub repo thumbnail / social card)

Outputs each tile as a 512px PNG plus a labeled montage grid (rows = style,
cols = palette).  Run: python generate_palettes.py   (needs cairosvg, Pillow)
"""
from __future__ import annotations

import io
from pathlib import Path

import cairosvg
from PIL import Image, ImageDraw, ImageFont

HERE = Path(__file__).resolve().parent

# Each palette: tile gradient (3 stops), homology-arm color(s), edit color(s),
# recruitment-arrow color, and the neutral "strand"/barcode-light color.
PALETTES = {
    "purple": {
        "label": "Purple (original)",
        "tile": ["#7c4dff", "#6a3fe0", "#4527a0"],
        "arm": ["#1de9b6", "#1de9b6"],
        "edit": ["#ff5d8f", "#ff5d8f"],
        "arrow": "#ffd54f",
        "light": "#ffffff",
    },
    "tidewave": {
        "label": "Tidewave (ocean)",
        "tile": ["#0bd3d3", "#1f6fe0", "#0a1f5c"],
        "arm": ["#affff0", "#23e6c8"],
        "edit": ["#ffd166", "#ff5da2", "#c531ff"],
        "arrow": "#ffe08a",
        "light": "#eafcff",
    },
    "emerald": {
        "label": "Emerald",
        "tile": ["#34e89e", "#15a986", "#0f3443"],
        "arm": ["#d6fff0", "#5ef0c0"],
        "edit": ["#ffe066", "#ff8a3d", "#ff5d5d"],
        "arrow": "#fff0a8",
        "light": "#eafff7",
    },
    "magma": {
        "label": "Magma",
        "tile": ["#ff7a18", "#dd2476", "#3a0ca3"],
        "arm": ["#ffe9c7", "#ffc06b"],
        "edit": ["#9bffe0", "#33e0ff", "#1a8bff"],
        "arrow": "#fff2c2",
        "light": "#fff4e8",
    },
    "cyan": {
        "label": "Indigo-Cyan",
        "tile": ["#00c6ff", "#0072ff", "#101046"],
        "arm": ["#e6fffb", "#5ff0e6"],
        "edit": ["#ffd166", "#ff4d8d", "#ff2e63"],
        "arrow": "#ffe08a",
        "light": "#eafdff",
    },
}


def _grad(id_, colors, vertical=False):
    x2, y2 = ("0", "1") if vertical else ("1", "1")
    stops = "".join(
        f'<stop offset="{i/(len(colors)-1):.3g}" stop-color="{c}"/>'
        for i, c in enumerate(colors)
    )
    return (
        f'<linearGradient id="{id_}" x1="0" y1="0" x2="{x2}" y2="{y2}">'
        f"{stops}</linearGradient>"
    )


def make_svg(p: dict, vibrant: bool) -> str:
    """Render the donor-bridge geometry for palette p in flat or vibrant style."""
    if vibrant:
        defs = (
            _grad("bg", p["tile"])
            + _grad("edit", p["edit"], vertical=True)
            + _grad("arm", p["arm"], vertical=True)
            + _grad("strand", [p["light"], p["arm"][-1]])
            + f'<radialGradient id="glow" cx="0.5" cy="0.42" r="0.6">'
            f'<stop offset="0" stop-color="{p["edit"][len(p["edit"])//2]}" stop-opacity="0.5"/>'
            f'<stop offset="1" stop-color="{p["edit"][len(p["edit"])//2]}" stop-opacity="0"/>'
            "</radialGradient>"
            '<filter id="soft" x="-30%" y="-30%" width="160%" height="160%">'
            '<feDropShadow dx="0" dy="3" stdDeviation="3.5" flood-color="#06122e" flood-opacity="0.5"/>'
            "</filter>"
        )
        arm_fill, edit_fill, strand_fill, light = (
            "url(#arm)",
            "url(#edit)",
            "url(#strand)",
            p["light"],
        )
        sheen = '<rect x="8" y="8" width="240" height="120" rx="54" ry="54" fill="#ffffff" opacity="0.06"/>'
        glow = '<circle cx="128" cy="118" r="92" fill="url(#glow)"/>'
        fx = ' filter="url(#soft)"'
    else:
        defs = _grad("bg", p["tile"])
        arm_fill, edit_fill, strand_fill, light = (
            p["arm"][-1],
            p["edit"][len(p["edit"]) // 2],
            light_strand(p),
            p["light"],
        )
        sheen = glow = ""
        fx = ""

    arrow = p["arrow"]
    return f"""<svg width="256" height="256" viewBox="0 0 256 256" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="MAGESTIC">
<defs>{defs}</defs>
<rect x="8" y="8" width="240" height="240" rx="54" ry="54" fill="url(#bg)"/>
{sheen}{glow}
<g{fx}>
  <rect x="26"  y="156" width="76" height="11" rx="5.5" fill="{strand_fill}" opacity="0.95"/>
  <rect x="154" y="156" width="76" height="11" rx="5.5" fill="{strand_fill}" opacity="0.95"/>
  <rect x="26"  y="174" width="76" height="11" rx="5.5" fill="{strand_fill}" opacity="0.6"/>
  <rect x="154" y="174" width="76" height="11" rx="5.5" fill="{strand_fill}" opacity="0.6"/>
</g>
<g{fx}>
  <rect x="86"  y="86" width="26" height="20" rx="6" fill="{arm_fill}"/>
  <rect x="115" y="80" width="26" height="26" rx="6" fill="{edit_fill}"/>
  <rect x="144" y="86" width="26" height="20" rx="6" fill="{arm_fill}"/>
</g>
<path d="M99 110 Q104 134 118 152" stroke="{arrow}" stroke-width="6" fill="none" stroke-linecap="round" opacity="0.95"/>
<path d="M157 110 Q152 134 138 152" stroke="{arrow}" stroke-width="6" fill="none" stroke-linecap="round" opacity="0.95"/>
<g{fx}>
  <rect x="104" y="196" width="9"  height="38" rx="3" fill="{arm_fill}"/>
  <rect x="117" y="196" width="6"  height="38" rx="3" fill="{light}"/>
  <rect x="127" y="196" width="11" height="38" rx="3" fill="{edit_fill}"/>
  <rect x="142" y="196" width="7"  height="38" rx="3" fill="{strand_fill}"/>
</g>
</svg>"""


def light_strand(p):
    """Flat style: a single readable strand color (use the light neutral)."""
    return p["light"]


def render(svg_str: str, px: int) -> Image.Image:
    png = cairosvg.svg2png(bytestring=svg_str.encode(), output_width=px, output_height=px)
    return Image.open(io.BytesIO(png)).convert("RGBA")


def main() -> None:
    styles = [("flat", False), ("vibrant", True)]
    size = 248
    pad = 22
    col_label_h = 30
    row_label_w = 86

    keys = list(PALETTES)
    ncol = len(keys)
    nrow = len(styles)

    W = row_label_w + ncol * size + (ncol + 1) * pad
    H = col_label_h + nrow * size + (nrow + 1) * pad
    sheet = Image.new("RGBA", (W, H), "#0f1117")
    draw = ImageDraw.Draw(sheet)

    def font(sz, bold=True):
        name = "Arial Bold.ttf" if bold else "Arial.ttf"
        try:
            return ImageFont.truetype(f"/System/Library/Fonts/Supplemental/{name}", sz)
        except OSError:
            return ImageFont.load_default()

    fcol, frow = font(18), font(16)

    # column headers (palette names)
    for c, key in enumerate(keys):
        x = row_label_w + pad + c * (size + pad)
        lbl = PALETTES[key]["label"]
        tb = draw.textbbox((0, 0), lbl, font=fcol)
        draw.text((x + (size - (tb[2] - tb[0])) / 2, 6), lbl, fill="#e6e6ef", font=fcol)

    for r, (sname, vib) in enumerate(styles):
        y = col_label_h + pad + r * (size + pad)
        # row label (rotated-ish: just place at left, vertically centered)
        draw.text((10, y + size / 2 - 8), sname.upper(), fill="#9aa0b5", font=frow)
        for c, key in enumerate(keys):
            p = PALETTES[key]
            s = make_svg(p, vib)
            img = render(s, size)
            x = row_label_w + pad + c * (size + pad)
            sheet.paste(img, (x, y), img)
            # save individual tile @512
            (HERE / f"{key}_{sname}.png").write_bytes(
                cairosvg.svg2png(bytestring=s.encode(), output_width=512, output_height=512)
            )
            (HERE / f"{key}_{sname}.svg").write_text(s)

    out = HERE / "palettes_montage.png"
    sheet.convert("RGB").save(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
