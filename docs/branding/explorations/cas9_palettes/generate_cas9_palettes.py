#!/usr/bin/env python3
"""Develop the Cas9 "Enzyme in Action" mark in clean, professional palettes.

Geometry (kept from the original exploration):
  - genomic DNA duplex with base-pair rungs, cleaved in the middle
  - bilobed Cas9 protein gripping the cut, active-site cleft pointing down
  - a cleavage flash at the break
  - the guide RNA (sgRNA) threading into the active site  -- REFINED here to a
    cleaner, continuous strand (the old tight dashes read as stray dots at small
    sizes)
  - an integrated, trackable barcode at the repaired site

Only COLOR changes across palettes. Each is tuned to be restrained / journal-
clean: one tile family, one enzyme hue with tonal lobes, ONE warm accent for the
cut, a calm strand neutral, and barcode stripes drawn from the same set.

Outputs each tile @512 PNG + SVG and a labeled montage. Needs cairosvg, Pillow.
"""
from __future__ import annotations

import io
from pathlib import Path

import cairosvg
from PIL import Image, ImageDraw, ImageFont

HERE = Path(__file__).resolve().parent

# tile = [top, bottom]; lobe_back/lobe_front = [top,bottom] tonal pair + stroke;
# cut = warm accent (spark outer); spark_core = hot center; guide = sgRNA strand;
# dna = [strand_top, strand_bottom, rung]; bc = 5 barcode colors; on_light = is
# the tile light (=> dark text-safe DNA already handled by dna colors).
PALETTES = {
    "slate_steel": {
        "label": "Slate / Steel-Blue",
        "tile": ["#1b2a38", "#0c151d"],
        "lobe_front": ["#5b9bd5", "#2e6aa6", "#1c4773"],
        "lobe_back": ["#7fb4e0", "#4a86c4", "#2a5a8f"],
        "cut": "#ffb454",
        "spark_core": "#fff0cf",
        "guide": "#86d2ff",
        "dna": ["#e6eef4", "#a9c2d3", "#6f8aa0"],
        "bc": ["#5b9bd5", "#e6eef4", "#86d2ff", "#e6eef4", "#ffb454"],
    },
    "forest": {
        "label": "Forest (refined green)",
        "tile": ["#173a2e", "#0a1c17"],
        "lobe_front": ["#3ddc84", "#16a35a", "#0a5c3c"],
        "lobe_back": ["#6ef0a8", "#2fce98", "#0f8f6e"],
        "cut": "#ffc24b",
        "spark_core": "#fff4d4",
        "guide": "#7fe6d0",
        "dna": ["#e3f1ec", "#a8cbc0", "#6f998c"],
        "bc": ["#3ddc84", "#e3f1ec", "#7fe6d0", "#e3f1ec", "#ffc24b"],
    },
    "clinical_light": {
        "label": "Clinical (light / journal)",
        "tile": ["#f3f6fa", "#e3e9f1"],
        "lobe_front": ["#3f6bb0", "#274a86", "#1b3463"],
        "lobe_back": ["#6f97cf", "#3f6bb0", "#2a4f8f"],
        "cut": "#e2502a",
        "spark_core": "#ffd9b0",
        "guide": "#2b8fd0",
        "dna": ["#46566b", "#7a8aa0", "#aab6c6"],
        "bc": ["#274a86", "#46566b", "#2b8fd0", "#46566b", "#e2502a"],
    },
    "plum_amber": {
        "label": "Plum / Amber",
        "tile": ["#2b2342", "#150f28"],
        "lobe_front": ["#9b7ce0", "#6d4dc0", "#472a8f"],
        "lobe_back": ["#bda3f0", "#8a6ad8", "#5a3eb0"],
        "cut": "#ffb454",
        "spark_core": "#fff0cf",
        "guide": "#7fd6e6",
        "dna": ["#ece8f6", "#bdb4d6", "#8a7fb0"],
        "bc": ["#9b7ce0", "#ece8f6", "#7fd6e6", "#ece8f6", "#ffb454"],
    },
    "teal_copper": {
        "label": "Teal / Copper",
        "tile": ["#0f2e33", "#071719"],
        "lobe_front": ["#2bc8b8", "#159a8e", "#0a5c54"],
        "lobe_back": ["#5ee6d6", "#2bb8a8", "#0f8f80"],
        "cut": "#e88a3c",
        "spark_core": "#ffe0bf",
        "guide": "#86e0ff",
        "dna": ["#e0f0f0", "#a8cccc", "#6f9999"],
        "bc": ["#2bc8b8", "#e0f0f0", "#86e0ff", "#e0f0f0", "#e88a3c"],
    },
}


def _vgrad(id_, top, bottom):
    return (
        f'<linearGradient id="{id_}" x1="0" y1="0" x2="0" y2="1">'
        f'<stop offset="0" stop-color="{top}"/>'
        f'<stop offset="1" stop-color="{bottom}"/></linearGradient>'
    )


def make_svg(p: dict) -> str:
    lf, lb = p["lobe_front"], p["lobe_back"]
    d_top, d_bot, rung = p["dna"]
    bc = p["bc"]
    defs = (
        _vgrad("tile", *p["tile"])
        + _vgrad("lobe", lf[0], lf[1])
        + _vgrad("lobe2", lb[0], lb[1])
        + f'<radialGradient id="spark" cx="0.5" cy="0.5" r="0.5">'
        f'<stop offset="0" stop-color="{p["spark_core"]}"/>'
        f'<stop offset="0.45" stop-color="{p["cut"]}"/>'
        f'<stop offset="1" stop-color="{p["cut"]}" stop-opacity="0"/></radialGradient>'
    )

    def rungs(x0):
        return "".join(
            f'<line x1="{x0+i*16}" y1="159" x2="{x0+i*16}" y2="167"/>' for i in range(5)
        )

    return f"""<svg width="256" height="256" viewBox="0 0 256 256" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="MAGESTIC Cas9 cleavage mark">
<defs>{defs}</defs>
<rect x="8" y="8" width="240" height="240" rx="52" ry="52" fill="url(#tile)"/>

<!-- DNA duplex, cleaved -->
<rect x="22"  y="150" width="86" height="9" rx="4.5" fill="{d_top}"/>
<rect x="22"  y="167" width="86" height="9" rx="4.5" fill="{d_bot}"/>
<g stroke="{rung}" stroke-width="3" stroke-linecap="round">{rungs(34)}</g>
<rect x="148" y="150" width="86" height="9" rx="4.5" fill="{d_top}"/>
<rect x="148" y="167" width="86" height="9" rx="4.5" fill="{d_bot}"/>
<g stroke="{rung}" stroke-width="3" stroke-linecap="round">{rungs(158)}</g>

<!-- guide RNA (sgRNA): clean continuous strand threading into the cleft -->
<path d="M44 56 Q92 78 120 128" stroke="{p['guide']}" stroke-width="5" fill="none"
      stroke-linecap="round" opacity="0.95" stroke-dasharray="11 7"/>

<!-- cleavage flash -->
<circle cx="128" cy="163" r="34" fill="url(#spark)"/>
<path d="M128 136 L121 156 L133 163 L121 170 L128 190"
      stroke="{p['spark_core']}" stroke-width="6" fill="none" stroke-linejoin="round" stroke-linecap="round"/>

<!-- bilobed Cas9 -->
<path d="M70 92 C58 70 78 50 104 54 C128 58 134 82 124 104 C118 118 96 124 84 116 C74 110 74 102 70 92 Z"
      fill="url(#lobe2)" stroke="{lb[2]}" stroke-width="3"/>
<path d="M118 92 C112 66 140 50 168 58 C196 66 200 100 184 124 C172 142 146 148 128 138 C112 130 116 110 118 92 Z"
      fill="url(#lobe)" stroke="{lf[2]}" stroke-width="3"/>
<path d="M120 122 L128 150 L140 124" fill="{p['tile'][1]}" opacity="0.85"/>

<!-- integrated barcode -->
<g>
  <rect x="104" y="198" width="8"  height="36" rx="2.5" fill="{bc[0]}"/>
  <rect x="116" y="198" width="4"  height="36" rx="2"   fill="{bc[1]}"/>
  <rect x="124" y="198" width="11" height="36" rx="2.5" fill="{bc[2]}"/>
  <rect x="139" y="198" width="4"  height="36" rx="2"   fill="{bc[3]}"/>
  <rect x="147" y="198" width="8"  height="36" rx="2.5" fill="{bc[4]}"/>
</g>
</svg>"""


def render(svg_str: str, px: int) -> Image.Image:
    png = cairosvg.svg2png(bytestring=svg_str.encode(), output_width=px, output_height=px)
    return Image.open(io.BytesIO(png)).convert("RGBA")


def main() -> None:
    keys = list(PALETTES)
    size = 256
    pad = 24
    label_h = 34
    ncol = len(keys)

    W = ncol * size + (ncol + 1) * pad
    H = size + 2 * pad + label_h
    sheet = Image.new("RGBA", (W, H), "#0f1117")
    draw = ImageDraw.Draw(sheet)
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 18)
    except OSError:
        font = ImageFont.load_default()

    for i, key in enumerate(keys):
        p = PALETTES[key]
        s = make_svg(p)
        img = render(s, size)
        x = pad + i * (size + pad)
        sheet.paste(img, (x, pad), img)
        lbl = p["label"]
        tb = draw.textbbox((0, 0), lbl, font=font)
        draw.text((x + (size - (tb[2] - tb[0])) / 2, pad + size + 8), lbl,
                  fill="#e6e6ef", font=font)
        (HERE / f"{key}.png").write_bytes(
            cairosvg.svg2png(bytestring=s.encode(), output_width=512, output_height=512)
        )
        (HERE / f"{key}.svg").write_text(s)

    out = HERE / "cas9_palettes_montage.png"
    sheet.convert("RGB").save(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
