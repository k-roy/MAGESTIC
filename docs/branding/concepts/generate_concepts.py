#!/usr/bin/env python3
"""Explore alternative MAGESTIC mark concepts.

Each concept renders the same deep-purple "cell" tile but a different inner
motif, all trying to capture two ideas the current mark under-states:

    1. EDITING  — a CRISPR double-strand break (DSB) at the target locus
    2. DONOR RECRUITMENT — a homology-arm donor brought to the cut site to
       template the precise edit, leaving an integrated barcode

Outputs: one PNG per concept + a labeled side-by-side montage for comparison.
Run:  python docs/branding/concepts/generate_concepts.py
Requires: cairosvg, Pillow.
"""
from __future__ import annotations

import io
from pathlib import Path

import cairosvg
from PIL import Image, ImageDraw, ImageFont

HERE = Path(__file__).resolve().parent

PURPLE_A = "#7c4dff"
PURPLE_B = "#4527a0"
TEAL = "#1de9b6"
AMBER = "#ffd54f"
WHITE = "#ffffff"
CORAL = "#ff5d8f"

TILE = """
  <defs>
    <linearGradient id="bg" x1="0" y1="0" x2="1" y2="1">
      <stop offset="0" stop-color="{a}"/>
      <stop offset="1" stop-color="{b}"/>
    </linearGradient>
  </defs>
  <rect x="8" y="8" width="240" height="240" rx="52" ry="52" fill="url(#bg)"/>
""".format(a=PURPLE_A, b=PURPLE_B)


def svg(body: str) -> str:
    return (
        '<svg width="256" height="256" viewBox="0 0 256 256" '
        'xmlns="http://www.w3.org/2000/svg">' + TILE + body + "</svg>"
    )


# ---------------------------------------------------------------- concept A
# "Donor bridge": a DSB gap in the genomic duplex; a donor with matching
# homology-arm ends descends to bridge the cut, carrying the edit; barcode below.
CONCEPT_A = svg(f"""
  <!-- genomic duplex with a clean DSB gap in the middle -->
  <rect x="30"  y="150" width="78" height="7" rx="3.5" fill="{WHITE}" opacity="0.9"/>
  <rect x="148" y="150" width="78" height="7" rx="3.5" fill="{WHITE}" opacity="0.9"/>
  <rect x="30"  y="163" width="78" height="7" rx="3.5" fill="{WHITE}" opacity="0.5"/>
  <rect x="148" y="163" width="78" height="7" rx="3.5" fill="{WHITE}" opacity="0.5"/>
  <!-- donor with homology arms (teal) + central edit (coral) recruited above -->
  <g>
    <rect x="92"  y="92" width="22" height="14" rx="4" fill="{TEAL}"/>
    <rect x="116" y="92" width="24" height="14" rx="4" fill="{CORAL}"/>
    <rect x="142" y="92" width="22" height="14" rx="4" fill="{TEAL}"/>
  </g>
  <!-- recruitment arrows from donor to cut -->
  <path d="M104 110 Q108 132 120 146" stroke="{AMBER}" stroke-width="5" fill="none" stroke-linecap="round"/>
  <path d="M152 110 Q148 132 136 146" stroke="{AMBER}" stroke-width="5" fill="none" stroke-linecap="round"/>
  <!-- integrated barcode below the repaired site -->
  <rect x="112" y="184" width="7"  height="34" rx="2" fill="{TEAL}"/>
  <rect x="123" y="184" width="4"  height="34" rx="2" fill="{WHITE}"/>
  <rect x="131" y="184" width="9"  height="34" rx="2" fill="{AMBER}"/>
""")


# ---------------------------------------------------------------- concept B
# "Cas9 wedge": a PacMan-style nuclease clamped on the duplex at the cut, guide
# RNA inside; donor recruited on a tether; barcode integrated.
CONCEPT_B = svg(f"""
  <!-- genomic duplex -->
  <rect x="30" y="150" width="196" height="7" rx="3.5" fill="{WHITE}" opacity="0.85"/>
  <rect x="30" y="163" width="196" height="7" rx="3.5" fill="{WHITE}" opacity="0.45"/>
  <!-- Cas9 wedge (PacMan) at cut site -->
  <path d="M128 96
           A40 40 0 1 0 162 144
           L128 136 Z" fill="{AMBER}"/>
  <!-- guide RNA seed inside -->
  <rect x="112" y="120" width="18" height="6" rx="3" fill="{PURPLE_B}"/>
  <!-- donor recruited on a tether from upper-left -->
  <rect x="36" y="58" width="20" height="12" rx="4" fill="{TEAL}"/>
  <rect x="56" y="58" width="16" height="12" rx="4" fill="{CORAL}"/>
  <rect x="72" y="58" width="20" height="12" rx="4" fill="{TEAL}"/>
  <path d="M92 70 Q120 96 124 124" stroke="{TEAL}" stroke-width="4" fill="none"
        stroke-dasharray="2 7" stroke-linecap="round"/>
  <!-- integrated barcode -->
  <rect x="150" y="184" width="7" height="32" rx="2" fill="{TEAL}"/>
  <rect x="161" y="184" width="4" height="32" rx="2" fill="{WHITE}"/>
  <rect x="169" y="184" width="9" height="32" rx="2" fill="{AMBER}"/>
""")


# ---------------------------------------------------------------- concept C
# "Scissors + donor": explicit scissors cutting the strand; donor arc above with
# homology arms; edit highlighted; barcode integrated.
CONCEPT_C = svg(f"""
  <rect x="30" y="156" width="196" height="8" rx="4" fill="{WHITE}" opacity="0.85"/>
  <!-- scissors: two blades + pivot -->
  <circle cx="128" cy="128" r="7" fill="{AMBER}"/>
  <path d="M128 128 L96 104"  stroke="{AMBER}" stroke-width="7" stroke-linecap="round"/>
  <path d="M128 128 L160 104" stroke="{AMBER}" stroke-width="7" stroke-linecap="round"/>
  <circle cx="92"  cy="100" r="9" fill="none" stroke="{AMBER}" stroke-width="6"/>
  <circle cx="164" cy="100" r="9" fill="none" stroke="{AMBER}" stroke-width="6"/>
  <path d="M128 135 L118 160" stroke="{AMBER}" stroke-width="7" stroke-linecap="round"/>
  <path d="M128 135 L138 160" stroke="{AMBER}" stroke-width="7" stroke-linecap="round"/>
  <!-- donor arc recruited above, homology arms + edit -->
  <path d="M70 70 Q128 40 186 70" stroke="{TEAL}" stroke-width="10" fill="none" stroke-linecap="round"/>
  <circle cx="128" cy="49" r="9" fill="{CORAL}"/>
  <!-- integrated barcode -->
  <rect x="112" y="186" width="7"  height="32" rx="2" fill="{TEAL}"/>
  <rect x="123" y="186" width="4"  height="32" rx="2" fill="{WHITE}"/>
  <rect x="131" y="186" width="9"  height="32" rx="2" fill="{AMBER}"/>
""")


# ---------------------------------------------------------------- concept D
# "Recruitment loop": the donor is pulled along a bold curved arrow to the cut
# gap; strong sense of recruitment; edit + barcode at the repaired locus.
CONCEPT_D = svg(f"""
  <!-- duplex with DSB gap -->
  <rect x="34"  y="158" width="74" height="8" rx="4" fill="{WHITE}" opacity="0.9"/>
  <rect x="148" y="158" width="74" height="8" rx="4" fill="{WHITE}" opacity="0.9"/>
  <!-- donor (homology arms + edit) at upper left -->
  <rect x="40" y="64" width="20" height="13" rx="4" fill="{TEAL}"/>
  <rect x="60" y="64" width="16" height="13" rx="4" fill="{CORAL}"/>
  <rect x="76" y="64" width="20" height="13" rx="4" fill="{TEAL}"/>
  <!-- bold recruitment arrow sweeping donor to the cut -->
  <path d="M100 74 C168 78 196 120 142 150" stroke="{AMBER}" stroke-width="9"
        fill="none" stroke-linecap="round"/>
  <path d="M150 140 L142 152 L130 145 Z" fill="{AMBER}"/>
  <!-- integrated barcode at repaired site -->
  <rect x="112" y="184" width="7"  height="34" rx="2" fill="{TEAL}"/>
  <rect x="123" y="184" width="4"  height="34" rx="2" fill="{WHITE}"/>
  <rect x="131" y="184" width="9"  height="34" rx="2" fill="{AMBER}"/>
""")


CONCEPTS = {
    "A_donor_bridge": ("A — Donor bridges the cut", CONCEPT_A),
    "B_cas9_wedge": ("B — Cas9 wedge + tethered donor", CONCEPT_B),
    "C_scissors_donor": ("C — Scissors + donor arc", CONCEPT_C),
    "D_recruitment_arrow": ("D — Donor recruitment sweep", CONCEPT_D),
}


def render(svg_str: str, px: int) -> Image.Image:
    png = cairosvg.svg2png(bytestring=svg_str.encode(), output_width=px, output_height=px)
    return Image.open(io.BytesIO(png)).convert("RGBA")


def main() -> None:
    size = 300
    pad = 28
    label_h = 40
    imgs = []
    for key, (label, s) in CONCEPTS.items():
        img = render(s, size)
        (HERE / f"concept_{key}.png").write_bytes(
            cairosvg.svg2png(bytestring=s.encode(), output_width=512, output_height=512)
        )
        imgs.append((label, img))

    n = len(imgs)
    W = n * size + (n + 1) * pad
    H = size + 2 * pad + label_h
    sheet = Image.new("RGBA", (W, H), "#f5f3fb")
    draw = ImageDraw.Draw(sheet)
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 20)
    except OSError:
        font = ImageFont.load_default()
    for i, (label, img) in enumerate(imgs):
        x = pad + i * (size + pad)
        sheet.paste(img, (x, pad), img)
        tb = draw.textbbox((0, 0), label, font=font)
        tw = tb[2] - tb[0]
        draw.text((x + (size - tw) / 2, pad + size + 8), label, fill="#2a2342", font=font)
    out = HERE / "concepts_montage.png"
    sheet.convert("RGB").save(out)
    print("wrote", out)
    for key in CONCEPTS:
        print("wrote", HERE / f"concept_{key}.png")


if __name__ == "__main__":
    main()
