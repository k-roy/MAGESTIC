#!/usr/bin/env python3
"""Generate MAGESTIC branding assets from the canonical SVG mark.

Outputs (all written next to this script in docs/branding/):
    magestic_favicon.png        64x64    site favicon
    magestic_apple_touch.png    180x180  apple-touch icon
    magestic_banner_light.png   wordmark banner for light backgrounds
    magestic_banner_dark.png    wordmark banner for dark backgrounds
    magestic_social_card.png    1280x640 GitHub Open-Graph / social-preview card

The mark itself is hand-authored vector art in ``magestic_mark.svg`` (edit that
file to change the logo); this script only rasterizes it and composes the
banners. Re-run after editing the SVG:

    python docs/branding/generate_branding.py

Requires: cairosvg, matplotlib, Pillow.
"""
from __future__ import annotations

import io
from pathlib import Path

import cairosvg
import matplotlib.pyplot as plt
from matplotlib import font_manager
from PIL import Image

HERE = Path(__file__).resolve().parent
MARK_SVG = HERE / "magestic_mark.svg"

# Indigo-Cyan brand palette (matches magestic_mark.svg)
CYAN = "#00c6ff"
BLUE = "#0072ff"
INDIGO = "#101046"
ARM = "#5ff0e6"          # mint homology-arm / barcode
EDIT = "#ff4d8d"         # central edit accent
ACCENT = EDIT            # underbar / rule color
INK = "#0a0f24"          # deep indigo ink for light-bg text
PAPER = "#ffffff"
CARD_TOP = "#0b1f3a"     # social-card gradient
CARD_BOT = "#06101f"


def _render_mark(px: int) -> Image.Image:
    """Rasterize the SVG mark to a square RGBA PIL image of side ``px``."""
    png_bytes = cairosvg.svg2png(
        url=str(MARK_SVG), output_width=px, output_height=px
    )
    return Image.open(io.BytesIO(png_bytes)).convert("RGBA")


def write_icons() -> None:
    _render_mark(64).save(HERE / "magestic_favicon.png")
    _render_mark(180).save(HERE / "magestic_apple_touch.png")


def write_banner(dark: bool) -> None:
    bg = INK if dark else PAPER
    word_color = PAPER if dark else INK
    tag_color = "#9fc6e0" if dark else "#52617a"

    fig = plt.figure(figsize=(8.2, 2.0), dpi=200)
    fig.patch.set_facecolor(bg)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # mark on the left
    mark = _render_mark(360)
    ax.imshow(mark, extent=(0.03, 0.225, 0.18, 0.92), aspect="auto", zorder=3)

    # wordmark
    ax.text(
        0.265, 0.60, "MAGESTIC",
        fontsize=46, fontweight="bold", color=word_color,
        family="DejaVu Sans", va="center", ha="left",
        zorder=3,
    )
    # accent underbar
    ax.add_patch(plt.Rectangle((0.267, 0.34), 0.40, 0.030, color=ACCENT, zorder=3))
    # tagline
    ax.text(
        0.268, 0.20,
        "Multiplexed Accurate Genome Editing with Short, Trackable,\n"
        "Integrated Cellular barcodes",
        fontsize=10.5, color=tag_color, family="DejaVu Sans",
        va="center", ha="left", linespacing=1.25, zorder=3,
    )

    suffix = "dark" if dark else "light"
    fig.savefig(
        HERE / f"magestic_banner_{suffix}.png",
        facecolor=bg, bbox_inches="tight", pad_inches=0.12,
    )
    plt.close(fig)


def write_social_card() -> None:
    """1280x640 Open-Graph card GitHub shows when the repo link is shared."""
    from PIL import ImageDraw, ImageFont, ImageFilter

    W, H = 1280, 640
    # vertical indigo gradient background
    col = Image.new("RGB", (1, H))
    tr, tg, tb = Image.new("RGB", (1, 1), CARD_TOP).getpixel((0, 0))
    br, bg_, bb = Image.new("RGB", (1, 1), CARD_BOT).getpixel((0, 0))
    for y in range(H):
        t = y / (H - 1)
        col.putpixel((0, y), (int(tr + (br - tr) * t), int(tg + (bg_ - tg) * t),
                              int(tb + (bb - tb) * t)))
    card = col.resize((W, H))

    glow = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    ImageDraw.Draw(glow).ellipse([110, 150, 490, 530], fill=EDIT + "55")
    card.paste(glow.filter(ImageFilter.GaussianBlur(60)),
               (0, 0), glow.filter(ImageFilter.GaussianBlur(60)))

    mark = _render_mark(360)
    card.paste(mark, (130, 140), mark)

    d = ImageDraw.Draw(card)

    def font(sz, bold=True):
        name = "Arial Bold.ttf" if bold else "Arial.ttf"
        try:
            return ImageFont.truetype(f"/System/Library/Fonts/Supplemental/{name}", sz)
        except OSError:
            return ImageFont.load_default()

    tx = 560
    d.text((tx, 210), "MAGESTIC", font=font(86), fill="#ffffff")
    d.rectangle([tx + 4, 318, tx + 150, 326], fill=ACCENT)
    tagline = "Barcoded CRISPR editing & analysis for yeast functional genomics"
    words, lines, cur = tagline.split(), [], ""
    f = font(30, bold=False)
    for w in words:
        test = (cur + " " + w).strip()
        if d.textbbox((0, 0), test, font=f)[2] > 620:
            lines.append(cur)
            cur = w
        else:
            cur = test
    lines.append(cur)
    yy = 350
    for ln in lines:
        d.text((tx, yy), ln, font=f, fill="#aebfd6")
        yy += 40
    d.text((tx, yy + 18), "github.com/k-roy/MAGESTIC", font=font(22, bold=False), fill=ACCENT)
    card.save(HERE / "magestic_social_card.png")


def main() -> None:
    assert MARK_SVG.exists(), f"missing {MARK_SVG}"
    write_icons()
    write_banner(dark=False)
    write_banner(dark=True)
    write_social_card()
    for f in sorted(HERE.glob("magestic_*.png")):
        print("wrote", f.relative_to(HERE.parent.parent))


if __name__ == "__main__":
    # font_manager import kept for environments that register custom fonts;
    # the banners fall back to DejaVu Sans which ships with matplotlib.
    _ = font_manager
    main()
