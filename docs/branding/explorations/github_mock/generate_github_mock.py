#!/usr/bin/env python3
"""Mock up the real GitHub touchpoints for the branding finalists.

For each chosen mark, render one composite sheet showing it exactly where it
will live on GitHub:

  1. AVATAR / repo-list row  — the org avatar at 20/40/96 px, plus a simulated
     "k-roy/MAGESTIC" repo row on GitHub's dark UI (#0d1117).  This is the size
     where a busy mark falls apart, so it's the real stress test.
  2. SOCIAL PREVIEW CARD      — the 1280x640 Open-Graph card GitHub shows when
     the repo link is shared (Slack/Twitter/etc).  Mark + wordmark + tagline.
  3. README BANNER            — the wide header image at the top of README.md.

Run: python generate_github_mock.py   (needs cairosvg, Pillow)
"""
from __future__ import annotations

import io
from pathlib import Path

import cairosvg
from PIL import Image, ImageDraw, ImageFont, ImageFilter

ROOT = Path(__file__).resolve().parent
EXPL = ROOT.parent

TAGLINE = "Barcoded CRISPR editing & analysis for yeast functional genomics"
REPO = "k-roy / MAGESTIC"
DESC = "Multiplexed Accurate Genome Editing with Short, Trackable, Integrated Cellular barcodes"

# theme per finalist: svg path, display name, banner background (dark/light),
# title color, tagline color, accent (rule/underline), lang-dot color.
FINALISTS = {
    "donor_cyan": {
        "svg": EXPL / "palettes/cyan_vibrant.svg",
        "name": "Donor-Bridge — Indigo-Cyan",
        "bg": ["#0b1f3a", "#06101f"],
        "title": "#ffffff",
        "tag": "#aebfd6",
        "accent": "#ff4d8d",
        "light": False,
    },
    "donor_tidewave": {
        "svg": EXPL / "palettes/tidewave_vibrant.svg",
        "name": "Donor-Bridge — Tidewave",
        "bg": ["#0a2740", "#06121f"],
        "title": "#ffffff",
        "tag": "#a9c6da",
        "accent": "#ff5da2",
        "light": False,
    },
    "cas9_slate": {
        "svg": EXPL / "cas9_palettes/slate_steel.svg",
        "name": "Cas9 — Slate / Steel-Blue",
        "bg": ["#141f2b", "#0a121b"],
        "title": "#ffffff",
        "tag": "#9fb2c6",
        "accent": "#ffb454",
        "light": False,
    },
    "cas9_clinical": {
        "svg": EXPL / "cas9_palettes/clinical_light.svg",
        "name": "Cas9 — Clinical (light)",
        "bg": ["#f5f7fb", "#e6ebf3"],
        "title": "#1b2a44",
        "tag": "#52617a",
        "accent": "#e2502a",
        "light": True,
    },
}

GH_BG = "#0d1117"      # GitHub dark canvas
GH_CARD = "#161b22"    # GitHub panel
GH_BLUE = "#2f81f7"    # GitHub link blue
GH_TEXT = "#c9d1d9"
GH_MUTED = "#8b949e"
GH_BORDER = "#30363d"


def font(sz, bold=True, mono=False):
    paths = []
    if mono:
        paths = ["/System/Library/Fonts/Menlo.ttc"]
    elif bold:
        paths = [
            "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
            "/System/Library/Fonts/HelveticaNeue.ttc",
        ]
    else:
        paths = [
            "/System/Library/Fonts/Supplemental/Arial.ttf",
            "/System/Library/Fonts/HelveticaNeue.ttc",
        ]
    for p in paths:
        try:
            return ImageFont.truetype(p, sz)
        except OSError:
            continue
    return ImageFont.load_default()


def render_mark(svg_path: Path, px: int) -> Image.Image:
    png = cairosvg.svg2png(bytestring=svg_path.read_bytes(), output_width=px, output_height=px)
    return Image.open(io.BytesIO(png)).convert("RGBA")


def rounded(img: Image.Image, radius: int) -> Image.Image:
    mask = Image.new("L", img.size, 0)
    ImageDraw.Draw(mask).rounded_rectangle([0, 0, img.size[0], img.size[1]], radius, fill=255)
    out = img.copy()
    out.putalpha(mask)
    return out


def circle(img: Image.Image) -> Image.Image:
    mask = Image.new("L", img.size, 0)
    ImageDraw.Draw(mask).ellipse([0, 0, img.size[0], img.size[1]], fill=255)
    out = img.copy()
    out.putalpha(mask)
    return out


def vgrad(size, top, bottom):
    w, h = size
    base = Image.new("RGB", (1, h))
    tr, tg, tb = Image.new("RGB", (1, 1), top).getpixel((0, 0))
    br, bg, bb = Image.new("RGB", (1, 1), bottom).getpixel((0, 0))
    for y in range(h):
        t = y / max(1, h - 1)
        base.putpixel((0, y), (int(tr + (br - tr) * t), int(tg + (bg - tg) * t), int(tb + (bb - tb) * t)))
    return base.resize((w, h))


def text_center(draw, cx, y, s, fnt, fill):
    tb = draw.textbbox((0, 0), s, font=fnt)
    draw.text((cx - (tb[2] - tb[0]) / 2, y), s, font=fnt, fill=fill)
    return tb[3] - tb[1]


# ----------------------------------------------------------- panels
def avatar_row(theme) -> Image.Image:
    """GitHub dark canvas: avatar at 3 sizes + a simulated repo row."""
    W, H = 1280, 300
    im = Image.new("RGB", (W, H), GH_BG)
    d = ImageDraw.Draw(im)
    big = render_mark(theme["svg"], 384)

    # three avatar sizes (orgs = rounded square)
    d.text((40, 28), "Org avatar / repo-list icon", font=font(20), fill=GH_MUTED)
    x = 40
    for px in (96, 40, 20):
        a = rounded(big.resize((px, px), Image.LANCZOS), int(px * 0.22))
        im.paste(a, (x, 78), a)
        d.text((x, 78 + px + 6), f"{px}px", font=font(13, bold=False), fill=GH_MUTED)
        x += px + 34
    # circle variant (user avatar)
    c = circle(big.resize((96, 96), Image.LANCZOS))
    im.paste(c, (x + 10, 78), c)
    d.text((x + 10, 78 + 96 + 6), "circle", font=font(13, bold=False), fill=GH_MUTED)

    # simulated repo-list row on the right
    rx, ry, rw, rh = 470, 70, 760, 150
    d.rounded_rectangle([rx, ry, rx + rw, ry + rh], 12, fill=GH_CARD, outline=GH_BORDER, width=1)
    av = rounded(big.resize((44, 44), Image.LANCZOS), 10)
    im.paste(av, (rx + 22, ry + 22), av)
    d.text((rx + 80, ry + 20), REPO, font=font(22), fill=GH_BLUE)
    d.rounded_rectangle([rx + 80 + 250, ry + 22, rx + 80 + 250 + 56, ry + 22 + 22], 10,
                        outline=GH_BORDER, width=1)
    d.text((rx + 80 + 258, ry + 24), "Public", font=font(13, bold=False), fill=GH_MUTED)
    # wrap description
    desc = "Barcoded CRISPR editing & analysis for yeast functional genomics"
    d.text((rx + 80, ry + 56), desc, font=font(15, bold=False), fill=GH_TEXT)
    # meta line: lang dot + stars
    my = ry + 100
    d.ellipse([rx + 80, my, rx + 80 + 14, my + 14], fill=theme["accent"])
    d.text((rx + 100, my - 2), "Python", font=font(14, bold=False), fill=GH_MUTED)
    d.text((rx + 190, my - 2), "★ 0", font=font(14, bold=False), fill=GH_MUTED)
    d.text((rx + 250, my - 2), "MIT", font=font(14, bold=False), fill=GH_MUTED)
    return im


def social_card(theme) -> Image.Image:
    """1280x640 OG card."""
    W, H = 1280, 640
    im = vgrad((W, H), theme["bg"][0], theme["bg"][1])
    d = ImageDraw.Draw(im)

    mark = render_mark(theme["svg"], 360)
    # soft glow behind the mark
    glow = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    gd = ImageDraw.Draw(glow)
    gd.ellipse([110, 150, 110 + 380, 150 + 380], fill=theme["accent"] + "55")
    glow = glow.filter(ImageFilter.GaussianBlur(60))
    im.paste(glow, (0, 0), glow)
    im.paste(mark, (130, 140), mark)

    tx = 560
    d.text((tx, 210), "MAGESTIC", font=font(86), fill=theme["title"])
    d.rectangle([tx + 4, 318, tx + 150, 326], fill=theme["accent"])
    # tagline (wrap to ~2 lines)
    words = TAGLINE.split()
    lines, cur = [], ""
    fnt = font(30, bold=False)
    for w in words:
        test = (cur + " " + w).strip()
        if d.textbbox((0, 0), test, font=fnt)[2] > 620:
            lines.append(cur)
            cur = w
        else:
            cur = test
    lines.append(cur)
    yy = 350
    for ln in lines:
        d.text((tx, yy), ln, font=fnt, fill=theme["tag"])
        yy += 40
    d.text((tx, yy + 18), "github.com/k-roy/MAGESTIC", font=font(22, bold=False),
           fill=theme["accent"])
    return im


def banner(theme) -> Image.Image:
    """Wide README header, ~1280x300."""
    W, H = 1280, 300
    im = vgrad((W, H), theme["bg"][0], theme["bg"][1])
    d = ImageDraw.Draw(im)
    mark = render_mark(theme["svg"], 200)
    im.paste(mark, (60, (H - 200) // 2), mark)
    tx = 300
    d.text((tx, 96), "MAGESTIC", font=font(64), fill=theme["title"])
    d.text((tx + 4, 176), DESC[:46] + "…", font=font(22, bold=False), fill=theme["tag"])
    d.rectangle([tx + 4, 168, tx + 120, 174], fill=theme["accent"])
    return im


def composite(key, theme) -> Path:
    pad = 26
    cap_h = 40
    panels = [avatar_row(theme), social_card(theme), banner(theme)]
    captions = ["1 · Avatar & repo-list row (GitHub dark UI)",
                "2 · Social preview card (1280×640)",
                "3 · README banner"]
    W = max(p.width for p in panels) + 2 * pad
    H = pad + cap_h + sum(p.height + cap_h + pad for p in panels)
    sheet = Image.new("RGB", (W, H), "#0a0c10")
    d = ImageDraw.Draw(sheet)
    d.text((pad, 14), theme["name"], font=font(26), fill="#ffffff")
    y = pad + cap_h
    for p, cap in zip(panels, captions):
        d.text((pad, y - 26), cap, font=font(16), fill="#8b949e")
        sheet.paste(p, (pad, y))
        y += p.height + cap_h + pad
    out = ROOT / f"github_mock_{key}.png"
    sheet.save(out)
    print("wrote", out)
    return out


def main():
    for key, theme in FINALISTS.items():
        composite(key, theme)


if __name__ == "__main__":
    main()
