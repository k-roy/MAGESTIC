#!/usr/bin/env python3
"""Barcoding figure: ~100% BC1 integration at the dedicated barcode locus.

Simplified / iconic. Distinct from editing (which happens at the target locus):
the barcode is integrated at a separate, dedicated barcode locus. Integration is
forced to ~100% because two orthogonal nucleases — I-SceI and SaCas9 (guide X) —
cut BOTH the plasmid and the barcode locus at the same engineered site; SaCas9
also removes the plasmid afterwards.

Engineered cut site: gTatccTAGGGATAACAGGGTAATgagtAC
(contains the 18-bp I-SceI site TAGGGATAACAGGGTAAT).

Run: python docs/figures/generate_barcode_integration.py   (needs cairosvg)
"""
from __future__ import annotations

import math
from pathlib import Path

import magfig as M

HERE = Path(__file__).resolve().parent
W, H = 1060, 430
CITE = "Roy et al., bioRxiv (2024)"

SEQ = "gTatccTAGGGATAACAGGGTAATgagtAC"
ISCE_START, ISCE_LEN = 6, 18          # 18-bp I-SceI site within SEQ
CHARW = 6.6                           # monospace advance at size 11 (uniform)


def duplex(x0, x1, y, color, sw=3, gap=14):
    return (M.line(x0, y, x1, y, color, sw=sw) + "\n" +
            M.line(x0, y + gap, x1, y + gap, color, sw=sw))


def cutmark(x, y, pal):
    return (f'<path d="M{x-7} {y-12} L{x+3} {y} L{x-5} {y+8} L{x+7} {y+20}" '
            f'stroke="{pal["amber"]}" stroke-width="2.6" fill="none" '
            f'stroke-linejoin="round" stroke-linecap="round"/>')


def tag(cx, cy, label, color, fill, pal, size=9.5):
    return M.protein(cx, cy, 14 + 7 * len(label), 19, label, color, fill, pal,
                     size=size)


def build(pal: dict) -> M.Canvas:
    c = M.Canvas(W, H, pal)

    c.add(M.text(40, 34, "100% barcode integration at the barcode locus",
                 size=18, weight="700", color=pal["title"]))
    c.add(M.text(40, 56, "I-SceI + SaCas9 cut both the plasmid and the barcode "
                 "locus — forcing BC1 integration and removing the plasmid",
                 size=12.5, color=pal["label"]))

    # ── 1 · barcoded plasmid (barcode flanked by cut sites)
    pcx, pcy, pr = 150, 252, 64
    c.add(f'<circle cx="{pcx}" cy="{pcy}" r="{pr}" fill="none" '
          f'stroke="{pal["divider"]}" stroke-width="8"/>')
    # barcode arc (striped) across the top
    a0, a1 = 205, 335
    c.add(M.arc(pcx, pcy, pr, a0, a1, pal["amber"], 8, cap="butt"))
    for i in range(7):
        a = math.radians(a0 + (a1 - a0) * i / 6)
        x0 = pcx + (pr - 4) * math.cos(a); y0 = pcy + (pr - 4) * math.sin(a)
        x1 = pcx + (pr + 4) * math.cos(a); y1 = pcy + (pr + 4) * math.sin(a)
        c.add(M.line(x0, y0, x1, y1, pal["bg"] if i % 2 else pal["amber"], sw=1.6))
    c.add(M.text(pcx, pcy - pr - 12, "barcode", size=10.5, color=pal["amber"],
                 weight="700", anchor="middle"))
    # cut marks flanking the barcode
    for a, nm, col, fill in ((200, "I-SceI", pal["edit"], pal["edit_l"]),
                             (340, "SaCas9", pal["cyan"], pal["cyan_l"])):
        ar = math.radians(a)
        cx = pcx + pr * math.cos(ar); cy = pcy + pr * math.sin(ar)
        c.add(cutmark(cx, cy - 4, pal))
        c.add(tag(cx + (28 if math.cos(ar) > 0 else -28), cy + 18, nm, col, fill, pal))
    c.add(M.text(pcx, pcy + pr + 30, "barcoded plasmid", size=11.5, weight="700",
                 color=pal["heading"], anchor="middle"))

    c.add(M.arrow(pcx + pr + 12, pcy, 330, pcy, pal["label"], sw=2))
    c.add(M.text(300, pcy - 10, "cut both", size=9.5,
                 italic=True, color=pal["muted"], anchor="middle"))

    # ── 2 · barcode locus cut at the same engineered site
    gx0, gx1, cutx, gy = 372, 720, 548, 250
    c.add(duplex(gx0, cutx - 12, gy, pal["dna"]))
    c.add(duplex(cutx + 12, gx1, gy, pal["dna"]))
    c.add(cutmark(cutx, gy, pal))
    c.add(tag(cutx - 34, gy - 26, "I-SceI", pal["edit"], pal["edit_l"], pal))
    c.add(tag(cutx + 36, gy - 26, "SaCas9", pal["cyan"], pal["cyan_l"], pal))
    c.add(M.text(cutx, gy + 36, "barcode locus", size=10.5, color=pal["muted"],
                 anchor="middle"))
    # engineered sequence (mono, uniform) with the I-SceI site underlined
    width = len(SEQ) * CHARW
    start_x = cutx - width / 2
    sy = gy + 58
    c.add(M.text(cutx, sy, SEQ, size=11, family=M.MONO, color=pal["heading"],
                 anchor="middle"))
    ux0 = start_x + ISCE_START * CHARW
    ux1 = ux0 + ISCE_LEN * CHARW
    c.add(M.line(ux0, sy + 5, ux1, sy + 5, pal["amber"], sw=2.4))
    c.add(M.text((ux0 + ux1) / 2, sy + 17, "18-bp I-SceI site (SaCas9 guide X overlaps)",
                 size=8.5, italic=True, color=pal["amber"], anchor="middle"))

    c.add(M.arrow(gx1 + 8, gy + 7, 812, gy + 7, pal["label"], sw=2.2))
    c.add(M.text((gx1 + 8 + 812) / 2, gy - 4, "integrate", size=11, weight="700",
                 color=pal["heading"], anchor="middle"))

    # ── 3 · BC1 integrated (~100%) at the barcode locus
    ex0, ex1, ecx = 818, 1030, 924
    c.add(duplex(ex0, ecx - 26, gy, pal["dna"]))
    c.add(duplex(ecx + 26, ex1, gy, pal["dna"]))
    bc, bw = M.barcode(ecx - 24, gy - 5, 24, pal,
                       widths=(4, 3, 5, 3, 4, 6, 3),
                       colors=[pal["amber"]] * 7, gap=2.5)
    c.add(bc)
    c.add(M.text(ecx, gy - 22, "BC1 integrated", size=11.5, color=pal["amber"],
                 weight="700", anchor="middle"))
    c.add(M.text(ecx, gy + 40, "~100% efficiency", size=10.5, italic=True,
                 color=pal["muted"], anchor="middle"))

    # plasmid-removal note (bottom-left of stage 3, clears the right edge)
    rcx = ex0 + 14
    c.add(f'<circle cx="{rcx}" cy="352" r="15" fill="none" stroke="{pal["divider"]}" '
          f'stroke-width="4"/>')
    c.add(M.line(rcx - 9, 343, rcx + 9, 361, pal["edit"], sw=2.4))
    c.add(M.line(rcx + 9, 343, rcx - 9, 361, pal["edit"], sw=2.4))
    c.add(M.text(rcx + 24, 356, "plasmid removed (SaCas9)", size=9.5,
                 color=pal["muted"], anchor="start"))

    c.add(M.text(40, H - 12, CITE, size=8.5, color=pal["muted"], anchor="start",
                 italic=True))
    return c


if __name__ == "__main__":
    M.render(build, "barcode_integration", HERE)
