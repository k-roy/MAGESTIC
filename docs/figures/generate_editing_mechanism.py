#!/usr/bin/env python3
"""Editing figure: the full MAGESTIC 3.0 editing mechanism (manuscript-grade).

Mirrors Kevin's own 3-panel schematic (MAGESTIC_3_and_REDI.pdf + his hand-drawn
reference): three routes deliver the donor / repair template to the Cas9 break,
converging on MAGESTIC 3.0:

  1. Plasmid recruitment  — the dsDNA plasmid donor, via LexA-DBD binding the
     LexO operator array.
  2. Retron recruitment   — the retron-derived ssDNA donor (MS2 stem-loops),
     via MCP binding the MS2 loops.
  3. Plasmid assembly     — immediate repair enhancement from a co-transformed
     linearized backbone + marker cassette assembling at the break.

Routes 1 + 2 share a SINGLE MCP-LexA(DBD)-FHA fusion protein: MCP binds MS2,
LexA binds LexO, and FHA binds pT (phospho-threonine) at the break — the engaged
adapter differs per route, FHA always docks the break. Iconic style (simple
shapes, no molecular realism) per Kevin.

NB pT = phospho-threonine (verbatim, the FHA-binding mark). Panel 3 geometry is
the least-resolved part of the reference — confirm with Kevin.

Run: python docs/figures/generate_editing_mechanism.py   (needs cairosvg)
"""
from __future__ import annotations

import math
from pathlib import Path

import magfig as M

HERE = Path(__file__).resolve().parent
W, H = 1200, 588
COLS = (260, 600, 940)
BRK_Y = 360                      # break top-strand y (shared across panels)


def duplex(x0, x1, y, color, sw=3, gap=12):
    return (M.line(x0, y, x1, y, color, sw=sw) + "\n" +
            M.line(x0, y + gap, x1, y + gap, color, sw=sw))


def jagged_cut(cx, y, pal, sw=2.6):
    return (f'<path d="M{cx-6} {y-12} L{cx+3} {y-1} L{cx-4} {y+6} '
            f'L{cx+6} {y+17}" stroke="{pal["amber"]}" stroke-width="{sw}" '
            f'fill="none" stroke-linejoin="round" stroke-linecap="round"/>')


def lobe(cx, cy, label, color, fill, pal, size=10):
    return M.protein(cx, cy, 14 + 7.5 * len(label), 21, label, color, fill, pal,
                     size=size)


def break_unit(cx, pal, y=BRK_Y, w=150):
    """Duplex with a jagged cut flanked by two pT (phospho-threonine) ovals."""
    half = w / 2
    out = []
    for px in (cx - 44, cx + 44):
        out.append(f'<ellipse cx="{px}" cy="{y+6}" rx="22" ry="15" '
                   f'fill="{pal["edit_l"]}" stroke="{pal["edit"]}" '
                   f'stroke-width="1.2"/>')
        out.append(M.text(px, y + 9.5, "pT", size=9, weight="700",
                          color=pal["edit"], anchor="middle"))
    out.append(duplex(cx - half, cx - 12, y, pal["dna"]))
    out.append(duplex(cx + 12, cx + half, y, pal["dna"]))
    out.append(jagged_cut(cx, y + 6, pal))
    return "\n".join(out)


def fusion(cx, pal, mode):
    """The single MCP-LexA(DBD)-FHA fusion protein, engaged per recruitment mode.

    mode='lexa' (plasmid route): LexA-DBD on top (binds LexO above), MCP to side.
    mode='mcp'  (retron route):  MCP on top (binds MS2 above), LexA to side.
    FHA is always the lower lobe, tethered to pT at the break.
    """
    top_y, fha_y = 296, 330
    out = []
    if mode == "lexa":
        out.append(lobe(cx, top_y, "LexA", pal["cyan"], pal["cyan_l"], pal))
        out.append(lobe(cx - 40, top_y, "MCP", pal["amber"], pal["amber_l"], pal,
                        size=9))
        out.append(M.line(cx - 22, top_y, cx - 14, top_y, pal["slate"], sw=2))
    else:
        out.append(lobe(cx, top_y, "MCP", pal["amber"], pal["amber_l"], pal))
        out.append(lobe(cx - 40, top_y, "LexA", pal["cyan"], pal["cyan_l"], pal,
                        size=9))
        out.append(M.line(cx - 22, top_y, cx - 14, top_y, pal["slate"], sw=2))
    out.append(lobe(cx, fha_y, "FHA", pal["blue"], pal["blue_l"], pal))
    out.append(M.line(cx, top_y + 11, cx, fha_y - 11, pal["slate"], sw=2))
    # FHA -> pT at the break
    out.append(M.line(cx, fha_y + 11, cx, BRK_Y + 4, pal["blue"], sw=1.6,
                      dash="2 4"))
    return "\n".join(out)


def header(c, cx, title, sub, pal):
    c.add(M.text(cx, 100, title, size=13, weight="700", color=pal["heading"],
                 anchor="middle"))
    c.add(M.text(cx, 117, sub, size=9, italic=True, color=pal["muted"],
                 anchor="middle"))


def build(pal: dict) -> M.Canvas:
    c = M.Canvas(W, H, pal)

    c.add(M.text(40, 36, "The MAGESTIC 3.0 editing mechanism",
                 size=19, weight="700", color=pal["title"]))
    c.add(M.text(40, 59, "three routes deliver the donor to the Cas9 break for "
                 "maximal HDR", size=12.5, color=pal["label"]))

    # ── Panel 1 · Plasmid recruitment (dsDNA donor, LexA-DBD : LexO) ──────────
    cx = COLS[0]
    header(c, cx, "Plasmid recruitment", "dsDNA plasmid donor", pal)
    pcx, pcy, pr = cx, 198, 52
    c.add(f'<circle cx="{pcx}" cy="{pcy}" r="{pr}" fill="none" '
          f'stroke="{pal["divider"]}" stroke-width="7"/>')
    c.add(M.arc(pcx, pcy, pr, 150, 210, pal["cyan"], 7, cap="round"))   # guide
    c.add(M.arc(pcx, pcy, pr, 330, 26, pal["blue"], 7, cap="round"))    # donor
    ea = math.radians(2)
    c.add(f'<circle cx="{pcx+pr*math.cos(ea):.0f}" cy="{pcy+pr*math.sin(ea):.0f}" '
          f'r="4" fill="{pal["edit"]}"/>')
    c.add(M.text(pcx - pr - 6, pcy - 2, "guide", size=9, color=pal["cyan"],
                 weight="600", anchor="end"))
    c.add(M.text(pcx + pr + 6, pcy + 6, "donor*", size=9, color=pal["blue"],
                 weight="600", anchor="start"))
    for a in (80, 92, 104):
        ar = math.radians(a)
        c.add(f'<circle cx="{pcx+pr*math.cos(ar):.0f}" '
              f'cy="{pcy+pr*math.sin(ar):.0f}" r="4" fill="{pal["teal"]}"/>')
    c.add(M.text(pcx + 40, pcy + pr + 4, "LexO", size=8.5, color=pal["teal"],
                 anchor="middle"))
    c.add(M.line(pcx, pcy + pr, pcx, 285, pal["teal"], sw=1.6))
    c.add(fusion(cx, pal, "lexa"))
    c.add(break_unit(cx, pal))

    # ── Panel 2 · Retron recruitment (ssDNA donor, MCP : MS2 loops) ───────────
    cx = COLS[1]
    header(c, cx, "Retron recruitment", "retron-derived ssDNA donor", pal)
    ry = 206
    c.add(M.line(cx - 78, ry, cx + 30, ry, pal["teal"], sw=2.6))         # transcript
    for i in range(3):
        c.add(M.stem_loop(cx - 64 + i * 18, ry, 22, pal["teal"], sw=2.2, w=8))
    # msd loop (terminal) + edit dot
    c.add(f'<path d="M{cx+30} {ry} C{cx+58} {ry-26}, {cx+92} {ry+18}, '
          f'{cx+62} {ry+30} C{cx+40} {ry+38}, {cx+28} {ry+12}, {cx+30} {ry}" '
          f'fill="none" stroke="{pal["blue"]}" stroke-width="2.4"/>')
    c.add(f'<circle cx="{cx-2}" cy="{ry}" r="4" fill="{pal["edit"]}"/>')
    c.add(M.text(cx - 50, ry + 24, "MS2 loops", size=8.5, color=pal["teal"],
                 anchor="middle"))
    c.add(M.line(cx, ry + 6, cx, 285, pal["teal"], sw=1.6))
    c.add(fusion(cx, pal, "mcp"))
    c.add(break_unit(cx, pal))

    # ── Panel 3 · Plasmid assembly (PROPOSED; mechanism not fully understood) ──
    #    The one linear donor is a dual-purpose template: it templates the
    #    genomic edit AND templates DNA synthesis from the marker fragment to
    #    complete a new circular plasmid — without consuming the backbone into
    #    the circle. (Prose in the footer; panel stays iconic.)
    cx = COLS[2]
    header(c, cx, "Plasmid assembly", "(proposed)", pal)
    dy = 158
    c.add(M.line(cx - 84, dy, cx + 18, dy, pal["dna"], sw=4))
    c.add(M.line(cx - 36, dy, cx + 14, dy, pal["blue"], sw=4))           # donor region
    c.add(f'<circle cx="{cx-2}" cy="{dy}" r="4" fill="{pal["edit"]}"/>')
    c.add(M.text(cx - 30, dy - 12, "linear donor", size=8.5,
                 color=pal["muted"], anchor="middle"))
    # (a) templates the genomic edit → down into the break
    c.add(M.arrow(cx - 2, dy + 8, cx, 300, pal["label"], sw=1.8, dash="3 4"))
    c.add(M.text(cx - 24, 244, "templates", size=8, italic=True,
                 color=pal["slate"], anchor="end"))
    c.add(M.text(cx - 24, 255, "the edit", size=8, italic=True,
                 color=pal["slate"], anchor="end"))
    # (b) templates synthesis from the marker fragment → a new circular plasmid
    ccx, ccy, cr = cx + 60, 216, 28
    c.add(M.arc(ccx, ccy, cr, 40, 240, pal["teal"], 5))                  # marker (existing)
    a0, a1 = math.radians(240), math.radians(40)
    x0, y0 = ccx + cr * math.cos(a0), ccy + cr * math.sin(a0)
    x1, y1 = ccx + cr * math.cos(a1), ccy + cr * math.sin(a1)
    c.add(f'<path d="M{x0:.1f} {y0:.1f} A{cr} {cr} 0 1 1 {x1:.1f} {y1:.1f}" '
          f'stroke="{pal["blue"]}" stroke-width="2.4" fill="none" '
          f'stroke-dasharray="3 3"/>')                                   # new synthesis
    c.add(M.text(ccx, ccy + 3, "marker", size=8, color=pal["teal"],
                 anchor="middle"))
    c.add(M.text(ccx, ccy + cr + 13, "new plasmid", size=8, color=pal["muted"],
                 anchor="middle"))
    c.add(M.arrow(cx + 20, dy + 2, ccx - cr - 3, ccy - 8, pal["label"], sw=1.6,
                  dash="3 4"))
    c.add(M.text(cx + 34, dy + 16, "templates", size=8, italic=True,
                 color=pal["slate"], anchor="start"))
    c.add(M.text(cx + 34, dy + 27, "synthesis", size=8, italic=True,
                 color=pal["slate"], anchor="start"))
    c.add(break_unit(cx, pal))

    # ── convergence → MAGESTIC 3.0 ────────────────────────────────────────────
    conv = (600, 432)
    for cx in COLS:
        c.add(M.arrow(cx, BRK_Y + 28, conv[0] + (cx - 600) * 0.18, conv[1] - 6,
                      pal["label"], sw=2))
    c.add(M.text(conv[0], conv[1] + 22, "MAGESTIC 3.0", size=17, weight="700",
                 color=pal["title"], anchor="middle"))
    # edited product
    ey = conv[1] + 44
    c.add(duplex(conv[0] - 64, conv[0] + 64, ey, pal["dna"]))
    c.add(f'<circle cx="{conv[0]}" cy="{ey+6}" r="5.5" fill="{pal["edit"]}"/>')
    c.add(M.text(conv[0], ey + 30, "maximal HDR  ·  precise edit", size=10.5,
                 italic=True, color=pal["muted"], anchor="middle"))

    # ── legend (one line each) ────────────────────────────────────────────────
    c.add(M.text(40, H - 46, "MCP–LexA(DBD)–FHA — a single fusion protein:  "
                 "MCP binds MS2 loops  ·  LexA binds the LexO array  ·  "
                 "FHA binds pT at the break", size=9, color=pal["label"],
                 anchor="start"))
    c.add(M.text(40, H - 32, "pT = phospho-threonine (the FHA-binding mark at the "
                 "Cas9 break)", size=9, italic=True, color=pal["muted"],
                 anchor="start"))
    c.add(M.text(40, H - 18, "Plasmid assembly (proposed, mechanism not fully "
                 "resolved): the one linear donor templates both the genomic edit "
                 "and synthesis from the marker fragment into a new circular "
                 "plasmid — without consuming the backbone.   ·   Roy et al., in "
                 "preparation", size=9, italic=True, color=pal["muted"],
                 anchor="start"))
    return c


if __name__ == "__main__":
    M.render(build, "editing_mechanism", HERE)
