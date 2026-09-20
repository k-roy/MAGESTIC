#!/usr/bin/env python3
"""Hero schematic: the end-to-end MAGESTIC workflow.

design -> build -> EDIT -> barcode -> phenotype -> analyze, as one left-to-right
flow. The headline step is the editing event (centered, emphasized): Cas9 cuts
the genomic target, a retron reverse-transcriptase makes the donor as ssDNA in
vivo, and LexA-FHA (+ MCP on the retron stem-loops) recruits that donor to the
cut to template the precise HDR edit. Barcoding follows.

Mechanism per MAGESTIC 3.0 manuscript figures (Roy et al., in prep). Package
facts: 20bp guide + 129bp donor in a 200-mer, BC0 (10bp) on the plasmid, BC1
(20bp) integrated and linked via shared BC0; phenotyping pooled (BC1 NGS) or
arrayed (REDI); scored with DESeq2 + tiered hit calling.

Run: python docs/figures/generate_pipeline_overview.py   (needs cairosvg)
"""
from __future__ import annotations

from pathlib import Path

import magfig as M

HERE = Path(__file__).resolve().parent
W, H = 1040, 384


def recruit_glyph(cx, y, pal):
    """Tiny donor-recruitment motif: a cut duplex with a donor arc + edit star."""
    out = []
    # duplex with a central break
    out.append(M.line(cx - 30, y, cx - 6, y, pal["dna"], sw=3))
    out.append(M.line(cx + 6, y, cx + 30, y, pal["dna"], sw=3))
    # donor arc descending into the cut, edit star at apex
    out.append(f'<path d="M{cx-18:.0f} {y-3:.0f} Q{cx:.0f} {y-22:.0f} '
               f'{cx+18:.0f} {y-3:.0f}" stroke="{pal["teal"]}" stroke-width="2.4" '
               f'fill="none" stroke-linecap="round"/>')
    out.append(f'<circle cx="{cx}" cy="{y-16}" r="3.2" fill="{pal["edit"]}"/>')
    return "\n".join(out)


def build(pal: dict) -> M.Canvas:
    c = M.Canvas(W, H, pal)

    # ── title
    c.add(M.text(30, 30, "The MAGESTIC workflow", size=17, weight="700",
                 color=pal["title"]))
    c.add(M.text(30, 49, "precise editing by retron ssDNA donor recruitment, "
                 "then trackable barcodes — from variant design to scored hits",
                 size=11.5, color=pal["label"]))

    # ── stage geometry (6 stages, stage index 2 = the editing headline)
    n = 6
    chip_w, chip_h = 142, 56
    edit_h = 72
    left, right = 30, W - 30
    cy = 176
    span = right - left
    step = (span - chip_w) / (n - 1)
    centers = [left + chip_w / 2 + i * step for i in range(n)]

    # ── phase bands (behind chips)
    phases = [
        ("DESIGN & BUILD", 0, 1, pal["indigo"]),
        ("EDIT & BARCODE", 2, 3, pal["edit"]),
        ("PHENOTYPE & ANALYZE", 4, 5, pal["teal"]),
    ]
    band_top, band_h = 78, 152
    for name, i0, i1, col in phases:
        x0 = centers[i0] - chip_w / 2 - 10
        x1 = centers[i1] + chip_w / 2 + 10
        core = (name == "EDIT & BARCODE")
        c.add(M.rrect(x0, band_top, x1 - x0, band_h, 12,
                      pal["edit_l"] if core else pal["divider"],
                      opacity=0.5 if core else 0.55))
        c.add(M.text((x0 + x1) / 2, band_top + 18, name, size=10.5, weight="700",
                     color=col, anchor="middle", spacing="0.5"))

    # ── stages: (label, sub, color, fill)
    stages = [
        ("Guide + donor", "20 + 129 bp in 200-mer", pal["indigo"], pal["indigo_l"]),
        ("Plasmid library", "+ BC0 (10 bp)", pal["blue"], pal["blue_l"]),
        ("Edit target site", None, pal["edit"], pal["edit_l"]),          # headline
        ("Genomic barcode", "BC1 (20 bp)", pal["cyan"], pal["cyan_l"]),
        ("Phenotype", "pooled NGS / REDI", pal["teal"], pal["teal_l"]),
        ("DESeq2 + hits", "tiered calling", pal["amber"], pal["amber_l"]),
    ]
    for i, (cx, (lab, sub, col, fill)) in enumerate(zip(centers, stages)):
        if i == 2:
            continue  # editing chip drawn separately (emphasized)
        c.add(M.chip(cx, cy, chip_w, chip_h, lab, col, fill, pal,
                     size=12.5, sub=sub))

    # ── editing headline chip: taller, thick accent border, two sub-lines
    ex = centers[2]
    ey0 = cy - edit_h / 2
    c.add(M.rrect(ex - chip_w / 2, ey0, chip_w, edit_h, 10, pal["edit_l"],
                  stroke=pal["edit"], sw=3))
    c.add(M.recruit_glyph(ex, ey0 - 12, pal) if hasattr(M, "recruit_glyph")
          else recruit_glyph(ex, ey0 - 12, pal))
    c.add(M.text(ex, cy - 12, "Edit target site", size=12.5, weight="700",
                 color=pal["edit"], anchor="middle"))
    c.add(M.text(ex, cy + 4, "Cas9 cut + retron", size=9.5, color=pal["label"],
                 anchor="middle"))
    c.add(M.text(ex, cy + 16, "ssDNA donor recruited", size=9.5, color=pal["label"],
                 anchor="middle"))
    c.add(M.text(ex, cy + edit_h / 2 + 13, "in vivo (LexA-FHA)", size=8.5,
                 italic=True, color=pal["muted"], anchor="middle"))

    # ── flow arrows between chips (account for the taller editing chip)
    def edge(i, side):
        half = (edit_h if i == 2 else chip_h) / 2  # noqa: unused but documents intent
        return centers[i] + (chip_w / 2 + 3) * (1 if side == "r" else -1)

    for a in range(n - 1):
        c.add(M.arrow(edge(a, "r"), cy, edge(a + 1, "l"), cy, pal["label"], sw=2))

    # ── barcode motif above the barcode chip (the trackable signature)
    bc, _ = M.barcode(centers[3] - 14, band_top + 30, 13, pal)
    c.add(bc)

    # ── CLI tags under the relevant stages
    cli_y = cy + chip_h / 2 + 20
    clis = {
        0: "magestic-design-vcf",
        1: "magestic-gdb-pipeline",
        3: "magestic-bc1-link",
    }
    for i, cmd in clis.items():
        c.add(M.cli(centers[i], cli_y, cmd, pal, anchor="middle"))
    # screen-analysis spans the pooled phenotype + DESeq2 stages (one pipeline)
    c.add(M.cli((centers[4] + centers[5]) / 2, cli_y, "magestic-screen-analysis",
                pal, anchor="middle"))

    # ── REDI arrayed-phenotyping branch off Phenotype (clears the CLI row)
    redi_cy = 322
    rcx = centers[4]
    c.add(M.arrow(rcx, cli_y + 12, rcx, redi_cy - 16, pal["muted"],
                  sw=1.8, dash="3 4"))
    c.add(M.chip(rcx, redi_cy, 196, 32,
                 "REDI: arrayed, indexed", pal["amber"], pal["amber_l"], pal,
                 size=11, r=9))
    c.add(M.cli(rcx, redi_cy + 28, "magestic-redi", pal, anchor="middle"))
    c.add(M.text(rcx, redi_cy + 42, "Cre lox66 × lox71 to known positions",
                 size=8.5, italic=True, color=pal["muted"], anchor="middle"))

    return c


if __name__ == "__main__":
    M.render(build, "pipeline_overview", HERE)
