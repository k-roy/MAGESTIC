"""
Panel drawing functions for locus plots.

Each module provides draw_*() functions for specific panel types.
Panels are designed to be composable for different figure layouts.

Available Panels:

VEP Panels (vep_panels.py):
  - draw_evo2_panel: Evo2 scores (colored by codon_variant_type)
  - draw_esm1v_panel: ESM1v scores (colored by codon_variant_type)
  - draw_shorkie_pergene_panel: Shorkie per-gene effects (colored by gene_type)
  - draw_yorzoi_pergene_panel: Yorzoi per-gene effects (colored by gene_type)
  - draw_vep_panel_styled: Generic styled VEP panel

Gene Track (gene_track.py):
  - draw_gene_track: Box-arrow gene structure diagrams
  - get_genes_in_region: Get genes overlapping a region

Note: Only complete variants are shown in VEP panels (sub-variants filtered out).

Author: Kevin R. Roy
Date: 2026-02-10
"""

# VEP panel functions
from .vep_panels import (
    draw_evo2_panel,
    draw_esm1v_panel,
    draw_shorkie_pergene_panel,
    draw_yorzoi_pergene_panel,
    draw_vep_panel_styled,
    extract_position_from_variants,
)

# Gene track functions
from .gene_track import (
    draw_gene_track,
    assign_feature_levels,
    get_genes_in_region,
)

__all__ = [
    # VEP panels (genomic/protein models - colored by codon_variant_type)
    "draw_evo2_panel",
    "draw_esm1v_panel",
    # VEP panels (RNA expression models - colored by gene_type)
    "draw_shorkie_pergene_panel",
    "draw_yorzoi_pergene_panel",
    # Generic VEP panel
    "draw_vep_panel_styled",
    "extract_position_from_variants",
    # Gene track
    "draw_gene_track",
    "assign_feature_levels",
    "get_genes_in_region",
]
