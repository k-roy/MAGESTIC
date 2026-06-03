"""
magestic.plotting — Modular visualization for MAGESTIC locus plots.

Ported from common/scripts/locus_plotting/.

Modules:
    config  — Color schemes, markers, constants
    panels  — Individual panel drawing functions
      gene_track  — Box-arrow gene schematic
      vep_panels  — VEP prediction panels (Evo2, Yorzoi, Shorkie, ESM1v)

Usage:
    from magestic.plotting import (
        CODON_VARIANT_COLORS,
        GENE_TYPE_COLORS,
        STRAIN_MARKERS,
    )

    from magestic.plotting.panels import (
        draw_gene_track,
        draw_vep_scatter,
    )

Author: Kevin R. Roy
"""

from .config import (
    # Color schemes
    CODON_VARIANT_COLORS,
    GENE_TYPE_COLORS,
    STRAIN_ORIGIN_COLORS,
    CONTROL_TYPE_COLORS,
    EFFECT_DIRECTION_COLORS,
    # Marker shapes
    SUBLIBRARY_MARKERS,
    STRAIN_ORIGIN_MARKERS,
    STRAIN_MARKERS,
    GENE_TYPE_MARKERS,
    SOURCE_MARKERS,
    # Column backgrounds
    COLUMN_COLORS,
    # Model info
    MODEL_NAMES,
    MODEL_SCORE_COLUMNS,
    # Pilot genes
    PILOT_GENES,
)

__all__ = [
    # Color schemes
    "CODON_VARIANT_COLORS",
    "GENE_TYPE_COLORS",
    "STRAIN_ORIGIN_COLORS",
    "CONTROL_TYPE_COLORS",
    "EFFECT_DIRECTION_COLORS",
    # Marker shapes
    "SUBLIBRARY_MARKERS",
    "STRAIN_ORIGIN_MARKERS",
    "STRAIN_MARKERS",
    "GENE_TYPE_MARKERS",
    "SOURCE_MARKERS",
    # Column backgrounds
    "COLUMN_COLORS",
    # Model info
    "MODEL_NAMES",
    "MODEL_SCORE_COLUMNS",
    # Pilot genes
    "PILOT_GENES",
]
