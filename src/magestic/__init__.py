"""
MAGESTIC — Multiplexed Accurate Genome Editing with Short, Trackable, Integrated Cellular barcodes

Tools and pipelines for designing and analyzing MAGESTIC 3.0 libraries.

Modules:
    magestic.genomics    — Core genomics utilities (sequence, genome I/O, GFF annotation)
    magestic.utils       — Cross-pipeline shared utilities
    magestic.design      — Guide-donor oligo design (VCF-based and saturation)
    magestic.pipelines   — Analysis pipelines (guide_donor_bc0, bc1_donor_bc0, bc1_from_screen)
    magestic.plotting    — Locus visualization (gene tracks, VEP panels)
    magestic.projects    — Project-specific code (NRD1, VEP)

Citation:
    Roy KR, et al. (2018). Multiplexed precision genome editing with trackable genomic barcodes
    in yeast. Nature Biotechnology 36:512-520.

Author: Kevin R. Roy
"""

__version__ = "3.0.0"
__author__ = "Kevin R. Roy"
__email__ = "kevinrjroy@gmail.com"
