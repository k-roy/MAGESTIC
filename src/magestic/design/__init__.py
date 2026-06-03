"""
magestic.design — Guide-donor oligo design for MAGESTIC 3.0 libraries.

Provides functions for designing guide RNA / donor DNA oligonucleotide pairs
for MAGESTIC libraries. Supports SpCas9 (NGG PAM), SpG (NGNN PAM), LbCas12a
(TTTV PAM), and impLbCas12a (expanded PAM set) nucleases.

Two design workflows:
    1. VCF-based design:  generate oligos for natural variants from a VCF file
    2. Saturation design: generate all amino acid substitutions for a target ORF

Modules:
    nuclease              — NucleaseConfig dataclass and pre-defined nuclease configs
    guide_donor_functions — Reference data loaders, sequence utils, codon utils,
                             guide disruption checks, donor assembly
    oligo_design          — ORF-level PAM finding, Azimuth score loading,
                             restriction site avoidance, donor shifting
    vcf_processing        — VCF loading, linked-variant clustering, gene-list filtering
    pam_search            — Genome-wide PAM site enumeration and off-target analysis
    guide_selection       — Allowable guide classification from off-target analysis
    oligo_generation      — VCF-based guide-donor oligo design orchestrator

CLI entry points (see magestic.design._cli):
    magestic-design-vcf        — VCF-based oligo design
    magestic-design-saturation — Saturation mutagenesis oligo design
    magestic-pam-search        — Genome-wide PAM site enumeration
"""

from .guide_donor_functions import (
    load_codon_table,
    load_processed_guides,
    load_guides_from_fasta,
    rev_comp,
    base_match,
    hamming_distance,
    restriction_site_in_seq,
    guide_disruption,
    assemble_donor,
)
from .nuclease import (
    NucleaseConfig,
    expand_pam_motif,
    get_nuclease,
    SPCAS9,
    SPG,
    LBCAS12A,
    IMP_LBCAS12A,
    NUCLEASE_REGISTRY,
)

__all__ = [
    # guide_donor_functions
    "load_codon_table",
    "load_processed_guides",
    "load_guides_from_fasta",
    "rev_comp",
    "base_match",
    "hamming_distance",
    "restriction_site_in_seq",
    "guide_disruption",
    "assemble_donor",
    # nuclease
    "NucleaseConfig",
    "expand_pam_motif",
    "get_nuclease",
    "SPCAS9",
    "SPG",
    "LBCAS12A",
    "IMP_LBCAS12A",
    "NUCLEASE_REGISTRY",
]
