"""
magestic.design.nuclease — Nuclease system configuration for MAGESTIC oligo design.

Each CRISPR nuclease used in MAGESTIC has a distinct PAM requirement, guide
orientation (upstream vs. downstream of PAM), restriction enzyme cloning site,
and poly-T tolerance (for U6 promoter-driven expression).

Typical usage::

    from magestic.design.nuclease import SPG, SPCAS9, LBCAS12A, expand_pam_motif

    pams = expand_pam_motif(SPG.pam_pattern)   # ['AGAG', 'AGAT', ..., 'TGTG']
    print(SPG.pam_length)                        # 4

Author: Kevin R. Roy <kevinrjroy@gmail.com>
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property
from typing import List


# ---------------------------------------------------------------------------
# IUPAC PAM expansion
# ---------------------------------------------------------------------------

_IUPAC_EXPAND: dict[str, str] = {
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}


def expand_pam_motif(motif: str) -> List[str]:
    """
    Expand an IUPAC PAM motif into all concrete 4-base sequences it matches.

    Parameters
    ----------
    motif : str
        IUPAC PAM motif (e.g. ``"NGG"``, ``"NGNN"``, ``"TTTV"``).

    Returns
    -------
    List[str]
        Sorted list of all unambiguous sequences compatible with *motif*.

    Examples
    --------
    >>> expand_pam_motif("NGG")
    ['AGG', 'CGG', 'GGG', 'TGG']
    >>> len(expand_pam_motif("NGNN"))
    64
    """
    # Start with one partial sequence per character slot
    expanded: list[list[str]] = [[c] for c in motif]
    # Expand each position that carries an IUPAC code
    result: list[str] = [""]
    for pos_chars in expanded:
        new_result: list[str] = []
        base = pos_chars[0]
        alternatives = _IUPAC_EXPAND.get(base, base)  # if not IUPAC, literal char
        for prefix in result:
            for alt in alternatives:
                new_result.append(prefix + alt)
        result = new_result
    return sorted(result)


# ---------------------------------------------------------------------------
# NucleaseConfig dataclass
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class NucleaseConfig:
    """
    Immutable configuration for a CRISPR nuclease used in MAGESTIC library design.

    Parameters
    ----------
    name : str
        Human-readable nuclease name (e.g. ``"SpCas9"``, ``"SpG"``, ``"LbCas12a"``).
    pam_pattern : str
        IUPAC PAM motif (e.g. ``"NGG"``, ``"NGNN"``, ``"TTTV"``).
    guide_length : int
        Length of the spacer sequence (bp). Default 20.
    pam_upstream : bool
        ``True`` if the PAM is 3′ of the guide (Cas9 family);
        ``False`` if the PAM is 5′ of the guide (Cas12a family).
    re_site : str
        Internal restriction enzyme recognition sequence used for Golden Gate
        cloning of the guide into the oligo backbone (e.g. ``"GTTTGAAGAGC"``
        for BspQI/SapI).
    max_poly_t : int
        Maximum run of consecutive T's tolerated in the guide sequence.
        Longer runs terminate U6-driven transcription. Default 4.
    fwd_priming_site : str
        Forward PCR priming sequence prepended to every oligo
        (lowercase in oligo; used for subpool amplification).
    internal_cloning_site : str
        Linker between guide and donor in the final oligo (lowercase).
        Contains the internal RE recognition site.
    """

    name: str
    pam_pattern: str
    guide_length: int = 20
    pam_upstream: bool = True
    re_site: str = "GTTTGAAGAGC"
    max_poly_t: int = 4
    fwd_priming_site: str = "ctgcgattggcaggcgcgcc"
    internal_cloning_site: str = "gtttgaagagc"

    # ------------------------------------------------------------------
    # Derived properties (computed once on first access)
    # ------------------------------------------------------------------

    @cached_property
    def pam_length(self) -> int:
        """Number of bases in the PAM motif."""
        return len(self.pam_pattern)

    @cached_property
    def all_pams(self) -> List[str]:
        """All concrete PAM sequences compatible with :attr:`pam_pattern`."""
        return expand_pam_motif(self.pam_pattern)

    @cached_property
    def re_site_rc(self) -> str:
        """Reverse complement of :attr:`re_site` (the RC cuts on the other strand)."""
        from magestic.design.guide_donor_functions import rev_comp
        return rev_comp(self.re_site)

    @cached_property
    def restriction_sites(self) -> tuple[str, str]:
        """Both orientations of the internal RE site that must be excluded from donors."""
        # Strip to the 7-base recognition core (GAAGAGC / GCTCTTC for BspQI)
        core = self.re_site[-7:] if len(self.re_site) >= 7 else self.re_site
        return core, self.re_site_rc[:7] if len(self.re_site_rc) >= 7 else self.re_site_rc

    @property
    def oligo_fixed_overhead(self) -> int:
        """Bases consumed by priming + guide + cloning site (not including donor or rev primer)."""
        return len(self.fwd_priming_site) + self.guide_length + len(self.internal_cloning_site)

    def donor_length_for_oligo(self, oligo_length: int, rev_primer_length: int) -> int:
        """
        Calculate the donor region length given a target total oligo length.

        Parameters
        ----------
        oligo_length : int
            Desired total oligo length in bases (e.g. 200).
        rev_primer_length : int
            Length of the reverse subpool priming sequence.

        Returns
        -------
        int
            Number of bases available for the donor (homology arms + variant).
        """
        return oligo_length - self.oligo_fixed_overhead - rev_primer_length

    def __str__(self) -> str:
        orientation = "upstream" if self.pam_upstream else "downstream"
        return (
            f"NucleaseConfig({self.name}, PAM={self.pam_pattern}, "
            f"guide={self.guide_length} nt, PAM {orientation} of guide)"
        )


# ---------------------------------------------------------------------------
# Pre-defined nuclease configs shipped with the package
# ---------------------------------------------------------------------------

#: SpCas9 with canonical NGG PAM (Jinek et al., 2012).
SPCAS9 = NucleaseConfig(
    name="SpCas9",
    pam_pattern="NGG",
    guide_length=20,
    pam_upstream=True,
    re_site="GTTTGAAGAGC",
    max_poly_t=4,
    fwd_priming_site="ctgcgattggcaggcgcgcc",
    internal_cloning_site="gtttgaagagc",
)

#: SpG, the relaxed-PAM SpCas9 variant targeting NGNN PAMs (Walton et al., 2020).
#: Enables ~4× more target sites than canonical SpCas9.
SPG = NucleaseConfig(
    name="SpG",
    pam_pattern="NGNN",
    guide_length=20,
    pam_upstream=True,
    re_site="GTTTGAAGAGC",
    max_poly_t=4,
    fwd_priming_site="ctgcgattggcaggcgcgcc",
    internal_cloning_site="gtttgaagagc",
)

#: LbCas12a with canonical TTTV PAM (Zetsche et al., 2015).
LBCAS12A = NucleaseConfig(
    name="LbCas12a",
    pam_pattern="TTTV",
    guide_length=20,
    pam_upstream=False,  # PAM is 5′ of the protospacer for Cas12a
    re_site="TTTCGAAGAGC",
    max_poly_t=5,
    fwd_priming_site="ctgcgattggcaggcgcgcc",
    internal_cloning_site="tttcgaagagc",
)

#: Improved LbCas12a (impLbCas12a) with expanded PAM set including TTTV,
#: TNTN, TACV, TTCV, TCCV, CTCV, CCCV, VTTV PAMs.
IMP_LBCAS12A = NucleaseConfig(
    name="impLbCas12a",
    pam_pattern="TNTN",  # Represents one of several accepted motifs; see all_pams
    guide_length=20,
    pam_upstream=False,
    re_site="TTTCGAAGAGC",
    max_poly_t=5,
    fwd_priming_site="ctgcgattggcaggcgcgcc",
    internal_cloning_site="tttcgaagagc",
)

#: Registry mapping CLI name strings to NucleaseConfig objects.
NUCLEASE_REGISTRY: dict[str, NucleaseConfig] = {
    "SpCas9": SPCAS9,
    "SpG": SPG,
    "LbCas12a": LBCAS12A,
    "impLbCas12a": IMP_LBCAS12A,
}


def get_nuclease(name: str) -> NucleaseConfig:
    """
    Look up a pre-defined nuclease by name.

    Parameters
    ----------
    name : str
        Nuclease name (case-sensitive). One of ``"SpCas9"``, ``"SpG"``,
        ``"LbCas12a"``, ``"impLbCas12a"``.

    Returns
    -------
    NucleaseConfig

    Raises
    ------
    ValueError
        If *name* is not found in the registry.
    """
    if name not in NUCLEASE_REGISTRY:
        raise ValueError(
            f"Unknown nuclease '{name}'. "
            f"Available: {sorted(NUCLEASE_REGISTRY)}"
        )
    return NUCLEASE_REGISTRY[name]
