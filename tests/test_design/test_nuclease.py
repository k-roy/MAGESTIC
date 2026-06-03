"""
Unit tests for magestic.design.nuclease.

Tests cover PAM expansion, NucleaseConfig properties, and the nuclease
registry lookup.
"""

import pytest
from magestic.design.nuclease import (
    expand_pam_motif,
    NucleaseConfig,
    SPCAS9,
    SPG,
    LBCAS12A,
    IMP_LBCAS12A,
    NUCLEASE_REGISTRY,
    get_nuclease,
)


# ---------------------------------------------------------------------------
# expand_pam_motif
# ---------------------------------------------------------------------------

class TestExpandPamMotif:
    def test_ngg(self):
        pams = expand_pam_motif("NGG")
        assert sorted(pams) == ["AGG", "CGG", "GGG", "TGG"]

    def test_ngnn_count(self):
        # N×G×N×N → 4×1×4×4 = 64 combinations
        pams = expand_pam_motif("NGNN")
        assert len(pams) == 64

    def test_tttv(self):
        # V = A, C, G → 3 combinations
        pams = expand_pam_motif("TTTV")
        assert sorted(pams) == ["TTTA", "TTTC", "TTTG"]

    def test_literal_bases_unchanged(self):
        pams = expand_pam_motif("AAGG")
        assert pams == ["AAGG"]

    def test_sorted_output(self):
        pams = expand_pam_motif("NGG")
        assert pams == sorted(pams)


# ---------------------------------------------------------------------------
# NucleaseConfig properties
# ---------------------------------------------------------------------------

class TestNucleaseConfig:
    def test_pam_length_cas9(self):
        assert SPCAS9.pam_length == 3

    def test_pam_length_spg(self):
        assert SPG.pam_length == 4

    def test_all_pams_spg_count(self):
        # NGNN = 64 PAMs
        assert len(SPG.all_pams) == 64

    def test_pam_upstream_cas9(self):
        assert SPCAS9.pam_upstream is True

    def test_pam_upstream_cas12a(self):
        assert LBCAS12A.pam_upstream is False

    def test_re_site_rc(self):
        # re_site = "GTTTGAAGAGC" → RC = "GCTCTTCAAAC"
        assert SPCAS9.re_site_rc == "GCTCTTCAAAC"

    def test_restriction_sites_tuple(self):
        core, core_rc = SPCAS9.restriction_sites
        assert "GAAGAGC" in core
        assert "GCTCTTC" in core_rc

    def test_donor_length_for_oligo(self):
        # oligo=200, fwd=20, guide=20, internal=11, rev=20 → donor = 200-20-20-11-20 = 129
        donor_len = SPCAS9.donor_length_for_oligo(oligo_length=200, rev_primer_length=20)
        assert donor_len == 200 - len(SPCAS9.fwd_priming_site) - 20 - len(SPCAS9.internal_cloning_site) - 20

    def test_frozen(self):
        with pytest.raises(Exception):  # frozen=True → AttributeError
            SPCAS9.name = "mutated"

    def test_oligo_fixed_overhead(self):
        overhead = SPCAS9.oligo_fixed_overhead
        assert overhead == len(SPCAS9.fwd_priming_site) + 20 + len(SPCAS9.internal_cloning_site)


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

class TestNucleaseRegistry:
    def test_all_nucleases_present(self):
        for name in ("SpCas9", "SpG", "LbCas12a", "impLbCas12a"):
            assert name in NUCLEASE_REGISTRY

    def test_get_nuclease_spg(self):
        n = get_nuclease("SpG")
        assert n is SPG

    def test_get_nuclease_unknown(self):
        with pytest.raises(ValueError, match="Unknown nuclease"):
            get_nuclease("UnknownNuclease")
