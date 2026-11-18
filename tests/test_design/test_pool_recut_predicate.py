"""The re-cut predicate, tested against the incident that produced it.

NN 581. A finished 24,000-oligo pool shipped 127 designs whose donor still carried a perfect
``guide + NGG`` target, so the repaired site would be cut again and the edit lost. Two independent
causes, both measured from the molecules rather than inferred:

  * 111 designs (25.40 % of 437) were products of the B1 codon unit rule, which REBUILDS a donor
    when it splits a linked variant. A unit installs only its own substituted positions, so a unit
    whose positions lie outside the protospacer no longer disrupts the guide its parent disrupted.
    The design-time guard ran before the rebuild. Failure rate among everything else: 0.06 %.
  * 16 ``syn_recoded`` helpers were never disruption-checked at all. For 15 of the 16 the recoding
    lands 1 to 6 bases OUTSIDE the retained protospacer+PAM; the 16th changes only the ``N`` of
    ``NGG``, which is not a disruption.

Both classes are regression-tested here against :func:`donor_retains_target`, which is now the one
implementation shared by the gate and by every design-time step that rebuilds a donor.
"""

from __future__ import annotations

import pytest

from magestic.design.pool import donor_retains_target, restriction_site_offsets

GUIDE = "TTCCTTCTTCTGCTGCGCAG"
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


# --------------------------------------------------------------------------- the basic contract
def test_perfect_target_is_retained():
    assert donor_retains_target(f"AAAA{GUIDE}TGGAAAA", GUIDE) is True


def test_a_mismatch_inside_the_protospacer_disrupts():
    broken = GUIDE[:-1] + "A"
    assert donor_retains_target(f"AAAA{broken}TGGAAAA", GUIDE) is False


def test_every_ngg_pam_is_accepted():
    for pam in ("AGG", "CGG", "GGG", "TGG"):
        assert donor_retains_target(f"AAAA{GUIDE}{pam}AAAA", GUIDE) is True, pam


def test_breaking_either_g_of_the_pam_disrupts():
    assert donor_retains_target(f"AAAA{GUIDE}TAGAAAA", GUIDE) is False
    assert donor_retains_target(f"AAAA{GUIDE}TGAAAAA", GUIDE) is False


def test_the_test_is_strand_symmetric_and_the_fixture_is_not_palindromic():
    # 🔴 a palindromic fixture makes a strand test vacuous, so assert it is not one FIRST
    donor = "AAAACCA" + rc(GUIDE) + "AAAA"
    assert donor != rc(donor)
    assert GUIDE not in donor, "the guide must be absent from the plus strand for this to test -"
    assert donor_retains_target(donor, GUIDE) is True


# ------------------------------------------------------- cause 1: the unit rule's rebuilt donors
def test_a_unit_that_installs_positions_outside_the_protospacer_still_recuts():
    """The parent MNV disrupted the guide; this unit's own positions are downstream of the PAM."""
    parent_disrupting = f"AAAA{GUIDE[:-2]}TT" + "TGG" + "GGGGAAAA"
    assert donor_retains_target(parent_disrupting, GUIDE) is False
    # the unit installs only the trailing positions, leaving protospacer and PAM intact
    unit_only_downstream = f"AAAA{GUIDE}" + "TGG" + "TTTTAAAA"
    assert donor_retains_target(unit_only_downstream, GUIDE) is True


# ------------------------------------------- cause 2: the syn_recoded helpers that never disrupt
def test_an_edit_one_base_past_the_pam_does_not_disrupt():
    """The measured YAL054C shape: the site ends at 45922 and the edit sits at 45923."""
    assert donor_retains_target(f"AAAA{GUIDE}TGG" + "T" + "AAAA", GUIDE) is True


@pytest.mark.parametrize("gap", [1, 2, 3, 4, 5, 6])
def test_an_edit_adjacent_to_the_site_never_disrupts_it(gap):
    """15 of the 16 helpers missed their own site by 1 to 6 bases. Every distance still re-cuts."""
    assert donor_retains_target("A" * gap + GUIDE + "TGG" + "A" * 10, GUIDE) is True
    assert donor_retains_target("A" * 10 + GUIDE + "TGG" + "A" * gap, GUIDE) is True


def test_changing_only_the_n_of_ngg_is_not_a_disruption():
    """The 16th helper. TGG -> AGG keeps a perfect NGG, so the site survives."""
    assert donor_retains_target(f"AAAA{GUIDE}AGGAAAA", GUIDE) is True


# ------------------------------------------------------------------------------ degenerate input
@pytest.mark.parametrize("donor,guide", [("", GUIDE), (f"AAAA{GUIDE}TGG", ""), ("", "")])
def test_missing_sequence_cannot_claim_a_retained_target(donor, guide):
    assert donor_retains_target(donor, guide) is False


def test_a_donor_shorter_than_the_site_cannot_retain_it():
    assert donor_retains_target(GUIDE[:10], GUIDE) is False


def test_lowercase_input_is_handled():
    assert donor_retains_target(f"aaaa{GUIDE.lower()}tggaaaa", GUIDE.lower()) is True


def test_an_empty_pam_cannot_match_everything():
    assert donor_retains_target(f"AAAA{GUIDE}TGGAAAA", GUIDE, pam="") is False


# ---------------------------------------------------------------- the restriction-site predicate
def test_offsets_find_both_orientations_of_the_sapi_site():
    assert restriction_site_offsets("AAGCTCTTCAAGAAGAGCTT", ("GCTCTTC", "GAAGAGC")) == [2, 11]


def test_offsets_are_sorted_and_count_overlapping_occurrences():
    seq = "GCTCTTCGCTCTTC"
    assert restriction_site_offsets(seq, ("GCTCTTC",)) == [0, 7]


def test_offsets_ignore_an_empty_site_rather_than_matching_everywhere():
    assert restriction_site_offsets("ACGTACGT", ("", "ACGT")) == [0, 4]


def test_offsets_on_empty_sequence_is_empty():
    assert restriction_site_offsets("", ("GCTCTTC",)) == []


# ------------------------------------------------------- IUPAC PAMs and PAM-before-guide layout
def test_an_iupac_pam_is_matched_as_a_pattern():
    assert donor_retains_target(f"AAAA{GUIDE}TGGAAAA", GUIDE, pam="NRG") is True   # TGG: R=A/G
    assert donor_retains_target(f"AAAA{GUIDE}TTGAAAA", GUIDE, pam="NRG") is False  # TTG: T not R


def test_pam_first_places_the_pam_before_the_guide():
    cas12a = f"AAAA{'TTTA'}{GUIDE}AAAA"
    assert donor_retains_target(cas12a, GUIDE, pam="TTTV", pam_first=True) is True
    # the same donor read with the Cas9 layout must NOT match
    assert donor_retains_target(cas12a, GUIDE, pam="TTTV", pam_first=False) is False


def test_ttt_t_does_not_satisfy_a_tttv_pam():
    assert donor_retains_target(
        f"AAAA{'TTTT'}{GUIDE}AAAA", GUIDE, pam="TTTV", pam_first=True
    ) is False


# -------------------------------------------- the off-by-one that both older copies carried
def test_a_site_flush_against_the_three_prime_end_is_retained():
    """The window count is len - span + 1. Dropping the +1 never examines the last window."""
    flush = "AAAA" + GUIDE + "TGG"
    assert flush.endswith("TGG")
    assert donor_retains_target(flush, GUIDE) is True


def test_a_site_filling_the_whole_donor_is_retained():
    assert donor_retains_target(GUIDE + "TGG", GUIDE) is True


def test_the_legacy_guide_disruption_now_agrees_with_the_predicate_at_the_boundary():
    """Both older copies reported a flush site as DISRUPTED, a false safe verdict."""
    from magestic.design.guide_donor_functions import guide_disruption as legacy
    from magestic.design.oligo_generation import guide_disruption as legacy2

    flush = "AAAA" + GUIDE + "TGG"
    # legacy returns True for "disrupted", so a retained site must give False
    assert legacy(GUIDE, "TGG", flush, True) is False
    assert legacy2(GUIDE, "NGG", flush, True, len(flush)) is False
    assert donor_retains_target(flush, GUIDE) is True


def test_legacy_and_predicate_agree_across_every_offset():
    """Sweep the site through the donor; the two implementations must never disagree."""
    from magestic.design.guide_donor_functions import guide_disruption as legacy

    for pad_left in range(0, 6):
        for pad_right in range(0, 6):
            donor = "A" * pad_left + GUIDE + "TGG" + "A" * pad_right
            assert legacy(GUIDE, "TGG", donor, True) is not donor_retains_target(donor, GUIDE), (
                f"disagreement at pad {pad_left}/{pad_right}"
            )
