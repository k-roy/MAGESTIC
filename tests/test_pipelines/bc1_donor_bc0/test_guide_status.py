"""guide_status separates the two things guide_match=False used to mean (planning 30_ step 6).

On the 2026-09-10 rebuild, guide_match=False covered 711,464 rows with no linked guide at all AND
223,954 rows whose linked guide disagreed with the donor. Worse, the matcher narrows by guide only
`if overlap:`, so a designed guide that named no design the donor allowed was silently discarded --
which is how ambiguous rows came to carry a guide matching none of their own candidates.
"""
import pandas as pd
import pytest

from magestic.pipelines.bc1_donor_bc0.core.matching import (
    GUIDE_CONFLICT, GUIDE_MATCH, GUIDE_NOT_DESIGNED, GUIDE_NO_GUIDE,
    build_lookup_tables, match_dataframe, match_sequence)

DONOR_A = "A" * 60 + "TTTT" + "A" * 65      # shared by two designs that install DIFFERENT variants
DONOR_B = "C" * 129
GUIDE_A1 = "AAAACCCCGGGGTTTTAAAA"
GUIDE_A2 = "CCCCGGGGTTTTAAAACCCC"
GUIDE_B = "GGGGTTTTAAAACCCCGGGG"
GUIDE_UNDESIGNED = "TTTTTTTTTTTTTTTTTTTT"

OLIGOS = pd.DataFrame({
    "oligo_name": ["d_a1", "d_a2", "d_b"],
    "guide": [GUIDE_A1, GUIDE_A2, GUIDE_B],
    "donor": [DONOR_A, DONOR_A, DONOR_B],
    "individual_variants": ["10_A_C_SNP", "20_G_T_SNP", "30_C_G_SNP"],
})
TABLES = build_lookup_tables(OLIGOS)


def test_no_guide():
    r = match_sequence(DONOR_A, None, TABLES)
    assert r.guide_status == GUIDE_NO_GUIDE and r.guide_match is False and r.num_guide_oligos == 0


def test_guide_not_designed():
    r = match_sequence(DONOR_A, GUIDE_UNDESIGNED, TABLES)
    assert r.guide_status == GUIDE_NOT_DESIGNED and r.guide_match is False and r.num_guide_oligos == 0


def test_guide_match_narrows_an_ambiguous_donor():
    ambiguous = match_sequence(DONOR_A, None, TABLES)
    assert ambiguous.num_candidate_variants == 2        # the donor alone cannot decide
    r = match_sequence(DONOR_A, GUIDE_A1, TABLES)
    assert r.guide_status == GUIDE_MATCH and r.guide_match is True
    assert r.oligo_name == "d_a1" and r.is_unambiguous and r.num_guide_oligos == 1


def test_conflict_is_recorded_instead_of_discarded():
    """GUIDE_B is a real designed guide, but it names only d_b, which DONOR_A does not allow. The donor
    answer is still kept -- what changes is that the disagreement is now on the row."""
    r = match_sequence(DONOR_A, GUIDE_B, TABLES)
    assert r.guide_status == GUIDE_CONFLICT
    assert r.guide_match is False                        # unchanged, back-compat
    assert r.num_guide_oligos == 1
    assert r.oligo_name in ("d_a1", "d_a2")              # donor answer kept
    assert r.num_candidate_variants == 2                 # and it is still reported as ambiguous


def test_guide_only_match_when_there_is_no_donor_answer():
    r = match_sequence("G" * 129, GUIDE_B, TABLES)       # donor matches nothing
    assert r.guide_status == GUIDE_MATCH and r.match_method == "guide_only" and r.oligo_name == "d_b"


def test_dataframe_reports_every_state_including_through_the_vectorized_fast_path():
    df = pd.DataFrame({
        "donor": [DONOR_A, DONOR_A, DONOR_A, DONOR_B, DONOR_A],
        "guide": [None, GUIDE_UNDESIGNED, GUIDE_A2, GUIDE_B, GUIDE_B],
    })
    out = match_dataframe(df, TABLES, donor_col="donor", guide_col="guide")
    assert out["guide_status"].tolist() == [
        GUIDE_NO_GUIDE, GUIDE_NOT_DESIGNED, GUIDE_MATCH, GUIDE_MATCH, GUIDE_CONFLICT]
    assert out["guide_match"].tolist() == [False, False, True, True, False]
    assert out["num_guide_oligos"].tolist() == [0, 0, 1, 1, 1]
    # guide_match=False no longer conflates: rows 0 and 1 have no usable guide, row 4 has a disagreement
    assert (out["guide_status"] == GUIDE_CONFLICT).sum() == 1


def test_guide_match_is_exactly_status_match():
    """The back-compat invariant: nothing that reads guide_match changes behaviour."""
    for donor, guide in [(DONOR_A, None), (DONOR_A, GUIDE_UNDESIGNED), (DONOR_A, GUIDE_A1),
                         (DONOR_A, GUIDE_B), (DONOR_B, GUIDE_B), ("G" * 129, GUIDE_B),
                         ("G" * 129, GUIDE_UNDESIGNED)]:
        r = match_sequence(donor, guide, TABLES)
        assert r.guide_match is (r.guide_status == GUIDE_MATCH), (donor[:5], guide)
