"""Ambiguity flags ship with the match (planning 30_ arm_bdb-F1, IO-XF1).

An ambiguous match names only the sorted-first candidate, so the flags and the candidate list are the only record
that other edits were equally consistent with the read. Both matcher paths are covered: the vectorized
perfect-donor fast path (anchor-bearing donors, no guide) and the per-row fallback (sliding window)."""
import random

import pandas as pd

from magestic.pipelines.bc1_donor_bc0.core.matching import (
    DONOR_VARIABLE_REGION_ANCHORS, build_lookup_tables, format_candidates, match_dataframe)
from magestic.utils.oligo_annotations import COLUMNS_BC1_DONOR_BC0

rng = random.Random(20260910)


def _seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


D_SHARED, D_DUP, D_WIN = _seq(129), _seq(129), _seq(129)
# A second design whose donor shares D_WIN's first 70 nt (a 60-mer window in common), installing another edit.
D_WIN2 = D_WIN[:70] + _seq(59)
OLIGOS = pd.DataFrame({
    "oligo_name": ["SpCas9_a", "SpCas9_b", "LbCas12a_c", "SpCas9_d", "SpCas9_e", "SpCas9_f"],
    "guide": [_seq(20) for _ in range(6)],
    "donor": [D_SHARED, D_SHARED, D_DUP, D_DUP, D_WIN, D_WIN2],
    # a/b: two guides, ONE edit.  c/d: identical donor, two different edits.  e/f: share a window, two edits.
    "individual_variants": ["100_A_C_SNP", "100_A_C_SNP", "200_G_T_SNP", "300_C_A_SNP", "400_T_G_SNP", "500_A_G_SNP"],
})
ANCHOR = DONOR_VARIABLE_REGION_ANCHORS[0]


def _match(df_oligos, donors):
    tables = build_lookup_tables(df_oligos)
    return match_dataframe(pd.DataFrame({"donor": donors}), tables).set_index("donor")


def test_loader_requests_variant_identity():
    assert "individual_variants" in COLUMNS_BC1_DONOR_BC0


def test_fast_path_flags_and_candidates():
    q_shared, q_dup = "ACGTACGT" + ANCHOR + D_SHARED, "ACGTACGT" + ANCHOR + D_DUP
    out = _match(OLIGOS, [q_shared, q_dup])
    r = out.loc[q_shared]
    assert (r.oligo_name, r.match_method, bool(r.is_unambiguous)) == ("SpCas9_a", "perfect", True)
    assert (r.num_candidates, r.num_candidate_variants, r.variant) == (2, 1, "100_A_C_SNP")
    assert r.candidate_oligos is None                     # same edit twice is not ambiguous
    r = out.loc[q_dup]
    assert (r.oligo_name, r.match_method, bool(r.is_unambiguous)) == ("LbCas12a_c", "perfect_ambiguous", False)
    assert (r.num_candidates, r.num_candidate_variants) == (2, 2)
    assert r.candidate_oligos == "LbCas12a_c;SpCas9_d"   # the other edit is visible


def test_one_edit_two_encodings_is_unambiguous_only_with_variant_column():
    # IO-XF1: one edit designed with two different donor encodings that share a 60-mer window. With individual_variants
    # the truncated read is unambiguous; without it the matcher falls back to donor SEQUENCE and calls 2 variants.
    enc1 = _seq(129)
    enc2 = enc1[:70] + _seq(59)
    olig = pd.DataFrame({"oligo_name": ["SpCas9_g", "SpCas9_h"], "guide": [_seq(20), _seq(20)],
                         "donor": [enc1, enc2], "individual_variants": ["600_C_T_SNP", "600_C_T_SNP"]})
    q = enc1[:65]
    r = _match(olig, [q]).loc[q]
    assert (r.match_method, bool(r.is_unambiguous), r.variant) == ("sliding_window", True, "600_C_T_SNP")
    assert r.candidate_oligos is None
    r = _match(olig.drop(columns=["individual_variants"]), [q]).loc[q]
    assert (r.match_method, bool(r.is_unambiguous), r.num_candidate_variants) == ("sliding_window_ambiguous", False, 2)


def test_fallback_sliding_window_flags_and_candidates():
    q = D_WIN[:65]                                        # truncated donor: only the shared 60-mer windows match
    r = _match(OLIGOS, [q]).loc[q]
    assert r.match_method == "sliding_window_ambiguous" and not bool(r.is_unambiguous)
    assert r.candidate_oligos == "SpCas9_e;SpCas9_f" and r.oligo_name == "SpCas9_e"


def test_candidate_cap():
    names = [f"o{i:03d}" for i in range(60)]
    s = format_candidates(names, 2)
    assert s.endswith(";+10_more") and s.count(";") == 50
    assert format_candidates(names, 1) is None
