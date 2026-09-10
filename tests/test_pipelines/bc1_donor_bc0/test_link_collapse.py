"""Tests for core/link_collapse.py and the label-not-filter semantics (planning 28_, 2026-09-10).

Synthetic only -- no Oak paths. Each test pins one rule that the pre-2026-09-10 code broke.
"""
import numpy as np
import pandas as pd
import pytest

from magestic.pipelines.bc1_donor_bc0.core.link_collapse import (
    aggregate_links, attach_short_donors, collapse_links, decode_masks, mask_columns,
    masks_for_codes, popcount_rows)

N_IDS = 131                      # spans three uint64 words on purpose
MC = mask_columns(N_IDS)


def _with_masks(df, codes):
    m = masks_for_codes(np.asarray(codes), N_IDS)
    for k, c in enumerate(MC):
        df[c] = m[:, k]
    return df


def _support(df):
    return popcount_rows(df[MC].to_numpy(dtype=np.uint64))


def test_bitmask_roundtrip_across_words():
    m = masks_for_codes(np.array([0, 63, 64, 130]), N_IDS)
    assert m.shape == (4, 3)
    assert list(popcount_rows(m)) == [1, 1, 1, 1]
    names = [f"r{i}" for i in range(N_IDS)]
    assert decode_masks(m, names) == ["r0", "r63", "r64", "r130"]
    assert decode_masks((m[0] | m[3])[None, :], names) == ["r0;r130"]


def test_aggregate_sums_ors_and_preserves_top_row():
    df = _with_masks(pd.DataFrame({"bc1": list("AABA"), "donor": list("ddde"), "bc0": list("xxyx"),
                                   "counts": [5, 9, 1, 2], "SPS": ["s1", "s2", "s3", "s4"]}),
                     [0, 1, 0, 0])
    df.loc[3, "donor"] = "d"
    a = aggregate_links(df, ["bc1", "donor", "bc0"], "counts", ["counts"], MC, ["SPS"])
    ra = a[a.bc1 == "A"]
    assert ra.counts.iloc[0] == 16 and ra.SPS.iloc[0] == "s2"
    assert _support(ra)[0] == 2


def test_aggregate_is_associative():
    df = _with_masks(pd.DataFrame({"bc1": list("AABA"), "donor": list("dded"), "bc0": list("xxyx"),
                                   "counts": [5, 9, 1, 2]}), [0, 1, 0, 3])
    k = ["bc1", "donor", "bc0"]
    whole = aggregate_links(df, k, "counts", ["counts"], MC)
    parts = aggregate_links(pd.concat([aggregate_links(df.iloc[:2], k, "counts", ["counts"], MC),
                                       aggregate_links(df.iloc[2:], k, "counts", ["counts"], MC)]),
                            k, "counts", ["counts"], MC)
    srt = lambda x: x.sort_values(k).reset_index(drop=True)
    pd.testing.assert_frame_equal(srt(whole), srt(parts))


def test_aggregate_refuses_nan_keys():
    df = pd.DataFrame({"bc1": ["A"], "donor": [None], "bc0": ["x"], "counts": [1]})
    with pytest.raises(ValueError):
        aggregate_links(df, ["bc1", "donor", "bc0"], "counts", ["counts"], [])


def test_collapse_unions_replicates_never_sums():
    base = pd.DataFrame({"bc1": ["A", "A"], "bc0": ["x", "x"],
                         "donor": ["ACGTACGTAC", "ACGTACGTAA"], "counts": [10, 1]})
    same = collapse_links(_with_masks(base.copy(), [3, 3]), ["bc1", "bc0"], "donor", "counts",
                          ["counts"], MC, threshold=5)
    assert len(same) == 1 and same.counts.iloc[0] == 11 and _support(same)[0] == 1
    diff = collapse_links(_with_masks(base.copy(), [3, 7]), ["bc1", "bc0"], "donor", "counts",
                          ["counts"], MC, threshold=5)
    assert _support(diff)[0] == 2


def test_collapse_ranks_by_counts_only():
    df = _with_masks(pd.DataFrame({"bc1": ["A", "A"], "bc0": ["x", "x"],
                                   "donor": ["AAAAAAAAAA", "AAAAAAAAAT"], "counts": [100, 1]}),
                     [0, 130])        # the minority row's mask is a huge integer in word 2
    out = collapse_links(df, ["bc1", "bc0"], "donor", "counts", ["counts"], MC, threshold=2)
    assert out.donor.iloc[0] == "AAAAAAAAAA"


def test_collapse_never_filters():
    df = _with_masks(pd.DataFrame({"bc1": ["A", "C"], "bc0": ["x", "y"],
                                   "donor": ["AAAAAAAAAA", "CCCCCCCCCC"], "counts": [1, 1]}), [0, 1])
    assert len(collapse_links(df, ["bc1", "bc0"], "donor", "counts", ["counts"], MC)) == 2


LONG = "G" * 150 + "ACGTACGTAC" * 8 + "TTTTTTTT"      # 238 nt, a 2x300-era donor
SHORT = LONG[-88:]                                   # the 2x150-era read of the same donor


def test_short_donor_attaches_to_its_unique_long_parent():
    df = _with_masks(pd.DataFrame({"bc1": ["A", "A", "A", "B"],
                                   "bc0": ["CCCCCCCCCC", "CCCCCCCCCC", "CCCCCCCCCA", "GGGGGGGGGG"],
                                   "donor": [LONG, SHORT, "T" * 88, SHORT], "counts": [50, 30, 4, 7]}),
                     [0, 1, 2, 3])
    out, st = attach_short_donors(df, "donor", "counts", ["counts"], MC)
    parent = out[(out.bc1 == "A") & (out.donor == LONG)]
    assert parent.counts.iloc[0] == 80 and _support(parent)[0] == 2
    assert parent.donor_lengths_observed.iloc[0] == "238;88"
    assert ((out.bc1 == "A") & (out.donor == "T" * 88)).sum() == 1      # not a suffix: stays
    assert (out.bc1 == "B").sum() == 1                                  # other bc1: untouched
    assert st["attached"] == 1 and st["no_parent"] == 1


def test_short_donor_with_two_parents_is_flagged_not_merged():
    long2 = "C" * 150 + LONG[-88:]
    df = _with_masks(pd.DataFrame({"bc1": ["A"] * 3, "bc0": ["CCCCCCCCCC"] * 3,
                                   "donor": [LONG, long2, SHORT], "counts": [50, 40, 5]}), [0, 1, 2])
    out, st = attach_short_donors(df, "donor", "counts", ["counts"], MC)
    assert len(out) == 3 and out.donor_suffix_ambiguous.sum() == 1 and st["ambiguous"] == 1


def test_suffix_chain_resolves_longest_first():
    df = _with_masks(pd.DataFrame({"bc1": ["A"] * 3, "bc0": ["CCCCCCCCCC"] * 3,
                                   "donor": [LONG, LONG[-150:], SHORT], "counts": [50, 20, 5]}),
                     [0, 1, 2])
    out, _ = attach_short_donors(df, "donor", "counts", ["counts"], MC)
    assert len(out) == 1 and out.counts.iloc[0] == 75 and _support(out)[0] == 3


@pytest.mark.parametrize("codes", [(0, 60), (0, 63), (1, 52, 53), (5, 64 + 60), (0, 127, 130)])
def test_collapse_preserves_every_mask_bit_alongside_singletons(codes):
    """Regression (review 2026-09-10): uint64 singletons + inferred-int64 merged rows
    concatenated to float64 and dropped low bits of any word that also held a bit >= 2**53."""
    n = len(codes)
    df = pd.DataFrame({"bc1": ["A"] * n + ["B"], "bc0": ["x"] * n + ["y"],
                       "donor": ["AAAAAAAAAA"] + ["AAAAAAAAAT"] * (n - 1) + ["CCCCCCCCCC"],
                       "counts": [10] + [1] * (n - 1) + [5]})
    df = _with_masks(df, list(codes) + [7])
    out = collapse_links(df, ["bc1", "bc0"], "donor", "counts", ["counts"], MC, threshold=2)
    assert all(out[c].dtype == np.uint64 for c in MC)
    assert _support(out[out.bc1 == "A"])[0] == n
    assert _support(out[out.bc1 == "B"])[0] == 1


def test_dominant_long_parent_takes_the_short_donor():
    long2 = "C" * 150 + LONG[-88:]
    base = pd.DataFrame({"bc1": ["A"] * 3, "bc0": ["CCCCCCCCCC"] * 3, "donor": [LONG, long2, SHORT]})
    dom = _with_masks(base.assign(counts=[100, 10, 5]), [0, 1, 2])
    out, st = attach_short_donors(dom, "donor", "counts", ["counts"], MC)          # 10x >= 3x -> attach
    assert len(out) == 2 and out[out.donor == LONG].counts.iloc[0] == 105 and st["attached_by_dominance"] == 1
    near = _with_masks(base.assign(counts=[100, 50, 5]), [0, 1, 2])
    out, st = attach_short_donors(near, "donor", "counts", ["counts"], MC)         # 2x < 3x -> ambiguous
    assert len(out) == 3 and out.donor_suffix_ambiguous.sum() == 1 and st["ambiguous"] == 1
    strict, st = attach_short_donors(dom, "donor", "counts", ["counts"], MC, dominance=float("inf"))
    assert len(strict) == 3 and st["ambiguous"] == 1


def test_nan_counts_are_refused_not_cast_to_int64_min():
    df = pd.DataFrame({"bc1": ["A"], "donor": ["d"], "bc0": ["x"], "counts": [np.nan]})
    with pytest.raises(ValueError):
        aggregate_links(df, ["bc1", "donor", "bc0"], "counts", ["counts"], [])
    with pytest.raises(ValueError):
        collapse_links(df, ["bc1"], "donor", "counts", ["counts"])
