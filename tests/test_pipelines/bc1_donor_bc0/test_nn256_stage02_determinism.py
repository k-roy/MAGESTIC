"""NN 256 -- bc1 stage 02 must not depend on the order of its input manifest.

NN 253 S8 measured the defect: the same 18 samples with a SHUFFLED manifest moved ~14,000 of
1,630,592 links (0.43%; 7,087 only in run a, 6,994 only in run b) while the net row count moved by
only 93. A row-count check therefore hides it, and so does a same-order rerun -- two runs of one
order were already byte-identical. **The control that matters is a DIFFERENT order.**

Cause: ``collapse_links`` / ``aggregate_links`` sorted stably on keys that are not TOTAL, so the
tie-winner was whichever row arrived first, which is exactly what shuffling the manifest changes.

Every test here RUNS the code (a source-text assertion would have passed while the thing it named
was unbound -- that happened earlier in the same session). Each one is a known-positive: it FAILS on
the pre-NN-256 module. Which hole each covers is named in its docstring, because a determinism test
that passes on both the broken and the fixed tree is testing nothing.
"""
import numpy as np
import pandas as pd
import pytest

from magestic.pipelines.bc1_donor_bc0.core import link_collapse as LC
from magestic.pipelines.bc1_donor_bc0.core.link_collapse import (
    aggregate_links, attach_short_donors, collapse_links, mask_columns, masks_for_codes,
    merge_orphan_bc1)

N_IDS = 70                        # spans two uint64 words on purpose
MC = mask_columns(N_IDS)
SUM = ["counts"]


def _with_masks(df, codes):
    m = masks_for_codes(np.asarray(codes), N_IDS)
    for k, c in enumerate(MC):
        df[c] = m[:, k]
    return df


def _canon(df, keys=("bc1", "donor", "bc0")):
    """Row-order-free view of a frame: sort on a TOTAL key, drop the index, order the columns."""
    k = [c for c in keys if c in df.columns]
    out = df.sort_values(k, kind="mergesort").reset_index(drop=True)
    return out[sorted(out.columns)]


def _perm(df, seed):
    return df.sample(frac=1.0, random_state=seed).reset_index(drop=True)


# ===========================================================================
# A. the ordering-key primitive
# ===========================================================================

def test_order_key_is_injective_where_str_is_not():
    """`str()` maps None, NaN and "" all to "", and 0 and "0" both to "0". A sort keyed on it
    therefore still lets the INPUT ROW ORDER decide between genuinely different values."""
    vals = [None, float("nan"), "", 0, "0", 1, "1", True, "True"]
    keys = [LC._order_key(v) for v in vals]
    # the two nulls are one class by design; everything else must be distinguishable
    non_null = [k for v, k in zip(vals, keys) if not (v is None or (isinstance(v, float) and v != v))]
    assert len(set(non_null)) == len(non_null), dict(zip(map(repr, vals), keys))
    # a missing value must never WIN a tie against a present one
    assert LC._order_key(None) > LC._order_key("")
    assert LC._order_key(None) > LC._order_key("zzz")
    # among plain strings the order is unchanged, so no non-tie behaviour moves
    assert (LC._order_key("abc") < LC._order_key("abd")) is ("abc" < "abd")


def test_order_key_survives_a_pandas_multicolumn_sort():
    """THE GUARD FOR THE TRAP THAT BIT THE FIRST VERSION OF THIS FIX.

    pandas' MULTI-column `sort_values` over object strings silently does NOT order them when the
    strings contain a NUL byte: it goes through lexsort_indexer -> Categorical -> factorize, whose
    hash table truncates at the NUL, so every key collapses into one category and the sort returns
    the INPUT ROW ORDER. A single-column sort on the same data orders correctly, which is why it is
    so easy to miss. The first NN 256 `_order_key` used "\x01s\x00" + v and so reintroduced, inside
    its own fix, the order dependence it was written to remove.

    Two assertions: the key must contain no NUL, and a pandas two-column sort keyed on it must
    agree with Python's own ordering.
    """
    raw = ["s9_0", "s5_0", "s0_0", "s10_0", "", "zz"]
    keys = [LC._order_key(v) for v in raw]
    assert not any("\x00" in k for k in keys), "an ordering key must not contain a NUL byte"
    want = [r for _, r in sorted(zip(keys, raw))]
    df = pd.DataFrame({"_seg": [0] * len(raw), "_k0": keys, "raw": raw})
    got = df.sort_values(["_seg", "_k0"], kind="mergesort")["raw"].tolist()
    assert got == want, f"pandas multi-column sort disagrees with Python: {got} != {want}"
    # and the vectorized path the module actually uses
    codes = LC._rank_codes(np.array(keys, dtype=object))
    assert [raw[i] for i in np.lexsort((codes,))] == want


def test_reserved_column_is_rejected():
    df = pd.DataFrame({"bc1": ["A"], "donor": ["d"], "bc0": ["x"], "counts": [1],
                       LC.REP_RANK_COL: [1]})
    with pytest.raises(ValueError):
        aggregate_links(df, ["bc1", "donor", "bc0"], "counts", ["counts", LC.REP_RANK_COL], [])


# ===========================================================================
# B. aggregate_links -- the fold's preserve pick must be ASSOCIATIVE
# ===========================================================================
#
# KNOWN-POSITIVE for the NN 256 `_reprank` fix. In 02 the cross-sample fold calls aggregate_links
# once per batch of `--accumulate-every` samples; the accumulated row's `counts` is a SUM, so a
# representative ranked on `counts` can be beaten by an accumulation of smaller rows. The batch
# boundaries follow the MANIFEST ROW ORDER, so the preserved columns -- SPS, donor_bc0_fragment,
# trafo_ID, V_library -- were manifest-order-dependent on any manifest longer than one batch.
# The 18-sample / default-24 end-to-end test cannot see this: it does exactly one fold.

def _fold(frames, every, key, preserve):
    """02's incremental fold, verbatim in structure."""
    acc, pending = None, []
    for f in frames:
        pending.append(f)
        if len(pending) >= every:
            acc = aggregate_links(pd.concat([*( [acc] if acc is not None else [] ), *pending],
                                            ignore_index=True),
                                  key, "counts", SUM, MC, preserve)
            pending = []
    if pending:
        acc = aggregate_links(pd.concat([*( [acc] if acc is not None else [] ), *pending],
                                        ignore_index=True), key, "counts", SUM, MC, preserve)
    return acc


def test_fold_preserve_pick_is_batch_size_invariant_minimal():
    """The 3-row counterexample. One key; per-sample rows (3,"Y") (3,"Y") (5,"X").

    One pass -> max rank 5 -> "X". Folded two-at-a-time -> the accumulator becomes (6,"Y") and
    beats (5,"X") -> "Y". Both must now say "X": the pick ranks on the max INDIVIDUAL rank.
    """
    rows = [("Y", 3), ("Y", 3), ("X", 5)]
    frames = [_with_masks(pd.DataFrame({"bc1": ["A"], "donor": ["d"], "bc0": ["x"],
                                        "counts": [c], "SPS": [s]}), [i])
              for i, (s, c) in enumerate(rows)]
    key, pres = ["bc1", "donor", "bc0"], ["SPS"]
    one = _fold(frames, every=99, key=key, preserve=pres)
    two = _fold(frames, every=2, key=key, preserve=pres)
    assert one["counts"].iloc[0] == two["counts"].iloc[0] == 11
    assert one["SPS"].iloc[0] == "X", "one-pass pick changed"
    assert two["SPS"].iloc[0] == one["SPS"].iloc[0], (
        "the fold's preserve pick depends on the batch size, and batch boundaries follow the "
        "manifest order")


@pytest.mark.parametrize("every", [1, 2, 3, 5, 99])
@pytest.mark.parametrize("seed", [0, 1, 7, 42])
def test_fold_is_invariant_to_batch_size_and_frame_order(every, seed):
    """Same property on a wider fixture, and over the frame ORDER as well as the batch size."""
    rng = np.random.default_rng(0)
    frames = []
    for i in range(12):
        n = 6
        frames.append(_with_masks(pd.DataFrame({
            "bc1": ["A" * 20, "C" * 20, "G" * 20] * 2,
            "donor": ["ACGT" * 10, "ACGT" * 10, "TTTT" * 10] * 2,
            "bc0": ["x" * 8, "x" * 8, "y" * 8] * 2,
            "counts": rng.integers(1, 6, n).tolist(),
            "SPS": [f"s{i}_{j}" for j in range(n)],
            "V_library": [f"V{j}" for j in range(n)],
        }), [i] * n))
    key, pres = ["bc1", "donor", "bc0"], ["SPS", "V_library"]
    ref = _fold(frames, every=99, key=key, preserve=pres)
    shuffled = list(frames)
    np.random.default_rng(seed).shuffle(shuffled)
    got = _fold(shuffled, every=every, key=key, preserve=pres)
    pd.testing.assert_frame_equal(_canon(ref.drop(columns=["_reprank"], errors="ignore")),
                                  _canon(got.drop(columns=["_reprank"], errors="ignore")))


def test_aggregate_preserve_tie_is_total_not_input_order():
    """KNOWN-POSITIVE for the type-tagged key. The old key was `str(v)`, which mapped None and ""
    to the same string, so which raw value was emitted fell back to the input row order."""
    base = _with_masks(pd.DataFrame({
        "bc1": ["A"] * 4, "donor": ["d"] * 4, "bc0": ["x"] * 4,
        "counts": [9, 9, 9, 9], "SPS": ["", None, float("nan"), ""]}), [0, 1, 2, 3])
    outs = set()
    for seed in range(12):
        o = aggregate_links(_perm(base, seed), ["bc1", "donor", "bc0"], "counts", SUM, MC, ["SPS"])
        outs.add(repr(o["SPS"].iloc[0]))
    assert len(outs) == 1, f"the preserved value depends on the input row order: {outs}"


# ===========================================================================
# C. collapse_links -- _peel's top selection must be TOTAL
# ===========================================================================
#
# KNOWN-POSITIVE for the `_peel` fix and the live cause of the measured 0.43% churn. `_peel` took
# "the first row whose rank equals the maximum", i.e. the first in INPUT ROW ORDER, and the
# representative's seq_col and preserve_cols are exactly what the function emits.

def _tie_group(n=24, count=7):
    """One (bc1, bc0) group, n donors ALL TIED on counts and all within one edit of each other,
    so every one of them is a legitimate representative and only the tie-break decides."""
    stem = "ACGTACGTACGTACGTACGT"
    donors = [stem[:19] + b for b in "ACGT"][:4]
    donors = [stem[:19] + c for c in "ACGT"] + [stem[:18] + c + "T" for c in "ACGT"]
    donors = (donors * ((n // len(donors)) + 1))[:n]
    donors = [d + f"{i:03d}" for i, d in enumerate(donors)]       # keep them distinct, still <=1+ edits apart
    return _with_masks(pd.DataFrame({
        "bc1": ["A" * 20] * n, "bc0": ["x" * 8] * n, "donor": donors,
        "counts": [count] * n, "SPS": [f"s{i:03d}" for i in range(n)]}),
        list(range(n)))


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 7, 13, 42, 99, 1234, 31337])
def test_collapse_links_is_permutation_invariant(seed):
    """n=24 tied donors: well past numpy's 16-element insertion-sort fallback, so the real sort
    path is exercised (the trap that produced a false 'NOT REPRODUCED' in NN 250 D9)."""
    base = _tie_group(n=24)
    a = collapse_links(base, ["bc1", "bc0"], "donor", "counts", SUM, MC, ["SPS"], threshold=5)
    b = collapse_links(_perm(base, seed), ["bc1", "bc0"], "donor", "counts", SUM, MC, ["SPS"],
                       threshold=5)
    pd.testing.assert_frame_equal(_canon(a), _canon(b))


@pytest.mark.parametrize("n_jobs", [1, 2])
def test_collapse_links_is_permutation_invariant_parallel(n_jobs):
    """Same property with the ProcessPoolExecutor branch forced."""
    base = pd.concat([_tie_group(n=20).assign(bc1=b * 20) for b in "ACGT"], ignore_index=True)
    a = collapse_links(base, ["bc1", "bc0"], "donor", "counts", SUM, MC, ["SPS"], threshold=5,
                       n_jobs=n_jobs)
    b = collapse_links(_perm(base, 5), ["bc1", "bc0"], "donor", "counts", SUM, MC, ["SPS"],
                       threshold=5, n_jobs=n_jobs)
    pd.testing.assert_frame_equal(_canon(a), _canon(b))


# ===========================================================================
# D. attach_short_donors -- the two sorted() calls
# ===========================================================================
#
# These two ties are LATENT at the default --attach-dominance 3.0: a tie at the top of `cands`
# fails `rank >= 3 * rank` either way, so no attach happens and the order is unobservable. They
# become live at --attach-dominance <= 1, which the CLI exposes. Tested there, honestly labelled:
# this hole is real and now closed, but it is NOT what produced NN 253's 0.43% churn.

def _suffix_tie(n_parents=6, parent_len=60, short_len=30, rank=10):
    parents, rows = [], []
    for i in range(n_parents):
        p = ("ACGT" * 20)[:parent_len - 3] + f"{i:03d}"
        p = f"{i:03d}" + ("ACGT" * 20)[:parent_len - 3]      # differ at the 5' end only
        parents.append(p)
        rows.append({"bc1": "A" * 20, "donor": p, "bc0": "x" * 8, "counts": rank})
    short = parents[0][-short_len:]                           # a suffix of EVERY parent
    rows.append({"bc1": "A" * 20, "donor": short, "bc0": "x" * 8, "counts": 1})
    return _with_masks(pd.DataFrame(rows), list(range(len(rows))))


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 7, 13, 42, 99])
def test_attach_short_donors_is_permutation_invariant_at_dominance_1(seed):
    """KNOWN-POSITIVE for the terminal keys on link_collapse :368 and :384. Several equal-length,
    equal-rank parents all qualify for one short donor; at dominance 1.0 the tie ATTACHES, so
    which parent `top2[0]` names is observable -- and it was the input row order that named it."""
    base = _suffix_tie()
    a, sa = attach_short_donors(base, "donor", "counts", SUM, MC, dominance=1.0)
    b, sb = attach_short_donors(_perm(base, seed), "donor", "counts", SUM, MC, dominance=1.0)
    assert sa == sb, f"attach stats depend on input row order: {sa} vs {sb}"
    pd.testing.assert_frame_equal(_canon(a), _canon(b))


@pytest.mark.parametrize("seed", [0, 3, 42])
def test_attach_short_donors_is_permutation_invariant_at_default_dominance(seed):
    base = _suffix_tie()
    a, sa = attach_short_donors(base, "donor", "counts", SUM, MC)
    b, sb = attach_short_donors(_perm(base, seed), "donor", "counts", SUM, MC)
    assert sa == sb
    pd.testing.assert_frame_equal(_canon(a), _canon(b))


# ===========================================================================
# E. merge_orphan_bc1 -- :571, on the FULL FRAME, not just the stats
# ===========================================================================
#
# NN 250 D9 already tested the stats and the bc1->root map. The stats can match while the emitted
# rows differ, which is the whole lesson of NN 253 (equal counts, different link set), so this
# compares the frame.

def _orphan_fixture(n_tied=40):
    root, orphan = "A" * 20, "A" * 19 + "C"
    rows = []
    for i in range(n_tied):
        d = f"{'ACGT' * 10}{i:03d}"
        rows.append({"bc1": root, "donor": d, "bc0": f"BC{i:08d}", "counts": 7})
        rows.append({"bc1": orphan, "donor": d, "bc0": f"BC{i:08d}", "counts": 7})
    rows.append({"bc1": root, "donor": "ACGT" * 10 + "999", "bc0": "BC99999999", "counts": 9})
    return _with_masks(pd.DataFrame(rows), [0] * len(rows))


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 7, 13, 42, 99, 1234, 31337])
def test_merge_orphan_bc1_full_frame_is_permutation_invariant(seed):
    base = _orphan_fixture()
    a, sa, ma = merge_orphan_bc1(base.copy(), "counts", SUM, MC, ["SPS"] if "SPS" in base else [])
    b, sb, mb = merge_orphan_bc1(_perm(base, seed), "counts", SUM, MC,
                                 ["SPS"] if "SPS" in base else [])
    assert sa == sb, f"orphan stats depend on input row order (seed={seed})"
    assert dict(zip(ma["bc1"], ma["bc1_root"])) == dict(zip(mb["bc1"], mb["bc1_root"]))
    pd.testing.assert_frame_equal(_canon(a), _canon(b))


# ===========================================================================
# F. the whole stage-02 collapse chain, in three orders
# ===========================================================================
#
# The unit-level mirror of the end-to-end acceptance test: the same calls 02 makes, in the same
# order, over a fixture built to be tie-rich at every step. Three input orders; the emitted frame
# must be identical. A same-order rerun is deliberately NOT the control -- it was already identical
# before the fix.

def _chain_fixture():
    rng = np.random.default_rng(11)
    rows = []
    bc1s = ["A" * 20, "A" * 19 + "C", "C" * 20, "G" * 20, "T" * 20]   # first two are Hamming 1
    for bi, bc1 in enumerate(bc1s):
        for g in range(5):
            parent = f"{g:03d}" + ("ACGT" * 20)[:57]                  # 60 nt
            variants = [parent] + [parent[:k] + "A" + parent[k + 1:] for k in (10, 20)]
            variants += [parent[-30:]]                                # the 2x150-era suffix read
            for v in variants:
                for bc0 in ("x" * 8, "x" * 7 + "y"):                  # within one edit
                    rows.append({"bc1": bc1, "donor": v, "bc0": bc0,
                                 "counts": int(rng.integers(3, 9)),   # many exact ties
                                 "SPS": f"sps_{bi}{g}", "V_library": f"V{g}"})
    df = pd.DataFrame(rows)
    return _with_masks(df, rng.integers(0, N_IDS, len(df)).tolist())


def _chain(df, accumulate_every=None):
    key, pres = ["bc1", "donor", "bc0"], ["SPS", "V_library"]
    if accumulate_every:                                 # exercise the fold, as 02 does
        parts = [df.iloc[i:i + 7] for i in range(0, len(df), 7)]
        acc = _fold(parts, accumulate_every, key, pres)
    else:
        acc = aggregate_links(df, key, "counts", SUM, MC, pres)
    acc = collapse_links(acc, ["bc1", "bc0"], "donor", "counts", SUM, MC, pres, threshold=5)
    acc = collapse_links(acc, ["bc1", "donor"], "bc0", "counts", SUM, MC, pres, threshold=1)
    acc, st = attach_short_donors(acc, "donor", "counts", SUM, MC, max_mismatch=5, bc0_max_edit=1)
    px = pres + ["donor_suffix_ambiguous", "donor_lengths_observed"]
    acc = collapse_links(acc, ["bc1", "bc0"], "donor", "counts", SUM, MC, px, threshold=5)
    acc = collapse_links(acc, ["donor", "bc0"], "bc1", "counts", SUM, MC, px, threshold=1)
    acc = collapse_links(acc, ["bc1", "donor"], "bc0", "counts", SUM, MC, px, threshold=1)
    acc, ost, _m = merge_orphan_bc1(acc, "counts", SUM, MC, px, donor_max_mismatch=5)
    return acc, st, ost


@pytest.mark.parametrize("seed", [1, 2, 3, 7, 42, 99])
def test_full_chain_is_permutation_invariant(seed):
    base = _chain_fixture()
    a, sta, osa = _chain(base)
    b, stb, osb = _chain(_perm(base, seed))
    assert sta == stb, f"suffix-attach stats depend on row order: {sta} vs {stb}"
    assert osa == osb, "orphan-merge stats depend on row order"
    pd.testing.assert_frame_equal(_canon(a.drop(columns=["_reprank"], errors="ignore")),
                                  _canon(b.drop(columns=["_reprank"], errors="ignore")))


@pytest.mark.parametrize("every", [1, 3, 8])
def test_full_chain_is_invariant_to_the_fold_batch_size(every):
    """Same chain, with the incremental fold at three batch sizes. 02's batch boundaries are a
    function of the manifest order, so a result that moves with the batch size is a result that
    moves with the manifest."""
    base = _chain_fixture()
    ref, _s, _o = _chain(base)
    got, _s2, _o2 = _chain(base, accumulate_every=every)
    pd.testing.assert_frame_equal(_canon(ref.drop(columns=["_reprank"], errors="ignore")),
                                  _canon(got.drop(columns=["_reprank"], errors="ignore")))


def test_the_comparison_would_detect_a_real_difference():
    """The positive control for the comparison itself: `_canon` + assert_frame_equal must FAIL on a
    pair that genuinely differs. A determinism assertion is meaningless without this -- a detector
    returning 0 may simply have been disabled."""
    base = _chain_fixture()
    with_merge, _s, _o = _chain(base)
    no_merge = _chain(base)[0]
    key, pres = ["bc1", "donor", "bc0"], ["SPS", "V_library"]
    acc = aggregate_links(base, key, "counts", SUM, MC, pres)
    acc = collapse_links(acc, ["bc1", "bc0"], "donor", "counts", SUM, MC, pres, threshold=5)
    assert len(acc) != len(with_merge), "fixture is not discriminating"
    with pytest.raises(AssertionError):
        pd.testing.assert_frame_equal(
            _canon(acc.drop(columns=["_reprank"], errors="ignore")),
            _canon(with_merge.drop(columns=["_reprank"], errors="ignore")))
