"""The orphan-bc1 merge and its star-only guard (Kevin 2026-09-10, decision C; planning 33_).

The three exact-key collapses in step 02 group by (donor, bc0), so a bc1 read error whose donor read
ALSO carries an error never lands in its parent's group and survives as its own barcode. This merges it,
inside step 02 rather than in a post-hoc script against the shipped table.
"""
import numpy as np
import pandas as pd
import pytest

from magestic.pipelines.bc1_donor_bc0.core.link_collapse import (
    merge_orphan_bc1, _bc1_hamming1_pairs, _donor_compatible, _edit_le1)

SUM = ["counts", "n_a", "n_b"]
OR = ["_repmask_0"]
PRESERVE = ["SPS"]
D1 = "ACGT" * 32 + "ACGTA"      # 133 nt
D2 = D1[:66] + "T" + D1[67:]    # 1 mismatch from D1


def _row(bc1, donor, bc0, counts, mask=1, n_a=None, sps="S"):
    return {"bc1": bc1, "donor": donor, "bc0": bc0, "counts": counts,
            "n_a": counts if n_a is None else n_a, "n_b": 0, "_repmask_0": np.uint64(mask), "SPS": sps}


def _df(rows):
    d = pd.DataFrame(rows)
    d["_repmask_0"] = d["_repmask_0"].astype("uint64")
    return d


def test_hamming1_pair_finder():
    bc = np.array(["A" * 20, "C" + "A" * 19, "CC" + "A" * 18, "G" * 20])
    pairs = {tuple(p) for p in _bc1_hamming1_pairs(bc)}
    assert pairs == {(0, 1), (1, 2)}          # 0 vs 2 is Hamming 2, 3 is far from everything
    assert len(_bc1_hamming1_pairs(np.array(["A" * 20]))) == 0


def test_edit_and_donor_predicates():
    assert _edit_le1("ACGTACGTAC", "ACGTACGTAC") and _edit_le1("ACGTACGTAC", "ACGTTCGTAC")
    assert _edit_le1("ACGTACGTAC", "ACGTACGTA") and not _edit_le1("ACGTACGTAC", "TTTTTTTTTT")
    assert _donor_compatible(D1, D1, 5) and _donor_compatible(D1, D2, 5)
    assert _donor_compatible(D1[-88:], D1, 5)            # a 2x150 read of a 2x300 parent
    assert not _donor_compatible(D1, "G" * 133, 5)
    assert not _donor_compatible(D1, None, 5) and not _donor_compatible(D1, "", 5)


def test_orphan_is_merged_counts_summed_and_replicates_unioned():
    parent, orphan = "A" * 20, "C" + "A" * 19
    df = _df([_row(parent, D1, "ACGTACGTAC", 100, mask=0b01),
              _row(orphan, D2, "ACGTACGTAG", 5, mask=0b10)])
    out, stats, m = merge_orphan_bc1(df, "counts", SUM, OR, PRESERVE)
    assert stats["accepted_pairs"] == 1 and stats["bc1_removed"] == 1 and stats["reads_moved"] == 5
    assert out["bc1"].nunique() == 1 and out["bc1"].iloc[0] == parent   # root = higher counts
    assert int(out["counts"].sum()) == 105 and int(out["n_a"].sum()) == 105
    # replicate membership is a UNION, never a sum: two distinct replicates -> both bits set
    assert int(out["_repmask_0"].max()) == 0b11
    assert m.to_dict("records") == [{"bc1": orphan, "bc1_root": parent}]
    assert stats["reads_in"] == stats["reads_out"] == 105


def test_incompatible_donor_or_bc0_is_not_merged():
    a, b = "A" * 20, "C" + "A" * 19
    far_donor = _df([_row(a, D1, "ACGTACGTAC", 100), _row(b, "G" * 133, "ACGTACGTAC", 5)])
    out, stats, m = merge_orphan_bc1(far_donor, "counts", SUM, OR, PRESERVE)
    assert stats["accepted_pairs"] == 0 and stats["bc1_removed"] == 0 and m.empty
    far_bc0 = _df([_row(a, D1, "ACGTACGTAC", 100), _row(b, D1, "TTTTTTTTTT", 5)])
    out, stats, m = merge_orphan_bc1(far_bc0, "counts", SUM, OR, PRESERVE)
    assert stats["accepted_pairs"] == 0 and stats["bc1_removed"] == 0


def test_star_guard_refuses_a_transitive_chain():
    """x - y - z: x and z are Hamming 2 apart, so only y has direct evidence against the root. Naive
    connected components would swallow z as well; that is how components spanning distance 8 arose."""
    x, y, z = "A" * 20, "C" + "A" * 19, "CC" + "A" * 18
    df = _df([_row(x, D1, "ACGTACGTAC", 100), _row(y, D1, "ACGTACGTAC", 10),
              _row(z, D1, "ACGTACGTAC", 3)])
    guarded, gs, gm = merge_orphan_bc1(df, "counts", SUM, OR, PRESERVE, star_guard=True)
    assert gs["accepted_pairs"] == 2
    assert set(guarded["bc1"]) == {x, z}          # y merged into x; z kept, it is Hamming 2 from x
    assert gm.to_dict("records") == [{"bc1": y, "bc1_root": x}]
    assert int(guarded.loc[guarded["bc1"] == x, "counts"].iloc[0]) == 110

    unguarded, us, um = merge_orphan_bc1(df, "counts", SUM, OR, PRESERVE, star_guard=False)
    assert set(unguarded["bc1"]) == {x}           # the chain swallows z
    assert int(unguarded["counts"].sum()) == 113
    assert us["bc1_removed"] == 2


def test_root_is_the_highest_count_member_with_a_deterministic_tie_break():
    hi, lo = "T" + "A" * 19, "A" * 20
    df = _df([_row(lo, D1, "ACGTACGTAC", 5), _row(hi, D1, "ACGTACGTAC", 500)])
    out, _s, m = merge_orphan_bc1(df, "counts", SUM, OR, PRESERVE)
    assert out["bc1"].tolist() == [hi] and m["bc1_root"].tolist() == [hi]
    tie = _df([_row("A" * 20, D1, "ACGTACGTAC", 7), _row("C" + "A" * 19, D1, "ACGTACGTAC", 7)])
    out2, _s2, _m2 = merge_orphan_bc1(tie, "counts", SUM, OR, PRESERVE)
    assert out2["bc1"].tolist() == ["A" * 20]     # lexicographic tie-break, reproducible


def test_reads_are_conserved_and_nothing_is_dropped():
    rows, rng = [], np.random.default_rng(11)
    for i in range(40):
        bc = "".join(rng.choice(list("ACGT"), 20))
        rows.append(_row(bc, D1, "ACGTACGTAC", int(rng.integers(1, 500))))
        if i % 4 == 0:                                    # give some of them a 1-nt neighbour
            nb = ("T" if bc[0] != "T" else "G") + bc[1:]
            rows.append(_row(nb, D2, "ACGTACGTAC", int(rng.integers(1, 20))))
    df = _df(rows)
    out, stats, m = merge_orphan_bc1(df, "counts", SUM, OR, PRESERVE)
    assert stats["reads_in"] == stats["reads_out"] == int(df["counts"].sum())
    assert stats["bc1_removed"] == len(m)
    assert out["bc1"].nunique() == df["bc1"].nunique() - stats["bc1_removed"]
    assert set(m["bc1"]).isdisjoint(set(out["bc1"]))      # a remapped bc1 never survives


def test_a_second_unrelated_link_of_an_orphan_is_re_keyed_but_not_force_merged():
    """Only the row that carries the accepted evidence is snapped onto the root's link. Another link of
    the same barcode is re-keyed to the root bc1 and keeps its own donor and bc0."""
    parent, orphan = "A" * 20, "C" + "A" * 19
    other = "T" * 133
    df = _df([_row(parent, D1, "ACGTACGTAC", 100),
              _row(orphan, D2, "ACGTACGTAG", 5),
              _row(orphan, other, "GGGGGGGGGG", 9)])
    out, stats, m = merge_orphan_bc1(df, "counts", SUM, OR, PRESERVE)
    assert stats["rows_snapped_to_root_link"] == 1
    assert set(out["bc1"]) == {parent}
    assert len(out) == 2
    merged_row = out[out["donor"] == D1]
    assert int(merged_row["counts"].iloc[0]) == 105
    kept = out[out["donor"] == other]
    assert int(kept["counts"].iloc[0]) == 9 and kept["bc0"].iloc[0] == "GGGGGGGGGG"
    assert stats["reads_in"] == stats["reads_out"] == 114


def test_dirty_and_empty_inputs():
    empty = _df([_row("A" * 20, D1, "ACGTACGTAC", 1)]).iloc[0:0]
    out, stats, m = merge_orphan_bc1(empty, "counts", SUM, OR, PRESERVE)
    assert len(out) == 0 and stats["accepted_pairs"] == 0 and m.empty
    dirty = _df([_row("A" * 19, D1, "ACGTACGTAC", 5), _row("NNNN" + "A" * 16, D1, "ACGTACGTAC", 5)])
    out, stats, m = merge_orphan_bc1(dirty, "counts", SUM, OR, PRESERVE)
    assert stats["bc1_clean_20nt"] == 0 and len(out) == 2   # non-20-nt ACGT barcodes are left alone
