"""A pool rebuild reports WHICH SEGMENT changed, not just how many rows did (NN 505 §4).

**The bug this pins.** The TRIUMPH rebuild changed 6,964 designs. Reported as one number
that looks catastrophic. Split by oligo segment it was two very different findings:
**6,603 differed only in the reverse priming site** — a relabelling, cheap and expected —
and **394 differed in the donor** — a real design change with a real cost. One pass over the
two pools separates them, and without that split a rebuild cannot be signed off.

**The contract.** :func:`diff_pools` joins on ``oligo_id`` and reports added / removed /
changed / unchanged, by category, plus for every changed row which declared segment
differs — with ``changed_only_in`` giving the NN 505 split directly.
"""

from __future__ import annotations

from magestic.design.pool import diff_pools

from . import _pool_fixtures as fx


def _swap_segment(seq: str, start: int, end: int | None, filler: str) -> str:
    """Overwrite one segment in place, keeping the sequence the same length."""
    start_i = start if start >= 0 else len(seq) + start
    end_i = len(seq) if end in (None, 0) else (end if end > 0 else len(seq) + end)
    width = end_i - start_i
    assert width > 0, (start, end, len(seq))
    return seq[:start_i] + (filler * width)[:width] + seq[end_i:]


def test_an_identical_pool_reports_no_change():
    old = fx.built_pool()
    new = [dict(r) for r in old]
    rep = diff_pools(old, new, fx.recipe())

    assert rep.changed == 0 and rep.added == 0 and rep.removed == 0
    assert rep.unchanged == 100


def test_the_nn505_split_is_recovered():
    """🔴 The decisive test: most rows differ only in the reverse priming site, a few only
    in the donor, and the report says so as two separate numbers."""
    old = fx.built_pool()
    new = [dict(r) for r in old]

    # 70 rows: reverse priming site only (the cheap relabelling)
    for r in new[:70]:
        r["oligo_seq"] = _swap_segment(r["oligo_seq"], -20, 0, "T")
    # 12 rows: donor only (the real design change)
    for r in new[70:82]:
        r["oligo_seq"] = _swap_segment(r["oligo_seq"], 51, -20, "A")

    rep = diff_pools(old, new, fx.recipe())

    assert rep.changed == 82
    assert rep.changed_only_in["reverse_priming"] == 70
    assert rep.changed_only_in["donor"] == 12
    assert rep.changed_in_multiple == 0
    assert rep.changed_unattributed == 0


def test_a_row_changing_in_two_segments_is_counted_separately():
    old = fx.built_pool()
    new = [dict(r) for r in old]
    r = new[0]
    seq = _swap_segment(r["oligo_seq"], -20, 0, "T")
    r["oligo_seq"] = _swap_segment(seq, 51, -20, "A")

    rep = diff_pools(old, new, fx.recipe())

    assert rep.changed == 1
    assert rep.changed_in_multiple == 1
    assert rep.changed_only_in == {}
    assert rep.changed_in["donor"] == 1 and rep.changed_in["reverse_priming"] == 1


def test_added_and_removed_rows_are_reported_by_category():
    old = fx.built_pool()
    new = [dict(r) for r in old if r["subpool"] != "syn_recoded"]
    for r in fx.make_rows("ptc", 3, seed=99):
        r["oligo_id"] = f"ptc_new_{r['oligo_id']}"
        r["subpool"] = "ptc"
        new.append(r)

    rep = diff_pools(old, new, fx.recipe())

    assert rep.removed == 5
    assert rep.added == 3
    assert rep.by_arm["syn_recoded"]["removed"] == 5
    assert rep.by_arm["ptc"]["added"] == 3
    assert rep.n_old == 100 and rep.n_new == 98


def test_a_whole_arm_disappearing_is_visible_in_the_per_arm_table():
    """The NN 505 ``syn_recoded`` loss, as the diff would have shown it."""
    old = fx.built_pool()
    new = [dict(r) for r in old if r["subpool"] != "syn_recoded"]

    rep = diff_pools(old, new, fx.recipe())

    assert rep.by_arm["syn_recoded"]["old"] == 5
    assert "new" not in rep.by_arm["syn_recoded"]
    assert rep.by_arm["syn_recoded"]["removed"] == 5
    assert "syn_recoded" in rep.summary()


def test_a_change_no_declared_segment_covers_is_reported_not_swallowed():
    """🔴 An unexplained change must be counted, never folded into 'unchanged'."""
    old = fx.built_pool()
    new = [dict(r) for r in old]
    new[0]["oligo_seq"] = new[0]["oligo_seq"] + "AAAA"  # length change

    rep = diff_pools(old, new, fx.recipe(), segments={"guide": (20, 40)})

    assert rep.changed == 1
    assert rep.changed_unattributed == 1


def test_the_diff_runs_without_a_recipe_for_a_legacy_pool():
    old = fx.built_pool()
    new = [dict(r) for r in old]
    for r in new[:5]:
        r["oligo_seq"] = _swap_segment(r["oligo_seq"], -20, 0, "T")

    rep = diff_pools(
        old, new,
        id_column="oligo_id", sequence_column="oligo_seq", arm_column="subpool",
        segments={"reverse_priming": (-20, 0)},
    )

    assert rep.changed == 5
    assert rep.changed_only_in["reverse_priming"] == 5


def test_the_summary_reads_as_two_numbers_with_different_meanings():
    old = fx.built_pool()
    new = [dict(r) for r in old]
    for r in new[:70]:
        r["oligo_seq"] = _swap_segment(r["oligo_seq"], -20, 0, "T")
    for r in new[70:82]:
        r["oligo_seq"] = _swap_segment(r["oligo_seq"], 51, -20, "A")

    text = diff_pools(old, new, fx.recipe()).summary()

    assert "changed in exactly one segment" in text
    assert "reverse_priming: 70" in text
    assert "donor: 12" in text
