"""Non-targeting controls match the GC DISTRIBUTION and come out at exactly N (NN 496/497).

**The bugs these pin**, both paid for by failed builds of the real pool:

1. 🔴 **Distribution, not mean.** A mean-matched scramble set passes a test of means and
   fails a two-sample KS test against the real pool. The real ``gc_frac`` is binned at 2 %
   and scrambles are drawn per bin in proportion, which is what Kevin's wording asked for.

2. 🔴 **Exactly N needs a top-up *and* a trim.** The per-bin targets round independently, so
   they sum to one or two **less** than N — every early build landed short and was chased by
   hand — *and* sometimes one **more**: the v5 pool's GC histogram differed from v4's by six
   designs, rounded to 1001, and the job died with the whole pack discarded. Both
   corrections draw from the modal GC bin, where the real designs are densest, so the
   distribution match is perturbed least.

**The contract.** ``make_nontargeting_controls`` returns exactly *n* rows or raises. A short
arm is never silently shipped, and the KS statistic against the reference is reported on
every build.
"""

from __future__ import annotations

import pytest

from magestic.design.pool import ks_statistic, make_nontargeting_controls
from magestic.design.pool.scramble import gc_fraction

from . import _pool_fixtures as fx


def _reference_gc(n: int = 200, seed: int = 11) -> list[float]:
    """GC fractions of *n* real-looking 200-mers."""
    rows = fx.make_rows("natural", n, seed=seed)
    return [gc_fraction(r["oligo_seq"]) for r in rows]


def test_exactly_n_controls_are_made():
    rows, report = make_nontargeting_controls(40, _reference_gc(), fx.ARCHITECTURE)
    assert len(rows) == 40
    assert report["requested"] == report["made"] == 40


def test_ids_are_unique_and_contiguous_after_any_correction():
    rows, _ = make_nontargeting_controls(25, _reference_gc(), fx.ARCHITECTURE)
    ids = [r["oligo_id"] for r in rows]
    assert len(set(ids)) == 25
    assert ids[0] == "scrambled_00001"
    assert ids[-1] == "scrambled_00025"


def test_every_control_has_the_designed_architecture_and_one_cloning_site():
    rows, _ = make_nontargeting_controls(20, _reference_gc(), fx.ARCHITECTURE)
    for r in rows:
        seq = r["oligo_seq"]
        assert len(seq) == 200
        assert seq.startswith(fx.FORWARD_PRIMING)
        assert seq.endswith(fx.REVERSE_TAIL)
        assert sum(seq.count(s) for s in ("GCTCTTC", "GAAGAGC")) == 1
        assert set(seq) <= set("ACGT")


def test_the_gc_distribution_is_matched_not_just_the_mean():
    """🔴 The decisive check: a small KS statistic against the real distribution."""
    reference = _reference_gc(300)
    rows, report = make_nontargeting_controls(150, reference, fx.ARCHITECTURE)

    made = [r["gc_frac"] for r in rows]
    assert report["ks_statistic_vs_reference"] == pytest.approx(
        ks_statistic(reference, made), abs=1e-5
    )
    assert report["ks_statistic_vs_reference"] < 0.15, (
        f"KS {report['ks_statistic_vs_reference']} — the scrambles do not track the real "
        f"GC distribution"
    )


def test_a_mean_matched_set_would_fail_this_check():
    """The positive control for the test above: matching the mean is not enough.

    Every oligo drawn at the reference MEAN gives a far tighter GC spread, which a test of
    means passes and the KS test rejects. Without this, "KS < 0.15" is not evidence that the
    distribution match is doing any work.
    """
    reference = _reference_gc(300)
    mean_gc = sum(reference) / len(reference)
    narrow = [mean_gc] * 150

    assert abs(sum(narrow) / len(narrow) - mean_gc) < 1e-9, "means agree by construction"
    assert ks_statistic(reference, narrow) > 0.4, (
        "a mean-matched set must be rejected by the KS statistic, or the distribution "
        "match is not being tested"
    )


def test_the_top_up_path_is_exercised_and_still_gives_exactly_n():
    """Bins rounding DOWN: 0.5/0.25/0.25 at n=10 rounds to 5+2+2 = 9, one short."""
    reference = [0.41] * 50 + [0.45] * 25 + [0.49] * 25
    rows, report = make_nontargeting_controls(10, reference, fx.ARCHITECTURE)

    assert report["topped_up"] >= 1, "this reference must undershoot before the top-up"
    assert report["trimmed"] == 0
    assert len(rows) == 10


def test_the_trim_path_is_exercised_and_still_gives_exactly_n():
    """Bins rounding UP: three equal bins at n=11 round to 4+4+4 = 12, one over.

    This is the mirror case that killed a real build after the top-up had been added.
    """
    reference = [0.41] * 30 + [0.45] * 30 + [0.49] * 30
    rows, report = make_nontargeting_controls(11, reference, fx.ARCHITECTURE)

    assert report["trimmed"] >= 1, "this reference must overshoot before the trim"
    assert len(rows) == 11
    ids = [r["oligo_id"] for r in rows]
    assert len(set(ids)) == 11, "ids must be renumbered contiguously after a trim"
    assert ids[-1] == "scrambled_00011"


def test_a_rebuild_from_the_same_seed_reproduces_the_same_controls():
    reference = _reference_gc()
    a, _ = make_nontargeting_controls(20, reference, fx.ARCHITECTURE, seed=4242)
    b, _ = make_nontargeting_controls(20, reference, fx.ARCHITECTURE, seed=4242)
    assert [r["oligo_seq"] for r in a] == [r["oligo_seq"] for r in b]


def test_a_different_seed_gives_different_controls():
    reference = _reference_gc()
    a, _ = make_nontargeting_controls(20, reference, fx.ARCHITECTURE, seed=1)
    b, _ = make_nontargeting_controls(20, reference, fx.ARCHITECTURE, seed=2)
    assert [r["oligo_seq"] for r in a] != [r["oligo_seq"] for r in b]


def test_zero_controls_is_a_valid_empty_arm():
    rows, report = make_nontargeting_controls(0, _reference_gc(), fx.ARCHITECTURE)
    assert rows == []
    assert report["made"] == 0


def test_an_empty_reference_raises_rather_than_inventing_a_distribution():
    with pytest.raises(ValueError, match="no distribution to match"):
        make_nontargeting_controls(10, [], fx.ARCHITECTURE)


def test_an_impossible_architecture_raises():
    from magestic.design.pool import OligoArchitecture

    too_long = OligoArchitecture(
        forward_priming="A" * 100,
        internal_cloning_site="C" * 60,
        reverse_tail="G" * 60,
        total_length=200,
    )
    with pytest.raises(ValueError, match="no variable region"):
        make_nontargeting_controls(1, [0.5], too_long)


def test_an_unreachable_target_raises_rather_than_returning_a_short_arm():
    """🔴 A short arm is never silently shipped."""
    from magestic.design.pool import OligoArchitecture

    # a GC bin the architecture can never reach: the fixed parts alone overshoot it
    arch = OligoArchitecture(
        forward_priming="GCGCGCGCGCGCGCGCGCGC",
        internal_cloning_site="GTTTGAAGAGC",
        reverse_tail="GCGCGCGCGCGCGCGCGCGC",
        total_length=200,
    )
    with pytest.raises(ValueError, match="could not be matched"):
        make_nontargeting_controls(5, [0.05] * 10, arch, max_tries_per_oligo=50)
