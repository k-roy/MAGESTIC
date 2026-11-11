"""G-COMPOSITION: the built pool must be the pool the recipe declared (NN 505, NN 506 §4).

**The bug this pins.** The TRIUMPH pool shipped with three composition regressions and not
one existing gate saw them, because every gate asked whether the pool was well-formed and
none asked whether it was the pool that was specified:

1. the ``syn_recoded`` control arm disappeared entirely — 252 oligos, a whole category;
2. an arm can also be partially lost, which is the same defect one size down;
3. 🔴 composition drifted, ``natural`` absorbing 454 slots from every other category, **and
   the total stayed at exactly 24,000** — so a matching total was read as evidence that the
   composition was right.

**The contract.** Every arm's built count equals its target, and the sum equals ``total``.
The remainder arm's target is derived from ``total`` minus the *declared* targets, so it
cannot absorb a loss: a shrunken arm still fails its own target. 🔴 **The gate never pads to
fit** — a mismatch is a failure, never something the gate corrects.

Each test below plants the exact defect the gate is meant to catch, because a gate that has
never been shown to fire is not evidence.
"""

from __future__ import annotations

import copy

import pytest

from magestic.design.pool import gate_composition
from magestic.design.pool.gates import FAIL, PASS

from . import _pool_fixtures as fx


def test_a_correct_pool_passes():
    result = gate_composition(fx.built_pool(), fx.recipe())
    assert result.status == PASS, result.detail
    # every arm accounted for, including the derived remainder
    by_arm = {row[0]: row for row in result.rows}
    assert by_arm["natural"][2] == 45, "remainder target should be 100 - 55"
    assert all(row[4] == 0 for row in result.rows), result.rows


def test_a_whole_category_vanishing_is_caught():
    """Defect 1, planted: the 252-oligo ``syn_recoded`` arm disappears."""
    pool = fx.built_pool(syn_recoded=0, natural=50)  # total still 100
    result = gate_composition(pool, fx.recipe())

    assert result.status == FAIL
    assert "syn_recoded" in result.detail
    assert "ABSENT" in result.detail
    assert "target 5" in result.detail


def test_a_partially_lost_category_is_caught():
    """Defect 2, planted: an arm survives but short."""
    pool = fx.built_pool(missense_recoded=14, natural=51)  # total still 100
    result = gate_composition(pool, fx.recipe())

    assert result.status == FAIL
    assert "missense_recoded" in result.detail
    assert "target 20, built 14" in result.detail
    assert "-6" in result.detail


def test_drift_is_caught_even_though_the_total_is_exactly_right():
    """🔴 Defect 3, planted — the one that escaped every existing gate.

    ``natural`` absorbs what the other arms lost, so ``len(pool)`` is exactly ``total``.
    A gate that checked only the total would pass this pool. This one must not.
    """
    pool = fx.built_pool(ptc=6, syn=7, missense_recoded=18, natural=54)
    assert len(pool) == 100, "the planted defect must keep the total exactly on target"

    result = gate_composition(pool, fx.recipe())

    assert result.status == FAIL, (
        "a pool whose total matches but whose composition drifted must FAIL — this is "
        "precisely the NN 505 regression"
    )
    for arm in ("ptc", "syn", "missense_recoded", "natural"):
        assert arm in result.detail, f"{arm} drifted but is not named in the failure"


def test_the_remainder_arm_cannot_absorb_a_loss():
    """The derived target moves with the declared ones, so absorbing is never silent.

    Here only ``ptc`` is short and ``natural`` takes the slack. ``natural``'s target is
    ``total`` minus the DECLARED targets, so it stays 45 and both arms are reported wrong.
    """
    pool = fx.built_pool(ptc=0, natural=55)
    result = gate_composition(pool, fx.recipe())

    assert result.status == FAIL
    by_arm = {row[0]: row for row in result.rows}
    assert by_arm["natural"][2] == 45, "the remainder target must not follow the built count"
    assert by_arm["natural"][3] == 55
    assert by_arm["ptc"][4] == -10


def test_an_undeclared_arm_in_the_pool_is_caught():
    """A category added to the pool but never written down is also a failure."""
    pool = fx.built_pool()
    for row in fx.make_rows("mystery", 5):
        row["subpool"] = "mystery"
        pool.append(row)

    result = gate_composition(pool, fx.recipe())

    assert result.status == FAIL
    assert "mystery" in result.detail
    assert "not declared" in result.detail


def test_a_wrong_total_fails_and_says_it_never_pads():
    pool = fx.built_pool(natural=40)  # 95 oligos
    result = gate_composition(pool, fx.recipe())

    assert result.status == FAIL
    assert "declared 100, built 95" in result.detail
    assert "never pads to fit" in result.detail


def test_the_gate_does_not_modify_the_pool():
    """🔴 It never pads to fit: the pool it was handed comes back untouched."""
    pool = fx.built_pool(ptc=3)
    before = copy.deepcopy(pool)

    result = gate_composition(pool, fx.recipe())

    assert result.status == FAIL
    assert pool == before, "G-COMPOSITION must report a mismatch, never correct one"


def test_every_arm_appears_in_the_report_rows_even_when_passing():
    """The per-arm table is the artifact a reader checks, so it is always populated."""
    result = gate_composition(fx.built_pool(), fx.recipe())
    assert result.header == ("arm", "kind", "target", "built", "delta")
    assert {row[0] for row in result.rows} == {
        "natural", "ptc", "syn", "missense_recoded", "syn_recoded", "scrambled",
    }


@pytest.mark.parametrize(
    "overrides",
    [
        {"syn_recoded": 0, "natural": 50},
        {"ptc": 6, "syn": 7, "missense_recoded": 18, "natural": 54},
        {"natural": 40},
    ],
)
def test_no_planted_defect_ever_passes(overrides):
    """A sweep: none of the planted shapes may come back PASS."""
    assert gate_composition(fx.built_pool(**overrides), fx.recipe()).status == FAIL
