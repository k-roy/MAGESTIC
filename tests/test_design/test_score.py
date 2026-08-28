"""The fitted score must reproduce `planning/648_fit_design_score.py` exactly, and must have the
shape properties the design relies on. If the fit is re-run, this test is what catches drift."""
import json
import math
import os

import pytest

from magestic.design.domain_insertion import score

FIT = os.environ.get(
    "MAGESTIC_648_FIT",
    os.path.expanduser("~/work/UCSB/wilson_lab/planning/648_design_score_fit.json"))


def _fit():
    if not os.path.exists(FIT):
        pytest.skip("648 fit JSON not available at %s" % FIT)
    return json.load(open(FIT))


def test_constants_match_the_fit():
    c = _fit()["SCORE_PY_CONSTANTS"]
    assert score.TRACT_KNOT == c["TRACT_KNOT"]
    for name, key in (("TRACT_INTERCEPT", "TRACT_INTERCEPT"), ("TRACT_BETA_N", "TRACT_BETA_N"),
                      ("TRACT_BETA_HINGE", "TRACT_BETA_HINGE"),
                      ("FIDELITY_NO_WM", "FIDELITY_NO_WM"), ("FIDELITY_WM", "FIDELITY_WM")):
        assert getattr(score, name) == pytest.approx(c[key], abs=1e-12), name


def test_p_tract_matches_the_fitted_table():
    tbl = _fit()["P_tract"]
    for n, want in tbl.items():
        assert score.p_tract(int(n)) == pytest.approx(want, abs=1e-5), n


def test_p_tract_is_flat_then_falls():
    """The whole point of the hinge: FLAT through the knot, then falling."""
    for n in range(0, score.TRACT_KNOT + 1):
        assert score.p_tract(n) == pytest.approx(score.p_tract(0), abs=1e-12), n
    for n in range(score.TRACT_KNOT, 20):
        assert score.p_tract(n + 1) < score.p_tract(n), n


def test_p_tract_never_rewards_recoding():
    """🔴 REGRESSION GUARD. An unconstrained fit on SEC14 gives P_tract(3) > P_tract(0), because
    n=0 there is entirely the `no_recoding` CONSTRUCT CLASS. That made the v2 designer add codons
    at MRD1 res554 and WBP1 res283 for no reason. P_tract must be monotone NON-INCREASING."""
    for n in range(0, 25):
        assert score.p_tract(n + 1) <= score.p_tract(n) + 1e-12, n
    for n in range(0, 25):
        assert score.implied_lambda(n) >= -1e-12, n


def test_implied_lambda_brackets_the_retired_constant():
    """`policy.LAMBDA_DEFAULT` was 0.02 flat. It is right at n=4 and wrong on both sides."""
    assert score.implied_lambda(0) == pytest.approx(0.0, abs=1e-9)   # the first codons are FREE
    assert score.implied_lambda(4) == pytest.approx(0.02, abs=0.01)
    assert score.implied_lambda(10) > 4 * 0.02    # badly undercharged in the tail


def test_wm_gap_costs_more_than_no_gap():
    assert score.p_fidelity(1) < score.p_fidelity(0)
    assert score.p_fidelity(0) == score.FIDELITY_NO_WM
    assert score.p_fidelity(3) == score.FIDELITY_WM       # any unrecodable codon, not a count


def test_unscored_guide_loses_without_poisoning_the_arithmetic():
    assert score.design_score(None, 0.0, 0, 0) == float("-inf")
    assert score.design_score(None, 0.0, 0, 0, log=False) == 0.0
    assert score.design_score(0.9, 0.0, 0, 0) > score.design_score(None, 0.0, 0, 0)


def test_score_is_monotone_in_every_argument():
    base = dict(efficacy=0.9, sv=0.01, net_recoded=4, wm_in_gap=0)
    s0 = score.design_score(**base)
    assert score.design_score(**{**base, "efficacy": 0.95}) > s0
    assert score.design_score(**{**base, "sv": 0.20}) < s0
    assert score.design_score(**{**base, "net_recoded": 9}) < s0
    assert score.design_score(**{**base, "wm_in_gap": 1}) < s0


def test_the_643_regression_cannot_recur():
    """YEL002C res283: a q=0.931 design with ONE mandatory strip codon must beat a q=0.476 design
    with five. Under the retired lexicographic rank_key the first was vetoed by `blocked`."""
    shipped = score.design_score(0.9332, 0.00168, 1, 0)
    fallthrough = score.design_score(0.4771, 0.00157, 5, 0)
    assert shipped > fallthrough


def test_log_and_linear_agree():
    a = score.design_score(0.8, 0.02, 5, 1, log=True)
    b = score.design_score(0.8, 0.02, 5, 1, log=False)
    assert math.exp(a) == pytest.approx(b, rel=1e-12)
