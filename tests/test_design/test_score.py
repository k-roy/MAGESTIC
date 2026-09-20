"""The fitted score must reproduce `planning/648_fit_design_score.py` exactly, and must have the
shape properties the design relies on. If the fit is re-run, this test is what catches drift."""
import json
import math
import os

import pytest

from magestic.design.domain_insertion import score

_P = os.path.expanduser("~/work/UCSB/wilson_lab/planning")
FIT = os.environ.get("MAGESTIC_648_FIT", os.path.join(_P, "648_design_score_fit.json"))
PDIST_FIT = os.environ.get("MAGESTIC_670_FIT", os.path.join(_P, "670_pdist_plateau.json"))
EFF_FIT = os.environ.get("MAGESTIC_671_FIT", os.path.join(_P, "671_efficacy_calibration.json"))


def _json(path, what):
    if not os.path.exists(path):
        pytest.skip("%s not available at %s" % (what, path))
    return json.load(open(path))


def _fit():
    return _json(FIT, "648 fit JSON")


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
    assert score.design_score(None, 0.0, 0, 0, cut_distance=0) == float("-inf")
    assert score.design_score(None, 0.0, 0, 0, cut_distance=0, log=False) == 0.0
    assert score.design_score(0.9, 0.0, 0, 0, cut_distance=0) > score.design_score(None, 0.0, 0, 0, cut_distance=0)


def test_score_is_monotone_in_every_argument():
    base = dict(efficacy=0.9, sv=0.01, net_recoded=4, wm_in_gap=0, cut_distance=4)
    s0 = score.design_score(**base)
    assert score.design_score(**{**base, "efficacy": 0.96}) > s0   # 0.90->0.96 crosses a calibration bin
    assert score.design_score(**{**base, "cut_distance": 14}) < s0
    assert score.design_score(**{**base, "sv": 0.20}) < s0
    assert score.design_score(**{**base, "net_recoded": 9}) < s0
    assert score.design_score(**{**base, "wm_in_gap": 1}) < s0


def test_the_643_regression_cannot_recur():
    """YEL002C res283: a q=0.931 design with ONE mandatory strip codon must beat a q=0.476 design
    with five. Under the retired lexicographic rank_key the first was vetoed by `blocked`."""
    shipped = score.design_score(0.9332, 0.00168, 1, 0, cut_distance=3)
    fallthrough = score.design_score(0.4771, 0.00157, 5, 0, cut_distance=3)
    assert shipped > fallthrough


def test_log_and_linear_agree():
    a = score.design_score(0.8, 0.02, 5, 1, cut_distance=7, log=True)
    b = score.design_score(0.8, 0.02, 5, 1, cut_distance=7, log=False)
    assert math.exp(a) == pytest.approx(b, rel=1e-12)


# --- P_dist and the efficacy calibration, added 2026-08-28 -------------------------------------
# Same contract as `test_constants_match_the_fit`: the shipped tables are EMITTED by their fit
# scripts, never transcribed, and these assert the module reproduces them.


def test_pdist_constants_match_the_fit():
    want = _json(PDIST_FIT, "670 P_dist fit")["SCORE_PY_CONSTANTS"]["PDIST_BINS"]
    got = [list(b) for b in score.PDIST_BINS]
    assert len(got) == len(want)
    for g, w in zip(got, want):
        assert g == pytest.approx(w, abs=1e-6)


def test_efficacy_calibration_matches_the_fit():
    want = _json(EFF_FIT, "671 calibration")["SCORE_PY_CONSTANTS"]["EFFICACY_CAL_BINS"]
    got = [list(b) for b in score.EFFICACY_CAL_BINS]
    assert len(got) == len(want)
    for g, w in zip(got, want):
        assert g[:2] == pytest.approx(w[:2], abs=1e-9)
        assert g[2] == pytest.approx(w[2], abs=1e-6)


def test_p_dist_is_flat_through_the_plateau():
    """🔴 REGRESSION GUARD, and the reason `662`'s interpolated table was NOT shipped. Per-nt CE
    rates across d = 0..5 show no trend (Cochran-Armitage z = -1.55, p = 0.121, n = 1,777). A
    slope here decides 78 positions the data cannot decide, and -- because `policy`'s
    fewest-codons tie-break fires only on EXACT score equality -- silently disables it (`670`)."""
    for d in range(0, 6):
        assert score.p_dist(d) == pytest.approx(score.p_dist(0), abs=1e-12), d
    assert score.p_dist(6) < score.p_dist(5)


def test_p_dist_never_rewards_a_farther_cut():
    for d in range(0, 80):
        assert score.p_dist(d + 1) <= score.p_dist(d) + 1e-12, d
    assert score.p_dist(200) > 0.0          # clamped, never log(0)


def test_p_dist_uses_the_magnitude_of_a_signed_distance():
    """🔴 `designer.cut_distance` is SIGNED -- `policy.wm_in_gap` branches on `>= 0`. Dropping the
    abs() would score an upstream cut at -12 as if it were the plateau."""
    for d in (1, 4, 7, 12, 25):
        assert score.p_dist(-d) == score.p_dist(d), d


def test_p_eff_is_monotone_and_quantised():
    """Isotonic by construction: the raw bin rates invert at [0.60,0.80) vs [0.80,0.90), and an
    un-pooled curve would pay a design for LOWER predicted efficacy."""
    xs = [i / 200.0 for i in range(0, 201)]
    vs = [score.p_eff(x) for x in xs]
    for a, b in zip(vs, vs[1:]):
        assert b >= a - 1e-12
    assert len(set(vs)) == len(score.EFFICACY_CAL_BINS)   # a step function, on purpose


def test_p_eff_shrinks_efficacy_below_the_codon_term():
    """Kevin, 2026-08-28: efficacy must weigh less than the syn-codon count. Raw efficacy spans
    12.1x against P_tract's 1.59x; calibrated it spans 1.51x."""
    cal = [b[2] for b in score.EFFICACY_CAL_BINS]
    eff_range = max(cal) / min(cal)
    codon_range = score.p_tract(0) / score.p_tract(10)
    assert eff_range < codon_range
    assert eff_range == pytest.approx(1.51, abs=0.02)


def test_calibration_creates_the_ties_the_tiebreak_needs():
    """The mechanism behind the demotion: two guides whose RAW efficacies differ but land in one
    calibration bin must score EQUAL, so `policy`'s fewest-codons tie-break decides them."""
    a = score.design_score(0.82, 0.0, 3, 0, cut_distance=2)
    b = score.design_score(0.88, 0.0, 1, 0, cut_distance=2)
    assert a == pytest.approx(b, abs=1e-12)


def test_distance_now_outranks_efficacy():
    """🔴 THE `659` DEFECT. A top-efficacy guide cutting 12 nt away must LOSE to a mid-efficacy
    guide cutting inside the plateau. Before P_dist the scorer could not see this at all."""
    far_good = score.design_score(0.99, 0.0, 0, 0, cut_distance=12)
    near_ok = score.design_score(0.85, 0.0, 0, 0, cut_distance=2)
    assert near_ok > far_good


def test_cut_distance_is_required():
    """No default, on purpose: v2 shipped a silently distance-blind scorer (`659`), and a default
    would let a caller that forgets to pass it reproduce that bug as a no-op."""
    with pytest.raises(TypeError):
        score.design_score(0.9, 0.0, 0, 0)


def test_wilson_bounds_bracket_the_point_estimate():
    for d in (0, 5, 8, 12, 18, 30):
        assert score.p_dist(d, "lo_ci") <= score.p_dist(d) <= score.p_dist(d, "hi_ci"), d
