"""
The design score for the AsLOV2 insertion library -- ONE scalar, no gates.

    S = efficacy x P_tract(net_recoded) x P_fidelity(wm_in_gap) x (1 - sv)

Fitted in `planning/648_fit_design_score.py`; the rationale, the data and the limits are in
`planning/648_design_score_fit.md` and `planning/646_design_score_spec.md`.

WHY THIS REPLACES `policy.rank_key`
-----------------------------------
`rank_key` was lexicographic -- `(blocked, net_recoded >= CLIFF, -(q - lam*net_recoded))`. Two
booleans sat in front of the score, so a candidate could be vetoed regardless of how good it was:
at YEL002C res283 a q = 0.931 design sorted BELOW a q = 0.476 one because a single mandatory
cloning-site strip codon tripped `is_blocked` (`planning/643`). A veto cannot express
"slightly riskier", and it is what made a designer change able to make a position worse.

Here every consideration is a FACTOR in one scalar. Blocking, the recoding cliff and the W/M gap
are priced, not vetoed. Combined with variant enumeration in `policy.design_position`, that makes
"a new capability can never make a position worse" a property of the construction rather than
something to be checked by rebuilding the library and diffing it.

UNITS: every factor is a probability, and they multiply. That is the point. `lam` could never be
fitted because it priced a probability of CUTTING against a probability of COMPLETE INCORPORATION
-- different conditionals with no experiment varying both. Each factor here is estimated on its own
data and the exchange rate falls out.

🔴 ASSUMES INDEPENDENCE between cutting and incorporation. No dataset tests that. It is a better
assumption than a hand-chosen lambda; it is not a measurement.
🔴 THE CONSTANTS BELOW ARE EMITTED BY THE FIT, NOT TRANSCRIBED. `648_fit_design_score.py` writes
them to `648_design_score_fit.json["SCORE_PY_CONSTANTS"]`, and
`tests/test_design/test_score.py::test_constants_match_the_fit` asserts this module reproduces the
fitted `P_tract` table to 1e-6. Hand-copying them is how a shipped module silently drifts from the
analysis that justified it.
"""
import math

# --- P_tract: SEC14, n = 7,121 variants / 313 tracts -------------------------------------------
# Continuous hinge in the logit, knot chosen by AIC over k in 2..9 (k=3 wins; dAIC 66.1 vs a plain
# line). Window fixed effects; SEs clustered on the recoding TRACT, which is the unit of
# replication (`298` T2 was retracted for treating clones as independent).
#
#   below the knot  slope CONSTRAINED TO ZERO -- P_tract is FLAT through 3 recoded codons
#   above the knot  beta = +0.370, p_clustered 5e-4   -- the hinge carries the whole effect
#
# 🔴 TWO CORRECTIONS THAT ARE LOAD-BEARING (found 2026-08-27 by running the v2 designer):
#   (i)  SEC14's `no_recoding` class is EXCLUDED. It is 341 of 341 variants at n = 0 -- there is no
#        such thing in that dataset as a RECODED design with zero recoded codons. Leaving it in
#        makes the fit read a gap between two CONSTRUCT CLASSES as a dose response.
#   (ii) the below-knot slope is FIXED AT ZERO. Unconstrained, (i) produced P_tract(3) > P_tract(0)
#        -- "add codons to improve survival" -- and the v2 designer duly moved MRD1 res554 and
#        WBP1 res283 from 1 recoded codon to 3 for no reason. Within the recoded classes the
#        empirical rate is 0.050 / 0.071 / 0.051 / 0.068 at n = 1..4: flat, with no evidence that
#        3 beats 1. P_tract is therefore MONOTONE NON-INCREASING by construction and cannot
#        manufacture a reason to recode.
TRACT_KNOT = 3
TRACT_INTERCEPT = -2.993546674982344   # base logit at the mean window effect
TRACT_BETA_N = 0.0   # CONSTRAINED FLAT -- see the class-confound note below
TRACT_BETA_HINGE = 0.3699948042189722

# --- P_fidelity: NNS WGS, CE vs IE, at DESIGN level ---------------------------------------------
# `102`'s mechanism, refitted: an unrecodable Trp/Met inside the wild-type gap beside the cassette
# drives PARTIAL incorporation. Clustered logit beta_wm = 2.585, p = 0.0040; codon count on this
# axis is dead once the flag is in (p = 0.60), reproducing `102`.
# 🔴 The W/M arm is only 15 designs. 95% CI on P(complete) is 0.55-0.93. Direction solid,
# magnitude not pinned -- do not quote 0.800 as if it were tight.
FIDELITY_NO_WM = 0.96794            # 267 designs, IE 3.2 %
FIDELITY_WM = 0.80000               # 15 designs,  IE 20.0 % (95 % CI 0.55-0.93)

# A design with NO scored guide. Left explicit rather than defaulting to 0.0, which would make every
# unscored guide look catastrophic, or to 1.0, which would make it look perfect. `597` hit exactly
# this and fell back to a codon-only rule rather than silently scoring q = 0.
EFFICACY_WHEN_UNSCORED = None


def p_tract(net_recoded: int) -> float:
    """P(the design is not lost) as a function of recoded codons. Monotone above the knot."""
    n = max(int(net_recoded), 0)
    logit = TRACT_INTERCEPT + TRACT_BETA_N * n + TRACT_BETA_HINGE * max(n - TRACT_KNOT, 0)
    return 1.0 - 1.0 / (1.0 + math.exp(-logit))


def p_fidelity(wm_in_gap: int) -> float:
    """P(complete incorporation | edited). `wm_in_gap` counts UNRECODABLE codons in the wild-type
    run between the cassette and the nearest edit -- the `102` mechanism, not merely a gap length."""
    return FIDELITY_WM if int(wm_in_gap or 0) > 0 else FIDELITY_NO_WM


def design_score(efficacy, sv, net_recoded, wm_in_gap, *, log=True):
    """The scalar. Returns -inf (or 0.0 with log=False) when efficacy is unavailable, so an
    unscored guide loses to any scored one without poisoning the arithmetic."""
    if efficacy is None:
        return float("-inf") if log else 0.0
    eff = min(max(float(efficacy), 1e-9), 1.0)
    sv_ok = min(max(1.0 - float(sv or 0.0), 1e-9), 1.0)
    pt = p_tract(net_recoded)
    pf = p_fidelity(wm_in_gap)
    if not log:
        return eff * pt * pf * sv_ok
    return math.log(eff) + math.log(pt) + math.log(pf) + math.log(sv_ok)


def implied_lambda(n: int, q_scale: float = 0.85) -> float:
    """The per-codon cost of going n -> n+1, in the additive efficacy units the retired
    `policy.LAMBDA_DEFAULT = 0.02` was expressed in. Provided so the new score stays auditable
    against the old one: it is 0.024 at n=4 (lam was about right there), NEGATIVE below 3
    (the first three codons are FREE), and 0.105 at n=9 -- five times what lam charged."""
    return q_scale * (1.0 - p_tract(n + 1) / p_tract(n))
