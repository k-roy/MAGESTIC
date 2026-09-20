"""
The design score for the AsLOV2 insertion library -- ONE scalar, no gates.

    S = P_eff(efficacy) x P_dist(cut_distance) x P_tract(net_recoded)
                        x P_fidelity(wm_in_gap) x (1 - sv)

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

THE 2026-08-28 REVISION -- two factors, both MEASUREMENTS, no new free parameter
--------------------------------------------------------------------------------
Kevin, 2026-08-28: *"informed by our previous data and first principles, but not too complicated
... we should weight the predicted efficacy lower than the number of syn mutations required."*

(1) `P_dist` ADDED. v2 kept cut-to-insertion distance as a hard gate (`policy.MAX_DISTANCE = 30`)
    but dropped it as a RANKING term (`planning/659`), leaving the score blind to the single
    largest measured effect in the data: 74.2 % correct edits at 0-5 nt falling to 30.0 % at
    10-14 nt (`662` on the cut axis, `670` pooled). It is FLAT through 5 nt -- see `p_dist`.

(2) `efficacy` CALIBRATED. It was multiplied in as though it were a probability. It is a MIN-MAX
    NORMALISED MODEL OUTPUT, so its numeric range (0.082-0.993, a 12.1x spread) is set by the
    extremes of SCORE's training set and not by any outcome rate. In our v1.3 / SNR52 regime the
    CE rate it actually separates spans 0.606-0.912, a 1.51x spread (`671`). `P_eff` replaces the
    raw value with the rate measured at that value, shrinking its authority 3.0x.

    ⭐ THIS IS THE DEMOTION KEVIN ASKED FOR, AND IT IS NOT A WEIGHT. Two consequences:
      * across full domains, calibrated efficacy (1.51x) now sits just BELOW the codon term
        (P_tract 0->10 codons, 1.59x), where raw efficacy sat 7.6x above it;
      * more importantly, a calibrated value takes only as many distinct levels as there were
        bins to measure it with -- we cannot rank on a difference finer than we measured. Two
        designs in the same bin TIE, and `policy.design_position`'s fewest-codons tie-break
        decides them. Below the knot (`P_tract` exactly flat, 1-3 recoded codons) that moves
        517 of the `659` hit positions from "raw efficacy decides" to "fewest codons decides"
        for 326 of them -- 63.1 %, against 1.5 % before (`671`).

🔴 ASSUMES INDEPENDENCE between cutting and incorporation. No dataset tests that. It is a better
assumption than a hand-chosen lambda; it is not a measurement.
🔴 P_dist RESTS ON AN EXTRAPOLATION. SCORE Data S1 measures SNV/MNV/INDEL variants; our designed
change set is a 408-nt cassette. The farthest-designed-change mapping (`662` section 1) is the
best available, not a validated one, and the bins past 6 nt rest on ~100 observations.
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

# --- P_dist: SCORE Data S1, n = 1,875 scored clones, indexed on the CUT axis --------------------
# `cut = PAM_location - 3` on '+' and `+ 3` on '-', verified against the data rather than assumed
# (`662` section 1). Bins are (midpoint, p, lo95, hi95); linear interpolation between midpoints,
# FLAT below the first, held at the last beyond it.
#
# 🔴 THE FIRST BIN'S MIDPOINT IS 5.0, NOT 1.0, AND THAT IS THE WHOLE POINT. `662` shipped separate
#    0-2 and 3-5 bins and interpolated between their midpoints, which puts a ~1.2 % SLOPE across
#    d = 1..4. There is no such slope in the data: per-nucleotide CE rates across d = 0..5 give a
#    Cochran-Armitage trend of z = -1.55, p = 0.121 at n = 1,777 -- flat by measurement (`670`).
#    Shipping the slope would have been the same artifact `662` section 4 retracted on the PAM axis,
#    and it would have SILENTLY DISABLED the fewest-codons tie-break, which fires only on exact
#    score equality: it decided 78 plateau positions that the data cannot decide. The two bins are
#    therefore POOLED (74.5 % / 73.6 % -> 74.2 %, CIs 72.0-76.8 and 69.5-77.3 overlapping), and the
#    curve is flat to 5 nt. Repairs the IDENTICAL 434/455 material distance losses either way.
PDIST_BINS = (
    (5.0, 0.7422622397298818, 0.7214181812093063, 0.762061090613634),
    (7.5, 0.5517241379310345, 0.4245185822161052, 0.6725034870705885),
    (12.0, 0.3, 0.14547527396899385, 0.5189767762289793),
    (17.0, 0.07692307692307693, 0.013710066379696323, 0.3331453392824736),
    (40.0, 0.0, 0.0, 0.35433884297520657),
)

# --- P_eff: the efficacy calibration, v1.3 / SNR52 only, n = 1,323 ------------------------------
# (lo, hi, P(CE)) over the RAW min-max efficacy. Bin edges chosen for SUFFICIENCY -- the fewest
# bins keeping every bin above n = 100, so each level carries a usable Wilson interval.
#
# 🔴 ISOTONIC (pool-adjacent-violators). The raw rates are NOT monotone -- [0.60,0.80) reads 84.7 %
#    against [0.80,0.90)'s 82.5 %, on heavily overlapping intervals -- and shipped as measured the
#    curve would PAY A DESIGN FOR HAVING LOWER PREDICTED EFFICACY. That is the same defect
#    `TRACT_BETA_N = 0.0` above records: a fit must not be able to manufacture a reason to prefer
#    the worse input. PAVA merges exactly the inverted pair; the 1.51x range is unchanged.
#
# 🔴 v1.3 ONLY, deliberately. Data S1 pools three MAGESTIC systems whose CE rates are 33.3 / 50.3 /
#    82.7 % (`664` section 2). `pS1380` is SNR52, i.e. the v1.3-like regime, and it is also the one
#    where the score is LEAST predictive (AUC 0.626 for CE vs 0.723 in v1.2) because much of what
#    SCORE.NE encodes is Pol III termination at T-runs -- a failure mode of the tRNA promoter we do
#    not use (`664` section 5).
EFFICACY_CAL_BINS = (
    (0.00, 0.60, 0.6058394160583942),
    (0.60, 0.90, 0.8304195804195804),
    (0.90, 0.95, 0.8369905956112853),
    (0.95, 1.00, 0.9118644067796611),
)


def p_tract(net_recoded: int) -> float:
    """P(the design is not lost) as a function of recoded codons. Monotone above the knot."""
    n = max(int(net_recoded), 0)
    logit = TRACT_INTERCEPT + TRACT_BETA_N * n + TRACT_BETA_HINGE * max(n - TRACT_KNOT, 0)
    return 1.0 - 1.0 / (1.0 + math.exp(-logit))


def p_dist(cut_distance, which: str = "p") -> float:
    """P(correct edit) as a function of cut-to-farthest-designed-change distance.

    🔴 `cut_distance` IS SIGNED in `designer.py` (`policy.wm_in_gap` branches on `>= 0` to find
    which side the cassette sits on). Only the magnitude is meaningful here, and dropping the
    `abs()` would silently score an upstream cut at -4 as if it were 0. `which` selects the point
    estimate or a Wilson bound, so a close call can be re-run at both (`662` section 3).
    """
    col = {"p": 1, "lo_ci": 2, "hi_ci": 3}[which]
    d = abs(float(cut_distance))
    xs = [b[0] for b in PDIST_BINS]
    ys = [b[col] for b in PDIST_BINS]
    if d <= xs[0]:
        v = ys[0]                       # the measured plateau -- flat, not interpolated
    elif d >= xs[-1]:
        v = ys[-1]
    else:
        i = max(j for j in range(len(xs) - 1) if d > xs[j])
        f = (d - xs[i]) / (xs[i + 1] - xs[i])
        v = ys[i] + f * (ys[i + 1] - ys[i])
    return min(max(v, 1e-9), 1.0)


def p_eff(efficacy) -> float:
    """The raw min-max efficacy mapped to the CE rate measured at that value in v1.3 / SNR52.

    Returns a STEP function on purpose: we cannot rank on a difference finer than we measured, so
    two guides inside one bin tie and `policy`'s fewest-codons tie-break decides them.
    """
    e = min(max(float(efficacy), 0.0), 1.0)
    for lo, hi, p in EFFICACY_CAL_BINS:
        if lo <= e < hi:
            return p
    return EFFICACY_CAL_BINS[-1][2]      # e == 1.0 exactly: the top bin is closed on the right


def p_fidelity(wm_in_gap: int) -> float:
    """P(complete incorporation | edited). `wm_in_gap` counts UNRECODABLE codons in the wild-type
    run between the cassette and the nearest edit -- the `102` mechanism, not merely a gap length."""
    return FIDELITY_WM if int(wm_in_gap or 0) > 0 else FIDELITY_NO_WM


def design_score(efficacy, sv, net_recoded, wm_in_gap, *, cut_distance, log=True,
                 which: str = "p"):
    """The scalar. Returns -inf (or 0.0 with log=False) when efficacy is unavailable, so an
    unscored guide loses to any scored one without poisoning the arithmetic.

    `cut_distance` is KEYWORD-ONLY AND REQUIRED. It has no default on purpose: v2 shipped a scorer
    that was silently blind to distance (`659`), and a default would let a caller that forgets to
    pass it reproduce that bug as a no-op instead of a TypeError.
    """
    if efficacy is None:
        return float("-inf") if log else 0.0
    eff = p_eff(efficacy)
    sv_ok = min(max(1.0 - float(sv or 0.0), 1e-9), 1.0)
    pd = p_dist(cut_distance, which)
    pt = p_tract(net_recoded)
    pf = p_fidelity(wm_in_gap)
    if not log:
        return eff * pd * pt * pf * sv_ok
    return math.log(eff) + math.log(pd) + math.log(pt) + math.log(pf) + math.log(sv_ok)


def implied_lambda(n: int, q_scale: float = 0.85) -> float:
    """The per-codon cost of going n -> n+1, in the additive efficacy units the retired
    `policy.LAMBDA_DEFAULT = 0.02` was expressed in. Provided so the new score stays auditable
    against the old one: it is 0.024 at n=4 (lam was about right there), NEGATIVE below 3
    (the first three codons are FREE), and 0.105 at n=9 -- five times what lam charged."""
    return q_scale * (1.0 - p_tract(n + 1) / p_tract(n))
