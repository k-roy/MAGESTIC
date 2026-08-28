"""
Guide-selection policy for the essentialome AsLOV2 library.

Installed as `magestic/design/domain_insertion/policy.py` (2026-07-30, Kevin: *"Apply arm F at
lambda 0.02 and regenerate the figures"*). The master copy lives in this workspace at
`planning/120_policy.py`; `120_apply_arm_f.py` installs it and patches the designer.

`design_oligo` returns the FIRST guide that yields any design, ordered by distance bucket then
predicted quality. `118` measured what that costs, and this module replaces it with the policy that
won:

  gate 1  Trp/Met-BLOCKED  (WT-identical run >= 3 nt between the cassette and the nearest SNV;
                            `102`/`103`: 22-25 % partial incorporation vs 2-4 % at <= 2)
  gate 2  >= 7 recoded codons  (`101` NE cliff, `108b` editing efficiency, `103` WGS replication)
  then    maximise  quality - lam * net_recoded,  quality = efficacy - SV   (`guide_quality`)

Both gates are PURITY failures -- an unedited or partially-edited locus reads as "no phenotype" and
no sequencing depth recovers it (HANDOFF, "abundance and PURITY are NOT one currency"). Everything
below them is an abundance trade, which is what `lam` prices.

`lam` = 0.02 efficacy-units per recoded codon, from Kevin's own bound (2026-07-30): *"If we need 10
codons of recoding to get a 0.95 guide, but another guide with 0.84 sits close by, seems better to
go with the 0.84"* => 10 codons must outweigh 0.11 efficacy (lam > 0.011), and one codon must not
flip a 0.11 gap (lam < 0.11). `118`'s sweep confirmed the lower bound empirically: at lam = 0.01 the
>= 7-codon fraction rose to 8.44 %, WORSE than the unfixed designer's 7.21 %.

Measured on 5,000 positions (`118`, arm F lam=0.02) against the stock designer:

    designable  95.40 -> 97.18 %      blocked   2.516 -> 0.412 %
    >= 7 codons  7.212 -> 4.384 %     codons     3.829 -> 2.672
    efficacy    0.8603 -> 0.8873

Better on every axis at once. It also strictly dominates the pure-scalar arm (E) and the banded arm
(D): same gates as D, yet fewer codons AND higher efficacy.
"""
from .designer import design_oligo
from .guides import guide_quality, guides_near

CLIFF = 7                    # `101`/`103`/`108b`: NE and editing efficiency break at >= 7 codons
BLOCK_GAP = 3                # `102`/`103`: WT gap >= 3 nt beside the cassette is the blocked state
LAMBDA_DEFAULT = 0.02        # efficacy-units per recoded codon (`118` arm F)
MIN_DISRUPTION_STEPS = (4, 3, 2)
MAX_DISTANCE = 30


def gap_beside_cassette(design, window, ins):
    """WT-identical nt between the cassette and the nearest recoded base.

    This is the `82`/`102` connectivity metric: a conversion tract can copy an SNV WITHOUT reaching
    the insertion, destroying the protospacer while leaving the cassette uninstalled -- a clone that
    looks quietly wild type. Returns a 99 sentinel when there is no SNV at all.
    """
    n_left, n_right = len(design.left_arm), len(design.right_arm)
    nat_l, nat_r = window[ins - n_left:ins], window[ins:ins + n_right]
    g_r = next((j for j, (a, b) in enumerate(zip(design.right_arm, nat_r)) if a != b), None)
    g_l = next((j for j, (a, b) in enumerate(zip(design.left_arm[::-1], nat_l[::-1])) if a != b), None)
    found = [g for g in (g_r, g_l) if g is not None]
    return min(found) if found else 99


def is_blocked(design, gap):
    """A design with ZERO recoded bases is NOT blocked -- it is the opposite.

    The hazard is a tract copying an SNV without the cassette; with no SNV to copy, the only way to
    destroy the site is to install the cassette, so coupling is perfect. `gap_beside_cassette`
    returns its 99 sentinel in that case and must not be read as a large gap.
    """
    return design.net_recoded > 0 and gap >= BLOCK_GAP


# --- `618` two-hurdle block rule (Kevin, 2026-08-21). DEFAULT OFF. -------------------------------
TWO_HURDLE_BLOCK = False     # flip to True to use hazard_blocked() instead of is_blocked()


def hazardous_edits(design, window, ins, guide_length):
    """Window-frame indices of the changed bases that clear BOTH hurdles.

    hurdle 1 -- NOT shielded: |edit - cut| <= |insertion - cut|. Incorporation is first-come,
                first-served outward from the cut, so an edit further from the cut than the
                insertion cannot be copied without the cassette, whatever the wild-type gap.
    hurdle 2 -- inside protospacer+PAM, so it can block re-cleavage. Outside that window the site
                stays cuttable and is simply re-cut until the insertion lands.
    """
    from .designer import protospacer_span
    n_left, n_right = len(design.left_arm), len(design.right_arm)
    nat_l = window[ins - n_left:ins]
    nat_r = window[ins:ins + n_right]
    changed = [ins - 1 - j for j, (a, b) in enumerate(zip(design.left_arm[::-1], nat_l[::-1])) if a != b]
    changed += [ins + j for j, (a, b) in enumerate(zip(design.right_arm, nat_r)) if a != b]
    if not changed:
        return []
    cut = ins + design.cut_distance                      # cut junction, window frame
    d_ins = abs(design.cut_distance)                     # insertion's distance from the cut
    p_lo, p_hi = protospacer_span(design.guide, guide_length)
    return [i for i in changed
            if abs(i - cut) <= d_ins                     # hurdle 1 failed: not shielded
            and p_lo <= i < p_hi]                        # hurdle 2 failed: blocks re-cleavage


def hazard_gap(design, window, ins, guide_length):
    """Distance from the cassette to the nearest HAZARDOUS edit; 99 sentinel when there is none.

    Same 99 convention as `gap_beside_cassette`, and the same trap: 99 means "no hazard", NOT
    "a very large gap". `hazard_blocked` reads it only together with the emptiness test.
    """
    haz = hazardous_edits(design, window, ins, guide_length)
    if not haz:
        return 99, 0
    n_left = len(design.left_arm)
    gaps = [(i - ins) if i >= ins else (ins - 1 - i) for i in haz]
    return min(gaps), len(haz)


def hazard_blocked(design, window, ins, guide_length):
    """Two-hurdle replacement for `is_blocked`.

    Blocked only when a hazardous edit exists AND it is disconnected from the cassette
    (>= BLOCK_GAP). An in-protospacer edit that `_connect_edits_to_insertion` chained back to the
    insertion has gap <= MAX_GAP_DEFAULT and stays unblocked -- so this keeps rejecting exactly the
    configuration the guard exists for, including the post-final-strip case.
    """
    gap, n = hazard_gap(design, window, ins, guide_length)
    return n > 0 and gap >= BLOCK_GAP


# --- `621` (B) guide-local fallback = Kevin's rule (c). DEFAULT OFF. -----------------------------
# "When an improved design is rejected for any reason, use that guide's PREVIOUS design rather than
#  a different one."  With (A) and (C) in place this should rarely fire; its value is as a general
#  backstop, so no change of this kind can silently make a position worse than what already ships.
#
# 🔴 Kevin asked for the FIRE COUNT (2026-08-22).  A module-level counter is useless in a 400-task
# array, so the runner MUST dump `fallback_report()` per chunk -- see `622`.
GUIDE_LOCAL_FALLBACK = False

import collections as _collections
FALLBACK_COUNTS = _collections.Counter()


def fallback_reset():
    FALLBACK_COUNTS.clear()


def fallback_report():
    """Plain dict, safe to json.dump per chunk."""
    return dict(FALLBACK_COUNTS)


def _glf_recover(rec, model, chrom_seq, residue, spec, usage, scores, slide_range,
                 max_distance, window, ins):
    """Replace `rec` in place with the guide's UNCREDITED design when that ranks better.

    Reasons counted separately, because they mean different things:
      no_design      the credited design does not exist at all (the `604` junction class)
      blocked        the credited design exists but `is_blocked`/`hazard_blocked` fires (`599`)
      outranked_check the credited design is fine; compared anyway (the `627` unconditional arm --
                     this is the one that catches HCA4-style rank reordering)
      taken          the uncredited design was actually adopted
      kept_credited  the uncredited design existed but did not rank better -- no change
      unavailable    neither state yields a design; nothing to fall back to
    """
    d_c = rec["design"]
    reason = ("no_design" if d_c is None
              else "blocked" if rec.get("blocked") else "outranked_check")
    FALLBACK_COUNTS["considered"] += 1
    FALLBACK_COUNTS["considered_" + reason] += 1

    d_u, md_u = stepped_design(model, chrom_seq, residue, spec, usage=usage,
                               allowed_guides={rec["protospacer"]}, scores=scores,
                               slide_range=slide_range, max_distance=max_distance,
                               credit_pam_separated=False)
    if d_u is None:
        FALLBACK_COUNTS["unavailable"] += 1
        return
    alt = {"guide": rec["guide"], "protospacer": rec["protospacer"], "design": d_u,
           "min_disruption_used": md_u, "quality": rec["quality"]}
    alt["gap"] = gap_beside_cassette(d_u, window, ins)
    if TWO_HURDLE_BLOCK:
        alt["hazard_gap"], alt["n_hazard"] = hazard_gap(d_u, window, ins, spec.guide_length)
        alt["blocked"] = hazard_blocked(d_u, window, ins, spec.guide_length)
    else:
        alt["blocked"] = is_blocked(d_u, alt["gap"])
    alt["net_recoded"] = d_u.net_recoded
    if d_c is None or rank_key(alt) < rank_key(rec):
        alt["credit_fallback"] = True
        rec.clear()
        rec.update(alt)
        FALLBACK_COUNTS["taken"] += 1
        FALLBACK_COUNTS["taken_" + reason] += 1
    else:
        FALLBACK_COUNTS["kept_credited"] += 1


def stepped_design(model, chrom_seq, residue, spec, *, usage, allowed_guides, scores,
                   slide_range, max_distance=MAX_DISTANCE, credit_pam_separated=False):
    """`112` clause 3: design at min_disruption 4; if the result needs >= CLIFF codons, retry at 3
    then 2 and take the first that clears the cliff. NOT a blanket relaxation -- `99`'s exchange
    rate (1 codon ~ 1.4 disruption points) makes that a net loss; it only pays where it crosses the
    cliff, which the linear rate does not price. Returns (design_or_None, min_disruption_used)."""
    best, used = None, None
    for md in MIN_DISRUPTION_STEPS:
        d = design_oligo(model, chrom_seq, residue, spec, usage=usage, max_distance=max_distance,
                         min_disruption=md, allowed_guides=allowed_guides, scores=scores,
                         slide_range=slide_range, credit_pam_separated=credit_pam_separated)
        if d is not None and best is None:
            best, used = d, md
        if d is not None and d.net_recoded < CLIFF:
            return d, md
    return best, used


def evaluate_candidates(model, chrom_seq, residue, spec, *, usage, allowed_guides, scores,
                        slide_range, max_distance=MAX_DISTANCE, flank=200,
                        credit_pam_separated=False):
    """Force each candidate guide in turn and return what each would actually produce.

    Yields one record per candidate within `max_distance`, in `guides_near` order (proximity, then
    strand, then PAM index) -- the order matters, because ranking ties are broken by it.
    """
    window, c_idx, _, _ = model.window_with_coords(
        chrom_seq, model.cds_to_genomic(residue * 3), flank)
    ins = c_idx + 1
    out = []
    for g in guides_near(window, ins, max_distance, guide_length=spec.guide_length):
        p = g.protospacer
        if allowed_guides is not None and p not in allowed_guides:
            continue                                  # off-target screen: the precomputed set
        d, md = stepped_design(model, chrom_seq, residue, spec, usage=usage,
                               allowed_guides={p}, scores=scores, slide_range=slide_range,
                               max_distance=max_distance,
                               credit_pam_separated=credit_pam_separated)
        rec = {"guide": g, "protospacer": p, "design": d, "min_disruption_used": md,
               "quality": guide_quality(p, scores)}
        if d is not None:
            rec["gap"] = gap_beside_cassette(d, window, ins)
            if TWO_HURDLE_BLOCK:
                rec["hazard_gap"], rec["n_hazard"] = hazard_gap(d, window, ins, spec.guide_length)
                rec["blocked"] = hazard_blocked(d, window, ins, spec.guide_length)
            else:
                rec["blocked"] = is_blocked(d, rec["gap"])
            rec["net_recoded"] = d.net_recoded
        # `627` (Kevin, 2026-08-22): UNCONDITIONAL both-ways comparison. The conditional form
        # (only when the credited design is None or blocked) covered every regression that had
        # actually OCCURRED, but `612` §6 was explicit that this shows the case did not occur, not
        # that it cannot -- and `623` then produced one: YJL033W_HCA4_res537, where the credited
        # design was neither None nor blocked, merely OUTRANKED after crediting reordered the
        # candidates (efficacy 0.624 -> 0.255). Comparing both states for EVERY candidate makes
        # `min(ON) <= min(OFF)` in `rank_key` hold by construction, which is the actual proof.
        if GUIDE_LOCAL_FALLBACK and credit_pam_separated:
            _glf_recover(rec, model, chrom_seq, residue, spec, usage, scores, slide_range,
                         max_distance, window, ins)
        out.append(rec)
    return out, window, ins


def rank_key(rec, lam=LAMBDA_DEFAULT):
    """Arm F. Lower sorts better. Only call on records whose design is not None."""
    q = rec["quality"]
    q = 0.0 if q is None else q
    return (bool(rec["blocked"]),
            1 if rec["net_recoded"] >= CLIFF else 0,
            -(q - lam * rec["net_recoded"]))


def design_position(model, chrom_seq, residue, spec, *, usage, allowed_guides, scores,
                    slide_range, max_distance=MAX_DISTANCE, lam=LAMBDA_DEFAULT, flank=200,
                    return_candidates=False, credit_pam_separated=False):
    """Design one position under the `118` arm-F policy.

    Returns (design, min_disruption_used) -- or (design, min_disruption_used, candidates, window,
    ins) when `return_candidates` is set, which is what the `119` figures use to label every
    candidate PAM with the design it would have given.
    """
    # 🔴 NOTE FOR POOLS 1/2: `rank_key` scores candidates partly by `net_recoded`, so crediting a
    # PAM-separated guide does not merely shrink recoding -- it can REORDER the candidates and pick a
    # different guide. That is a larger blast radius than the heterologous pools, where scores=None
    # and ordering is pure proximity. Every changed oligo must therefore be audited, not assumed.
    cands, window, ins = evaluate_candidates(
        model, chrom_seq, residue, spec, usage=usage, allowed_guides=allowed_guides,
        scores=scores, slide_range=slide_range, max_distance=max_distance, flank=flank,
        credit_pam_separated=credit_pam_separated)
    usable = [c for c in cands if c["design"] is not None]
    # `sorted` is stable, so ties fall back to guides_near order -- deterministic, and the same
    # order `118` measured.
    usable.sort(key=lambda c: rank_key(c, lam))
    best = usable[0] if usable else None
    if best is not None:
        best["chosen"] = True
    design = None if best is None else best["design"]
    used = None if best is None else best["min_disruption_used"]
    if return_candidates:
        for c in cands:
            c.setdefault("chosen", False)
        return design, used, cands, window, ins
    return design, used
