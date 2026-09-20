"""
Guide-donor design policy for the AsLOV2 insertion library -- SINGLE PASS.

Kevin, 2026-08-27: *"a single, clean, updated version of the magestic guide-donor design code for
the insertion module arm of the package that will produce the best version of the library on one
go, without needing multiple passes."*

THE CONSTRUCTION, and the property it buys
------------------------------------------
    every capability ENLARGES the candidate set, and the choice is argmax of ONE scalar over it

so adding a capability cannot lower the maximum. "v2 is never worse than v1" is therefore a
property of the construction, checkable inside a single run
(`design_position(..., return_candidates=True)` exposes every variant it considered), rather than
something established by building the library twice and diffing it.

🔴 THE GUARANTEE IS IN THE NEW SCORE, NOT THE OLD ONE. v2 will produce designs that are worse under
the retired `q - 0.02*net_recoded` at some positions. That is intended: `648` measured that
lam = 0.02 is right at 4 recoded codons and wrong on both sides -- the first three codons are FREE
(implied lam is NEGATIVE) and ten codons cost 0.13/codon, six times what lam charged.

WHAT V1 DID AND WHY IT NEEDED TWO PASSES
----------------------------------------
Four places SUBSTITUTED a design instead of ADDING a candidate, and the ranking could veto on a
boolean (`planning/653` §2):

  1. `rank_key` was lexicographic -- (blocked, net_recoded >= CLIFF, -(q - lam*nr)). A q = 0.931
     design sorted BELOW a q = 0.476 one at YEL002C res283 because one mandatory cloning-site strip
     codon tripped `is_blocked` (`planning/643`).
  2. `credit_pam_separated` REPLACED a guide's design, so the uncredited design for that same guide
     became unreachable -- the shipped design was not in the new candidate set at all.
  3. `stepped_design` tried min_disruption 4 -> 3 -> 2 and RETURNED THE FIRST under the cliff.
  4. the donor slide returned the first that worked.

`GUIDE_LOCAL_FALLBACK` was invented to patch (2) by re-deriving and comparing -- the two-pass
pattern inside the designer. Here (1) is fixed by `score.design_score`, and (2) and (3) become
enumeration axes. (4) is NOT enumerated; see `SLIDE_IS_NOT_ENUMERATED` below.

Retired from v1, deliberately: `rank_key`, `LAMBDA_DEFAULT`, `CLIFF`, `is_blocked`,
`TWO_HURDLE_BLOCK`, `GUIDE_LOCAL_FALLBACK`, `stepped_design`. v1 is preserved in git at tag
`insertion-v1`; `git show insertion-v1:src/magestic/design/domain_insertion/policy.py` is the
whole of it, which is why no shadow copy is kept beside this file.
"""
import math
from dataclasses import dataclass, field

from .designer import (design_oligo, protospacer_span, spans_insertion,
                       insertion_separates_pam, alternatives, codons_in_window)
from .guides import guides_near
from . import score
from .score import design_score

MAX_DISTANCE = 30

# --- the indistinguishability band (Kevin, 2026-08-28) -------------------------------------------
# Two designs whose scores differ by less than this are treated as TIED and ordered by the explicit
# preference in `design_position` instead. `677`/`685` measured why this is necessary: `P_eff` takes
# only 4 distinct values and `P_dist` is flat inside the plateau -- both deliberately, so the
# fewest-codons rule can decide -- but `(1 - sv)` takes 4,635 values and therefore never ties. On
# 17.9 % of contested positions every informative factor tied and `sv` alone picked the winner,
# which produced 5,483 STRICTLY DOMINATED designs (`682`) and ~17,000 that took a lower-efficacy
# guide for no stated reason (`680`).
#
# 🔴 THE VALUE IS BOUNDED BY THE CALIBRATION, NOT CHOSEN FOR EFFECT. The smallest gap between
# adjacent `P_eff` levels is 0.8304 -> 0.8370, i.e. 0.79 %. Staying strictly below that guarantees
# the band can never swallow a distinction the calibration actually draws; it only catches designs
# the score genuinely cannot separate.
SCORE_TOLERANCE = 0.005          # relative, in LINEAR score

# --- the enumeration axes ------------------------------------------------------------------------
# Every combination is generated for every guide and scored. These are NOT settings.
CREDIT_VARIANTS = (False, True)          # `576`: credit a PAM-separated guide as self-disrupting
MIN_DISRUPTION_VARIANTS = (4, 3, 2)      # `112` clause 3 -- all of them, not first-success
JUNCTION_RECODE_VARIANTS = (False, True)  # `621`(C): recode a cloning site across an oligo junction

# 🔴 NOT ENUMERATED, stated so it is not mistaken for covered. The donor slide (17 offsets) is
# chosen inside `design_oligo` by `_slide_anchor_score`, which is a principled criterion (maximise
# uninterrupted distal homology), not a first-success accident. Enumerating it would multiply the
# search by ~17 for a axis that is already optimised. The argmax below therefore covers
# guide x credit x min_disruption x junction_recode, and NOT slide.
SLIDE_IS_NOT_ENUMERATED = True

# The Trp/Met missense bridge (`627`(E)) is a LAST RESORT and a DIFFERENT CONSTRUCT CLASS, never a
# silent winner: a lethal W->Y otherwise reads as "AsLOV2 not tolerated at this residue"
# (`621` §21, `151` §5). It is tried only when nothing else designs, and the result is labelled.
NO_SYN_CODONS = ("TGG", "ATG")


@dataclass
class Variant:
    """One reachable design for one guide, with how it was produced."""
    design: object
    protospacer: str
    guide: object
    credit: bool
    min_disruption: int
    junction_recode: bool
    wm_missense_bridge: bool = False
    score: float = float("-inf")
    wm_in_gap: int = 0
    chosen: bool = False

    @property
    def label(self):
        bits = []
        if self.credit:
            bits.append("credit")
        if self.junction_recode:
            bits.append("jr")
        if self.wm_missense_bridge:
            bits.append("wm")
        bits.append("md%d" % self.min_disruption)
        return "+".join(bits)


def wm_in_gap(design, window, ins, codons, usage):
    """UNRECODABLE codons in the wild-type run between the cassette and the nearest edit.

    This is the `102` mechanism and the input to `score.p_fidelity`: a Trp/Met that CANNOT be
    recoded forces a wild-type gap a conversion tract can stop inside, giving partial incorporation
    -- the cassette absent while the recoding is present. A gap the designer merely CHOSE not to
    fill is a different object; this counts only the ones it could not.
    """
    n_left, n_right = len(design.left_arm), len(design.right_arm)
    nat_l, nat_r = window[ins - n_left:ins], window[ins:ins + n_right]
    chg = [ins - 1 - j for j, (a, b) in enumerate(zip(design.left_arm[::-1], nat_l[::-1])) if a != b]
    chg += [ins + j for j, (a, b) in enumerate(zip(design.right_arm, nat_r)) if a != b]
    if not chg:
        return 0
    right = [i for i in chg if i >= ins]
    left = [i for i in chg if i < ins]
    if design.cut_distance >= 0:
        near = min(right) if right else None
        span = set(range(ins, near)) if near is not None else set()
    else:
        near = max(left) if left else None
        span = set(range(near + 1, ins)) if near is not None else set()
    if not span:
        return 0
    n = 0
    for c in codons:
        if c.is_split or not any(i in span for i in c.indices):
            continue
        cur = "".join(window[i] for i in c.indices).upper()
        if cur in NO_SYN_CODONS or not [a for a in alternatives(cur, usage) if a != cur]:
            n += 1
    return n


def design_variants(model, chrom_seq, residue, spec, *, usage, protospacer, scores, slide_range,
                    max_distance=MAX_DISTANCE, allow_wm_bridge=False):
    """EVERY design reachable for ONE guide, deduped by oligo.

    Replaces v1's `stepped_design`, which returned the FIRST min_disruption under the cliff and
    threw the rest away. Nothing is filtered here -- filtering is what made a capability able to
    make a position worse. Scoring happens in `design_position`.
    """
    out, seen = [], set()
    combos = [(c, md, jr)
              for c in CREDIT_VARIANTS
              for md in MIN_DISRUPTION_VARIANTS
              for jr in JUNCTION_RECODE_VARIANTS]
    if allow_wm_bridge:
        combos = [(c, md, jr) for c in CREDIT_VARIANTS
                  for md in MIN_DISRUPTION_VARIANTS for jr in (True,)]
    for credit, md, jr in combos:
        d = design_oligo(model, chrom_seq, residue, spec, usage=usage,
                         max_distance=max_distance, min_disruption=md,
                         allowed_guides={protospacer}, scores=scores, slide_range=slide_range,
                         credit_pam_separated=credit, junction_recode=jr,
                         wm_missense_bridge=(True if allow_wm_bridge else None))
        if d is None or d.oligo in seen:
            continue
        seen.add(d.oligo)
        out.append(Variant(design=d, protospacer=protospacer, guide=d.guide, credit=credit,
                           min_disruption=md, junction_recode=jr,
                           wm_missense_bridge=bool(allow_wm_bridge)))
    return out


def design_position(model, chrom_seq, residue, spec, *, usage, allowed_guides, scores,
                    slide_range, max_distance=MAX_DISTANCE, flank=200,
                    return_candidates=False, allow_wm_bridge=True):
    """Design one position: argmax of `score.design_score` over guides x variants.

    Returns (design, variant) -- or (design, variant, variants, window, ins) with
    `return_candidates`, which is what the property test and the figures read.
    """
    window, c_idx, _, cds_pos = model.window_with_coords(
        chrom_seq, model.cds_to_genomic(residue * 3), flank)
    ins = c_idx + 1
    codons = codons_in_window(cds_pos)

    variants = []
    for g in guides_near(window, ins, max_distance, guide_length=spec.guide_length):
        if allowed_guides is not None and g.protospacer not in allowed_guides:
            continue                                   # off-target screen: the precomputed set
        variants.extend(design_variants(
            model, chrom_seq, residue, spec, usage=usage, protospacer=g.protospacer,
            scores=scores, slide_range=slide_range, max_distance=max_distance))

    # LAST RESORT ONLY, and labelled: the W/M missense bridge is a different construct class.
    if not variants and allow_wm_bridge:
        for g in guides_near(window, ins, max_distance, guide_length=spec.guide_length):
            if allowed_guides is not None and g.protospacer not in allowed_guides:
                continue
            variants.extend(design_variants(
                model, chrom_seq, residue, spec, usage=usage, protospacer=g.protospacer,
                scores=scores, slide_range=slide_range, max_distance=max_distance,
                allow_wm_bridge=True))

    for v in variants:
        d = v.design
        v.wm_in_gap = wm_in_gap(d, window, ins, codons, usage)
        v.score = design_score(d.efficacy, d.sv_score, d.net_recoded, v.wm_in_gap,
                               cut_distance=d.cut_distance)

    best = None
    if variants:
        # TIE-BREAK: score first, then FEWEST recoded codons, then enumeration order.
        #
        # 🔴 The codon tie-break is not cosmetic, and 2026-08-28 made it LOAD-BEARING rather than
        # a corner case. `P_tract` is FLAT through the knot, so a 1-codon and a 3-codon design at
        # the same guide score IDENTICALLY -- and enumeration order alone would hand the win to
        # whichever axis happened to come first (it picked the 3-codon one).
        #
        # `score.p_eff` now quantises efficacy to the four levels its calibration can actually
        # resolve, and `score.p_dist` is FLAT through 5 nt by measurement. Both changes CREATE
        # ties on purpose, so this line decides 63.1 % of the below-knot positions where raw
        # efficacy previously decided 98.5 % of them on differences finer than the data
        # (`671`). Anything that reintroduces a spurious slope in either factor silently
        # disables it -- that is why `662`'s interpolated P_dist was not shipped (`670`).
        # Every recoded codon is a real edit whose risks the score does not fully see: off-target
        # consequences, a partially-copied multi-SNV codon that need not stay synonymous (`86`),
        # and simply more that can go wrong. Where the model is indifferent, do less.
        # This is Kevin's call at MRD1 res554 / WBP1 res283, 2026-08-27: "the insertion alone is
        # enough to disrupt this guide ... remove the syn codons altogether."
        #
        # `variants` is built in a fixed order (guides_near order, then the fixed combo order), so
        # the final fallback is deterministic.
        order = {id(v): i for i, v in enumerate(variants)}
        top = max(v.score for v in variants)
        # scores are logs, so a relative linear tolerance is an additive offset here
        floor = top + math.log(1.0 - SCORE_TOLERANCE)
        band = [v for v in variants if v.score >= floor]

        # 🔴 PLATEAU-AWARE, and the `max(..., PLATEAU)` is the whole trick. `score.p_dist` is FLAT
        # through PLATEAU nt BY MEASUREMENT (`670`: Cochran-Armitage z = -1.55, p = 0.12,
        # n = 1,777), so ranking on a cut difference INSIDE that band would assert a difference we
        # tested for and did not find -- the same error as `662`'s interpolated curve, just moved
        # into the tie-break. Clamping every sub-plateau cut to one value makes them tie, so they
        # fall through to raw efficacy; cuts that genuinely leave the plateau still rank.
        #
        # Kevin, 2026-08-28, choosing this over a plain "closest cut" second key: it "honors the
        # ordering wherever distance is a measured difference, and declines to rank on it where we
        # tested and found none." Measured, the two differ at 13,192 positions, 99.9 % of them
        # inside the plateau (`686`).
        #
        # RAW efficacy, not `P_eff`: inside the band the calibrated value carries no information --
        # that is what put the design in the band. Calibration strips efficacy's false SCALE, not
        # its ranking power (AUC 0.626 for a correct edit -- weak, but real and non-zero), so the
        # raw score is the right last resort once codons and real distance are spent.
        plateau = score.PDIST_BINS[0][0]

        def _prefer(v):
            d = v.design
            cut = abs(d.cut_distance)
            eff = -1.0 if d.efficacy is None else float(d.efficacy)
            return (d.net_recoded,          # 1. fewest recoded codons
                    max(cut, plateau),      # 2. closest cut, but only outside the flat plateau
                    -eff,                   # 3. highest RAW predicted efficacy
                    cut,                    # 4. closest cut, now including inside the plateau
                    order[id(v)])           # 5. deterministic fallback

        best = min(band, key=_prefer)
        best.chosen = True
    design = None if best is None else best.design
    if return_candidates:
        return design, best, variants, window, ins
    return design, best
