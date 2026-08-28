"""
Guide-donor oligo designer for AsLOV2 domain insertion.

Implements the five design considerations from the project brief:

  1. use a guide that passes off-target screening, favouring high-efficacy / low-SV guides
  2. synonymous recoding from the cut site toward the insertion site, to destroy the protospacer
     and remove microhomology that would let the donor strand-hop back to the genome
  3. synonymous removal of BspQI / PaqCI sites from the donor arms, or a +/-3 nt donor slide when
     the site sits at the donor-backbone junction
  4. the donor always starts on a codon boundary
  5. arms are GENOMIC (intron-containing) while codons are read from the SPLICED CDS

The output of `design_oligo` is intended to be handed straight to
`magestic.design.domain_insertion.validate.validate_oligo`, which re-derives everything
independently.
"""

from dataclasses import dataclass, field

from .genome import rc
from .guides import Guide, guides_near
from .spec import DISALLOWED_IN_DONOR, OligoSpec
from .validate import CODON_TABLE, translate

# aa -> synonymous codons
SYNONYMS = {}
for _c, _a in CODON_TABLE.items():
    SYNONYMS.setdefault(_a, []).append(_c)


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def build_codon_usage(cds_sequences) -> dict:
    """Codon -> relative frequency within its amino-acid family, counted from real CDSs."""
    counts = {}
    for seq in cds_sequences:
        s = seq.upper()
        for i in range(0, len(s) - 2, 3):
            c = s[i:i + 3]
            if c in CODON_TABLE:
                counts[c] = counts.get(c, 0) + 1
    usage = {}
    for aa, codons in SYNONYMS.items():
        tot = sum(counts.get(c, 0) for c in codons) or 1
        for c in codons:
            usage[c] = counts.get(c, 0) / tot
    return usage


@dataclass
class Codon:
    """A codon inside the design window."""

    residue: int              # 1-based residue number in the protein
    indices: tuple            # window indices of its three nucleotides (transcript order)

    @property
    def is_split(self) -> bool:
        return not (self.indices[1] == self.indices[0] + 1
                    and self.indices[2] == self.indices[1] + 1)

    def seq(self, window: str) -> str:
        return "".join(window[i] for i in self.indices)


def codons_in_window(cds_positions) -> list:
    """Group window indices into whole codons using their CDS positions."""
    by_residue = {}
    for idx, cp in enumerate(cds_positions):
        if cp is None:
            continue
        res = (cp - 1) // 3 + 1
        by_residue.setdefault(res, {})[(cp - 1) % 3] = idx
    out = []
    for res, slots in sorted(by_residue.items()):
        if len(slots) == 3:
            out.append(Codon(residue=res, indices=(slots[0], slots[1], slots[2])))
    return out


def alternatives(codon: str, usage: dict = None) -> list:
    """Synonymous codons, most-different-first, breaking ties by codon usage."""
    aa = CODON_TABLE.get(codon.upper())
    if aa is None:
        return []
    alts = [c for c in SYNONYMS.get(aa, []) if c != codon.upper()]
    return sorted(alts, key=lambda c: (-hamming(c, codon.upper()),
                                       -(usage or {}).get(c, 0.0)))


@dataclass
class DesignResult:
    orf: str
    residue: int                      # insertion goes AFTER this residue
    oligo: str
    guide: Guide
    left_arm: str
    right_arm: str
    recoded: list = field(default_factory=list)   # (residue, from_codon, to_codon)
    donor_shift: int = 0
    notes: list = field(default_factory=list)
    offtargets: int = None
    disruption: int = 0
    efficacy: float = None
    sv_score: float = None
    donor_starts_in_noncoding: bool = False
    left_cds_positions: list = field(default_factory=list)
    right_cds_positions: list = field(default_factory=list)
    # `recoded` is an OP LOG -- it records every apply(), including round-trips, so it overstates
    # how many codons actually differ from the genome (22.4% of rows did, pre-fix). `net_recoded`
    # is the count that actually differ. Panel finding M4 (2026-07-26).
    net_recoded: int = 0

    @property
    def cut_distance(self) -> int:
        return self._cut_distance

    _cut_distance: int = 0


def _guide_site_present(seq: str, guide: Guide) -> bool:
    """Is the protospacer + a functional NGG PAM still present on either strand?"""
    p = guide.protospacer
    for hay in (seq.upper(), rc(seq.upper())):
        i = hay.find(p)
        while i != -1:
            pam = hay[i + len(p):i + len(p) + 3]
            if len(pam) == 3 and pam[1:] == "GG":
                return True
            i = hay.find(p, i + 1)
    return False



def near_cognate_recut(seq: str, guide: Guide, max_mm: int = 3, seed_len: int = 12,
                       max_seed_mm: int = 1) -> bool:
    """
    Is there a site in `seq` that SpCas9 would still plausibly cut with this guide?

    `_guide_site_present` is exact-match only, so a repaired allele carrying 1-3 mismatches but a
    functional NGG and an intact seed slips through both it and validator C8. Cas9 tolerates
    PAM-distal mismatches; what it does not tolerate is mismatches in the PAM-proximal ~12 nt.

    Adversarial panel finding M5 (2026-07-26): 12/4,795 designs left a non-self-disrupting repaired
    allele with a functional NGG and <=1 mismatch in the seed (5 with ZERO), and a further
    **142 (~3%)** had the RIGHT_LINKER's terminal A reconstitute the guide's PAM-distal base across
    the `TCCGGA|right_arm` junction. Consequence is a silent false negative: the correctly edited
    chromosome is re-cut beside the fresh cassette, so the barcode is recovered from the pool but the
    intended optogenetic allele is gone.

    Returns True if any 20-mer + NGG lies within `max_mm` total mismatches and at most
    `max_seed_mm` within the PAM-proximal `seed_len` bases.
    """
    p = guide.protospacer.upper()
    L = len(p)
    for hay in (seq.upper(), rc(seq.upper())):
        for i in range(0, len(hay) - L - 2):
            pam = hay[i + L:i + L + 3]
            if len(pam) < 3 or pam[1:] != "GG":
                continue
            cand = hay[i:i + L]
            mm = [k for k in range(L) if cand[k] != p[k]]
            if len(mm) > max_mm:
                continue
            # seed = the PAM-proximal end, i.e. the LAST `seed_len` bases of the protospacer
            if sum(1 for k in mm if k >= L - seed_len) > max_seed_mm:
                continue
            return True
    return False

def _oligo_sites_are_canonical(oligo: str, spec: OligoSpec) -> bool:
    """
    The assembled oligo must carry EXACTLY the intended cloning sites:
    one BspQI (in the cloning site) and one inward PaqCI pair (in the stuffer). Anything else
    means a site was created across a junction.
    """
    from .spec import ASCI_SITE, BSPQI_SITE, NOTI_SITE, PAQCI_SITE

    o = oligo.upper()
    n_bspqi = o.count(BSPQI_SITE) + o.count(rc(BSPQI_SITE))
    if n_bspqi != spec.expected_bspqi_sites:
        return False
    if o.count(PAQCI_SITE) != 1 or o.count(rc(PAQCI_SITE)) != 1:
        return False
    # AscI appears once in the constant 5' flank by design; NotI never.
    if o.count(NOTI_SITE) or o.count(rc(NOTI_SITE)):
        return False
    if (o.count(ASCI_SITE) + o.count(rc(ASCI_SITE))) > (
            spec.fwd_const.upper().count(ASCI_SITE)
            + spec.fwd_const.upper().count(rc(ASCI_SITE))):
        return False
    return True


# --- `627` (E) W/M conservative-missense bridge (Kevin, 2026-08-22). DEFAULT OFF. ----------------
# The 2,151 shipped pool-1 W/M rescue oligos come from `154_wm_rescue_production.py`, whose work
# list is `designed==yes && trp_met_blocked==1` -- an ADDITIVE PASS OVER DESIGNS THAT ALREADY EXIST.
# A position whose connectivity gate fails outright returns None, is written `designed==no`, and the
# rescue pass never sees it. The escape hatch is real but wired one stage too late (`617`: 153 of
# 159 DESIGN_GAINED are connectivity failures, `codons_no_synonym > 0` in all 153).
#
# This puts the same substitutions INSIDE the gate. Conservative, single-codon, and the only two
# single-codon amino acids there are:
#     Met ATG -> Leu (TTG / CTG)      Trp TGG -> Tyr (TAT / TAC)
#
# 🔴 LABELLING IS LOAD-BEARING (`151` §5). A design that used this carries a DELIBERATE MISSENSE and
# must NEVER be emitted as `construct_class = internal`, or it is analysed as a clean insertion
# allele. `WM_MISSENSE_NOTE` is appended to `notes` so the builder can route it, and every such
# design should ship paired with its syn-only control (`550`).
WM_MISSENSE_BRIDGE = False
WM_MISSENSE_NOTE = "bridge closed by conservative missense (627E)"
WM_MISSENSE_ALTS = {"ATG": ("TTG", "CTG"), "TGG": ("TAT", "TAC")}


def _wm_bridge_alts(cur, usage):
    """Synonymous alternatives first; conservative missense ONLY for Met/Trp and ONLY when the flag
    is on and no synonymous option exists at all."""
    alts = [a for a in alternatives(cur, usage) if a != cur]
    if alts or not WM_MISSENSE_BRIDGE:
        return alts, False
    return list(WM_MISSENSE_ALTS.get(cur.upper(), ())), True


# --- `621` (C) junction-site recoding (Kevin, 2026-08-22). DEFAULT OFF. --------------------------
JUNCTION_RECODE = False      # flip to True to recode a junction cloning site instead of refusing

# How far into an arm a junction-spanning site can reach. The longest DISALLOWED_IN_DONOR site is
# 8 nt (AscI/NotI), so a site straddling a junction can occupy at most its first/last 7 nt.
JR_REACH = 8


def _jr_assemble(spec, guide, left_arm, right_arm):
    return (spec.fwd_const + guide.protospacer + spec.cloning_site
            + left_arm + spec.left_linker + spec.stuffer + spec.right_linker
            + right_arm + spec.rev_primer)


def _recode_junction_site(wl, codons, lo, hi, ins, spec, guide, usage, recoded, wl0):
    """Remove a cloning site created ACROSS an oligo junction by recoding ONE arm codon.

    Returns (left_arm, right_arm, oligo) on success, or None.  Mutates `wl` and appends to
    `recoded` ONLY on success -- a failed attempt leaves both untouched.

    Scope, and why each restriction is there:
      * only codons lying wholly inside an arm and within JR_REACH of an arm terminus -- those are
        the only ones whose bases can participate in a junction-spanning site;
      * never a split codon (recoding one would need both exon halves);
      * a codon overlapping the protospacer+PAM is allowed ONLY when it is CONTIGUOUS with the
        insertion (gap <= MAX_GAP_DEFAULT) -- Kevin's "with us" rule, 2026-08-22. 🔴 `622` measured
        why the exemption is REQUIRED, not optional: for a PAM-separated minus-strand guide the
        insertion abuts the protospacer, so the arm-start codon -- the `AGA` Arg that is the whole
        fix -- lies INSIDE the protospacer. Excluding it outright made the fix fire at other guides
        and never at the one that needs it (3 of 6 DESIGN_LOST recovered nothing; the 3 that did
        cost 11/6/10 codons against a shipped 3/3/2). An edit at gap 0 is connected to the cassette
        by construction, so it clears the two-hurdle test -- it cannot be copied without the
        cassette. A protospacer codon FURTHER out stays excluded: it would rescore disruption after
        the M4 re-gate and land in the hazard window;
      * never a codon already carrying an introduced edit -- the WT guard from `_strip_disallowed`
        and `apply()`: reverting one would silently undo a counted disruption;
      * single-SNV alternatives first (Kevin, 2026-07-29: a partially-copied multi-base codon need
        not be synonymous), then the rest in `alternatives()` order;
      * the result must leave BOTH arms free of DISALLOWED_IN_DONOR and pass
        `_oligo_sites_are_canonical`.
    """
    p_lo, p_hi = protospacer_span(guide, spec.guide_length)

    def _gap_from_cassette(c):
        """Wild-type bases between the insertion junction and this codon's nearest base."""
        return min((i - ins) if i >= ins else (ins - 1 - i) for i in c.indices)

    def _protospacer_ok(c):
        if not any(p_lo <= i < p_hi for i in c.indices):
            return True
        return _gap_from_cassette(c) <= MAX_GAP_DEFAULT      # "with us": connected by construction

    def _touches_terminus(c):
        idx = c.indices
        return (min(idx) < lo + JR_REACH or max(idx) >= ins - JR_REACH
                or min(idx) < ins + JR_REACH or max(idx) >= hi - JR_REACH)

    cands = [c for c in codons
             if not c.is_split
             and all(lo <= i < hi for i in c.indices)
             and _protospacer_ok(c)
             and all(wl[i] == wl0[i] for i in c.indices)
             and _touches_terminus(c)]
    # nearest-to-a-terminus first, so the shortest possible reach into the arm is tried first
    cands.sort(key=lambda c: min(min(abs(i - b) for b in (lo, ins - 1, ins, hi - 1))
                                 for i in c.indices))

    for codon in cands:
        cur = codon.seq("".join(wl))
        alts = [a for a in alternatives(cur, usage) if a != cur]
        alts.sort(key=lambda a: sum(1 for x, y in zip(a, cur) if x != y))   # single SNV first
        for alt in alts:
            for k, idx in enumerate(codon.indices):
                wl[idx] = alt[k]
            left_arm = "".join(wl[lo:ins])
            right_arm = "".join(wl[ins:hi])
            clean = all(site not in seg and rc(site) not in seg
                        for seg in (left_arm, right_arm) for site in DISALLOWED_IN_DONOR)
            if clean:
                oligo = _jr_assemble(spec, guide, left_arm, right_arm)
                if _oligo_sites_are_canonical(oligo, spec):
                    recoded.append((codon.residue, cur, alt))
                    return left_arm, right_arm, oligo
            for k, idx in enumerate(codon.indices):     # revert
                wl[idx] = cur[k]
    return None


def spans_insertion(guide: Guide, guide_length: int, insertion: int) -> bool:
    """
    Does the insertion point fall INSIDE the protospacer + PAM?

    If so the guide is **self-disrupting**: a 420-nt cassette dropped into the middle of the
    protospacer destroys the site outright, in the repaired allele and in the donor alike (the
    donor's copy is split across the PaqCI stuffer). No synonymous recoding is needed to stop
    re-cutting, which matters most exactly where recoding is impossible -- protospacers lying in
    5' UTR or intron, where there are no codons to change.

    Kevin's stated priority (2026-07-25): "favour guides that span both donor arms (i.e. the
    insertion directly disrupts it), then look at guides on each side. If one side has a blocker,
    we look to the other side."
    """
    # The insertion must split the PROTOSPACER, not merely the PAM. Splitting only the PAM is
    # not disruptive here: the G-S linker `GGTTCT` begins with GG, so a sense guide whose
    # protospacer ends at the arm/linker boundary has a functional NGG PAM RECONSTITUTED by the
    # linker, and the donor would be cut. Caught by validator check C8.
    lo, hi = protospacer_span(guide, guide_length)
    if guide.strand == "+":
        p_lo, p_hi = lo, hi - 3           # drop the PAM, which sits 3' of the protospacer
    else:
        p_lo, p_hi = lo + 3, hi           # antisense: PAM sits 5' in window coordinates
    return p_lo < insertion < p_hi


def insertion_separates_pam(guide: Guide, guide_length: int, insertion: int) -> bool:
    """Does the insertion fall BETWEEN an intact protospacer and its PAM?

    `spans_insertion` requires the insertion to split the PROTOSPACER. But a guide whose protospacer
    lies wholly on one side of the insertion with its PAM wholly on the other is *also* destroyed
    structurally: 408 nt of AsLOV2 (or the 38-nt stuffer in the donor) lands between spacer and PAM,
    and Cas9 cannot cut a spacer whose PAM is a cassette away. Kevin identified this on 2026-08-19
    from the bovine TdT residue-271 alignment.

    ONLY the exact-abutment case is reported, and that is not a stylistic choice. The LEFT linker is
    `GGTTCT`, so it begins `GG`:
        insertion at the protospacer 3' end -> the 3 nt after the spacer become `GGT` -> not NGG
        insertion 1 nt into the PAM         -> they become  N + `GG`                  -> IS an NGG
        insertion 2 nt into the PAM         -> they become  N,G + `G`                 -> IS an NGG
    One nucleotide further and the linker RECONSTITUTES a functional PAM (the hazard `spec.py`
    documents as `creates_junction_pam`). Crediting those would licence skipping the very recoding
    that then prevents a re-cut.

    NECESSARY, NOT SUFFICIENT: callers must still clear `_guide_site_present` on both the repaired
    allele and the donor, and `near_cognate_recut` on the repaired allele. Those gates are unchanged
    and run regardless of this predicate.
    """
    lo, hi = protospacer_span(guide, guide_length)
    if guide.strand == "+":
        # spacer [lo, hi-3), PAM [hi-3, hi). Safe iff the insertion abuts the spacer 3' end.
        return insertion == hi - 3
    # antisense: PAM [lo, lo+3), spacer [lo+3, hi). Safe iff the insertion abuts the spacer 5' end.
    return insertion == lo + 3


def protospacer_span(guide: Guide, guide_length: int) -> tuple:
    """[start, end) window indices covered by the protospacer AND its PAM."""
    if guide.strand == "+":
        return guide.pam_index - guide_length, guide.pam_index + 3
    return guide.pam_index, guide.pam_index + 3 + guide_length


def disruption_score(original: str, edited: str, degenerate_offset: int = None) -> int:
    """
    Mismatches introduced across the protospacer + PAM.

    The legacy MAGESTIC saturation-editing code requires
    MIN_GUIDE_DONOR_DISRUPTION_SCORE = 6. This project uses **4** by default (Kevin, 2026-07-25):

        "we may be able to reduce the total disruption score by only requiring ~4 positions in the
         guide disrupted, especially if we are close to the PAM and cut site (i.e. the seed region).
         Adding additional positions decreases the length of homology to mediate repair, and
         increases the chances of the edit only incorporating the guide-disruptive changes but not
         the changes closer to the insertion."

    That last clause is the real hazard: guide-disruptive changes sit next to the cut, the cassette
    sits further away, so a donor loaded with changes near the cut offers an attractive partial
    repair product that carries the disruption but NOT the insert.

    **`degenerate_offset` masks SpCas9's degenerate PAM base (the N of NGG).** A synonymous change
    there cannot stop re-cutting, so counting it inflates the score and lets phase 2 stop early with
    the seed untouched. Adversarial panel finding M1 (2026-07-26): 604/4,795 designs reported
    `disruption = 4` with a fully intact NGG and only 3 informative mismatches. Use
    `degenerate_offset_in_span()` to compute it -- the offset is strand-dependent.
    """
    n = sum(a != b for a, b in zip(original.upper(), edited.upper()))
    if degenerate_offset is not None and 0 <= degenerate_offset < min(len(original), len(edited)):
        if original[degenerate_offset].upper() != edited[degenerate_offset].upper():
            n -= 1
    return n


def degenerate_offset_in_span(guide: Guide, guide_length: int) -> int:
    """
    Offset, WITHIN `protospacer_span`, of the degenerate PAM base (the N of NGG).

    `protospacer_span` returns the protospacer followed by its PAM for a sense guide, and the PAM
    followed by the protospacer for an antisense one (window coordinates always run left-to-right
    in gene sense). So the N sits at:
      '+'  offset `guide_length`  -- first base of the trailing NGG
      '-'  offset 2               -- the antisense PAM reads CCN on this strand, N last
    """
    return guide_length if guide.strand == "+" else 2


def seed_anchor(guide: Guide) -> int:
    """
    Window index of the PAM-proximal END of the protospacer -- the seed.

    Phase 2 sorts candidate codons by distance to this. The previous anchor was
    `guide.pam_index + 3`, i.e. one base PAST the PAM, which let a codon lying entirely outside the
    protospacer rank first and let phase 2 spend its whole budget on non-contributing positions
    (adversarial panel finding M2, 2026-07-26; YDR236C res 8 recoded codons 16 and 15, contributing
    0 and 0 informative disruption, then stopped with the seed untouched).
    """
    return guide.pam_index - 1 if guide.strand == "+" else guide.pam_index + 3


# Toggle for the single-SNV recoding preference; see planning/86 for the measured trade-off.
PREFER_SINGLE_SNV = False

# Largest wild-type run tolerated between adjacent introduced mutations (and between the cassette
# and the nearest edit). Artifact 10's co-incorporation curve: P(complete) = 94.7% at 0 nt, 87.4%
# at 1, 81.9% at 2, 77.1% at 3, 75.3% at 4, then a cliff to 25.9% at 5. A value of 3 therefore
# ACCEPTS ~23% incomplete co-incorporation at the worst permitted gap -- exactly where a
# conversion-tract boundary can drop the cassette while keeping the cut-destroying mutations.
# (Kevin, 2026-07-29, on SPP381 codon 17: "why did we not recode the E ... This would eliminate the
# dangerous 3 bp homology." It was left because the gap equals max_gap, not because E is unrecodable
# -- GAA->GAG is a single SNV.) Swept in planning/89.
MAX_GAP_DEFAULT = 3


def _recode_to_destroy(window_list, codons, guide, lo, hi, usage, recoded, spec,
                       min_disruption: int = 4, recode_cut_to_insertion: bool = True,
                       insertion: int = None, max_gap: int = None, seed_depth: int = 12,
                       protect_boundary_codons: bool = False,
                       spread_disruption_target: int = 5,
                       credit_pam_separated: bool = False):
    """
    Tile synonymous substitutions ANCHORED AT THE INSERTION POINT, working outward through the
    cut and a short way into the guide.

    Policy (Kevin, 2026-07-25). The region that matters most is BETWEEN the domain-insertion site
    and the cut, not the protospacer:

      * The reconstructed pair table (artifact 10) gives P(complete co-incorporation) vs the
        WT-identical gap between adjacent introduced mutations: 94.7% at 0 nt, 87.4% at 1,
        81.9% at 2, 77.1% at 3, 75.3% at 4, then a cliff to 25.9% at 5 and ~17% at 6+.
      * The reference POL3 design recodes only codons lying under the protospacer, leaving a
        6-nt perfect-homology gap between the AsLOV2 cassette and the nearest edit. That block of
        synonymous changes sits beside the cut and will be copied; the cassette is 6 nt further
        away across perfect homology and is the part at risk of being dropped -- P(complete)
        ~17.6%. The failure mode is a clone carrying the recoding but NOT the domain.
      * Recoding every codon under the guide is unnecessary; what is needed is enough disruption
        to stop re-cutting.

    So: phase 1 fills the insertion->cut(+seed) window so no gap exceeds `max_gap`; phase 2 adds
    protospacer edits, PAM-proximal first, only until the disruption threshold is met.
    """
    p_lo, p_hi = protospacer_span(guide, spec.guide_length)
    self_disrupting = (insertion is not None
                       and spans_insertion(guide, spec.guide_length, insertion))
    # OFF BY DEFAULT so every shipped build reproduces byte-identically (`320`/`312` rule:
    # parameterise the tested code path, never fork it). When on, a PAM-separated guide counts as
    # self-disrupting, which skips PHASE 2 ONLY; the site-presence and near-cognate gates below are
    # untouched and still decide.
    if credit_pam_separated and insertion is not None and not self_disrupting:
        self_disrupting = insertion_separates_pam(guide, spec.guide_length, insertion)
    p_lo, p_hi = max(p_lo, 0), min(p_hi, len(window_list))
    original = "".join(window_list[p_lo:p_hi])
    window0 = list(window_list)          # pristine copy, to locate introduced mutations

    def usable(c):
        if protect_boundary_codons and insertion is not None:
            # The two codons flanking the insertion. Leaving them native keeps aligners from
            # opening the insertion gap mid-codon (cosmetic), at the cost of a larger
            # WT-identical gap beside the cassette (real -- artifact 10).
            if any(abs(i - insertion) <= 3 for i in c.indices):
                return False
        if c.is_split or not all(lo <= i < hi for i in c.indices):
            return False
        # Same WT guard as `apply()`. Without it a codon whose ONLY synonym is the wild-type codon
        # passes `usable` and then no-ops in `apply` -- harmless but it makes the phase loops spend
        # iterations on candidates that can never contribute (panel note on M3).
        wt = "".join(window0[i] for i in c.indices)
        return bool([a for a in alternatives(c.seq("".join(window_list)), usage) if a != wt])

    def apply(c):
        cur = c.seq("".join(window_list))
        # WILD-TYPE GUARD. `alternatives()` excludes only the CURRENT codon, so the pristine genomic
        # codon is always a legal return -- and for any family whose synonyms are all Hamming-1 it is
        # frequently the top pick. Phase 2 would then re-visit a codon phase 1 had already recoded and
        # put it straight back (`alternatives('CCT')[0] == 'CCA'` and `alternatives('CCA')[0] == 'CCT'`
        # is a deterministic two-cycle), restoring perfect homology beside the cassette -- the exact
        # harm `protect_boundary_codons=False` exists to avoid. Invisible downstream: a reverted codon
        # is synonymous, in-frame and codon-phased, so C1-C8 all pass.
        # Adversarial panel finding M3 (2026-07-26): 437/4,795 designs ended with a >=5 nt WT-identical
        # stretch beside the cassette; with this guard that falls to 52, mean gap 1.36 -> 0.95, and
        # ZERO designable positions are lost.
        wt = "".join(window0[i] for i in c.indices)
        alts = [a for a in alternatives(cur, usage) if a != wt]
        if not alts:
            return False
        # PREFER_SINGLE_SNV (Kevin, 2026-07-29). A codon changed at 2-3 positions can be PARTIALLY
        # copied by a conversion tract, and a partial copy need not be synonymous -- measured, 74%
        # of multi-SNV codons have a non-synonymous partial state. A 1-SNV change has no partial
        # state. Off by default because it TRADES AGAINST disruption density: fewer mismatches per
        # codon can force recoding more codons to reach `min_disruption`, pushing edits further from
        # the insertion. Measure before enabling (see planning/86).
        if PREFER_SINGLE_SNV:
            alts = sorted(alts, key=lambda a: (sum(1 for k in range(3) if a[k] != cur[k]),
                                               alts.index(a)))
        for k, idx in enumerate(c.indices):
            window_list[idx] = alts[0][k]
        recoded.append((c.residue, cur, alts[0]))
        return True

    def edit_marks(cassette_mark=None):
        """
        Indices of introduced mutations so far, plus a mark standing for the cassette itself.

        `insertion` is a JUNCTION index (between W[insertion-1] and W[insertion]), not a base, so it
        cannot be used as a mark on both sides. Using it raw put the mark at `w_lo` for direction>0
        (counted) but at `w_hi+1` for direction<0 (filtered out by `_max_gap_in`), so the LEFT side
        never measured the cassette-to-nearest-edit gap while the right side did -- the second half
        of adversarial panel finding M6 (2026-07-26). Callers pass the adjacent BASE on the cut side.
        """
        marks = {i for i, (a, b) in enumerate(zip(window_list, window0)) if a != b}
        if cassette_mark is not None:
            marks.add(cassette_mark)
        elif insertion is not None:
            marks.add(insertion)
        return marks

    ins = insertion if insertion is not None else guide.cut_junction
    cut = guide.cut_junction
    if max_gap is None:
        max_gap = MAX_GAP_DEFAULT
    direction = 1 if cut >= ins else -1

    # Window: from the insertion point, through the cut, `seed_depth` nt into the guide -- and
    # STRICTLY ON THE CUT SIDE. Recoding the far side buys nothing: the invading strand travels
    # from the cut toward the insertion, so only homology on that path can let it switch back.
    # Changes beyond the insertion just shorten the homology available for repair.
    # (Kevin, 2026-07-25, on NBL1: "I don't see any reason to recode the Lys codon on the right
    # side of the insertion given that the guide + PAM are entirely confined to the left side.")
    # The bounds are half-open at `ins` so the first far-side codon cannot sneak in -- an earlier
    # version used <= on both ends and recoded exactly that codon.
    #
    # The gap-fill window STOPS AT THE CUT and never runs past it (`seed_depth` extends the
    # *phase 2* protospacer search, not this one).
    # (Kevin, 2026-07-25, on YOS1: "When we are already obtaining >4 disruptions on the domain
    # insertion side of the guide, it does not help us to add additional synonymous recodings on
    # the other side. In fact, it could be detrimental as we can have microhomologies that might
    # only include the recodings.")
    # The failure mode is specific and silent: a recoding sitting BEYOND the cut can be copied by
    # a conversion tract that never reaches the insertion point. That tract destroys the
    # protospacer -- so the locus is never re-cut -- while the cassette is not installed. The clone
    # then looks quietly WT at the protein level and is a false negative in the screen, not a
    # visible failure. Edits between the insertion and the cut cannot do this: any tract covering
    # them lies on the path the invading strand takes to the insertion.
    #
    # M6 (adversarial panel, 2026-07-26): these two bounds were ASYMMETRIC. Bases strictly between
    # the insertion and cut junctions are [ins, cut-1] going right and [cut, ins-1] going left, so
    # the old right bound `cut` included one base PAST the cut while the left bound did not. With
    # max_gap=3 that made phase-1 gap-fill fire at cut_distance >= +4 on the right but only <= -5 on
    # the left -- identical geometry, different design, affecting 7.6% of the background stratum.
    #
    # DECISION (2026-07-27): make them symmetric by WIDENING THE LEFT, not narrowing the right.
    # Narrowing the right to `cut-1` loses designable positions (3 flipped True->False in the panel's
    # test) AND enlarges the WT-identical gap beside the cassette, which artifact 10's
    # co-incorporation curve says is the harmful direction. The right bound was arguably correct all
    # along and the left was under-recoding, so both sides now include the base immediately beyond
    # the cut. Revisit with Kevin if the strict "strictly between the junctions" convention is wanted.
    if max_gap is None:
        max_gap = MAX_GAP_DEFAULT
    if direction < 0:
        w_lo, w_hi = cut - 1, ins - 1                # cut side is to the LEFT of the insertion
    else:
        w_lo, w_hi = ins, cut                        # cut side is to the RIGHT
    w_lo, w_hi = max(w_lo, lo), min(w_hi, hi)

    _deg = degenerate_offset_in_span(guide, spec.guide_length)

    def score_now():
        # `_deg` masks the N of NGG -- a change there cannot prevent re-cutting (M1).
        return disruption_score(original, "".join(window_list[p_lo:p_hi]), degenerate_offset=_deg)

    # --- phase 1: SPREAD from the insertion TOWARDS the cut -----------------
    #
    # ⛔ BUG FIXED 2026-07-29 (Kevin, from the SPP381/YBR152W spot check):
    #   "There should always be codon spreading from the insertion site TOWARDS the guide cut
    #    site, or until a disruption score of ~5 is achieved, NOT the other way around."
    #
    # The old break condition was `_max_gap_in(marks, w_lo, w_hi) <= max_gap`. `_max_gap_in`
    # measures gaps only BETWEEN adjacent marks and uses w_lo/w_hi purely as FILTERS, never as
    # sentinels -- so the distance from the outermost edit to the CUT was invisible to it. One
    # codon of three adjacent edits therefore gave max-gap 0 and broke the loop immediately,
    # leaving everything from ins+3 to the cut wild-type. Measured on SPP381 (cut at +16): edits
    # at offsets 0-2 only, then a 9-nt WT run to the cut.
    #
    # Why that is dangerous, not merely untidy: a conversion tract whose boundary falls in that
    # WT run copies the cut-proximal mutations WITHOUT the cassette. Those mutations destroy the
    # protospacer, so the locus is never re-cut, and the clone is a SILENT false negative -- it
    # survives, resists re-cutting, and scores as edited by barcode while carrying no domain.
    # This is the same failure mode the docstring already warns about for edits beyond the cut.
    #
    # The fix adds the cut-side window edge as a SENTINEL, so the loop keeps spreading until the
    # outermost edit is within `max_gap` of the cut, with an early stop once disruption reaches
    # `spread_disruption_target` (Kevin's "~5").
    if recode_cut_to_insertion:
        cut_edge = w_hi if direction > 0 else w_lo

        def _gap_to_cut():
            marks = edit_marks(cassette_mark=(ins if direction > 0 else ins - 1))
            return _max_gap_in(set(marks) | {cut_edge}, w_lo, w_hi)

        span = sorted((c for c in codons if usable(c)
                       and all(w_lo <= i <= w_hi for i in c.indices)),
                      key=lambda c: min(abs(i - ins) for i in c.indices))
        for c in span:
            if score_now() >= spread_disruption_target:
                break
            if _gap_to_cut() <= max_gap:
                break
            apply(c)


    # --- phase 2: top up disruption, PAM-proximal first ---------------------
    if not self_disrupting and (
            score_now() < min_disruption or _guide_site_present("".join(window_list), guide)):
        anchor = seed_anchor(guide)

        def on_cut_side(c):
            """
            Phase 2 is restricted to the cut side of the insertion, exactly as phase 1 is.

            Without this, a protospacer that straddles the insertion lets an edit land on the FAR
            side of the cassette from the cut. Such an edit is reachable by a conversion tract that
            never crosses the insertion: the tract destroys the protospacer -- so the locus is never
            re-cut -- while the cassette is not installed, giving a clone that looks quietly WT and
            scores as a false negative rather than a visible failure.
            (Kevin, 2026-07-25, YOS1: "it does not help us to add additional synonymous recodings on
            the other side. In fact, it could be detrimental as we can have microhomologies that
            might only include the recodings.")
            Caught at YER074W-A residue 26, where codon 27 was recoded despite the cut lying 3 nt to
            the LEFT of the insertion.
            """
            return all(i < ins for i in c.indices) if direction < 0 \
                else all(i >= ins for i in c.indices)

        proto = sorted((c for c in codons if usable(c) and on_cut_side(c)
                        and any(p_lo <= i < p_hi for i in c.indices)),
                       key=lambda c: min(abs(i - anchor) for i in c.indices))
        for c in proto:
            if score_now() >= min_disruption and not _guide_site_present(
                    "".join(window_list), guide):
                break
            apply(c)

    score = score_now()
    if self_disrupting and not _guide_site_present("".join(window_list), guide):
        # The cassette interrupts the protospacer, so disruption is structural rather than
        # sequence-based and the score threshold does not apply. The site-presence check is kept
        # as a guard: a linker or junction can occasionally reconstitute a functional PAM.
        return True, score
    still = _guide_site_present("".join(window_list), guide)
    return (not still) and score >= min_disruption, score


def _max_gap_in(marks, w_lo, w_hi):
    """Largest run of WT-identical positions between adjacent introduced mutations in [w_lo, w_hi]."""
    pts = sorted(m for m in marks if w_lo <= m <= w_hi)
    if len(pts) < 2:
        return w_hi - w_lo
    return max(b - a - 1 for a, b in zip(pts, pts[1:]))



def _slide_anchor_score(window, lo, ins, hi):
    """
    How much UNINTERRUPTED homology would survive site-stripping for this arm placement?

    A disallowed cloning site inside an arm has to be recoded away. If it sits near the arm's OUTER
    end, that recoding destroys the distal anchor -- the uninterrupted match the invading strand
    needs. Artifact 36 measured the consequence: **100% of `site_EDGE` designs ended below 30 nt of
    effective homology**, versus ~1% of background positions.

    But the +/-3 nt donor slide moves the ARM BOUNDARY. A site that sits just inside an arm at one
    shift often falls just OUTSIDE it at another -- no recoding needed, anchor intact. The slide loop
    used to try shifts in a fixed order (0, 3, -3, 6, -6) and take the first that merely *worked*,
    so it would happily recode a site at the arm edge when a one-codon slide would have excluded it.

    Returns min(left_anchor, right_anchor), where each anchor is the distance from that arm's outer
    end to the first disallowed site inside it (or the full arm length if the arm is clean). Higher
    is better; callers try shifts in descending order.
    """
    seq = "".join(window).upper()
    left, right = seq[lo:ins], seq[ins:hi]

    def first_site(arm, from_outer_end):
        hits = []
        for site in DISALLOWED_IN_DONOR:
            for probe in (site, rc(site)):
                j = arm.find(probe)
                while j != -1:
                    hits.append((j, j + len(probe)))
                    j = arm.find(probe, j + 1)
        if not hits:
            return len(arm)
        # left arm: outer end is index 0, so the anchor is the distance to the nearest site start.
        # right arm: outer end is the last base, so measure back from len(arm).
        if from_outer_end == "start":
            return min(a for a, _ in hits)
        return min(len(arm) - b for _, b in hits)

    return min(first_site(left, "start"), first_site(right, "end"))


def _connect_edits_to_insertion(window_list, window0, codons, lo, hi, usage, recoded,
                                ins, cut, guide, spec, max_gap=None, diag=None):
    """
    Guarantee that every mutation INSIDE the protospacer/PAM is connected back to the cassette.

    Kevin, 2026-07-29: *"if an RE site exists within the guide or PAM, [check] that all the codons
    between the RE-site-stripping SNV and the insertion are recoded."*

    Why this is needed even though phase 1 already spreads to the cut: `_strip_disallowed` runs
    AFTER the recoding phases and can introduce a fresh mutation anywhere in an arm, including deep
    inside the protospacer. SPP381 (YBR152W res 13) is the worked example -- two BspQI sites in the
    right arm, one of them at protospacer offsets 13-20, stripped by recoding codon 18 while codons
    15-17 stayed wild-type.

    The hazard is specific: a mutation inside the protospacer destroys the cut site. If a gene
    conversion tract can copy it WITHOUT reaching the cassette -- which a long wild-type run between
    the two permits -- the locus becomes uncuttable while carrying no domain. That clone survives,
    is never re-cut, and scores as edited by barcode: a SILENT false negative.

    Mutations OUTSIDE the protospacer+PAM are deliberately NOT connected. They cannot render the
    site uncuttable on their own, so copying them without the cassette leaves the locus still
    cuttable -- a visible failure, not a silent one. (SPP381's second strip, at codon 23 past the
    PAM, is exactly this benign class.)

    Returns True if the invariant holds on exit.
    """
    if max_gap is None:
        max_gap = MAX_GAP_DEFAULT
    direction = 1 if cut >= ins else -1
    # `protospacer_span` ALREADY returns [start, end) over the protospacer AND its PAM, so no
    # padding is needed. Scope is exactly Kevin's rule (2026-07-29): "every introduced mutation
    # THAT COULD DISRUPT THE GUIDE AND/OR PAM must be connected back to the insertion by recoded
    # codons" -- mutations outside this span cannot make the locus uncuttable and are left alone.
    p_lo, p_hi = protospacer_span(guide, spec.guide_length)
    p_lo, p_hi = max(p_lo, lo), min(p_hi, hi)

    def cut_side_marks(in_protospacer_only):
        """All introduced mutations on the cut side; optionally only those inside the guide/PAM.

        The DISTINCTION MATTERS and getting it wrong is expensive: the protospacer filter selects
        WHICH mutation must be connected (only one that can make the locus uncuttable), but the
        CONNECTIVITY GAP must be measured over EVERY mutation on the cut side. Measuring the gap
        with the protospacer filter applied ignores the phase-1 edits sitting between the cassette
        and the protospacer, so an already-connected design looks disconnected and is rejected.
        That mistake cost 21% of designability (97.3% -> 76.0%) before it was caught.
        """
        out = []
        for i in range(lo, hi):
            if window_list[i] == window0[i]:
                continue
            if in_protospacer_only and not (p_lo <= i < p_hi):
                continue
            if (direction > 0 and i >= ins) or (direction < 0 and i < ins):
                out.append(i)
        return out

    for _ in range(24):
        marks = cut_side_marks(True)
        if not marks:
            return True
        outer = max(marks) if direction > 0 else min(marks)
        w_lo, w_hi = (ins, outer) if direction > 0 else (outer, ins - 1)
        pts = sorted({m for m in cut_side_marks(False) if w_lo <= m <= w_hi}
                     | {ins if direction > 0 else ins - 1})
        gap = max((b - a - 1 for a, b in zip(pts, pts[1:])), default=0)
        if gap <= max_gap:
            return True
        # fill the codon nearest the insertion that sits in the offending stretch
        cands = sorted((c for c in codons
                        if not c.is_split
                        and all(w_lo <= i <= w_hi for i in c.indices)
                        and all(window_list[i] == window0[i] for i in c.indices)),
                       key=lambda c: min(abs(i - ins) for i in c.indices))
        applied = False
        for codon in cands:
            cur = codon.seq("".join(window_list))
            wt = "".join(window0[i] for i in codon.indices)
            _wm_alts, _wm_used = _wm_bridge_alts(cur, usage)
            for alt in _wm_alts:
                if alt == wt:
                    continue
                for k, idx in enumerate(codon.indices):
                    window_list[idx] = alt[k]
                seg_l = "".join(window_list[lo:ins])
                seg_r = "".join(window_list[ins:hi])
                clean = all(site not in seg and rc(site) not in seg
                            for seg in (seg_l, seg_r) for site in DISALLOWED_IN_DONOR)
                if clean:
                    recoded.append((codon.residue, cur, alt))
                    if _wm_used and diag is not None:
                        diag["wm_missense"] = diag.get("wm_missense", 0) + 1
                    if _wm_used:
                        notes_wm = globals().get("_WM_NOTES_SINK")
                        if notes_wm is not None:
                            notes_wm.append(WM_MISSENSE_NOTE)
                    applied = True
                    break
                for k, idx in enumerate(codon.indices):      # revert
                    window_list[idx] = cur[k]
            if applied:
                break
        if not applied:
            if diag is not None:
                # Classify WHY the stretch cannot be closed. Every base in the offending run is
                # one of: non-coding (no codon exists to recode), inside a SPLIT codon (excluded
                # -- recoding one would need both exon halves), a single-codon amino acid
                # (Met/Trp: no synonymous alternative exists at all), or blocked because every
                # synonymous alternative re-creates a disallowed cloning site.
                worst = None
                for a, b in zip(pts, pts[1:]):
                    if b - a - 1 > max_gap:
                        worst = (a + 1, b - 1)
                        break
                g_lo, g_hi = worst if worst else (w_lo, w_hi)
                noncoding = sum(1 for c2 in codons for i in c2.indices if False)  # placeholder
                idx_in = list(range(g_lo, g_hi + 1))
                cod_in = [c2 for c2 in codons if any(g_lo <= i <= g_hi for i in c2.indices)]
                covered = {i for c2 in cod_in for i in c2.indices}
                noncoding = sum(1 for i in idx_in if i not in covered)
                nosyn = sum(1 for c2 in cod_in
                            if len([x for x in alternatives(c2.seq("".join(window_list)), usage)
                                    if x != "".join(window0[i] for i in c2.indices)]) == 0)
                split = sum(1 for c2 in cod_in if c2.is_split)
                diag.update(stretch_nt=g_hi - g_lo + 1, codons_in_stretch=len(cod_in),
                            noncoding_nt=noncoding, codons_no_synonym=nosyn,
                            codons_split=split,
                            codons_blocked_by_site=max(0, len(cod_in) - nosyn - split))
            return False          # cannot close the stretch synonymously
    if diag is not None:
        diag.update(stretch_nt=-1, codons_in_stretch=-1, noncoding_nt=-1,
                    codons_no_synonym=-1, codons_split=-1, codons_blocked_by_site=-1)
    return False


def _strip_disallowed(window_list, codons, lo, hi, usage, recoded, ins=None, window0=None):
    """
    Synonymously remove BspQI/PaqCI/AscI/NotI sites from the donor arms.

    `ins` splits the region into the two arms. This matters: in the assembled oligo the arms are
    NOT adjacent -- the linkers plus the ~408-nt AsLOV2 stuffer sit between them -- so scanning
    `window_list[lo:hi]` as one block invents **phantom sites straddling the insertion point**
    that cannot exist in the real construct. Recoding those away burns synonymous changes on a
    constraint that isn't there, shortening the homology arms for nothing. That is the same class
    of unnecessary recoding removed from phase 2 (Kevin, 2026-07-25, YOS1).

    Sites created at a real junction (arm|linker, cloning_site|arm, arm|rev_primer) are NOT
    handled here -- they only exist after assembly and are caught by
    `_oligo_sites_are_canonical`, whose remedy is the 3-nt donor slide in the caller.

    `ins=None` preserves the old single-block behaviour for callers that have no insertion point.
    """
    segments = [(lo, hi)] if ins is None else [(lo, ins), (ins, hi)]
    for _ in range(12):
        offending = None
        for seg_lo, seg_hi in segments:
            arm = "".join(window_list[seg_lo:seg_hi])
            for site in DISALLOWED_IN_DONOR:
                for probe in (site, rc(site)):
                    j = arm.find(probe)
                    if j != -1:
                        offending = (seg_lo + j, seg_lo + j + len(probe))
                        break
                if offending:
                    break
            if offending:
                break
        if not offending:
            return True
        a, b = offending
        fixed = False
        for codon in codons:
            if codon.is_split or not any(a <= i < b for i in codon.indices):
                continue
            cur = codon.seq("".join(window_list))
            # WT guard, as in `apply()` (M3/M4). Without it the strip pass can pick the pristine
            # genomic codon as its "alternative" for a codon `_recode_to_destroy` already moved,
            # silently undoing a disruption that was counted toward `min_disruption`.
            _wt = "".join(window0[i] for i in codon.indices) if window0 is not None else None

            # PREFER A SINGLE SNV (Kevin, 2026-07-29): "For RE site stripping ... it is better to do
            # a single SNV where possible. First, it gives us better homology, and second, if for
            # some reason the HR machinery only incorporates one SNV, it is the correct syn codon."
            #
            # The second reason is the load-bearing one and it is not obvious. A codon changed at 2-3
            # positions can be PARTIALLY copied by a conversion tract, and a partial copy need not be
            # synonymous: COF1 codon 109 AGA(R) -> CGT(R) yields CGA(R) if only the first base is
            # taken but AGT(**Ser**) if only the last is. Measured over 15,386 recoded codons:
            # 20.6% carry >1 SNV, and **74.0% of those have a partial copy that is NON-synonymous**.
            # A single-SNV change cannot fail this way -- there is no partial state.
            #
            # Ordering is (n_snv_vs_current, changed base lands INSIDE the offending site, usage
            # rank). One base anywhere in a 7-bp recognition site destroys it, so a 1-SNV option is
            # sufficient whenever its changed base falls in the site; the tie-break keeps the
            # original usage preference among equals.
            _alts = list(alternatives(cur, usage))
            def _rank(alt, _cur=cur, _a=a, _b=b, _ci=codon.indices):
                d = [k for k in range(3) if alt[k] != _cur[k]]
                inside = any(_a <= _ci[k] < _b for k in d)
                return (len(d), 0 if inside else 1, _alts.index(alt))
            for alt in sorted(_alts, key=_rank):
                if _wt is not None and alt == _wt:
                    continue
                for k, idx in enumerate(codon.indices):
                    window_list[idx] = alt[k]
                new_segs = ["".join(window_list[s_lo:s_hi]) for s_lo, s_hi in segments]
                if all(s not in seg and rc(s) not in seg
                       for seg in new_segs for s in DISALLOWED_IN_DONOR):
                    recoded.append((codon.residue, cur, alt))
                    return True
                if all(seg.find(probe) == -1 for seg in new_segs):
                    recoded.append((codon.residue, cur, alt))
                    fixed = True
                    break
                for k, idx in enumerate(codon.indices):      # revert
                    window_list[idx] = cur[k]
            if fixed:
                break
        if not fixed:
            return False
    return False


def design_oligo(model, chrom_seq, residue, spec: OligoSpec, genome=None,
                 max_distance: int = 30, usage: dict = None, flank: int = 400,
                 max_offtargets: int = 1, min_disruption: int = 4,
                 recode_cut_to_insertion: bool = True, allowed_guides: set = None,
                 scores: dict = None, sv_weight: float = 1.0, sv_max: float = None,
                 require_score: bool = False, distance_buckets: tuple = (6, 10, 15, 20),
                 reasons: list = None, prefer_spanning: str = "within_bucket",
                 # Widened from (0,3,-3,6,-6) on 2026-07-27. Ordered by `_slide_anchor_score`,
                 # the extra +/-9/+/-12 codons rescue most `site_EDGE` designs (26.1% -> 4.9%
                 # below 30 nt of effective homology). Cost: the largest left arm grows 72 -> 78,
                 # so read 1 needs 129 in-oligo cycles rather than 123 -- still inside 2x150
                 # (artifact 37). Do NOT widen further without re-checking that budget.
                 slide_range: tuple = (0, 3, -3, 6, -6, 9, -9, 12, -12),
                 protect_boundary_codons: bool = False,
                 credit_pam_separated: bool = False):
    """
    Design one guide-donor oligo inserting the domain AFTER `residue`.

    Returns a DesignResult, or None if no guide within `max_distance` yields a clean donor.
    `genome` (dict of chromosome sequences) enables the exact off-target screen.
    """
    center = model.cds_to_genomic(residue * 3)
    W, c_idx, coords, cds_pos = model.window_with_coords(chrom_seq, center, flank)
    ins = c_idx + 1                                  # junction just after the codon's last nt
    codons = codons_in_window(cds_pos)

    left_len, right_len = spec.split_arms()
    if ins - left_len < 0 or ins + right_len > len(W):
        return None

    candidates = guides_near(W, ins, max_distance, guide_length=spec.guide_length)

    # Consideration 1, ranking half: prefer high predicted efficacy and low SV-formation
    # potential. `guides_near` returns proximity order; when SCORE predictions are supplied we
    # re-rank by quality and use proximity only to break ties. Efficacy is 1 - SCORE.NE (the
    # inverse of predicted UNEDITED likelihood) -- NOT Azimuth, which was trained on NHEJ indel
    # propensity in human cells and does not describe HDR donor editing in yeast.
    if scores is not None:
        from .guides import guide_quality
        ranked = []
        for g in candidates:
            q = guide_quality(g.protospacer, scores, sv_weight=sv_weight)
            if q is None:
                if require_score:
                    continue
                q = float("-inf")          # unscored guides fall to the back, but stay usable
            sv = (scores.get(g.protospacer.upper()) or {}).get("sv")
            if sv_max is not None and sv is not None and sv > sv_max:
                continue
            ranked.append((q, g))
        # PROXIMITY DOMINATES QUALITY. SCORE paper Data S1: the correct-edit rate is ~75% when
        # the edit is within 6 nt of the PAM, 66% at 7-10, 55% at 11-15, and collapses to 8% at
        # 16-20 and 0% beyond 20. A better guide 16 nt away is far worse than an average guide
        # 4 nt away. So candidates are bucketed by distance FIRST and ranked by SCORE quality
        # only WITHIN a bucket. Ranking purely by quality (an earlier version) roughly halved
        # the expected yield of correct edits.
        def _bucket(g):
            # A self-disrupting guide (insertion inside the protospacer) needs no recoding to
            # stop re-cutting and works even where the protospacer sits in non-coding sequence.
            # `prefer_spanning="first"` puts it above every distance bucket; the default
            # "within_bucket" keeps DISTANCE dominant (the one axis we have measured) and uses
            # self-disruption only to break ties inside a bucket -- see artifact 13.
            if prefer_spanning == "first" and spans_insertion(g, spec.guide_length, ins):
                return -1
            d = abs(g.cut_junction - ins)
            for i, edge in enumerate(distance_buckets):
                if d <= edge:
                    return i
            return len(distance_buckets)

        def _span_rank(g):
            """0 for a self-disrupting guide, 1 otherwise -- so spanning sorts first."""
            if prefer_spanning != "within_bucket":
                return 0
            return 0 if spans_insertion(g, spec.guide_length, ins) else 1

        candidates = [g for _, _, _, g in sorted(
            ((_bucket(g), _span_rank(g), -q, g) for q, g in ranked),
            key=lambda t: (t[0], t[1], t[2], abs(t[3].cut_junction - ins)))]

    for guide in candidates:
        # Preferred off-target screen: membership in the precomputed allowed set for the
        # MAGESTIC background strain. This is NOT an S288C screen -- the background carries
        # designed fixes plus WGS-called SNVs/indels (yKR831), and the allowed set was
        # recomputed against it.
        if allowed_guides is not None:
            if guide.protospacer not in allowed_guides:
                continue
            n_off = 0
        elif genome is not None:
            from .guides import count_offtargets
            n_off = count_offtargets(guide.protospacer, genome)
            if n_off > max_offtargets:
                continue
        else:
            n_off = None

        _rej = None
        # Try slides in descending order of the uninterrupted homology they would leave, so a shift
        # that moves a cloning site OUT of an arm beats one that leaves it in to be recoded at the
        # arm's outer end. Ties (including the common case of no sites at all) keep the original
        # 0, 3, -3, 6, -6 preference, so unaffected designs are untouched.
        _order = tuple(slide_range)
        def _score(sh):
            l_s, r_s = left_len + sh, right_len - sh
            if l_s < spec.min_arm or r_s < spec.min_arm:
                return -1
            lo_s, hi_s = ins - l_s, ins + r_s
            if lo_s < 0 or hi_s > len(W):
                return -1
            return _slide_anchor_score(W, lo_s, ins, hi_s)
        # `621`(C): pass 0 is exactly today's behaviour. Pass 1 re-walks the SAME slides with
        # junction recoding enabled, and only runs because pass 0 produced nothing -- which is what
        # makes the fix purely additive (Kevin: "only after every donor slide has failed").
        _jr_shifts = sorted(_order, key=lambda sh: (-_score(sh), _order.index(sh)))
        _jr_plan = [(0, _sh) for _sh in _jr_shifts]
        if JUNCTION_RECODE:
            _jr_plan += [(1, _sh) for _sh in _jr_shifts]
        for _jr_pass, shift in _jr_plan:
            wl = list(W)
            recoded, notes = [], []

            # Sliding the donor moves the LEFT/RIGHT ARM BOUNDARY, keeping the insertion point
            # pinned to the intended residue and the total donor length constant. An earlier
            # version shifted `ins` itself, which silently relocated the insertion by up to two
            # residues -- caught by the validator's C5 arm-vs-native-protein check.
            # The shift must be a multiple of 3 so the left arm stays codon-phased, which is
            # exactly why consideration 3 specifies +/-3 nt.
            left_len_s, right_len_s = left_len + shift, right_len - shift
            if left_len_s < spec.min_arm or right_len_s < spec.min_arm:
                _rej = "arm shorter than min_arm after donor slide"
                continue
            lo, hi = ins - left_len_s, ins + right_len_s
            if lo < 0 or hi > len(wl):
                _rej = "donor runs off the design window"
                continue

            # Pristine copy of the window BEFORE any recoding, so the strip pass can tell the
            # wild-type codon from a codon phase 1/2 deliberately moved (M3/M4).
            wl0 = list(wl)

            # consideration 2: synonymous recoding across the protospacer + PAM
            destroyed, score = _recode_to_destroy(
                wl, codons, guide, lo, hi, usage, recoded, spec,
                min_disruption=min_disruption,
                recode_cut_to_insertion=recode_cut_to_insertion, insertion=ins,
                protect_boundary_codons=protect_boundary_codons,
                credit_pam_separated=credit_pam_separated)
            _spans_ins = spans_insertion(guide, spec.guide_length, ins)
            # 🔴 MUST mirror the extension made inside `_recode_to_destroy`. Crediting a
            # PAM-separated guide there (skipping phase 2) WITHOUT extending this rollback branch
            # leaves `destroyed=False` with phase-1 edits still applied and no rollback, and the
            # design can then die in `_connect_edits_to_insertion` -- so a position that previously
            # designed cleanly falls through to a worse, further candidate. Measured 2026-08-19 on
            # murine core 368: cut -3 / 4 codons became cut -11 / 6 codons, i.e. the "improvement"
            # made that oligo strictly worse. Both sites must move together.
            if credit_pam_separated and not _spans_ins:
                _spans_ins = insertion_separates_pam(guide, spec.guide_length, ins)
            if not destroyed:
                notes.append(f"guide disruption score {score} < {min_disruption}, "
                             "or the site survives recoding")
                _rej = (f"disruption score {score} < {min_disruption}, or guide site "
                        "survives synonymous recoding")
                if _spans_ins:
                    # `117` G1: the insertion splits the protospacer, so 408 nt of AsLOV2 (and the
                    # 38-nt stuffer in the donor) destroy the site structurally. Recoding to
                    # min_disruption is not required here -- and the codons just spent chasing an
                    # unreachable score are pure cost (`99`), so restore the pristine window and
                    # take the free design. near_cognate_recut below remains the real safety gate.
                    wl[:] = list(wl0)
                    del recoded[:]
                    notes.append("spanning guide: cassette disrupts the site structurally; "
                                 "synonymous recoding rolled back")

            # consideration 3: no cloning sites inside the donor
            if not _strip_disallowed(wl, codons, lo, hi, usage, recoded, ins=ins, window0=wl0):
                _rej = "cloning site inside a donor arm could not be removed synonymously"
                continue

            # ORDER MATTERS (Kevin, 2026-07-29): recode first, THEN strip, THEN reconnect. A strip
            # can drop a fresh mutation deep inside the protospacer with wild-type sequence between
            # it and the cassette -- copyable by a tract that never installs the domain (silent
            # false negative). Reconnect, then re-strip because the new codons can themselves
            # create a site.
            _cdiag = {}
            if not _connect_edits_to_insertion(wl, wl0, codons, lo, hi, usage, recoded,
                                               ins, guide.cut_junction, guide, spec,
                                               diag=_cdiag):
                _rej = ("edits inside the protospacer could not be connected to the insertion "
                        "synonymously | " + ";".join(f"{k}={v}" for k, v in sorted(_cdiag.items())))
                continue
            if not _strip_disallowed(wl, codons, lo, hi, usage, recoded, ins=ins, window0=wl0):
                _rej = "cloning site created while reconnecting edits to the insertion"
                continue

            # M4 (adversarial panel, 2026-07-26): `_strip_disallowed` mutates the SAME window after
            # the disruption gate has already been decided, so the score attached to the design was
            # the stale pre-strip value. 60/4,795 shipped rows disagreed with their reported
            # `disruption`, and 3 non-self-disrupting rows shipped BELOW min_disruption. `destroyed`
            # is a gate, not a label -- so re-score and RE-GATE here, after the last mutation.
            _plo, _phi = protospacer_span(guide, spec.guide_length)
            _plo, _phi = max(_plo, 0), min(_phi, len(wl))
            score = disruption_score(
                "".join(wl0[_plo:_phi]), "".join(wl[_plo:_phi]),
                degenerate_offset=degenerate_offset_in_span(guide, spec.guide_length))
            _self_disrupting = spans_insertion(guide, spec.guide_length, ins)
            if credit_pam_separated and not _self_disrupting:
                # must mirror `_recode_to_destroy`, or this re-gate would demand a disruption score
                # the design was deliberately excused from earning.
                _self_disrupting = insertion_separates_pam(guide, spec.guide_length, ins)
            # `117` G2 (the `116` defect): the old line tested the arms JOINED CONTIGUOUSLY --
            # a molecule that is never made. Test the two that are: the repaired allele (arms +
            # 408 nt AsLOV2) and the donor oligo (arms + the 38-nt PaqCI stuffer).
            from .spec import ASLOV2_CDS as _LOV117
            _la117 = "".join(wl[lo:ins]); _ra117 = "".join(wl[ins:hi])
            _rep117 = _la117 + spec.left_linker + _LOV117 + spec.right_linker + _ra117
            _don117 = _la117 + spec.left_linker + spec.stuffer + spec.right_linker + _ra117
            _site_left = (_guide_site_present(_rep117, guide)
                          or _guide_site_present(_don117, guide))
            if _site_left or (not _self_disrupting and score < min_disruption):
                _rej = (f"post-strip re-gate failed: score {score} < {min_disruption}"
                        if not _site_left else "guide site survives after cloning-site strip")
                continue

            left_arm = "".join(wl[lo:ins])
            right_arm = "".join(wl[ins:hi])

            # M5: mismatch-tolerant re-cut check on the REPAIRED ALLELE (arms + linkers + cassette),
            # which also covers the arm|linker junctions where the linker can restore a guide base.
            from .spec import ASLOV2_CDS as _LOV
            _repaired = (left_arm + spec.left_linker + _LOV + spec.right_linker + right_arm)
            if near_cognate_recut(_repaired, guide):
                _rej = "repaired allele retains a near-cognate site (functional NGG, intact seed)"
                continue

            # Consideration 4: the donor should begin on a codon boundary -- but ONLY when the
            # donor start actually lies in coding sequence. Kevin's phrasing is explicit:
            # "this only applies when we are already ~20+ codons into the CDS, otherwise we will
            # of course have part of the 5' UTR sequence there."
            #
            # An earlier version rejected any design whose donor started outside the CDS
            # (cds_pos[lo] is None), which silently killed EVERY insertion within ~22 codons of
            # the start codon or ~21 of the stop, plus everything whose arm crosses an intron.
            # For YOS1 (85 aa, 3 exons) that was 56 of 84 positions -- with 4-7 usable guides
            # available at every one of them. Frame is set by the insertion point sitting on a
            # codon boundary, not by where the arm happens to start.
            cp = cds_pos[lo]
            donor_starts_in_noncoding = cp is None
            if cp is not None and (cp - 1) % 3 != 0:
                _rej = "donor start inside the CDS is not codon-phased"
                continue

            oligo = (spec.fwd_const + guide.protospacer + spec.cloning_site
                     + left_arm + spec.left_linker + spec.stuffer + spec.right_linker
                     + right_arm + spec.rev_primer)

            # `117` G3 -- THE DECISIVE ONE. This was an unconditional reject on the same
            # `destroyed` flag as G1, with no self-disrupting exemption and no fresh rejection
            # reason, so diagnostics blamed G1 and the real gate stayed invisible.
            if not destroyed and not _spans_ins:
                continue

            # Consideration 3, the part that only shows up after assembly: a cloning site can be
            # created ACROSS a junction (arm|linker, linker|stuffer, cloning_site|arm) even when
            # neither piece contains one. Kevin's fix for exactly this case is to slide the donor
            # by 3 nt, which is what the enclosing loop does.
            if _jr_pass == 1 and not _oligo_sites_are_canonical(oligo, spec):
                _jr_fix = _recode_junction_site(wl, codons, lo, hi, ins, spec, guide,
                                                usage, recoded, wl0)
                if _jr_fix is not None:
                    left_arm, right_arm, oligo = _jr_fix
                    notes.append("junction cloning site removed by synonymous recoding (621C)")
            if not _oligo_sites_are_canonical(oligo, spec):
                notes.append(f"cloning site created at a junction; donor slid {shift} nt")
                _rej = "cloning site created ACROSS a junction in the assembled oligo"
                continue

            res = DesignResult(orf=model.orf, residue=residue, oligo=oligo, guide=guide,
                               left_arm=left_arm, right_arm=right_arm, recoded=recoded,
                               donor_shift=shift, notes=notes, offtargets=n_off,
                               left_cds_positions=cds_pos[lo:ins],
                               right_cds_positions=cds_pos[ins:hi])
            res.donor_starts_in_noncoding = donor_starts_in_noncoding
            res._cut_distance = guide.cut_junction - ins
            res.disruption = score
            res.net_recoded = sum(
                1 for c in codons
                if all(lo <= i < hi for i in c.indices)
                and any(wl[i] != wl0[i] for i in c.indices))
            if scores is not None:
                _s = scores.get(guide.protospacer.upper())
                if _s:
                    res.efficacy, res.sv_score = _s["efficacy"], _s["sv"]
            return res

        # Every donor slide for this guide failed. Record why, mirroring the PTC script's
        # `reasons_not_targetable`, so an undesignable position can be diagnosed rather than
        # merely counted.
        if reasons is not None and _rej:
            reasons.append((guide.protospacer, guide.cut_junction - ins, _rej))
    return None
