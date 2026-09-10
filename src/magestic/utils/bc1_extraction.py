"""
Canonical bc1 identification and extraction — ONE definition, every dataset.

Replaces two divergent implementations:
  * ``bc1_from_screen/snakemake/02_parse_bc1_reads.extract_bc1_from_sequence``
    -- exact ``str.find`` prefix (zero error tolerance) + a suffix check that
    degrades to no check at all, returning an unlabelled guess.
  * ``utils/sequence_anchors.extract_bc1_between_anchors`` -- the right shape
    (both-anchor, length-tolerant) but never called for bc1, and carrying the
    same silent fixed-length fallback.

Definition of a bc1 (see planning/pipeline_provenance_rebuild/22_bc1_extraction_review.md):

    [ upstream anchor ][ bc1: N20 ][ layout-specific tail ][ UNIVERSAL anchor ]
      CAGGTCATGCTC                    0-18 nt              CGTGCATAAC

The universal downstream anchor is present in ALL FOUR known layouts at
96.9-99.5% (18_), at a layout-specific offset -- so THE OBSERVED OFFSET IS THE
LAYOUT CALL, per read. Nothing has to be configured per run or per sample.

    offset +18 -> A (1.0/2.5, with the 18-nt spacer)
    offset  +1 -> B / D (no spacer)
    offset +10 -> C (3.0 / yL406 era)

Error tolerance is BOUNDED AND COUNTED, never a similarity ratio: the upstream
anchor is found by pigeonhole k-mer seeding (a pattern with <= k mismatches must
contain at least one of its k+1 non-overlapping seeds intact), then verified by
Hamming distance and accepted iff <= k. ``k`` is an explicit resolved parameter
recorded on every row -- k=1 for PCR amplicons, k=0 for WGS (13_).

Nothing is ever silently dropped, normalized or padded. Every read gets a tier.
"""

from dataclasses import dataclass, asdict
from typing import Optional, List, Tuple, Dict

# --- the canonical anchors -------------------------------------------------
UPSTREAM_ANCHOR = "CAGGTCATGCTC"        # 12-mer; 7 hits genome-wide
UNIVERSAL_DOWNSTREAM = "CGTGCATAAC"     # 10-mer; 8 hits genome-wide (~= the 12-mer)
DOWNSTREAM_FALLBACK = "CGTGCATA"        # 8-mer; 185 hits (26x) -- OFFSET-CONSTRAINED USE ONLY
BC1_LENGTH = 20

#: The anchor-tier vocabulary, in decreasing order of evidence. SINGLE SOURCE
#: OF TRUTH -- every consumer must iterate this rather than an inline literal.
#: A hardcoded 5-tuple in 02_parse_bc1_reads.py made the four tiers added on
#: 2026-09-09 invisible in the log, with printed percentages that did not sum
#: to 100%. Adding a tier must not require finding every list that names them.
ANCHOR_TIERS = (
    "confirmed",                     # full 10-mer, exact, in the window
    "confirmed_mismatch",            # full 10-mer, <=1 substitution, at the offset
    "confirmed_truncated",           # read ends mid-anchor; fitting prefix exact
    "confirmed_truncated_mismatch",  # truncated AND <=1 substitution
    "confirmed_indel",               # <=1 indel, offset-constrained
    "offset_fallback",               # 8-mer at the established offset
    "prefix_only",                   # upstream verified, no downstream confirmation
    "length_mismatch",               # read too short for a full barcode
    "none",                          # no upstream anchor within k on either strand
)

#: Tiers that yield a usable bc1. The other two return bc1=None by construction,
#: so they can never appear as a row in a bc1-keyed table.
EMITTING_TIERS = ANCHOR_TIERS[:6]

#: observed universal-anchor offset -> layout label
LAYOUT_BY_OFFSET: Dict[int, str] = {18: "A_1.0/2.5", 1: "B_no_spacer", 10: "C_3.0"}
#: widest window in which the universal anchor may sit (max offset + its length)
MAX_DOWNSTREAM_OFFSET = 40
#: shortest anchor prefix accepted at an ESTABLISHED offset when the read is
#: truncated mid-anchor. 6 nt at ONE fixed position = 1/4096 by chance, and it
#: sits 20 nt from an already-verified 12-mer upstream anchor.
MIN_TRUNCATED_ANCHOR = 6
#: floor for accepting a truncated anchor prefix that ALSO carries a mismatch.
#: At one fixed offset that is C(8,1)*3/4^8 ~= 3.7e-4 by chance.
MIN_TRUNCATED_ANCHOR_MM = 8

_COMP = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(s: str) -> str:
    return s.translate(_COMP)[::-1]


@dataclass(frozen=True)
class Bc1Call:
    """One read's bc1 verdict. ``bc1`` is None unless a barcode was extracted."""
    bc1: Optional[str]
    anchor_tier: str          # confirmed | offset_fallback | prefix_only | length_mismatch | none
    orientation: Optional[str]  # fwd | rev  (of the READ relative to the anchor pair)
    anchor_offset: Optional[int]
    layout: Optional[str]
    k_used: Optional[int]     # mismatches actually found in the upstream anchor
    bc1_length: Optional[int]
    prefix_observed: Optional[str] = None   # the upstream anchor AS READ
    prefix_start: Optional[int] = None      # its index in the strand actually used
    strand_seq: Optional[str] = None        # that strand (fwd read, or its rev-comp)
    reason: str = ""

    def as_row(self) -> dict:
        return asdict(self)


def _hamming_at_most(a: str, b: str, k: int) -> Optional[int]:
    """Hamming distance between equal-length strings, or None if it exceeds k."""
    d = 0
    for x, y in zip(a, b):
        if x != y:
            d += 1
            if d > k:
                return None
    return d


def _edit_at_most_1(obs: str, exp: str) -> Optional[int]:
    """Levenshtein distance between ``obs`` and ``exp`` if <= 1, else None.

    Unlike Hamming this sees INDELS, which is the point: a 1-nt insertion in the
    tail shifts every downstream base by one, so a Hamming comparison against a
    shifted string looks like ~75% mismatches and is silently discarded. Length
    difference must be 0 or 1; anything larger cannot be within edit distance 1.
    """
    la, lb = len(obs), len(exp)
    if abs(la - lb) > 1:
        return None
    if la == lb:
        d = 0
        for x, y in zip(obs, exp):
            if x != y:
                d += 1
                if d > 1:
                    return None
        return d
    # one is longer by exactly 1: find the first divergence, skip that base
    lo, sh = (obs, exp) if la > lb else (exp, obs)
    i = 0
    while i < len(sh) and lo[i] == sh[i]:
        i += 1
    return 1 if lo[i + 1:] == sh[i:] else None


def _seed_positions(seq: str, pattern: str, k: int) -> List[int]:
    """Candidate start positions via pigeonhole seeding.

    A match with <= k mismatches must contain at least one of the k+1
    non-overlapping seeds intact, so seeding on them cannot miss such a match.
    Exhaustive scan is used when the seeds would be too short to be selective.
    """
    n, m = len(seq), len(pattern)
    if m > n:
        return []
    if k == 0:
        out, start = [], 0
        while True:
            i = seq.find(pattern, start)
            if i < 0:
                return out
            out.append(i)
            start = i + 1
    seed_len = m // (k + 1)
    if seed_len < 5:                       # too short to seed usefully
        return list(range(n - m + 1))
    cands = set()
    for s in range(k + 1):
        off = s * seed_len
        seed = pattern[off:off + seed_len]
        start = 0
        while True:
            i = seq.find(seed, start)
            if i < 0:
                break
            p = i - off
            if 0 <= p <= n - m:
                cands.add(p)
            start = i + 1
    return sorted(cands)


def _find_upstream(seq: str, anchor: str, k: int, protect: int) -> Optional[Tuple[int, int]]:
    """Locate the upstream anchor. Returns (start, mismatches) or None.

    ``protect`` bases at the 3' end of the anchor (those adjacent to the
    barcode) must match EXACTLY: a mismatch there can shift the barcode start,
    which a mismatch mid-anchor cannot (13_).
    """
    i = seq.find(anchor)                   # fast path: the overwhelming majority
    if i >= 0:
        return i, 0
    if k == 0:
        return None
    best = None
    tail_ok = len(anchor) - protect
    for p in _seed_positions(seq, anchor, k):
        cand = seq[p:p + len(anchor)]
        if protect and cand[tail_ok:] != anchor[tail_ok:]:
            continue
        d = _hamming_at_most(cand, anchor, k)
        if d is not None and (best is None or d < best[1]):
            best = (p, d)
            if d == 0:
                break
    return best


def _scan_downstream(tail: str, established_offset: Optional[int]) -> Tuple[Optional[int], str]:
    """Find the universal anchor in the post-barcode tail.

    Returns (offset, tier). Everything beyond the first exact scan is
    OFFSET-CONSTRAINED -- tested at one already-established position, never
    scanned -- so the specificity cost of loosening is a single-position cost,
    not a window-wide one.

    Tiers, in the order they are attempted:
      confirmed            full 10-mer, exact, anywhere in the window
      confirmed_mismatch   full 10-mer with <=1 mismatch AT the established offset
      confirmed_truncated  read ends mid-anchor; longest fitting prefix (>=6 nt)
                           matches exactly at the established offset
      offset_fallback      8-mer at the established offset
      prefix_only          none of the above
    """
    j = tail.find(UNIVERSAL_DOWNSTREAM)
    if 0 <= j <= MAX_DOWNSTREAM_OFFSET:
        return j, "confirmed"
    if established_offset is None:
        return None, "prefix_only"
    o = established_offset
    seg = tail[o:o + len(UNIVERSAL_DOWNSTREAM)]

    # (a) full-length anchor with a single substitution
    if len(seg) == len(UNIVERSAL_DOWNSTREAM):
        if _hamming_at_most(seg, UNIVERSAL_DOWNSTREAM, 1) is not None:
            return o, "confirmed_mismatch"
    # (b) read truncated mid-anchor: longest fitting prefix must match exactly
    elif len(seg) >= MIN_TRUNCATED_ANCHOR:
        exp = UNIVERSAL_DOWNSTREAM[:len(seg)]
        if seg == exp:
            return o, "confirmed_truncated"
        # truncated AND carrying a substitution -- the two degradations co-occur
        # far more often than either alone (measured: it is the single largest
        # residual class on layout A). Needs a longer floor than the exact case.
        if len(seg) >= MIN_TRUNCATED_ANCHOR_MM and _hamming_at_most(seg, exp, 1) is not None:
            return o, "confirmed_truncated_mismatch"
    # (c) the 8-mer fallback
    if tail[o:o + len(DOWNSTREAM_FALLBACK)] == DOWNSTREAM_FALLBACK:
        return o, "offset_fallback"
    # (d) INDEL-aware last resort. A 1-nt insertion or deletion upstream of the
    #     anchor shifts it by one, so the anchor is intact but at o+/-1; an indel
    #     INSIDE the anchor is caught by the edit-distance test. Still
    #     offset-constrained: only o-1, o, o+1 are ever examined.
    for cand in (o, o - 1, o + 1):
        if cand < 0:
            continue
        seg2 = tail[cand:cand + len(UNIVERSAL_DOWNSTREAM)]
        if len(seg2) < MIN_TRUNCATED_ANCHOR_MM:
            continue
        exp2 = UNIVERSAL_DOWNSTREAM[:len(seg2)]
        d = _edit_at_most_1(seg2, exp2)
        if d is not None:
            return cand, "confirmed_indel"
        # DELETION inside the anchor: the observed run is one SHORTER than the
        # anchor, so it must be compared against the FULL anchor, not against a
        # same-length prefix (which is the truncation case, handled at (b)).
        # Require the read to actually extend far enough that this is a deletion
        # and not a truncated read.
        if len(tail) - cand >= len(UNIVERSAL_DOWNSTREAM):
            seg3 = tail[cand:cand + len(UNIVERSAL_DOWNSTREAM) - 1]
            if (len(seg3) >= MIN_TRUNCATED_ANCHOR_MM
                    and _edit_at_most_1(seg3, UNIVERSAL_DOWNSTREAM) is not None):
                return cand, "confirmed_indel"
    return None, "prefix_only"


def extract_bc1(
    seq: str,
    k: int = 1,
    protect: int = 0,
    bc1_length: int = BC1_LENGTH,
    established_offset: Optional[int] = None,
    search_both_strands: bool = True,
) -> Bc1Call:
    """Extract the bc1 from one read. Never returns an unlabelled guess.

    Args:
        k: mismatch budget for the upstream anchor. 1 = PCR amplicons,
           0 = WGS (no amplification error to recover). Recorded on the row.
        protect: 3'-most anchor bases required to match EXACTLY. Default 0.
           This was 1 on the reasoning that "a mismatch adjacent to the barcode
           can shift its start". That is true of an INDEL and false of a
           SUBSTITUTION -- and k only ever admits substitutions, so the guard
           was aimed at a failure mode it could not encounter, while rejecting
           real reads. Measured 2026-09-09: dropping it recovers +0.05 to
           +0.27% of reads, of which 94-100% are DOWNSTREAM-CONFIRMED (i.e. the
           barcode start was independently verified by the universal anchor at
           its established offset), and ZERO barcodes changed for reads both
           settings extracted. The downstream anchor is the real guard against a
           shifted start; ``protect`` was a weaker proxy for it. Set protect=1
           when running WITHOUT downstream confirmation.
        established_offset: this sample's modal universal-anchor offset, if
           known. Enables the 8-mer fallback for reads whose tail is too short
           for the 10-mer (layout D). Omit to disable the fallback entirely.
        search_both_strands: 36.6% of files are reverse-dominant (09_).
    """
    seq = seq.upper()
    strands = ((seq, "fwd"),) if not search_both_strands else (
        (seq, "fwd"), (reverse_complement(seq), "rev"))

    best_miss = None
    for s, orient in strands:
        hit = _find_upstream(s, UPSTREAM_ANCHOR, k, protect)
        if hit is None:
            continue
        start, mism = hit
        b0 = start + len(UPSTREAM_ANCHOR)
        tail = s[b0 + bc1_length:]
        offset, tier = _scan_downstream(tail, established_offset)

        if offset is not None:
            bc1 = s[b0:b0 + bc1_length]
            obs = s[start:start + len(UPSTREAM_ANCHOR)]
            if len(bc1) != bc1_length:
                return Bc1Call(None, "length_mismatch", orient, offset,
                               LAYOUT_BY_OFFSET.get(offset), mism, len(bc1), obs, start, s,
                               f"read ends {bc1_length - len(bc1)} nt short of a full barcode")
            return Bc1Call(bc1, tier, orient, offset,
                           LAYOUT_BY_OFFSET.get(offset), mism, bc1_length, obs, start, s)

        # prefix found, no downstream confirmation -- keep it, but SAY SO
        bc1 = s[b0:b0 + bc1_length]
        obs = s[start:start + len(UPSTREAM_ANCHOR)]
        if len(bc1) != bc1_length:
            best_miss = Bc1Call(None, "length_mismatch", orient, None, None, mism,
                                len(bc1), obs, start, s, "read too short for a full barcode")
            continue
        best_miss = Bc1Call(bc1, "prefix_only", orient, None, None, mism, bc1_length, obs,
                            start, s, "upstream anchor verified; no downstream anchor in window")

    if best_miss is not None:
        return best_miss
    return Bc1Call(None, "none", None, None, None, None, None, None, None, None,
                   "upstream anchor not found within k on either strand")


# =============================================================================
# Per-barcode prefix catalog
# =============================================================================
#
# Kevin, 2026-09-09: "a catalog of barcode prefixes identified, so that as we are
# matching barcodes, we can distinguish between seq errors and true prefix
# differences for each barcode."
#
# The separating principle: a TRUE structural difference is REPRODUCIBLE across
# independent reads of the same barcode. A sequencing error is not. A barcode
# seen 500x with the canonical anchor and twice with a variant has two errors; a
# barcode seen 290 of 300 times with that same variant has a real difference --
# a subpool cloned on a different backbone, say. The canonical anchor is
# therefore a PRIOR, not a rule, and each barcode gets its own empirical answer.

from collections import defaultdict, Counter as _Counter


class PrefixCatalog:
    """Observed upstream-anchor context per bc1, and the error-vs-real verdict.

    Two passes: ``observe`` every confirmed read, then ``classify`` the doubtful
    ones against what that barcode actually looks like.
    """

    #: a context must clear BOTH to count as a real difference rather than noise
    def __init__(self, min_support: int = 3, min_fraction: float = 0.10):
        self.min_support = min_support
        self.min_fraction = min_fraction
        self._ctx: Dict[str, "_Counter[Tuple[str, Optional[int]]]"] = defaultdict(_Counter)

    def observe(self, bc1: str, prefix_observed: str,
                anchor_offset: Optional[int], count: int = 1) -> None:
        if bc1 and prefix_observed:
            self._ctx[bc1][(prefix_observed, anchor_offset)] += count

    def observe_call(self, call: "Bc1Call", count: int = 1) -> None:
        if call.bc1 and call.prefix_observed:
            self.observe(call.bc1, call.prefix_observed, call.anchor_offset, count)

    def total(self, bc1: str) -> int:
        return sum(self._ctx.get(bc1, {}).values())

    def consensus(self, bc1: str) -> Optional[Tuple[str, Optional[int]]]:
        c = self._ctx.get(bc1)
        return c.most_common(1)[0][0] if c else None

    def classify(self, bc1: str, prefix_observed: str,
                 anchor_offset: Optional[int] = None) -> Dict[str, object]:
        """Is this context a real property of this barcode, or a read error?

        Returns a dict with ``verdict``, ``support``, ``total``, ``fraction``
        and the barcode's ``consensus``. Verdicts:
          novel_barcode      -- bc1 absent from the catalog; no opinion possible
          canonical          -- context IS this barcode's consensus
          variant_confirmed  -- non-consensus, but reproducible for THIS barcode
                                (>= min_support reads AND >= min_fraction)
          error_likely       -- non-consensus and rare for this barcode
        """
        c = self._ctx.get(bc1)
        if not c:
            return {"verdict": "novel_barcode", "support": 0, "total": 0,
                    "fraction": 0.0, "consensus": None}
        key = (prefix_observed, anchor_offset)
        tot = sum(c.values())
        sup = c.get(key, 0)
        cons = c.most_common(1)[0][0]
        frac = sup / tot if tot else 0.0
        if tot < self.MIN_COVERAGE_FOR_CONSENSUS and key != cons:
            # not enough reads of THIS barcode to tell a variant from a misread
            return {"verdict": "insufficient_coverage", "support": sup, "total": tot,
                    "fraction": frac, "consensus": cons}
        if key == cons:
            verdict = "canonical"
        elif sup >= self.min_support and frac >= self.min_fraction:
            verdict = "variant_confirmed"
        else:
            verdict = "error_likely"
        return {"verdict": verdict, "support": sup, "total": tot,
                "fraction": frac, "consensus": cons}

    #: a consensus drawn from fewer reads than this is not a consensus. Measured
    #: 2026-09-09: at 1.3 reads/barcode EVERY "non-canonical consensus" was a
    #: low-coverage artifact (2/2, 1/2, 4/4 reads) or an N base-call failure.
    #: Two reads cannot separate a real variant from one misread base.
    MIN_COVERAGE_FOR_CONSENSUS = 10

    def variant_barcodes(self, min_total: Optional[int] = None
                         ) -> List[Tuple[str, Tuple[str, Optional[int]], int, int]]:
        """Barcodes whose CONSENSUS context is genuinely not the canonical anchor.

        Only barcodes with >= ``min_total`` observations are considered: below
        that the consensus is noise, not signal. Contexts containing N are
        excluded outright -- an uncalled base is an instrument failure, not a
        construct difference.
        """
        floor = self.MIN_COVERAGE_FOR_CONSENSUS if min_total is None else min_total
        out = []
        for bc1, c in self._ctx.items():
            tot = sum(c.values())
            if tot < floor:
                continue
            cons = c.most_common(1)[0][0]
            if cons[0] != UPSTREAM_ANCHOR and "N" not in cons[0]:
                out.append((bc1, cons, c[cons], tot))
        return sorted(out, key=lambda r: -r[3])

    def coverage_report(self) -> Dict[str, object]:
        """Is this catalog deep enough to be believed?

        Per-barcode coverage is the limiting resource: the error-vs-real call is
        only possible for barcodes observed enough times for a minority context
        to be distinguishable from a single misread base. Build the catalog at
        SCREEN scale, not per file.
        """
        tots = sorted(sum(c.values()) for c in self._ctx.values())
        if not tots:
            return {"barcodes": 0, "usable": 0, "usable_fraction": 0.0}
        n_ok = sum(1 for t in tots if t >= self.MIN_COVERAGE_FOR_CONSENSUS)
        return {"barcodes": len(tots),
                "observations": sum(tots),
                "mean_coverage": sum(tots) / len(tots),
                "median_coverage": tots[len(tots) // 2],
                "usable": n_ok,
                "usable_fraction": n_ok / len(tots),
                "min_coverage_required": self.MIN_COVERAGE_FOR_CONSENSUS}

    def summary(self) -> Dict[str, int]:
        n_multi = sum(1 for c in self._ctx.values() if len(c) > 1)
        return {"barcodes": len(self._ctx),
                "barcodes_with_multiple_contexts": n_multi,
                "barcodes_with_noncanonical_consensus": len(self.variant_barcodes()),
                "observations": sum(sum(c.values()) for c in self._ctx.values())}
