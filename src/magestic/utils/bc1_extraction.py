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
    return None, "prefix_only"


def extract_bc1(
    seq: str,
    k: int = 1,
    protect: int = 1,
    bc1_length: int = BC1_LENGTH,
    established_offset: Optional[int] = None,
    search_both_strands: bool = True,
) -> Bc1Call:
    """Extract the bc1 from one read. Never returns an unlabelled guess.

    Args:
        k: mismatch budget for the upstream anchor. 1 = PCR amplicons,
           0 = WGS (no amplification error to recover). Recorded on the row.
        protect: 3'-most anchor bases required to match exactly.
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
            if len(bc1) != bc1_length:
                return Bc1Call(None, "length_mismatch", orient, offset,
                               LAYOUT_BY_OFFSET.get(offset), mism, len(bc1),
                               f"read ends {bc1_length - len(bc1)} nt short of a full barcode")
            return Bc1Call(bc1, tier, orient, offset,
                           LAYOUT_BY_OFFSET.get(offset), mism, bc1_length)

        # prefix found, no downstream confirmation -- keep it, but SAY SO
        bc1 = s[b0:b0 + bc1_length]
        if len(bc1) != bc1_length:
            best_miss = Bc1Call(None, "length_mismatch", orient, None, None, mism,
                                len(bc1), "read too short for a full barcode")
            continue
        best_miss = Bc1Call(bc1, "prefix_only", orient, None, None, mism, bc1_length,
                            "upstream anchor verified; no downstream anchor in window")

    if best_miss is not None:
        return best_miss
    return Bc1Call(None, "none", None, None, None, None, None,
                   "upstream anchor not found within k on either strand")
