"""
Unified sliding window matching for donor sequences.

This module implements consistent donor matching for both gdb and bdb data:
1. Build lookup tables by sliding a window across DESIGNED oligos
2. Match observed donors by sliding a window across the OBSERVED sequence
3. This handles truncations, extensions, and internal mutations consistently

The key insight is that we slide the window on BOTH sides:
- When building the lookup: slide across designed donor
- When matching: slide across observed donor

This fixes the previous discrepancy where gdb used sliding window but bdb
used fixed extraction.

Uses the pre-parsed harmonized QTL oligo table for annotations, which
already contains extracted guide and donor sequences, eliminating the
need for runtime extraction from oligo_seq.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Set, Tuple
from pathlib import Path
import pandas as pd
import logging
import re
import sys

from ..config import (
    GUIDE_START, GUIDE_LENGTH, DONOR_START, SPCAS9_CLONING_SITE,
    MIDDLE_DONOR_LENGTH,
)

# Import harmonized oligo annotations
from magestic.utils.oligo_annotations import (
    load_oligo_annotations,
    get_annotation_columns_for_bc1_donor_bc0,
    HARMONIZED_OLIGO_TABLE,
)

logger = logging.getLogger(__name__)


# Default window size for middle donor matching
DEFAULT_WINDOW_SIZE = 60

# Anchors used to pre-trim observed bc1_donor_bc0 reads down to the variable
# donor region. The observed donor extracted at the parse step carries a
# constant 5' flank (SaCas9/SceI fragment + scaffold leader) terminating in a
# cloning-enzyme site; everything AFTER that anchor is the designed donor
# (~129bp / ~117bp). Trimming up to and including this anchor lets the
# perfect-donor lookup fire.
#
# Two anchor systems are in use across the MAGESTIC plasmid generations:
#   GTTTGAAGAGC -- SpCas9 / SpG / SaCas9 (MAGESTIC 1.0 SpCas9, 2.5 SpCas9, 3.0 SpG)
#   TTTCGAAGAGC -- LbCas12a / impLbCas12a (MAGESTIC 1.0/2.5 LbCas12a libraries:
#                  V450, V454, V455; legacy 20210822 / 20220607 / 20211206 screens)
# Without the second anchor every V450/V454/V455 read returns an untrimmed
# 238-bp donor that the perfect lookup cannot hit -- they fall through to the
# 60-bp sliding-window fallback instead. Forward-ported from the standalone
# `04_map_to_oligos.py` script (HANDOFF_20260527).
DONOR_VARIABLE_REGION_ANCHOR = "GTTTGAAGAGC"
DONOR_VARIABLE_REGION_ANCHORS = [
    "GTTTGAAGAGC",  # SpCas9, SpG, SaCas9
    "TTTCGAAGAGC",  # LbCas12a, impLbCas12a
]

# Canonical position of the anchor in the parsed donor field for current
# MAGESTIC 3.0 amplicons. 88% of well-formed reads place it at exactly 98;
# Q-drop / single-base sequencing errors in the 11-bp anchor account for
# another ~10.5% (12% of all reads lose the exact-match anchor, of which
# ~89% rescue with Hamming-1 at the canonical position).
DEFAULT_ANCHOR_POS = 98
DEFAULT_ANCHOR_POS_SLACK = 10
DEFAULT_ANCHOR_MAX_MISMATCH = 1


def _hamming(a: str, b: str) -> int:
    """Count mismatches between two equal-length strings. N counts as mismatch."""
    return sum(x != y for x, y in zip(a, b))


def _find_anchor_tolerant(
    donor: str,
    anchor: str,
    expected_pos: int,
    slack: int,
    max_mismatch: int,
) -> int:
    """
    Return position of the best-matching anchor window in [expected_pos-slack,
    expected_pos+slack], allowing up to max_mismatch substitutions. -1 if no
    window meets the threshold. Ties broken by closeness to expected_pos.
    """
    L = len(anchor)
    if len(donor) < L:
        return -1
    start = max(0, expected_pos - slack)
    end = min(len(donor) - L + 1, expected_pos + slack + 1)
    best_pos = -1
    best_mm = max_mismatch + 1
    best_dist = None
    for i in range(start, end):
        mm = _hamming(donor[i:i + L], anchor)
        if mm > max_mismatch:
            continue
        dist = abs(i - expected_pos)
        if mm < best_mm or (mm == best_mm and (best_dist is None or dist < best_dist)):
            best_mm = mm
            best_pos = i
            best_dist = dist
            if mm == 0 and dist == 0:
                return i
    return best_pos


def trim_to_variable_region(
    donor: Optional[str],
    anchors: Sequence[str] = DONOR_VARIABLE_REGION_ANCHORS,
    expected_pos: int = DEFAULT_ANCHOR_POS,
    pos_slack: int = DEFAULT_ANCHOR_POS_SLACK,
    max_mismatch: int = DEFAULT_ANCHOR_MAX_MISMATCH,
    anchor: Optional[str] = None,
) -> Optional[str]:
    """
    Trim an observed donor down to the variable region following the anchor.

    Tries each anchor in `anchors` in order. The SpCas9/SpG anchor
    (GTTTGAAGAGC) hits in ~88% of MAGESTIC 3.0 reads; the LbCas12a anchor
    (TTTCGAAGAGC) is needed for V450/V454/V455 reads in the legacy MAGESTIC
    1.0 / 2.5 screens. Anchor-conditional + idempotent:

    1. None/NaN: returned unchanged.
    2. For each anchor in `anchors`, in order:
       a. Exact match anywhere in the donor: trim at first occurrence (fast
          path, ~88% of yL406-like reads on the first anchor).
       b. No exact match: search a Hamming-up-to-`max_mismatch` window in
          [expected_pos - pos_slack, expected_pos + pos_slack]. This rescues
          the ~12% of reads where a single-base sequencer Q-drop destroyed
          the exact match. Ties broken by closeness to `expected_pos`.
    3. No anchor in `anchors` matched: return donor unchanged (already-
       trimmed donors, clean gdb donors, reads with corrupted flank all
       flow to sliding-window fallback).

    Args:
        donor: Observed donor sequence (may be None).
        anchors: Anchors to try, in order. Default is the SpCas9/SpG +
            LbCas12a pair covering all MAGESTIC plasmid generations.
        expected_pos: Canonical 0-indexed start of the anchor in the donor
            field (98 for MAGESTIC 3.0 bc1_donor_bc0 amplicons).
        pos_slack: How far either side of expected_pos to scan for the
            tolerant match.
        max_mismatch: Hamming budget for the tolerant search. Set to 0 to
            disable the tolerant path entirely and reproduce exact-only
            behavior.
        anchor: Back-compat singular parameter. If provided, overrides
            `anchors` with a single-element list. Preserves the previous API
            for callers passing one anchor by keyword.

    Returns:
        The trimmed donor, or the original donor if no anchor is found.
    """
    if not donor:
        return donor
    if anchor is not None:
        anchors = (anchor,)
    # 1) Exact-match pass across all anchors (fastest).
    for a in anchors:
        idx = donor.find(a)
        if idx >= 0:
            return donor[idx + len(a):]
    if max_mismatch <= 0:
        return donor
    # 2) Tolerant pass across all anchors; first hit wins.
    for a in anchors:
        pos = _find_anchor_tolerant(donor, a, expected_pos, pos_slack, max_mismatch)
        if pos >= 0:
            return donor[pos + len(a):]
    return donor


GUIDE_NO_GUIDE = "no_guide"
GUIDE_NOT_DESIGNED = "guide_not_designed"
GUIDE_MATCH = "match"
GUIDE_CONFLICT = "conflict"
GUIDE_STATUSES = (GUIDE_NO_GUIDE, GUIDE_NOT_DESIGNED, GUIDE_MATCH, GUIDE_CONFLICT)


@dataclass
class OligoMatch:
    """Result of matching a sequence to designed oligos."""
    oligo_name: Optional[str] = None

    # Match quality indicators
    perfect_donor: bool = False
    perfect_middle_region: bool = False

    # `guide_match` is a BOOLEAN and therefore conflates states that mean opposite things: False was
    # both "this row has no linked guide at all" and "the linked guide disagrees with the donor".
    # Measured on the 2026-09-10 rebuild, False covered 711,464 rows of the first kind and 223,954 of
    # the second. `guide_status` separates them (planning 30_ step 6):
    #   no_guide            no guide on this row
    #   guide_not_designed  a guide is present but is not any designed guide in scope
    #   match               the guide's designs intersect the donor's candidates; it narrowed the match
    #   conflict            the guide IS designed, but names NO design the donor allows
    # `conflict` is the state that used to vanish: the matcher narrows by guide only `if overlap:`, so
    # when the guide's designs did not intersect the donor's it kept the donor answer and said nothing.
    # That is how ambiguous rows ended up carrying a guide matching none of their own candidates.
    guide_match: bool = False
    guide_status: str = GUIDE_NO_GUIDE
    num_guide_oligos: int = 0

    # NN 250 / D13. Do the candidate oligos agree on the PAM TUPLE (PAM, PAM_coord, PAM_strand)?
    # None = unknown (the oligo table carried no PAM columns), which is a THIRD state and must
    # never be read as agreement. False does NOT make the match ambiguous at the EDIT level -- see
    # _resolve_candidates and 250_bc1_pipeline_fixes.md §9 -- it says the PAM ATTRIBUTION is
    # undetermined, which is what the per-well PAM-class tables depend on.
    candidates_agree_on_pam: Optional[bool] = None
    num_candidate_pams: int = 0

    # Confidence
    #
    # `is_unambiguous` is VARIANT-level, not oligo_name-level. Many designed
    # oligos share an identical donor (the HDR repair template) but differ in
    # which guide cuts -- they install the SAME edit. Collapsing candidate
    # oligos to their `donor_changes` (the installed variant) shows a single
    # variant in ~100% of perfect-donor matches even though the oligo count is
    # 4-6. A match is therefore "unambiguous" when all candidate oligos agree
    # on the installed variant, regardless of how many redundant guides exist.
    is_unambiguous: bool = False
    num_candidate_oligos: int = 0     # distinct candidate oligo_names
    num_candidate_variants: int = 0   # distinct installed variants among candidates
    variant: Optional[str] = None     # donor_changes (installed edit) when unambiguous

    # Details
    match_method: str = "none"  # "perfect", "middle_region", "guide_only", "none"
    candidate_oligos: List[str] = field(default_factory=list)


@dataclass
class OligoLookupTables:
    """
    Precomputed lookup tables for fast oligo matching.

    Contains multiple lookup strategies:
    - perfect_donor: exact donor sequence match
    - sliding_window: 60bp windows from designed donors
    - guide: guide sequence to oligos
    - guide_donor: combined guide+donor for perfect match
    """
    # Exact matches
    perfect_donor_to_oligos: Dict[str, List[str]] = field(default_factory=dict)
    guide_to_oligos: Dict[str, List[str]] = field(default_factory=dict)
    guide_donor_to_oligo: Dict[str, str] = field(default_factory=dict)

    # Sliding window matches (key = 60bp window)
    window_to_oligos: Dict[str, List[str]] = field(default_factory=dict)

    # oligo_name -> installed variant (donor_changes). Lets the matcher collapse
    # candidate oligos that share a donor to the single edit they all install,
    # so variant-level unambiguity can be assessed in every match branch.
    oligo_to_variant: Dict[str, str] = field(default_factory=dict)

    # oligo_name -> (PAM, PAM_coord, PAM_strand). Empty when the oligo table carries no PAM
    # columns. NN 250 / D13: lets a match report whether its candidates agree on the PAM TUPLE,
    # which is what the per-well PAM-class consumers key on, without touching edit identity.
    oligo_to_pam: Dict[str, tuple] = field(default_factory=dict)

    # Metadata
    window_size: int = DEFAULT_WINDOW_SIZE
    # NN 250 / D10. The anchor-trim thresholds ride on the lookup rather than being read from
    # module defaults inside match_sequence, so a CLI override reaches EVERY match path -- the
    # vectorized fast path, the per-row fallback, and each per-subpool lookup the SPS-scoped
    # closure builds. Threading them through call arguments instead is how a flag ends up working
    # unscoped and being silently ignored once --sps-mapping is given.
    anchor_pos_slack: int = DEFAULT_ANCHOR_POS_SLACK
    anchor_max_mismatch: int = DEFAULT_ANCHOR_MAX_MISMATCH
    num_oligos: int = 0
    #: NN 250 / C20: oligos sharing one guide+donor key, which guide_donor_to_oligo collapses.
    guide_donor_key_collisions: int = 0

    def get_stats(self) -> Dict:
        """Get statistics about the lookup tables."""
        return {
            "num_oligos": self.num_oligos,
            "num_unique_donors": len(self.perfect_donor_to_oligos),
            "num_unique_guides": len(self.guide_to_oligos),
            "num_unique_windows": len(self.window_to_oligos),
            "window_size": self.window_size,
            "guide_donor_key_collisions": self.guide_donor_key_collisions,
            "num_oligos_with_pam": len(self.oligo_to_pam),
        }


def load_oligo_design(filepath: Path) -> pd.DataFrame:
    """
    Load and parse an oligo design file.

    The oligo file has columns: oligo_name, oligo_seq
    The oligo_seq is a 200mer with structure:
    - positions 0-19: 5' adapter (lowercase)
    - positions 20-39: guide (uppercase)
    - positions 40-50: scaffold/cloning site (lowercase)
    - positions 51-179: donor (uppercase)
    - positions 180-199: 3' adapter (lowercase)

    DEPRECATED: Use load_oligos_from_harmonized_table() instead, which
    provides pre-parsed guide and donor sequences with full annotations.

    Args:
        filepath: Path to the oligo design TSV file

    Returns:
        DataFrame with extracted guide and donor columns
    """
    logger.info(f"Loading oligo design from {filepath}")
    df = pd.read_csv(filepath, sep='\t')

    # Check if guide and donor columns already exist (pre-parsed)
    if 'guide' in df.columns and 'donor' in df.columns:
        logger.info("  Using pre-parsed guide and donor columns")
    else:
        # Extract guide and donor using slice positions
        oligo_col = 'oligo_seq' if 'oligo_seq' in df.columns else 'oligo_sequence'
        if oligo_col not in df.columns:
            raise ValueError(f"No oligo sequence column found. Expected 'oligo_seq' or 'oligo_sequence'")

        df['guide'] = df[oligo_col].str[GUIDE_START:GUIDE_START + GUIDE_LENGTH]
        df['donor'] = df[oligo_col].str[DONOR_START:].str.extract(r'([A-Z]+)')[0]

    # Clean up - ensure uppercase
    df['guide'] = df['guide'].str.upper()
    df['donor'] = df['donor'].str.upper()

    logger.info(f"  Loaded {len(df)} oligos")
    logger.info(f"  Guide lengths: {df['guide'].str.len().describe()}")
    logger.info(f"  Donor lengths: {df['donor'].str.len().describe()}")

    return df


def load_oligos_from_harmonized_table(
    v_libraries: Optional[List[str]] = None,
    nuclease_types: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Load oligo designs from the harmonized QTL oligo table.

    This is the PREFERRED method for loading oligo designs. The harmonized
    table already has pre-parsed guide and donor sequences plus full
    annotations.

    Args:
        v_libraries: Filter by V library (V416, V419, V520, V521). None = all.
        nuclease_types: Filter by nuclease/PAM type (SpG_NGNG, SpCas9_NGG, etc.).
            None = all.

    Returns:
        DataFrame with oligo_name, guide, donor, and annotation columns
    """
    logger.info("Loading oligos from harmonized table...")

    # Load all columns needed for bc1_donor_bc0 pipeline
    columns = get_annotation_columns_for_bc1_donor_bc0()
    df = load_oligo_annotations(columns=columns, add_is_control=True)

    # Filter by V library if specified
    if v_libraries:
        df = df[df['guide_donor_bc0_plasmid_library'].isin(v_libraries)]
        logger.info(f"  Filtered to V libraries {v_libraries}: {len(df)} oligos")

    # Filter by nuclease type if specified
    if nuclease_types:
        df = df[df['guide_nuclease_PAM'].isin(nuclease_types)]
        logger.info(f"  Filtered to nuclease types {nuclease_types}: {len(df)} oligos")

    # Ensure uppercase for matching
    df['guide'] = df['guide'].str.upper()
    df['donor'] = df['donor'].str.upper()

    logger.info(f"  Loaded {len(df)} oligos from harmonized table")
    logger.info(f"  Guide lengths: {df['guide'].str.len().describe()}")
    logger.info(f"  Donor lengths: {df['donor'].str.len().describe()}")

    return df


# Column used as the "installed variant" identity for collapsing guide-redundant
# candidate oligos. This MUST be the canonical per-edit identifier:
#   individual_variants = "<pos>_<ref>_<alt>_<SNP|INDEL>" (e.g. "25340_A_C_SNP")
#   26,126 distinct values, 87.4% populated, exactly 1.00 per donor, and it
#   correctly MERGES the ~2.67 distinct donor encodings (synonymous / strand /
#   shift variants) that install the same edit.
#
# NOT donor_changes: that column (242 distinct, 28.3% populated, mostly "[]") is
# donor-CONSTRUCTION metadata ("restriction site GCTCTTC ... removed with
# synonymous codons"), shared by up to 37,400 distinct donors -- collapsing by
# it conflates thousands of different edits and badly inflates is_unambiguous.
# NOT the raw donor sequence: 69,799 distinct, it over-splits each edit into
# ~2.67 encodings, deflating is_unambiguous. donor-seq is only the per-row
# FALLBACK below for the ~12.6% (controls / unannotated) with no variant id.
VARIANT_COLUMN = "individual_variants"

#: An `individual_variants` value that names a genomic position: `<pos>_<ref>_<alt>_SNP`.
POSITIONAL_VARIANT_RE = re.compile(r"^\d+_[ACGT]+_[ACGT]+_(SNP|INDEL|DEL|INS)$", re.IGNORECASE)


def _nonblank(x) -> bool:
    if x is None:
        return False
    if isinstance(x, float) and x != x:          # NaN
        return False
    s = str(x).strip()
    return s != "" and s.upper() not in ("NA", "NAN", "NONE")


def _norm_coord(x) -> str:
    """Normalise a coordinate to a STABLE string: an integral value becomes an integer string.

    🔴 Found by the NN 250 real-data run, not by reading the code. The harmonized table stores
    ``variant_coord`` as ``706445.0``, so the variant id came out ``chrII:706445.0:deletion``.
    That is stable only while pandas infers the same dtype for the column; read as int in one run
    and float in another, one edit would get TWO different ids -- and a join on them returns
    nothing, which reads as a biological null.

    A genuinely fractional coordinate is left as-is rather than truncated: truncating would merge
    two distinct coordinates, which is the opposite of the defect being fixed.
    """
    s = str(x).strip()
    try:
        f = float(s)
    except (TypeError, ValueError):
        return s
    return str(int(f)) if f == int(f) else s


def qualified_variant_id(chrom, individual_variants, variant_coord=None, *,
                         allow_unqualified=False):
    """A CHROMOSOME-QUALIFIED, POSITIONAL variant id, or ``None`` when one cannot be built.

    🔴 NN 250 / D15. `individual_variants` was used directly as the variant identity that decides
    `is_unambiguous`, and MEASURED on the real shipped harmonized table (146,250 rows) it is
    neither chromosome-qualified nor reliably positional:

      * 101,681 positional rows (`<pos>_<ref>_<alt>_SNP`) over 22,252 distinct ids, **50 of which
        appear on MORE THAN ONE CHROMOSOME** -- so two genuinely different edits merge into one
        variant and an ambiguous match reports `is_unambiguous=True`;
      * 26,101 rows carry a **CLASS LABEL**, not an identity, across 3,874 labels -- `PTC` alone
        covers 3,054 rows. Treating `PTC` as a variant collapses 3,054 different premature-stop
        designs into ONE variant, so a match among any of them looks unambiguous;
      * 18,468 blank.

    🔴 **CORRECTED 2026-10-02 (NN 253 §S2; Kevin's D15 ruling: "Remove the key").** An earlier
    version of this function ALSO keyed on `variant_coord`, returning ``f"{c}:{vc}:{iv}"``. That was
    backwards and it corrupted production output. All 327 groups it "corrected" have IDENTICAL
    `individual_variants` and a **byte-identical donor**, with ``abs(variant_coord - PAM_coord) <= 3``
    on 1,586/1,586 rows and `codon_variant_type = complete_CDS_deletion` on 1,568: one donor and
    ~5 redundant guides whose coordinates span up to 2,383 bp. Keying on the coordinate therefore
    SPLIT one real edit into several ids and manufactured false ambiguity -- **1,492 production rows,
    708 of them `perfect_ambiguous`, downgraded on a single plate.** Removing the branch measured
    327 -> 0 false splits and recoveries 104 -> 151, left the positional id space byte-identical
    (108,880 rows / 24,241 distinct ids), kept `test_known_answer_roundtrip` green, and broke only
    the 10 tests that asserted the removed branch.

    `variant_coord` is still ACCEPTED and is now deliberately IGNORED, so existing callers keep
    working. A class label is routed to the donor-sequence fallback that :func:`build_lookup_tables`
    already has for unannotated rows; a coordinate never turns a label into an identity.

    THE TWO HALVES ARE DECOUPLED, deliberately:

      * **A class label is never an identity.** That holds whether or not a chromosome is
        available, so it is enforced unconditionally.
      * **Chromosome qualification** needs a `chrom` column. With ``allow_unqualified=True`` a
        POSITIONAL id with no chromosome is returned unqualified -- the pre-NN-250 behaviour --
        because demanding a chromosome nobody supplied would send every row to the donor fallback
        and destroy variant-level unambiguity outright. (It did: 17 tests, including the whole NN
        240 roundtrip, whose fixtures carry valid positional ids and no `chrom`.) The caller warns
        when it takes that path, since cross-chromosome collisions stay possible there.
    """
    iv = str(individual_variants).strip() if _nonblank(individual_variants) else None
    if iv is None:
        return None
    if not _nonblank(chrom):
        # No chromosome. A positional id is still an identity, just an unqualified one.
        if allow_unqualified and POSITIONAL_VARIANT_RE.match(iv):
            return iv
        return None
    c = str(chrom).strip()
    if POSITIONAL_VARIANT_RE.match(iv):
        return f"{c}:{iv}"
    # A class label names nothing. It goes to the donor-sequence fallback rather than being
    # given a coordinate, because a coordinate does not make a label an identity -- see the
    # D15 note above for why coordinate-keying was removed.
    return None


def build_lookup_tables(
    oligo_df: pd.DataFrame,
    window_size: int = DEFAULT_WINDOW_SIZE,
    variant_col: str = VARIANT_COLUMN,
    anchor_pos_slack: int = DEFAULT_ANCHOR_POS_SLACK,
    anchor_max_mismatch: int = DEFAULT_ANCHOR_MAX_MISMATCH,
) -> OligoLookupTables:
    """
    Build lookup tables from oligo design DataFrame.

    This function creates:
    1. Perfect donor lookup (exact match)
    2. Sliding window lookup (60bp windows across designed donor)
    3. Guide lookup (for guide-based matching)
    4. Combined guide+donor lookup (for perfect match)

    PERFORMANCE: Uses vectorized pandas operations (~10-50x faster than iterrows).

    Args:
        oligo_df: DataFrame with columns: oligo_name, guide, donor
        window_size: Size of sliding window (default 60bp)

    Returns:
        OligoLookupTables with all lookup dictionaries
    """
    logger.info(f"Building lookup tables (window_size={window_size})...")

    tables = OligoLookupTables(window_size=window_size,
                               anchor_pos_slack=anchor_pos_slack,
                               anchor_max_mismatch=anchor_max_mismatch)

    # Filter valid rows
    valid_df = oligo_df.dropna(subset=['guide', 'donor']).copy()
    tables.num_oligos = len(valid_df)

    if tables.num_oligos == 0:
        logger.warning("No valid oligos found after filtering NaN values")
        return tables

    # 1. Perfect donor lookup (vectorized groupby)
    tables.perfect_donor_to_oligos = (
        valid_df.groupby('donor')['oligo_name']
        .apply(lambda x: x.drop_duplicates().tolist())
        .to_dict()
    )

    # 2. Guide lookup (vectorized groupby)
    tables.guide_to_oligos = (
        valid_df.groupby('guide')['oligo_name']
        .apply(lambda x: x.drop_duplicates().tolist())
        .to_dict()
    )

    # 3. Combined guide+donor lookup (vectorized)
    valid_df['guide_donor_key'] = (
        valid_df['guide'] + SPCAS9_CLONING_SITE + valid_df['donor']
    )
    # NN 250 / C20: `dict(zip(...))` is LAST-WINS and uncounted -- the same defect class as
    # 04's oligo-attribute map. The map is kept (consumers rely on it) but the collisions are
    # COUNTED and surfaced in get_stats(), so "two designs share one guide+donor key" can never
    # again be resolved silently by row order.
    tables.guide_donor_key_collisions = int(
        len(valid_df) - valid_df['guide_donor_key'].nunique()
    )
    tables.guide_donor_to_oligo = dict(
        zip(valid_df['guide_donor_key'], valid_df['oligo_name'])
    )
    if tables.guide_donor_key_collisions:
        logger.warning(
            f"  {tables.guide_donor_key_collisions:,} guide+donor keys are shared by more than "
            f"one oligo_name; guide_donor_to_oligo keeps whichever row came LAST. Counted in "
            f"get_stats()['guide_donor_key_collisions'] (NN 250 / C20)."
        )

    # 4. Sliding window lookup - need to iterate but optimize
    # Build using a set accumulator then convert to dict
    window_accumulator = {}

    # Vectorized length check
    valid_df['donor_len'] = valid_df['donor'].str.len()
    windowed_df = valid_df[valid_df['donor_len'] >= window_size]

    # Process using itertuples (faster than iterrows)
    for row in windowed_df[['oligo_name', 'donor']].itertuples(index=False):
        oligo_name = row.oligo_name
        donor = row.donor
        for i in range(len(donor) - window_size + 1):
            window = donor[i:i + window_size]
            if window not in window_accumulator:
                window_accumulator[window] = set()
            window_accumulator[window].add(oligo_name)

    # Convert sets to lists
    tables.window_to_oligos = {k: list(v) for k, v in window_accumulator.items()}

    # 5. oligo_name -> installed variant (individual_variants). PER-ROW fallback
    # to the donor sequence where the variant id is missing -- the variant
    # column is only ~87% populated, and a plain .astype(str) would map every
    # NaN to the single string "nan", collapsing genuinely-distinct unannotated
    # edits into one fake variant (inflating is_unambiguous in the sliding-
    # window branch). Two oligos sharing a donor install the same edit, so the
    # donor sequence is the correct conservative identity when the explicit
    # annotation is absent. If the whole column is missing (legacy design
    # files), fall back to donor sequence for every row.
    if variant_col in valid_df.columns:
        # 🔴 NN 250 / D15: qualify the identity by CHROMOSOME and require it to be POSITIONAL.
        # A bare position collides across chromosomes (50 ids measured on the real table) and a
        # class label such as `PTC` is not an identity at all (3,054 rows). Rows that cannot be
        # made positional fall through to the donor-sequence fallback below, which is what this
        # builder already did for unannotated rows -- the label is never kept as an identity.
        _chrom = valid_df['chrom'] if 'chrom' in valid_df.columns else pd.Series(
            [None] * len(valid_df), index=valid_df.index)
        _vcoord = valid_df['variant_coord'] if 'variant_coord' in valid_df.columns else pd.Series(
            [None] * len(valid_df), index=valid_df.index)
        _has_chrom = 'chrom' in valid_df.columns
        qualified = [
            qualified_variant_id(c, v, vc, allow_unqualified=not _has_chrom)
            for c, v, vc in zip(_chrom, valid_df[variant_col], _vcoord)
        ]
        qualified = pd.Series(qualified, index=valid_df.index, dtype=object)
        variant_series = qualified.where(qualified.notna(), valid_df['donor']).astype(str)
        n_fallback = int(qualified.isna().sum())
        n_raw_missing = int(valid_df[variant_col].isna().sum())
        if not _has_chrom:
            logger.warning(
                "  No 'chrom' column: variant ids CANNOT be chromosome-qualified, so positional "
                "ids are used unqualified and two edits at the same coordinate on different "
                "chromosomes still merge (50 such ids measured on the real harmonized table). "
                "Class labels are still routed to the donor fallback. Load 'chrom' and "
                "'variant_coord' with the oligo table to close the collision (NN 250 / D15)."
            )
        if n_fallback:
            logger.info(
                f"  Variant id '{variant_col}' is not a chromosome-qualified position for "
                f"{n_fallback:,} ({n_fallback/len(valid_df):.1%}) oligos "
                f"({n_raw_missing:,} of those were blank; the rest are class labels with no "
                f"coordinate); using donor sequence as the per-row fallback identity for those."
            )
    else:
        logger.info(
            f"  Variant column '{variant_col}' not found; using donor sequence "
            f"as the variant identity for unambiguity assessment."
        )
        variant_series = valid_df['donor']
    tables.oligo_to_variant = dict(zip(valid_df['oligo_name'], variant_series))

    # NN 250 / D13: the PAM TUPLE, which is what NN 241's attack differed on (same gene, same
    # chromosome, same coordinate, same ref->alt; different PAM / PAM_strand / PAM_coord).
    _pam_cols = [c for c in ("PAM", "PAM_coord", "PAM_strand") if c in valid_df.columns]
    if _pam_cols:
        tables.oligo_to_pam = dict(zip(
            valid_df['oligo_name'],
            zip(*(valid_df[c].astype(str) for c in _pam_cols)),
        ))

    stats = tables.get_stats()
    logger.info(f"  Built lookup tables:")
    logger.info(f"    Unique donors: {stats['num_unique_donors']}")
    logger.info(f"    Unique guides: {stats['num_unique_guides']}")
    logger.info(f"    Unique windows: {stats['num_unique_windows']}")

    return tables


def match_donor_sliding_window(
    donor: str,
    tables: OligoLookupTables
) -> Set[str]:
    """
    Match an observed donor using sliding window.

    Slides a window across the OBSERVED donor and looks up each window
    in the precomputed lookup table. This handles:
    - Truncated donors (5' or 3' deletions)
    - Extended donors (chimeras)
    - Internal mutations (as long as 60bp window is intact)

    Args:
        donor: The observed donor sequence
        tables: Precomputed lookup tables

    Returns:
        Set of matching oligo names
    """
    matches = set()
    window_size = tables.window_size

    if len(donor) < window_size:
        # Donor too short for sliding window
        return matches

    # Slide window across OBSERVED donor
    for i in range(len(donor) - window_size + 1):
        window = donor[i:i + window_size]
        if window in tables.window_to_oligos:
            matches.update(tables.window_to_oligos[window])

    return matches


def _variant_set(oligos, tables: OligoLookupTables) -> Set[str]:
    """Collapse a collection of candidate oligo_names to their installed variants."""
    if not tables.oligo_to_variant:
        # No variant table -> treat each oligo as its own variant (legacy).
        return set(oligos)
    return {tables.oligo_to_variant.get(o, o) for o in oligos}


def _resolve_candidates(
    result: OligoMatch,
    candidates,
    tables: OligoLookupTables,
    method_unique: str,
    method_ambiguous: str,
) -> None:
    """
    Populate `result` from a candidate oligo set using VARIANT-level unambiguity.

    Candidate oligos that share the same installed variant (donor_changes) are
    redundant guides for one edit; a match is unambiguous when they all agree on
    the variant, even if >1 oligo_name is present. The reported `oligo_name` is
    chosen deterministically (sorted-first) among the candidates so reruns are
    reproducible.

    🔴 A CONFLICTING GUIDE VETOES UNAMBIGUITY (NN 250 / D12, from the NN 241 attack). Every
    match branch sets ``result.guide_status`` BEFORE calling this, so the veto lives here once
    rather than in each branch. ``guide_status == conflict`` means the read carries a DESIGNED
    guide that names no design the donor allows: the donor answer is still kept and reported, but
    the call is internally contradicted and must not tell a downstream ``is_unambiguous`` filter
    that it is clean. The attack it closes: one sequencing error turns a read into another
    design's exact donor while the read still carries its own design's guide, and the row shipped
    ``match_method='perfect', is_unambiguous=True``. ``no_guide`` and ``guide_not_designed`` are
    NOT evidence against the donor and are deliberately left alone.
    """
    candidates = sorted(set(candidates))
    variants = _variant_set(candidates, tables)
    result.candidate_oligos = candidates
    result.num_candidate_oligos = len(candidates)
    result.num_candidate_variants = len(variants)
    result.oligo_name = candidates[0]

    # 🔴 NN 250 / D13 -- a FLAG, deliberately NOT a veto on is_unambiguous. See
    # planning/222_triumph_embl_library_design/250_bc1_pipeline_fixes.md §9: within one plasmid
    # library, 35,004 of 35,885 multi-oligo donor groups (97.5%) disagree on this tuple, because
    # the library is DESIGNED with redundant guides at different PAM sites installing the same
    # edit. Vetoing would strip unambiguity from the library's own design. The damage NN 241
    # measured is to PAM ATTRIBUTION (the per-well PAM-class tables, all derived from oligo_name)
    # at ~17% while the EDIT is right -- so the PAM-class consumers get a flag to read, and edit
    # identity is left alone. None (not True) when the table carries no PAM columns: unknown is a
    # third state and must not read as agreement.
    if tables.oligo_to_pam:
        pams = {tables.oligo_to_pam.get(o) for o in candidates}
        pams.discard(None)
        result.num_candidate_pams = len(pams)
        result.candidates_agree_on_pam = (len(pams) == 1) if pams else None

    if len(variants) == 1:
        result.is_unambiguous = result.guide_status != GUIDE_CONFLICT
        result.variant = next(iter(variants))
        result.match_method = method_unique
    else:
        result.match_method = method_ambiguous


def match_sequence(
    donor: Optional[str],
    guide: Optional[str],
    tables: OligoLookupTables,
    pre_trim: bool = True,
) -> OligoMatch:
    """
    Match a donor (and optionally guide) to designed oligos.

    Matching priority:
    1. Perfect donor match (highest confidence)
    2. Sliding window match on donor
    3. Guide-only match (if donor fails but guide is available)

    Unambiguity is assessed at the VARIANT level (installed edit / donor_changes)
    rather than oligo_name: oligos that share a donor install the same edit and
    differ only by which guide cut, so they are collapsed. Guide, when present,
    is still used to narrow candidates before the variant collapse.

    Args:
        donor: The observed donor sequence (can be None)
        guide: The guide sequence if available (can be None)
        tables: Precomputed lookup tables
        pre_trim: If True (default), anchor-trim the observed donor down to the
            variable region before lookup. Anchor-conditional and idempotent:
            already-trimmed / anchorless donors pass through untouched, so the
            gdb path is unaffected.

    Returns:
        OligoMatch with results and confidence indicators
    """
    result = OligoMatch()

    if donor is None and guide is None:
        return result

    # Normalize: STRIP, then uppercase.
    #
    # 🔴 NN 250 / D16. Stripping was missing, and one TRAILING space destroyed the perfect-donor
    # route: measured on the two-design fixture, `perfect_donor` went True -> False and the row
    # fell through to the 60-nt sliding window. The NN 241 attack measured the consequence on real
    # data -- perfect-donor matching 100% -> 0% and accuracy 68.6% -> 40.9% -- which presents as a
    # FAILED SYNTHESIS BATCH rather than as a parsing bug. A LEADING space survived, because the
    # anchor trim discards everything before the anchor, so the defect was asymmetric and the
    # easier half of it was invisible.
    if donor:
        donor = donor.strip().upper()
        if pre_trim:
            donor = trim_to_variable_region(
                donor,
                pos_slack=tables.anchor_pos_slack,
                max_mismatch=tables.anchor_max_mismatch,
            )
    if guide:
        # A padded guide would otherwise read as an UNDESIGNED guide, which silently suppresses
        # the D12 conflict guard instead of firing it.
        guide = guide.strip().upper()

    # Guide state, independent of whether the guide ends up narrowing anything.
    guide_oligos = tables.guide_to_oligos.get(guide) if guide else None
    if not guide:
        result.guide_status = GUIDE_NO_GUIDE
    elif guide_oligos is None:
        result.guide_status = GUIDE_NOT_DESIGNED
    else:
        # provisional: a designed guide that never intersects the donor's candidates is a CONFLICT
        result.guide_status = GUIDE_CONFLICT
        result.num_guide_oligos = len(guide_oligos)

    # 1. Try perfect donor match
    if donor and donor in tables.perfect_donor_to_oligos:
        candidates = tables.perfect_donor_to_oligos[donor]
        result.perfect_donor = True
        result.perfect_middle_region = True

        # Narrow by guide first if it disambiguates further (rarely available).
        if guide_oligos is not None:
            overlap = set(candidates) & set(guide_oligos)
            if overlap:
                result.guide_match = True
                result.guide_status = GUIDE_MATCH
                candidates = list(overlap)
            # else: guide_status stays GUIDE_CONFLICT -- the donor answer is kept, but the
            # disagreement is now recorded instead of discarded.

        _resolve_candidates(
            result, candidates, tables,
            method_unique="perfect",
            method_ambiguous="perfect_ambiguous",
        )
        return result

    # 2. Try sliding window match on donor
    if donor:
        window_matches = match_donor_sliding_window(donor, tables)

        if window_matches:
            result.perfect_middle_region = True

            if guide_oligos is not None:
                overlap = window_matches & set(guide_oligos)
                if overlap:
                    result.guide_match = True
                    result.guide_status = GUIDE_MATCH
                    window_matches = overlap
                # else: guide_status stays GUIDE_CONFLICT

            _resolve_candidates(
                result, window_matches, tables,
                method_unique="sliding_window",
                method_ambiguous="sliding_window_ambiguous",
            )
            return result

    # 3. Guide-only match (fallback): nothing to conflict WITH, so a designed guide is a match.
    if guide_oligos is not None:
        candidates = guide_oligos
        result.guide_match = True
        result.guide_status = GUIDE_MATCH
        _resolve_candidates(
            result, candidates, tables,
            method_unique="guide_only",
            method_ambiguous="guide_only_ambiguous",
        )
        return result

    return result


def _match_row(args):
    """Helper function for parallel matching."""
    donor, guide, tables = args
    if pd.isna(donor):
        return None
    return match_sequence(donor, guide if pd.notna(guide) else None, tables)


CANDIDATE_LIST_CAP = 50
CAS12A_OLIGO_PREFIXES = ("LbCas12a_", "impLbCas12a_")


def repair_gdb01_shifted_guides(guides: pd.Series, tables: "OligoLookupTables"):
    """BRIDGE, retired 2026-09-11: undo gdb Step 02's 1-nt-early guide extraction on Cas12a reads (planning 30_).

    🔴 THE ROOT FIX IS IN gdb STEP 02, not here. ``magestic.utils.sequence_anchors.find_re_site`` now tries every
    cloning-site pattern exactly before trying any approximately, and scores approximate matches by substitutions
    only, so the shift cannot happen. Tables built from 2026-09-11 onward need no repair and this function will
    restore ~0 rows against them.

    This is kept, and reachable only via bdb 04's ``--gdb01-repair`` (OFF by default), so guide-donor-bc0 tables
    built BEFORE that fix stay readable. It is a bridge in the sense of PIPELINE_CORRECTNESS_POLICY.md §1: it is
    labelled in the output (``guide_repaired_gdb01`` / ``guide_gdb_raw``), switchable, and it protects exactly one
    consumer while leaving the wrong guides on disk for every other reader. Re-run gdb Step 02 instead.

    gdb Step 02 searches the SpCas9 RE site first; on a Cas12a read the fuzzy fallback ties the true TTTCGAAGAGC window
    with the window one base earlier whenever the guide ends in G, and keeps the earlier one. The extracted guide is
    then 'T' + guide[:19] (the base before a Cas12a guide is T). Measured on the rescued Oak tables bdb 03 consumes
    (V438 LbCas12a / V439 impLbCas12a): for G-ending designs 99.3% / 98.8% of reads carry the shifted string; A/C/T-ending
    designs are >99% exact.

    Repair rule (deterministic, conservative): a guide that is NOT a designed guide, starts with 'T', and whose
    guide[1:] + 'G' IS a designed guide belonging only to Cas12a oligos, is replaced by that designed guide.
    Returns (guides_out, repaired_flag). Guides that are themselves designed are never touched."""
    known = tables.guide_to_oligos
    up = guides.where(guides.isna(), guides.astype(str).str.upper())

    def _fix(x):
        if not isinstance(x, str) or not x or x in known or x[0] != "T":
            return None
        r = x[1:] + "G"
        oligos = known.get(r)
        if oligos and all(str(o).startswith(CAS12A_OLIGO_PREFIXES) for o in oligos):
            return r
        return None

    uniq = up.dropna().unique()
    fixmap = {u: _fix(u) for u in uniq}
    rep = up.map(fixmap)
    flag = rep.notna()
    return guides.where(~flag, rep), flag


def format_candidates(candidates, n_variants: int) -> Optional[str]:
    """';'-joined candidate oligo_names for a VARIANT-ambiguous match (n_variants > 1), else None.

    The reported oligo_name of an ambiguous match is only the sorted-first candidate (_resolve_candidates), so
    the candidates must ship with the row for a consumer to see that other edits were equally consistent with
    the read (planning 30_ arm_bdb-F1). Capped at CANDIDATE_LIST_CAP names plus a '+N_more' token."""
    if n_variants <= 1 or not candidates:
        return None
    c = sorted(set(candidates))
    if len(c) > CANDIDATE_LIST_CAP:
        return ";".join(c[:CANDIDATE_LIST_CAP]) + f";+{len(c) - CANDIDATE_LIST_CAP}_more"
    return ";".join(c)


def _get_perfect_result_map(tables: OligoLookupTables) -> Dict[str, tuple]:
    """donor -> (oligo_name, is_unambiguous, variant, num_candidate_oligos,
    num_candidate_variants, match_method, candidate_oligos) for the perfect-donor, NO-guide case.

    Precomputes _resolve_candidates' variant-level collapse ONCE per `tables`
    (cached on the object) so match_dataframe's vectorized fast path can .map()
    it. Encodes exactly match_sequence's perfect-donor branch without guide
    narrowing (perfect_donor=True, perfect_middle_region=True, guide_match=False).
    """
    cached = getattr(tables, "_perfect_result_map", None)
    if cached is not None:
        return cached
    o2v = tables.oligo_to_variant
    m: Dict[str, tuple] = {}
    for donor, oligos in tables.perfect_donor_to_oligos.items():
        cands = sorted(set(oligos))
        variants = {o2v.get(o, o) for o in cands} if o2v else set(cands)
        nv = len(variants)
        # NN 250 / D13: the PAM verdict is precomputed here too, or the perfect-donor majority
        # (which takes the vectorized fast path) would silently carry None for every row.
        if tables.oligo_to_pam:
            _pams = {tables.oligo_to_pam.get(o) for o in cands}
            _pams.discard(None)
            _pam_ok = (len(_pams) == 1) if _pams else None
            _npam = len(_pams)
        else:
            _pam_ok, _npam = None, 0
        m[donor] = (
            cands[0],
            nv == 1,
            (next(iter(variants)) if nv == 1 else None),
            len(cands),
            nv,
            "perfect" if nv == 1 else "perfect_ambiguous",
            format_candidates(cands, nv),
            _pam_ok,
            _npam,
        )
    try:
        tables._perfect_result_map = m
    except Exception:
        pass
    return m


#: Fraction of GUIDE-BEARING rows whose guide is not any designed guide in scope, above which the
#: run is warned. A guide that matches no design in scope came from another library, so a high
#: fraction means the lookup scope and the sample do not belong to each other.
#:
#: 🔴 PROVISIONAL, and labelled as such: 0.5 is a structural bound, not a measured base rate. The
#: NN 241 cross-library-contamination attack absorbed at 46.7% assigned, so a run whose guides are
#: majority-undesigned is certainly wrong; the true per-screen base rate has NOT been measured
#: here, and the real-data run in §10 is where it should be replaced with a measured value.
GUIDE_NOT_DESIGNED_WARN_FRACTION = 0.5


def guide_provenance_report(status_counts) -> Dict[str, object]:
    """Run-level guide provenance: what fraction of GUIDE-BEARING rows carry an undesigned guide.

    🔴 NN 250 / D14. `guide_not_designed` was already computed per row and never checked in
    aggregate, so cross-library contamination had nothing to trip: the NN 241 attack absorbed at
    46.7% assigned with 47 rows still flagged unambiguous.

    Rows with NO guide are excluded from the denominator -- they are not evidence either way, and
    leaving them in would dilute the fraction to nothing on the production bc1-donor-bc0 input,
    which mostly has no guide column at all.
    """
    nd = int(status_counts.get(GUIDE_NOT_DESIGNED, 0))
    bearing = nd + int(status_counts.get(GUIDE_MATCH, 0)) + int(status_counts.get(GUIDE_CONFLICT, 0))
    frac = (nd / bearing) if bearing else 0.0
    return {
        "guide_bearing_rows": bearing,
        "guide_not_designed": nd,
        "guide_not_designed_fraction": frac,
        "threshold": GUIDE_NOT_DESIGNED_WARN_FRACTION,
        "exceeds_threshold": bool(bearing and frac > GUIDE_NOT_DESIGNED_WARN_FRACTION),
    }


def donor_length_report(observed_lengths, designed_lengths) -> Dict[str, object]:
    """How many observed trimmed-donor lengths are OFF the designed set, and the histogram.

    🔴 NN 250 / D16, second half. A systematic length problem -- stray whitespace, a truncation,
    an off-by-one trim -- removes the perfect-donor route for EVERY read at once, and the run then
    presents as a failed synthesis batch rather than as a parsing bug (the NN 241 attack measured
    perfect-donor matching 100% -> 0% and accuracy 68.6% -> 40.9% from a single trailing space).
    The lengths are the cheapest place to see it, so they are reported rather than assumed.
    """
    designed = sorted({int(x) for x in designed_lengths})
    hist: Dict[int, int] = {}
    for L in observed_lengths:
        try:
            Li = int(L)
        except (TypeError, ValueError):
            continue
        hist[Li] = hist.get(Li, 0) + 1
    n = sum(hist.values())
    off = sum(c for L, c in hist.items() if L not in set(designed))
    return {
        "designed_lengths": designed,
        "length_histogram": hist,
        "n_donors": n,
        "off_design_lengths": off,
        "off_design_fraction": (off / n) if n else 0.0,
    }


#: Fraction of off-design donor lengths above which the run is warned about. Not a filter and not
#: a threshold on any result: nothing is dropped, the run is only told. Set deliberately high
#: (half the donors) because legitimate inputs carry truncated and chimeric donors; the condition
#: this is aimed at -- a systematic trim/whitespace fault -- puts it at ~1.0.
DONOR_LENGTH_WARN_FRACTION = 0.5


def match_dataframe(
    df: pd.DataFrame,
    tables: OligoLookupTables,
    donor_col: str = "donor",
    guide_col: str = "guide",
    prefix: str = "",
    n_jobs: int = 1,
    chunk_size: int = 10000,
    pre_trim: bool = True,
) -> pd.DataFrame:
    """
    Match all rows in a DataFrame to designed oligos.

    Adds columns for match results:
    - {prefix}oligo_name: Matched oligo name
    - {prefix}perfect_donor: True if perfect donor match
    - {prefix}perfect_middle_region: True if middle region matched
    - {prefix}guide_match: True if guide contributed to match (kept for back-compat)
    - {prefix}guide_status: no_guide / guide_not_designed / match / conflict -- see OligoMatch
    - {prefix}num_guide_oligos: how many designs carry this guide (0 when it is not a designed guide)
    - {prefix}is_unambiguous: True if match is unambiguous
    - {prefix}match_method: Description of how match was made
    - {prefix}num_candidates / {prefix}num_candidate_variants: candidate oligos / distinct installed variants
    - {prefix}variant: installed variant id when unambiguous
    - {prefix}candidate_oligos: ';'-joined candidates when num_candidate_variants > 1, else None

    PERFORMANCE: Supports parallel processing with n_jobs > 1 and
    processes in chunks for memory efficiency.

    Args:
        df: DataFrame with donor and optionally guide columns
        tables: Precomputed lookup tables
        donor_col: Name of donor column
        guide_col: Name of guide column (can be missing)
        prefix: Prefix for output column names (e.g., "bdb_")
        n_jobs: Number of parallel workers (default: 1 = sequential)
        chunk_size: Size of chunks for processing (default: 10000)

    Returns:
        DataFrame with added match columns
    """
    from tqdm import tqdm
    import numpy as np

    logger.info(f"Matching {len(df)} sequences to oligos (n_jobs={n_jobs})...")

    df = df.copy()
    n = len(df)
    has_guide = guide_col in df.columns

    # Result-column arrays for a single vectorized writeback at the end
    # (replaces the old per-row df.iloc x 9 assignment).
    oligo = [None] * n
    pdon = [False] * n
    pmid = [False] * n
    gmatch = [False] * n
    gstatus = [GUIDE_NO_GUIDE] * n
    ngoligo = [0] * n
    unamb = [False] * n
    method = ["none"] * n
    ncand = [0] * n
    ncandv = [0] * n
    var = [None] * n
    cands = [None] * n
    pam_ok = [None] * n          # NN 250 / D13
    npam = [0] * n

    # ---------------------------------------------------------------------
    # VECTORIZED FAST PATH (2026-06-02) — the perfect-donor majority.
    # Replicates match_sequence's perfect-donor branch WITHOUT guide-narrowing
    # (the dominant case: production bc1_donor_bc0 input has no guide column).
    # Exact-anchor trim is vectorized (str.split == trim_to_variable_region's
    # exact pass); the precomputed donor->result map encodes _resolve_candidates'
    # variant-level collapse. EVERY other row (no exact anchor / tolerant-trim /
    # sliding-window / guide-disambiguating / no perfect hit) falls back to the
    # UNCHANGED per-row match_sequence below, so output is identical to the
    # all-per-row path (validated old==new on real data).
    # ---------------------------------------------------------------------
    pm = _get_perfect_result_map(tables)
    fast_mask = np.zeros(n, dtype=bool)
    mapped_vals = None
    if pre_trim and n:
        donor_s = df[donor_col]
        # NN 250 / D16: strip BEFORE upper, matching match_sequence exactly. Without it one
        # trailing space moved a row off the perfect-donor fast path.
        upu = donor_s.astype(str).str.strip().str.upper().where(donor_s.notna())
        a1, a2 = DONOR_VARIABLE_REGION_ANCHORS[0], DONOR_VARIABLE_REGION_ANCHORS[1]
        t1 = upu.str.split(a1, n=1).str[1]   # part after first anchor1, else NaN
        t2 = upu.str.split(a2, n=1).str[1]   # part after first anchor2, else NaN
        trimmed = t1.where(t1.notna(), t2)   # anchor1 priority (matches trim())
        mapped = trimmed.map(pm)             # result-tuple or NaN

        # NN 250 / D16: report the trimmed-donor length distribution against the designed set.
        # A systematic trim/whitespace fault takes the perfect route away from EVERY read at
        # once, and that presents as failed synthesis unless something says otherwise. Nothing
        # is dropped or filtered here -- the run is only told.
        _designed = {len(d) for d in tables.perfect_donor_to_oligos}
        _trimmed_ok = trimmed.dropna()
        # `.str` raises on a series that carries no strings at all (an all-NaN column comes back
        # as float64), so the emptiness check has to come FIRST -- it is not merely an
        # optimization. This cost 7 tests before it was guarded.
        if _designed and len(_trimmed_ok):
            _rep = donor_length_report(_trimmed_ok.astype(str).str.len().tolist(), _designed)
            if _rep["n_donors"] and _rep["off_design_fraction"] > DONOR_LENGTH_WARN_FRACTION:
                _top = sorted(_rep["length_histogram"].items(), key=lambda kv: -kv[1])[:5]
                logger.warning(
                    "donor length: %.1f%% of %d trimmed donors have an OFF-DESIGN length "
                    "(designed: %s; most common observed: %s). A systematic trim or whitespace "
                    "fault looks exactly like a failed synthesis batch -- check the parse before "
                    "reading this as biology (NN 250 / D16).",
                    100 * _rep["off_design_fraction"], _rep["n_donors"],
                    _rep["designed_lengths"], _top,
                )
        if has_guide:
            guide_s = df[guide_col]
            g_present = guide_s.notna().to_numpy(dtype=bool)
            g_block = guide_s.astype(str).str.strip().str.upper().where(guide_s.notna()).map(
                lambda x: isinstance(x, str) and x in tables.guide_to_oligos
            ).fillna(False).to_numpy(dtype=bool)
        else:
            g_present = np.zeros(n, dtype=bool)
            g_block = np.zeros(n, dtype=bool)
        fast_mask = mapped.notna().to_numpy(dtype=bool) & (~g_block)
        mapped_vals = mapped.to_numpy()
        for i in np.flatnonzero(fast_mask):
            t = mapped_vals[i]  # (oligo, is_unambiguous, variant, n_oligo, n_var, method, candidates)
            oligo[i] = t[0]; pdon[i] = True; pmid[i] = True
            unamb[i] = bool(t[1]); var[i] = t[2]
            ncand[i] = int(t[3]); ncandv[i] = int(t[4]); method[i] = t[5]; cands[i] = t[6]
            pam_ok[i] = t[7]; npam[i] = int(t[8])          # NN 250 / D13
            # the fast path is only taken where the guide is NOT a designed guide (~g_block), so the
            # guide can never have narrowed here: it is absent, or present and undesigned.
            gstatus[i] = GUIDE_NOT_DESIGNED if g_present[i] else GUIDE_NO_GUIDE

    # ---------------------------------------------------------------------
    # FALLBACK — unchanged per-row match_sequence for every non-fast row.
    # NaN-donor rows are skipped (None result), matching the original which did
    # `match_sequence(d, g) if pd.notna(d) else None` (guide-only never fires on
    # a NaN donor).
    # ---------------------------------------------------------------------
    donor_vals = df[donor_col].to_numpy()
    guide_vals = df[guide_col].to_numpy() if has_guide else None
    rest = np.flatnonzero(~fast_mask)
    for i in tqdm(rest, desc="Matching (fallback)", disable=(len(rest) == 0)):
        d = donor_vals[i]
        d = d if isinstance(d, str) else None
        if d is None:
            continue
        g = guide_vals[i] if has_guide else None
        g = g if (has_guide and isinstance(g, str)) else None
        r = match_sequence(d, g, tables, pre_trim=pre_trim)
        oligo[i] = r.oligo_name; pdon[i] = r.perfect_donor; pmid[i] = r.perfect_middle_region
        gmatch[i] = r.guide_match; unamb[i] = r.is_unambiguous; method[i] = r.match_method
        gstatus[i] = r.guide_status; ngoligo[i] = r.num_guide_oligos
        ncand[i] = r.num_candidate_oligos; ncandv[i] = r.num_candidate_variants; var[i] = r.variant
        cands[i] = format_candidates(r.candidate_oligos, r.num_candidate_variants)
        pam_ok[i] = r.candidates_agree_on_pam; npam[i] = r.num_candidate_pams   # NN 250 / D13

    # Single vectorized writeback
    df[f'{prefix}oligo_name'] = oligo
    df[f'{prefix}perfect_donor'] = pdon
    df[f'{prefix}perfect_middle_region'] = pmid
    df[f'{prefix}guide_match'] = gmatch
    df[f'{prefix}guide_status'] = gstatus
    df[f'{prefix}num_guide_oligos'] = ngoligo
    df[f'{prefix}is_unambiguous'] = unamb
    df[f'{prefix}match_method'] = method
    df[f'{prefix}num_candidates'] = ncand
    df[f'{prefix}num_candidate_variants'] = ncandv
    df[f'{prefix}variant'] = var
    df[f'{prefix}candidate_oligos'] = cands
    # 🔴 NN 250 / D13: emit these ONLY when the lookup actually carries PAM data. Writing them
    # unconditionally put `None` / `0` on every row of any run whose oligo table lacks the PAM
    # columns -- a CONSTANT COLUMN, which names nothing, and precisely the failure D17 exists to
    # avoid. `object` dtype deliberately: None (unknown) must survive as None rather than being
    # coerced to False, which would read as "the PAMs disagree".
    if tables.oligo_to_pam:
        df[f'{prefix}candidates_agree_on_pam'] = pd.Series(pam_ok, index=df.index, dtype=object)
        df[f'{prefix}num_candidate_pams'] = npam
    else:
        logger.info("  no PAM columns in the oligo lookup: candidates_agree_on_pam NOT emitted "
                    "(an all-None column would name nothing; NN 250 / D13)")

    n_matched = sum(1 for o in oligo if o is not None)
    n_perfect = sum(1 for x in pdon if x)
    logger.info(f"  Matched {n_matched}/{n} (perfect_donor={n_perfect}; "
                f"vectorized fast-path={int(fast_mask.sum())}, fallback={int((~fast_mask).sum())})")
    status_counts = {k: 0 for k in GUIDE_STATUSES}
    for st in gstatus:
        status_counts[st] = status_counts.get(st, 0) + 1
    logger.info(f"  guide_status: {status_counts}")
    if status_counts[GUIDE_CONFLICT]:
        logger.warning(
            f"  {status_counts[GUIDE_CONFLICT]:,} rows have guide_status='conflict': a DESIGNED guide "
            f"naming no design the donor allows. The donor answer was kept. These were invisible before "
            f"(guide_match=False meant this and 'no guide' alike). Since NN 250 / D12 these rows also "
            f"carry is_unambiguous=False."
        )
    # NN 250 / D14: the run-level guide-provenance check.
    _gp = guide_provenance_report(status_counts)
    if _gp["exceeds_threshold"]:
        logger.warning(
            "  🔴 %.1f%% of the %d guide-bearing rows have guide_status='guide_not_designed' "
            "(threshold %.0f%%). A guide that is not ANY designed guide in scope is a guide from "
            "ANOTHER library: this is the signature of cross-library contamination, which absorbed "
            "silently at 46.7%% assigned with 47 rows flagged unambiguous in the NN 241 attack. "
            "Check --v-libraries and the sample sheet before reading any assignment from this run "
            "(NN 250 / D14).",
            100 * _gp["guide_not_designed_fraction"], _gp["guide_bearing_rows"],
            100 * _gp["threshold"],
        )
    return df


def load_and_build_lookup(
    oligo_files: Optional[List[Path]] = None,
    window_size: int = DEFAULT_WINDOW_SIZE,
    use_harmonized_table: bool = True,
    v_libraries: Optional[List[str]] = None,
    nuclease_types: Optional[List[str]] = None,
    anchor_pos_slack: int = DEFAULT_ANCHOR_POS_SLACK,
    anchor_max_mismatch: int = DEFAULT_ANCHOR_MAX_MISMATCH,
) -> OligoLookupTables:
    """
    Load oligo designs and build unified lookup tables.

    PREFERRED: Use use_harmonized_table=True (default) to load from the
    pre-parsed harmonized table with full annotations.

    LEGACY: Provide oligo_files to load from individual design files.

    Args:
        oligo_files: List of paths to oligo design files (legacy mode)
        window_size: Size of sliding window
        use_harmonized_table: Use harmonized table (default True)
        v_libraries: Filter by V library when using harmonized table
        nuclease_types: Filter by nuclease type when using harmonized table

    Returns:
        Combined OligoLookupTables
    """
    if use_harmonized_table and oligo_files is None:
        # PREFERRED: Load from harmonized table
        combined_df = load_oligos_from_harmonized_table(
            v_libraries=v_libraries,
            nuclease_types=nuclease_types,
        )
    elif oligo_files:
        # LEGACY: Load from individual files
        if use_harmonized_table:
            logger.info("Note: oligo_files provided, ignoring use_harmonized_table")

        all_oligos = []
        for filepath in oligo_files:
            if filepath.exists():
                df = load_oligo_design(filepath)
                all_oligos.append(df)
            else:
                logger.warning(f"Oligo file not found: {filepath}")

        if not all_oligos:
            raise ValueError("No oligo files loaded")

        combined_df = pd.concat(all_oligos, ignore_index=True)
        logger.info(f"Combined oligos from {len(oligo_files)} files")
    else:
        raise ValueError(
            "Either provide oligo_files or use use_harmonized_table=True"
        )

    # Remove duplicates by oligo_name
    combined_df = combined_df.drop_duplicates(subset=['oligo_name'])
    logger.info(f"Building lookup from {len(combined_df)} unique oligos")

    return build_lookup_tables(combined_df, window_size=window_size,
                               anchor_pos_slack=anchor_pos_slack,
                               anchor_max_mismatch=anchor_max_mismatch)


def build_lookup_from_harmonized(
    v_libraries: Optional[List[str]] = None,
    nuclease_types: Optional[List[str]] = None,
    window_size: int = DEFAULT_WINDOW_SIZE,
) -> OligoLookupTables:
    """
    Convenience function to build lookup directly from harmonized table.

    This is the simplest way to get started with oligo matching.

    Args:
        v_libraries: Filter by V library (V416, V419, V520, V521). None = all.
        nuclease_types: Filter by nuclease type. None = all.
        window_size: Size of sliding window

    Returns:
        OligoLookupTables ready for matching
    """
    return load_and_build_lookup(
        oligo_files=None,
        window_size=window_size,
        use_harmonized_table=True,
        v_libraries=v_libraries,
        nuclease_types=nuclease_types,
    )


# ---------------------------------------------------------------------------
# SPS scoping: match a row only against the subpool its own molecule came from
# ---------------------------------------------------------------------------
#
# Every bc1-donor-bc0 read carries the subpool priming site (SPS) between the donor and bc0, and the SPS
# names which design SUBPOOL the molecule was amplified from. Matching every row against every subpool at
# once lets a donor tie be broken by alphabetical order across subpools instead: measured on the
# 2026-09-10 rebuild, 313,570 of 4,458,480 named rows with a readable SPS (7.0%) carried a name from a
# different subpool than the read itself came from, and the transformation independently agreed with the
# SPS on 99.3% of them.
#
# 🔴 A POST-HOC CANDIDATE FILTER CANNOT FIX THIS. An UNAMBIGUOUS row carries no candidate list -- the
# right design was never a candidate to begin with. The SEARCH has to be scoped, which is why this is a
# root fix and not a bridge. Ordered by Kevin 2026-09-10: "SPS scoping in 04: implement and re-run."
#
# 🔴 SPS NAMES THE SUBPOOL ONLY (policy §4). V450 and V454 share one SPS; V451 and V455 share the other.
# A scaffold claim can NEVER be read out of an SPS -- scaffold comes from the transformation (the cut
# vector) or from the guide-donor-bc0 linkage.


#: Accepted names for the library column of an SPS mapping, in priority order. The SHIPPED file
#: (/oak/.../common/keyfiles/sps_library_mapping.tsv) uses ``V_library``; ``library_id`` is the
#: canonical name this loader was written against. Checked BY NAME, never positionally.
SPS_LIBRARY_COLUMNS = ("library_id", "V_library")


def load_sps_mapping(path: str, v_libraries: Optional[Sequence[str]] = None) -> Dict[str, str]:
    """Read an ``sps -> library_id`` TSV into an uppercased dict, optionally restricted to a scope.

    Raises rather than silently returning an empty map: a mistyped path that yields ``{}`` would turn
    SPS scoping off without saying so, and the run would look like it had scoped when it had not.
    """
    m = pd.read_csv(path, sep="\t", dtype=str)
    if "sps" not in m.columns:
        raise ValueError(f"{path}: SPS mapping lacks 'sps'; expected columns sps, "
                         f"{' or '.join(SPS_LIBRARY_COLUMNS)}")
    # NN 250 / D11a: the SHIPPED mapping (/oak/.../common/keyfiles/sps_library_mapping.tsv) names
    # this column `V_library`, so demanding `library_id` rejected the only real input there is.
    # Named explicitly, never taken positionally.
    lib_col = next((c for c in SPS_LIBRARY_COLUMNS if c in m.columns), None)
    if lib_col is None:
        raise ValueError(f"{path}: SPS mapping lacks a library column; expected one of "
                         f"{SPS_LIBRARY_COLUMNS} alongside 'sps'")

    pairs = [(str(k).upper(), v) for k, v in zip(m["sps"], m[lib_col])
             if isinstance(k, str) and isinstance(v, str) and k and v]

    # 🔴 NN 250 / D11c: FILTER BEFORE COLLAPSING. One SPS legitimately maps to several subpools
    # (V629/V545 share one SPS, V631/V546 the other -- policy §4: an SPS names the SUBPOOL, and
    # a scaffold claim can never be read out of it). The old code built a dict first, so a
    # duplicate SPS collapsed LAST-WINS and dropped V629/V631 outright; scoping to V629 then
    # raised "empty after filtering", a legitimate request failing with a misleading message.
    if v_libraries:
        keep = set(v_libraries)
        pairs = [(k, v) for k, v in pairs if v in keep]

    if not pairs:
        raise ValueError(f"{path}: SPS mapping is empty after filtering to {v_libraries}; refusing to "
                         f"run unscoped while claiming to be scoped")

    # 🔴 NN 250 / D11b: an SPS that still names more than one library IN SCOPE names neither
    # uniquely. Picking one silently is precisely the mechanism behind the 7.0% cross-subpool
    # mis-naming that SPS scoping exists to remove, so refuse and say which.
    by_sps: Dict[str, set] = {}
    for k, v in pairs:
        by_sps.setdefault(k, set()).add(v)
    ambiguous = {k: sorted(v) for k, v in by_sps.items() if len(v) > 1}
    usable = {k: next(iter(v)) for k, v in by_sps.items() if len(v) == 1}
    if ambiguous:
        # DROPPED, LOUDLY AND COUNTED -- not silently resolved, and not fatal on its own. An SPS
        # that names two libraries in scope cannot scope anything, so it is removed from the map
        # and its rows fall back to the full scope via assign_sps_scope's NaN path (which is
        # already reported per run). Refusing the whole file instead would reject the SHIPPED
        # mapping outright, since every SPS in it is shared by a 2024 and a 2025 subpool.
        logger.warning(
            "%s: %d of %d SPS map to MORE THAN ONE library within scope %s and are DROPPED from "
            "the scoping map (their rows fall back to the full scope): %s. An SPS names the "
            "design SUBPOOL only, so two libraries sharing one SPS cannot be told apart by it "
            "(policy §4). Narrow --v-libraries to scope those rows (NN 250 / D11).",
            path, len(ambiguous), len(by_sps), v_libraries or "ALL", ambiguous,
        )
    if not usable:
        raise ValueError(
            f"{path}: no SPS uniquely names a library within scope {v_libraries or 'ALL'} -- "
            f"{len(ambiguous)} of {len(by_sps)} are shared by several libraries ({ambiguous}). "
            f"Refusing to run unscoped while claiming to be scoped; narrow --v-libraries to one "
            f"of the libraries sharing each SPS (NN 250 / D11)."
        )
    return usable


def assign_sps_scope(df: pd.DataFrame, sps2v: Dict[str, str], sps_col: str = "SPS") -> pd.Series:
    """Per row, the subpool its own SPS names, or NaN when the SPS is absent or unrecognised."""
    if sps_col not in df.columns:
        return pd.Series([None] * len(df), index=df.index, dtype=object)
    return df[sps_col].astype(str).str.upper().map(sps2v)


UNSCOPED = "__UNSCOPED__"


def match_dataframe_scoped(df: pd.DataFrame, scope: pd.Series, build_tables, match_fn=None, **kw):
    """Match each group of rows against ONLY the lookup its ``scope`` names.

    ``build_tables(v)`` returns the lookup tables for subpool ``v``; it is called with ``None`` for rows
    whose scope is NaN, which must return the caller's full unscoped lookup. Row ORDER and the original
    index are restored with ``sort_index()``, so this is a drop-in for ``match_dataframe``.

    ``build_tables`` is injected rather than called directly so the grouping can be tested without the
    harmonized oligo table.
    """
    if match_fn is None:
        match_fn = match_dataframe
    if df.empty:
        return match_fn(df, build_tables(None), **kw)
    key = scope.fillna(UNSCOPED)
    parts = []
    for v, sub in df.groupby(key, sort=False):
        tables = build_tables(None if v == UNSCOPED else v)
        parts.append(match_fn(sub, tables, **kw))
    out = pd.concat(parts).sort_index()
    if len(out) != len(df):
        raise RuntimeError(f"SPS-scoped match changed row count: {len(df)} in, {len(out)} out")
    return out


# ---------------------------------------------------------------------------
# bc1 PROVENANCE AND RECURRENCE (NN 250 / D17, from NN 249's mate-pair test)
# ---------------------------------------------------------------------------
#
# 🔴 WHY NO FIELD COULD EVER HAVE NAMED THE TEMPLATE FROM SEQUENCE. NN 249 established that the
# genome-integrated bc1 cassette and the guide-donor-PLASMID bc1 cassette are THE SAME CONSTRUCT
# in opposite orientations, sequence-identical for 80 bp to the left and 2,421 bp to the right of
# the N20. With 2x150 reads, almost every read pair at the locus is genuinely uninformative:
# unambiguous genomic support requires a read reaching reference position <= 10189 or >= 12711.
# So the "matcher cannot distinguish plasmid-derived from genome-derived reads" gap is NOT an
# oversight in the matcher -- the discrimination is MATE GEOMETRY (left-reach), which lives in the
# WGS lane, not in sequence the matcher ever sees.
#
# The existing `bc1_from_wgs_source` column is a CONSTANT (`G` x 615) and therefore names nothing.
# What follows is PLUMBING ONLY: it carries NN 249's per-bc1 call and refuses to carry a constant.
# This NN does not compute left-reach.
#
# ⚠️ Also from NN 249, recorded so nobody reads the old column as evidence: the SpG run's
# "plasmid" route interrogated plasmid positions 1768-1787, which is the GUIDE SPACER (AscI +
# sgRNA-scaffold flanks). So `P_barcode == NA` across all 768 wells says NOTHING about plasmid
# bc1, and must not be cited as evidence of absence.

#: The provenance column NN 249 fills. Values are its call vocabulary, not this module's.
BC1_PROVENANCE_COLUMN = "bc1_from_wgs_source"


def attach_bc1_provenance(df, calls, *, bc1_col: str = "bc1",
                          provenance_col: str = BC1_PROVENANCE_COLUMN,
                          require_variation: bool = True, return_stats: bool = False):
    """Join NN 249's per-bc1 provenance call onto ``df``. Refuses to attach a CONSTANT column.

    Args:
        df: any bc1-keyed frame.
        calls: a frame with ``bc1`` and ``provenance_col``, or None when no call was produced.
        require_variation: raise if every called row carries the same value. 🔴 A constant column
            names nothing -- that is the exact failure this replaces, and re-creating it with a
            new name would be worse than leaving the gap open.
        return_stats: also return the join statistics.

    When ``calls`` is None, NO column is added. An always-NaN column is as empty a claim as a
    constant one, and a downstream reader cannot tell "not computed" from "computed as unknown".
    """
    stats = {"supplied": calls is not None, "n_rows": int(len(df)),
             "n_called": 0, "n_uncalled": int(len(df)), "distinct_calls": 0, "call_counts": {}}
    if calls is None or len(calls) == 0:
        return (df, stats) if return_stats else df
    if provenance_col not in calls.columns or bc1_col not in calls.columns:
        raise ValueError(
            f"provenance calls must carry '{bc1_col}' and '{provenance_col}'; got "
            f"{list(calls.columns)}"
        )

    vals = calls[provenance_col].dropna()
    distinct = sorted(set(vals.astype(str)))
    if require_variation and len(distinct) <= 1:
        raise ValueError(
            f"refusing to attach '{provenance_col}': it does not vary across the supplied calls "
            f"(every value is {distinct or ['<all null>']}). A constant column names nothing -- "
            f"the column this replaces was 'G' on all 615 rows. Supply a real per-bc1 call "
            f"(NN 250 / D17)."
        )

    m = calls.drop_duplicates(subset=[bc1_col]).set_index(bc1_col)[provenance_col]
    out = df.copy()
    out[provenance_col] = out[bc1_col].map(m)
    n_called = int(out[provenance_col].notna().sum())
    stats.update(
        n_called=n_called,
        n_uncalled=int(len(out) - n_called),
        distinct_calls=len(distinct),
        call_counts={k: int(v) for k, v in out[provenance_col].value_counts().items()},
    )
    logger.info(
        f"bc1 provenance '{provenance_col}': {n_called:,}/{len(out):,} rows carry a call "
        f"({stats['n_uncalled']:,} uncalled, counted not blanked); values {stats['call_counts']}"
    )
    return (out, stats) if return_stats else out


def flag_bc1_recurrence(df, *, bc1_col: str = "bc1", well_col: str = "well"):
    """``(df + bc1_n_wells/bc1_recurs_across_wells, stats)`` -- a bc1 seen in more than one well.

    🔴 NN 250 / NN 249-3. The decisive evidence that the residual mis-assignment is a MAPPING
    ERROR and not contamination was that 3 of 3 recurring barcodes mis-map to the SAME wrong gene
    in both clones carrying them -- contamination cannot be reproducible. Nothing in the pipeline
    surfaced that class, so it had to be found by hand. ``bc1_n_wells`` counts DISTINCT wells, so
    a barcode appearing twice in one well is not recurrence.
    """
    out = df.copy()
    nw = out.groupby(bc1_col, dropna=False, observed=True)[well_col].transform("nunique")
    out["bc1_n_wells"] = nw.astype("int64")
    out["bc1_recurs_across_wells"] = out["bc1_n_wells"] > 1
    per_bc1 = out.drop_duplicates(subset=[bc1_col])
    stats = {
        "distinct_bc1": int(out[bc1_col].nunique(dropna=False)),
        "recurring_bc1": int((per_bc1["bc1_n_wells"] > 1).sum()),
        "rows_flagged": int(out["bc1_recurs_across_wells"].sum()),
        "max_wells_per_bc1": int(out["bc1_n_wells"].max()) if len(out) else 0,
    }
    if stats["recurring_bc1"]:
        logger.warning(
            f"  {stats['recurring_bc1']:,} of {stats['distinct_bc1']:,} bc1 appear in MORE THAN "
            f"ONE well (max {stats['max_wells_per_bc1']}). A barcode that mis-maps to the SAME "
            f"wrong target in every well carrying it is a reproducible MAPPING error, not "
            f"contamination -- check those before reading any per-well assignment (NN 250 / "
            f"NN 249-3)."
        )
    return out, stats


#: Screen-side read-count floor below which a bc1 -> oligo assignment is LABELLED low-support.
#:
#: MEASURED, from NN 249: the 21 mis-assigned wells' barcodes are 13-20x more likely to be rare in
#: the screen reference (OR 13.12 [2.64, 65.13], p = 0.0019 against NE-other, which matches on
#: outcome class and so removes the CE-selection objection). The four cases it classed "rare in
#: screen" had screen-side counts of **1, 2, 2 and 7** -- "the oligo call rests on a handful of
#: reads". The floor is set at 8, the smallest integer above all four, so it captures the measured
#: class without reaching into the well-supported ones (the same table's abundant call is 28,228).
#: 🔶 It is a floor fitted to FOUR observations; treat it as provisional and re-derive it per
#: screen when a larger labelled set exists.
SCREEN_SUPPORT_FLOOR = 8


def label_low_screen_support(df, *, count_col: str = "screen_counts",
                             floor: int = None, out_col: str = "bc1_low_screen_support"):
    """LABEL (never filter) bc1 whose screen-side read support is below ``floor``.

    🔴 NN 250 / NN 249-2. Labelled rather than filtered, because this codebase's standing rule is
    "LABEL, never filter" (see ``02_combine_and_collapse``): dropping these rows would move a
    denominator with nothing to show for it. Missing support is a THIRD state (NaN), counted
    separately -- it is not evidence of low support. No count column means no column at all,
    rather than a column of False.
    """
    f = SCREEN_SUPPORT_FLOOR if floor is None else int(floor)
    stats = {"supplied": count_col in df.columns, "floor": f, "n_rows": int(len(df)),
             "n_below_floor": 0, "n_missing_support": 0, "fraction_below_floor": 0.0}
    if count_col not in df.columns:
        return df, stats
    out = df.copy()
    c = pd.to_numeric(out[count_col], errors="coerce")
    out[out_col] = (c < f).where(c.notna())
    stats["n_below_floor"] = int((c < f).sum())
    stats["n_missing_support"] = int(c.isna().sum())
    evaluable = int(c.notna().sum())
    stats["fraction_below_floor"] = (stats["n_below_floor"] / evaluable) if evaluable else 0.0
    if stats["n_below_floor"]:
        logger.warning(
            f"  {stats['n_below_floor']:,} of {evaluable:,} bc1 have screen-side support < {f} "
            f"reads and are LABELLED low-support (not dropped). NN 249 measured such barcodes "
            f"13-20x over-represented among mis-assigned wells (OR 13.12, p=0.0019) "
            f"(NN 250 / NN 249-2)."
        )
    return out, stats
