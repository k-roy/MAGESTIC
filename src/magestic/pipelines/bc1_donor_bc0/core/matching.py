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

# ---------------------------------------------------------------------------
# Edit certification (planning 218, ledger item 10; BDB-H05 / DONOR-WINDOW-R1)
# ---------------------------------------------------------------------------
# A window hit says a 60-mer of the read equals a 60-mer of a design. It does NOT say the read carries the
# edit that design installs: for a 129-nt donor a window avoiding any single offset always exists, so the
# UNEDITED donor of a design still hits that design. 98 real registry rows keep oligo, perfect_middle_region
# and closest-middle ED0 after their annotated SNVs are reverted to reference. Certification therefore asks a
# separate question -- does the read carry this design's intended allele at every intended position, inside an
# exactly matched window? -- and answers with one of these statuses. No new threshold: the only length used is
# the lookup's existing window_size.
EDIT_FULL_DONOR = "full_donor"
EDIT_POSITIONS_MATCHED = "edit_positions_matched"
EDIT_POSITIONS_MISMATCHED = "edit_positions_mismatched"
EDIT_POSITIONS_NOT_COVERED = "edit_positions_not_covered"
EDIT_POSITIONS_UNCONFIRMED = "edit_positions_unconfirmed"
EDIT_NOT_ASSESSABLE = "not_assessable"
EDIT_NO_DONOR_EVIDENCE = "no_donor_evidence"
# Ordered by how much the row tells you about the edit, certifying first; used to summarise a candidate SET.
EDIT_EVIDENCE_STATUSES = (
    EDIT_FULL_DONOR, EDIT_POSITIONS_MATCHED, EDIT_POSITIONS_MISMATCHED, EDIT_POSITIONS_NOT_COVERED,
    EDIT_POSITIONS_UNCONFIRMED, EDIT_NOT_ASSESSABLE, EDIT_NO_DONOR_EVIDENCE,
)
EDIT_CERTIFYING = frozenset({EDIT_FULL_DONOR, EDIT_POSITIONS_MATCHED})

# Registry columns that place a designed donor on the genome. The convention is INFERRED from the 2026-09
# geometry scan of the shipped registry and validated per row in design_edit_positions: the donor occupies the
# zero-based half-open interval [donor_start_coord, donor_end_coord) on the forward strand and variant
# coordinates are one-based. Of the four conventions tested, that one passed 94,586 of 146,250 rows and the
# zero-based-variant alternative never passed alone (216/preflight/donor_geometry/summaries.json).
# 35.3% of shipped rows lack donor_strand and so stay not_assessable until the registry generator backfills
# it (Kevin, decision D-1, 2026-09-23: accept fail-closed, then backfill).
DONOR_GEOMETRY_COLUMNS = ("donor_start_coord", "donor_end_coord", "donor_strand")
EDIT_POSITION_VARIANT_COLUMN = "individual_variants"
_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_SNV_TOKEN = re.compile(r"(?:chr[^_]+_)?(\d+)_([ACGT])_([ACGT])(?:_(?:SNP|SNV))?")


def _hamming(a: str, b: str) -> int:
    """Count mismatches between two equal-length strings. N counts as mismatch."""
    return sum(x != y for x, y in zip(a, b))


def _coord(value) -> Optional[int]:
    """An integer coordinate, or None. Accepts the registry's float spelling ('344980.0')."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if number != number or number != int(number):      # NaN, or not an integer
        return None
    return int(number)


def design_edit_positions(donor, individual_variants, donor_start, donor_end, donor_strand):
    """0-based offsets, in the DESIGNED donor, of every allele this design intends to install.

    Returns ``(offsets, None)`` or ``(None, reason)``. The row is accepted only when the designed donor itself
    carries the alternate allele (complemented on '-') at every projected offset, so a wrong convention, a
    gapped design or a stale annotation fails closed instead of pointing at the wrong base. Only single-base
    substitutions are projected: multi-base and indel alleles need a gapped transform this registry does not
    carry, and they are reported not assessable rather than guessed at.
    """
    if not isinstance(donor, str) or not donor:
        return None, "no_donor"
    if donor_strand not in ("+", "-"):
        return None, "no_donor_strand"
    start, end = _coord(donor_start), _coord(donor_end)
    if start is None or end is None:
        return None, "no_donor_coordinates"
    if end - start != len(donor):
        return None, "interval_length_disagrees_with_donor"
    if not isinstance(individual_variants, str) or not individual_variants.strip():
        return None, "no_individual_variants"
    offsets = []
    for raw in re.split(r"[,;]", individual_variants):
        token = raw.strip().strip("[]").strip().strip("'\"")
        parsed = _SNV_TOKEN.fullmatch(token)
        if not parsed:
            return None, "not_a_single_base_substitution"
        position, ref, alt = int(parsed.group(1)), parsed.group(2), parsed.group(3)
        if ref == alt:
            return None, "ref_equals_alt"
        offset = (position - 1 - start) if donor_strand == "+" else (end - position)
        expected = alt if donor_strand == "+" else alt.translate(_COMPLEMENT)
        if not (0 <= offset < len(donor)) or donor[offset] != expected:
            return None, "alternate_allele_not_at_projected_offset"
        offsets.append(offset)
    if len(set(offsets)) != len(offsets):
        return None, "two_alleles_at_one_offset"
    return tuple(sorted(offsets)), None


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
    # True only when the EDIT-bearing region matched: edit_evidence in EDIT_CERTIFYING (planning 218).
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

    # Confidence
    #
    # `is_unambiguous` is VARIANT-level, not oligo_name-level. Many designed
    # oligos share an identical donor (the HDR repair template) but differ in
    # which guide cuts -- they install the SAME edit. Collapsing candidate
    # oligos to their `donor_changes` (the installed variant) shows a single
    # variant in ~100% of perfect-donor matches even though the oligo count is
    # 4-6. A match is therefore "unambiguous" when all candidate oligos agree
    # on the installed variant, regardless of how many redundant guides exist.
    # Since planning 218 (ledger item 10) it ALSO requires certifying edit evidence: a window-only
    # match that cannot show the design's intended allele stays a candidate.
    is_unambiguous: bool = False
    num_candidate_oligos: int = 0     # distinct candidate oligo_names
    num_candidate_variants: int = 0   # distinct installed variants among candidates (the OLD is_unambiguous)
    variant: Optional[str] = None     # installed edit when unambiguous AND certified

    # What the observed donor SHOWS about the candidates' designed edit (planning 218, ledger item 10).
    # `edit_evidence` summarises the candidate SET -- the strongest status any candidate reaches -- and
    # certifies `oligo_name` only when `is_unambiguous` is True, because an ambiguous row's `oligo_name` is
    # merely the sorted-first candidate. A window hit alone is a CANDIDATE, never a certification.
    edit_evidence: str = EDIT_NO_DONOR_EVIDENCE
    num_edit_certified_variants: int = 0

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

    # oligo_name -> its designed donor, and -> the 0-based donor offsets of the alleles it installs
    # (None when the design's geometry cannot be projected; see design_edit_positions). Both are needed to
    # ask whether a read actually shows a candidate's edit.
    oligo_to_donor: Dict[str, Optional[str]] = field(default_factory=dict)
    oligo_to_edit_positions: Dict[str, Optional[Tuple[int, ...]]] = field(default_factory=dict)

    # Metadata
    window_size: int = DEFAULT_WINDOW_SIZE
    num_oligos: int = 0

    def get_stats(self) -> Dict:
        """Get statistics about the lookup tables."""
        return {
            "num_oligos": self.num_oligos,
            "num_unique_donors": len(self.perfect_donor_to_oligos),
            "num_unique_guides": len(self.guide_to_oligos),
            "num_unique_windows": len(self.window_to_oligos),
            "window_size": self.window_size,
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


def build_lookup_tables(
    oligo_df: pd.DataFrame,
    window_size: int = DEFAULT_WINDOW_SIZE,
    variant_col: str = VARIANT_COLUMN,
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

    tables = OligoLookupTables(window_size=window_size)

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
    tables.guide_donor_to_oligo = dict(
        zip(valid_df['guide_donor_key'], valid_df['oligo_name'])
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
        variant_series = (
            valid_df[variant_col]
            .where(valid_df[variant_col].notna(), valid_df['donor'])
            .astype(str)
        )
        n_fallback = valid_df[variant_col].isna().sum()
        if n_fallback:
            logger.info(
                f"  Variant id '{variant_col}' missing for {n_fallback:,} "
                f"({n_fallback/len(valid_df):.1%}) oligos; using donor sequence "
                f"as the per-row fallback identity for those."
            )
    else:
        logger.info(
            f"  Variant column '{variant_col}' not found; using donor sequence "
            f"as the variant identity for unambiguity assessment."
        )
        variant_series = valid_df['donor']
    tables.oligo_to_variant = dict(zip(valid_df['oligo_name'], variant_series))

    # 6. Edit geometry: the donor and the donor-coordinate offsets of the alleles each design installs.
    # One extra pass over the designs (no extra pass over windows). A design whose geometry cannot be
    # projected is not assessable: such rows can only be certified by the full donor being present.
    conflicting = set()
    for name, donor in zip(valid_df['oligo_name'], valid_df['donor']):
        previous = tables.oligo_to_donor.get(name)
        if previous is not None and previous != donor:
            conflicting.add(name)
        tables.oligo_to_donor[name] = donor
    for name in conflicting:
        # One name, two donors: the dict keeps only the last, so neither the full-donor test nor a projected
        # position can be trusted for it. Blanking the donor makes design_edit_evidence return not_assessable
        # (it guards `if not design`). load_and_build_lookup dedupes by name, so production never gets here;
        # without this line a duplicated name could certify full_donor against the wrong donor.
        tables.oligo_to_donor[name] = None
    reasons = {}
    have_geometry = all(c in valid_df.columns for c in DONOR_GEOMETRY_COLUMNS)
    if have_geometry and EDIT_POSITION_VARIANT_COLUMN in valid_df.columns:
        cols = ['oligo_name', 'donor', EDIT_POSITION_VARIANT_COLUMN, *DONOR_GEOMETRY_COLUMNS]
        for name, donor, variants, start, end, strand in valid_df[cols].itertuples(index=False, name=None):
            positions, reason = design_edit_positions(donor, variants, start, end, strand)
            if name in conflicting:
                positions, reason = None, "oligo_name_has_conflicting_donors"
            tables.oligo_to_edit_positions[name] = positions
            reasons[reason or "assessable"] = reasons.get(reason or "assessable", 0) + 1
        logger.info(f"  Edit positions by design: {reasons}")
    else:
        missing = [c for c in (*DONOR_GEOMETRY_COLUMNS, EDIT_POSITION_VARIANT_COLUMN)
                   if c not in valid_df.columns]
        tables.oligo_to_edit_positions = {name: None for name in tables.oligo_to_donor}
        logger.warning(
            f"  Design table lacks {missing}: no design's intended alleles can be placed in donor "
            f"coordinates, so a window-only match can never be certified (edit_evidence=not_assessable). "
            f"Only reads containing the FULL designed donor will be certified.")

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
    """
    candidates = sorted(set(candidates))
    variants = _variant_set(candidates, tables)
    result.candidate_oligos = candidates
    result.num_candidate_oligos = len(candidates)
    result.num_candidate_variants = len(variants)
    result.oligo_name = candidates[0]

    if len(variants) == 1:
        result.is_unambiguous = True
        result.variant = next(iter(variants))
        result.match_method = method_unique
    else:
        result.match_method = method_ambiguous


def _observed_windows(donor: str, window_size: int) -> Dict[str, List[int]]:
    """Every window of the OBSERVED donor -> the offsets it starts at. Built once per read."""
    windows: Dict[str, List[int]] = {}
    for i in range(len(donor) - window_size + 1):
        windows.setdefault(donor[i:i + window_size], []).append(i)
    return windows


def design_edit_evidence(observed, observed_windows, design, positions, window_size) -> str:
    """What the observed donor shows about ONE design's intended alleles (the EDIT_* statuses).

    An intended position counts as SEEN only inside a window of the design that occurs exactly in the read;
    inside such a window the read's base at that position is the design's allele by construction. Positions
    outside every exact window are classified from the diagonals of the matched windows, which is diagnosis,
    not certification.
    """
    if not design:
        return EDIT_NOT_ASSESSABLE
    if design in observed:
        return EDIT_FULL_DONOR
    if not positions:
        return EDIT_NOT_ASSESSABLE
    starts, diagonals = set(), set()
    for j in range(len(design) - window_size + 1):
        hits = observed_windows.get(design[j:j + window_size])
        if hits:
            starts.add(j)
            diagonals.update(i - j for i in hits)
    uncovered = [p for p in positions if not any((p - k) in starts for k in range(window_size))]
    if not uncovered:
        return EDIT_POSITIONS_MATCHED
    mismatched = truncated = False
    for p in uncovered:
        spanning = [d for d in diagonals if 0 <= p + d < len(observed)]
        if not spanning:
            truncated = True
        elif all(observed[p + d] != design[p] for d in spanning):
            mismatched = True
    if mismatched:
        return EDIT_POSITIONS_MISMATCHED
    if truncated:
        return EDIT_POSITIONS_NOT_COVERED
    return EDIT_POSITIONS_UNCONFIRMED


def _certify_candidates(result: OligoMatch, observed: str, tables: OligoLookupTables) -> None:
    """Replace the window branch's assumed certification with what the read actually shows.

    Candidates, counts, `oligo_name` and `match_method` are untouched -- a window-only match remains a
    CANDIDATE. Only `perfect_middle_region`, `is_unambiguous` and `variant` change, and only downwards.
    """
    windows = _observed_windows(observed, tables.window_size)
    cache, statuses = {}, {}
    for name in result.candidate_oligos:
        design = tables.oligo_to_donor.get(name)
        positions = tables.oligo_to_edit_positions.get(name)
        key = (design, positions)
        if key not in cache:
            cache[key] = design_edit_evidence(observed, windows, design, positions, tables.window_size)
        statuses[name] = cache[key]
    certified = [n for n, s in statuses.items() if s in EDIT_CERTIFYING]
    result.num_edit_certified_variants = len(_variant_set(certified, tables))
    if statuses:
        result.edit_evidence = min(statuses.values(), key=EDIT_EVIDENCE_STATUSES.index)
    result.perfect_middle_region = result.edit_evidence in EDIT_CERTIFYING
    if not result.perfect_middle_region:
        result.is_unambiguous = False
        result.variant = None


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

    # Normalize to uppercase
    if donor:
        donor = donor.upper()
        if pre_trim:
            donor = trim_to_variable_region(donor)
    if guide:
        guide = guide.upper()

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
        result.edit_evidence = EDIT_FULL_DONOR
        result.num_edit_certified_variants = result.num_candidate_variants
        return result

    # 2. Try sliding window match on donor
    if donor:
        window_matches = match_donor_sliding_window(donor, tables)

        if window_matches:
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
            # A matched 60-mer names CANDIDATES. Whether the read shows their edit is a separate question.
            _certify_candidates(result, donor, tables)
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
        # A designed guide names a construct; it is no evidence that this plasmid's donor carries that
        # construct's edit -- the donor on this row matched nothing. Certification needs donor evidence
        # (reviews/bdb_attribution_contract/PLAN_AND_REVIEW.md, `variant_candidates` row). Decision D-2,
        # implemented as recommended; reversing it means deleting these three lines and test_guideonly_*.
        result.edit_evidence = EDIT_NO_DONOR_EVIDENCE
        result.is_unambiguous = False
        result.variant = None
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
        m[donor] = (
            cands[0],
            nv == 1,
            (next(iter(variants)) if nv == 1 else None),
            len(cands),
            nv,
            "perfect" if nv == 1 else "perfect_ambiguous",
            format_candidates(cands, nv),
        )
    try:
        tables._perfect_result_map = m
    except Exception:
        pass
    return m


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
    - {prefix}perfect_middle_region: True if the EDIT-bearing region matched (edit_evidence certifies)
    - {prefix}guide_match: True if guide contributed to match (kept for back-compat)
    - {prefix}guide_status: no_guide / guide_not_designed / match / conflict -- see OligoMatch
    - {prefix}num_guide_oligos: how many designs carry this guide (0 when it is not a designed guide)
    - {prefix}is_unambiguous: True if match is unambiguous and its edit is certified
    - {prefix}match_method: Description of how match was made
    - {prefix}num_candidates / {prefix}num_candidate_variants: candidate oligos / distinct installed variants
    - {prefix}variant: installed variant id when unambiguous and certified
    - {prefix}candidate_oligos: ';'-joined candidates when num_candidate_variants > 1, else None
    - {prefix}edit_evidence: what the donor shows about the candidates' designed edit (EDIT_* statuses)
    - {prefix}num_edit_certified_variants: distinct installed variants whose edit the donor certifies

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
    editev = [EDIT_NO_DONOR_EVIDENCE] * n
    ncertv = [0] * n

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
        upu = donor_s.astype(str).str.upper().where(donor_s.notna())
        a1, a2 = DONOR_VARIABLE_REGION_ANCHORS[0], DONOR_VARIABLE_REGION_ANCHORS[1]
        t1 = upu.str.split(a1, n=1).str[1]   # part after first anchor1, else NaN
        t2 = upu.str.split(a2, n=1).str[1]   # part after first anchor2, else NaN
        trimmed = t1.where(t1.notna(), t2)   # anchor1 priority (matches trim())
        mapped = trimmed.map(pm)             # result-tuple or NaN
        if has_guide:
            guide_s = df[guide_col]
            g_present = guide_s.notna().to_numpy(dtype=bool)
            g_block = guide_s.astype(str).str.upper().where(guide_s.notna()).map(
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
            # the fast path is only taken where the guide is NOT a designed guide (~g_block), so the
            # guide can never have narrowed here: it is absent, or present and undesigned.
            gstatus[i] = GUIDE_NOT_DESIGNED if g_present[i] else GUIDE_NO_GUIDE
            editev[i] = EDIT_FULL_DONOR          # the fast path is exact full-donor equality
            ncertv[i] = int(t[4])                # every candidate of a perfect donor is certified

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
        editev[i] = r.edit_evidence; ncertv[i] = r.num_edit_certified_variants

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
    df[f'{prefix}edit_evidence'] = editev
    df[f'{prefix}num_edit_certified_variants'] = ncertv

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
            f"(guide_match=False meant this and 'no guide' alike)."
        )
    evidence_counts = {}
    for status in editev:
        evidence_counts[status] = evidence_counts.get(status, 0) + 1
    logger.info(f"  edit_evidence: {evidence_counts}")
    uncertified = sum(1 for status, how in zip(editev, method)
                      if how.startswith("sliding_window") and status not in EDIT_CERTIFYING)
    if uncertified:
        logger.warning(
            f"  {uncertified:,} rows matched a designed 60-mer window but do NOT show that design's intended "
            f"allele(s). They keep their candidate names and are NOT certified (is_unambiguous=False, "
            f"perfect_middle_region=False). Before this change they were reported perfect_middle_region=True "
            f"and could be unambiguous (ledger item 10).")
    return df


def load_and_build_lookup(
    oligo_files: Optional[List[Path]] = None,
    window_size: int = DEFAULT_WINDOW_SIZE,
    use_harmonized_table: bool = True,
    v_libraries: Optional[List[str]] = None,
    nuclease_types: Optional[List[str]] = None,
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

    return build_lookup_tables(combined_df, window_size=window_size)


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


def load_sps_mapping(path: str, v_libraries: Optional[Sequence[str]] = None) -> Dict[str, str]:
    """Read an ``sps -> library_id`` TSV into an uppercased dict, optionally restricted to a scope.

    Raises rather than silently returning an empty map: a mistyped path that yields ``{}`` would turn
    SPS scoping off without saying so, and the run would look like it had scoped when it had not.
    """
    m = pd.read_csv(path, sep="\t", dtype=str)
    missing = [c for c in ("sps", "library_id") if c not in m.columns]
    if missing:
        raise ValueError(f"{path}: SPS mapping lacks {missing}; expected columns sps, library_id")
    out = {str(k).upper(): v for k, v in zip(m["sps"], m["library_id"])
           if isinstance(k, str) and isinstance(v, str) and k and v}
    if v_libraries:
        keep = set(v_libraries)
        out = {k: v for k, v in out.items() if v in keep}
    if not out:
        raise ValueError(f"{path}: SPS mapping is empty after filtering to {v_libraries}; refusing to "
                         f"run unscoped while claiming to be scoped")
    return out


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
