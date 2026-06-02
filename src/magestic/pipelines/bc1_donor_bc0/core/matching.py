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
from typing import Dict, List, Optional, Set, Tuple
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

# Anchor used to pre-trim observed bc1_donor_bc0 reads down to the variable
# donor region. The observed donor extracted at the parse step carries a
# constant 5' flank (SaCas9/SceI fragment + scaffold leader) terminating in the
# SpCas9 cloning site `GTTTGAAGAGC`; everything AFTER that anchor is the
# designed donor (~129bp / ~117bp). Trimming up to and including this anchor
# lets the perfect-donor lookup fire and removes the flanking-driven ambiguity.
DONOR_VARIABLE_REGION_ANCHOR = "GTTTGAAGAGC"

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
    anchor: str = DONOR_VARIABLE_REGION_ANCHOR,
    expected_pos: int = DEFAULT_ANCHOR_POS,
    pos_slack: int = DEFAULT_ANCHOR_POS_SLACK,
    max_mismatch: int = DEFAULT_ANCHOR_MAX_MISMATCH,
) -> Optional[str]:
    """
    Trim an observed donor down to the variable region following the anchor.

    Tolerant + anchor-conditional + idempotent:
    1. None/NaN: returned unchanged.
    2. Exact match anywhere in the donor: trim at first occurrence (fast path,
       covers ~88% of yL406-like reads).
    3. No exact match: search a Hamming-up-to-`max_mismatch` window in
       [expected_pos - pos_slack, expected_pos + pos_slack]. This rescues the
       ~12% of reads where a single-base sequencer Q-drop (most commonly
       N at idx 0 or idx 10) destroyed the exact match. Ties broken by
       closeness to `expected_pos`.
    4. Still no match: return donor unchanged (already-trimmed donors, clean
       gdb donors, reads with corrupted flank all flow to sliding-window
       fallback).

    Args:
        donor: Observed donor sequence (may be None).
        anchor: Constant 3' boundary of the flank (default SpCas9 cloning site).
        expected_pos: Canonical 0-indexed start of the anchor in the donor
            field (98 for MAGESTIC 3.0 bc1_donor_bc0 amplicons).
        pos_slack: How far either side of expected_pos to scan for the
            tolerant match.
        max_mismatch: Hamming budget for the tolerant search. Set to 0 to
            disable the tolerant path entirely and reproduce the prior
            exact-only behavior.

    Returns:
        The trimmed donor, or the original donor if no anchor is found.
    """
    if not donor:
        return donor
    L = len(anchor)
    idx = donor.find(anchor)
    if idx >= 0:
        return donor[idx + L:]
    if max_mismatch <= 0:
        return donor
    pos = _find_anchor_tolerant(donor, anchor, expected_pos, pos_slack, max_mismatch)
    if pos < 0:
        return donor
    return donor[pos + L:]


@dataclass
class OligoMatch:
    """Result of matching a sequence to designed oligos."""
    oligo_name: Optional[str] = None

    # Match quality indicators
    perfect_donor: bool = False
    perfect_middle_region: bool = False
    guide_match: bool = False

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

    # 1. Try perfect donor match
    if donor and donor in tables.perfect_donor_to_oligos:
        candidates = tables.perfect_donor_to_oligos[donor]
        result.perfect_donor = True
        result.perfect_middle_region = True

        # Narrow by guide first if it disambiguates further (rarely available).
        if guide and guide in tables.guide_to_oligos:
            overlap = set(candidates) & set(tables.guide_to_oligos[guide])
            if overlap:
                result.guide_match = True
                candidates = list(overlap)

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

            if guide and guide in tables.guide_to_oligos:
                overlap = window_matches & set(tables.guide_to_oligos[guide])
                if overlap:
                    result.guide_match = True
                    window_matches = overlap

            _resolve_candidates(
                result, window_matches, tables,
                method_unique="sliding_window",
                method_ambiguous="sliding_window_ambiguous",
            )
            return result

    # 3. Guide-only match (fallback)
    if guide and guide in tables.guide_to_oligos:
        candidates = tables.guide_to_oligos[guide]
        result.guide_match = True
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
    - {prefix}guide_match: True if guide contributed to match
    - {prefix}is_unambiguous: True if match is unambiguous
    - {prefix}match_method: Description of how match was made

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
    from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
    from tqdm import tqdm

    logger.info(f"Matching {len(df)} sequences to oligos (n_jobs={n_jobs})...")

    # Initialize result columns
    df = df.copy()
    df[f'{prefix}oligo_name'] = None
    df[f'{prefix}perfect_donor'] = False
    df[f'{prefix}perfect_middle_region'] = False
    df[f'{prefix}guide_match'] = False
    df[f'{prefix}is_unambiguous'] = False
    df[f'{prefix}match_method'] = "none"
    df[f'{prefix}num_candidates'] = 0
    df[f'{prefix}num_candidate_variants'] = 0
    df[f'{prefix}variant'] = None

    # Check if guide column exists
    has_guide = guide_col in df.columns

    # Prepare data for matching
    donors = df[donor_col].values
    guides = df[guide_col].values if has_guide else [None] * len(df)

    # Match using vectorized approach for lookups where possible
    match_counts = {"perfect": 0, "sliding_window": 0, "guide_only": 0, "none": 0}

    # Process in chunks to manage memory and show progress
    results = []
    total_rows = len(df)

    for start_idx in tqdm(range(0, total_rows, chunk_size), desc="Matching chunks"):
        end_idx = min(start_idx + chunk_size, total_rows)
        chunk_donors = donors[start_idx:end_idx]
        chunk_guides = guides[start_idx:end_idx]

        if n_jobs > 1:
            # Parallel processing within chunk
            # Note: ThreadPoolExecutor is better here because tables is shared
            with ThreadPoolExecutor(max_workers=n_jobs) as executor:
                chunk_results = list(executor.map(
                    lambda args: match_sequence(
                        args[0], args[1] if pd.notna(args[1]) else None, tables,
                        pre_trim=pre_trim,
                    ) if pd.notna(args[0]) else None,
                    zip(chunk_donors, chunk_guides)
                ))
        else:
            # Sequential processing
            chunk_results = [
                match_sequence(d, g if pd.notna(g) else None, tables, pre_trim=pre_trim)
                if pd.notna(d) else None
                for d, g in zip(chunk_donors, chunk_guides)
            ]

        results.extend(chunk_results)

    # Apply results to DataFrame
    for idx, result in enumerate(results):
        if result is None:
            continue

        df.iloc[idx, df.columns.get_loc(f'{prefix}oligo_name')] = result.oligo_name
        df.iloc[idx, df.columns.get_loc(f'{prefix}perfect_donor')] = result.perfect_donor
        df.iloc[idx, df.columns.get_loc(f'{prefix}perfect_middle_region')] = result.perfect_middle_region
        df.iloc[idx, df.columns.get_loc(f'{prefix}guide_match')] = result.guide_match
        df.iloc[idx, df.columns.get_loc(f'{prefix}is_unambiguous')] = result.is_unambiguous
        df.iloc[idx, df.columns.get_loc(f'{prefix}match_method')] = result.match_method
        df.iloc[idx, df.columns.get_loc(f'{prefix}num_candidates')] = result.num_candidate_oligos
        df.iloc[idx, df.columns.get_loc(f'{prefix}num_candidate_variants')] = result.num_candidate_variants
        df.iloc[idx, df.columns.get_loc(f'{prefix}variant')] = result.variant

        # Count match types
        if result.perfect_donor:
            match_counts["perfect"] += 1
        elif result.perfect_middle_region:
            match_counts["sliding_window"] += 1
        elif result.guide_match:
            match_counts["guide_only"] += 1
        else:
            match_counts["none"] += 1

    logger.info(f"  Match results:")
    logger.info(f"    Perfect donor: {match_counts['perfect']}")
    logger.info(f"    Sliding window: {match_counts['sliding_window']}")
    logger.info(f"    Guide only: {match_counts['guide_only']}")
    logger.info(f"    No match: {match_counts['none']}")

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
