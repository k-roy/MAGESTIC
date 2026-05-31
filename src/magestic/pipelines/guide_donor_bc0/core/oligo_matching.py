"""
Oligo matching functions for the Guide-Donor-BC0 pipeline.

Provides functions to match guide-donor sequences to designed oligos
using exact and partial matching strategies.
"""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd

from ..config import (
    GUIDE_START, GUIDE_END,
    DONOR_START, DONOR_END,
    MIDDLE_DONOR_START, MIDDLE_DONOR_END,
    MIDDLE_DONOR_LENGTH,
    SPCAS9_CLONING_SITE,
    PERFECT_LEADER,
)

logger = logging.getLogger(__name__)


@dataclass
class OligoMatchResult:
    """Result of oligo matching."""
    oligo_name: Optional[str]
    match_type: str  # 'exact', 'guide_only', 'donor_only', 'partial', 'unmatched'
    confidence: float  # 0.0 to 1.0


def load_oligo_pool(
    filename: Path,
    guide_slice: Tuple[int, int] = (GUIDE_START, GUIDE_END),
    donor_slice: Tuple[int, int] = (DONOR_START, DONOR_END),
) -> pd.DataFrame:
    """
    Load and process oligo pool design file.

    Args:
        filename: Path to oligo design TSV file
        guide_slice: (start, end) positions for guide in oligo_seq
        donor_slice: (start, end) positions for donor in oligo_seq

    Returns:
        DataFrame with added guide, donor, and middle_donor columns
    """
    df = pd.read_csv(filename, sep="\t")

    # Extract guide and donor sequences
    df["guide"] = df["oligo_seq"].str.slice(*guide_slice)
    df["donor"] = df["oligo_seq"].str.slice(*donor_slice).str.upper()

    # Extract middle donor region for partial matching
    middle_start = donor_slice[0] + MIDDLE_DONOR_START
    middle_end = donor_slice[0] + MIDDLE_DONOR_END
    df["middle_donor"] = df["oligo_seq"].str.slice(middle_start, middle_end).str.upper()

    return df


def build_lookup_dictionaries(oligo_pool_df: pd.DataFrame) -> Dict[str, Dict]:
    """
    Build lookup dictionaries for matching sequences to oligos.

    Uses vectorized pandas operations for improved performance (~10-50x faster
    than iterrows).

    Args:
        oligo_pool_df: DataFrame from load_oligo_pool()

    Returns:
        Dictionary containing:
        - guide_to_oligos: guide -> [oligo_names]
        - perfect_donor_to_oligos: donor -> [oligo_names]
        - partial_donor_to_oligos: partial_donor -> [oligo_names]
        - guide_donor_to_oligo: guide+cloning_site+donor -> oligo_name
    """
    # 1. Exact guide+donor match (vectorized)
    oligo_pool_df = oligo_pool_df.copy()
    oligo_pool_df["guide_donor"] = (
        oligo_pool_df["guide"] + SPCAS9_CLONING_SITE + oligo_pool_df["donor"]
    )
    guide_donor_to_oligo = dict(
        zip(oligo_pool_df["guide_donor"], oligo_pool_df["oligo_name"])
    )

    # 2. Guide -> oligos (vectorized groupby)
    guide_to_oligos = (
        oligo_pool_df.groupby("guide")["oligo_name"]
        .apply(lambda x: x.drop_duplicates().tolist())
        .to_dict()
    )

    # 3. Donor -> oligos (vectorized groupby)
    perfect_donor_to_oligos = (
        oligo_pool_df.groupby("donor")["oligo_name"]
        .apply(lambda x: x.drop_duplicates().tolist())
        .to_dict()
    )

    # 4. Middle donor -> oligos (vectorized groupby)
    partial_donor_to_oligos = (
        oligo_pool_df.groupby("middle_donor")["oligo_name"]
        .apply(lambda x: x.drop_duplicates().tolist())
        .to_dict()
    )

    return {
        "guide_to_oligos": guide_to_oligos,
        "perfect_donor_to_oligos": perfect_donor_to_oligos,
        "partial_donor_to_oligos": partial_donor_to_oligos,
        "guide_donor_to_oligo": guide_donor_to_oligo,
    }


def resolve_oligo_name(
    guide: str,
    donor: str,
    lookups: Dict[str, Dict],
    middle_donor_start: int = MIDDLE_DONOR_START,
    middle_donor_end: int = MIDDLE_DONOR_END,
) -> OligoMatchResult:
    """
    Resolve oligo name from guide and donor sequences.

    Matching priority:
    1. Exact guide+donor match (highest confidence)
    2. Guide-only match (if unique)
    3. Donor-only match (if unique)
    4. Partial donor match (if unique)

    Args:
        guide: Guide sequence (20bp)
        donor: Donor sequence (~129bp)
        lookups: Dictionary from build_lookup_dictionaries()
        middle_donor_start: Start position for partial donor matching
        middle_donor_end: End position for partial donor matching

    Returns:
        OligoMatchResult with oligo_name, match_type, and confidence
    """
    donor = donor.upper()
    guide_donor = guide + SPCAS9_CLONING_SITE + donor

    # Priority 1: Exact guide+donor match
    if guide_donor in lookups["guide_donor_to_oligo"]:
        return OligoMatchResult(
            oligo_name=lookups["guide_donor_to_oligo"][guide_donor],
            match_type="exact",
            confidence=1.0,
        )

    # Priority 2: Check guide matches
    guide_matches = lookups["guide_to_oligos"].get(guide, [])

    # Priority 3: Check perfect donor matches
    donor_matches = lookups["perfect_donor_to_oligos"].get(donor, [])

    # If both match uniquely to same oligo
    if len(guide_matches) == 1 and len(donor_matches) == 1:
        if guide_matches[0] == donor_matches[0]:
            return OligoMatchResult(
                oligo_name=guide_matches[0],
                match_type="exact",
                confidence=0.95,
            )

    # If guide matches uniquely
    if len(guide_matches) == 1:
        return OligoMatchResult(
            oligo_name=guide_matches[0],
            match_type="guide_only",
            confidence=0.7,
        )

    # If donor matches uniquely
    if len(donor_matches) == 1:
        return OligoMatchResult(
            oligo_name=donor_matches[0],
            match_type="donor_only",
            confidence=0.7,
        )

    # Priority 4: Try partial donor matching
    middle_donor = donor[middle_donor_start:middle_donor_end]
    partial_matches = lookups["partial_donor_to_oligos"].get(middle_donor, [])

    if len(partial_matches) == 1:
        return OligoMatchResult(
            oligo_name=partial_matches[0],
            match_type="partial",
            confidence=0.5,
        )

    # No unique match found
    return OligoMatchResult(
        oligo_name=None,
        match_type="unmatched",
        confidence=0.0,
    )


def process_counts_file(
    input_path: Path,
    output_path: Path,
    lookups: Dict[str, Dict],
    chunk_size: int = 100_000,
) -> Dict[str, int]:
    """
    Match every row in a Step-02 counts file to the designed oligo pool.

    Reads *input_path* (produced by :func:`fastq_processing.process_library`)
    in chunks, calls :func:`resolve_oligo_name` for each row, and writes the
    result to *output_path*.

    Parameters
    ----------
    input_path : Path
        ``{library_ID}_guide_donor_bc0_counts.tsv`` from Step 02.
        Expected columns: ``library_ID, total_counts, ..., leader, guide,
        donor, SPS, bc0, msd, donor_bc0_fragment, seq``
    output_path : Path
        Destination for the matched file
        (``{library_ID}_matched_guide_donor_bc0_counts.tsv``).
    lookups : dict
        Output of :func:`build_lookup_dictionaries`.
    chunk_size : int
        Rows per pandas chunk (default 100 000).

    Returns
    -------
    dict
        ``total``, ``matched``, ``perfect_match`` counts.
    """
    stats: Dict[str, int] = {"total": 0, "matched": 0, "perfect_match": 0}
    first_chunk = True

    for chunk in pd.read_csv(input_path, sep="\t", chunksize=chunk_size):
        stats["total"] += len(chunk)
        rows = []
        for _, row in chunk.iterrows():
            guide = row.get("guide", "NA")
            donor = row.get("donor", "NA")
            leader = row.get("leader", "NA")
            bc0 = row.get("bc0", "NA")
            donor_bc0_fragment = row.get("donor_bc0_fragment", "NA")
            library_id = row.get("library_ID", "NA")
            counts = row.get("total_counts", 0)

            if pd.isna(guide) or guide == "NA":
                rows.append({
                    "library_ID": library_id,
                    "scaffold": "NA",
                    "bc0": bc0,
                    "counts": counts,
                    "perfect_leader": False,
                    "leader": leader,
                    "perfect_guide": False,
                    "guide": guide,
                    "perfect_donor": False,
                    "perfect_middle_donor_region": False,
                    "donor": donor,
                    "donor_bc0_fragment": donor_bc0_fragment,
                    "unambiguous_guide": False,
                    "unambiguous_donor": False,
                    "oligo_name": "NA",
                    "match_type": "unmatched",
                })
                continue

            result = resolve_oligo_name(
                str(guide), str(donor) if not pd.isna(donor) else "", lookups
            )
            guide_matches = lookups["guide_to_oligos"].get(str(guide), [])
            donor_key = str(donor).upper() if not pd.isna(donor) else ""
            donor_matches = lookups["perfect_donor_to_oligos"].get(donor_key, [])
            middle_donor = donor_key[MIDDLE_DONOR_START:MIDDLE_DONOR_END]
            middle_matches = lookups["partial_donor_to_oligos"].get(middle_donor, [])

            matched = result.oligo_name is not None
            if matched:
                stats["matched"] += 1
            if result.match_type == "exact":
                stats["perfect_match"] += 1

            rows.append({
                "library_ID": library_id,
                "scaffold": "WT_SpCas9",
                "bc0": bc0,
                "counts": counts,
                "perfect_leader": str(leader) == PERFECT_LEADER,
                "leader": leader,
                "perfect_guide": len(guide_matches) > 0,
                "guide": guide,
                "perfect_donor": len(donor_matches) > 0,
                "perfect_middle_donor_region": len(middle_matches) > 0,
                "donor": donor,
                "donor_bc0_fragment": donor_bc0_fragment,
                "unambiguous_guide": len(guide_matches) == 1,
                "unambiguous_donor": len(donor_matches) == 1,
                "oligo_name": result.oligo_name if matched else "NA",
                "match_type": result.match_type,
            })

        out_df = pd.DataFrame(rows)
        out_df.to_csv(
            output_path,
            sep="\t",
            index=False,
            header=first_chunk,
            mode="w" if first_chunk else "a",
        )
        first_chunk = False
        logger.info(
            "  %s: %d sequences processed (%d matched)",
            input_path.stem, stats["total"], stats["matched"],
        )

    return stats


def calculate_bc0_purity(
    bc0_counts: pd.DataFrame,
    oligo_col: str = "oligo_name",
    bc0_col: str = "bc0",
    count_col: str = "count",
) -> pd.DataFrame:
    """
    Calculate BC0 purity (fraction of reads from dominant oligo).

    Args:
        bc0_counts: DataFrame with bc0, oligo_name, count columns
        oligo_col: Name of oligo column
        bc0_col: Name of BC0 column
        count_col: Name of count column

    Returns:
        DataFrame with bc0, dominant_oligo, purity, total_reads columns
    """
    # Group by BC0 and calculate stats
    bc0_stats = bc0_counts.groupby(bc0_col).agg({
        count_col: ['sum', 'max'],
        oligo_col: lambda x: x.loc[x.idxmax()] if len(x) > 0 else None,
    }).reset_index()

    bc0_stats.columns = [bc0_col, 'total_reads', 'max_reads', 'dominant_oligo']
    bc0_stats['purity'] = bc0_stats['max_reads'] / bc0_stats['total_reads']

    return bc0_stats[[bc0_col, 'dominant_oligo', 'purity', 'total_reads']]
