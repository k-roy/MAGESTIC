"""
Oligo matching functions for the Guide-Donor-BC0 pipeline.

Provides functions to match guide-donor sequences to designed oligos
using exact and partial matching strategies.
"""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np
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


def load_oligo_pool_harmonized(
    v_libraries: Optional[List[str]] = None,
    nuclease_types: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Load designed oligos from the harmonized QTL table instead of slicing a raw
    200mer order file.

    The harmonized table carries ``guide`` / ``donor`` already extracted *per
    nuclease*. This fixes LbCas12a matching: the fixed [20:40]/[51:180] slice in
    :func:`load_oligo_pool` is correct only for the SpCas9/SpG 200mer layout and
    mis-extracts the LbCas12a guide/donor (~0%% matches against a raw LbCas12a
    design). Optional ``v_libraries`` / ``nuclease_types`` scoping restricts the
    lookup to one library's nuclease/sublibrary, reducing cross-sublibrary
    donor/guide ambiguity on the dense 2024 SpG design.

    Returns the same schema as :func:`load_oligo_pool`
    (``oligo_name``, ``guide``, ``donor``, ``middle_donor``) so
    :func:`build_lookup_dictionaries` consumes it unchanged.
    """
    from magestic.utils.oligo_annotations import load_oligo_annotations

    cols = [
        "oligo_name", "guide", "donor",
        "guide_nuclease_PAM", "guide_donor_bc0_plasmid_library",
    ]
    df = load_oligo_annotations(columns=cols, add_is_control=False)

    if v_libraries:
        df = df[df["guide_donor_bc0_plasmid_library"].isin(v_libraries)]
    if nuclease_types:
        df = df[df["guide_nuclease_PAM"].isin(nuclease_types)]

    df = df.dropna(subset=["guide", "donor"]).copy()
    df["guide"] = df["guide"].astype(str).str.upper()
    df["donor"] = df["donor"].astype(str).str.upper()
    df["middle_donor"] = df["donor"].str.slice(MIDDLE_DONOR_START, MIDDLE_DONOR_END)
    return df[["oligo_name", "guide", "donor", "middle_donor"]]


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


def _derive_match_maps(lookups: Dict[str, Dict]) -> Dict[str, Dict]:
    """Precompute count + single-oligo maps from the list-valued lookups.

    Lets the per-row :func:`resolve_oligo_name` cascade be applied vectorized:
    ``*_count`` give the candidate count (for the perfect_/unambiguous_ flags and
    the cascade gates) and ``*_single`` give the lone oligo when a key is unique.
    Built once per call, O(#oligos).
    """
    g2o = lookups["guide_to_oligos"]
    d2o = lookups["perfect_donor_to_oligos"]
    m2o = lookups["partial_donor_to_oligos"]
    return {
        "guide_count": {k: len(v) for k, v in g2o.items()},
        "donor_count": {k: len(v) for k, v in d2o.items()},
        "middle_count": {k: len(v) for k, v in m2o.items()},
        "guide_single": {k: v[0] for k, v in g2o.items() if len(v) == 1},
        "donor_single": {k: v[0] for k, v in d2o.items() if len(v) == 1},
        "middle_single": {k: v[0] for k, v in m2o.items() if len(v) == 1},
    }


def process_counts_file(
    input_path: Path,
    output_path: Path,
    lookups: Dict[str, Dict],
    chunk_size: int = 100_000,
) -> Dict[str, int]:
    """
    Match every row in a Step-02 counts file to the designed oligo pool.

    Reads *input_path* (produced by :func:`fastq_processing.process_library`)
    in chunks and writes the matched result to *output_path*.

    PERFORMANCE (2026-06-12, R4): fully vectorized — replaces the prior
    row-by-row ``iterrows`` (which called :func:`resolve_oligo_name` per row) with
    a vectorized exact guide+donor ``Series.map`` plus an ``np.select`` cascade
    that reproduces ``resolve_oligo_name``'s branch order exactly
    (exact → both-unique-same → guide-only → donor-only → partial → unmatched).
    The flag columns (perfect_/unambiguous_) are computed by vectorized count
    lookups. Output is byte-identical to the legacy per-row version (verified by
    all-branch parity test); no per-row Python loop, so no joblib needed.

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
    maps = _derive_match_maps(lookups)
    gd2o = lookups["guide_donor_to_oligo"]
    first_chunk = True

    for chunk in pd.read_csv(input_path, sep="\t", chunksize=chunk_size):
        n = len(chunk)
        stats["total"] += n
        idx = chunk.index

        def col(name, default):
            return chunk[name] if name in chunk.columns else pd.Series([default] * n, index=idx)

        guide = col("guide", "NA")
        donor_raw = col("donor", "NA")
        leader = col("leader", "NA")
        guide_s = guide.astype(str)
        # Uppercased donor for lookups; NaN -> "" (matches the legacy
        # `str(donor) if not pd.isna(donor) else ""`). Output keeps original case.
        donor_up = donor_raw.where(donor_raw.notna(), "").astype(str).str.upper()
        nag = (guide.isna() | (guide_s == "NA")).to_numpy()

        # Vectorized candidate counts -> flag columns + cascade gates.
        ng = guide_s.map(maps["guide_count"]).fillna(0).to_numpy()
        nd = donor_up.map(maps["donor_count"]).fillna(0).to_numpy()
        middle = donor_up.str.slice(MIDDLE_DONOR_START, MIDDLE_DONOR_END)
        nm = middle.map(maps["middle_count"]).fillna(0).to_numpy()
        gsv = guide_s.map(maps["guide_single"]).to_numpy(dtype=object)
        dsv = donor_up.map(maps["donor_single"]).to_numpy(dtype=object)
        msv = middle.map(maps["middle_single"]).to_numpy(dtype=object)

        # Priority 1: exact guide+cloning_site+donor.
        oligo_exact = (guide_s + SPCAS9_CLONING_SITE + donor_up).map(gd2o)
        has_exact = oligo_exact.notna().to_numpy() & ~nag

        # Remaining cascade (resolve_oligo_name order), all vectorized.
        both_same = (~has_exact) & ~nag & (ng == 1) & (nd == 1) & (gsv == dsv)
        guide_only = (~has_exact) & (~both_same) & ~nag & (ng == 1)
        donor_only = (~has_exact) & (~both_same) & (~guide_only) & ~nag & (nd == 1)
        partial = (~has_exact) & (~both_same) & (~guide_only) & (~donor_only) & ~nag & (nm == 1)

        oligo = np.full(n, "NA", dtype=object)
        oligo[has_exact] = oligo_exact.to_numpy(dtype=object)[has_exact]
        oligo[both_same] = gsv[both_same]
        oligo[guide_only] = gsv[guide_only]
        oligo[donor_only] = dsv[donor_only]
        oligo[partial] = msv[partial]

        mtype = np.full(n, "unmatched", dtype=object)
        mtype[has_exact] = "exact"
        mtype[both_same] = "exact"
        mtype[guide_only] = "guide_only"
        mtype[donor_only] = "donor_only"
        mtype[partial] = "partial"

        stats["matched"] += int((oligo != "NA").sum())
        stats["perfect_match"] += int((mtype == "exact").sum())

        # Flag columns (legacy hardcodes all False for the NA-guide rows).
        perfect_leader = (leader.astype(str).to_numpy() == PERFECT_LEADER)
        perfect_guide = ng > 0
        perfect_donor = nd > 0
        perfect_middle = nm > 0
        unamb_guide = ng == 1
        unamb_donor = nd == 1
        for arr in (perfect_leader, perfect_guide, perfect_donor, perfect_middle,
                    unamb_guide, unamb_donor):
            arr[nag] = False

        out_df = pd.DataFrame({
            "library_ID": col("library_ID", "NA").to_numpy(),
            "scaffold": np.where(nag, "NA", "WT_SpCas9"),
            "bc0": col("bc0", "NA").to_numpy(),
            "counts": col("total_counts", 0).to_numpy(),
            "perfect_leader": perfect_leader,
            "leader": leader.to_numpy(),
            "perfect_guide": perfect_guide,
            "guide": guide.to_numpy(),
            "perfect_donor": perfect_donor,
            "perfect_middle_donor_region": perfect_middle,
            "donor": donor_raw.to_numpy(),
            "donor_bc0_fragment": col("donor_bc0_fragment", "NA").to_numpy(),
            "unambiguous_guide": unamb_guide,
            "unambiguous_donor": unamb_donor,
            "oligo_name": oligo,
            "match_type": mtype,
        })
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
