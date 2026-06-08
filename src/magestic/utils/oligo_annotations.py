"""
Oligo Annotations Utility Module

Provides unified access to the pre-parsed harmonized QTL oligo table,
eliminating the need for runtime oligo name parsing.

The harmonized table contains 70,656 oligos from:
- 20240411_SpG_QTL (V520, V521): SpG NGNG/NGNH
- 20210422_SpCas9_LbCas12a_QTL (V416, V419): SpCas9 NGG, LbCas12a

TERMINOLOGY (see common/scripts/TERMINOLOGY.md for full reference):
-----------
Gene name columns have TWO distinct concepts - do not conflate:

1. LOCUS columns (QTL being dissected):
   - systematic_gene_name_locus, common_gene_name_locus
   - Use for: hit aggregation, QTL mapping attribution

2. OVERLAPPING columns (CDS where variant resides):
   - overlapping_systematic_gene_name, overlapping_common_gene_name
   - Use for: variant consequence, CDS-level analysis

~14% of variants have locus != overlapping (e.g., RPS19A locus -> SMF1 CDS).

Usage:
    from magestic.utils.oligo_annotations import (
        load_oligo_annotations,
        get_annotation_columns_for_bc1_from_screen,
        get_annotation_columns_for_bc1_donor_bc0,
    )

    # Load full table (cached after first load)
    oligo_df = load_oligo_annotations()

    # Get specific columns for bc1_from_screen pipeline
    annot_df = load_oligo_annotations(
        columns=get_annotation_columns_for_bc1_from_screen()
    )

    # Merge with DESeq2 results
    results_annotated = results_df.merge(annot_df, on="oligo_name", how="left")

Author: Kevin R. Roy
Date: 2026-02-20
Updated: 2026-02-21 (added terminology reference)
"""

import logging
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from .path_utils import HARMONIZED_VARIANTS_FILE

logger = logging.getLogger(__name__)

# Canonical path to harmonized oligo table
HARMONIZED_OLIGO_TABLE = HARMONIZED_VARIANTS_FILE

# Column name mappings: harmonized table -> legacy annotation.py names
# Use these to maintain backward compatibility with existing code
COLUMN_MAPPING_BC1_FROM_SCREEN = {
    # Harmonized table column -> bc1_from_screen annotation.py expected name
    "oligo_name": "oligo_name",
    "guide_nuclease_PAM": "cas_type",
    "systematic_gene_name_locus": "gene_systematic",
    "common_gene_name_locus": "gene_common",
    "first_aa_num": "codon_position",
    "ref_aas": "original_aa",
    "alt_aas": "new_aa",
    "ref_codons": "original_codon",
    "alt_codons": "new_codon",
    "num_US_syn_changes": "upstream_changes",
    "num_DS_syn_changes": "downstream_changes",
    "donor_strand": "strand",
    "chrom": "chromosome",
    "PAM_strand": "guide_strand",
    "PAM_coord": "guide_position",
    "PAM": "pam_sequence",
    "leftmost_codon_coord": "edit_start",
    "rightmost_codon_coord": "edit_end",
    "codon_variant_type": "variant_type",
}

# Additional useful columns from harmonized table (not in legacy annotation.py)
EXTRA_COLUMNS_BC1_FROM_SCREEN = [
    "aa_change",                     # Human-readable AA change string
    "shared_aa_missense_change",     # For AA-level aggregation
    "SV_prediction_score",           # SCORE-SV editing prediction
    "individual_variants",           # Raw variant coordinates
    "category",                      # natural_variants vs codon_changes
    "is_sub_variant",                # Linkage relationship flag
    "sub_variant_of",                # Parent oligo names
    "is_incomplete_frameshift",      # Artifact flag
    "SPS",                           # Sublibrary primer sequence
    "guide_donor_bc0_plasmid_library",  # V library (V520, V521, V416, V419)
    "oligo_pool",                    # Twist pool identifier
]

# Columns needed for bc1_donor_bc0 pipeline
COLUMNS_BC1_DONOR_BC0 = [
    "oligo_name",
    "oligo_sequence",  # Full 200bp sequence
    "guide",           # 20bp guide
    "donor",           # ~129bp donor
    "guide_nuclease_PAM",
    "SPS",
    "guide_donor_bc0_plasmid_library",
    "codon_variant_type",
    "systematic_gene_name_locus",
    "common_gene_name_locus",
]


def get_annotation_columns_for_bc1_from_screen(
    include_extras: bool = True,
    use_legacy_names: bool = True,
) -> List[str]:
    """
    Get list of columns needed for bc1_from_screen pipeline.

    Args:
        include_extras: Include additional useful columns beyond legacy annotation.py
        use_legacy_names: Return legacy column names (for backward compatibility)

    Returns:
        List of column names to load from harmonized table
    """
    if use_legacy_names:
        # Return harmonized table column names (will be renamed after load)
        cols = list(COLUMN_MAPPING_BC1_FROM_SCREEN.keys())
    else:
        cols = list(COLUMN_MAPPING_BC1_FROM_SCREEN.keys())

    if include_extras:
        cols.extend(EXTRA_COLUMNS_BC1_FROM_SCREEN)

    return cols


def get_annotation_columns_for_bc1_donor_bc0() -> List[str]:
    """
    Get list of columns needed for bc1_donor_bc0 pipeline.

    Returns:
        List of column names to load from harmonized table
    """
    return COLUMNS_BC1_DONOR_BC0.copy()


@lru_cache(maxsize=4)
def _load_cached(
    filepath: str,
    columns_tuple: Optional[tuple] = None,
) -> pd.DataFrame:
    """
    Internal cached loader. Uses tuple for columns since lists aren't hashable.
    """
    columns = list(columns_tuple) if columns_tuple else None

    logger.info(f"Loading oligo annotations from {filepath}")
    df = pd.read_csv(filepath, sep="\t", usecols=columns)
    logger.info(f"  Loaded {len(df):,} oligos with {len(df.columns)} columns")

    return df


def load_oligo_annotations(
    columns: Optional[List[str]] = None,
    rename_to_legacy: bool = False,
    add_is_control: bool = True,
    filepath: Optional[Path] = None,
) -> pd.DataFrame:
    """
    Load harmonized oligo annotations table.

    This is the primary function for accessing oligo annotations.
    Results are cached for performance.

    Args:
        columns: Specific columns to load. If None, loads all columns.
        rename_to_legacy: Rename columns to match legacy annotation.py names
        add_is_control: Add 'is_control' column (True for background_variant oligos)
        filepath: Override default harmonized table path

    Returns:
        DataFrame with oligo annotations
    """
    if filepath is None:
        filepath = HARMONIZED_OLIGO_TABLE

    if not filepath.exists():
        raise FileNotFoundError(
            f"Harmonized oligo table not found: {filepath}\n"
            "Run the harmonization script first."
        )

    # Convert columns to tuple for caching
    columns_tuple = tuple(columns) if columns else None

    df = _load_cached(str(filepath), columns_tuple).copy()

    # Add is_control column
    if add_is_control:
        if "systematic_gene_name_locus" in df.columns:
            df["is_control"] = (
                df["systematic_gene_name_locus"] == "background_variant"
            )
        elif "oligo_name" in df.columns:
            # Fallback: check oligo_name for control patterns
            df["is_control"] = (
                df["oligo_name"].str.lower().str.contains("control|background", na=False)
            )
        else:
            df["is_control"] = False

    # Rename to legacy column names if requested
    if rename_to_legacy:
        rename_map = {
            old: new for old, new in COLUMN_MAPPING_BC1_FROM_SCREEN.items()
            if old in df.columns and old != new
        }
        df = df.rename(columns=rename_map)

    return df


def get_oligo_annotation(
    oligo_name: str,
    annotation_df: Optional[pd.DataFrame] = None,
) -> Optional[Dict]:
    """
    Get annotation for a single oligo by name.

    Args:
        oligo_name: The oligo name to look up
        annotation_df: Pre-loaded annotation DataFrame (for efficiency)

    Returns:
        Dictionary with annotation fields, or None if not found
    """
    if annotation_df is None:
        annotation_df = load_oligo_annotations()

    match = annotation_df[annotation_df["oligo_name"] == oligo_name]

    if len(match) == 0:
        return None

    return match.iloc[0].to_dict()


def annotate_dataframe(
    df: pd.DataFrame,
    oligo_col: str = "oligo_name",
    columns: Optional[List[str]] = None,
    rename_to_legacy: bool = False,
) -> pd.DataFrame:
    """
    Annotate a DataFrame by merging with oligo annotations.

    This is the recommended way to add annotations to DESeq2 results
    or other analysis outputs.

    Args:
        df: DataFrame with oligo names
        oligo_col: Column containing oligo names
        columns: Specific annotation columns to add. If None, adds all.
        rename_to_legacy: Use legacy column names from annotation.py

    Returns:
        DataFrame with annotation columns added
    """
    # Load annotations
    annotation_df = load_oligo_annotations(
        columns=columns,
        rename_to_legacy=rename_to_legacy,
        add_is_control=True,
    )

    # Merge
    result = df.merge(
        annotation_df,
        left_on=oligo_col,
        right_on="oligo_name",
        how="left",
    )

    # Report match rate
    n_matched = result["oligo_name"].notna().sum() if "oligo_name" in result.columns else 0
    n_total = len(df)
    match_pct = 100 * n_matched / n_total if n_total > 0 else 0

    logger.info(f"Annotated {n_matched:,}/{n_total:,} rows ({match_pct:.1f}%)")

    return result


# ============================================================================
# Convenience functions for specific use cases
# ============================================================================

def get_oligos_by_gene(
    gene: str,
    annotation_df: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Get all oligos targeting a specific gene."""
    if annotation_df is None:
        annotation_df = load_oligo_annotations()

    return annotation_df[
        (annotation_df["systematic_gene_name_locus"] == gene) |
        (annotation_df["common_gene_name_locus"] == gene)
    ]


def get_oligos_by_library(
    library: str,
    annotation_df: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Get all oligos from a specific V library (V416, V419, V520, V521)."""
    if annotation_df is None:
        annotation_df = load_oligo_annotations()

    return annotation_df[
        annotation_df["guide_donor_bc0_plasmid_library"] == library
    ]


def get_synonymous_oligos(
    annotation_df: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Get all synonymous variant oligos (used as null controls)."""
    if annotation_df is None:
        annotation_df = load_oligo_annotations()

    return annotation_df[
        annotation_df["codon_variant_type"] == "synonymous"
    ]


def get_missense_oligos(
    annotation_df: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Get all missense variant oligos."""
    if annotation_df is None:
        annotation_df = load_oligo_annotations()

    return annotation_df[
        annotation_df["codon_variant_type"] == "missense"
    ]


# ============================================================================
# Validation
# ============================================================================

def validate_harmonized_table(filepath: Optional[Path] = None) -> bool:
    """
    Validate that the harmonized table has expected structure.

    Returns:
        True if valid, raises ValueError otherwise
    """
    if filepath is None:
        filepath = HARMONIZED_OLIGO_TABLE

    df = pd.read_csv(filepath, sep="\t", nrows=100)

    # Check required columns
    required_cols = [
        "oligo_name", "oligo_sequence", "guide", "donor",
        "guide_nuclease_PAM", "codon_variant_type",
        "systematic_gene_name_locus", "common_gene_name_locus",
    ]

    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Harmonized table missing required columns: {missing}")

    # Check oligo_name uniqueness
    if df["oligo_name"].duplicated().any():
        raise ValueError("Harmonized table has duplicate oligo_name values")

    logger.info("Harmonized table validation passed")
    return True


if __name__ == "__main__":
    # Quick validation when run directly
    logging.basicConfig(level=logging.INFO)
    validate_harmonized_table()

    df = load_oligo_annotations()
    print(f"\nLoaded {len(df):,} oligos")
    print(f"Columns: {list(df.columns)}")
    print(f"\nVariant types:")
    print(df["codon_variant_type"].value_counts())
    print(f"\nLibraries:")
    print(df["guide_donor_bc0_plasmid_library"].value_counts())
