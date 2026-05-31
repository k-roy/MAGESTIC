#!/usr/bin/env python3
"""
Hit Aggregation Module
======================

This module provides functions for aggregating hit calling results across
multiple dimensions: tiers, nuclease libraries, amino acid changes, and
variant categories.

These aggregations happen AFTER the main hit calling pipeline (WLS, IHW,
tier-aware thresholds) and provide additional validation and summarization.

Key Features:
    - Cross-tier validation: Variants significant in both abundance_high AND abundance_medium
    - Nuclease-library aggregation: Concordance across SpG_NGG, SpCas9_NGG, etc.
    - Amino acid aggregation: Group variants producing same AA change
    - Category stratification: Breakdown by category and codon_variant_type

NOTE: Abundance tiers here refer to DESeq2 analysis strategy based on t0 read counts,
distinct from confidence tiers in bc1_donor_bc0_pipeline for oligo attribution.

Pipeline Position:
    This module works on the OUTPUT of hit calling, not during it:

    1. DESeq2 (per-comparison) → BC1-level results
    2. Variance calibration → weights from synonymous
    3. WLS aggregation → variant-level (robust M-estimators)
    4. IHW FDR control → covariate-adjusted p-values
    5. Tier-aware hit calling → z-score thresholds by tier
    6. [THIS MODULE] Aggregation refinements:
       a. aggregate_by_tier() → cross-tier validation
       b. aggregate_across_nuclease_libraries() → library concordance
       c. aggregate_by_amino_acid_change() → AA-level effects
       d. stratify_hits_by_category() → category breakdown

Usage:
    from bc1_from_screen_pipeline.core.aggregation import (
        aggregate_by_tier,
        aggregate_across_nuclease_libraries,
        aggregate_by_amino_acid_change,
        stratify_hits_by_category,
    )

    # Load combined hit calling results
    hits_df = pd.concat([pd.read_csv(f, sep='\t') for f in hit_files])

    # Run aggregations
    tier_agg = aggregate_by_tier(hits_df)
    nuclease_agg = aggregate_across_nuclease_libraries(hits_df)
    aa_agg = aggregate_by_amino_acid_change(hits_df)
    category_strata = stratify_hits_by_category(hits_df)

Author: Kevin R. Roy
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
from tqdm import tqdm
import logging

from .sublibrary import NUCLEASE_PAM_PREFIXES

logger = logging.getLogger(__name__)


# ============================================================================
# Category Stratification
# ============================================================================

def stratify_hits_by_category(
    hits_df: pd.DataFrame,
    category_col: str = "category",
    codon_type_col: str = "codon_variant_type",
    sig_col: str = "is_significant",
) -> Dict[str, pd.DataFrame]:
    """Stratify hits by category and codon_variant_type.

    Creates a dictionary of DataFrames, one per category__codon_variant_type
    combination. Useful for:
    - Category-specific hit rates
    - Comparing natural vs designed variants
    - Identifying category-specific effects

    Args:
        hits_df: DataFrame with hit calling results.
        category_col: Column name for variant category.
        codon_type_col: Column name for codon variant type.
        sig_col: Column name for significance flag.

    Returns:
        Dictionary mapping strata keys to filtered DataFrames.
        Keys are in format "category__codon_variant_type".
    """
    hits_df = hits_df.copy()

    # Create stratification key
    hits_df["strat_key"] = (
        hits_df[category_col].fillna("unknown").astype(str)
        + "__"
        + hits_df[codon_type_col].fillna("unknown").astype(str)
    )

    strata = {}
    for key in sorted(hits_df["strat_key"].unique()):
        strata[key] = hits_df[hits_df["strat_key"] == key].copy()

    # Log summary
    logger.info("Hits by category__codon_variant_type:")
    for key, df in strata.items():
        n_sig = df[sig_col].sum()
        n_total = len(df)
        pct = 100 * n_sig / n_total if n_total > 0 else 0
        logger.info(f"  {key}: {n_sig}/{n_total} significant ({pct:.1f}%)")

    return strata


# ============================================================================
# Nuclease-Library Aggregation
# ============================================================================

def aggregate_across_nuclease_libraries(
    hits_df: pd.DataFrame,
    variant_col: str = "sublibrary_variant_combo",
    comparison_col: str = "comparison",
    lfc_col: str = "mean_log2fc",
    concordance_col: str = "concordance_fraction",
    cv_col: str = "within_variant_cv",
    sig_col: str = "is_significant",
    zscore_col: str = "z_score",
    sv_score_col: str = "sv_score",
    gene_col: str = "common_gene_name_locus",
    category_col: str = "category",
    n_bc1_col: str = "n_bc1",
    nuclease_prefixes: Optional[List[str]] = None,
    min_concordance: float = 0.5,
) -> pd.DataFrame:
    """Aggregate evidence across nuclease-library combinations.

    For variants tested in multiple nuclease libraries (SpG_NGG, SpG_NGNG,
    SpCas9_NGG, etc.), aggregates results to determine overall significance.

    Includes NGNH efficiency exception: SpG_NGNH has known lower editing
    efficiency, so if a variant is significant in all non-NGNH libraries
    but not NGNH, it still passes aggregation.

    Args:
        hits_df: DataFrame with hit calling results.
        variant_col: Column name for variant IDs (includes nuclease prefix).
        comparison_col: Column name for comparison IDs.
        lfc_col: Column name for log2 fold change.
        concordance_col: Column name for concordance fraction.
        cv_col: Column name for coefficient of variation.
        sig_col: Column name for significance flag.
        zscore_col: Column name for z-scores.
        sv_score_col: Column name for sub-variant scores.
        gene_col: Column name for gene annotation.
        category_col: Column name for variant category.
        n_bc1_col: Column name for BC1 counts.
        nuclease_prefixes: List of nuclease prefixes to strip. If None, uses
            NUCLEASE_PAM_PREFIXES from sublibrary.py for consistency.
        min_concordance: Minimum library concordance for aggregated significance.

    Returns:
        DataFrame with one row per (comparison, variant_identity) with
        aggregated statistics across nuclease libraries.
    """
    if nuclease_prefixes is None:
        nuclease_prefixes = NUCLEASE_PAM_PREFIXES

    hits_df = hits_df.copy()

    # Extract variant identity (remove nuclease prefix)
    pattern = "^(" + "|".join(nuclease_prefixes) + ")_"
    hits_df["variant_identity"] = hits_df[variant_col].str.replace(
        pattern, "", regex=True
    )

    # Aggregate by comparison and variant_identity
    agg_dict = {
        variant_col: ["nunique", lambda x: sorted(x.unique())],
        n_bc1_col: "sum",
        lfc_col: ["mean", "median", "std"],
        concordance_col: "mean",
        cv_col: "mean",
        sig_col: ["any", "all", "sum"],
        zscore_col: ["mean", lambda x: max(x.abs()) if len(x) > 0 else np.nan],
    }

    # Add optional columns
    if sv_score_col in hits_df.columns:
        agg_dict[sv_score_col] = "first"
    if gene_col in hits_df.columns:
        agg_dict[gene_col] = lambda x: list(x.unique())
    if category_col in hits_df.columns:
        agg_dict[category_col] = "first"

    # Group and aggregate
    agg_df = hits_df.groupby([comparison_col, "variant_identity"]).agg(agg_dict)

    # Flatten column names
    agg_df.columns = [
        f"{a}_{b}" if b != "<lambda_0>" and b != "<lambda_1>" else
        ("nuclease_libraries" if a == variant_col else "max_abs_z_score")
        for a, b in agg_df.columns
    ]
    agg_df = agg_df.reset_index()

    # Rename columns for clarity
    rename_map = {
        f"{variant_col}_nunique": "n_nuclease_libraries",
        f"{n_bc1_col}_sum": "n_bc1_total",
        f"{lfc_col}_mean": "mean_log2fc",
        f"{lfc_col}_median": "median_log2fc",
        f"{lfc_col}_std": "std_log2fc_across_libraries",
        f"{concordance_col}_mean": "mean_concordance",
        f"{cv_col}_mean": "mean_cv",
        f"{sig_col}_any": "any_significant",
        f"{sig_col}_all": "all_significant",
        f"{sig_col}_sum": "n_significant",
        f"{zscore_col}_mean": "mean_z_score",
    }
    agg_df.rename(columns=rename_map, inplace=True)

    # Library concordance (protect against division by zero)
    agg_df["library_concordance"] = np.where(
        agg_df["n_nuclease_libraries"] > 0,
        agg_df["n_significant"] / agg_df["n_nuclease_libraries"],
        0.0
    )

    # NGNH efficiency exception - vectorized for performance
    # Build the logic using vectorized operations where possible
    # For rows with >1 library, check if NGNH is failing while others pass

    # First, identify which rows have NGNH in their libraries (need to iterate for list checking)
    has_ngnh = agg_df["nuclease_libraries"].apply(
        lambda libs: any("NGNH" in lib for lib in libs) if isinstance(libs, (list, set)) else False
    )

    # Count non-NGNH libraries for each row
    n_non_ngnh = agg_df["nuclease_libraries"].apply(
        lambda libs: len([lib for lib in libs if "NGNH" not in lib]) if isinstance(libs, (list, set)) else 0
    )

    # NGNH exception: more than 1 library, has NGNH, and significant count equals non-NGNH count
    agg_df["ngnh_low_efficiency_exception"] = (
        (agg_df["n_nuclease_libraries"] > 1) &
        has_ngnh &
        (agg_df["n_significant"] == n_non_ngnh)
    )

    # Aggregated significance
    agg_df["is_significant_aggregated"] = (
        agg_df["any_significant"] & (agg_df["library_concordance"] >= min_concordance)
    ) | agg_df["ngnh_low_efficiency_exception"]

    # Log summary
    logger.info("Nuclease-library aggregation:")
    logger.info(f"  Unique variants tested: {agg_df['variant_identity'].nunique()}")
    logger.info(f"  Significant (aggregated): {agg_df['is_significant_aggregated'].sum()}")
    logger.info(f"  NGNH efficiency exceptions: {agg_df['ngnh_low_efficiency_exception'].sum()}")

    return agg_df


# ============================================================================
# Amino Acid Aggregation
# ============================================================================

def aggregate_by_amino_acid_change(
    hits_df: pd.DataFrame,
    comparison_col: str = "comparison",
    variant_col: str = "sublibrary_variant_combo",
    aa_change_col: str = "aa_change",
    shared_aa_col: str = "shared_aa_missense_change",
    category_col: str = "category",
    lfc_col: str = "mean_log2fc",
    sig_col: str = "is_significant",
    zscore_col: str = "z_score",
    gene_col: str = "common_gene_name_locus",
    n_bc1_col: str = "n_bc1",
    use_shared_aa: bool = True,
    min_concordance: float = 0.5,
) -> pd.DataFrame:
    """Aggregate variants that produce the same amino acid change.

    Groups variants by their amino acid change (either shared_aa_missense_change
    or aa_change) and checks for concordance. This is useful for identifying:
    - Robust effects seen across multiple oligo designs
    - Interesting codon effects (when natural variant and designed codon differ)

    Args:
        hits_df: DataFrame with hit calling results.
        comparison_col: Column name for comparison IDs.
        variant_col: Column name for variant IDs.
        aa_change_col: Column name for amino acid change.
        shared_aa_col: Column name for shared amino acid change (preferred).
        category_col: Column name for variant category.
        lfc_col: Column name for log2 fold change.
        sig_col: Column name for significance flag.
        zscore_col: Column name for z-scores.
        gene_col: Column name for gene annotation.
        n_bc1_col: Column name for BC1 counts.
        use_shared_aa: If True and shared_aa_col exists, use it. Otherwise aa_change.
        min_concordance: Minimum design concordance for aggregated significance.

    Returns:
        DataFrame with one row per (comparison, amino_acid_change) with
        aggregated statistics.
    """
    hits_df = hits_df.copy()

    # Determine grouping column
    if use_shared_aa and shared_aa_col in hits_df.columns:
        aa_col = shared_aa_col
        logger.info(f"Using {shared_aa_col} for amino acid aggregation...")
    else:
        aa_col = aa_change_col
        logger.info(f"Using {aa_change_col} for amino acid aggregation...")

    # Check column exists
    if aa_col not in hits_df.columns:
        logger.warning(f"Column {aa_col} not found - skipping AA aggregation")
        return pd.DataFrame()

    # Filter for variants with AA changes
    aa_variants = hits_df[hits_df[aa_col].notna()].copy()

    if len(aa_variants) == 0:
        logger.info(f"  No variants with {aa_col} found")
        return pd.DataFrame()

    # Build aggregation
    agg_dict = {
        variant_col: ["nunique", lambda x: list(x.unique())],
        category_col: lambda x: sorted(x.unique()),
        n_bc1_col: "sum",
        lfc_col: ["mean", "std"],
        sig_col: ["any", "all", "sum"],
        zscore_col: "mean",
    }

    if gene_col in aa_variants.columns:
        agg_dict[gene_col] = lambda x: list(x.unique())
    if aa_change_col in aa_variants.columns and aa_change_col != aa_col:
        agg_dict[aa_change_col] = lambda x: list(x.dropna().unique())

    # Group and aggregate
    agg_df = aa_variants.groupby([comparison_col, aa_col]).agg(agg_dict)

    # Flatten column names
    agg_df.columns = [
        f"{a}_{b}" if b not in ["<lambda_0>", "<lambda_1>"] else
        ("variant_designs" if a == variant_col else
         "categories" if a == category_col else
         "genes" if a == gene_col else
         "aa_changes" if a == aa_change_col else a)
        for a, b in agg_df.columns
    ]
    agg_df = agg_df.reset_index()

    # Rename columns
    rename_map = {
        f"{variant_col}_nunique": "n_variant_designs",
        f"{n_bc1_col}_sum": "n_bc1_total",
        f"{lfc_col}_mean": "mean_log2fc",
        f"{lfc_col}_std": "std_log2fc",
        f"{sig_col}_any": "any_significant",
        f"{sig_col}_all": "all_significant",
        f"{sig_col}_sum": "n_significant",
        f"{zscore_col}_mean": "mean_z_score",
    }
    agg_df.rename(columns=rename_map, inplace=True)

    # Category concordance check - vectorized for performance
    # Check if both natural_variants and codon_changes are present in categories
    has_natural = agg_df["categories"].apply(
        lambda cats: "natural_variants" in cats if isinstance(cats, (list, set)) else False
    )
    has_codon = agg_df["categories"].apply(
        lambda cats: "codon_changes" in cats if isinstance(cats, (list, set)) else False
    )

    # Both categories present: they agree if all significant or none significant
    both_present = has_natural & has_codon
    agreement = agg_df["all_significant"] | ~agg_df["any_significant"]

    # Concordant if only one category (trivially true) or both agree
    agg_df["category_concordant"] = ~both_present | agreement
    agg_df["interesting_codon_effect"] = ~agg_df["category_concordant"]

    # Aggregated significance (protect against division by zero)
    design_concordance = np.where(
        agg_df["n_variant_designs"] > 0,
        agg_df["n_significant"] / agg_df["n_variant_designs"],
        0.0
    )
    agg_df["is_significant_aggregated"] = agg_df["any_significant"] & (
        design_concordance >= min_concordance
    )

    # Log summary
    logger.info(f"  Grouping column: {aa_col}")
    logger.info(f"  Unique AA changes tested: {len(agg_df)}")
    logger.info(f"  Significant (aggregated): {agg_df['is_significant_aggregated'].sum()}")
    logger.info(f"  Interesting codon effects flagged: {agg_df['interesting_codon_effect'].sum()}")

    return agg_df


# ============================================================================
# Tier Aggregation (Cross-Tier Validation)
# ============================================================================

def aggregate_by_tier(
    hits_df: pd.DataFrame,
    variant_col: str = "sublibrary_variant_combo",
    comparison_col: str = "comparison",
    tier_col: str = "tier",
    sig_col: str = "is_significant",
    lfc_col: str = "mean_log2fc",
    zscore_col: str = "z_score",
    category_col: str = "category",
    gene_col: str = "common_gene_name_locus",
) -> pd.DataFrame:
    """Aggregate significance counts by abundance tier and identify cross-tier validated variants.

    Cross-tier validation identifies variants that are significant in BOTH:
    - abundance_high: High-abundance individual BC1s (≥100 t0 reads, most reliable)
    - abundance_medium: Aggregated pseudo-BC1s from lower-abundance BC1s (20-99 t0 reads)

    Variants passing in both tiers are more robust, as they show the effect
    in two independent subsets of the data.

    Args:
        hits_df: DataFrame with hit calling results.
        variant_col: Column name for variant IDs.
        comparison_col: Column name for comparison IDs.
        tier_col: Column name for tier assignments.
        sig_col: Column name for significance flag.
        lfc_col: Column name for log2 fold change.
        zscore_col: Column name for z-scores.
        category_col: Column name for variant category.
        gene_col: Column name for gene annotation.

    Returns:
        DataFrame with one row per variant, containing tier-specific
        significance counts and cross-tier validation flag.
    """
    logger.info("Aggregating by tier...")

    if tier_col not in hits_df.columns:
        logger.warning(f"No {tier_col} column found - skipping tier aggregation")
        return pd.DataFrame()

    hits_df = hits_df.copy()

    # Group by variant across all comparisons and tiers
    def count_tier_significant(x, tier_name):
        """Count significant hits in a specific tier."""
        return (
            (hits_df.loc[x.index, tier_col] == tier_name) &
            hits_df.loc[x.index, sig_col]
        ).sum()

    agg_dict = {
        comparison_col: "nunique",
        tier_col: [
            lambda x: (x == "abundance_high").sum(),
            lambda x: (x == "abundance_medium").sum(),
        ],
        sig_col: [
            lambda x: count_tier_significant(x, "abundance_high"),
            lambda x: count_tier_significant(x, "abundance_medium"),
        ],
        lfc_col: "mean",
        zscore_col: "mean",
    }

    if category_col in hits_df.columns:
        agg_dict[category_col] = "first"
    if gene_col in hits_df.columns:
        agg_dict[gene_col] = "first"

    tier_agg = hits_df.groupby(variant_col).agg(agg_dict)

    # Flatten column names
    tier_agg.columns = [
        "n_comparisons" if a == comparison_col else
        "n_abundance_high_rows" if a == tier_col and "abundance_high" in str(b) else
        "n_abundance_medium_rows" if a == tier_col and "abundance_medium" in str(b) else
        "n_abundance_high_significant" if a == sig_col and "abundance_high" in str(b) else
        "n_abundance_medium_significant" if a == sig_col and "abundance_medium" in str(b) else
        f"{a}_{b}" if b and b != "first" else a
        for a, b in tier_agg.columns
    ]

    tier_agg = tier_agg.reset_index()

    # Rename remaining columns
    rename_map = {
        f"{lfc_col}_mean": "mean_log2fc",
        f"{zscore_col}_mean": "mean_z_score",
    }
    tier_agg.rename(columns=rename_map, inplace=True, errors="ignore")

    # Cross-tier validation flags
    tier_agg["abundance_high_significant"] = tier_agg["n_abundance_high_significant"] > 0
    tier_agg["abundance_medium_significant"] = tier_agg["n_abundance_medium_significant"] > 0
    tier_agg["cross_tier_validated"] = (
        tier_agg["abundance_high_significant"] & tier_agg["abundance_medium_significant"]
    )

    # Log summary
    n_cross_tier = tier_agg["cross_tier_validated"].sum()
    n_high_only = (tier_agg["abundance_high_significant"] & ~tier_agg["abundance_medium_significant"]).sum()
    n_medium_only = (tier_agg["abundance_medium_significant"] & ~tier_agg["abundance_high_significant"]).sum()

    logger.info(f"  abundance_high-only significant: {n_high_only:,}")
    logger.info(f"  abundance_medium-only significant: {n_medium_only:,}")
    logger.info(f"  Cross-tier validated: {n_cross_tier:,}")

    return tier_agg


# ============================================================================
# Convenience Function: Run All Aggregations
# ============================================================================

def run_all_aggregations(
    hits_df: pd.DataFrame,
    output_dir: Optional[str] = None,
) -> Dict[str, pd.DataFrame]:
    """Run all aggregation analyses on hit calling results.

    Convenience function that runs all four aggregation functions and
    optionally saves results to files.

    Args:
        hits_df: DataFrame with hit calling results.
        output_dir: Optional directory to save results. If provided, saves
            each aggregation to a separate TSV file.

    Returns:
        Dictionary with keys:
        - "category_strata": Dict of DataFrames by category
        - "nuclease_library": Nuclease-library aggregation
        - "amino_acid": Amino acid aggregation
        - "tier": Tier aggregation
    """
    from pathlib import Path

    results = {}

    # 1. Category stratification
    logger.info("\n1. Category stratification...")
    results["category_strata"] = stratify_hits_by_category(hits_df)

    # 2. Nuclease-library aggregation
    logger.info("\n2. Nuclease-library aggregation...")
    results["nuclease_library"] = aggregate_across_nuclease_libraries(hits_df)

    # 3. Amino acid aggregation
    logger.info("\n3. Amino acid level aggregation...")
    results["amino_acid"] = aggregate_by_amino_acid_change(hits_df)

    # 4. Tier aggregation
    logger.info("\n4. Tier aggregation...")
    results["tier"] = aggregate_by_tier(hits_df)

    # Save if output_dir provided
    if output_dir:
        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)

        # Save main aggregations
        if len(results["nuclease_library"]) > 0:
            results["nuclease_library"].to_csv(
                out_path / "hits_aggregated_by_nuclease_library.tsv",
                sep="\t", index=False
            )

        if len(results["amino_acid"]) > 0:
            results["amino_acid"].to_csv(
                out_path / "hits_aggregated_by_amino_acid.tsv",
                sep="\t", index=False
            )

        if len(results["tier"]) > 0:
            results["tier"].to_csv(
                out_path / "hits_summary_by_tier.tsv",
                sep="\t", index=False
            )
            # Also save cross-tier validated subset
            cross_tier = results["tier"][results["tier"]["cross_tier_validated"]]
            cross_tier.to_csv(
                out_path / "hits_cross_tier_validated.tsv",
                sep="\t", index=False
            )

        # Save category strata
        for cat, df in results["category_strata"].items():
            safe_name = cat.replace("/", "_")
            df.to_csv(
                out_path / f"hits_category_{safe_name}.tsv",
                sep="\t", index=False
            )

        logger.info(f"\nResults saved to: {out_path}")

    return results
