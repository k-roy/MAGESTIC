#!/usr/bin/env python3
"""
Time Series Validation Module
=============================

This module provides time-series analysis functions for multi-timepoint screens.
It complements the slope_fitness.py module by adding validation metrics that
assess whether fitness effects follow expected temporal patterns.

Key Features:
    - TimeSeriesConfig: Configuration for validation thresholds
    - monotonic_score: Fraction of timepoints following expected trend (0-1)
    - time_series_validated: Combined criterion (R² + monotonic + direction)
    - Condition-specific analysis for multi-condition screens

Relationship to slope_fitness.py:
    - slope_fitness.py: Calculates fitness slopes from raw counts (BC1-level)
    - time_series.py: Validates hit calls across timepoints (variant-level)

    These modules work together but at different stages:
    1. slope_fitness: counts → BC1 fitness slopes → variant aggregation
    2. time_series: hits across timepoints → validation metrics

Usage:
    from bc1_from_screen_pipeline.core.time_series import (
        TimeSeriesConfig,
        analyze_fitness_slopes_multi_timepoint,
    )

    ts_config = TimeSeriesConfig(min_timepoints=3, min_r2_for_validation=0.7)
    results = analyze_fitness_slopes_multi_timepoint(hits_df, ts_config)

Author: Kevin R. Roy

References:
    - FitSeq: Li et al. (2018) Nat Methods 15(4):267-270
    - CRISPEY: Mullis et al. (2023) Cell 186(6):1307
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from tqdm import tqdm
import logging

logger = logging.getLogger(__name__)


# ============================================================================
# Configuration
# ============================================================================

@dataclass
class TimeSeriesConfig:
    """Configuration for time series / fitness slope validation.

    Attributes:
        min_timepoints: Minimum timepoints required for slope calculation.
        min_r2_for_validation: Minimum R² for time_series_validated flag.
        min_monotonic_score: Minimum fraction of timepoints following expected trend.
        require_direction_consistent: All non-zero timepoints must have same sign.
        screen_dates: List of screen date prefixes for multi-timepoint screens
            (e.g., ["20251014", "20251018"]). Used to identify which comparisons
            belong to multi-timepoint screens.
        generations: List of generation values for multi-timepoint screens
            (e.g., [0, 6, 12, 21]).
    """
    min_timepoints: int = 3
    min_r2_for_validation: float = 0.7
    min_monotonic_score: float = 0.75
    require_direction_consistent: bool = True
    screen_dates: List[str] = field(default_factory=lambda: [])
    generations: List[int] = field(default_factory=lambda: [0, 6, 12, 21])


# ============================================================================
# Utility Functions
# ============================================================================

def extract_generation_from_comparison(comparison: str) -> Optional[int]:
    """Extract generation number from comparison ID.

    Assumes comparison IDs end with the generation number, e.g.:
        '20251014_20251018_Caffeine_6' -> 6
        '20251014_20251018_YPD_21' -> 21
        'Caffeine_20' -> 20

    Args:
        comparison: Comparison ID string.

    Returns:
        Generation number if found, None otherwise.
    """
    if pd.isna(comparison):
        return None

    parts = str(comparison).split("_")
    # Last part should be the generation number
    try:
        return int(parts[-1])
    except (ValueError, IndexError):
        return None


def extract_condition_from_comparison(comparison: str) -> Optional[str]:
    """Extract condition name from comparison ID.

    Handles multi-part condition names:
        '20251014_20251018_Caffeine_6' -> 'Caffeine'
        '20251014_20251018_Cadmium_chloride_12' -> 'Cadmium_chloride'
        'Caffeine_20' -> 'Caffeine'

    Args:
        comparison: Comparison ID string.

    Returns:
        Condition name if found, None otherwise.
    """
    if pd.isna(comparison):
        return None

    parts = str(comparison).split("_")
    # Find the screen date parts (8 digits) and condition parts
    condition_parts = []
    for i, part in enumerate(parts):
        # Skip leading date parts (8 digits)
        if part.isdigit() and len(part) == 8:
            continue
        # Skip trailing generation number
        if i == len(parts) - 1 and part.isdigit():
            continue
        condition_parts.append(part)

    return "_".join(condition_parts) if condition_parts else None


def is_multi_timepoint_screen(
    comparison: str,
    screen_dates: Optional[List[str]] = None,
) -> bool:
    """Check if comparison is from a multi-timepoint screen.

    Args:
        comparison: Comparison ID string.
        screen_dates: List of screen date prefixes that indicate multi-timepoint
            screens (e.g., ["20251014", "20251018", "20251014_20251018"]).

    Returns:
        True if comparison is from a multi-timepoint screen.
    """
    if pd.isna(comparison):
        return False

    if screen_dates is None or len(screen_dates) == 0:
        # No screen dates specified - can't determine
        return False

    comp_str = str(comparison)
    return any(date in comp_str for date in screen_dates)


# ============================================================================
# Time Series Analysis
# ============================================================================

def analyze_time_series_for_variant(
    variant_df: pd.DataFrame,
    ts_config: TimeSeriesConfig,
    lfc_col: str = "mean_log2fc",
    comparison_col: str = "comparison",
) -> Dict:
    """Analyze time series for a single variant across timepoints.

    For a variant that appears in multiple timepoints within the same condition,
    calculates fitness metrics including:
    - fitness_slope: log2FC change per generation (from linear regression)
    - fitness_slope_r2: R² of the linear fit
    - monotonic_score: Fraction of timepoints following expected trend
    - direction_consistent: Whether all non-zero timepoints have same sign
    - time_series_validated: Combined criterion passes all thresholds

    Args:
        variant_df: DataFrame with rows for this variant across timepoints.
            Must have columns for comparison, mean_log2fc, and generation.
        ts_config: Time series configuration.
        lfc_col: Column name for log2 fold change values.
        comparison_col: Column name for comparison IDs.

    Returns:
        Dictionary with time series metrics.
    """
    # Extract generation from comparison
    variant_df = variant_df.copy()
    variant_df["generation"] = variant_df[comparison_col].apply(
        extract_generation_from_comparison
    )
    variant_df = variant_df[variant_df["generation"].notna()].copy()

    if len(variant_df) < ts_config.min_timepoints:
        return {
            "fitness_slope": np.nan,
            "fitness_slope_r2": np.nan,
            "monotonic_score": np.nan,
            "direction_consistent": False,
            "time_series_validated": False,
            "n_timepoints": len(variant_df),
        }

    # Sort by generation
    variant_df = variant_df.sort_values("generation")
    generations = variant_df["generation"].values.astype(float)
    log2fc_values = variant_df[lfc_col].values

    # Linear regression for fitness slope
    if len(generations) >= 2:
        slope, intercept, r_value, p_value, std_err = scipy_stats.linregress(
            generations, log2fc_values
        )
        r2 = r_value ** 2
    else:
        slope, r2 = np.nan, np.nan

    # Monotonicity score: fraction of timepoints following expected trend
    # This is more informative than a simple boolean
    if len(log2fc_values) >= 2 and pd.notna(slope):
        # For negative effects (fitness cost), expect decreasing log2FC
        # For positive effects (fitness benefit), expect increasing log2FC
        # Note: slope > 0 with NaN slope returns False, which would incorrectly set direction = -1
        expected_direction = 1 if slope > 0 else -1

        n_following = 0
        for i in range(1, len(log2fc_values)):
            actual_direction = 1 if log2fc_values[i] >= log2fc_values[i-1] else -1
            # Count as following if matches expected OR no change
            if actual_direction == expected_direction or log2fc_values[i] == log2fc_values[i-1]:
                n_following += 1

        monotonic_score = n_following / (len(log2fc_values) - 1)
    else:
        monotonic_score = np.nan

    # Direction consistency: all non-zero timepoints same sign?
    non_zero = log2fc_values[log2fc_values != 0]
    if len(non_zero) >= 2:
        direction_consistent = (
            (non_zero > 0).all() or (non_zero < 0).all()
        )
    else:
        direction_consistent = True  # Only one non-zero point, trivially consistent

    # Time series validated if passes all criteria
    time_series_validated = (
        (not np.isnan(r2)) and
        (r2 >= ts_config.min_r2_for_validation) and
        (not np.isnan(monotonic_score)) and
        (monotonic_score >= ts_config.min_monotonic_score) and
        (not ts_config.require_direction_consistent or direction_consistent)
    )

    return {
        "fitness_slope": slope,
        "fitness_slope_r2": r2,
        "monotonic_score": monotonic_score,
        "direction_consistent": direction_consistent,
        "time_series_validated": time_series_validated,
        "n_timepoints": len(variant_df),
    }


def analyze_fitness_slopes_multi_timepoint(
    hits_df: pd.DataFrame,
    ts_config: Optional[TimeSeriesConfig] = None,
    variant_col: str = "sublibrary_variant_combo",
    comparison_col: str = "comparison",
    lfc_col: str = "mean_log2fc",
    sig_col: str = "is_significant",
    gene_col: str = "common_gene_name_locus",
    category_col: str = "category",
    show_progress: bool = True,
) -> pd.DataFrame:
    """Analyze fitness slopes for multi-timepoint screen variants.

    For each significant variant in multi-timepoint screens, calculates:
    - fitness_slope: log2FC change per generation
    - fitness_slope_r2: goodness of linear fit
    - monotonic_score: fraction following expected trend
    - time_series_validated: passes all QC criteria

    The analysis is performed per variant PER CONDITION. This means a variant
    tested in both Caffeine and Cobalt will have separate time-series results
    for each condition.

    Args:
        hits_df: DataFrame with hit calling results. Must have columns for
            variant ID, comparison, log2FC, significance, gene, and category.
        ts_config: Time series configuration. If None, uses default.
        variant_col: Column name for variant IDs.
        comparison_col: Column name for comparison IDs.
        lfc_col: Column name for log2 fold change values.
        sig_col: Column name for significance flag.
        gene_col: Column name for gene annotation.
        category_col: Column name for variant category.
        show_progress: Show progress bar.

    Returns:
        DataFrame with one row per variant-condition combination, containing
        time series validation metrics. Empty if no multi-timepoint data.
    """
    logger.info("Analyzing fitness slopes for multi-timepoint screens...")

    if ts_config is None:
        ts_config = TimeSeriesConfig()

    # Check for required columns
    required_cols = [variant_col, comparison_col, lfc_col, sig_col]
    missing = [c for c in required_cols if c not in hits_df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Filter to multi-timepoint screens
    hits_df = hits_df.copy()
    hits_df["is_multi_timepoint"] = hits_df[comparison_col].apply(
        lambda x: is_multi_timepoint_screen(x, ts_config.screen_dates)
    )
    multi_timepoint_df = hits_df[hits_df["is_multi_timepoint"]].copy()

    if len(multi_timepoint_df) == 0:
        logger.warning("No multi-timepoint screen data found - skipping time series analysis")
        return pd.DataFrame()

    logger.info(f"  Multi-timepoint screen rows: {len(multi_timepoint_df):,}")

    # Extract condition for grouping
    multi_timepoint_df["condition"] = multi_timepoint_df[comparison_col].apply(
        extract_condition_from_comparison
    )
    n_conditions = multi_timepoint_df["condition"].nunique()
    logger.info(f"  Conditions: {n_conditions}")

    # Get significant variants in multi-timepoint screens
    sig_multi_timepoint = multi_timepoint_df[multi_timepoint_df[sig_col]]
    sig_variants = sig_multi_timepoint[variant_col].unique()
    logger.info(f"  Significant variants in multi-timepoint screens: {len(sig_variants):,}")

    # Analyze time series for each variant within each condition
    results = []
    iterator = tqdm(sig_variants, desc="  Analyzing time series") if show_progress else sig_variants

    for variant in iterator:
        variant_df = multi_timepoint_df[multi_timepoint_df[variant_col] == variant]

        for condition in variant_df["condition"].unique():
            cond_df = variant_df[variant_df["condition"] == condition]

            if len(cond_df) < 2:
                continue

            ts_metrics = analyze_time_series_for_variant(
                cond_df,
                ts_config,
                lfc_col=lfc_col,
                comparison_col=comparison_col,
            )

            result = {
                variant_col: variant,
                "condition": condition,
                "mean_log2fc": cond_df[lfc_col].mean(),
                **ts_metrics,
            }

            # Add gene and category if available
            if gene_col in cond_df.columns:
                result["gene"] = cond_df[gene_col].iloc[0]
            if category_col in cond_df.columns:
                result["category"] = cond_df[category_col].iloc[0]

            results.append(result)

    if not results:
        logger.warning("No variants with sufficient timepoints for time series analysis")
        return pd.DataFrame()

    results_df = pd.DataFrame(results)

    # Summary stats
    n_validated = results_df["time_series_validated"].sum()
    logger.info(f"  Time series analysis results:")
    logger.info(f"    Variants analyzed: {len(results_df):,}")
    logger.info(f"    Time series validated: {n_validated:,} ({100*n_validated/len(results_df):.1f}%)")

    if len(results_df) > 0:
        mean_r2 = results_df["fitness_slope_r2"].mean()
        mean_mono = results_df["monotonic_score"].mean()
        logger.info(f"    Mean R²: {mean_r2:.3f}")
        logger.info(f"    Mean monotonic score: {mean_mono:.3f}")

    return results_df


# ============================================================================
# Integration Helpers
# ============================================================================

def add_time_series_validation_to_hits(
    hits_df: pd.DataFrame,
    time_series_df: pd.DataFrame,
    variant_col: str = "sublibrary_variant_combo",
) -> pd.DataFrame:
    """Add time series validation columns to hits DataFrame.

    Merges time series results back to the hits DataFrame, adding columns:
    - fitness_slope, fitness_slope_r2, monotonic_score
    - direction_consistent, time_series_validated

    Args:
        hits_df: DataFrame with hit calling results.
        time_series_df: DataFrame from analyze_fitness_slopes_multi_timepoint().
        variant_col: Column name for variant IDs.

    Returns:
        hits_df with additional time series validation columns.
    """
    if len(time_series_df) == 0:
        # Add empty columns
        for col in ["fitness_slope", "fitness_slope_r2", "monotonic_score",
                    "direction_consistent", "time_series_validated"]:
            hits_df[col] = np.nan
        return hits_df

    # Extract condition from comparison in hits_df
    hits_df = hits_df.copy()
    hits_df["_condition"] = hits_df["comparison"].apply(extract_condition_from_comparison)

    # Merge on variant + condition
    merge_cols = [variant_col, "condition", "fitness_slope", "fitness_slope_r2",
                  "monotonic_score", "direction_consistent", "time_series_validated"]
    ts_subset = time_series_df[[c for c in merge_cols if c in time_series_df.columns]].copy()

    hits_df = hits_df.merge(
        ts_subset,
        left_on=[variant_col, "_condition"],
        right_on=[variant_col, "condition"],
        how="left",
        suffixes=("", "_ts"),
    )

    # Clean up
    hits_df.drop(columns=["_condition", "condition"], errors="ignore", inplace=True)

    return hits_df
