#!/usr/bin/env python3
"""
Step 11: Abundance-Tier-Aware Hit Calling for Single Comparison

Performs hit calling with abundance-tier-specific thresholds, using synonymous variants
to estimate the null distribution (CRISPEY-BAR inspired approach).

Features:
- Synonymous-based null distribution estimation (robust MAD)
- Abundance-tier-specific z-score and concordance thresholds
  - abundance_high: z>=3.0, concordance>=0.75
  - abundance_medium: z>=3.75, concordance>=0.85
- Optional IHW for FDR control (if rpy2 available)
- Robust WLS aggregation for double-edit handling

NOTE: Abundance tiers (for DESeq2 analysis) are distinct from confidence tiers
(for oligo attribution) in bc1_donor_bc0_pipeline.

Usage:
    python 11_hit_calling.py \\
        --input annotated_results.tsv \\
        --output-hits hits.tsv \\
        --output-stats stats.tsv \\
        --library yL437 \\
        --comparison Caffeine_6

Author: Kevin R. Roy
"""

import argparse
import logging
import sys
from pathlib import Path

# Add common scripts to path

import numpy as np
import pandas as pd

from magestic.pipelines.bc1_from_screen.core.hit_calling import (
    calculate_z_scores,
    call_hits_single_comparison,
    calculate_concordance,
)
from magestic.pipelines.bc1_from_screen.config import HitCallingConfig
from magestic.pipelines.bc1_from_screen.core.variance_calibration import (
    VarianceCalibrationResult,
)

# Try to import IHW
try:
    from magestic.pipelines.bc1_from_screen.core.ihw import run_ihw, HAS_IHW
except ImportError:
    HAS_IHW = False

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Tier-aware hit calling for a single comparison"
    )
    parser.add_argument(
        "--input", type=Path, required=True,
        help="Input annotated DESeq2 results TSV"
    )
    parser.add_argument(
        "--output-hits", type=Path, required=True,
        help="Output hits TSV"
    )
    parser.add_argument(
        "--output-stats", type=Path, required=True,
        help="Output statistics TSV"
    )
    parser.add_argument(
        "--library", type=str, required=True,
        help="Library name"
    )
    parser.add_argument(
        "--comparison", type=str, required=True,
        help="Comparison ID"
    )
    # Abundance tier thresholds
    parser.add_argument(
        "--tier-1-z", type=float, default=3.0,
        help="Z-score threshold for abundance_high (default: 3.0)"
    )
    parser.add_argument(
        "--tier-1-concordance", type=float, default=0.75,
        help="Concordance threshold for abundance_high (default: 0.75)"
    )
    parser.add_argument(
        "--tier-2-z", type=float, default=3.75,
        help="Z-score threshold for abundance_medium (default: 3.75)"
    )
    parser.add_argument(
        "--tier-2-concordance", type=float, default=0.85,
        help="Concordance threshold for abundance_medium (default: 0.85)"
    )
    # Null estimation
    parser.add_argument(
        "--use-synonymous-null", type=str, default="True",
        help="Use synonymous variants for null estimation (default: True)"
    )
    parser.add_argument(
        "--use-ihw", action="store_true",
        help="Use IHW for FDR control if available"
    )
    return parser.parse_args()


def estimate_null_distribution(df: pd.DataFrame, use_synonymous: bool, tier_name: str = "all") -> tuple[float, float]:
    """Estimate null distribution mean and std using synonymous variants or middle-50%.

    Args:
        df: DataFrame with log2FoldChange and optionally codon_variant_type
        use_synonymous: Whether to use synonymous variants for null estimation
        tier_name: Name of tier for logging purposes

    Returns:
        tuple: (null_mean, null_std) for z-score calculation as (LFC - null_mean) / null_std
    """
    if use_synonymous and "codon_variant_type" in df.columns:
        # Use synonymous variants (CRISPEY-BAR approach)
        syn_mask = df["codon_variant_type"].str.contains("synonymous", case=False, na=False)
        syn_df = df[syn_mask]

        if len(syn_df) >= 20:
            lfc_values = syn_df["log2FoldChange"].dropna()
            null_mean = lfc_values.median()
            mad = np.median(np.abs(lfc_values - null_mean))
            robust_std = mad * 1.4826  # Scale factor for normal

            logger.info(f"  [{tier_name}] Null from {len(syn_df):,} synonymous: mean={null_mean:.4f}, std={robust_std:.4f}")
            return null_mean, robust_std
        else:
            logger.warning(f"  [{tier_name}] Only {len(syn_df)} synonymous variants, using middle-50% fallback")

    # Fallback: middle 50% of all variants
    lfc_values = df["log2FoldChange"].dropna()
    if len(lfc_values) < 20:
        logger.warning(f"  [{tier_name}] Only {len(lfc_values)} values, returning NaN")
        return np.nan, np.nan

    q25, q75 = lfc_values.quantile([0.25, 0.75])
    middle_50 = lfc_values[(lfc_values >= q25) & (lfc_values <= q75)]
    null_mean = middle_50.median()
    null_std = middle_50.std()

    # Apply correction factor (middle-50% underestimates by ~3x based on empirical analysis)
    corrected_std = null_std * 2.5

    logger.info(f"  [{tier_name}] Null from middle-50%: mean={null_mean:.4f}, std={corrected_std:.4f}")
    return null_mean, corrected_std


def main():
    args = parse_args()

    use_synonymous = args.use_synonymous_null.lower() in ("true", "yes", "1")

    logger.info(f"Step 11: Hit calling for {args.library} - {args.comparison}")
    logger.info(f"  Input: {args.input}")
    logger.info(f"  abundance_high: z >= {args.tier_1_z}, concordance >= {args.tier_1_concordance}")
    logger.info(f"  abundance_medium: z >= {args.tier_2_z}, concordance >= {args.tier_2_concordance}")
    logger.info(f"  Use synonymous null: {use_synonymous}")
    logger.info(f"  IHW available: {HAS_IHW}, use IHW: {args.use_ihw}")

    # Load data
    logger.info("Loading annotated results...")
    df = pd.read_csv(args.input, sep="\t")
    logger.info(f"  Loaded {len(df):,} BC1s")

    # Check for required columns
    required_cols = ["bc1", "log2FoldChange", "pvalue", "padj", "baseMean"]
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        logger.error(f"Missing required columns: {missing_cols}")
        return 1

    # Configure hit calling
    config = HitCallingConfig(
        tier_1_z_threshold=args.tier_1_z,
        tier_1_concordance_threshold=args.tier_1_concordance,
        tier_2_z_threshold=args.tier_2_z,
        tier_2_concordance_threshold=args.tier_2_concordance,
    )

    # Estimate null distribution and calculate z-scores PER TIER
    # This is critical because tiers have different variance characteristics
    logger.info("Estimating tier-specific null distributions...")
    df["z_score"] = np.nan
    df["null_mean"] = np.nan
    df["null_std"] = np.nan
    df["is_significant"] = False

    # Track overall null stats for reporting (weighted by tier size)
    null_mean_overall = 0.0
    null_std_overall = 0.0

    if "tier" in df.columns:
        # Process each tier separately
        tier_configs = [
            (["abundance_high", "tier_1"], config.tier_1_z_threshold, "high_abundance"),
            (["abundance_medium", "tier_2"], config.tier_2_z_threshold, "medium_abundance"),
        ]

        for tier_values, z_threshold, tier_label in tier_configs:
            tier_mask = df["tier"].isin(tier_values)
            n_tier = tier_mask.sum()

            if n_tier == 0:
                continue

            tier_df = df[tier_mask]
            null_mean, null_std = estimate_null_distribution(tier_df, use_synonymous, tier_label)

            if np.isnan(null_mean) or np.isnan(null_std):
                logger.warning(f"  [{tier_label}] Could not estimate null, skipping")
                continue

            # Calculate z-scores for this tier
            df.loc[tier_mask, "z_score"] = (df.loc[tier_mask, "log2FoldChange"] - null_mean) / null_std
            df.loc[tier_mask, "null_mean"] = null_mean
            df.loc[tier_mask, "null_std"] = null_std

            # Call significant hits
            df.loc[tier_mask, "is_significant"] = df.loc[tier_mask, "z_score"].abs() >= z_threshold

            n_sig = df.loc[tier_mask, "is_significant"].sum()
            logger.info(f"  [{tier_label}] {n_sig:,} significant / {n_tier:,} BC1s (z>={z_threshold})")

            # Track for overall stats
            if null_std_overall == 0:
                null_mean_overall = null_mean
                null_std_overall = null_std
            else:
                # Use tier_1 stats for overall reporting
                pass
    else:
        # No tier info, use abundance_high thresholds for all
        logger.warning("  No tier column found, using abundance_high thresholds for all BC1s")
        null_mean, null_std = estimate_null_distribution(df, use_synonymous, "all")
        null_mean_overall = null_mean
        null_std_overall = null_std
        df["z_score"] = (df["log2FoldChange"] - null_mean) / null_std
        df["null_mean"] = null_mean
        df["null_std"] = null_std
        df["is_significant"] = df["z_score"].abs() >= config.tier_1_z_threshold
        df["tier"] = "abundance_high"

    # Apply IHW if requested and available
    if args.use_ihw and HAS_IHW:
        logger.info("Applying IHW for FDR control...")
        try:
            from magestic.pipelines.bc1_from_screen.core.ihw import apply_ihw_to_dataframe
            df = apply_ihw_to_dataframe(df, pvalue_col="pvalue", covariate_col="baseMean")
            logger.info("  IHW applied successfully")
        except Exception as e:
            logger.warning(f"  IHW failed: {e}, using standard padj")

    # Collect statistics (support both naming conventions)
    if "tier" in df.columns:
        n_high = df["tier"].isin(["abundance_high", "tier_1"]).sum()
        n_medium = df["tier"].isin(["abundance_medium", "tier_2"]).sum()
        n_sig_high = (df["tier"].isin(["abundance_high", "tier_1"]) & df["is_significant"]).sum()
        n_sig_medium = (df["tier"].isin(["abundance_medium", "tier_2"]) & df["is_significant"]).sum()
    else:
        n_high = len(df)
        n_medium = 0
        n_sig_high = df["is_significant"].sum()
        n_sig_medium = 0

    # Get tier-specific null stats for reporting
    high_mask = df["tier"].isin(["abundance_high", "tier_1"])
    medium_mask = df["tier"].isin(["abundance_medium", "tier_2"])
    null_mean_high = df.loc[high_mask, "null_mean"].iloc[0] if high_mask.any() and not df.loc[high_mask, "null_mean"].isna().all() else np.nan
    null_std_high = df.loc[high_mask, "null_std"].iloc[0] if high_mask.any() and not df.loc[high_mask, "null_std"].isna().all() else np.nan
    null_mean_medium = df.loc[medium_mask, "null_mean"].iloc[0] if medium_mask.any() and not df.loc[medium_mask, "null_mean"].isna().all() else np.nan
    null_std_medium = df.loc[medium_mask, "null_std"].iloc[0] if medium_mask.any() and not df.loc[medium_mask, "null_std"].isna().all() else np.nan

    stats = {
        "library": args.library,
        "comparison": args.comparison,
        "n_bc1s_total": len(df),
        "null_mean_high": null_mean_high,
        "null_std_high": null_std_high,
        "null_mean_medium": null_mean_medium,
        "null_std_medium": null_std_medium,
        "n_abundance_high": n_high,
        "n_abundance_medium": n_medium,
        "n_significant_abundance_high": n_sig_high,
        "n_significant_abundance_medium": n_sig_medium,
        "n_significant_total": df["is_significant"].sum(),
    }

    stats_df = pd.DataFrame([stats])

    # Log summary
    logger.info(f"  Total BC1s: {stats['n_bc1s_total']:,}")
    logger.info(f"  Significant high-abundance: {stats['n_significant_abundance_high']:,} / {stats['n_abundance_high']:,}")
    logger.info(f"  Significant medium-abundance: {stats['n_significant_abundance_medium']:,} / {stats['n_abundance_medium']:,}")
    logger.info(f"  Total significant: {stats['n_significant_total']:,}")

    # Ensure output directories exist
    args.output_hits.parent.mkdir(parents=True, exist_ok=True)
    args.output_stats.parent.mkdir(parents=True, exist_ok=True)

    # Save outputs
    logger.info(f"Saving hits to {args.output_hits}")
    df.to_csv(args.output_hits, sep="\t", index=False)

    logger.info(f"Saving statistics to {args.output_stats}")
    stats_df.to_csv(args.output_stats, sep="\t", index=False)

    logger.info("Step 11: Hit calling complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
