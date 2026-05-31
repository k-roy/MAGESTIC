#!/usr/bin/env python3
"""
Step G0: Tier-Aware Hit Calling for Single Comparison

Performs hit calling with tier-specific thresholds, using synonymous variants
to estimate the null distribution (CRISPEY-BAR inspired approach).

Features:
- Synonymous-based null distribution estimation (robust MAD)
- Tier-specific z-score and concordance thresholds
- Optional IHW for FDR control (if rpy2 available)
- Robust WLS aggregation for double-edit handling

Usage:
    python step_g0_hit_calling.py \\
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
sys.path.insert(0, str(Path(__file__).parents[2]))

import numpy as np
import pandas as pd

from bc1_from_screen_pipeline.core.hit_calling import (
    calculate_z_scores,
    call_hits_by_tier,
    calculate_concordance,
    HitCallingConfig,
)
from bc1_from_screen_pipeline.core.variance_calibration import (
    estimate_null_from_synonymous,
    VarianceCalibrationResult,
)

# Try to import IHW
try:
    from bc1_from_screen_pipeline.core.ihw import run_ihw, HAS_IHW
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
    # Tier thresholds
    parser.add_argument(
        "--tier-1-z", type=float, default=3.0,
        help="Z-score threshold for tier_1 (default: 3.0)"
    )
    parser.add_argument(
        "--tier-1-concordance", type=float, default=0.75,
        help="Concordance threshold for tier_1 (default: 0.75)"
    )
    parser.add_argument(
        "--tier-2-z", type=float, default=3.75,
        help="Z-score threshold for tier_2 (default: 3.75)"
    )
    parser.add_argument(
        "--tier-2-concordance", type=float, default=0.85,
        help="Concordance threshold for tier_2 (default: 0.85)"
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


def estimate_null_distribution(df: pd.DataFrame, use_synonymous: bool) -> float:
    """Estimate null distribution std using synonymous variants or middle-50%."""
    if use_synonymous and "codon_variant_type" in df.columns:
        # Use synonymous variants (CRISPEY-BAR approach)
        syn_mask = df["codon_variant_type"].str.contains("synonymous", case=False, na=False)
        syn_df = df[syn_mask]

        if len(syn_df) >= 20:
            lfc_values = syn_df["log2FoldChange"].dropna()
            median_lfc = lfc_values.median()
            mad = np.median(np.abs(lfc_values - median_lfc))
            robust_std = mad * 1.4826  # Scale factor for normal

            logger.info(f"  Null estimation from {len(syn_df):,} synonymous variants")
            logger.info(f"    MAD-based robust std: {robust_std:.4f}")
            return robust_std
        else:
            logger.warning(f"  Only {len(syn_df)} synonymous variants, using middle-50% fallback")

    # Fallback: middle 50% of all variants
    lfc_values = df["log2FoldChange"].dropna()
    q25, q75 = lfc_values.quantile([0.25, 0.75])
    middle_50 = lfc_values[(lfc_values >= q25) & (lfc_values <= q75)]
    null_std = middle_50.std()

    # Apply correction factor (middle-50% underestimates by ~3x based on empirical analysis)
    corrected_std = null_std * 2.5

    logger.info(f"  Null estimation from middle-50% (corrected): {corrected_std:.4f}")
    return corrected_std


def main():
    args = parse_args()

    use_synonymous = args.use_synonymous_null.lower() in ("true", "yes", "1")

    logger.info(f"Step G0: Hit calling for {args.library} - {args.comparison}")
    logger.info(f"  Input: {args.input}")
    logger.info(f"  Tier 1: z >= {args.tier_1_z}, concordance >= {args.tier_1_concordance}")
    logger.info(f"  Tier 2: z >= {args.tier_2_z}, concordance >= {args.tier_2_concordance}")
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

    # Estimate null distribution
    logger.info("Estimating null distribution...")
    null_std = estimate_null_distribution(df, use_synonymous)

    # Calculate z-scores
    logger.info("Calculating z-scores...")
    df["z_score"] = df["log2FoldChange"] / null_std

    # Configure hit calling
    config = HitCallingConfig(
        tier_1_z_threshold=args.tier_1_z,
        tier_1_concordance_threshold=args.tier_1_concordance,
        tier_2_z_threshold=args.tier_2_z,
        tier_2_concordance_threshold=args.tier_2_concordance,
    )

    # Call hits by tier
    logger.info("Calling hits by tier...")
    if "tier" in df.columns:
        df = call_hits_by_tier(df, config)
    else:
        # No tier info, use tier_1 thresholds for all
        logger.warning("  No tier column found, using tier_1 thresholds for all BC1s")
        df["is_significant"] = df["z_score"].abs() >= config.tier_1_z_threshold
        df["tier"] = "tier_1"

    # Apply IHW if requested and available
    if args.use_ihw and HAS_IHW:
        logger.info("Applying IHW for FDR control...")
        try:
            from bc1_from_screen_pipeline.core.ihw import apply_ihw_to_dataframe
            df = apply_ihw_to_dataframe(df, pvalue_col="pvalue", covariate_col="baseMean")
            logger.info("  IHW applied successfully")
        except Exception as e:
            logger.warning(f"  IHW failed: {e}, using standard padj")

    # Collect statistics
    stats = {
        "library": args.library,
        "comparison": args.comparison,
        "n_bc1s_total": len(df),
        "null_std": null_std,
        "n_tier_1": (df["tier"] == "tier_1").sum() if "tier" in df.columns else len(df),
        "n_tier_2": (df["tier"] == "tier_2").sum() if "tier" in df.columns else 0,
        "n_significant_tier_1": ((df["tier"] == "tier_1") & df["is_significant"]).sum() if "tier" in df.columns else df["is_significant"].sum(),
        "n_significant_tier_2": ((df["tier"] == "tier_2") & df["is_significant"]).sum() if "tier" in df.columns else 0,
        "n_significant_total": df["is_significant"].sum(),
    }

    stats_df = pd.DataFrame([stats])

    # Log summary
    logger.info(f"  Total BC1s: {stats['n_bc1s_total']:,}")
    logger.info(f"  Significant tier_1: {stats['n_significant_tier_1']:,}")
    logger.info(f"  Significant tier_2: {stats['n_significant_tier_2']:,}")
    logger.info(f"  Total significant: {stats['n_significant_total']:,}")

    # Ensure output directories exist
    args.output_hits.parent.mkdir(parents=True, exist_ok=True)
    args.output_stats.parent.mkdir(parents=True, exist_ok=True)

    # Save outputs
    logger.info(f"Saving hits to {args.output_hits}")
    df.to_csv(args.output_hits, sep="\t", index=False)

    logger.info(f"Saving statistics to {args.output_stats}")
    stats_df.to_csv(args.output_stats, sep="\t", index=False)

    logger.info("Step G0: Hit calling complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
