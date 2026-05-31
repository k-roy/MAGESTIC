#!/usr/bin/env python3
"""
Step 15: Aggregate Hit Calling Results Across Comparisons

Combines per-comparison hit calling results from Step 14 and runs:
1. Cross-tier validation: Identify variants significant in both tier_1 AND tier_2
2. Nuclease-library aggregation: Concordance across SpG_NGG, SpCas9_NGG, etc.
3. Amino acid aggregation: Group variants producing the same AA change
4. Category stratification: Breakdown by category and codon_variant_type
5. Time series validation: Fitness slope analysis for multi-timepoint screens

This step runs AFTER individual comparison hit calling (Step 14) and provides
additional validation and summarization across the entire screen.

Usage:
    python 15_aggregate_hits.py \\
        --input-dir hit_calling/per_comparison/ \\
        --output-dir hit_calling/aggregated/ \\
        --screen-name 20250912_SpG_QTL_screen

    # For multi-timepoint screens:
    python 15_aggregate_hits.py \\
        --input-dir hit_calling/per_comparison/ \\
        --output-dir hit_calling/aggregated/ \\
        --screen-name 20250912_SpG_QTL_screen \\
        --screen-dates 20251014 20251018 \\
        --multi-timepoint

Author: Kevin R. Roy
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from datetime import timedelta
import time

# Add common scripts to path

import pandas as pd
from tqdm import tqdm

from magestic.pipelines.bc1_from_screen.core.aggregation import (
    stratify_hits_by_category,
    aggregate_across_nuclease_libraries,
    aggregate_by_amino_acid_change,
    aggregate_by_tier,
    run_all_aggregations,
)
from magestic.pipelines.bc1_from_screen.core.time_series import (
    TimeSeriesConfig,
    analyze_fitness_slopes_multi_timepoint,
    add_time_series_validation_to_hits,
)

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Aggregate hit calling results across comparisons"
    )
    parser.add_argument(
        "--input-dir", type=Path, required=True,
        help="Directory containing per-comparison hit files (hits_*.tsv)"
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True,
        help="Output directory for aggregated results"
    )
    parser.add_argument(
        "--screen-name", type=str, default="screen",
        help="Screen name for output file prefixes"
    )

    # Multi-timepoint options
    parser.add_argument(
        "--multi-timepoint", action="store_true",
        help="Enable time series analysis for multi-timepoint screens"
    )
    parser.add_argument(
        "--screen-dates", type=str, nargs="+", default=[],
        help="Screen date prefixes for multi-timepoint screens (e.g., 20251014 20251018)"
    )
    parser.add_argument(
        "--min-timepoints", type=int, default=3,
        help="Minimum timepoints for time series validation (default: 3)"
    )
    parser.add_argument(
        "--min-r2", type=float, default=0.7,
        help="Minimum R² for time series validation (default: 0.7)"
    )
    parser.add_argument(
        "--min-monotonic-score", type=float, default=0.75,
        help="Minimum monotonic score for time series validation (default: 0.75)"
    )

    # Column name overrides
    parser.add_argument(
        "--variant-col", type=str, default="sublibrary_variant_combo",
        help="Column name for variant identifiers"
    )
    parser.add_argument(
        "--comparison-col", type=str, default="comparison",
        help="Column name for comparison identifiers"
    )

    return parser.parse_args()


def load_all_hit_files(input_dir: Path) -> pd.DataFrame:
    """Load all per-comparison hit files from directory."""
    hit_files = sorted(input_dir.glob("hits_*.tsv"))

    if len(hit_files) == 0:
        raise FileNotFoundError(
            f"No hit files found in {input_dir}. "
            "Run step_g0 for all comparisons first."
        )

    logger.info(f"Found {len(hit_files)} hit files in {input_dir}")

    all_hits = []
    for f in tqdm(hit_files, desc="Loading hits"):
        df = pd.read_csv(f, sep="\t", low_memory=False)
        all_hits.append(df)

    hits_df = pd.concat(all_hits, ignore_index=True)
    logger.info(f"Combined: {len(hits_df):,} rows from {len(hit_files)} comparisons")

    # Report tier breakdown if available
    if "tier" in hits_df.columns:
        tier_counts = hits_df["tier"].value_counts()
        logger.info("Tier breakdown:")
        for tier, count in sorted(tier_counts.items()):
            n_sig = (hits_df[hits_df["tier"] == tier]["is_significant"]).sum()
            logger.info(f"  {tier}: {count:,} rows, {n_sig:,} significant")

    return hits_df


def main():
    args = parse_args()
    start_time = time.time()

    logger.info("=" * 80)
    logger.info("STEP G1: AGGREGATE HIT CALLING RESULTS")
    logger.info("=" * 80)

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load all per-comparison hit files
    hits_df = load_all_hit_files(args.input_dir)

    # Summary statistics
    n_comparisons = hits_df[args.comparison_col].nunique()
    n_variants = hits_df[args.variant_col].nunique()
    n_significant = hits_df["is_significant"].sum()

    logger.info(f"\nSummary:")
    logger.info(f"  Comparisons: {n_comparisons}")
    logger.info(f"  Unique variants: {n_variants:,}")
    logger.info(f"  Significant hits: {n_significant:,}")

    # Run aggregations
    logger.info("\n" + "=" * 80)
    logger.info("RUNNING AGGREGATIONS")
    logger.info("=" * 80)

    results = run_all_aggregations(hits_df, output_dir=str(args.output_dir))

    # Time series analysis (if enabled)
    time_series_df = pd.DataFrame()
    if args.multi_timepoint and args.screen_dates:
        logger.info("\n5. Time series / fitness slope analysis...")

        ts_config = TimeSeriesConfig(
            min_timepoints=args.min_timepoints,
            min_r2_for_validation=args.min_r2,
            min_monotonic_score=args.min_monotonic_score,
            require_direction_consistent=True,
            screen_dates=args.screen_dates,
        )

        time_series_df = analyze_fitness_slopes_multi_timepoint(
            hits_df,
            ts_config,
            variant_col=args.variant_col,
            comparison_col=args.comparison_col,
        )

        if len(time_series_df) > 0:
            # Save time series results
            time_series_df.to_csv(
                args.output_dir / "fitness_slope_analysis_multi_timepoint.tsv",
                sep="\t", index=False
            )

            # Save time series validated subset
            ts_validated = time_series_df[time_series_df["time_series_validated"]]
            ts_validated.to_csv(
                args.output_dir / "time_series_validated_hits.tsv",
                sep="\t", index=False
            )

            logger.info(f"  Time series validated: {len(ts_validated):,}")

            # Add time series validation to hits
            hits_df = add_time_series_validation_to_hits(
                hits_df, time_series_df, variant_col=args.variant_col
            )

    # Save main combined results
    logger.info("\nSaving combined results...")

    # All hits
    hits_df.to_csv(
        args.output_dir / f"{args.screen_name}_variant_hits_all.tsv",
        sep="\t", index=False
    )

    # Significant hits only
    sig_hits = hits_df[hits_df["is_significant"]].copy()
    if "mean_log2fc" in sig_hits.columns:
        sig_hits = sig_hits.sort_values("mean_log2fc", key=abs, ascending=False)
    sig_hits.to_csv(
        args.output_dir / f"{args.screen_name}_variant_hits_significant.tsv",
        sep="\t", index=False
    )

    # Generate summary statistics
    stats = {
        "screen_name": args.screen_name,
        "n_comparisons": n_comparisons,
        "n_variants": n_variants,
        "n_significant_total": int(n_significant),
    }

    # Add tier-specific stats
    if "tier" in hits_df.columns:
        for tier in hits_df["tier"].unique():
            tier_df = hits_df[hits_df["tier"] == tier]
            stats[f"n_{tier}_rows"] = len(tier_df)
            stats[f"n_{tier}_significant"] = int(tier_df["is_significant"].sum())

    # Add cross-tier validation stats
    if len(results.get("tier", pd.DataFrame())) > 0:
        tier_agg = results["tier"]
        stats["n_cross_tier_validated"] = int(tier_agg["cross_tier_validated"].sum())
        stats["n_tier1_only"] = int(
            (tier_agg["tier1_significant"] & ~tier_agg["tier2_significant"]).sum()
        )
        stats["n_tier2_only"] = int(
            (tier_agg["tier2_significant"] & ~tier_agg["tier1_significant"]).sum()
        )

    # Add nuclease-library aggregation stats
    if len(results.get("nuclease_library", pd.DataFrame())) > 0:
        nuc_agg = results["nuclease_library"]
        stats["n_nuclease_agg_significant"] = int(
            nuc_agg["is_significant_aggregated"].sum()
        )
        stats["n_ngnh_exceptions"] = int(
            nuc_agg["ngnh_low_efficiency_exception"].sum()
        )

    # Add time series stats
    if len(time_series_df) > 0:
        stats["n_time_series_analyzed"] = len(time_series_df)
        stats["n_time_series_validated"] = int(
            time_series_df["time_series_validated"].sum()
        )
        stats["mean_fitness_slope_r2"] = float(
            time_series_df["fitness_slope_r2"].mean()
        )

    # Save stats
    with open(args.output_dir / "aggregation_stats.json", "w") as f:
        json.dump(stats, f, indent=2, default=str)

    # Final summary
    elapsed = time.time() - start_time
    logger.info("\n" + "=" * 80)
    logger.info(f"COMPLETE in {timedelta(seconds=int(elapsed))}")
    logger.info("=" * 80)
    logger.info(f"\nOutput directory: {args.output_dir}")
    logger.info("Key files:")
    logger.info(f"  {args.screen_name}_variant_hits_significant.tsv - Final significant hits")
    logger.info(f"  hits_summary_by_tier.tsv - Tier breakdown")
    logger.info(f"  hits_cross_tier_validated.tsv - Cross-tier validated hits")
    logger.info(f"  hits_aggregated_by_nuclease_library.tsv - Library concordance")
    logger.info(f"  hits_aggregated_by_amino_acid.tsv - AA-level aggregation")
    if args.multi_timepoint:
        logger.info(f"  fitness_slope_analysis_multi_timepoint.tsv - Time series analysis")
        logger.info(f"  time_series_validated_hits.tsv - Time-validated hits")


if __name__ == "__main__":
    main()
