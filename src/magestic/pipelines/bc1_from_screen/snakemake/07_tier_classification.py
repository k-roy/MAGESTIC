#!/usr/bin/env python3
"""
Step 04: Abundance Tier Classification for BC1 Barcodes

Classifies BC1 barcodes into abundance tiers based on t0 read counts within each library.

Abundance Tier Definitions:
- abundance_high (≥100 reads): Individual BC1 analysis in DESeq2
- abundance_medium (20-99 reads): Aggregate by oligo before DESeq2
- abundance_low (<20 reads): Excluded from analysis (low confidence)

NOTE: These "abundance tiers" are distinct from the "confidence_tier" column
in bc1_donor_bc0_pipeline which classifies oligo attribution confidence.

Usage:
    python 04_tier_classification.py \\
        --bc1-counts bc1_counts_long.tsv \\
        --bc1-reference bc1_reference_table.tsv \\
        --sample-key sample_key.tsv \\
        --library yL406 \\
        --output-tier-mapping tier_mapping.tsv \\
        --output-stats tier_stats.tsv

Author: Kevin R. Roy
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

# Add common scripts to path

from magestic.pipelines.bc1_from_screen.config import TierConfig

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Classify BC1s into abundance tiers based on t0 counts"
    )
    parser.add_argument(
        "--bc1-counts", type=Path, required=True,
        help="Input BC1 counts long-format TSV"
    )
    parser.add_argument(
        "--bc1-reference", type=Path, required=True,
        help="BC1 reference table with oligo assignments"
    )
    parser.add_argument(
        "--sample-key", type=Path, required=True,
        help="Sample key TSV with condition, generations, sample_name columns"
    )
    parser.add_argument(
        "--library", type=str, required=True,
        help="Library name (e.g., yL437, yL442, yL406)"
    )
    parser.add_argument(
        "--output-tier-mapping", type=Path, required=True,
        help="Output BC1 tier mapping TSV"
    )
    parser.add_argument(
        "--output-stats", type=Path, required=True,
        help="Output tier statistics TSV"
    )
    parser.add_argument(
        "--tier-1-min", type=int, default=100,
        help="Minimum t0 reads for abundance_high (default: 100)"
    )
    parser.add_argument(
        "--tier-2-min", type=int, default=20,
        help="Minimum t0 reads for abundance_medium (default: 20)"
    )
    parser.add_argument(
        "--t0-condition", type=str, default="YPD",
        help="Condition to identify t0 samples (default: YPD)"
    )
    parser.add_argument(
        "--screen-dates", type=str, nargs="+",
        help="Screen dates to filter for this library (optional)"
    )
    return parser.parse_args()


def calculate_t0_counts(
    bc1_counts: pd.DataFrame,
    sample_key: pd.DataFrame,
    t0_condition: str = "YPD",
    screen_dates: list = None,
    bc1_col: str = "bc1",
    sample_col: str = "sample_name",
    count_col: str = "counts",
) -> pd.DataFrame:
    """
    Calculate total t0 counts per BC1.

    Identifies t0 samples as those with:
    - condition matching t0_condition (default: YPD)
    - generations == 0
    - Optionally filtered by screen_dates

    Args:
        bc1_counts: Long-format BC1 counts (bc1, sample_name, counts)
        sample_key: Sample metadata with condition, generations, sample_name
        t0_condition: Condition value identifying t0 samples
        screen_dates: Optional list of screen dates to filter
        bc1_col: Column name for BC1 barcodes
        sample_col: Column name for sample identifiers
        count_col: Column name for read counts

    Returns:
        DataFrame with bc1 and t0_total columns
    """
    # Identify t0 samples from sample key
    # Handle generations as numeric, converting strings/floats to int for comparison
    generations_numeric = pd.to_numeric(sample_key["generations"], errors="coerce").fillna(-1).astype(int)
    t0_mask = (
        (sample_key["condition"] == t0_condition) &
        (generations_numeric == 0)
    )

    # Optionally filter by screen dates
    if screen_dates:
        t0_mask = t0_mask & (sample_key["screen_date"].isin(screen_dates))

    t0_samples = sample_key.loc[t0_mask, sample_col].unique().tolist()

    if not t0_samples:
        raise ValueError(
            f"No t0 samples found for condition='{t0_condition}', generations=0"
            + (f", screen_dates={screen_dates}" if screen_dates else "")
        )

    logger.info(f"  Found {len(t0_samples)} t0 samples")

    # Filter BC1 counts to t0 samples
    t0_counts = bc1_counts[bc1_counts[sample_col].isin(t0_samples)]

    # Sum counts per BC1 across all t0 samples
    bc1_t0_totals = (
        t0_counts
        .groupby(bc1_col)[count_col]
        .sum()
        .reset_index()
        .rename(columns={count_col: "t0_total"})
    )

    return bc1_t0_totals


def classify_bc1_tiers(
    bc1_t0_counts: pd.DataFrame,
    tier_config: TierConfig,
    bc1_col: str = "bc1",
    t0_col: str = "t0_total",
) -> pd.DataFrame:
    """
    Classify BC1s into tiers based on t0 read counts.

    Args:
        bc1_t0_counts: DataFrame with bc1 and t0_total columns
        tier_config: TierConfig with tier thresholds
        bc1_col: Column name for BC1 barcodes
        t0_col: Column name for t0 counts

    Returns:
        DataFrame with bc1, t0_total, and tier columns
    """
    result = bc1_t0_counts.copy()
    result["tier"] = result[t0_col].apply(tier_config.get_tier)
    return result


def create_tier_summary(tier_mapping: pd.DataFrame, library: str) -> pd.DataFrame:
    """
    Create summary statistics for abundance tier classification.

    Args:
        tier_mapping: DataFrame with tier assignments
        library: Library name

    Returns:
        DataFrame with tier summary statistics
    """
    summary_data = []

    for tier in ["abundance_high", "abundance_medium", "abundance_low"]:
        tier_data = tier_mapping[tier_mapping["tier"] == tier]
        summary_data.append({
            "library": library,
            "tier": tier,
            "n_bc1s": len(tier_data),
            "total_t0_reads": tier_data["t0_total"].sum() if len(tier_data) > 0 else 0,
            "mean_t0_reads": tier_data["t0_total"].mean() if len(tier_data) > 0 else 0,
            "median_t0_reads": tier_data["t0_total"].median() if len(tier_data) > 0 else 0,
        })

    return pd.DataFrame(summary_data)


def main():
    args = parse_args()

    logger.info(f"Step 04: Abundance tier classification for library {args.library}")
    logger.info(f"  Input BC1 counts: {args.bc1_counts}")
    logger.info(f"  BC1 reference: {args.bc1_reference}")
    logger.info(f"  Sample key: {args.sample_key}")
    logger.info(f"  t0 condition: {args.t0_condition}")
    logger.info(f"  Tier thresholds: abundance_high >= {args.tier_1_min}, abundance_medium >= {args.tier_2_min}")

    # Load BC1 counts
    logger.info("Loading BC1 counts...")
    bc1_counts = pd.read_csv(args.bc1_counts, sep="\t")
    logger.info(f"  Loaded {len(bc1_counts):,} rows")

    # Detect count column name (handles both 'counts' and 'sample_counts')
    if "counts" in bc1_counts.columns:
        count_col = "counts"
    elif "sample_counts" in bc1_counts.columns:
        count_col = "sample_counts"
    else:
        raise ValueError("BC1 counts file must have 'counts' or 'sample_counts' column")
    logger.info(f"  Using count column: {count_col}")

    # Load BC1 reference table
    logger.info("Loading BC1 reference table...")
    ref_table = pd.read_csv(args.bc1_reference, sep="\t")
    logger.info(f"  Loaded {len(ref_table):,} BC1s")

    # Load sample key
    logger.info("Loading sample key...")
    sample_key = pd.read_csv(args.sample_key, sep="\t")
    logger.info(f"  Loaded {len(sample_key):,} samples")

    # Verify required columns
    required_cols = ["sample_name", "condition", "generations"]
    missing = [c for c in required_cols if c not in sample_key.columns]
    if missing:
        raise ValueError(f"Sample key missing required columns: {missing}")

    # Configure tier thresholds
    tier_config = TierConfig(
        tier_1_min_reads=args.tier_1_min,
        tier_2_min_reads=args.tier_2_min,
    )

    # Calculate t0 counts
    logger.info(f"Calculating t0 counts (condition={args.t0_condition}, generations=0)...")
    t0_counts = calculate_t0_counts(
        bc1_counts,
        sample_key,
        t0_condition=args.t0_condition,
        screen_dates=args.screen_dates,
        count_col=count_col,
    )
    logger.info(f"  Found {len(t0_counts):,} BC1s with t0 counts")

    # Classify into tiers
    logger.info("Classifying BC1s into tiers...")
    tier_mapping = classify_bc1_tiers(t0_counts, tier_config)

    # Add oligo information from reference table
    if "bc1" in ref_table.columns and "oligo_name" in ref_table.columns:
        tier_mapping = tier_mapping.merge(
            ref_table[["bc1", "oligo_name"]].drop_duplicates(),
            on="bc1",
            how="left"
        )
        n_with_oligo = tier_mapping["oligo_name"].notna().sum()
        logger.info(f"  {n_with_oligo:,} BC1s mapped to oligos")

    # Create summary statistics
    tier_stats = create_tier_summary(tier_mapping, args.library)

    # Log summary
    for tier in ["abundance_high", "abundance_medium", "abundance_low"]:
        n_bc1s = len(tier_mapping[tier_mapping["tier"] == tier])
        logger.info(f"  {tier}: {n_bc1s:,} BC1s")

    # Ensure output directories exist
    args.output_tier_mapping.parent.mkdir(parents=True, exist_ok=True)
    args.output_stats.parent.mkdir(parents=True, exist_ok=True)

    # Save outputs
    logger.info(f"Saving tier mapping to {args.output_tier_mapping}")
    tier_mapping.to_csv(args.output_tier_mapping, sep="\t", index=False)

    logger.info(f"Saving tier statistics to {args.output_stats}")
    tier_stats.to_csv(args.output_stats, sep="\t", index=False)

    logger.info("Step 04: Abundance tier classification complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
