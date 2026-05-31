#!/usr/bin/env python3
"""
Step D: Create DESeq2 Count Matrix with Consolidated T0

Creates a wide-format count matrix combining tier_1 and tier_2 BC1s for a library,
with optional t0 sample consolidation for robust normalization.

Features:
- Library-specific processing
- Consolidated t0 (combines technical replicates)
- Includes both tier_1 (individual) and tier_2 (aggregated) BC1s

Usage:
    python step_d_create_count_matrix.py \\
        --bc1-counts bc1_counts_long.tsv \\
        --tier-mapping tier_mapping.tsv \\
        --sample-key sample_key.tsv \\
        --library yL437 \\
        --output count_matrix.tsv

Author: Kevin R. Roy
"""

import argparse
import logging
import sys
from pathlib import Path

# Add common scripts to path
sys.path.insert(0, str(Path(__file__).parents[2]))

import pandas as pd
import numpy as np

from bc1_from_screen_pipeline.core.count_matrix import (
    create_wide_count_matrix,
    consolidate_t0_samples,
)

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create wide-format count matrix for DESeq2"
    )
    parser.add_argument(
        "--bc1-counts", type=Path, required=True,
        help="Input BC1 counts long-format TSV"
    )
    parser.add_argument(
        "--tier-mapping", type=Path, required=True,
        help="BC1 tier mapping TSV"
    )
    parser.add_argument(
        "--sample-key", type=Path, required=True,
        help="Sample key TSV with condition/generation info"
    )
    parser.add_argument(
        "--library", type=str, required=True,
        help="Library name to process"
    )
    parser.add_argument(
        "--output", type=Path, required=True,
        help="Output wide-format count matrix TSV"
    )
    parser.add_argument(
        "--consolidate-t0", type=str, default="True",
        help="Consolidate t0 samples (default: True)"
    )
    parser.add_argument(
        "--t0-groups", type=int, default=4,
        help="Number of t0 groups after consolidation (default: 4)"
    )
    parser.add_argument(
        "--screen-dates", type=str, nargs="+",
        help="Screen dates for this library"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    consolidate = args.consolidate_t0.lower() in ("true", "yes", "1")

    logger.info(f"Step D: Creating count matrix for library {args.library}")
    logger.info(f"  BC1 counts: {args.bc1_counts}")
    logger.info(f"  Tier mapping: {args.tier_mapping}")
    logger.info(f"  Sample key: {args.sample_key}")
    logger.info(f"  Consolidate t0: {consolidate} (groups: {args.t0_groups})")

    # Load data
    logger.info("Loading BC1 counts...")
    bc1_counts = pd.read_csv(args.bc1_counts, sep="\t")
    logger.info(f"  Loaded {len(bc1_counts):,} rows")

    logger.info("Loading tier mapping...")
    tier_mapping = pd.read_csv(args.tier_mapping, sep="\t")
    logger.info(f"  Loaded {len(tier_mapping):,} BC1 tier assignments")

    logger.info("Loading sample key...")
    sample_key = pd.read_csv(args.sample_key, sep="\t")
    logger.info(f"  Loaded {len(sample_key):,} sample entries")

    # Filter to relevant tiers (tier_1 and tier_2)
    relevant_tiers = ["tier_1", "tier_2"]
    tier_bc1s = tier_mapping[tier_mapping["tier"].isin(relevant_tiers)]["bc1"].unique()
    logger.info(f"  BC1s in tier_1/tier_2: {len(tier_bc1s):,}")

    # Filter BC1 counts to relevant BC1s
    bc1_counts = bc1_counts[bc1_counts["bc1"].isin(tier_bc1s)]
    logger.info(f"  Filtered counts: {len(bc1_counts):,} rows")

    # Filter to library-specific samples if screen_dates provided
    if args.screen_dates:
        # Filter sample_key and bc1_counts to these dates
        date_pattern = "|".join(args.screen_dates)
        if "sample_name" in bc1_counts.columns:
            bc1_counts = bc1_counts[bc1_counts["sample_name"].str.contains(date_pattern, na=False)]
        logger.info(f"  After screen date filter: {len(bc1_counts):,} rows")

    # Create wide-format matrix
    logger.info("Creating wide-format count matrix...")
    if "sample_name" in bc1_counts.columns and "counts" in bc1_counts.columns:
        wide_matrix = bc1_counts.pivot_table(
            index="bc1",
            columns="sample_name",
            values="counts",
            aggfunc="sum",
            fill_value=0
        )
    else:
        # Assume already wide format or use core function
        wide_matrix = create_wide_count_matrix(bc1_counts, tier_mapping)

    logger.info(f"  Matrix shape: {wide_matrix.shape[0]:,} BC1s x {wide_matrix.shape[1]:,} samples")

    # Consolidate t0 samples if requested
    if consolidate:
        logger.info(f"Consolidating t0 samples into {args.t0_groups} groups...")
        # Identify t0 samples (generation 0)
        t0_samples = [col for col in wide_matrix.columns if "_t0_" in col.lower() or "_gen0_" in col.lower() or col.endswith("_0")]

        if t0_samples:
            wide_matrix = consolidate_t0_samples(wide_matrix, t0_samples, args.t0_groups)
            logger.info(f"  Consolidated matrix: {wide_matrix.shape[0]:,} BC1s x {wide_matrix.shape[1]:,} samples")
        else:
            logger.warning("  No t0 samples found to consolidate")

    # Add tier information to output
    tier_info = tier_mapping[tier_mapping["bc1"].isin(wide_matrix.index)][["bc1", "tier"]].drop_duplicates()

    # Reset index for output
    wide_matrix = wide_matrix.reset_index()

    # Merge tier info
    wide_matrix = wide_matrix.merge(tier_info, on="bc1", how="left")

    # Reorder columns: bc1, tier, then samples
    cols = ["bc1", "tier"] + [c for c in wide_matrix.columns if c not in ["bc1", "tier"]]
    wide_matrix = wide_matrix[cols]

    # Ensure output directory exists
    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Save output
    logger.info(f"Saving count matrix to {args.output}")
    wide_matrix.to_csv(args.output, sep="\t", index=False)

    # Log summary by tier
    for tier in ["tier_1", "tier_2"]:
        n_tier = (wide_matrix["tier"] == tier).sum()
        logger.info(f"  {tier}: {n_tier:,} BC1s")

    logger.info("Step D: Count matrix creation complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
