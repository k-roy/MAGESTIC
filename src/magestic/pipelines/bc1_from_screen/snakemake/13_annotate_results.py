#!/usr/bin/env python3
"""
Step 13: Annotate DESeq2 Results with Variant Information

Merges DESeq2 results with BC1 reference table and oligo design to add:
- Tier classification
- Oligo name and variant information
- Gene annotations
- Codon variant type (missense, synonymous, nonsense)

Usage:
    python 13_annotate_results.py \\
        --deseq2-results deseq2_results.tsv \\
        --ref-table bc1_reference_table.tsv \\
        --tier-mapping tier_mapping.tsv \\
        --oligo-design oligo_design.tsv \\
        --output annotated.tsv \\
        --library yL437

Author: Kevin R. Roy
"""

import argparse
import logging
import sys
from pathlib import Path

# Add common scripts to path

import pandas as pd

from magestic.pipelines.bc1_from_screen.core.annotation import (
    create_annotation_lookup,
    annotate_deseq2_results,
)

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Annotate DESeq2 results with variant information"
    )
    parser.add_argument(
        "--deseq2-results", type=Path, required=True,
        help="Input DESeq2 results TSV"
    )
    parser.add_argument(
        "--ref-table", type=Path, required=True,
        help="BC1 reference table with oligo assignments"
    )
    parser.add_argument(
        "--tier-mapping", type=Path, required=True,
        help="BC1 tier mapping TSV"
    )
    parser.add_argument(
        "--oligo-design", type=Path, required=True,
        help="Oligo design TSV with variant annotations"
    )
    parser.add_argument(
        "--output", type=Path, required=True,
        help="Output annotated TSV"
    )
    parser.add_argument(
        "--library", type=str, required=True,
        help="Library name"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    logger.info(f"Step 13: Annotating results for library {args.library}")
    logger.info(f"  DESeq2 results: {args.deseq2_results}")
    logger.info(f"  Reference table: {args.ref_table}")
    logger.info(f"  Tier mapping: {args.tier_mapping}")
    logger.info(f"  Oligo design: {args.oligo_design}")

    # Load DESeq2 results
    logger.info("Loading DESeq2 results...")
    deseq2_df = pd.read_csv(args.deseq2_results, sep="\t")
    if deseq2_df.columns[0] == "Unnamed: 0" or "bc1" not in deseq2_df.columns:
        # First column is likely bc1 index
        deseq2_df = pd.read_csv(args.deseq2_results, sep="\t", index_col=0)
        deseq2_df = deseq2_df.reset_index()
        deseq2_df = deseq2_df.rename(columns={deseq2_df.columns[0]: "bc1"})
    logger.info(f"  Loaded {len(deseq2_df):,} BC1 results")

    # Load reference table
    logger.info("Loading reference table...")
    ref_table = pd.read_csv(args.ref_table, sep="\t")
    logger.info(f"  Loaded {len(ref_table):,} BC1 entries")

    # Load tier mapping
    logger.info("Loading tier mapping...")
    tier_mapping = pd.read_csv(args.tier_mapping, sep="\t")
    logger.info(f"  Loaded {len(tier_mapping):,} tier assignments")

    # Load oligo design
    logger.info("Loading oligo design...")
    oligo_design = pd.read_csv(args.oligo_design, sep="\t")
    logger.info(f"  Loaded {len(oligo_design):,} oligo designs")

    # Create annotation lookup
    logger.info("Creating annotation lookup...")
    annotation_lookup = create_annotation_lookup(ref_table, oligo_design)

    # Merge tier information
    logger.info("Merging tier information...")
    if "tier" in tier_mapping.columns:
        tier_cols = ["bc1", "tier"]
        if "t0_counts" in tier_mapping.columns:
            tier_cols.append("t0_counts")
        deseq2_df = deseq2_df.merge(
            tier_mapping[tier_cols].drop_duplicates(),
            on="bc1",
            how="left"
        )
        deseq2_df["tier"] = deseq2_df["tier"].fillna("unknown")

    # Annotate results
    logger.info("Annotating results...")
    annotated_df = annotate_deseq2_results(deseq2_df, annotation_lookup)

    # Add library column
    annotated_df["library"] = args.library

    # Log summary
    n_annotated = annotated_df["oligo_name"].notna().sum()
    logger.info(f"  Annotated {n_annotated:,} / {len(annotated_df):,} BC1s ({100*n_annotated/len(annotated_df):.1f}%)")

    if "tier" in annotated_df.columns:
        for tier in annotated_df["tier"].unique():
            n_tier = (annotated_df["tier"] == tier).sum()
            logger.info(f"    {tier}: {n_tier:,} BC1s")

    if "codon_variant_type" in annotated_df.columns:
        for vtype in annotated_df["codon_variant_type"].dropna().unique()[:5]:
            n_vtype = (annotated_df["codon_variant_type"] == vtype).sum()
            logger.info(f"    {vtype}: {n_vtype:,}")

    # Ensure output directory exists
    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Save output
    logger.info(f"Saving annotated results to {args.output}")
    annotated_df.to_csv(args.output, sep="\t", index=False)

    logger.info("Step 13: Annotation complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
