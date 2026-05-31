#!/usr/bin/env python3
"""
Step 11: Run DESeq2 for a Single Comparison

Runs PyDESeq2 differential expression analysis for one condition_generation comparison.
Each library is processed separately to ensure proper normalization.

Usage:
    python 11_deseq2_single.py \\
        --reference reference_counts.tsv \\
        --treatment treatment_counts.tsv \\
        --output results.tsv \\
        --comparison Caffeine_6 \\
        --library yL437

Author: Kevin R. Roy
"""

import argparse
import logging
import os
import sys
import warnings
from pathlib import Path

# CRITICAL: Limit threads BEFORE importing numpy/scipy/pydeseq2
# This prevents joblib from spawning too many workers on SLURM
n_threads = os.environ.get("SLURM_CPUS_PER_TASK", "4")
os.environ["OMP_NUM_THREADS"] = n_threads
os.environ["OPENBLAS_NUM_THREADS"] = n_threads
os.environ["MKL_NUM_THREADS"] = n_threads
os.environ["LOKY_MAX_CPU_COUNT"] = n_threads

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=UserWarning)

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run DESeq2 for a single comparison"
    )
    parser.add_argument(
        "--reference", type=Path, required=True,
        help="Reference (t0) count table"
    )
    parser.add_argument(
        "--treatment", type=Path, required=True,
        help="Treatment count table"
    )
    parser.add_argument(
        "--output", type=Path, required=True,
        help="Output DESeq2 results TSV"
    )
    parser.add_argument(
        "--comparison", type=str, required=True,
        help="Comparison ID (e.g., Caffeine_6)"
    )
    parser.add_argument(
        "--library", type=str, required=True,
        help="Library name (e.g., yL437)"
    )
    parser.add_argument(
        "--min-basemean", type=float, default=0,
        help="Minimum baseMean filter (default: 0, no filtering)"
    )
    return parser.parse_args()


def load_and_combine_counts(reference_file: Path, treatment_file: Path) -> tuple:
    """
    Load and combine reference + treatment count data.

    Returns:
        (counts_df, metadata_df)
    """
    ref_df = pd.read_csv(reference_file, sep="\t", index_col=0)
    treat_df = pd.read_csv(treatment_file, sep="\t", index_col=0)

    # Combine counts (ensure same BC1s)
    common_bc1s = ref_df.index.intersection(treat_df.index)
    ref_df = ref_df.loc[common_bc1s]
    treat_df = treat_df.loc[common_bc1s]

    counts_df = pd.concat([ref_df, treat_df], axis=1)

    # Create metadata
    ref_samples = ref_df.columns.tolist()
    treat_samples = treat_df.columns.tolist()

    metadata_df = pd.DataFrame({
        "sample_ID": ref_samples + treat_samples,
        "condition": ["reference"] * len(ref_samples) + ["treatment"] * len(treat_samples)
    }).set_index("sample_ID")

    return counts_df, metadata_df


def run_deseq2(counts_df: pd.DataFrame, metadata_df: pd.DataFrame) -> pd.DataFrame:
    """Run DESeq2 analysis."""
    # Ensure counts are integers
    counts_df = counts_df.round().astype(int)

    # Filter genes with zero counts across all samples
    counts_df = counts_df.loc[counts_df.sum(axis=1) > 0]

    logger.info(f"  Running DESeq2 on {len(counts_df):,} BC1s, {len(counts_df.columns)} samples")

    # Create DESeq2 dataset
    # Note: PyDESeq2/AnnData expects counts with rows=samples, columns=features
    # Our input has rows=BC1s (features), columns=samples, so we transpose
    inference = DefaultInference(n_cpus=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))

    dds = DeseqDataSet(
        counts=counts_df.T,  # Transpose: samples as rows, BC1s as columns
        metadata=metadata_df,
        design_factors="condition",
        refit_cooks=True,
        inference=inference,
    )

    # Run DESeq2 workflow
    dds.deseq2()

    # Get results (treatment vs reference)
    stats = DeseqStats(dds, contrast=["condition", "treatment", "reference"], inference=inference)
    stats.summary()

    # Extract results
    results_df = stats.results_df.copy()
    results_df.index.name = "bc1"

    return results_df


def main():
    args = parse_args()

    logger.info(f"Step 11: DESeq2 for {args.library} - {args.comparison}")
    logger.info(f"  Reference: {args.reference}")
    logger.info(f"  Treatment: {args.treatment}")
    logger.info(f"  Thread limit: {n_threads}")

    # Load data
    logger.info("Loading count data...")
    counts_df, metadata_df = load_and_combine_counts(args.reference, args.treatment)
    logger.info(f"  Loaded {len(counts_df):,} BC1s")
    logger.info(f"  Reference samples: {(metadata_df['condition'] == 'reference').sum()}")
    logger.info(f"  Treatment samples: {(metadata_df['condition'] == 'treatment').sum()}")

    # Run DESeq2
    logger.info("Running DESeq2...")
    results_df = run_deseq2(counts_df, metadata_df)

    # Add metadata columns
    results_df["library"] = args.library
    results_df["comparison"] = args.comparison

    # Apply baseMean filter if specified
    if args.min_basemean > 0:
        n_before = len(results_df)
        results_df = results_df[results_df["baseMean"] >= args.min_basemean]
        logger.info(f"  Filtered to baseMean >= {args.min_basemean}: {n_before:,} -> {len(results_df):,}")

    # Log summary statistics
    n_sig = (results_df["padj"] < 0.1).sum()
    logger.info(f"  Results: {len(results_df):,} BC1s, {n_sig:,} significant (padj < 0.1)")

    # Ensure output directory exists
    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Save results
    logger.info(f"Saving results to {args.output}")
    results_df.to_csv(args.output, sep="\t")

    logger.info("Step 11: DESeq2 complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
