#!/usr/bin/env python3
"""
Step 13b: DESeq2 Rarefaction Analysis

Performs rarefaction analysis to assess how DESeq2 hit calling is affected by
sequencing depth. Randomly subsamples read counts at different fractions and
runs DESeq2 analysis at each depth.

This helps determine:
1. Sequencing saturation for hit calling
2. Minimum depth needed for robust hit detection
3. How hit sensitivity/specificity changes with depth

Usage:
    python 13b_deseq_rarefaction.py \
        --reference reference_counts.tsv \
        --treatment treatment_counts.tsv \
        --output-dir rarefaction_results/ \
        --comparison Caffeine_21 \
        --library yL442 \
        --fractions 0.1 0.25 0.5 0.75 1.0 \
        --replicates 5

Author: Kevin R. Roy
"""

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple

# CRITICAL: Limit threads BEFORE importing numpy/scipy/pydeseq2
n_threads = os.environ.get("SLURM_CPUS_PER_TASK", "4")
os.environ["OMP_NUM_THREADS"] = n_threads
os.environ["OPENBLAS_NUM_THREADS"] = n_threads
os.environ["MKL_NUM_THREADS"] = n_threads
os.environ["LOKY_MAX_CPU_COUNT"] = n_threads

import numpy as np
import pandas as pd
import warnings

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
        description="DESeq2 rarefaction analysis for sequencing depth assessment"
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
        "--output-dir", type=Path, required=True,
        help="Output directory for rarefaction results"
    )
    parser.add_argument(
        "--comparison", type=str, required=True,
        help="Comparison name (e.g., Caffeine_21)"
    )
    parser.add_argument(
        "--library", type=str, required=True,
        help="Library name (e.g., yL442)"
    )
    parser.add_argument(
        "--fractions", type=float, nargs="+",
        default=[0.1, 0.25, 0.5, 0.75, 1.0],
        help="Subsampling fractions (default: 0.1 0.25 0.5 0.75 1.0)"
    )
    parser.add_argument(
        "--replicates", type=int, default=5,
        help="Number of replicates per fraction (default: 5)"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for reproducibility"
    )
    parser.add_argument(
        "--padj-threshold", type=float, default=0.05,
        help="Adjusted p-value threshold for hit calling (default: 0.05)"
    )
    parser.add_argument(
        "--lfc-threshold", type=float, default=0.5,
        help="Log2 fold change threshold for hit calling (default: 0.5)"
    )
    return parser.parse_args()


def subsample_counts(counts: pd.DataFrame, fraction: float, rng: np.random.Generator) -> pd.DataFrame:
    """
    Subsample counts by randomly removing reads.

    For each cell, if count is n, we sample from Binomial(n, fraction) to get
    the subsampled count. This preserves the expected proportion while adding
    sampling noise.
    """
    # Get numeric columns only
    bc1_col = counts.columns[0]
    sample_cols = counts.columns[1:]

    # Subsample each cell
    count_values = counts[sample_cols].values.astype(float)
    subsampled = rng.binomial(count_values.astype(int), fraction)

    result = counts.copy()
    result[sample_cols] = subsampled

    return result


def run_deseq2(ref_counts: pd.DataFrame, treat_counts: pd.DataFrame,
               comparison: str) -> pd.DataFrame:
    """
    Run DESeq2 analysis on subsampled counts.
    Returns results DataFrame with log2FoldChange, pvalue, padj.
    """
    # Combine reference and treatment
    bc1_col = ref_counts.columns[0]
    ref_samples = list(ref_counts.columns[1:])
    treat_samples = list(treat_counts.columns[1:])

    # Merge on BC1
    merged = ref_counts.merge(treat_counts, on=bc1_col)
    bc1_index = merged[bc1_col].values

    # Create count matrix (samples as rows, BC1s as columns for pydeseq2)
    all_samples = ref_samples + treat_samples
    count_matrix = merged[all_samples].T
    count_matrix.columns = bc1_index

    # Filter low-count BC1s (need at least 10 reads total)
    min_counts = count_matrix.sum(axis=0) >= 10
    count_matrix = count_matrix.loc[:, min_counts]

    if count_matrix.shape[1] < 100:
        logger.warning(f"  Too few BC1s after filtering: {count_matrix.shape[1]}")
        return None

    # Create metadata
    metadata = pd.DataFrame({
        "sample": all_samples,
        "condition": ["reference"] * len(ref_samples) + ["treatment"] * len(treat_samples)
    })
    metadata = metadata.set_index("sample")

    # Run DESeq2
    try:
        inference = DefaultInference(n_cpus=int(n_threads))
        dds = DeseqDataSet(
            counts=count_matrix,
            metadata=metadata,
            design_factors="condition",
            refit_cooks=True,
            inference=inference
        )
        dds.deseq2()

        stat_res = DeseqStats(dds, contrast=["condition", "treatment", "reference"], inference=inference)
        stat_res.summary()

        results = stat_res.results_df
        results["bc1"] = results.index
        results = results.reset_index(drop=True)

        return results

    except Exception as e:
        logger.warning(f"  DESeq2 error: {e}")
        return None


def call_hits(results: pd.DataFrame, padj_threshold: float, lfc_threshold: float) -> Dict:
    """
    Call hits based on adjusted p-value and log2 fold change thresholds.
    Returns summary statistics.
    """
    if results is None:
        return {
            "n_tested": 0,
            "n_significant": 0,
            "n_up": 0,
            "n_down": 0
        }

    n_tested = len(results)

    # Significant hits
    sig_mask = (results["padj"] < padj_threshold) & (results["padj"].notna())
    sig = results[sig_mask]

    n_significant = len(sig)
    n_up = len(sig[sig["log2FoldChange"] > lfc_threshold])
    n_down = len(sig[sig["log2FoldChange"] < -lfc_threshold])

    return {
        "n_tested": n_tested,
        "n_significant": n_significant,
        "n_up": n_up,
        "n_down": n_down
    }


def main():
    args = parse_args()

    logger.info(f"DESeq2 Rarefaction Analysis")
    logger.info(f"  Comparison: {args.comparison}")
    logger.info(f"  Library: {args.library}")
    logger.info(f"  Fractions: {args.fractions}")
    logger.info(f"  Replicates per fraction: {args.replicates}")

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    logger.info("Loading count tables...")
    ref_counts = pd.read_csv(args.reference, sep="\t")
    treat_counts = pd.read_csv(args.treatment, sep="\t")

    logger.info(f"  Reference: {ref_counts.shape[1]-1} samples")
    logger.info(f"  Treatment: {treat_counts.shape[1]-1} samples")

    # Calculate total reads
    total_ref_reads = ref_counts.iloc[:, 1:].sum().sum()
    total_treat_reads = treat_counts.iloc[:, 1:].sum().sum()
    total_reads = total_ref_reads + total_treat_reads
    logger.info(f"  Total reads: {total_reads:,}")

    # Initialize RNG
    rng = np.random.default_rng(args.seed)

    # Run rarefaction
    results_rows = []

    for fraction in sorted(args.fractions):
        depth = int(total_reads * fraction)
        logger.info(f"Fraction {fraction} ({depth:,} reads):")

        for rep in range(args.replicates):
            # Subsample
            ref_sub = subsample_counts(ref_counts, fraction, rng)
            treat_sub = subsample_counts(treat_counts, fraction, rng)

            actual_reads = ref_sub.iloc[:, 1:].sum().sum() + treat_sub.iloc[:, 1:].sum().sum()

            # Run DESeq2
            deseq_results = run_deseq2(ref_sub, treat_sub, args.comparison)

            # Call hits
            hits = call_hits(deseq_results, args.padj_threshold, args.lfc_threshold)

            logger.info(f"  Rep {rep+1}: {hits['n_significant']} hits "
                       f"({hits['n_up']} up, {hits['n_down']} down)")

            results_rows.append({
                "library": args.library,
                "comparison": args.comparison,
                "fraction": fraction,
                "target_depth": depth,
                "actual_depth": int(actual_reads),
                "replicate": rep + 1,
                "n_tested": hits["n_tested"],
                "n_significant": hits["n_significant"],
                "n_up": hits["n_up"],
                "n_down": hits["n_down"],
                "padj_threshold": args.padj_threshold,
                "lfc_threshold": args.lfc_threshold
            })

            # Save individual DESeq results for full depth
            if fraction == 1.0 and rep == 0 and deseq_results is not None:
                full_results_path = args.output_dir / f"{args.comparison}_full_depth_results.tsv"
                deseq_results.to_csv(full_results_path, sep="\t", index=False)

    # Save rarefaction summary
    results_df = pd.DataFrame(results_rows)
    summary_path = args.output_dir / f"{args.comparison}_rarefaction_summary.tsv"
    results_df.to_csv(summary_path, sep="\t", index=False)
    logger.info(f"Saved rarefaction summary to {summary_path}")

    # Calculate summary statistics across replicates
    summary_stats = results_df.groupby("fraction").agg({
        "n_significant": ["mean", "std"],
        "n_up": ["mean", "std"],
        "n_down": ["mean", "std"],
        "actual_depth": "mean"
    }).reset_index()

    summary_stats.columns = ["fraction", "mean_significant", "std_significant",
                             "mean_up", "std_up", "mean_down", "std_down", "mean_depth"]

    stats_path = args.output_dir / f"{args.comparison}_rarefaction_stats.tsv"
    summary_stats.to_csv(stats_path, sep="\t", index=False)
    logger.info(f"Saved summary statistics to {stats_path}")

    logger.info("Rarefaction analysis complete")


if __name__ == "__main__":
    main()
