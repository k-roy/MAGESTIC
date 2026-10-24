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

# Shared with step 10 (sibling module free of magestic imports, like config.py).
try:
    from magestic.pipelines.bc1_from_screen.transformation_scope import is_design_pool_column
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from transformation_scope import is_design_pool_column

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
    parser.add_argument(
        "--sample-metadata", type=Path, default=None,
        help="Sample key giving each count column's transformation (checked before fitting)"
    )
    parser.add_argument(
        "--feature-metadata", type=Path, default=None,
        help="bc1 table giving each feature's transformation (checked before fitting)"
    )
    parser.add_argument(
        "--transformation-column", type=str, default=None,
        help="Construction-derived transformation column in the metadata, e.g. trafo_ID (no default)"
    )
    parser.add_argument(
        "--sample-name-column", type=str, default="sample_name",
        help="Sample-name column in --sample-metadata (default: sample_name)"
    )
    parser.add_argument(
        "--feature-id-column", type=str, default="bc1",
        help="Feature column in --feature-metadata (default: bc1)"
    )
    parser.add_argument(
        "--pool-transformations", action="store_true", default=False,
        help="Fit columns or features from several transformations together (explicit pooling)"
    )
    parser.add_argument(
        "--size-factors-fit-type", choices=("ratio", "poscounts", "iterative"), default="ratio",
        help=("pydeseq2 size-factor estimator: ratio (median of ratios; the pydeseq2 default, kept here), poscounts (DESeq2's estimator for zero-heavy counts) or iterative. Under ratio, pydeseq2 switches to the iterative estimator when every bc1 has a zero in some sample, which a screen matrix with a consolidated t0 always does; that fit ran for hours on one yL406 comparison (planning/272).")
    )
    return parser.parse_args()


def validate_transformation_scope(counts_df: pd.DataFrame, args) -> None:
    """Refuse a matrix whose columns or features span transformations, unless pooling is explicit.

    Without --sample-metadata/--feature-metadata the behaviour is unchanged (step 10 is
    the gate). The transformation column has no default, and a design-pool column
    (oligo_name, SPS, V library IDs...) is refused. Every count column (or fitted
    feature) must have exactly one transformation value; duplicate rows must agree.
    """
    axes = [
        (getattr(args, "sample_metadata", None), getattr(args, "sample_name_column", "sample_name"),
         [str(c) for c in counts_df.columns], "count column"),
        (getattr(args, "feature_metadata", None), getattr(args, "feature_id_column", "bc1"),
         [str(i) for i in counts_df.index], "feature"),
    ]
    if not any(path for path, *_ in axes):
        return
    column = getattr(args, "transformation_column", None)
    if not column:
        raise ValueError("--transformation-column is required with --sample-metadata or --feature-metadata; "
                         "name a construction-derived column such as trafo_ID.")
    if is_design_pool_column(column):
        raise ValueError(f"--transformation-column '{column}' identifies a design pool, not a transformation; "
                         f"name a construction-derived column such as trafo_ID.")
    for path, key, ids, label in axes:
        if not path:
            continue
        table = pd.read_csv(path, sep="\t")
        for needed in (key, column):
            if needed not in table.columns:
                raise ValueError(f"{path} lacks column '{needed}'; available: {list(table.columns)}")
        values = {}
        for ident, value in zip(table[key].astype(str), table[column]):
            if pd.isna(value) or str(value).strip() == "":
                continue
            values.setdefault(ident, set()).add(str(value).strip())
        conflicting = sorted(i for i in ids if len(values.get(i, ())) > 1)
        if conflicting:
            raise ValueError(f"{len(conflicting)} {label}(s) have conflicting '{column}' values in {path}: "
                             f"{conflicting[:10]}")
        uncovered = sorted(i for i in ids if i not in values)
        if uncovered:
            raise ValueError(f"{len(uncovered)} {label}(s) have no '{column}' value in {path}: {uncovered[:10]}")
        present = sorted({next(iter(values[i])) for i in ids})
        if len(present) > 1 and not getattr(args, "pool_transformations", False):
            raise ValueError(
                f"The {label}s span {len(present)} transformations in '{column}': {present}. Fit each "
                f"transformation separately, or pass --pool-transformations to fit them together explicitly.")


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


def size_factor_note(counts_df: pd.DataFrame, fit_type: str) -> str:
    """Say which size-factor estimator pydeseq2 will actually use (ratio silently becomes iterative)."""
    zero_free = int((counts_df > 0).all(axis=1).sum())
    if fit_type == "ratio" and zero_free == 0:
        return (f"size factors: ratio requested, but 0 of {len(counts_df):,} BC1s are zero-free, so pydeseq2 "
                f"switches to the ITERATIVE estimator (slow at screen scale; see --size-factors-fit-type)")
    return f"size factors: {fit_type} ({zero_free:,} of {len(counts_df):,} BC1s are zero-free)"


def run_deseq2(counts_df: pd.DataFrame, metadata_df: pd.DataFrame,
               size_factors_fit_type: str = "ratio") -> pd.DataFrame:
    """Run DESeq2 analysis."""
    # Ensure counts are integers
    counts_df = counts_df.round().astype(int)

    # Filter genes with zero counts across all samples
    counts_df = counts_df.loc[counts_df.sum(axis=1) > 0]

    logger.info(f"  Running DESeq2 on {len(counts_df):,} BC1s, {len(counts_df.columns)} samples")
    logger.info(f"  {size_factor_note(counts_df, size_factors_fit_type)}")

    # Create DESeq2 dataset
    # Note: PyDESeq2/AnnData expects counts with rows=samples, columns=features
    # Our input has rows=BC1s (features), columns=samples, so we transpose
    inference = DefaultInference(n_cpus=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))

    dds = DeseqDataSet(
        counts=counts_df.T,  # Transpose: samples as rows, BC1s as columns
        metadata=metadata_df,
        design_factors="condition",
        refit_cooks=True,
        size_factors_fit_type=size_factors_fit_type,
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

    # Refuse a matrix spanning transformations before any model is built.
    try:
        validate_transformation_scope(counts_df, args)
    except ValueError as exc:
        logger.error(str(exc))
        return 1

    # Run DESeq2
    logger.info("Running DESeq2...")
    results_df = run_deseq2(counts_df, metadata_df, args.size_factors_fit_type)

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
