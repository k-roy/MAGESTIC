#!/usr/bin/env python3
"""
Step 16: Rarefaction Analysis for Sequencing Depth Assessment

Performs rarefaction analysis to assess sequencing saturation by subsampling
reads at different depths and counting unique barcodes recovered.

Supports four data types:
- guide_donor_bc0: Guide-donor-BC0 demultiplexing results
- bc1_donor_bc0: BC1-donor-BC0 demultiplexing results
- bc1_t0: BC1 counts from t0 samples (baseline)
- bc1_timepoints: BC1 counts from treatment timepoints

Usage:
    python 16_rarefaction.py \\
        --input counts.tsv \\
        --output rarefaction_results.tsv \\
        --data-type bc1_t0 \\
        --sample-name my_sample \\
        --depths 10000000 50000000 100000000 200000000

    # Predict saturation at higher depths
    python 16_rarefaction.py \\
        --input counts.tsv \\
        --output rarefaction_results.tsv \\
        --data-type bc1_t0 \\
        --sample-name my_sample \\
        --predict-depths 200000000 500000000 1000000000

Author: Kevin R. Roy
"""

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import List, Tuple

# CRITICAL: Limit threads BEFORE importing numpy
n_threads = os.environ.get("SLURM_CPUS_PER_TASK", "4")
os.environ["OMP_NUM_THREADS"] = n_threads
os.environ["OPENBLAS_NUM_THREADS"] = n_threads

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Default rarefaction depths (in reads)
DEFAULT_DEPTHS = [
    1_000_000,     # 1M
    5_000_000,     # 5M
    10_000_000,    # 10M
    25_000_000,    # 25M
    50_000_000,    # 50M
    100_000_000,   # 100M
    200_000_000,   # 200M
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Rarefaction analysis for sequencing depth assessment"
    )
    parser.add_argument(
        "--input", type=Path, required=True,
        help="Input counts TSV (columns: barcode, count or sample columns)"
    )
    parser.add_argument(
        "--output", type=Path, required=True,
        help="Output rarefaction results TSV"
    )
    parser.add_argument(
        "--data-type", type=str, required=True,
        choices=["guide_donor_bc0", "bc1_donor_bc0", "bc1_t0", "bc1_timepoints"],
        help="Type of data being analyzed"
    )
    parser.add_argument(
        "--sample-name", type=str, required=True,
        help="Sample name for output labeling"
    )
    parser.add_argument(
        "--barcode-column", type=str, default=None,
        help="Column name for barcode (auto-detected if not specified)"
    )
    parser.add_argument(
        "--count-column", type=str, default=None,
        help="Column name for counts (auto-detected if not specified)"
    )
    parser.add_argument(
        "--depths", type=int, nargs="+", default=None,
        help="Depths to subsample (default: 1M to 200M)"
    )
    parser.add_argument(
        "--predict-depths", type=int, nargs="+", default=None,
        help="Additional depths to predict using extrapolation"
    )
    parser.add_argument(
        "--n-iterations", type=int, default=10,
        help="Number of subsampling iterations per depth (default: 10)"
    )
    parser.add_argument(
        "--min-count-threshold", type=int, default=1,
        help="Minimum count to consider a barcode detected (default: 1)"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for reproducibility (default: 42)"
    )
    return parser.parse_args()


def auto_detect_columns(df: pd.DataFrame, data_type: str) -> Tuple[str, str]:
    """Auto-detect barcode and count columns based on data type."""
    # Common barcode column names by data type
    barcode_names = {
        "guide_donor_bc0": ["guide_donor_bc0", "gRNA", "guide", "barcode"],
        "bc1_donor_bc0": ["bc1_donor_bc0", "bc1", "barcode"],
        "bc1_t0": ["bc1", "barcode"],
        "bc1_timepoints": ["bc1", "barcode"],
    }

    # Common count column names
    count_names = ["count", "counts", "n_reads", "reads", "total"]

    # Try to find barcode column
    barcode_col = None
    for name in barcode_names.get(data_type, ["barcode"]):
        if name in df.columns:
            barcode_col = name
            break

    # If still not found, use first column if it looks like a barcode
    if barcode_col is None and len(df.columns) > 0:
        first_col = df.columns[0]
        if df[first_col].dtype == object:
            barcode_col = first_col

    # Try to find count column
    count_col = None
    for name in count_names:
        if name in df.columns:
            count_col = name
            break

    # If count column not found, look for numeric column
    if count_col is None:
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) > 0:
            # Prefer column that isn't the barcode
            for col in numeric_cols:
                if col != barcode_col:
                    count_col = col
                    break

    return barcode_col, count_col


def load_counts(
    input_file: Path,
    data_type: str,
    barcode_column: str = None,
    count_column: str = None
) -> pd.Series:
    """
    Load counts from input file.

    Returns:
        Series with barcode as index and counts as values
    """
    df = pd.read_csv(input_file, sep="\t")

    # Auto-detect columns if not specified
    if barcode_column is None or count_column is None:
        auto_bc, auto_count = auto_detect_columns(df, data_type)
        barcode_column = barcode_column or auto_bc
        count_column = count_column or auto_count

    if barcode_column is None:
        raise ValueError(f"Could not detect barcode column. Columns: {list(df.columns)}")
    if count_column is None:
        raise ValueError(f"Could not detect count column. Columns: {list(df.columns)}")

    logger.info(f"  Barcode column: {barcode_column}")
    logger.info(f"  Count column: {count_column}")

    # Create Series
    counts = df.set_index(barcode_column)[count_column]
    return counts


def subsample_counts(counts: pd.Series, target_depth: int, rng: np.random.Generator) -> pd.Series:
    """
    Subsample counts to a target depth using multinomial sampling.

    Args:
        counts: Series with barcode counts
        target_depth: Target total read count
        rng: Random number generator

    Returns:
        Subsampled counts series
    """
    total_reads = counts.sum()

    if target_depth >= total_reads:
        return counts.copy()

    # Calculate probabilities for each barcode
    probs = counts.values / total_reads

    # Multinomial sampling
    subsampled = rng.multinomial(target_depth, probs)

    return pd.Series(subsampled, index=counts.index)


def calculate_rarefaction_metrics(
    counts: pd.Series,
    min_count_threshold: int = 1
) -> dict:
    """Calculate metrics for a count distribution."""
    detected = (counts >= min_count_threshold).sum()
    total_reads = counts.sum()

    # Count distribution metrics
    counts_above_threshold = counts[counts >= min_count_threshold]

    return {
        "n_unique_detected": detected,
        "total_reads": total_reads,
        "mean_count": counts_above_threshold.mean() if len(counts_above_threshold) > 0 else 0,
        "median_count": counts_above_threshold.median() if len(counts_above_threshold) > 0 else 0,
        "max_count": counts_above_threshold.max() if len(counts_above_threshold) > 0 else 0,
    }


def run_rarefaction(
    counts: pd.Series,
    depths: List[int],
    n_iterations: int,
    min_count_threshold: int,
    seed: int
) -> pd.DataFrame:
    """
    Run rarefaction analysis at multiple depths.

    Returns:
        DataFrame with rarefaction results
    """
    rng = np.random.default_rng(seed)
    total_reads = counts.sum()

    results = []

    # Add actual depth first
    actual_metrics = calculate_rarefaction_metrics(counts, min_count_threshold)
    results.append({
        "depth": total_reads,
        "iteration": 0,
        "is_actual": True,
        **actual_metrics,
    })

    # Run subsampling at each depth
    for depth in depths:
        if depth >= total_reads:
            logger.info(f"  Depth {depth:,}: exceeds total reads ({total_reads:,}), skipping")
            continue

        logger.info(f"  Depth {depth:,}: running {n_iterations} iterations")

        for iteration in range(n_iterations):
            subsampled = subsample_counts(counts, depth, rng)
            metrics = calculate_rarefaction_metrics(subsampled, min_count_threshold)

            results.append({
                "depth": depth,
                "iteration": iteration + 1,
                "is_actual": False,
                **metrics,
            })

    return pd.DataFrame(results)


def fit_saturation_model(rarefaction_df: pd.DataFrame) -> dict:
    """
    Fit a saturation model to predict unique barcodes at higher depths.

    Uses the Michaelis-Menten-like model: y = (a * x) / (b + x)
    where y = unique barcodes, x = depth

    Returns:
        Model parameters dict
    """
    # Use mean values at each depth
    summary = rarefaction_df.groupby("depth")["n_unique_detected"].mean().reset_index()
    x = summary["depth"].values
    y = summary["n_unique_detected"].values

    if len(x) < 3:
        logger.warning("  Not enough data points to fit saturation model")
        return None

    try:
        from scipy.optimize import curve_fit

        def michaelis_menten(x, a, b):
            return (a * x) / (b + x)

        # Initial guesses
        y_max_est = y.max() * 2
        k_est = x[np.argmax(y > y.max() * 0.5)] if len(x) > 1 else x[0]

        popt, pcov = curve_fit(
            michaelis_menten,
            x,
            y,
            p0=[y_max_est, k_est],
            bounds=([0, 0], [np.inf, np.inf]),
            maxfev=10000
        )

        return {
            "a_max_unique": popt[0],
            "b_half_saturation": popt[1],
            "model": "michaelis_menten",
        }

    except Exception as e:
        logger.warning(f"  Could not fit saturation model: {e}")
        return None


def predict_at_depth(model_params: dict, depth: int) -> float:
    """Predict unique barcodes at a given depth using fitted model."""
    if model_params is None or model_params.get("model") != "michaelis_menten":
        return np.nan

    a = model_params["a_max_unique"]
    b = model_params["b_half_saturation"]

    return (a * depth) / (b + depth)


def main():
    args = parse_args()

    logger.info(f"Step 16: Rarefaction analysis for {args.sample_name}")
    logger.info(f"  Data type: {args.data_type}")
    logger.info(f"  Input: {args.input}")

    # Load counts
    logger.info("Loading counts...")
    counts = load_counts(
        args.input,
        args.data_type,
        args.barcode_column,
        args.count_column
    )

    total_reads = counts.sum()
    n_unique = len(counts)
    logger.info(f"  Total reads: {total_reads:,}")
    logger.info(f"  Unique barcodes: {n_unique:,}")

    # Determine depths to analyze
    depths = args.depths or DEFAULT_DEPTHS
    depths = sorted([d for d in depths if d < total_reads])

    if not depths:
        logger.warning("No depths below total reads - adding 50% and 25% depths")
        depths = [int(total_reads * 0.25), int(total_reads * 0.5)]

    logger.info(f"  Analyzing depths: {[f'{d:,}' for d in depths]}")

    # Run rarefaction
    logger.info("Running rarefaction analysis...")
    rarefaction_df = run_rarefaction(
        counts,
        depths,
        args.n_iterations,
        args.min_count_threshold,
        args.seed
    )

    # Add metadata
    rarefaction_df["sample_name"] = args.sample_name
    rarefaction_df["data_type"] = args.data_type

    # Fit saturation model
    logger.info("Fitting saturation model...")
    model_params = fit_saturation_model(rarefaction_df)

    if model_params:
        logger.info(f"  Max unique (asymptote): {model_params['a_max_unique']:,.0f}")
        logger.info(f"  Half-saturation depth: {model_params['b_half_saturation']:,.0f}")

        # Calculate percent saturation at actual depth
        predicted_at_actual = predict_at_depth(model_params, total_reads)
        percent_saturation = 100 * (counts >= args.min_count_threshold).sum() / model_params['a_max_unique']
        logger.info(f"  Percent saturation at current depth: {percent_saturation:.1f}%")

        # Add model predictions
        rarefaction_df["predicted_max_unique"] = model_params["a_max_unique"]
        rarefaction_df["predicted_half_saturation"] = model_params["b_half_saturation"]
        rarefaction_df["percent_saturation"] = 100 * rarefaction_df["n_unique_detected"] / model_params["a_max_unique"]

    # Predict at additional depths if requested
    if args.predict_depths and model_params:
        logger.info("Predicting at additional depths...")
        for depth in args.predict_depths:
            predicted = predict_at_depth(model_params, depth)
            logger.info(f"  At {depth:,} reads: {predicted:,.0f} unique barcodes predicted")

            # Add prediction row
            pred_row = {
                "depth": depth,
                "iteration": 0,
                "is_actual": False,
                "n_unique_detected": predicted,
                "total_reads": depth,
                "mean_count": np.nan,
                "median_count": np.nan,
                "max_count": np.nan,
                "sample_name": args.sample_name,
                "data_type": args.data_type,
                "predicted_max_unique": model_params["a_max_unique"],
                "predicted_half_saturation": model_params["b_half_saturation"],
                "percent_saturation": 100 * predicted / model_params["a_max_unique"],
                "is_prediction": True,
            }
            rarefaction_df = pd.concat([rarefaction_df, pd.DataFrame([pred_row])], ignore_index=True)

    # Fill is_prediction column
    if "is_prediction" not in rarefaction_df.columns:
        rarefaction_df["is_prediction"] = False

    # Create summary statistics
    logger.info("Creating summary...")
    summary = rarefaction_df[~rarefaction_df["is_prediction"]].groupby("depth").agg({
        "n_unique_detected": ["mean", "std", "min", "max"],
        "total_reads": "first",
    }).round(2)
    summary.columns = ["_".join(col).strip("_") for col in summary.columns]

    for _, row in summary.iterrows():
        depth = int(row.name)
        mean_unique = row["n_unique_detected_mean"]
        std_unique = row.get("n_unique_detected_std", 0)
        logger.info(f"  Depth {depth:,}: {mean_unique:,.0f} ± {std_unique:,.0f} unique barcodes")

    # Ensure output directory exists
    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Save results
    logger.info(f"Saving results to {args.output}")
    rarefaction_df.to_csv(args.output, sep="\t", index=False)

    # Also save summary
    summary_file = args.output.with_suffix(".summary.tsv")
    summary.to_csv(summary_file, sep="\t")
    logger.info(f"Saved summary to {summary_file}")

    logger.info("Step 16: Rarefaction analysis complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
