#!/usr/bin/env python3
"""
Step 04: Calculate BC0 Purity for Guide-Donor-BC0 Combinations

For each unique guide-bc0 combination, calculates the purity ratio:
the fraction of reads supporting the dominant donor sequence vs all reads
with that guide-bc0 pair.

This identifies guide-bc0 combinations that have contamination from multiple
donor sequences (low purity) vs clean single-donor combinations (high purity).

Input: Matched guide-donor-bc0 counts from Step 03
Output: Purity-filtered files with guide_bc0_purity column

Usage:
    python 04_calculate_bc0_purity.py
"""

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import List
import sys

import pandas as pd
from tqdm import tqdm

# ============================================================================
# Configuration - Read from environment variables set by SLURM script
# ============================================================================

import os

def get_config_from_env():
    """Get project configuration from environment variables."""
    project_dir_str = os.environ.get("PROJECT_DIR")
    if not project_dir_str:
        raise ValueError(
            "PROJECT_DIR environment variable not set. "
            "Run via SLURM script or set: export PROJECT_DIR=/path/to/project"
        )

    project_dir = Path(project_dir_str)
    datatype = os.environ.get("DATATYPE", "guide_donor_bc0")

    processed_data_dir = project_dir / "processed_data" / datatype

    return {
        "project_dir": project_dir,
        "processed_data_dir": processed_data_dir,
        "counts_dir": processed_data_dir / "guide_donor_bc0",
        "purity_dir": processed_data_dir / "guide_donor_bc0_purity_filtered",
        "log_dir": processed_data_dir / "log",
    }

# File patterns
INPUT_SUFFIX = "_matched_guide_donor_bc0_counts.tsv"
OUTPUT_SUFFIX = "_guide_donor_bc0_purity_filtered.tsv"


@dataclass
class PurityConfig:
    """Configuration for purity calculation."""
    # Minimum counts to include in analysis
    min_counts: int = 1
    # Library ID prefixes to process (e.g., "V" for V629, V631, etc.)
    library_prefix: str = "V"


def setup_directories(env_config: dict) -> None:
    """Create necessary output directories."""
    for directory in [env_config["purity_dir"], env_config["log_dir"]]:
        directory.mkdir(parents=True, exist_ok=True)


def calculate_bc0_purity(
    input_file: Path,
    output_file: Path,
    config: PurityConfig = PurityConfig()
) -> dict:
    """
    Calculate BC0 purity for a single sample file.

    For each guide-bc0 combination:
    1. Sum all counts across different donors
    2. Find the top (highest count) donor
    3. Calculate purity = top_count / total_count

    Args:
        input_file: Path to matched guide-donor-bc0 counts file
        output_file: Path to write purity-filtered output
        config: Purity calculation configuration

    Returns:
        Dict with processing statistics
    """
    df = pd.read_csv(input_file, sep="\t")

    # Filter out low count entries if needed
    if config.min_counts > 1:
        df = df[df["counts"] >= config.min_counts]

    if df.empty:
        return {
            "total_rows": 0,
            "unique_guide_bc0": 0,
            "mean_purity": None,
        }

    # Sum counts for each unique guide-bc0 combination
    summed_guide_bc0 = (
        df.groupby(["guide", "bc0"])["counts"]
        .sum()
        .reset_index()
        .rename(columns={"counts": "grouped_guide_bc0_counts"})
    )

    # Retain only the top guide-bc0 combination (highest count per group)
    idx = df.groupby(["guide", "bc0"])["counts"].idxmax()
    top_guide_bc0 = df.loc[idx].copy()

    # Merge summed counts with top combinations
    merged = top_guide_bc0.merge(
        summed_guide_bc0,
        on=["guide", "bc0"]
    )

    # Calculate purity: ratio of top count to total count for guide-bc0
    merged["guide_bc0_purity"] = (
        merged["counts"] / merged["grouped_guide_bc0_counts"]
    )

    # Save results
    merged.to_csv(output_file, sep="\t", index=False)

    return {
        "total_rows": len(df),
        "unique_guide_bc0": len(merged),
        "mean_purity": merged["guide_bc0_purity"].mean(),
        "median_purity": merged["guide_bc0_purity"].median(),
        "high_purity_count": (merged["guide_bc0_purity"] >= 0.9).sum(),
        "low_purity_count": (merged["guide_bc0_purity"] < 0.5).sum(),
    }


def find_input_files(input_dir: Path, suffix: str, prefix: str = "") -> List[Path]:
    """
    Find all input files matching the pattern.

    Args:
        input_dir: Directory to search
        suffix: File suffix to match
        prefix: Optional library prefix filter (e.g., "V")

    Returns:
        List of matching file paths
    """
    pattern = f"{prefix}*{suffix}" if prefix else f"*{suffix}"
    return sorted(input_dir.glob(pattern))


def main():
    """Main processing pipeline."""
    # Load configuration from environment
    env_config = get_config_from_env()
    counts_dir = env_config["counts_dir"]
    purity_dir = env_config["purity_dir"]

    start_time = datetime.now()
    print("=" * 70)
    print("Step 04: Calculate BC0 Purity")
    print(f"Started: {start_time}")
    print(f"Project: {env_config['project_dir']}")
    print("=" * 70)

    # Setup
    config = PurityConfig()
    setup_directories(env_config)

    # Check input directory
    if not counts_dir.exists():
        print(f"ERROR: Input directory not found: {counts_dir}")
        sys.exit(1)

    # Find input files
    input_files = find_input_files(counts_dir, INPUT_SUFFIX, config.library_prefix)

    if not input_files:
        print(f"ERROR: No input files found matching pattern: {config.library_prefix}*{INPUT_SUFFIX}")
        print(f"  Directory: {counts_dir}")
        sys.exit(1)

    print(f"\nFound {len(input_files)} files to process")
    print(f"Output directory: {purity_dir}")

    # Process each file
    all_stats = []

    for input_file in tqdm(input_files, desc="Processing samples"):
        library_id = input_file.stem.replace(INPUT_SUFFIX.replace(".tsv", ""), "")
        output_file = purity_dir / f"{library_id}{OUTPUT_SUFFIX}"

        print(f"\n  Processing: {library_id}")
        stats = calculate_bc0_purity(input_file, output_file, config)
        stats["library_id"] = library_id
        all_stats.append(stats)

        if stats["mean_purity"] is not None:
            print(f"    Unique guide-bc0 pairs: {stats['unique_guide_bc0']:,}")
            print(f"    Mean purity: {stats['mean_purity']:.3f}")
            print(f"    High purity (>=0.9): {stats['high_purity_count']:,}")
            print(f"    Low purity (<0.5): {stats['low_purity_count']:,}")

    # Write summary statistics
    if all_stats:
        summary_df = pd.DataFrame(all_stats)
        summary_file = purity_dir / "purity_summary.tsv"
        summary_df.to_csv(summary_file, sep="\t", index=False)
        print(f"\nSummary written to: {summary_file}")

    # Final summary
    elapsed = datetime.now() - start_time
    print(f"\n{'=' * 70}")
    print("Step D Complete")
    print(f"  Processed {len(input_files)} files")
    print(f"  Output directory: {purity_dir}")
    print(f"  Elapsed time: {elapsed}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
