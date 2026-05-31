"""
Step 04 core logic: Calculate BC0 purity for guide-donor-BC0 combinations.

For each unique (guide, BC0) pair, purity is defined as the fraction of reads
that support the single most common donor sequence.  High-purity pairs represent
clean, single-donor BC0 assignments; low-purity pairs indicate contamination
from multiple donor sequences sharing the same BC0.

Author: Kevin R. Roy
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


INPUT_SUFFIX = "_matched_guide_donor_bc0_counts.tsv"
OUTPUT_SUFFIX = "_guide_donor_bc0_purity_filtered.tsv"


@dataclass
class PurityConfig:
    """Configuration for BC0 purity calculation."""
    min_counts: int = 1
    """Minimum read count for a row to be included in the analysis."""
    library_prefix: str = "V"
    """Library ID prefix used when searching for input files (e.g. ``"V"`` for V629)."""


def find_input_files(
    input_dir: Path,
    suffix: str = INPUT_SUFFIX,
    prefix: str = "",
) -> List[Path]:
    """
    Return sorted paths matching ``{prefix}*{suffix}`` in *input_dir*.

    Parameters
    ----------
    input_dir : Path
        Directory to search.
    suffix : str
        File suffix pattern.
    prefix : str
        Optional library-ID prefix filter.

    Returns
    -------
    list[Path]
        Sorted matching file paths.
    """
    pattern = f"{prefix}*{suffix}" if prefix else f"*{suffix}"
    return sorted(input_dir.glob(pattern))


def calculate_bc0_purity(
    input_file: Path,
    output_file: Path,
    config: Optional[PurityConfig] = None,
) -> Dict[str, object]:
    """
    Calculate purity for every guide-BC0 pair in one matched-counts file.

    For each (guide, bc0) group:

    1. Sum counts across all donor variants.
    2. Retain the row with the highest count (dominant donor).
    3. ``purity = dominant_count / total_count``

    Parameters
    ----------
    input_file : Path
        Matched guide-donor-BC0 counts file from Step 03.
        Must have columns: ``guide``, ``bc0``, ``counts``.
    output_file : Path
        Destination TSV file (created or overwritten).
    config : PurityConfig, optional
        Purity configuration.  Defaults to ``PurityConfig()``.

    Returns
    -------
    dict
        Statistics: total_rows, unique_guide_bc0, mean_purity, median_purity,
        high_purity_count (>= 0.9), low_purity_count (< 0.5).
    """
    if config is None:
        config = PurityConfig()

    df = pd.read_csv(input_file, sep="\t")

    if config.min_counts > 1:
        df = df[df["counts"] >= config.min_counts]

    if df.empty:
        return {"total_rows": 0, "unique_guide_bc0": 0, "mean_purity": None}

    # Total counts per guide-bc0 pair
    summed = (
        df.groupby(["guide", "bc0"])["counts"]
        .sum()
        .reset_index()
        .rename(columns={"counts": "grouped_guide_bc0_counts"})
    )

    # Dominant donor for each guide-bc0 pair
    idx = df.groupby(["guide", "bc0"])["counts"].idxmax()
    top = df.loc[idx].copy()

    merged = top.merge(summed, on=["guide", "bc0"])
    merged["guide_bc0_purity"] = merged["counts"] / merged["grouped_guide_bc0_counts"]

    merged.to_csv(output_file, sep="\t", index=False)

    return {
        "total_rows": len(df),
        "unique_guide_bc0": len(merged),
        "mean_purity": merged["guide_bc0_purity"].mean(),
        "median_purity": merged["guide_bc0_purity"].median(),
        "high_purity_count": (merged["guide_bc0_purity"] >= 0.9).sum(),
        "low_purity_count": (merged["guide_bc0_purity"] < 0.5).sum(),
    }
