"""
magestic.pipelines.redi.core.purity — colony purity scoring.

Ported from `20250312_NNS_complex_satmut/scripts/step_d_*`. For each
(sample_name, REDI_bc) group:

1.  Identify the top bc1 (highest counts) → `top_bc1`.
2.  Compute Levenshtein distance from every observed bc1 to `top_bc1`; drop
    rows where the distance is in `edit_distance_drop` (default {1, 2}) —
    these are sequencing-error daughters of the dominant bc1.
3.  Sum the top-N (default 4) bc1 read counts → `total_counts`.
4.  Colony purity = top_bc1.counts / total_counts.

Returned columns: `top_bc1`, `edit_distance`, `total_counts`, `colony_purity`.

The engine version is config-driven (no hard-coded screen path) and accepts
either an in-memory DataFrame or a path to a TSV.
"""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd

try:  # rapidfuzz is faster, python-Levenshtein is the legacy choice
    from Levenshtein import distance as lev_distance
except ImportError:  # pragma: no cover
    from rapidfuzz.distance.Levenshtein import distance as lev_distance


def top_bc1_per_redi_bc(df: pd.DataFrame, group_cols: Iterable[str]) -> pd.DataFrame:
    """For each group, return the (group_cols, top_bc1) where top_bc1 is the
    bc1 with the highest `counts`.

    Equivalent to legacy step_d:

        top_bc1_df = df.sort_values('counts', ascending=False) \
            .groupby([...], as_index=False).first()[[..., 'bc1']]
    """
    group_cols = list(group_cols)
    top = (
        df.sort_values("counts", ascending=False)
        .groupby(group_cols, as_index=False)
        .first()[group_cols + ["bc1"]]
        .rename(columns={"bc1": "top_bc1"})
    )
    return top


def _edit_distance_vec(top: np.ndarray, observed: np.ndarray) -> list[float]:
    """Pairwise Levenshtein distance with NaN passthrough.

    Equivalent of step_d's `fast_edit_distance`.
    """
    out: list[float] = []
    for t, o in zip(top, observed):
        if isinstance(t, float) and np.isnan(t):
            out.append(np.nan)
        elif isinstance(o, float) and np.isnan(o):
            out.append(np.nan)
        elif t is None or o is None:
            out.append(np.nan)
        else:
            out.append(lev_distance(t, o))
    return out


def filter_edit_distance(
    df: pd.DataFrame,
    drop_distances: Iterable[int] = (1, 2),
) -> pd.DataFrame:
    """Drop rows where `edit_distance` is in `drop_distances`.

    Requires that `top_bc1` and `edit_distance` are already populated.
    """
    drop_set = set(drop_distances)
    return df[~df["edit_distance"].isin(drop_set)].copy()


def compute_purity(
    df: pd.DataFrame,
    group_cols: Iterable[str] = ("sample_name", "REDI_bc"),
    top_n: int = 4,
    edit_distance_drop: Iterable[int] = (1, 2),
) -> pd.DataFrame:
    """Run the full step_d purity calculation.

    Parameters
    ----------
    df : DataFrame with columns at least
        [bc1, REDI_bc, sample_name, counts, REDI_plate_color_384]
        plus the metadata columns from `combined_REDI_bc_MAGESTIC_bc_merged_df`.
    group_cols : grouping columns for the per-clone aggregation.
    top_n : number of top bc1's summed for the purity denominator (legacy 4).
    edit_distance_drop : edit-distance values to drop as sequencing-error
        daughters of the dominant bc1 (legacy {1, 2}).

    Returns
    -------
    A `top_bc1_df` with one row per (group_cols), columns including
    `top_bc1`, `edit_distance`, `colony_purity`, `total_counts`, and all
    metadata columns from the input.
    """
    group_cols = list(group_cols)

    # Step 1: top bc1 per group
    top = top_bc1_per_redi_bc(df, group_cols)
    df = df.merge(top, on=group_cols, how="left")

    # Step 2: edit distance from top_bc1
    df["edit_distance"] = _edit_distance_vec(
        df["top_bc1"].values, df["bc1"].values
    )

    # Step 3: drop near-duplicate bc1's (sequencing errors)
    df = filter_edit_distance(df, edit_distance_drop)

    # Drop rows without plate annotation (legacy step_d explicit filter)
    if "REDI_plate_color_384" in df.columns:
        df = df[df["REDI_plate_color_384"].notna()].copy()

    # Step 4: top-N totals
    df = df.sort_values(by="counts", ascending=False)
    top_n_df = df.groupby(group_cols).head(top_n).reset_index(drop=True)

    total_counts = (
        top_n_df.groupby(group_cols)["counts"]
        .sum()
        .reset_index()
        .rename(columns={"counts": "total_counts"})
    )
    top_n_df = top_n_df.merge(total_counts, on=group_cols, how="left")
    top_n_df["colony_purity"] = top_n_df["counts"] / top_n_df["total_counts"]

    # Step 5: keep only the top bc1 per group (one row per clone)
    top_n_df = top_n_df.sort_values(
        by=group_cols + ["counts"],
        ascending=[True] * len(group_cols) + [False],
    )
    top_clone = top_n_df.groupby(group_cols).first().reset_index()
    return top_clone
