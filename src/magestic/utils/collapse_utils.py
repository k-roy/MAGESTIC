#!/usr/bin/env python3
"""
Optimized utilities for collapsing sequences by edit distance.

This module provides high-performance sequence collapsing functions using:
- rapidfuzz for faster edit distance calculations (3-10x faster than python-Levenshtein)
- Multiprocessing for parallel processing across groups
- Efficient memory management with iterators

Author: Kevin R. Roy

Usage:
    from utils.collapse_utils import collapse_sequences, collapse_sequences_by_edit_distance

    collapsed_df = collapse_sequences(
        df,
        group_cols=["sample_id", "bc1"],
        seq_col="donor",
        count_cols="counts",
        edit_distance_threshold=5,
        n_jobs=16
    )
"""

import multiprocessing
import pandas as pd
import numpy as np
from typing import List, Union, Tuple, Dict, Any, Generator

# Try rapidfuzz first (faster), fall back to Levenshtein
try:
    from rapidfuzz.distance import Levenshtein as rf_levenshtein
    def edit_distance(s1: str, s2: str) -> int:
        """Calculate Levenshtein edit distance using rapidfuzz."""
        return rf_levenshtein.distance(s1, s2)
    DISTANCE_BACKEND = "rapidfuzz"
except ImportError:
    try:
        from Levenshtein import distance as lev_distance
        def edit_distance(s1: str, s2: str) -> int:
            """Calculate Levenshtein edit distance using python-Levenshtein."""
            return lev_distance(s1, s2)
        DISTANCE_BACKEND = "python-Levenshtein"
    except ImportError:
        raise ImportError(
            "Neither rapidfuzz nor python-Levenshtein is available. "
            "Install one with: pip install rapidfuzz"
        )


def _process_group_worker(args: Tuple) -> List[Dict[str, Any]]:
    """
    Worker function for parallel group processing.

    Collapses similar sequences within a group by edit distance threshold.
    Uses a greedy algorithm that iteratively selects the highest-count sequence
    and merges all sequences within the edit distance threshold.

    Args:
        args: Tuple of (group_df, group_cols, seq_col, count_cols, edit_distance_threshold)

    Returns:
        List of dicts representing collapsed rows
    """
    group, group_cols, seq_col, count_cols, edit_distance_threshold = args
    group = group.copy()
    result_rows = []

    # Pre-compute whether count_cols is a list
    is_multi_count = isinstance(count_cols, list)

    while not group.empty:
        # Find sequence with highest count
        if is_multi_count:
            top_idx = group[count_cols].sum(axis=1).idxmax()
        else:
            top_idx = group[count_cols].idxmax()

        top_seq = group.loc[top_idx, seq_col]

        # Skip if top sequence is NA
        if pd.isna(top_seq):
            group = group.drop(top_idx)
            continue

        # Calculate edit distance to top sequence
        # Using a list comprehension is faster than .apply() for small groups
        sequences = group[seq_col].values
        distances = np.array([
            edit_distance(str(seq), str(top_seq)) if pd.notna(seq) else 999
            for seq in sequences
        ])

        collapse_mask = distances <= edit_distance_threshold

        # Aggregate counts for similar sequences
        if is_multi_count:
            collapsed_counts = {
                col: group.loc[collapse_mask, col].sum()
                for col in count_cols
            }
        else:
            collapsed_counts = {count_cols: group.loc[collapse_mask, count_cols].sum()}

        # Build result row
        result_row = {col: group.iloc[0][col] for col in group_cols}
        result_row[seq_col] = top_seq
        result_row.update(collapsed_counts)
        result_rows.append(result_row)

        # Remove collapsed sequences
        group = group.loc[~collapse_mask]

    return result_rows


def _generate_group_args(
    grouped,
    group_cols: List[str],
    seq_col: str,
    count_cols: Union[str, List[str]],
    edit_distance_threshold: int
) -> Generator[Tuple, None, None]:
    """Generate arguments for parallel processing."""
    for _, group in grouped:
        yield (group, group_cols, seq_col, count_cols, edit_distance_threshold)


def collapse_sequences(
    df: pd.DataFrame,
    group_cols: List[str],
    seq_col: str,
    count_cols: Union[str, List[str]],
    edit_distance_threshold: int,
    show_progress: bool = True,
    n_jobs: int = 1
) -> pd.DataFrame:
    """
    Collapse similar sequences within each group based on edit distance threshold.

    This function groups the input DataFrame by `group_cols`, then within each group,
    iteratively merges sequences that are within `edit_distance_threshold` of each other.
    The sequence with the highest count is kept as the representative, and counts from
    similar sequences are summed.

    Args:
        df: Input DataFrame containing sequences and counts
        group_cols: Columns to group by before collapsing
        seq_col: Column containing sequences to collapse
        count_cols: Column(s) containing counts to aggregate (str or list of str)
        edit_distance_threshold: Maximum edit distance for collapsing sequences
        show_progress: Show progress bar (requires tqdm)
        n_jobs: Number of parallel jobs (1 for serial, >1 for parallel)

    Returns:
        DataFrame with collapsed sequences and aggregated counts

    Example:
        >>> collapsed = collapse_sequences(
        ...     df,
        ...     group_cols=["sample_id", "bc1"],
        ...     seq_col="donor",
        ...     count_cols="counts",
        ...     edit_distance_threshold=5,
        ...     n_jobs=16
        ... )
    """
    grouped = df.groupby(group_cols)
    n_groups = grouped.ngroups

    if show_progress:
        try:
            from tqdm import tqdm
            use_tqdm = True
        except ImportError:
            use_tqdm = False
            print(f"Collapsing {seq_col} ({n_groups} groups)...")
    else:
        use_tqdm = False

    if n_jobs > 1:
        # Parallel processing
        with multiprocessing.Pool(n_jobs) as pool:
            chunksize = max(1, n_groups // (n_jobs * 4))
            args_gen = _generate_group_args(
                grouped, group_cols, seq_col, count_cols, edit_distance_threshold
            )
            iterator = pool.imap_unordered(_process_group_worker, args_gen, chunksize=chunksize)

            if use_tqdm and show_progress:
                iterator = tqdm(iterator, total=n_groups, desc=f"Collapsing {seq_col}")

            result_rows = [row for group_result in iterator for row in group_result]
    else:
        # Serial processing
        groups = [group for _, group in grouped]
        if use_tqdm and show_progress:
            groups = tqdm(groups, desc=f"Collapsing {seq_col}")

        result_rows = []
        for group in groups:
            args = (group, group_cols, seq_col, count_cols, edit_distance_threshold)
            result_rows.extend(_process_group_worker(args))

    return pd.DataFrame(result_rows)


def collapse_sequences_by_edit_distance(
    df: pd.DataFrame,
    seq_col: str,
    count_col: str,
    edit_distance_threshold: int,
    additional_cols: List[str] = None
) -> pd.DataFrame:
    """
    Collapse sequences globally (no grouping) by edit distance.

    This is a simpler interface for collapsing all sequences in a DataFrame
    without prior grouping. Useful for final aggregation steps.

    Args:
        df: Input DataFrame
        seq_col: Column containing sequences
        count_col: Column containing counts
        edit_distance_threshold: Maximum edit distance for collapsing
        additional_cols: Additional columns to preserve (takes value from top sequence)

    Returns:
        DataFrame with collapsed sequences
    """
    if additional_cols is None:
        additional_cols = []

    df = df.sort_values(count_col, ascending=False).reset_index(drop=True)
    result_rows = []
    used = set()

    for idx, row in df.iterrows():
        if idx in used:
            continue

        seq = row[seq_col]
        if pd.isna(seq):
            continue

        # Find all sequences within threshold
        total_count = row[count_col]
        used.add(idx)

        for other_idx, other_row in df.iloc[idx+1:].iterrows():
            if other_idx in used:
                continue
            other_seq = other_row[seq_col]
            if pd.isna(other_seq):
                continue
            if edit_distance(str(seq), str(other_seq)) <= edit_distance_threshold:
                total_count += other_row[count_col]
                used.add(other_idx)

        result = {seq_col: seq, count_col: total_count}
        for col in additional_cols:
            result[col] = row[col]
        result_rows.append(result)

    return pd.DataFrame(result_rows)


def get_distance_backend() -> str:
    """Return the name of the edit distance backend being used."""
    return DISTANCE_BACKEND


# Module-level info
__all__ = [
    "collapse_sequences",
    "collapse_sequences_by_edit_distance",
    "edit_distance",
    "get_distance_backend"
]

if __name__ == "__main__":
    # Quick test
    print(f"Edit distance backend: {get_distance_backend()}")
    print(f"Test distance('hello', 'hallo'): {edit_distance('hello', 'hallo')}")
