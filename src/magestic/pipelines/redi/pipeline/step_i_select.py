"""
Step I — generic clone-selection helpers.

The legacy `step_i_select_clones_for_NNS_collection.py` is per-screen. The
engine exposes reusable selection primitives. Collection-specific eligibility,
caps, physical transfer registration and robot worklist generation still need
a separate implementation; no screen driver is bundled here.

Reusable primitives:
    select_top_n_per_group  — generic "take N highest-purity clones per group"
"""

from __future__ import annotations

import pandas as pd


def select_top_n_per_group(
    df: pd.DataFrame,
    group_cols: list[str],
    n: int = 3,
    sort_by: tuple[str, ...] = ("mean_purity",),
    ascending: tuple[bool, ...] = (False,),
) -> pd.DataFrame:
    """Take the top N rows per group, sorted descending by `sort_by`."""
    if type(n) is not int or n < 1:
        raise ValueError("n must be a positive integer")
    if not group_cols or len(group_cols) != len(set(group_cols)):
        raise ValueError("group columns must be nonempty and unique")
    if df.columns.duplicated().any() or any(c not in df for c in [*group_cols,*sort_by]):
        raise ValueError("selection requires unique, present group and sort columns")
    if df[group_cols].isna().any().any() or any(df[c].astype('string').str.strip().eq('').any() for c in group_cols):
        raise ValueError("selection group identities must not be missing")
    if not sort_by or len(set(sort_by)) != len(sort_by):
        raise ValueError("sort columns must be nonempty and unique")
    if len(sort_by) != len(ascending):
        raise ValueError("len(sort_by) must equal len(ascending)")
    return (
        df.sort_values(by=list(sort_by), ascending=list(ascending), kind="stable")
        .groupby(group_cols)
        .head(n)
        .reset_index(drop=True)
    )
