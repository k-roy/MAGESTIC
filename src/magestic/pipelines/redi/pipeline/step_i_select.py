"""
Step I — generic clone-selection helpers.

The legacy `step_i_select_clones_for_NNS_collection.py` is per-screen. The
engine exposes the reusable selection primitives; project-specific selection
logic (NRD1/NAB3/SEN1 ORF maps, sublib_to_gene dicts, 3-clones-per-site cap)
lives in `REDI/screens/<screen>/driver.py`.

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
    if len(sort_by) != len(ascending):
        raise ValueError("len(sort_by) must equal len(ascending)")
    return (
        df.sort_values(by=list(sort_by), ascending=list(ascending))
        .groupby(group_cols)
        .head(n)
        .reset_index(drop=True)
    )
