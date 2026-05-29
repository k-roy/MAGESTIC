"""
magestic.pipelines.redi.core.annotate — join REDI top-clone tables with the
bc1-donor-bc0 reference table.

Ported from `step_e_annotate_REDI_bc_MAGESTIC_bc_with_bc1_reference_table.py`.

Key behavior: overlapping column names (other than `bc1` and `counts`) are
prefixed `REDI_bc_bc1_` on the REDI side so the merged table is unambiguous.
`counts` is renamed to `REDI_bc_bc1_counts`. `amino_acid_position` is coerced
to nullable Int (legacy: fill NaN with -1 then int).
"""

from __future__ import annotations

import pandas as pd


def annotate_with_bc1_reference(
    top_clone_df: pd.DataFrame,
    bc1_reference_df: pd.DataFrame,
    redi_prefix: str = "REDI_bc_bc1_",
    aa_position_col: str = "amino_acid_position",
) -> pd.DataFrame:
    """Left-merge `top_clone_df` onto `bc1_reference_df` on `bc1`.

    Columns appearing in both (besides `bc1` and `counts`) are renamed on the
    REDI side with `redi_prefix`. `counts` → `redi_prefix + "counts"`.
    """
    overlap = set(top_clone_df.columns).intersection(set(bc1_reference_df.columns))
    rename_map = {
        col: f"{redi_prefix}{col}"
        for col in overlap
        if col not in {"bc1", "counts"}
    }
    if "counts" in top_clone_df.columns:
        rename_map["counts"] = f"{redi_prefix}counts"
    renamed = top_clone_df.rename(columns=rename_map)

    merged = renamed.merge(bc1_reference_df, on="bc1", how="left")

    if aa_position_col in merged.columns:
        merged[aa_position_col] = (
            pd.to_numeric(merged[aa_position_col], errors="coerce")
            .fillna(-1)
            .astype(int)
        )
    return merged
