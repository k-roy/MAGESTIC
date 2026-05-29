"""
magestic.pipelines.redi.core.cross_array — agreement between two REDI arrays.

Generalizes legacy `step_f_check_agreement_between_SHA435_and_yT182_REDI.py`,
which was hard-coded to "SHA345 REDI" vs "everything else". Engine version
takes a pair of selector predicates and returns a merged comparison table
with a `bc1_agreement` boolean and per-array purity/count columns prefixed by
array name.

Output columns include `<a>_bc1`, `<b>_bc1`, `<a>_colony_purity`,
`<b>_colony_purity`, `<a>_REDI_bc_bc1_counts`, `<b>_REDI_bc_bc1_counts`,
`<a>_REDI_bc_bc1_total_counts`, `<b>_REDI_bc_bc1_total_counts`,
`<a>_bc1_and_<b>_bc1_both_present`, and `bc1_agreement`.
"""

from __future__ import annotations

import pandas as pd


PLATE_GROUP_COLS = [
    "bc1",
    "plate_and_position",
    "MAT_alpha_array_1536_barcode",
    "REDI_plate_color_384",
    "REDI_row_384",
    "REDI_col_384",
]

PLATE_AGG = {
    "REDI_bc_bc1_counts": "sum",
    "REDI_bc_bc1_total_counts": "sum",
}


def _attach_plate_and_position(df: pd.DataFrame) -> pd.DataFrame:
    """Build the `plate_and_position` composite key used for cross-array joins."""
    df = df[
        df["MAT_alpha_array_1536_barcode"].notnull()
        & df["REDI_plate_color_384"].notnull()
    ].copy()
    df["plate_and_position"] = (
        df["MAT_alpha_array_1536_barcode"]
        .astype(str)
        .str.replace(r"\.0$", "", regex=True)
        + "_"
        + df["REDI_plate_color_384"]
        + "_"
        + df["REDI_row_384"].astype(str).str.replace(r"\.0$", "", regex=True)
        + "_"
        + df["REDI_col_384"].astype(str).str.replace(r"\.0$", "", regex=True)
    )
    return df


def _collapse_array(df: pd.DataFrame, array_name: str) -> pd.DataFrame:
    """Group, aggregate, dedupe-by-position, and prefix columns for one array."""
    g = df.groupby(PLATE_GROUP_COLS, as_index=False).agg(PLATE_AGG)
    g = g.sort_values(by="REDI_bc_bc1_counts", ascending=False).drop_duplicates(
        subset=["plate_and_position"]
    )
    g = g[g["REDI_bc_bc1_counts"] > 0].copy()
    g["colony_purity"] = (
        g["REDI_bc_bc1_counts"] / g["REDI_bc_bc1_total_counts"]
    )
    return g.rename(
        columns={
            "bc1": f"{array_name}_bc1",
            "REDI_bc_bc1_counts": f"{array_name}_REDI_bc_bc1_counts",
            "REDI_bc_bc1_total_counts": f"{array_name}_REDI_bc_bc1_total_counts",
            "colony_purity": f"{array_name}_colony_purity",
        }
    )


def compare_arrays(
    annotated_top_bc1_df: pd.DataFrame,
    array_a_mask: pd.Series,
    array_b_mask: pd.Series,
    array_a_name: str,
    array_b_name: str,
    min_counts: int = 10,
) -> pd.DataFrame:
    """Compare two arrays' REDI clone calls on shared plate-and-position keys.

    Parameters
    ----------
    annotated_top_bc1_df : DataFrame with the columns from
        `annotated_top_REDI_bc_MAGESTIC_bc1_df` (output of step_e).
    array_a_mask, array_b_mask : boolean Series selecting array A / B rows.
        Computed by the caller from `purpose` / validated flags / etc., so the
        engine doesn't hard-code SHA345 semantics.
    array_a_name, array_b_name : column-prefix names (e.g., "SHA345", "yT182").
    min_counts : drop rows with `REDI_bc_bc1_counts` below this floor
        (legacy step_f: 10).

    Returns
    -------
    Outer-merged dataframe with one row per `plate_and_position`, plus
    `bc1_agreement` and `<a>_bc1_and_<b>_bc1_both_present` flags.
    """
    df = annotated_top_bc1_df.query(f"REDI_bc_bc1_counts >= {min_counts}").copy()
    df = _attach_plate_and_position(df)
    df_a = _collapse_array(df[array_a_mask.loc[df.index]], array_a_name)
    df_b = _collapse_array(df[array_b_mask.loc[df.index]], array_b_name)
    merged = df_a.merge(df_b, how="outer", suffixes=(f"_{array_a_name}", f"_{array_b_name}"))
    merged[f"{array_a_name}_bc1_and_{array_b_name}_bc1_both_present"] = (
        merged[f"{array_a_name}_bc1"].notna()
        & merged[f"{array_b_name}_bc1"].notna()
    )
    merged["bc1_agreement"] = (
        merged[f"{array_a_name}_bc1"] == merged[f"{array_b_name}_bc1"]
    )
    return merged
