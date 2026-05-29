"""
magestic.pipelines.redi.qc.pipeline_qc — per-stage QC summaries.

Reads the standard step_c/d/f outputs and produces small summary tables for
the per-screen README and the pipeline_qc report.
"""

from __future__ import annotations

import pandas as pd


def summarize_read_attribution(combined_df: pd.DataFrame) -> pd.DataFrame:
    """For each sample_name, count rows attributed to bc1+REDI_bc vs bc1-only."""
    cols = ["sample_name"]
    if "purpose" in combined_df.columns:
        cols.append("purpose")
    return combined_df.groupby(cols, dropna=False).size().reset_index(name="n_rows")


def summarize_purity_distribution(top_clone_df: pd.DataFrame) -> pd.DataFrame:
    """Bin colony_purity into deciles and report row counts + median count."""
    if "colony_purity" not in top_clone_df.columns:
        raise KeyError("top_clone_df must contain 'colony_purity'")
    df = top_clone_df.copy()
    df["purity_bin"] = pd.cut(
        df["colony_purity"],
        bins=[0, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0001],
        labels=["<0.5", "0.5-0.7", "0.7-0.8", "0.8-0.9", "0.9-0.95", "0.95-0.99", ">=0.99"],
        include_lowest=True,
    )
    return (
        df.groupby("purity_bin", observed=True)
        .agg(n=("colony_purity", "size"), median_counts=("counts", "median"))
        .reset_index()
    )


def summarize_cross_array_agreement(
    merged_df: pd.DataFrame,
    array_a_name: str = "SHA345",
    array_b_name: str = "yT182",
) -> pd.DataFrame:
    """Counts of {both_present, A_only, B_only, agree, disagree}."""
    both_col = f"{array_a_name}_bc1_and_{array_b_name}_bc1_both_present"
    if both_col not in merged_df.columns:
        raise KeyError(f"missing column: {both_col}")
    both = merged_df[merged_df[both_col]]
    rows = {
        "both_present": int(both.shape[0]),
        f"{array_a_name}_only": int(
            (
                merged_df[f"{array_a_name}_bc1"].notna()
                & merged_df[f"{array_b_name}_bc1"].isna()
            ).sum()
        ),
        f"{array_b_name}_only": int(
            (
                merged_df[f"{array_b_name}_bc1"].notna()
                & merged_df[f"{array_a_name}_bc1"].isna()
            ).sum()
        ),
        "agree": int(both["bc1_agreement"].sum()),
        "disagree": int((~both["bc1_agreement"]).sum()),
    }
    return pd.DataFrame([rows]).T.reset_index()
