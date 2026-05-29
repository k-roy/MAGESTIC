"""
magestic.pipelines.redi.core.wgs_check — compare REDI-called bc1 to WGS bc1.

Ported from `step_g_check_REDI_bc1_against_WGS_bc1.py`. WGS is the
low-throughput ground truth (type-(i) correspondence per MASTER_PLAN §9d).

Behavior: WGS bc1 calls are truncated/padded to `bc1_length` (legacy: 20, with
'C' as the right-pad fill character) before comparison, matching the legacy
script.
"""

from __future__ import annotations

import pandas as pd


def normalize_wgs_bc1(series: pd.Series, length: int = 20, pad_char: str = "C") -> pd.Series:
    """Truncate/pad WGS bc1 to a fixed length, matching legacy step_g."""
    out = series.str.slice(0, length)
    return out.str.pad(length, side="right", fillchar=pad_char)


def compare_to_wgs(
    redi_df: pd.DataFrame,
    wgs_bc1_call_df: pd.DataFrame,
    wgs_plate_key_df: pd.DataFrame,
    array_a_name: str = "SHA345",
    array_b_name: str = "yT182",
    bc1_length: int = 20,
) -> pd.DataFrame:
    """Compare WGS bc1 calls against REDI-called bc1.

    Parameters
    ----------
    redi_df : merged cross-array dataframe (output of `compare_arrays`).
    wgs_bc1_call_df : DataFrame with columns `sample_name`, `bc1_from_wgs`.
    wgs_plate_key_df : DataFrame mapping sample_name to plate/well.
    array_a_name, array_b_name : column prefix names matching the
        `compare_arrays` output.
    bc1_length : length to normalize WGS bc1 to (default 20).

    Returns
    -------
    Outer merge with agreement booleans:
        `wgs_and_<a>_bc1_agreement`, `wgs_and_<b>_bc1_agreement`,
        `wgs_and_<a>_<b>_agreement`, `REDI_called_bc1`.
    """
    wgs = wgs_bc1_call_df.copy()
    wgs["bc1_from_wgs"] = normalize_wgs_bc1(wgs["bc1_from_wgs"], length=bc1_length)
    merged = wgs_plate_key_df.merge(wgs, on="sample_name", how="left").merge(
        redi_df, how="left"
    )
    merged[f"wgs_and_{array_a_name}_bc1_agreement"] = (
        merged["bc1_from_wgs"] == merged[f"{array_a_name}_bc1"]
    )
    merged[f"wgs_and_{array_b_name}_bc1_agreement"] = (
        merged["bc1_from_wgs"] == merged[f"{array_b_name}_bc1"]
    )
    merged[f"wgs_and_{array_a_name}_{array_b_name}_agreement"] = (
        merged[f"wgs_and_{array_a_name}_bc1_agreement"]
        & merged[f"wgs_and_{array_b_name}_bc1_agreement"]
    )
    merged["REDI_called_bc1"] = merged[f"{array_b_name}_bc1"].fillna(
        merged[f"{array_a_name}_bc1"]
    )
    return merged
