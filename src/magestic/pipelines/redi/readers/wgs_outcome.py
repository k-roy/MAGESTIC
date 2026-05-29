"""
magestic.pipelines.redi.readers.wgs_outcome — wrapper around Shengdi's WGS
outcome table.

Canonical input:
    /oak/.../collaborator/projects/WGS_MAGESTIC_QTL_bc_calling/output3/summaries/
    wgs_outcome_full_table.v0.tsv

A newer `v20250701.tsv` exists alongside; configs choose which version.
Schema documented in MAGESTIC-SCORE repo.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def load_wgs_outcome(
    path: str | Path,
    bc1_col: str = "bc1_from_wgs",
    keep_cols: list[str] | None = None,
) -> pd.DataFrame:
    """Read Shengdi's WGS outcome table.

    Parameters
    ----------
    path : path to `wgs_outcome_full_table.v?.tsv`.
    bc1_col : column name to rename → `bc1` (legacy is `bc1_from_wgs`).
    keep_cols : optional whitelist; if provided, only these columns are kept
        (besides the renamed `bc1`).
    """
    df = pd.read_csv(path, sep="\t")
    if bc1_col in df.columns and bc1_col != "bc1":
        df = df.rename(columns={bc1_col: "bc1"})
    if keep_cols is not None:
        cols = [c for c in keep_cols if c in df.columns]
        if "bc1" not in cols:
            cols.append("bc1")
        df = df[cols].drop_duplicates().reset_index(drop=True)
    return df


def load_wgs_plate_key(
    wgs_plate_key_path: str | Path,
    plate_384_to_96_format_path: str | Path,
) -> pd.DataFrame:
    """Load WGS_plate_key.tsv and attach 96-well sample_number via 384→96 format.

    Ports step_g's `WGS_plate_key_df` construction. The result has columns:

        plate_384, row_384, col_384, row_96, col_96, sample_number,
        sequencing_date, "Nextera outer primer set", sample_name
    """
    plate_384_to_96 = pd.read_csv(plate_384_to_96_format_path, sep="\t")
    plate_96 = plate_384_to_96[["row_96", "col_96"]].drop_duplicates().reset_index(drop=True)
    plate_96["sample_number"] = "sample_" + (plate_96.index + 1).astype(str)
    plate_384_to_96 = plate_384_to_96.merge(
        plate_96[["row_96", "col_96", "sample_number"]],
        on=["row_96", "col_96"],
        how="left",
    )
    wgs_plate_key = pd.read_csv(wgs_plate_key_path, sep="\t")
    wgs_plate_key = wgs_plate_key.merge(plate_384_to_96)
    wgs_plate_key["sample_name"] = (
        wgs_plate_key["sequencing_date"].astype(str)
        + "_"
        + wgs_plate_key["Nextera outer primer set"]
        + "_"
        + wgs_plate_key["sample_number"]
    )
    return wgs_plate_key
