"""
magestic.pipelines.redi.core.correspondence — unified two-type bc1 correspondence API.

REDI maintains two complementary correspondences between barcodes (MASTER_PLAN
§9d, J.3):

  Type (i)  bc1 ↔ edit   from WGS         (low-throughput ground truth)
  Type (ii) REDI-bc ↔ MAGESTIC-bc1 from re-REDI (high-throughput, same array)

This module exposes a single API for both and emits a joined publishable
artifact `data/REDI_bc1_correspondence_master.tsv` with a `concordance` flag.

Caller responsibility: type-(i) and type-(ii) tables are produced by upstream
steps (step_g + step_e respectively) and passed in here as DataFrames.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass
class CorrespondenceResult:
    """Joined type-(i) + type-(ii) table."""

    master: pd.DataFrame
    n_type_i_only: int
    n_type_ii_only: int
    n_both_agree: int
    n_both_disagree: int

    def summary(self) -> dict[str, int]:
        return {
            "WGS_only": self.n_type_i_only,
            "reREDI_only": self.n_type_ii_only,
            "both_agree": self.n_both_agree,
            "both_disagree": self.n_both_disagree,
            "total": len(self.master),
        }


def load_type_i_wgs(
    wgs_outcome_path: str,
    bc1_col: str = "bc1_from_wgs",
    edit_col: str | None = "edit_called_from_wgs",
) -> pd.DataFrame:
    """Load Shengdi's WGS outcome table → `bc1_to_edit_from_WGS.tsv` schema.

    Output schema:
        bc1, edit_from_wgs, sample_name (and any other passthrough cols
        present in the source).
    """
    df = pd.read_csv(wgs_outcome_path, sep="\t")
    rename = {bc1_col: "bc1"}
    if edit_col and edit_col in df.columns:
        rename[edit_col] = "edit_from_wgs"
    return df.rename(columns=rename)


def build_type_ii_redi_to_magestic(
    cross_array_df: pd.DataFrame,
    redi_bc_col: str = "REDI_bc",
    magestic_bc1_col: str = "REDI_called_bc1",
) -> pd.DataFrame:
    """Project cross-array output into the type-(ii) `REDI_bc ↔ MAGESTIC_bc1` schema.

    Caller passes the cross-array dataframe (output of `compare_arrays` or
    `compare_to_wgs`); this returns a deduplicated mapping table.
    """
    cols = [redi_bc_col, magestic_bc1_col]
    keep = [c for c in cols if c in cross_array_df.columns]
    if len(keep) != len(cols):
        raise KeyError(
            f"Cross-array dataframe missing required columns: "
            f"need {cols}, found {keep}"
        )
    out = cross_array_df[cols].drop_duplicates()
    out = out.dropna(subset=[magestic_bc1_col])
    out = out.rename(
        columns={
            redi_bc_col: "REDI_bc",
            magestic_bc1_col: "MAGESTIC_bc1",
        }
    )
    return out


def join_correspondence(
    type_i: pd.DataFrame,
    type_ii: pd.DataFrame,
    bc1_col_i: str = "bc1",
    bc1_col_ii: str = "MAGESTIC_bc1",
) -> CorrespondenceResult:
    """Outer-join type-(i) + type-(ii) on bc1 and emit master + summary.

    Concordance flag values:
        "WGS_only"     row present in type_i only
        "reREDI_only"  row present in type_ii only
        "both_agree"   both present and edit_from_wgs / MAGESTIC linkage consistent
        "both_disagree" both present but disagree
    """
    left = type_i.copy()
    right = type_ii.copy()
    left["_in_wgs"] = True
    right["_in_redi"] = True
    merged = left.merge(
        right, how="outer", left_on=bc1_col_i, right_on=bc1_col_ii
    )
    merged["_in_wgs"] = merged["_in_wgs"].fillna(False)
    merged["_in_redi"] = merged["_in_redi"].fillna(False)

    def _flag(row: pd.Series) -> str:
        if row["_in_wgs"] and not row["_in_redi"]:
            return "WGS_only"
        if row["_in_redi"] and not row["_in_wgs"]:
            return "reREDI_only"
        # both present. If the bc1 keys themselves agree (left = right) it's
        # agreement; otherwise disagreement.
        if row[bc1_col_i] == row[bc1_col_ii]:
            return "both_agree"
        return "both_disagree"

    merged["concordance"] = merged.apply(_flag, axis=1)
    merged = merged.drop(columns=["_in_wgs", "_in_redi"])

    counts = merged["concordance"].value_counts().to_dict()
    return CorrespondenceResult(
        master=merged,
        n_type_i_only=int(counts.get("WGS_only", 0)),
        n_type_ii_only=int(counts.get("reREDI_only", 0)),
        n_both_agree=int(counts.get("both_agree", 0)),
        n_both_disagree=int(counts.get("both_disagree", 0)),
    )


def write_correspondence_master(
    result: CorrespondenceResult,
    out_path: str,
) -> None:
    """Write the joined master TSV (`data/REDI_bc1_correspondence_master.tsv`)."""
    result.master.to_csv(out_path, sep="\t", index=False)
