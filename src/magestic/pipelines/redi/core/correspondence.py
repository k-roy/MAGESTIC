"""
magestic.pipelines.redi.core.correspondence — unified two-type bc1 correspondence API.

REDI maintains two complementary correspondences between barcodes (MASTER_PLAN
§9d, J.3):

  Type (i)  bc1 ↔ edit   from WGS         (registered sequencing observations)
  Type (ii) REDI-bc ↔ MAGESTIC-bc1 from re-REDI (high-throughput, same array)

This module exposes a single API for both and emits a joined publishable
artifact `data/REDI_bc1_correspondence_master.tsv` with a `concordance` flag.

Caller responsibility: type-(i) and type-(ii) tables are produced by upstream
steps (step_g + step_e respectively) and passed in here as DataFrames.
Agreement is self-consistency unless an independent expected identity exists.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd
from collections.abc import Sequence
from .identity import keyed_frame, valid_barcode


@dataclass
class CorrespondenceResult:
    """Joined type-(i) + type-(ii) table."""

    master: pd.DataFrame
    n_type_i_only: int
    n_type_ii_only: int
    n_both_agree: int
    n_both_disagree: int
    n_both_unassessed: int = 0

    def summary(self) -> dict[str, int]:
        return {
            "WGS_only": self.n_type_i_only,
            "reREDI_only": self.n_type_ii_only,
            "both_agree": self.n_both_agree,
            "both_disagree": self.n_both_disagree,
            "both_unassessed": self.n_both_unassessed,
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
    from ..readers.wgs_outcome import load_wgs_outcome
    df = load_wgs_outcome(wgs_outcome_path, bc1_col=bc1_col)
    rename = {}
    if edit_col and edit_col in df.columns:
        if edit_col != 'edit_from_wgs' and 'edit_from_wgs' in df:
            raise ValueError('edit observation rename collision with existing edit_from_wgs')
        rename[edit_col] = "edit_from_wgs"
    return df.rename(columns=rename)


def build_type_ii_redi_to_magestic(
    cross_array_df: pd.DataFrame,
    redi_bc_col: str = "REDI_bc",
    magestic_bc1_col: str = "REDI_called_bc1",
    *, position_keys: Sequence[str] = ("plate_and_position",),
) -> pd.DataFrame:
    """Project cross-array output into the type-(ii) `REDI_bc ↔ MAGESTIC_bc1` schema.

    Caller passes the cross-array dataframe (output of `compare_arrays` or
    `compare_to_wgs`); this returns a deduplicated mapping table.
    """
    cross_array_df = keyed_frame(cross_array_df, position_keys, label="type-II correspondence")
    cols = list(position_keys) + [redi_bc_col, magestic_bc1_col]
    keep = [c for c in cols if c in cross_array_df.columns]
    if len(keep) != len(cols):
        raise KeyError(
            f"Cross-array dataframe missing required columns: "
            f"need {cols}, found {keep}"
        )
    out = cross_array_df[cols].drop_duplicates()
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
    *, position_keys: Sequence[str] = ("plate_and_position",),
    bc1_length: int = 20,
) -> CorrespondenceResult:
    """Outer-join independently registered physical positions, then compare calls.

    Concordance flag values:
        "WGS_only"     row present in type_i only
        "reREDI_only"  row present in type_ii only
        "both_agree"   both present and edit_from_wgs / MAGESTIC linkage consistent
        "both_disagree" both present but disagree
    """
    if bc1_col_i in position_keys or bc1_col_ii in position_keys:
        raise ValueError("position keys cannot be observed barcodes")
    reserved = {'_in_wgs','_in_redi','_valid_wgs','_valid_redi','concordance'}
    collisions = reserved & (set(type_i) | set(type_ii))
    if collisions:
        raise ValueError(f'reserved correspondence column collision: {sorted(collisions)}')
    left = keyed_frame(type_i, position_keys, label="WGS correspondence")
    right = keyed_frame(type_ii, position_keys, label="REDI correspondence")
    if bc1_col_i == bc1_col_ii:
        raise ValueError("distinct observation column names are required")
    if bc1_col_i not in left or bc1_col_ii not in right:
        raise ValueError("barcode observation columns are required")
    if bc1_col_i in right or bc1_col_ii in left:
        raise ValueError("observation columns collide across correspondence tables")
    left["_in_wgs"] = True
    right["_in_redi"] = True
    merged = left.merge(
        right, how="outer", on=list(position_keys), validate="one_to_one"
    )
    merged["_in_wgs"] = merged["_in_wgs"].fillna(False)
    merged["_in_redi"] = merged["_in_redi"].fillna(False)

    def _flag(row: pd.Series) -> str:
        if row["_in_wgs"] and not row["_in_redi"]:
            return "WGS_only"
        if row["_in_redi"] and not row["_in_wgs"]:
            return "reREDI_only"
        if not row["_valid_wgs"] or not row["_valid_redi"]:
            return "both_unassessed"
        if row[bc1_col_i] == row[bc1_col_ii]:
            return "both_agree"
        return "both_disagree"

    merged["_valid_wgs"] = valid_barcode(merged[bc1_col_i], bc1_length)
    merged["_valid_redi"] = valid_barcode(merged[bc1_col_ii], bc1_length)
    merged["concordance"] = pd.Series([_flag(row) for _, row in merged.iterrows()],index=merged.index,dtype="string")
    merged = merged.drop(columns=["_in_wgs", "_in_redi", "_valid_wgs", "_valid_redi"])

    counts = merged["concordance"].value_counts().to_dict()
    return CorrespondenceResult(
        master=merged,
        n_type_i_only=int(counts.get("WGS_only", 0)),
        n_type_ii_only=int(counts.get("reREDI_only", 0)),
        n_both_agree=int(counts.get("both_agree", 0)),
        n_both_disagree=int(counts.get("both_disagree", 0)),
        n_both_unassessed=int(counts.get("both_unassessed", 0)),
    )


def write_correspondence_master(
    result: CorrespondenceResult,
    out_path: str,
) -> None:
    """Write the joined master TSV (`data/REDI_bc1_correspondence_master.tsv`)."""
    result.master.to_csv(out_path, sep="\t", index=False)
