"""
magestic.pipelines.redi.readers.keyfiles — keyfile loaders for REDI screens.

Ports `step_b_combine_keyfiles.py` and the REDI_bc key consolidation logic
from `step_c`.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def combine_keyfiles(
    sample_key_path: str | Path,
    gdna_plate_key_path: str | Path,
    redi_plate_barcode_key_path: str | Path,
) -> pd.DataFrame:
    """Merge sample_key + gDNA_plate_key + REDI_plate_barcode_key into a single
    `combined_key.tsv` analog, exactly as legacy step_b does.

    Critical: `REDI_plate_barcode` is coerced to numeric on both sides before
    the merge — both can otherwise contain non-integer placeholder values that
    break the join. (Legacy step_b does this with `pd.to_numeric(errors='coerce')`.)
    """
    sample_key = pd.read_csv(sample_key_path, sep="\t")
    gdna_plate_key = pd.read_csv(gdna_plate_key_path, sep="\t")
    redi_plate_barcode_key = pd.read_csv(redi_plate_barcode_key_path, sep="\t")

    combined = sample_key.merge(gdna_plate_key, on="gDNA_plate_name", how="left")
    combined["REDI_plate_barcode"] = pd.to_numeric(
        combined["REDI_plate_barcode"], errors="coerce"
    )
    redi_plate_barcode_key["REDI_plate_barcode"] = pd.to_numeric(
        redi_plate_barcode_key["REDI_plate_barcode"], errors="coerce"
    )
    combined = combined.merge(
        redi_plate_barcode_key, on="REDI_plate_barcode", how="left"
    )
    return combined


def load_redi_bc_keys(
    v5_key_path: str | Path | None = None,
    yt182_yt183_key_path: str | Path | None = None,
) -> pd.DataFrame:
    """Load + consolidate the REDI_bc keyfile(s).

    Ports the V5 + yT182/yT183 keyfile consolidation from step_c.
    The output schema is a uniform 6-column table:

        REDI_library, REDI_bc, REDI_plate_color_384, REDI_row_384,
        REDI_col_384, validated_by_yT100_mating (optional),
        confirmed_contaminated_by_yT100_mating (optional),
        clone_removed_in_RBYG_signature_set (optional)
    """
    parts: list[pd.DataFrame] = []

    if v5_key_path is not None:
        v5 = pd.read_csv(v5_key_path, sep="\t")
        v5["REDI_library"] = "SHA345"
        v5["REDI_bc"] = v5["REDI_bc_seq"].str.upper()
        v5_minimal = [
            "REDI_library",
            "REDI_bc",
            "REDI_plate_color_384",
            "REDI_row_384",
            "REDI_col_384",
        ]
        for opt in ["validated_by_yT100_mating", "confirmed_contaminated_by_yT100_mating"]:
            if opt in v5.columns:
                v5_minimal.append(opt)
        parts.append(v5[v5_minimal])

    if yt182_yt183_key_path is not None:
        y = pd.read_csv(yt182_yt183_key_path, sep="\t")
        y["REDI_bc"] = y["REDI_bc1_seq"].str.upper()
        y_minimal = ["REDI_library", "REDI_bc", "plate_384", "row_384", "col_384"]
        if "clone_removed_in_RBYG_signature_set" in y.columns:
            y_minimal.append("clone_removed_in_RBYG_signature_set")
        y = y[y_minimal]
        y = y.rename(
            columns={
                "plate_384": "REDI_plate_color_384",
                "row_384": "REDI_row_384",
                "col_384": "REDI_col_384",
            }
        )
        parts.append(y)

    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
