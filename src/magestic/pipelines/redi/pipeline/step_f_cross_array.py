"""
Step F — cross-array agreement (legacy `step_f_check_agreement_between_SHA435_and_yT182_REDI.py`).

Engine version reads array-pair definitions from `config.array_pairs`. For
each pair, calls `core.cross_array.compare_arrays` and writes per-pair output.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..core.config import REDIConfig
from ..core.cross_array import compare_arrays


def run_step_f(
    config: REDIConfig,
    annotated_path: Path | None = None,
) -> list[Path]:
    """Run step_f for every configured array pair. Returns list of output paths."""
    combined_data = config.paths.combined_data
    if annotated_path is None:
        annotated_path = combined_data / "annotated_top_REDI_bc_MAGESTIC_bc1_df.tsv"

    df = pd.read_csv(annotated_path, sep="\t")

    # Attach 1536 → 384 format if present in keyfiles
    plate_1536_path = config.paths.keyfiles / "1536_PIXL_rearray_format.tsv"
    plate_1536_df = (
        pd.read_csv(plate_1536_path, sep="\t") if plate_1536_path.exists() else None
    )

    out_paths: list[Path] = []
    for pair in config.array_pairs:
        # array A: purpose == purpose_a, validated_by_yT100_mating filter
        # applied when those columns exist
        mask_a = df["purpose"] == pair.purpose_a
        if "validated_by_yT100_mating" in df.columns:
            mask_a = mask_a & (df["validated_by_yT100_mating"] == True)  # noqa: E712
        if "confirmed_contaminated_by_yT100_mating" in df.columns:
            mask_a = mask_a & (
                df["confirmed_contaminated_by_yT100_mating"] == False  # noqa: E712
            )
        # array B: purpose != purpose_a (legacy SHA345-vs-everything-else
        # behavior) unless a more specific `purpose_b_not` is supplied
        mask_b = df["purpose"] != pair.purpose_a

        merged = compare_arrays(
            df,
            array_a_mask=mask_a,
            array_b_mask=mask_b,
            array_a_name=pair.name_a,
            array_b_name=pair.name_b,
            min_counts=config.purity.min_counts,
        )

        if plate_1536_df is not None:
            merged["plate_384"] = merged["REDI_plate_color_384"]
            merged["row_384"] = merged["REDI_row_384"]
            merged["col_384"] = merged["REDI_col_384"]
            merged = merged.merge(plate_1536_df, how="left")

        out = combined_data / f"{pair.name_a}_{pair.name_b}_merged_REDI_df.tsv"
        merged.to_csv(out, sep="\t", index=False)
        out_paths.append(out)
    return out_paths
