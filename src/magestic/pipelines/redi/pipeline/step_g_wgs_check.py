"""
Step G — compare REDI-called bc1 against WGS bc1 (type-(i) ground truth).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..core.config import REDIConfig
from ..core.wgs_check import compare_to_wgs
from ..readers.wgs_outcome import load_wgs_plate_key


def run_step_g(
    config: REDIConfig,
    cross_array_path: Path | None = None,
    mat_alpha_array_barcode_filter: int | None = 4996,
) -> Path:
    """Run step_g. Returns path to `WGS_bc1_call_outer_merge_with_REDI_df_and_plate_key.tsv`.

    `mat_alpha_array_barcode_filter`: legacy step_g filters to
    `MAT_alpha_array_1536_barcode == 4996` (the WGS subset). Pass None to skip.
    """
    if config.wgs_outcome_table is None:
        raise ValueError("config.wgs_outcome_table is required for step_g")
    combined_data = config.paths.combined_data

    # Default to the first array pair's output file
    if cross_array_path is None:
        if not config.array_pairs:
            raise ValueError(
                "No array_pairs configured; step_g needs a cross-array file or "
                "an explicit cross_array_path."
            )
        pair = config.array_pairs[0]
        cross_array_path = (
            combined_data / f"{pair.name_a}_{pair.name_b}_merged_REDI_df.tsv"
        )

    cross = pd.read_csv(cross_array_path, sep="\t")
    if mat_alpha_array_barcode_filter is not None:
        cross = cross.query(
            f"MAT_alpha_array_1536_barcode == {mat_alpha_array_barcode_filter}"
        )

    wgs_plate_key_path = config.paths.keyfiles / "WGS_plate_key.tsv"
    plate_384_to_96 = config.paths.keyfiles / "384_PIXL_rearray_format.tsv"
    plate_key = load_wgs_plate_key(wgs_plate_key_path, plate_384_to_96)

    wgs_table = pd.read_csv(config.wgs_outcome_table, sep="\t")
    wgs_bc1_call = (
        wgs_table[["sample_name", "bc1_from_wgs"]].drop_duplicates().reset_index(drop=True)
    )

    # Rename plate_key columns to match cross-array output
    plate_key = plate_key.rename(
        columns={
            "MAT_alpha_array_plate_384": "REDI_plate_color_384",
            "row_384": "REDI_row_384",
            "col_384": "REDI_col_384",
        }
    )
    pair = config.array_pairs[0]
    result = compare_to_wgs(
        redi_df=cross,
        wgs_bc1_call_df=wgs_bc1_call,
        wgs_plate_key_df=plate_key,
        array_a_name=pair.name_a,
        array_b_name=pair.name_b,
        bc1_length=config.barcode.bc1_length,
    )
    out = (
        combined_data
        / "WGS_bc1_call_outer_merge_with_REDI_df_and_plate_key.tsv"
    )
    result.to_csv(out, sep="\t", index=False)
    return out
