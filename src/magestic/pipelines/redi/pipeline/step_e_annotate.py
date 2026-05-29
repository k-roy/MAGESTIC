"""
Step E — annotate top-clone table with the bc1-donor-bc0 reference table.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..core.annotate import annotate_with_bc1_reference
from ..core.config import REDIConfig


def run_step_e(
    config: REDIConfig,
    top_clone_path: Path | None = None,
) -> Path:
    """Run step_e. Returns path to `annotated_top_REDI_bc_MAGESTIC_bc1_df.tsv`."""
    if config.bc1_reference_table is None:
        raise ValueError(
            "config.bc1_reference_table is None — step_e requires the "
            "bc1_donor_bc0 reference table path."
        )

    combined_data = config.paths.combined_data
    if top_clone_path is None:
        top_clone_path = combined_data / "top_REDI_bc_MAGESTIC_bc1_df.tsv"

    top_clone = pd.read_csv(top_clone_path, sep="\t")
    bc1_ref = pd.read_csv(config.bc1_reference_table, sep="\t")
    annotated = annotate_with_bc1_reference(top_clone, bc1_ref)

    out = combined_data / "annotated_top_REDI_bc_MAGESTIC_bc1_df.tsv"
    annotated.to_csv(out, sep="\t", index=False)
    annotated.head(100).to_csv(
        combined_data / "annotated_top_REDI_bc_MAGESTIC_bc1_df_head_100.tsv",
        sep="\t",
        index=False,
    )
    return out
