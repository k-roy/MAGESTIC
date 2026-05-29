"""
Step D — compute REDI clone purity and emit `top_REDI_bc_MAGESTIC_bc1_df.tsv`.

Thin orchestrator: reads the step_c output, calls
`magestic.pipelines.redi.core.purity.compute_purity`, writes the result.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..core.config import REDIConfig
from ..core.purity import compute_purity


def run_step_d(
    config: REDIConfig,
    combined_data_path: Path | None = None,
) -> Path:
    """Run step_d. Returns path to written `top_REDI_bc_MAGESTIC_bc1_df.tsv`."""
    combined_data = config.paths.combined_data
    if combined_data_path is None:
        combined_data_path = combined_data / "combined_REDI_bc_MAGESTIC_bc_merged_df.tsv"

    annotated = pd.read_csv(combined_data_path, sep="\t")
    top_clone = compute_purity(
        annotated,
        group_cols=("sample_name", "REDI_bc"),
        top_n=config.purity.top_n_for_purity,
        edit_distance_drop=config.purity.edit_distance_drop,
    )
    out = combined_data / "top_REDI_bc_MAGESTIC_bc1_df.tsv"
    top_clone.to_csv(out, sep="\t", index=False)
    # 100-row inspection slice (legacy convenience)
    top_clone.head(100).to_csv(
        combined_data / "top_REDI_bc_MAGESTIC_bc1_df_inspection.tsv",
        sep="\t",
        index=False,
    )
    return out
