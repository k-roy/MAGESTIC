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
from ..core.io import provenance, publish_tables
from ..readers.keyfiles import _read_key


def run_step_d(
    config: REDIConfig,
    combined_data_path: Path | None = None,
) -> Path:
    """Run step_d. Returns path to written `top_REDI_bc_MAGESTIC_bc1_df.tsv`."""
    combined_data = config.paths.combined_data
    if combined_data_path is None:
        combined_data_path = combined_data / "combined_REDI_bc_MAGESTIC_bc_merged_df.tsv"

    out = combined_data / "top_REDI_bc_MAGESTIC_bc1_df.tsv"
    prov = provenance('d_purity', out)
    prov.add_input(combined_data_path, role='registered_barcode_counts')
    annotated = _read_key(combined_data_path)
    top_clone = compute_purity(
        annotated,
        group_cols=("sample_name", "REDI_bc"),
        top_n=config.purity.top_n_for_purity,
        edit_distance_drop=config.purity.edit_distance_drop,
    )
    prov.set_parameters(group_cols=['sample_name', 'REDI_bc'],
        top_n=config.purity.top_n_for_purity, edit_distance_drop=list(config.purity.edit_distance_drop),
        denominator='all valid observed bc1 reads; not cell purity',
        minimum_count_filter='not applied here; all distributions retained for downstream pooling')
    prov.set_counts(input_rows=len(annotated), output_rows=len(top_clone),
        input_reads=sum(int(x) for x in annotated.counts),
        output_reads=sum(int(x) for x in top_clone.total_counts))
    publish_tables([(out, top_clone),
        (combined_data / 'top_REDI_bc_MAGESTIC_bc1_df_inspection.tsv', top_clone.head(100))], prov)
    return out
