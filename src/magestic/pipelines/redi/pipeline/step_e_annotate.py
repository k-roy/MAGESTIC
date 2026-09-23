"""
Step E — annotate top-clone table with the bc1-donor-bc0 reference table.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..core.annotate import annotate_with_bc1_reference
from ..core.config import REDIConfig
from ..core.io import provenance, publish_tables
from ..readers.keyfiles import _read_key


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

    out = combined_data / "annotated_top_REDI_bc_MAGESTIC_bc1_df.tsv"
    prov = provenance('e_annotate', out)
    prov.add_input(top_clone_path, role='observed_clone_distributions')
    prov.add_input(config.bc1_reference_table, role='candidate_design_reference')
    top_clone = _read_key(top_clone_path)
    bc1_ref = _read_key(config.bc1_reference_table)
    annotated = annotate_with_bc1_reference(top_clone, bc1_ref)

    prov.set_counts(input_rows=len(top_clone), reference_rows=len(bc1_ref), output_rows=len(annotated),
        ambiguous_rows=int(annotated.bc1_reference_ambiguous.sum()),
        unmatched_rows=int(annotated.bc1_reference_match_count.eq(0).sum()))
    publish_tables([(out, annotated),
        (combined_data / 'annotated_top_REDI_bc_MAGESTIC_bc1_df_head_100.tsv', annotated.head(100))], prov)
    return out
