"""
Step B — combine sample_key + gDNA_plate_key + REDI_plate_barcode_key.

Writes `combined_key.tsv` next to the source keyfiles (legacy behavior).
"""

from __future__ import annotations

from pathlib import Path

from ..core.config import REDIConfig
from ..readers.keyfiles import combine_keyfiles
from ..core.io import provenance, publish_tables


def run_step_b(config: REDIConfig, out_name: str = "combined_key.tsv") -> Path:
    """Run step_b and return path to the written `combined_key.tsv`."""
    k = config.paths.keyfiles
    if Path(out_name).name != out_name or out_name in ('', '.', '..'):
        raise ValueError('key output name must be a safe basename')
    out = k / out_name
    prov = provenance('b_combine_keys', out)
    for name in (config.sample_key, config.gdna_plate_key, config.redi_plate_barcode_key):
        prov.add_input(k / name, role='physical_assay_registration')
    combined = combine_keyfiles(
        k / config.sample_key,
        k / config.gdna_plate_key,
        k / config.redi_plate_barcode_key,
    )
    prov.set_counts(output_rows=len(combined))
    publish_tables([(out, combined)], prov)
    return out
