"""
Step B — combine sample_key + gDNA_plate_key + REDI_plate_barcode_key.

Writes `combined_key.tsv` next to the source keyfiles (legacy behavior).
"""

from __future__ import annotations

from pathlib import Path

from ..core.config import REDIConfig
from ..readers.keyfiles import combine_keyfiles


def run_step_b(config: REDIConfig, out_name: str = "combined_key.tsv") -> Path:
    """Run step_b and return path to the written `combined_key.tsv`."""
    k = config.paths.keyfiles
    combined = combine_keyfiles(
        k / config.sample_key,
        k / config.gdna_plate_key,
        k / config.redi_plate_barcode_key,
    )
    out = k / out_name
    combined.to_csv(out, sep="\t", index=False)
    return out
