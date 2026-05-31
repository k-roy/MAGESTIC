"""
V-library reclone alias mapping helpers.

Per HANDOFF-F.4 Strategy A, the harmonized QTL oligo table labels SPS-derived
plasmid libraries by 2024 originals (V545/V546/V429). The 2025 reclones
(V629/V631/V634) are content-equivalent: same SPS, same oligo content, just
re-prepped in the pS1439 backbone. The harmonized tables carry them in a
``V_aliases`` column as ``V545|V629``, ``V546|V631``, ``V429|V634``.

Screens that use 2025 reclones (yL437/yL442) need to remap to the 2024 V
before joining harmonized annotations.

Mapping (2025 reclone -> 2024 original):
    V629 -> V545   (SpG_NGNG)
    V631 -> V546   (SpG_NGNH)
    V634 -> V429   (NGG; SpG-paired = trafo_3, SpCas9-paired = trafo_4)

Author: Agent D (HANDOFF-D yL437/yL442 reprocessing, 2026-05-30)
"""
from __future__ import annotations

from typing import Dict
import pandas as pd

V_RECLONE_TO_ORIGINAL: Dict[str, str] = {
    "V629": "V545",
    "V631": "V546",
    "V634": "V429",
}

V_ORIGINAL_TO_RECLONES: Dict[str, str] = {
    "V545": "V629",
    "V546": "V631",
    "V429": "V634",
}


def harmonize_v_label(v: str) -> str:
    """Map a single V label (reclone or original) to the 2024 canonical V."""
    return V_RECLONE_TO_ORIGINAL.get(v, v)


def add_v_harmonized_column(
    df: pd.DataFrame,
    source_col: str = "V_library",
    target_col: str = "V_harmonized",
    inplace: bool = False,
) -> pd.DataFrame:
    """Add ``V_harmonized`` column mapping reclones -> 2024 originals.

    Idempotent: if ``target_col`` already exists, returns ``df`` unchanged.
    """
    if target_col in df.columns:
        return df
    out = df if inplace else df.copy()
    out[target_col] = out[source_col].map(V_RECLONE_TO_ORIGINAL).fillna(out[source_col])
    return out
