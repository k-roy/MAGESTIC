"""
magestic.pipelines.redi.readers.plate_db — agar plate barcode DB reader.

Consumes the snapshot artifacts produced by
`REDI/plate_database/export_plate_db.py` (TSV + per-stage parquet).
**Never** reads the live xlsx (`Rectangular_plate_barcode_database_for_PIXL_Spinnaker.xlsx`)
to avoid multi-user write-contention. If snapshot is missing, the loaders raise.

Sheet → plate ID range:
    Source                1000-1999
    Target                2000-2999
    REDI_mating           3000-3999
    Misc                  4000-4999
    Rearray_Target        5000-5999
    Consolidated_384      (no contiguous range; consolidation IDs)

API:
    load_plate_db(snapshot_dir)        → dict[stage_name → DataFrame]
    resolve_plate(plate_id, db)        → (stage, DataFrame row)
    trace_lineage(plate_id, db)        → list[parent plate IDs walking back]
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

STAGE_NAMES = (
    "Source",
    "Target",
    "REDI_mating",
    "Misc",
    "Rearray_Target",
    "Consolidated_384",
)

PLATE_ID_RANGES: dict[str, tuple[int, int]] = {
    "Source": (1000, 1999),
    "Target": (2000, 2999),
    "REDI_mating": (3000, 3999),
    "Misc": (4000, 4999),
    "Rearray_Target": (5000, 5999),
    # Consolidated_384 has no contiguous range — fall through to explicit lookup
}


def load_plate_db(snapshot_dir: str | Path) -> dict[str, pd.DataFrame]:
    """Load all snapshot parquet files into a dict keyed by stage name."""
    snapshot_dir = Path(snapshot_dir)
    if not snapshot_dir.is_dir():
        raise FileNotFoundError(
            f"Plate-DB snapshot directory not found: {snapshot_dir}. "
            f"Run REDI/plate_database/export_plate_db.py first."
        )
    out: dict[str, pd.DataFrame] = {}
    for stage in STAGE_NAMES:
        candidate = snapshot_dir / f"{stage}.parquet"
        if not candidate.exists():
            # The Consolidated_384 sheet may legitimately be absent in older
            # snapshots; tolerate missing optional sheets but warn.
            if stage == "Consolidated_384":
                continue
            raise FileNotFoundError(
                f"Plate-DB snapshot missing {candidate}. "
                f"Re-run export_plate_db.py."
            )
        out[stage] = pd.read_parquet(candidate)
    return out


def _stage_for_plate_id(plate_id: int) -> str | None:
    for stage, (lo, hi) in PLATE_ID_RANGES.items():
        if lo <= plate_id <= hi:
            return stage
    return None


def resolve_plate(
    plate_id: int,
    db: dict[str, pd.DataFrame],
    plate_id_col: str = "plate_id",
) -> tuple[str, pd.Series] | None:
    """Look up a plate by ID. Tries the range-based stage first; falls back to
    a full-table scan if not found (handles Consolidated_384 and edge cases).

    Returns (stage_name, row) or None if the ID isn't in any sheet.
    """
    stage_guess = _stage_for_plate_id(plate_id)
    if stage_guess and stage_guess in db:
        df = db[stage_guess]
        if plate_id_col in df.columns:
            hit = df[df[plate_id_col] == plate_id]
            if not hit.empty:
                return stage_guess, hit.iloc[0]
    # Fallback: scan all sheets
    for stage, df in db.items():
        if plate_id_col not in df.columns:
            continue
        hit = df[df[plate_id_col] == plate_id]
        if not hit.empty:
            return stage, hit.iloc[0]
    return None


def trace_lineage(
    plate_id: int,
    db: dict[str, pd.DataFrame],
    plate_id_col: str = "plate_id",
    parent_col: str = "source_plate_id",
    max_depth: int = 16,
) -> list[int]:
    """Walk back via `source_plate_id` (or whatever column the snapshot uses)
    to find this plate's lineage. Returns [plate_id, parent_id, grandparent_id, …].

    Stops at `max_depth` to avoid infinite loops on broken DBs.
    """
    out = [plate_id]
    seen = {plate_id}
    for _ in range(max_depth):
        hit = resolve_plate(plate_id, db, plate_id_col=plate_id_col)
        if hit is None:
            break
        _stage, row = hit
        if parent_col not in row.index:
            break
        parent = row.get(parent_col)
        if pd.isna(parent):
            break
        parent_int = int(parent)
        if parent_int in seen:
            break  # cycle guard
        out.append(parent_int)
        seen.add(parent_int)
        plate_id = parent_int
    return out
