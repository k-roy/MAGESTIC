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

import fcntl
import hashlib
import io
import json
from contextlib import nullcontext
from pathlib import Path
from decimal import Decimal, InvalidOperation

import pandas as pd

STAGE_NAMES = (
    "Source",
    "Target",
    "REDI_mating",
    "Misc",
    "Rearray_Target",
    "Consolidated_384",
    "Misc_7000",
)

# These are primary A1 headers in the archived workbooks/exporter. Rearray's
# side-by-side Target_plate_ID/first_source_plate table is deliberately excluded.
PRIMARY_ID_HEADERS = {
    "Source": ("rectangular plate barcode labels for PIXL-Spinnaker",),
    "Target": ("target_plate_ID",),
    "REDI_mating": ("barcode",),
    "Misc": ("PlusPlate barcode label", "rectangular plate barcode labels for PIXL-Spinnaker"),
    "Rearray_Target": ("384 YPD+G418 agar plate barcode",),
    "Consolidated_384": ("PlusPlate barcode label", "rectangular plate barcode labels for PIXL-Spinnaker"),
    "Misc_7000": ("PlusPlate_barcode_label",),
}


def _plate_integer(value) -> int:
    try:
        if isinstance(value, bool):
            raise ValueError
        number = Decimal(str(value))
        if not number.is_finite() or number != number.to_integral_value() or number <= 0:
            raise ValueError
        return int(number)
    except (InvalidOperation, ValueError, TypeError):
        raise ValueError(f"plate identifier must be a positive exact integer: {value!r}") from None


def _missing(value) -> bool:
    return pd.isna(value) or str(value).strip() in ("", "nan", "None", "<NA>")


def _plate_frame(stage: str, frame: pd.DataFrame, id_col: str) -> pd.DataFrame:
    if frame.columns.duplicated().any():
        raise ValueError(f"{stage}: duplicate plate snapshot column names")
    candidates = [id_col] + (list(PRIMARY_ID_HEADERS.get(stage, ())) if id_col == "plate_id" else [])
    candidates = list(dict.fromkeys(c for c in candidates if c in frame))
    if not candidates:
        raise ValueError(f"{stage}: missing primary plate identifier column {id_col!r}")
    identities = []
    for _, row in frame.iterrows():
        values = {_plate_integer(row[c]) for c in candidates if not _missing(row[c])}
        if len(values) > 1:
            raise ValueError(f"{stage}: conflicting primary plate identifier columns")
        identities.append(next(iter(values)) if values else pd.NA)
    out = frame.copy()
    out[id_col] = pd.Series(identities, index=out.index, dtype="Int64")
    if out.loc[out[id_col].notna(), id_col].duplicated().any():
        raise ValueError(f"{stage}: duplicate plate identifier; physical lookup is ambiguous")
    return out

PLATE_ID_RANGES: dict[str, tuple[int, int]] = {
    "Source": (1000, 1999),
    "Target": (2000, 2999),
    "REDI_mating": (3000, 3999),
    "Misc": (4000, 4999),
    "Rearray_Target": (5000, 5999),
    # Consolidated_384 has no contiguous range — fall through to explicit lookup
}


def _attested_parquets(snapshot_dir: Path) -> dict[str, str] | None:
    """`{parquet file name: md5}` from a modern export's completion record.

    The exporter publishes under an exclusive `plates.tsv.lock` and writes
    `plates.tsv.provenance.json` last (`core.io.publish_tables`). Either marker
    beside `plates/` identifies that protocol, so a missing record means an
    interrupted or damaged export, never a legacy snapshot. A snapshot with
    neither marker predates the protocol and returns None (loaded as before).
    """
    flat = snapshot_dir.parent / "plates.tsv"
    sidecar = Path(str(flat) + ".provenance.json")
    if not sidecar.exists():
        if Path(str(flat) + ".lock").exists():
            raise ValueError(
                f"Plate-DB snapshot {snapshot_dir} uses the attested export protocol but has no "
                f"completion provenance ({sidecar}); the export is incomplete. Re-run export_plate_db.py."
            )
        return None
    try:
        record = json.loads(sidecar.read_text())
        rows = record["outputs"]
        attested = {Path(row["path"]).name: row["md5"] for row in rows
                    if Path(row["path"]).parent.name == snapshot_dir.name and row["path"].endswith(".parquet")}
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise ValueError(f"Plate-DB completion provenance {sidecar} is unreadable: {exc}") from None
    if not attested:
        raise ValueError(f"Plate-DB completion provenance {sidecar} attests no parquet in {snapshot_dir}")
    return attested


def load_plate_db(snapshot_dir: str | Path) -> dict[str, pd.DataFrame]:
    """Load all snapshot parquet files into a dict keyed by stage name.

    A modern export is read under the exporter's lock (shared) and every parquet is
    checked against the completion record's MD5 before use; the frame is parsed
    from the same bytes that were hashed. A parquet the record does not attest, or
    an attested one that is missing or changed, is refused. Legacy snapshots with
    no protocol markers load unchanged.
    """
    snapshot_dir = Path(snapshot_dir)
    if not snapshot_dir.is_dir():
        raise FileNotFoundError(
            f"Plate-DB snapshot directory not found: {snapshot_dir}. "
            f"Run REDI/plate_database/export_plate_db.py first."
        )
    lock_path = snapshot_dir.parent / "plates.tsv.lock"
    with (lock_path.open("r") if lock_path.exists() else nullcontext()) as lock:
        if lock is not None:
            fcntl.flock(lock, fcntl.LOCK_SH)
        attested = _attested_parquets(snapshot_dir)
        if attested is not None:
            present = {p.name for p in snapshot_dir.glob("*.parquet")}
            unattested = sorted(present - set(attested))
            missing = sorted(set(attested) - present)
            if unattested or missing:
                raise ValueError(
                    f"Plate-DB snapshot {snapshot_dir} does not match its completion provenance: "
                    f"unattested={unattested}, missing={missing}"
                )
        out: dict[str, pd.DataFrame] = {}
        for stage in STAGE_NAMES:
            candidate = snapshot_dir / f"{stage}.parquet"
            if not candidate.exists():
                # These sheets may legitimately be absent in older snapshots.
                if stage in ("Consolidated_384", "Misc_7000"):
                    continue
                raise FileNotFoundError(
                    f"Plate-DB snapshot missing {candidate}. "
                    f"Re-run export_plate_db.py."
                )
            data = candidate.read_bytes()
            if attested is not None and hashlib.md5(data).hexdigest() != attested[candidate.name]:
                raise ValueError(
                    f"Plate-DB snapshot content changed after its attested export: {candidate} "
                    "no longer matches the completion provenance hash"
                )
            out[stage] = _plate_frame(stage, pd.read_parquet(io.BytesIO(data)), "plate_id")
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
    """Look up an exact primary plate ID across all registered snapshot sheets.

    Ranges never suppress a competing hit in another sheet. Duplicate identity
    is ambiguous even when one of its records lies in the nominal range.

    Returns (stage_name, row) or None if the ID isn't in any sheet.
    """
    plate_id = _plate_integer(plate_id)
    hits = []
    for stage, df in db.items():
        normalized = _plate_frame(stage, df, plate_id_col)
        hit = normalized[normalized[plate_id_col] == plate_id]
        if not hit.empty:
            hits.append((stage, hit.iloc[0]))
    if len(hits) > 1:
        raise ValueError(f"duplicate plate identifier {plate_id} across sheets; lookup is ambiguous")
    return hits[0] if hits else None


def trace_lineage(
    plate_id: int,
    db: dict[str, pd.DataFrame],
    plate_id_col: str = "plate_id",
    parent_col: str = "source_plate_id",
    max_depth: int = 16,
) -> list[int]:
    """Walk back via `source_plate_id` (or whatever column the snapshot uses)
    to find this plate's lineage. Returns [plate_id, parent_id, grandparent_id, …].

    Only an explicit single-parent column is followed. Raw mating/quadruplet
    descriptions require a separately registered graph, not row-adjacent fields.
    Cycles, unresolved parent records and a depth bound raise instead of returning
    a truncated chain that appears complete.
    """
    if type(max_depth) is not int or max_depth < 1:
        raise ValueError("max_depth must be a positive integer")
    plate_id = _plate_integer(plate_id)
    out = [plate_id]
    seen = {plate_id}
    while True:
        hit = resolve_plate(plate_id, db, plate_id_col=plate_id_col)
        if hit is None:
            raise ValueError(f"unresolved plate in lineage: {plate_id}")
        _stage, row = hit
        if parent_col not in row.index:
            return out
        parent = row.get(parent_col)
        if _missing(parent):
            return out
        parent_int = _plate_integer(parent)
        if parent_int in seen:
            raise ValueError(f"cycle in physical plate lineage at {parent_int}")
        if len(out) - 1 >= max_depth:
            raise ValueError("physical plate lineage exceeds max_depth; chain is incomplete")
        out.append(parent_int)
        seen.add(parent_int)
        plate_id = parent_int
