"""Tests for the plate-DB snapshot reader."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest


def test_plate_db_resolve_in_range(tmp_path: Path):
    from magestic.pipelines.redi.readers.plate_db import (
        load_plate_db,
        resolve_plate,
        trace_lineage,
    )

    # Build a tiny synthetic snapshot
    source = pd.DataFrame(
        {"plate_id": [1001, 1002], "source_plate_id": [None, 1001]}
    )
    target = pd.DataFrame(
        {"plate_id": [2001], "source_plate_id": [1001]}
    )
    redi = pd.DataFrame(
        {"plate_id": [3001], "source_plate_id": [2001]}
    )
    misc = pd.DataFrame({"plate_id": [4001], "source_plate_id": [None]})
    rearray = pd.DataFrame({"plate_id": [5001], "source_plate_id": [3001]})

    snap = tmp_path
    source.to_parquet(snap / "Source.parquet", index=False)
    target.to_parquet(snap / "Target.parquet", index=False)
    redi.to_parquet(snap / "REDI_mating.parquet", index=False)
    misc.to_parquet(snap / "Misc.parquet", index=False)
    rearray.to_parquet(snap / "Rearray_Target.parquet", index=False)

    db = load_plate_db(snap)
    assert set(db.keys()) == {"Source", "Target", "REDI_mating", "Misc", "Rearray_Target"}

    hit = resolve_plate(2001, db)
    assert hit is not None
    stage, row = hit
    assert stage == "Target"
    assert int(row["plate_id"]) == 2001

    # Lineage: 5001 (Rearray) → 3001 (REDI) → 2001 (Target) → 1001 (Source) → None
    lineage = trace_lineage(5001, db)
    assert lineage == [5001, 3001, 2001, 1001]


def test_plate_db_missing_dir(tmp_path: Path):
    from magestic.pipelines.redi.readers.plate_db import load_plate_db

    with pytest.raises(FileNotFoundError):
        load_plate_db(tmp_path / "does_not_exist")
