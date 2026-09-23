"""Modern snapshot completion must be checked without rejecting legacy data.

Two rejection cases are deliberately separate from the frozen GHI after suite.
Their five primary plate records are declared here, before export; no observed
barcode table supplies identity. Neither test edits the shared input workbook.
"""
import importlib.util
from pathlib import Path
import shutil

import pandas as pd
import pytest

from magestic.pipelines.redi.readers.plate_db import load_plate_db, resolve_plate


PRIMARY_SHEETS = [
    ('Source 1000-1999', 'rectangular plate barcode labels for PIXL-Spinnaker', 1001),
    ('Target 2000-2999', 'target_plate_ID', 2001),
    ('REDI_mating 3000-3999', 'barcode', 3001),
    ('misc 4000-4999', 'PlusPlate barcode label', 4001),
    ('Rearray_Target 5000-5999', '384 YPD+G418 agar plate barcode', 5001),
]
EXPECTED_IDS = {'Source': 1001, 'Target': 2001, 'REDI_mating': 3001,
                'Misc': 4001, 'Rearray_Target': 5001}


def declared_export(tmp_path):
    module_path = Path(__file__).resolve().parents[3]/'REDI/plate_database/export_plate_db.py'
    spec = importlib.util.spec_from_file_location('redi_plate_export_attestation_test', module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    workbook = tmp_path/'declared_primary.xlsx'
    with pd.ExcelWriter(workbook, engine='openpyxl') as writer:
        for sheet, header, identifier in PRIMARY_SHEETS:
            pd.DataFrame({header: [identifier]}).to_excel(writer, sheet_name=sheet, index=False)
    destination = tmp_path/'modern_snapshot'
    exported = module.export(workbook, destination)
    assert set(exported)-{'__union_tsv'} == set(EXPECTED_IDS)
    assert (destination/'plates.tsv.provenance.json').is_file()
    assert (destination/'plates.tsv.lock').is_file()

    # A historic, unattested snapshot has neither modern completion marker.
    # Copy exactly the five declared tables; do not invent a provenance file.
    legacy = tmp_path/'legacy_snapshot'/'plates'
    legacy.mkdir(parents=True)
    for stage in EXPECTED_IDS:
        shutil.copyfile(destination/'plates'/f'{stage}.parquet', legacy/f'{stage}.parquet')
    assert not (legacy.parent/'plates.tsv.provenance.json').exists()
    assert not (legacy.parent/'plates.tsv.lock').exists()
    old = load_plate_db(legacy)
    assert set(old) == set(EXPECTED_IDS)
    for stage, identifier in EXPECTED_IDS.items():
        assert len(old[stage]) == 1
        assert resolve_plate(identifier, old)[0] == stage
    return destination


def test_modern_plate_snapshot_without_completion_record_is_rejected(tmp_path):
    snapshot = declared_export(tmp_path)
    # Persistent writer lock identifies the modern publication protocol. A
    # missing sidecar can be an interrupted generation, never legacy evidence.
    (snapshot/'plates.tsv.provenance.json').unlink()
    with pytest.raises(ValueError, match='attest|provenance|incomplete|completion'):
        load_plate_db(snapshot/'plates')


def test_modern_plate_snapshot_with_modified_parquet_is_rejected(tmp_path):
    snapshot = declared_export(tmp_path)
    target = snapshot/'plates'/'Target.parquet'
    frame = pd.read_parquet(target)
    # Identity remains syntactically valid. Only content-attestation checking
    # can detect this changed record; ordinary ID validation must not mask it.
    frame['descriptive_note'] = 'changed after completed export'
    frame.to_parquet(target, index=False)
    with pytest.raises(ValueError, match='attest|provenance|hash|changed|content'):
        load_plate_db(snapshot/'plates')
