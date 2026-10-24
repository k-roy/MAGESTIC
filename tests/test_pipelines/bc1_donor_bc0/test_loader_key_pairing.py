"""Plate x sample key pairing after review of 2d52a83 (patch 29 v2).

File-naming columns shared by both keys are join keys; identity conflicts are refused
for every file they touch, however rows are ordered or de-duplicated; resolved values
keep a clean text spelling; a declared amplicon_level wins over derivation in every
producer; numbered step 02 legacy mode reads keys the way the 2x300 loader does.
"""
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import unquote

import pandas as pd
import pytest

from magestic.pipelines.bc1_donor_bc0.readers import data_loader
from magestic.pipelines.bc1_donor_bc0.readers.data_loader import (
    load_2x150_bc1_donor_bc0_data, load_2x300_data)

PARSED = '20250101_plate_1A_sample_1_bc1_donor_bc0.tsv'
PLATE = dict(description='bc1-donor-bc0', sequencing_date='20250101', outer_primers='plate_1A',
             inner_primers='inner_A', gDNA_plate_name='prep_A', polymerase='Q5')


def _write(tmp_path, plate_rows, sample_rows, table_extra=None):
    data, keys = tmp_path / 'data', tmp_path / 'keys'
    data.mkdir()
    keys.mkdir()
    pd.DataFrame(plate_rows).to_csv(keys / 'plate_key.tsv', sep='\t', index=False)
    pd.DataFrame(sample_rows).to_csv(keys / 'sample_key.tsv', sep='\t', index=False)
    # SPS added at the 218 merge (planning/272): master's NN 250 / C17 refuses to rebuild
    # donor_bc0_fragment without it (the SPS-less key is 40 nt and joins nothing). Fixture only;
    # no assertion changed.
    pd.DataFrame([dict(counts=7, bc1='G' * 20, donor='A' * 129, bc0='C' * 10, SPS='T' * 20,
                       **(table_extra or {}))]).to_csv(data / PARSED, sep='\t', index=False)
    return data, keys


def _decode(identity):
    return json.loads(unquote(identity[len('pcr:v1:'):]))


def _numbered02():
    path = Path(data_loader.__file__).parents[1] / 'snakemake/02_combine_and_collapse.py'
    spec = importlib.util.spec_from_file_location('test_key_pairing_numbered02', path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize('conflicting_first', [True, False])
def test_conflict_is_refused_even_when_a_sibling_row_agrees(tmp_path, conflicting_first):
    rows = [dict(sample_number=1, trafo_ID='trafo_4', yeast_library='yL_B'),
            dict(sample_number=1, trafo_ID='trafo_3', yeast_library='yL_A')]
    if not conflicting_first:
        rows.reverse()
    data, keys = _write(tmp_path, [dict(PLATE, trafo_ID='trafo_3', yeast_library='yL_A')], rows)
    with pytest.raises(ValueError, match='conflicting'):
        load_2x300_data(data, keys)


def test_date_mismatch_without_a_file_does_not_block_the_agreeing_file(tmp_path):
    data, keys = _write(tmp_path, [dict(PLATE, sequencing_date='20250101,20250102')],
                        [dict(sample_number=1, sequencing_date=20250101)])
    result, _ = load_2x300_data(data, keys)
    assert result['counts'].tolist() == [7]


def test_resolved_integer_date_keeps_its_spelling():
    combined = data_loader._cross_merge_keys(
        pd.DataFrame([dict(sequencing_date='20250101'), dict(sequencing_date='20250102')]),
        pd.DataFrame([dict(sample_number=1, sequencing_date=20250101)]))
    assert combined['sequencing_date'].tolist() == ['20250101']


def test_file_named_by_no_key_pair_fails_default_coverage(tmp_path):
    """The only sample row says 20250102; the file on disk is 20250101. No pair names it."""
    data, keys = _write(tmp_path, [PLATE], [dict(sample_number=1, sequencing_date=20250102)])
    with pytest.raises(ValueError):
        load_2x300_data(data, keys)
    result, _ = load_2x300_data(data, keys, require_full_coverage=False)
    assert result.empty


@pytest.mark.parametrize('plate_value,sample_value,conflict', [
    ('1_2', '12', True), ('1e1', '10', True), ('inf', 'Infinity', True),
    ('1', '1.0', False), ('prep_A ', 'prep_A', False)])
def test_what_counts_as_the_same_declaration(plate_value, sample_value, conflict):
    combined = data_loader._cross_merge_keys(
        pd.DataFrame([dict(gDNA_plate_name=plate_value)]),
        pd.DataFrame([dict(sample_number=1, gDNA_plate_name=sample_value)]))
    assert (combined.iloc[0]['_key_conflicts'] == 'gDNA_plate_name') is conflict


@pytest.mark.parametrize('spaced_key', ['plate', 'sample'])
def test_whitespace_in_either_key_still_matches_the_table(tmp_path, spaced_key):
    data, keys = _write(tmp_path, [dict(PLATE, trafo_ID='trafo_3 ' if spaced_key == 'plate' else 'trafo_3')],
                        [dict(sample_number=1, trafo_ID='trafo_3 ' if spaced_key == 'sample' else 'trafo_3')],
                        table_extra={'trafo_ID': 'trafo_3'})
    result, _ = load_2x300_data(data, keys)
    assert result['counts'].tolist() == [7]


def test_declared_amplicon_level_wins_in_every_producer(tmp_path):
    data, keys = _write(tmp_path, [PLATE], [dict(sample_number=1, amplicon_level='pooled_yL')])
    loaded, _ = load_2x300_data(data, keys)
    record = dict(PLATE, sample_number='1', amplicon_level='pooled_yL', sample_id='s1',
                  parsed_path=str(data / PARSED))
    record['plate_polymerase'] = record.pop('polymerase')
    manifest = tmp_path / 'manifest.tsv'
    pd.DataFrame([record]).to_csv(manifest, sep='\t', index=False)
    numbered = _numbered02().build_manifest(
        SimpleNamespace(manifest=manifest, trafo_key=None, allow_missing=False))
    assert _decode(loaded.iloc[0]['pcr_replicate_id']) == _decode(numbered.iloc[0]['pcr_replicate_id'])
    assert _decode(loaded.iloc[0]['pcr_replicate_id'])[5] == 'pooled_yL'

    short_data, short_keys = tmp_path / 'short', tmp_path / 'short_keys'
    short_data.mkdir()
    short_keys.mkdir()
    pd.DataFrame([dict(description='bc1-donor-bc0', sequencing_date=20250101, outer_primers='plate_1A',
                       sample_number=1, inner_primers='inner_A', gDNA_plate_name='prep_A',
                       polymerase='Q5', amplicon_level='pooled_yL')]).to_csv(
        short_keys / 'bc1_donor_bc0_combined_key.tsv', sep='\t', index=False)
    pd.DataFrame([dict(counts=7, bc1='G' * 20, donor_bc0_fragment='A' * 30 + 'C' * 10)]).to_csv(
        short_data / '20250101_plate_1A_1_bc1_donor_bc0_fragment_counts.tsv', sep='\t', index=False)
    short, _ = load_2x150_bc1_donor_bc0_data(short_data, short_keys)
    assert _decode(short.iloc[0]['pcr_replicate_id'])[5] == 'pooled_yL'


@pytest.mark.parametrize('shared_trafo', [False, True])
def test_numbered02_legacy_mode_matches_the_loader(tmp_path, shared_trafo):
    plate = dict(PLATE, trafo_ID='trafo_3') if shared_trafo else dict(PLATE)
    data, keys = _write(tmp_path, [plate], [dict(sample_number=1, trafo_ID='trafo_3')])
    loaded, _ = load_2x300_data(data, keys)
    legacy = _numbered02().build_manifest(SimpleNamespace(
        manifest=None, parsed_dir=str(data), sample_key=str(keys / 'sample_key.tsv'),
        plate_key=str(keys / 'plate_key.tsv'), run_date='20250101', trafo_key=None,
        allow_missing=True))
    assert _decode(loaded.iloc[0]['pcr_replicate_id']) == _decode(legacy.iloc[0]['pcr_replicate_id'])
