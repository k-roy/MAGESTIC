"""Declared scope through the plate x sample key merge, and one PCR identity rule.

Covers the 2x300 loader's resolution of columns that both keyfiles declare, and the
single missing-value rule inside build_pcr_replicate_id that the loaders and numbered
step 02 share. The BDB amplicon reads the guide-donor plasmid; nothing here concerns
genomic integration.
"""
import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from magestic.pipelines.bc1_donor_bc0.readers import data_loader
from magestic.pipelines.bc1_donor_bc0.readers.data_loader import (
    build_pcr_replicate_id, load_2x300_data)

PARSED = '20250101_plate_1A_sample_1_bc1_donor_bc0.tsv'


def write_keys(tmp_path, plate_extra=None, sample_extra=None, table_extra=None,
               sample_rows=None):
    data, keys = tmp_path / 'data', tmp_path / 'keys'
    data.mkdir()
    keys.mkdir()
    plate = dict(description='bc1-donor-bc0', sequencing_date='20250101',
                 outer_primers='plate_1A', inner_primers='inner_A',
                 gDNA_plate_name='prep_A', polymerase='Q5')
    plate.update(plate_extra or {})
    sample = dict(sample_number=1, yeast_library='yL_fixture', trafo_ID='trafo_3')
    sample.update(sample_extra or {})
    samples = sample_rows or [sample]
    pd.DataFrame([plate]).to_csv(keys / 'plate_key.tsv', sep='\t', index=False)
    pd.DataFrame(samples).to_csv(keys / 'sample_key.tsv', sep='\t', index=False)
    pd.DataFrame([dict(counts=7, bc1='G' * 20, donor='A' * 129, bc0='C' * 10,
                       **(table_extra or {}))]).to_csv(data / PARSED, sep='\t', index=False)
    return data, keys


def test_trafo_declared_by_both_keys_and_table_loads(tmp_path):
    """Equal declarations in both keys survive; the table's equal value is accepted."""
    data, keys = write_keys(tmp_path, plate_extra={'trafo_ID': 'trafo_3'},
                            table_extra={'trafo_ID': 'trafo_3'})
    result, _ = load_2x300_data(data, keys)
    assert result['counts'].tolist() == [7]
    assert result.iloc[0]['trafo_ID'] == 'trafo_3'
    assert result.iloc[0]['amplicon_level'] == 'trafo'


@pytest.mark.parametrize('column', ['trafo_ID', 'yeast_library', 'gDNA_plate_name', 'inner_primers'])
def test_conflicting_key_declarations_are_refused(tmp_path, column):
    sample_value = {'gDNA_plate_name': 'prep_A', 'inner_primers': 'inner_A'}.get(column)
    data, keys = write_keys(tmp_path, plate_extra={column: 'conflicting'},
                            sample_extra={column: sample_value} if sample_value else None)
    with pytest.raises(ValueError, match='conflicting'):
        load_2x300_data(data, keys)


@pytest.mark.parametrize('declaring_key', ['plate', 'sample'])
def test_one_sided_declaration_is_kept(tmp_path, declaring_key):
    data, keys = write_keys(tmp_path)
    for side in ('plate', 'sample'):
        path = keys / f'{side}_key.tsv'
        table = pd.read_csv(path, sep='\t')
        table['yeast_sublibrary'] = 'trafo_3' if side == declaring_key else None
        table.to_csv(path, sep='\t', index=False)
    result, _ = load_2x300_data(data, keys)
    assert result.iloc[0]['yeast_sublibrary'] == 'trafo_3'


def test_conflict_on_a_row_without_a_file_does_not_fail_the_load(tmp_path):
    """The cross product also pairs the plate with sample 2, whose file does not exist."""
    data, keys = write_keys(tmp_path, plate_extra={'trafo_ID': 'trafo_3'}, sample_rows=[
        dict(sample_number=1, trafo_ID='trafo_3'), dict(sample_number=2, trafo_ID='trafo_4')])
    result, _ = load_2x300_data(data, keys, require_full_coverage=False)
    assert result['counts'].tolist() == [7]


def test_disagreeing_descriptive_column_keeps_both_values(tmp_path):
    data, keys = write_keys(tmp_path, plate_extra={'notes': 'plate note'},
                            sample_extra={'notes': 'sample note'})
    result, _ = load_2x300_data(data, keys)
    assert result['counts'].tolist() == [7]
    combined = data_loader._cross_merge_keys(pd.read_csv(keys / 'plate_key.tsv', sep='\t'),
                                             pd.read_csv(keys / 'sample_key.tsv', sep='\t'))
    row = combined.iloc[0]
    assert pd.isna(row['notes'])
    assert (row['notes__plate_key'], row['notes__sample_key']) == ('plate note', 'sample note')


def test_pcr_identity_rules():
    base = dict(outer_primers='outer', plate_polymerase='Q5', condition='before',
                amplicon_level='trafo', sample_number=1)
    # Field boundaries are unambiguous even when values contain underscores.
    assert build_pcr_replicate_id(dict(base, gDNA_plate_name='A_B', inner_primers='C')) != \
        build_pcr_replicate_id(dict(base, gDNA_plate_name='A', inner_primers='B_C'))
    # Declared rule: absent, blank and NA spellings are all missing.
    missing = {build_pcr_replicate_id(dict(base, gDNA_plate_name=v, inner_primers='C'))
               for v in (None, float('nan'), '', ' ', 'NA', 'nan')}
    assert len(missing) == 1
    # Formatting differences for one declared PCR encode identically.
    spellings = [dict(base, gDNA_plate_name='A', inner_primers='C', sample_number=n)
                 for n in (1, 1.0, '1', ' 1 ', '01', 'sample_1', 'sample_01')]
    assert len({build_pcr_replicate_id(row) for row in spellings}) == 1
    # A missing amplicon_level is derived from trafo_ID, as the loaders derive it.
    derived = dict(base, gDNA_plate_name='A', inner_primers='C', trafo_ID='trafo_3')
    del derived['amplicon_level']
    assert build_pcr_replicate_id(derived) == build_pcr_replicate_id(
        dict(derived, amplicon_level='trafo'))


def _numbered02():
    path = Path(data_loader.__file__).parents[1] / 'snakemake/02_combine_and_collapse.py'
    spec = importlib.util.spec_from_file_location('test_loader_scope_numbered02', path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_numbered02_and_loader_encode_the_same_declared_pcr(tmp_path):
    data, keys = write_keys(tmp_path)
    loaded, _ = load_2x300_data(data, keys)
    plate = pd.read_csv(keys / 'plate_key.tsv', sep='\t').iloc[0].to_dict()
    sample = pd.read_csv(keys / 'sample_key.tsv', sep='\t').iloc[0].to_dict()
    record = {**plate, **sample, 'sample_id': 's1', 'parsed_path': str(data / PARSED)}
    record['plate_polymerase'] = record.pop('polymerase')
    manifest = tmp_path / 'manifest.tsv'
    pd.DataFrame([record]).to_csv(manifest, sep='\t', index=False)
    numbered = _numbered02().build_manifest(
        SimpleNamespace(manifest=manifest, trafo_key=None, allow_missing=False))
    assert numbered.iloc[0]['pcr_replicate_id'] == loaded.iloc[0]['pcr_replicate_id']


def test_numbered02_missing_fields_match_the_direct_builder(tmp_path):
    parsed = tmp_path / 'parsed.tsv'
    pd.DataFrame([dict(counts=7, bc1='G' * 20, donor='A' * 129, bc0='C' * 10)]).to_csv(
        parsed, sep='\t', index=False)
    row = dict(sample_id='s1', parsed_path=str(parsed), sequencing_date='20250101',
               sample_number='1', outer_primers='plate_1A')
    manifest = tmp_path / 'manifest.tsv'
    pd.DataFrame([row]).to_csv(manifest, sep='\t', index=False)
    numbered = _numbered02().build_manifest(
        SimpleNamespace(manifest=manifest, trafo_key=None, allow_missing=False))
    assert numbered.iloc[0]['pcr_replicate_id'] == build_pcr_replicate_id(row)
