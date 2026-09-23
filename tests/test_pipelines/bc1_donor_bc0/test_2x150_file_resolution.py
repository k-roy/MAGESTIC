"""A combined-key row loads its own file or none (218 patch 42).

The 2x150 loader and the donor-bc0 lookup tried the exact ``{date}_{outer}_{sample}_...`` name
and then a loose ``*_{sample}_*...`` glob, which accepts a single match from any other date or
plate. On the 20250912 yL437/yL442 layout that fallback resolves 577 of the 673 bc1-donor-bc0
key rows to a file named for another row, so each of the 96 files on disk loads under up to
seven rows with different metadata and different PCR identities: reads multiply and PCR support
inflates, silently (probe: planning/218_magestic_production_corrections/probes/
fallback_20250912.json). The loose glob now runs only when a row cannot name its exact file,
exactly as in the 2x300 loader. One physical file claimed by two rows with different metadata is
refused, and the 2x150 loader checks that every fragment table on disk was loaded.
"""
import os

import pandas as pd
import pytest

from magestic.pipelines.bc1_donor_bc0.readers.data_loader import (
    load_2x150_bc1_donor_bc0_data, load_donor_bc0_lookup)

SPS = 'CGATAGACGACTGGACAGCA'
BC1, BC0 = 'G' * 20, 'A' * 10
FRAGMENT = 'T' * 30 + SPS + BC0


def _key_row(date, sample, trafo='trafo_1', description='bc1-donor-bc0', outer='plate_2A', **extra):
    return dict(sequencing_date=date, outer_primers=outer, description=description,
                sample_number=sample, condition='plate wash', trafo_ID=trafo, **extra)


def _fragment_table(path, count):
    pd.DataFrame([dict(counts=count, bc1=BC1, bc0=BC0, SPS=SPS, donor_bc0_fragment=FRAGMENT)]).to_csv(
        path, sep='\t', index=False)


def _layout(tmp_path, rows, files):
    data, keys = tmp_path / 'data', tmp_path / 'keys'
    data.mkdir()
    keys.mkdir()
    pd.DataFrame(rows).to_csv(keys / 'bc1_donor_bc0_combined_key.tsv', sep='\t', index=False)
    for name, count in files.items():
        _fragment_table(data / name, count)
    return data, keys


def test_a_row_whose_file_is_absent_loads_nothing(tmp_path):
    """Only 20251125 was processed into this folder; the 20251212 re-sequencing rows name no file."""
    data, keys = _layout(tmp_path,
                         [_key_row('20251125', 'sample_1'), _key_row('20251212', 'sample_1', trafo='yL442'),
                          _key_row('20251212', 'sample_1', outer='plate_3B')],
                         {'20251125_plate_2A_sample_1_bc1_donor_bc0_fragment_counts.tsv': 7})
    frame, info = load_2x150_bc1_donor_bc0_data(data, keys)
    assert frame['counts'].sum() == 7
    assert len(frame) == 1 and frame.iloc[0]['sequencing_date'] == 20251125
    assert len(info.file_paths) == 1


def test_a_key_without_dates_still_uses_the_sample_glob(tmp_path):
    """Control: with no sequencing_date or outer_primers there is no exact name to try."""
    data, keys = tmp_path / 'data', tmp_path / 'keys'
    data.mkdir()
    keys.mkdir()
    pd.DataFrame([dict(description='bc1-donor-bc0', sample_number=4, trafo_ID='trafo_1')]).to_csv(
        keys / 'bc1_donor_bc0_combined_key.tsv', sep='\t', index=False)
    _fragment_table(data / 'run1_plate_2A_4_bc1_donor_bc0_fragment_counts.tsv', 9)
    frame, _ = load_2x150_bc1_donor_bc0_data(data, keys)
    assert frame['counts'].sum() == 9


def test_one_file_under_two_conflicting_rows_is_refused(tmp_path):
    data, keys = _layout(tmp_path,
                         [_key_row('20251125', 'sample_1'), _key_row('20251125', 'sample_1', trafo='yL442')],
                         {'20251125_plate_2A_sample_1_bc1_donor_bc0_fragment_counts.tsv': 7})
    with pytest.raises(ValueError, match='one physical file'):
        load_2x150_bc1_donor_bc0_data(data, keys)


def test_a_repeated_identical_row_loads_its_file_once(tmp_path):
    data, keys = _layout(tmp_path, [_key_row('20251125', 'sample_1'), _key_row('20251125', 'sample_1')],
                         {'20251125_plate_2A_sample_1_bc1_donor_bc0_fragment_counts.tsv': 7})
    frame, info = load_2x150_bc1_donor_bc0_data(data, keys)
    assert frame['counts'].sum() == 7 and len(info.file_paths) == 1


def test_a_fragment_table_no_row_names_fails_coverage(tmp_path):
    data, keys = _layout(tmp_path, [_key_row('20251125', 'sample_1')],
                         {'20251125_plate_2A_sample_1_bc1_donor_bc0_fragment_counts.tsv': 7,
                          '20251125_plate_2A_sample_2_bc1_donor_bc0_fragment_counts.tsv': 3})
    with pytest.raises(ValueError, match='COVERAGE FAILURE'):
        load_2x150_bc1_donor_bc0_data(data, keys)
    frame, _ = load_2x150_bc1_donor_bc0_data(data, keys, require_full_coverage=False)
    assert frame['counts'].sum() == 7


def test_donor_lookup_row_whose_file_is_absent_loads_nothing(tmp_path):
    """A donor-bc0 row for an unprocessed date must not borrow another date's table."""
    data, keys = tmp_path / 'data', tmp_path / 'keys'
    data.mkdir()
    keys.mkdir()
    pd.DataFrame([_key_row('20251125', 'sample_1', description='donor-bc0', outer='plate_3A'),
                  _key_row('20251212', 'sample_1', trafo='trafo_2', description='donor-bc0', outer='plate_3A')]).to_csv(
        keys / 'bc1_donor_bc0_combined_key.tsv', sep='\t', index=False)
    pd.DataFrame([dict(counts=11, donor='T' * 130, donor_bc0_fragment=FRAGMENT, bc0=BC0)]).to_csv(
        data / '20251125_plate_3A_sample_1_complete_donor_bc0_counts.tsv', sep='\t', index=False)
    lookup = load_donor_bc0_lookup(data, keys)
    assert int(lookup['counts'].sum()) == 11
    assert set(lookup['trafo_ID']) == {'trafo_1'}
    assert os.path.basename(lookup.iloc[0]['donor_recovery_source_files'].strip('[]"')) == \
        '20251125_plate_3A_sample_1_complete_donor_bc0_counts.tsv'
