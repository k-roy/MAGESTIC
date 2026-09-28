"""Stage I worklist writer: deterministic targets, PIXL run files, and a round trip that catches drift."""
import pandas as pd
import pytest

from magestic.pipelines.redi.core.worklist import (WorklistPlan, assign, edge_wells, parse_run_file,
                                                     run_file_lines, wells)
from magestic.pipelines.redi.pipeline.step_i_worklist import run_step_i_worklist, validate_generation


def picks(spec, density='384'):
    """spec: {plate: n} -> the first n wells of each plate in row order."""
    rows = []
    for plate, n in spec.items():
        for r, c in wells(int(density))[:n]:
            rows.append(dict(source_plate=plate, source_row=r, source_col=str(c), source_density=density,
                             bc1='ACGT' * 5))
    return pd.DataFrame(rows)


def test_row_major_fill_spills_to_next_target():
    t = assign(picks({'2351': 6, '2350': 384}), WorklistPlan(target_plates=['5169', '5170']))
    assert list(t.source_plate[:2]) == ['2350', '2350']            # sorted by source plate first
    assert (t.target_plate.iloc[0], t.target_row.iloc[0], t.target_col.iloc[0]) == ('5169', 'A', 1)
    assert (t.target_row.iloc[1], t.target_col.iloc[1]) == ('A', 2)
    assert (t.target_plate.iloc[383], t.target_row.iloc[383], t.target_col.iloc[383]) == ('5169', 'P', 24)
    assert (t.target_plate.iloc[384], t.target_row.iloc[384], t.target_col.iloc[384]) == ('5170', 'A', 1)


def test_runs_hold_at_most_two_sources_and_fill_the_empty_nest():
    t = assign(picks({'2350': 10, '2351': 10, '2352': 10}), WorklistPlan(target_plates=['5169']))
    assert t.groupby('run').source_plate.nunique().max() == 2 and t.run.nunique() == 2
    last = run_file_lines(t, 2, 2)
    assert last[:3] == ['5169,SBS,384,Target', '2352,SBS,384,Source', '0000,SBS,384,Source']
    assert last[3] == '2352,A,1,5169,A,21'


def test_legacy_header_carries_the_source_density():
    t = assign(picks({'4860': 3}, density='1536'), WorklistPlan(target_plates=['5102']))
    assert run_file_lines(t, 1, 2)[:2] == ['5102,SBS,384,Target', '4860,SBS,1536,Source']


def test_edge_exclusion_starts_at_B2():
    t = assign(picks({'2350': 5}), WorklistPlan(target_plates=['5169'], exclude_target_wells=edge_wells(384)))
    assert (t.target_row.iloc[0], t.target_col.iloc[0]) == ('B', 2)
    assert not any((r, c) in edge_wells(384) for r, c in zip(t.target_row, t.target_col))


@pytest.mark.parametrize('bad', [dict(source_row='Q'), dict(source_col='25'), dict(source_plate='0000'),
                                 dict(source_density='500')])
def test_invalid_source_coordinates_are_refused(bad):
    p = picks({'2350': 2})
    for k, v in bad.items():
        p.loc[0, k] = v
    with pytest.raises(ValueError):
        assign(p, WorklistPlan(target_plates=['5169']))


def test_1536_source_coordinates_are_accepted():
    p = pd.DataFrame([dict(source_plate='4800', source_row='AF', source_col='48', source_density='1536')])
    assert assign(p, WorklistPlan(target_plates=['5169'])).target_row.iloc[0] == 'A'


def test_duplicate_source_well_and_overfull_targets_are_refused():
    p = picks({'2350': 2})
    with pytest.raises(ValueError, match='same source well'):
        assign(pd.concat([p, p.iloc[:1]]), WorklistPlan(target_plates=['5169']))
    with pytest.raises(ValueError, match='exceed'):
        assign(picks({'2350': 384, '2351': 1}), WorklistPlan(target_plates=['5169']))


def test_restamped_source_is_renamed_and_the_original_kept():
    t = assign(picks({'2350': 2}), WorklistPlan(target_plates=['5169'], source_rename={'2350': '7201'}))
    assert set(t.source_plate) == {'7201'} and set(t.original_source_plate) == {'2350'}


def write(tmp_path, spec, **plan):
    src = tmp_path / 'picks.tsv'
    picks(spec).to_csv(src, sep='\t', index=False)
    out = tmp_path / 'gen1'
    run_step_i_worklist(src, WorklistPlan(target_plates=['5169', '5170'], **plan), out, 'TEST')
    return src, out


def test_round_trip_passes_and_every_run_file_matches_the_master(tmp_path):
    _, out = write(tmp_path, {'2350': 200, '2351': 200, '2352': 50})
    report = validate_generation(out)
    assert report['ok'] and report['transfers_in_files'] == report['transfers_in_master'] == 450
    assert (out / 'TEST_transfers.tsv.provenance.json').exists()
    header, body = parse_run_file(out / 'TEST_batch_1_1.csv')
    assert header[0] == ('5169', '384', 'Target') and len(body) > 0


def test_a_tampered_run_file_is_caught(tmp_path):
    _, out = write(tmp_path, {'2350': 20})
    f = out / 'TEST_batch_1_1.csv'
    lines = f.read_text().splitlines()
    lines[3] = lines[3].rsplit(',', 1)[0] + ',24'                  # move one transfer to another target well
    f.write_text('\n'.join(lines) + '\n')
    assert not validate_generation(out)['ok']


def test_a_generation_directory_is_never_rewritten(tmp_path):
    src, out = write(tmp_path, {'2350': 5})
    with pytest.raises(FileExistsError):
        run_step_i_worklist(src, WorklistPlan(target_plates=['5169']), out, 'TEST')
