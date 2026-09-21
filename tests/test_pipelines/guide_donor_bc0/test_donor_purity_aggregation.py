"""Donor purity must be invariant to splitting reads over leader variants."""
import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest

from magestic.pipelines.guide_donor_bc0.core import bc0_purity


@pytest.fixture(params=['core', 'numbered'])
def purity_route(request):
    if request.param == 'core':
        return bc0_purity
    path = Path(bc0_purity.__file__).parents[1]/'snakemake/04_calculate_bc0_purity.py'
    spec = importlib.util.spec_from_file_location('gdb_purity_numbered', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_purity(module, tmp_path, rows, **config):
    source, output = tmp_path/'matched.tsv', tmp_path/'purity.tsv'
    pd.DataFrame(rows, columns=['guide', 'bc0', 'donor', 'oligo_name', 'leader', 'counts']).to_csv(source, sep='\t', index=False)
    stats = module.calculate_bc0_purity(source, output, module.PurityConfig(**config))
    return pd.read_csv(output, sep='\t'), stats, output


def test_split_donor_support(purity_route, tmp_path):
    result, stats, _ = run_purity(purity_route, tmp_path, [
        ('G'*20, 'A'*10, 'C'*129, 'A', 'leader1', 6),
        ('G'*20, 'A'*10, 'C'*129, 'A', 'leader2', 6),
        ('G'*20, 'A'*10, 'T'*129, 'B', 'leader1', 10)])
    row = result.iloc[0]
    assert row.oligo_name == 'A'
    assert row.dominant_donor_counts == 12
    assert row.counts == 6  # Representative sequence-row count remains distinguishable.
    assert row.grouped_guide_bc0_counts == 22
    assert row.guide_bc0_purity == pytest.approx(12/22)
    assert stats['unique_guide_bc0'] == 1


def test_groups_and_representative_evidence(purity_route, tmp_path):
    result, _, _ = run_purity(purity_route, tmp_path, [
        ('G'*20, 'A'*10, 'C'*129, 'A', 'minor_leader', 2),
        ('G'*20, 'A'*10, 'C'*129, 'A', 'major_leader', 8),
        ('G'*20, 'A'*10, 'T'*129, 'B', 'other', 1),
        ('G'*20, 'C'*10, 'T'*129, 'B', 'other_bc0', 4),
        ('T'*20, 'A'*10, 'T'*129, 'B', 'other_guide', 5)])
    assert len(result) == 3
    row = result[(result.guide == 'G'*20) & (result.bc0 == 'A'*10)].iloc[0]
    assert row.leader == 'major_leader'
    assert row.counts == 8
    assert row.dominant_donor_counts == 10
    assert row.guide_bc0_purity == pytest.approx(10/11)
    assert sorted(result.grouped_guide_bc0_counts) == [4, 5, 11]


@pytest.mark.parametrize('min_counts', [1, 2])
def test_empty_replaces_old_output(purity_route, tmp_path, min_counts):
    (tmp_path/'purity.tsv').write_text('old stale data\n')
    rows = [] if min_counts == 1 else [('G'*20, 'A'*10, 'C'*129, 'A', 'leader', 1)]
    result, stats, _ = run_purity(purity_route, tmp_path, rows, min_counts=min_counts)
    assert result.empty
    assert 'dominant_donor_counts' in result.columns
    assert stats['unique_guide_bc0'] == 0


def test_purity_output_records_actual_input_and_thresholds(purity_route, tmp_path):
    _, _, output = run_purity(purity_route, tmp_path, [('G'*20, 'A'*10, 'C'*129, 'A', 'leader', 9)])
    sidecar = json.loads(Path(str(output)+'.provenance.json').read_text())
    assert sidecar['inputs'][0]['md5']
    assert sidecar['outputs'][0]['md5']
    assert sidecar['counts']['rows_out'] == 1
    assert sidecar['counts']['reads_in'] == 9
    assert sidecar['parameters']['min_counts'] == 1
    assert sidecar['parameters']['high_purity_threshold'] == 0.9
    assert sidecar['parameters']['low_purity_threshold'] == 0.5


def test_missing_donor_stays_in_denominator(purity_route, tmp_path):
    result, _, output = run_purity(purity_route, tmp_path, [
        ('G'*20, 'A'*10, 'C'*129, 'A', 'leader', 6),
        ('G'*20, 'A'*10, None, None, 'unknown', 4)])
    assert result.iloc[0].guide_bc0_purity == pytest.approx(0.6)
    assert result.iloc[0].grouped_guide_bc0_counts == 10
    sidecar = json.loads(Path(str(output)+'.provenance.json').read_text())
    assert sidecar['counts']['eligible_group_reads'] == 10


def test_unknown_donors_cannot_outvote_observed_sequence(purity_route, tmp_path):
    result, _, output = run_purity(purity_route, tmp_path, [
        ('G'*20, 'A'*10, None, 'guide_only', 'unknown1', 6),
        ('G'*20, 'A'*10, None, 'guide_only', 'unknown2', 6),
        ('G'*20, 'A'*10, 'C'*129, 'known', 'observed', 10)])
    row = result.iloc[0]
    assert row.donor == 'C'*129
    assert row.oligo_name == 'known'
    assert row.dominant_donor_counts == 10
    assert row.grouped_guide_bc0_counts == 22
    assert row.guide_bc0_purity == pytest.approx(10/22)
    sidecar = json.loads(Path(str(output)+'.provenance.json').read_text())
    assert sidecar['counts']['missing_donor_reads'] == 12


def test_all_unknown_donors_are_unassessable(purity_route, tmp_path):
    result, stats, output = run_purity(purity_route, tmp_path, [
        ('G'*20, 'A'*10, None, 'guide_only', 'unknown1', 6),
        ('G'*20, 'A'*10, None, 'guide_only', 'unknown2', 6)])
    assert len(result) == 1  # Preserve the group and its read denominator.
    row = result.iloc[0]
    assert pd.isna(row.donor) and pd.isna(row.guide_bc0_purity)
    assert row.dominant_donor_counts == 0
    assert row.grouped_guide_bc0_counts == 12
    assert not row.donor_observed
    assert stats['unassessable_groups'] == 1
    assert stats['high_purity_count'] == stats['low_purity_count'] == 0
    assert stats['mean_purity'] is None and stats['median_purity'] is None
    sidecar = json.loads(Path(str(output)+'.provenance.json').read_text())
    assert sidecar['counts']['eligible_group_reads'] == 12


def test_equal_donor_tie_and_sequence_row_threshold(purity_route, tmp_path):
    result, _, output = run_purity(purity_route, tmp_path, [
        ('G'*20, 'A'*10, 'C'*129, 'A', 'first', 4),
        ('G'*20, 'A'*10, 'T'*129, 'B', 'second', 4),
        ('G'*20, 'A'*10, 'T'*129, 'B', 'filtered', 1)], min_counts=2)
    assert result.iloc[0].oligo_name == 'A'
    assert result.iloc[0].guide_bc0_purity == pytest.approx(0.5)
    sidecar = json.loads(Path(str(output)+'.provenance.json').read_text())
    assert sidecar['counts']['filtered_reads'] == 1
    assert sidecar['counts']['reads_in'] == 9


def test_packaged_runner_retains_shared_sidecars(tmp_path):
    from magestic.pipelines.guide_donor_bc0.config import GuideDonorBC0Config
    from magestic.pipelines.guide_donor_bc0.pipeline.run_pipeline import run_bc0_purity

    config = GuideDonorBC0Config(project_name='tiny', project_dir=tmp_path,
                                  output_dir=tmp_path/'output')
    counts_dir = config.get_working_dir('guide_donor_bc0')
    counts_dir.mkdir(parents=True, exist_ok=True)
    source = counts_dir/('V_test'+bc0_purity.INPUT_SUFFIX)
    pd.DataFrame([{'guide': 'G'*20, 'bc0': 'A'*10, 'donor': 'C'*129,
                   'oligo_name': 'A', 'counts': 7}]).to_csv(source, sep='\t', index=False)
    assert run_bc0_purity(config)
    outdir = config.output_dir/'guide_donor_bc0_purity_filtered'
    table = outdir/('V_test'+bc0_purity.OUTPUT_SUFFIX)
    sidecar = json.loads(Path(str(table)+'.provenance.json').read_text())
    assert sidecar['parameters']['donor_key'] == ['guide', 'bc0', 'donor']
    assert sidecar['counts']['reads_in'] == 7
    summary = outdir/'purity_summary.tsv'
    assert pd.read_csv(summary, sep='\t').iloc[0].unique_guide_bc0 == 1
    summary_prov = json.loads(Path(str(summary)+'.provenance.json').read_text())
    assert summary_prov['inputs'][0]['path'] == str(source)
    assert summary_prov['counts']['guide_bc0_groups'] == 1
