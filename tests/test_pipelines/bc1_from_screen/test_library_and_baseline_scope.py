"""Baselines are scoped to one yeast library and to generation exactly 0.

Register findings CT01 (step 07 --library did not scope the t0 baseline), UP-C008
(an unknown library name disabled step 08's library filter even when the key lists
other libraries) and UP-C009 (generation 0.5 was truncated to 0 and treated as t0).
Behaviour is checked through the numbered scripts themselves.
"""
import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

import magestic.pipelines.bc1_from_screen as screen
from magestic.pipelines.bc1_from_screen.core.tier_classification import calculate_t0_counts

SNAKEMAKE = Path(screen.__file__).parent / 'snakemake'
BARCODE = 'A' * 20


def _step(name):
    spec = importlib.util.spec_from_file_location(f'test_scope_{name[:2]}', SNAKEMAKE / name)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _mixed_key():
    return pd.DataFrame([
        dict(sample_name='t0_437', condition='YPD', generations=0, screen_date='d1', library_ID='yL437'),
        dict(sample_name='t0_442', condition='YPD', generations=0, screen_date='d1', library_ID='yL442'),
    ])


def test_step07_baseline_uses_only_the_requested_library(tmp_path, monkeypatch):
    """10 reads in the yL437 baseline plus 100 in yL442's must stay 10 reads (tier low) for yL437."""
    counts = tmp_path / 'counts.tsv'
    pd.DataFrame([dict(bc1=BARCODE, sample_name='t0_437', counts=10),
                  dict(bc1=BARCODE, sample_name='t0_442', counts=100)]).to_csv(counts, sep='\t', index=False)
    reference = tmp_path / 'reference.tsv'
    pd.DataFrame([dict(bc1=BARCODE, oligo_name='o1')]).to_csv(reference, sep='\t', index=False)
    key = tmp_path / 'key.tsv'
    _mixed_key().to_csv(key, sep='\t', index=False)
    mapping, stats = tmp_path / 'tiers.tsv', tmp_path / 'stats.tsv'
    step07 = _step('07_tier_classification.py')
    monkeypatch.setattr(sys, 'argv', ['07', '--bc1-counts', str(counts), '--bc1-reference', str(reference),
                                      '--sample-key', str(key), '--library', 'yL437',
                                      '--output-tier-mapping', str(mapping), '--output-stats', str(stats)])
    assert step07.main() == 0
    row = pd.read_csv(mapping, sep='\t').iloc[0]
    assert (int(row['t0_total']), row['tier']) == (10, 'abundance_low')


def test_step08_refuses_a_library_the_key_does_not_list():
    """A typo must not silently disable the filter when the key lists yL437 and yL442."""
    with pytest.raises(ValueError, match='yL473'):
        _step('08_create_count_matrix.py').resolve_library_samples(_mixed_key(), 'yL473')


@pytest.mark.parametrize('key,expected', [
    (_mixed_key(), ['t0_437']),
    (_mixed_key().drop(columns='library_ID'), None),              # single-library screen
    (_mixed_key().assign(library_ID=['L1', 'L2']), None),         # column holds no yeast libraries
], ids=['declared_library', 'no_library_column', 'non_yeast_library_values'])
def test_step08_library_resolution_controls(key, expected):
    assert _step('08_create_count_matrix.py').resolve_library_samples(key, 'yL437') == expected


def test_steps_07_and_08_share_one_library_rule():
    assert _step('07_tier_classification.py').resolve_library_samples is \
        _step('08_create_count_matrix.py').resolve_library_samples


def test_fractional_generation_is_not_a_baseline():
    key = pd.DataFrame([
        dict(sample_name='g0', condition='YPD', generations=0, screen_date='d1'),
        dict(sample_name='g0_text', condition='YPD', generations='0', screen_date='d1'),
        dict(sample_name='g0.5', condition='YPD', generations=0.5, screen_date='d1'),
        dict(sample_name='g_bad', condition='YPD', generations='none', screen_date='d1'),
    ])
    assert sorted(_step('07_tier_classification.py').identify_t0_samples(key)) == ['g0', 'g0_text']
    long_table = pd.DataFrame([dict(bc1=BARCODE, sample_name=name, sample_counts=5)
                               for name in key['sample_name']])
    t0 = calculate_t0_counts(long_table, key, SimpleNamespace(name='yLtest', screen_dates=['d1']))
    assert int(t0['t0_total'].sum()) == 10
