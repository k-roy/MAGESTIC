"""Step 07 keeps one tier row per bc1 (register CA01, CA02, UPCL-07-DUPREF).

A bc1 that the reference table links to two oligos used to become two tier rows, so
its t0 reads were counted twice in the tier summary (120 reads reported as 240). It
now abstains: one row, no oligo_name, and both candidates kept.
"""
import importlib.util
import sys
from pathlib import Path

import pandas as pd

import magestic.pipelines.bc1_from_screen as screen

SNAKEMAKE = Path(screen.__file__).parent / 'snakemake'
A, B = 'A' * 20, 'C' * 20


def _step07():
    spec = importlib.util.spec_from_file_location('test_projection_07', SNAKEMAKE / '07_tier_classification.py')
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_multi_link_bc1_is_one_tier_row_and_counted_once(tmp_path, monkeypatch):
    counts = tmp_path / 'counts.tsv'
    pd.DataFrame([dict(bc1=A, sample_name='t0', counts=120),
                  dict(bc1=B, sample_name='t0', counts=50)]).to_csv(counts, sep='\t', index=False)
    reference = tmp_path / 'reference.tsv'
    pd.DataFrame([dict(bc1=A, oligo_name='o1'), dict(bc1=A, oligo_name='o2'),
                  dict(bc1=B, oligo_name='o3'), dict(bc1=B, oligo_name='o3')]).to_csv(
        reference, sep='\t', index=False)
    key = tmp_path / 'key.tsv'
    pd.DataFrame([dict(sample_name='t0', condition='YPD', generations=0)]).to_csv(key, sep='\t', index=False)
    mapping, stats = tmp_path / 'tiers.tsv', tmp_path / 'stats.tsv'
    monkeypatch.setattr(sys, 'argv', ['07', '--bc1-counts', str(counts), '--bc1-reference', str(reference),
                                      '--sample-key', str(key), '--library', 'yLtest',
                                      '--output-tier-mapping', str(mapping), '--output-stats', str(stats)])
    assert _step07().main() == 0
    tiers = pd.read_csv(mapping, sep='\t').set_index('bc1')
    assert len(tiers) == 2
    assert pd.isna(tiers.loc[A, 'oligo_name']) and tiers.loc[A, 'bc1_reference_oligo_candidates'] == 'o1;o2'
    assert tiers.loc[B, 'oligo_name'] == 'o3'
    assert int(pd.read_csv(stats, sep='\t')['total_t0_reads'].sum()) == 170
