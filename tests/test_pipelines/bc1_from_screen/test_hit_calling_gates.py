"""Hit-calling gates and library concordance (register SCR-D05a/b, SCR-D08a/b, UP-C009).

A null with no spread made every z infinite, so any nonzero effect passed any
threshold; the advertised concordance gate was parsed and never applied; library
concordance counted rows and could exceed 1; and the NGNH exception fired when NGNH
was the only library calling a variant. A 0.5 generation also still joined the
DESeq2 reference in core/count_matrix.py.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import magestic.pipelines.bc1_from_screen as screen
from magestic.pipelines.bc1_from_screen.core.aggregation import aggregate_across_nuclease_libraries
from magestic.pipelines.bc1_from_screen.core.count_matrix import split_by_comparison

SNAKEMAKE = Path(screen.__file__).parent / 'snakemake'


def _step14():
    spec = importlib.util.spec_from_file_location('test_gates_14', SNAKEMAKE / '14_hit_calling.py')
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _hits(tmp_path, monkeypatch, rows, **extra):
    table = tmp_path / 'annotated.tsv'
    pd.DataFrame(rows).to_csv(table, sep='\t', index=False)
    hits, stats = tmp_path / 'hits.tsv', tmp_path / 'stats.tsv'
    argv = ['14', '--input', str(table), '--output-hits', str(hits), '--output-stats', str(stats),
            '--library', 'yLtest', '--comparison', 'cmpA']
    for flag, value in extra.items():
        argv += [f'--{flag.replace("_", "-")}', str(value)]
    monkeypatch.setattr(sys, 'argv', argv)
    assert _step14().main() == 0
    return pd.read_csv(hits, sep='\t')


def _null_rows(lfc_values, **extra):
    """Synonymous rows that define the tier's null."""
    return [dict(bc1=f'S{i:03d}', log2FoldChange=value, pvalue=0.9, padj=0.9, baseMean=100,
                 tier='abundance_high', codon_variant_type='synonymous_codon_variant', **extra)
            for i, value in enumerate(lfc_values)]


def test_null_without_spread_calls_no_hits(tmp_path, monkeypatch):
    """25 identical synonymous LFCs give MAD 0; a 0.001 effect must not become an infinite z."""
    rows = _null_rows([0.0] * 25)
    rows.append(dict(bc1='T' * 20, log2FoldChange=0.001, pvalue=1e-9, padj=1e-9, baseMean=100,
                     tier='abundance_high', codon_variant_type='missense_codon_variant'))
    result = _hits(tmp_path, monkeypatch, rows).set_index('bc1')
    target = result.loc['T' * 20]
    assert not bool(target['is_significant'])
    assert not np.isfinite(pd.to_numeric(pd.Series([target['z_score']]), errors='coerce').iloc[0]) or \
        pd.isna(target['z_score'])


def test_concordance_gate_filters_where_the_column_exists(tmp_path, monkeypatch):
    rows = _null_rows(np.linspace(-1.0, 1.0, 25), concordance_fraction=0.9)
    rows.append(dict(bc1='L' * 20, log2FoldChange=8.0, pvalue=1e-12, padj=1e-12, baseMean=100,
                     tier='abundance_high', codon_variant_type='missense_codon_variant',
                     concordance_fraction=0.0))
    rows.append(dict(bc1='H' * 20, log2FoldChange=8.0, pvalue=1e-12, padj=1e-12, baseMean=100,
                     tier='abundance_high', codon_variant_type='missense_codon_variant',
                     concordance_fraction=0.9))
    result = _hits(tmp_path, monkeypatch, rows).set_index('bc1')
    assert bool(result.loc['H' * 20, 'is_significant'])
    assert not bool(result.loc['L' * 20, 'is_significant'])


def test_hits_are_still_called_without_a_concordance_column(tmp_path, monkeypatch):
    """Control: with no concordance data every BC1 passes that gate (and the run says so)."""
    rows = _null_rows(np.linspace(-1.0, 1.0, 25))
    rows.append(dict(bc1='H' * 20, log2FoldChange=8.0, pvalue=1e-12, padj=1e-12, baseMean=100,
                     tier='abundance_high', codon_variant_type='missense_codon_variant'))
    result = _hits(tmp_path, monkeypatch, rows).set_index('bc1')
    assert bool(result.loc['H' * 20, 'is_significant'])


def _nuclease_rows(designs):
    return pd.DataFrame([
        dict(sublibrary_variant_combo=name, comparison='cmpA', mean_log2fc=1.0,
             concordance_fraction=1.0, within_variant_cv=0.1, is_significant=sig,
             z_score=4.0, n_bc1=5, tier=tier)
        for name, sig, tier in designs])


def test_library_concordance_counts_designs_not_rows():
    """One design in two tiers, both significant: concordance is 1.0, not 2.0."""
    hits = _nuclease_rows([('SpG_NGG_varA', True, 'abundance_high'),
                           ('SpG_NGG_varA', True, 'abundance_medium')])
    agg = aggregate_across_nuclease_libraries(hits)
    assert agg['n_nuclease_libraries'].tolist() == [1]
    assert agg['library_concordance'].tolist() == [1.0]


def test_ngnh_exception_needs_ngnh_to_be_the_failing_library():
    only_ngnh = aggregate_across_nuclease_libraries(_nuclease_rows([
        ('SpG_NGNH_varA', True, 'abundance_high'), ('SpG_NGG_varA', False, 'abundance_high')]))
    assert not bool(only_ngnh['ngnh_low_efficiency_exception'].iloc[0])
    ngnh_fails = aggregate_across_nuclease_libraries(_nuclease_rows([
        ('SpG_NGNH_varA', False, 'abundance_high'), ('SpG_NGG_varA', True, 'abundance_high')]))
    assert bool(ngnh_fails['ngnh_low_efficiency_exception'].iloc[0])


def test_fractional_generation_is_not_the_deseq2_reference():
    matrix = pd.DataFrame({'ref': [10], 'half': [10], 'late': [10]}, index=['bc1_1'])
    metadata = pd.DataFrame([
        dict(sample_name='ref', condition='YPD', generations=0),
        dict(sample_name='half', condition='YPD', generations=0.5),
        dict(sample_name='late', condition='YPD', generations=6),
    ])
    comparisons = split_by_comparison(matrix, metadata)
    assert all(c['reference_samples'] == ['ref'] for c in comparisons.values())
    assert 'YPD_0.5g' in comparisons and 'YPD_6g' in comparisons
