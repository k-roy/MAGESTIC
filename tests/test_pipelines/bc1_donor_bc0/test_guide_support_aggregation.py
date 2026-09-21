"""Guide bridge support is summed over sequence rows within a source library."""
import importlib.util
from pathlib import Path
import sys

import magestic
import pandas as pd
import pytest


def loader():
    path = Path(magestic.__file__).parent/'pipelines/bc1_donor_bc0/snakemake/03_link_guide_donor_bc0.py'
    spec = importlib.util.spec_from_file_location('guide_support_test', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def guide_row(guide, count, library='V_test', sps='A'*20, **extra):
    return dict(library_ID=library, guide=guide, total_counts=count, SPS=sps,
                donor_bc0_fragment='C'*30+sps+'G'*10, **extra)


def run_loader(tmp_path, rows):
    source = tmp_path/'V_fixture_guide_donor_bc0_counts.tsv'
    pd.DataFrame(rows).to_csv(source, sep='\t', index=False)
    output, _ = loader().load_guide_table(tmp_path, '*_guide_donor_bc0_counts.tsv',
                                         '*_guide_donor_bc0_purity_filtered.tsv', str(tmp_path/'linked.tsv'))
    return output


@pytest.mark.parametrize('count_column', ['total_counts', 'counts'])
def test_actual_cli_selects_aggregated_guide_support(tmp_path, monkeypatch, count_column):
    guide_dir = tmp_path/'gdb'
    guide_dir.mkdir()
    guides = pd.DataFrame([guide_row('A'*20, 60, leader='first'),
                           guide_row('A'*20, 60, leader='second'),
                           guide_row('T'*20, 100, leader='other')])
    guides.rename(columns={'total_counts': count_column}).to_csv(
        guide_dir/'V_test_guide_donor_bc0_counts.tsv', sep='\t', index=False)
    source, dest = tmp_path/'bc1.tsv', tmp_path/'linked.tsv'
    pd.DataFrame([dict(bc1='A'*20, donor='C'*129, bc0='G'*20, total_counts=13)]).to_csv(source, sep='\t', index=False)
    monkeypatch.setattr(sys, 'argv', ['link', '--input', str(source), '--output', str(dest), '--guide-dir', str(guide_dir)])
    loader().main()
    result = pd.read_csv(dest, sep='\t')
    assert len(result) == 1
    row = result.iloc[0]
    assert row.guide == 'A'*20
    assert row.guide_link_counts == 120
    assert row.guide_link_n_guides == 2  # Winning support is not unique evidence.
    assert row.total_counts == 13  # Barcode-assay counts are not guide-assay reads.


def test_representative_metadata_comes_from_real_best_row(tmp_path):
    result = run_loader(tmp_path, [guide_row('A'*20, 7, sps='A'*20),
                                  guide_row('A'*20, 5, sps='T'*20),
                                  guide_row('C'*20, 10, sps='G'*20)])
    row = result.iloc[0]
    assert row.guide == 'A'*20
    assert row.guide_link_counts == 12
    assert row.SPS == 'A'*20
    assert row.library_ID == 'V_test'


def test_named_source_libraries_are_not_implicitly_pooled(tmp_path):
    result = run_loader(tmp_path, [guide_row('A'*20, 60, library='V_first'),
                                  guide_row('A'*20, 60, library='V_second'),
                                  guide_row('C'*20, 100, library='V_third')])
    row = result.iloc[0]
    assert row.guide == 'C'*20
    assert row.guide_link_counts == 100
    assert row.library_ID == 'V_third'
    assert row.guide_link_n_guides == 2
