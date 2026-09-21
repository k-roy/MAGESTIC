"""Missing inputs must fail; valid empty inputs must replace stale count tables."""
import importlib
import json
from pathlib import Path

import pandas as pd
import pytest

from magestic.pipelines.guide_donor_bc0.config import GuideDonorBC0Config
from magestic.pipelines.guide_donor_bc0.core import fastq_processing as core

runner = importlib.import_module('magestic.pipelines.guide_donor_bc0.pipeline.run_pipeline')
DATE = '20200102'


def setup_case(tmp_path, rows=(('s1', 'V_test'),)):
    fastqs = tmp_path/'fastqs'
    fastqs.mkdir()
    key = tmp_path/'key.tsv'
    key.write_text('sample_number\tlibrary_ID\tinner_primers\touter_primers\tsequencing_dates\tpolymerase\n' +
                   ''.join(f'{sample}\t{library}\tKR1952-KR1590\tplate\t{DATE}\tQ5\n'
                           for sample, library in rows))
    config = GuideDonorBC0Config(project_name='fixture', project_dir=tmp_path,
                                collapsed_fastq_dir=fastqs, sample_keyfile=key,
                                output_dir=tmp_path/'output')
    return config, {sample: fastqs/f'{DATE}_plate_{sample}_collapsed.fastq'
                    for sample, _ in rows}


def write_fastq(path, count=7):
    seq = 'GCAGGCGCGCC' + 'G'*20 + 'GTTTGAAGAGC' + 'C'*129 + 'ACGT'*5 + 'AT'*5 + 'AGGAAAACAGAC'
    path.write_text(f'@{count}\n{seq}\n+\n' + 'J'*len(seq) + '\n')


def table_path(config, library='V_test'):
    return config.get_working_dir('guide_donor_bc0')/f'{library}_guide_donor_bc0_counts.tsv'


@pytest.mark.parametrize('partial', [False, True])
def test_missing_rerun_preserves_old_table_and_sidecar_but_fails(tmp_path, partial):
    rows = (('s1', 'V_test'), ('s2', 'V_test')) if partial else (('s1', 'V_test'),)
    config, fastqs = setup_case(tmp_path, rows)
    for path in fastqs.values():
        write_fastq(path)
    assert runner.run_combine_fastq(config) is True
    table = table_path(config)
    sidecar = Path(str(table)+'.provenance.json')
    before = (table.read_bytes(), sidecar.read_bytes())
    fastqs['s1'].unlink()
    assert runner.run_combine_fastq(config) is False
    assert (table.read_bytes(), sidecar.read_bytes()) == before


@pytest.mark.parametrize('prior_output', [False, True])
def test_valid_empty_fastq_writes_empty_table_and_truthful_sidecar(tmp_path, prior_output):
    config, fastqs = setup_case(tmp_path)
    if prior_output:
        write_fastq(fastqs['s1'])
        assert runner.run_combine_fastq(config) is True
    fastqs['s1'].write_text('')
    assert runner.run_combine_fastq(config) is True
    table = table_path(config)
    frame = pd.read_csv(table, sep='\t')
    assert frame.empty
    assert {'library_ID', 'total_counts', 'fwd_0102', 'rev_0102', 'guide', 'donor', 'bc0'} <= set(frame)
    sidecar = json.loads(Path(str(table)+'.provenance.json').read_text())
    assert sidecar['counts']['reads_in'] == sidecar['counts']['rows_out'] == 0
    assert sidecar['outputs'][0]['rows'] == 0
    assert len(sidecar['inputs']) == 2  # key and the actual empty FASTQ


@pytest.mark.parametrize('partial', [False, True])
def test_core_missing_input_refuses_before_replacing_output(tmp_path, partial):
    rows = (('s1', 'V_test'), ('s2', 'V_test')) if partial else (('s1', 'V_test'),)
    config, fastqs = setup_case(tmp_path, rows)
    if partial:
        write_fastq(fastqs['s2'])
    samples, dates = core.load_sample_keyfile(config.sample_keyfile)
    counts_dir = config.get_working_dir('guide_donor_bc0')
    counts_dir.mkdir(parents=True)
    table = table_path(config)
    table.write_text('prior output must remain untouched\n')
    with pytest.raises(FileNotFoundError, match='collapsed FASTQ'):
        core.process_library('V_test', samples['V_test'], config.collapsed_fastq_dir,
                             counts_dir, dates)
    assert table.read_text() == 'prior output must remain untouched\n'


def test_runner_preflights_all_selected_libraries_before_publishing(tmp_path):
    config, fastqs = setup_case(tmp_path, (('s1', 'V_first'), ('s2', 'V_missing')))
    write_fastq(fastqs['s1'])
    assert runner.run_combine_fastq(config) is False
    assert not table_path(config, 'V_first').exists()
    assert not table_path(config, 'V_missing').exists()


def test_empty_and_nonempty_schemas_agree_for_repeated_month_day(tmp_path):
    config, fastqs = setup_case(tmp_path)
    config.sample_keyfile.write_text(config.sample_keyfile.read_text().replace(DATE, DATE+',20210102'))
    second = config.collapsed_fastq_dir/'20210102_plate_s1_collapsed.fastq'
    fastqs['s1'].write_text('')
    second.write_text('')
    assert runner.run_combine_fastq(config) is True
    empty_header = table_path(config).read_text().splitlines()[0].split('\t')
    assert len(empty_header) == len(set(empty_header))
    write_fastq(fastqs['s1'])
    assert runner.run_combine_fastq(config) is True
    populated_header = table_path(config).read_text().splitlines()[0].split('\t')
    assert empty_header == populated_header
    # This only preserves the existing shared MMDD schema. Distinguishing years
    # requires a separately reviewed migration of all producers and consumers.
