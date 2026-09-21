"""Count conservation across every GDB step02 collapsed-FASTQ reader."""
import importlib.util
from pathlib import Path
import sys

import magestic
import pandas as pd
import pytest

from magestic.pipelines.guide_donor_bc0.core import fastq_processing as core
from magestic.pipelines.guide_donor_bc0.core import sharded_combine as sharded
from magestic.utils.collapsed_fastq import CollapsedHeaderError


def load_script(name):
    path = Path(magestic.__file__).parent / 'pipelines/guide_donor_bc0/snakemake' / name
    spec = importlib.util.spec_from_file_location('count_test_' + path.stem, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(params=['core', 'numbered', 'optimized', 'sharded'])
def reader(request, tmp_path):
    def run(fastq):
        if request.param in ('core', 'numbered'):
            module = core if request.param == 'core' else load_script('02_combine_step_1_collapsed_fastq.py')
            counts = {}
            records, reads = module.process_collapsed_fastq(
                fastq, 'KR1952-KR1590', '20260920', 'V_test', counts, {'20260920'})
            return records, reads, sum(row['fwd_0920'] for row in counts.values())
        if request.param == 'optimized':
            module = load_script('02_combine_step_1_collapsed_fastq_optimized.py')
            frame = module.process_single_fastq((fastq, 'KR1952-KR1590', '20260920', 'V_test'))
            return len(frame), int(frame['counts'].sum()), int(frame['counts'].sum())
        row = dict(shard_id='1', fastq=str(fastq), library_ID='V_test', sample_number='sample_1',
                   sequencing_date='20260920', inner_primers='KR1952-KR1590', counts_col='fwd_0920', direction='fwd')
        result = sharded.parse_shard(row, tmp_path/'shards')
        frame = pd.read_csv(result['tsv'], sep='\t')
        return result['unique_reads_in'], result['reads_in'], int(frame['count'].sum())
    return run


def write_fastq(path, names):
    sequence = 'GACTAC' + 'ACGT'*5 + 'GTTTGAAGAGC' + 'A'*129 + 'C'*20 + 'G'*10 + core.SEQ_CONST.msd
    path.write_text(''.join(f'@{name}\n{sequence}\n+\n'+ 'J'*len(sequence)+'\n' for name in names))
    return path


def test_mixed_producer_headers_conserve_abundance(reader, tmp_path):
    # Both formats coexist in one file; per-file format detection cannot pass.
    path = write_fastq(tmp_path/'mixed.fastq', ['12', 'seq_1_count=208'])
    assert reader(path) == (2, 220, 220)


def test_legacy_headers_still_conserve_abundance(reader, tmp_path):
    assert reader(write_fastq(tmp_path/'legacy.fastq', ['12', '208'])) == (2, 220, 220)


@pytest.mark.parametrize('name', ['not_a_count', '0', 'seq_1_count=0'])
def test_unknown_or_zero_counts_fail_explicitly(reader, tmp_path, name):
    with pytest.raises(CollapsedHeaderError, match='read name'):
        reader(write_fastq(tmp_path/'invalid.fastq', [name]))
