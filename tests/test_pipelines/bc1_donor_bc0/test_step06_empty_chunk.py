"""An empty step-06 chunk merges with the others (218 patch 50, register BDB-H24).

split_bc1_reference_table writes a header-only chunk when there are fewer unmatched donors
than chunks. process_chunk then wrote pd.DataFrame(), a file with no header line, and
main's merge failed reading it (pandas EmptyDataError), so one empty chunk failed the run.
"""
import importlib.util
import sys
from pathlib import Path

import pandas as pd

from magestic.pipelines.bc1_donor_bc0.readers import data_loader


def _step06():
    path = Path(data_loader.__file__).parents[1] / 'snakemake/06_annotate_unmatched_donors.py'
    spec = importlib.util.spec_from_file_location('test_step06_empty_chunk', path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def run_chunk(step06, tmp_path, rows, index):
    chunk = tmp_path / f'chunk_{index}.tsv'
    pd.DataFrame(rows, columns=['bc1', 'donor', 'yeast_sublibrary', 'V_library']).to_csv(chunk, sep='\t', index=False)
    info = step06.OligoInfo('oligo_1', 'G' * 20, 'A' * 129, 'A' * 60, '')
    lookup = {info.oligo_name: info}
    index_ = step06.build_kmer_index(lookup)
    # No reference FASTA, so the chunk is matched against the designs but not aligned.
    step06.init_worker(lookup, lookup, index_, index_, {}, tmp_path / 'absent.fa', {},
                       tmp_path / 'output', tmp_path / 'progress')
    _, result = step06.process_chunk(chunk, index)
    return pd.read_csv(result, sep='\t')


def test_an_empty_chunk_merges_with_a_full_one(tmp_path):
    step06 = _step06()
    (tmp_path / 'progress').mkdir()
    full = run_chunk(step06, tmp_path, [dict(bc1='A' * 20, donor='A' * 129, yeast_sublibrary='SpG_NGNG',
                                             V_library='V629')], 0)
    empty = run_chunk(step06, tmp_path, [], 1)
    assert empty.empty and list(empty.columns) == list(full.columns)
    assert len(pd.concat([full, empty], ignore_index=True)) == 1
