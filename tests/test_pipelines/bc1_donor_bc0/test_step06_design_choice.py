"""Step 06 takes the nearer of the two design tables (218 patch 52, register BDB-H12).

A row labelled trafo_3 or trafo_4 (both V634, the 2021 NGG pool) was sent to the 2024 SpG table,
because any unrecognized label defaulted to 'spg', and a row with no label became 'SpG_NGNG'. On
the real references each pool's unmatched donors match their own table at a median edit distance
of 1 and the other at 2-14, and a third of the matched V634 rows of 20250912 carry names found only
in the 2024 table, so both tables are searched, the nearer design wins, and the row's own pool (its
SPS through sps_library_mapping.tsv, else its V library) breaks an exact tie.
"""
import importlib.util
import random
import sys
from pathlib import Path

import pandas as pd

from magestic.pipelines.bc1_donor_bc0.readers import data_loader

rng = random.Random(52)
NGG = ''.join(rng.choice('ACGT') for _ in range(129))
SPS_2021, SPS_2024 = 'CGATAGACGACTGGACAGCA', 'CTTGACTCGACAGTGAGAGC'


def mutate(sequence, positions):
    letters = list(sequence)
    for position in positions:
        letters[position] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[letters[position]]
    return ''.join(letters)


def _step06():
    path = Path(data_loader.__file__).parents[1] / 'snakemake/06_annotate_unmatched_donors.py'
    spec = importlib.util.spec_from_file_location('test_step06_design_choice', path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def annotate(tmp_path, rows, spg_donor, spcas9_donor=NGG, sps_pools=None):
    step06 = _step06()
    chunk = tmp_path / 'chunk.tsv'
    pd.DataFrame(rows).to_csv(chunk, sep='\t', index=False)
    spg = {'SpG_oligo': step06.OligoInfo('SpG_oligo', 'G' * 20, spg_donor, spg_donor[30:90], '')}
    spcas9 = {'SpCas9_oligo': step06.OligoInfo('SpCas9_oligo', 'C' * 20, spcas9_donor, spcas9_donor[30:90], '')}
    # The SPS pools are passed only when a test has them, so the other tests also run at a parent
    # commit whose init_worker has no such parameter and fail there on their routing, not the call.
    step06.init_worker(spg, spcas9, step06.build_kmer_index(spg), step06.build_kmer_index(spcas9), {},
                       tmp_path / 'absent.fa', {}, tmp_path / 'output', tmp_path / 'progress',
                       **({'sps_pools': sps_pools} if sps_pools else {}))
    _, result = step06.process_chunk(chunk, 0)
    return pd.read_csv(result, sep='\t')


def row(label=None, v_library=None, sps=None, donor=None):
    return dict(bc1='A' * 20, donor=donor or mutate(NGG, [50]), yeast_sublibrary=label, V_library=v_library, SPS=sps)


def test_a_trafo_labelled_v634_donor_takes_the_nearer_2021_design(tmp_path):
    """The 2024 table holds a design 7 edits away, as the real SpG table did for V634 donors."""
    result = annotate(tmp_path, [row('trafo_4', 'V634', SPS_2021)], spg_donor=mutate(NGG, range(5, 110, 16)))
    assert result.loc[0, 'closest_oligo_name'] == 'SpCas9_oligo'
    assert (result.loc[0, 'design_table'], result.loc[0, 'design_choice']) == ('spcas9', 'lower_edit_distance')
    assert result.loc[0, 'closest_donor_edit_distance'] == 1


def test_an_unlabelled_row_is_not_sent_to_the_spg_table(tmp_path):
    result = annotate(tmp_path, [row()], spg_donor=mutate(NGG, range(5, 110, 16)))
    assert result.loc[0, 'closest_oligo_name'] == 'SpCas9_oligo'
    assert pd.isna(result.loc[0, 'yeast_sublibrary'])  # no SpG_NGNG is filled in


def test_an_exact_tie_follows_the_rows_sps_pool(tmp_path):
    """The same donor designed in both pools: the SPS mapping names the 2024 pool."""
    result = annotate(tmp_path, [row(sps=SPS_2024)], spg_donor=NGG, sps_pools={SPS_2024: 'spg'})
    assert (result.loc[0, 'design_table'], result.loc[0, 'design_choice']) == ('spg', 'tie_pool')
    assert result.loc[0, 'closest_oligo_name'] == 'SpG_oligo'


def test_an_exact_tie_follows_the_v_library_without_an_sps_mapping(tmp_path):
    result = annotate(tmp_path, [row(v_library='V634')], spg_donor=NGG)
    assert (result.loc[0, 'design_table'], result.loc[0, 'design_choice']) == ('spcas9', 'tie_pool')


def test_a_tie_nothing_breaks_chooses_no_design(tmp_path):
    result = annotate(tmp_path, [row(label='trafo_4')], spg_donor=NGG)
    assert pd.isna(result.loc[0, 'closest_oligo_name']) and result.loc[0, 'design_choice'] == 'ambiguous'
    assert result.loc[0, 'confidence_tier'] != 'tier_1_low_ed_attributable'
    assert set(result.loc[0, 'closest_candidates'].split(';')) == {'spcas9:SpCas9_oligo:1', 'spg:SpG_oligo:1'}


def test_a_donor_found_in_one_table_only_takes_it(tmp_path):
    other = ''.join(random.Random(7).choice('ACGT') for _ in range(129))
    result = annotate(tmp_path, [row(label='SpG_NGNG')], spg_donor=other)
    assert (result.loc[0, 'design_table'], result.loc[0, 'design_choice']) == ('spcas9', 'only_spcas9')
