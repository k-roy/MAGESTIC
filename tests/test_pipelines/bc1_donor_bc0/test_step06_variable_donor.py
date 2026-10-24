"""Step 06 compares the variable donor, not the raw read (218 patch 51, register BDB-H11-R1).

Unmatched donors keep the raw read representation: a 5' plasmid prefix and cloning motif ahead of
the 129-nt variable donor (median 238 nt on the real 20250912 and yL406 references). Compared
whole, no raw donor was within MAX_EDIT_DISTANCE of a design, so the March 2026 step-06 runs
classed every unmatched donor tier_5_no_match (208,649 of 208,649 and 190,953 of 190,953). The
last 129 nt matched at a median edit distance of 1 (probes/step06_ground_truth.py).
"""
import importlib.util
import random
import sys
from pathlib import Path

import pandas as pd
import pytest

from magestic.pipelines.bc1_donor_bc0.readers import data_loader

rng = random.Random(51)
DESIGNED = ''.join(rng.choice('ACGT') for _ in range(129))
PREFIX = ''.join(rng.choice('ACGT') for _ in range(98)) + 'GTTTGAAGAGC'


def mutate(sequence, positions):
    letters = list(sequence)
    for position in positions:
        letters[position] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[letters[position]]
    return ''.join(letters)


def _step06():
    path = Path(data_loader.__file__).parents[1] / 'snakemake/06_annotate_unmatched_donors.py'
    spec = importlib.util.spec_from_file_location('test_step06_variable_donor', path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def annotate(tmp_path, donors):
    step06 = _step06()
    chunk = tmp_path / 'chunk.tsv'
    pd.DataFrame([dict(bc1=f'{i:020d}', donor=donor, yeast_sublibrary='SpG_NGNG', V_library='V629')
                  for i, donor in enumerate(donors)]).to_csv(chunk, sep='\t', index=False)
    info = step06.OligoInfo('oligo_1', 'G' * 20, DESIGNED, DESIGNED[30:90], '')
    lookup = {info.oligo_name: info}
    index = step06.build_kmer_index(lookup)
    step06.init_worker(lookup, lookup, index, index, {}, tmp_path / 'absent.fa', {},
                       tmp_path / 'output', tmp_path / 'progress')
    _, result = step06.process_chunk(chunk, 0)
    return pd.read_csv(result, sep='\t')


def test_a_raw_prefixed_donor_finds_its_design(tmp_path):
    raw = PREFIX + mutate(DESIGNED, [40, 100])
    assert len(raw) == 238
    row = annotate(tmp_path, [raw]).iloc[0]
    assert (row['closest_oligo_name'], row['closest_donor_edit_distance']) == ('oligo_1', 2)
    assert row['observed_donor_length'] == 238  # the raw read length is still reported


def test_a_designed_length_donor_is_compared_as_it_is(tmp_path):
    row = annotate(tmp_path, [mutate(DESIGNED, [40])]).iloc[0]
    assert (row['closest_oligo_name'], row['closest_donor_edit_distance']) == ('oligo_1', 1)


@pytest.mark.parametrize('donor, expected', [
    (PREFIX + DESIGNED, DESIGNED), (DESIGNED, DESIGNED), (DESIGNED[:100], DESIGNED[:100]), (None, None)])
def test_variable_donor_keeps_the_last_designed_length(donor, expected):
    assert _step06().variable_donor(donor) == expected
