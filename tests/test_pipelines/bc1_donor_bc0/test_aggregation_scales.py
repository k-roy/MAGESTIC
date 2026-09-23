"""BDB aggregation works on records, not on one frame per group (218 patch 48).

On the 20250912 screen (3.1 million 2x150 rows), aggregate_link_evidence ran for hours: it converted
every link group, and every step of the greedy sequence collapse, into its own frame and records.
Patch 48 converts each frame once and groups row positions (in the order frame.groupby(..., sort=False)
visits them), so the output is unchanged. The reference below is the pre-patch group-by-group code,
kept verbatim except that it calls the module's _combine wrapper.
"""
import random

import pandas as pd
import pytest
from rapidfuzz.distance.Levenshtein import distance

from magestic.pipelines.bc1_donor_bc0.core import aggregation
from magestic.pipelines.bc1_donor_bc0.core.aggregation import aggregate_link_evidence


def reference_combine_groups(frame, group_columns, count_columns):
    return [aggregation._combine(group, count_columns)
            for _, group in frame.groupby(group_columns, dropna=False, sort=False)]


def reference_collapse(frame, group_columns, sequence_column, count_columns, threshold):
    rows = []
    for _, group in frame.groupby(group_columns, dropna=False, sort=False):
        group = group.reset_index(drop=True)
        while len(group):
            top_index = group['counts'].idxmax()
            representative = group.loc[top_index].to_dict()
            sequence = aggregation._text(representative[sequence_column])
            if sequence is None:
                mask = group.index == top_index
            else:
                mask = group[sequence_column].map(
                    lambda value: aggregation._text(value) is not None and distance(str(value), sequence) <= threshold)
            rows.append(aggregation._combine(group.loc[mask], count_columns, representative))
            group = group.loc[~mask].reset_index(drop=True)
    return pd.DataFrame(rows, columns=frame.columns.union(pd.Index(['guide_candidates', 'guide_evidence_conflict']), sort=False))


def mutate(sequence, rng, edits):
    letters = list(sequence)
    for _ in range(edits):
        position = rng.randrange(len(letters))
        letters[position] = rng.choice('ACGT'.replace(letters[position], ''))
    return ''.join(letters)


def observations(seed, n=400):
    """Raw link evidence with ties, repeated links, unknown components, scopes and near-duplicates."""
    rng = random.Random(seed)
    bc1s = [''.join(rng.choice('ACGT') for _ in range(20)) for _ in range(10)]
    bc1s += [mutate(bc1, rng, 1) for bc1 in bc1s[:5]]
    donors = [''.join(rng.choice('ACGT') for _ in range(130)) for _ in range(5)]
    donors += [mutate(donor, rng, 3) for donor in donors[:3]]
    rows = []
    for _ in range(n):
        donor = rng.choice(donors + [None])
        bc0 = rng.choice(['A' * 10, 'C' * 10, 'G' * 10])
        rows.append(dict(
            bc1=rng.choice(bc1s + [None]), donor=donor, bc0=bc0,
            donor_bc0_fragment=None if donor else rng.choice(['T' * 60, 'A' * 60]),
            counts=rng.choice([0, 1, 2, 5, 5, 8, 13, 40]),
            pcr_replicate_id=f'p{rng.randrange(6)}',
            trafo_ID=rng.choice(['trafo_1', 'trafo_2', 'empty', None]),
            SPS=rng.choice(['CGATAGACGACTGGACAGCA', None]),
            guide=rng.choice([None, None, 'G' * 20, 'C' * 20])))
    return pd.DataFrame(rows)


@pytest.mark.parametrize('collapse', [False, True])
@pytest.mark.parametrize('seed', [1, 2, 3])
def test_matches_the_group_by_group_reference(monkeypatch, seed, collapse):
    frame = observations(seed)
    fast = aggregate_link_evidence(frame, min_reads=1, min_pcr_replicates=2, collapse_sequences=collapse,
                                   show_progress=False)
    monkeypatch.setattr(aggregation, '_combine_groups', reference_combine_groups)
    monkeypatch.setattr(aggregation, '_collapse', reference_collapse)
    reference = aggregate_link_evidence(frame, min_reads=1, min_pcr_replicates=2, collapse_sequences=collapse,
                                        show_progress=False)
    assert len(fast) > 20
    pd.testing.assert_frame_equal(fast, reference)


def test_grouping_converts_each_frame_once_not_once_per_group(monkeypatch):
    """2,000 links of one plasmid in 2,000 lineages: one bc1-collapse group with 2,000 rows."""
    rng = random.Random(7)
    frame = pd.DataFrame([dict(bc1=''.join(rng.choice('ACGT') for _ in range(20)), donor='T' * 130,
                               bc0='A' * 10, counts=rng.randrange(1, 50), pcr_replicate_id=f'p{i % 3}',
                               trafo_ID='trafo_1') for i in range(2000)])
    calls = []
    original = pd.DataFrame.to_dict

    def counting(self, *args, **kwargs):
        calls.append(len(self))
        return original(self, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, 'to_dict', counting)
    result = aggregate_link_evidence(frame, min_reads=1, min_pcr_replicates=2, collapse_sequences=True,
                                     show_progress=False)
    assert len(result) > 1900
    assert len(calls) <= 12, f'{len(calls)} frame-to-records conversions for 2,000 links'
