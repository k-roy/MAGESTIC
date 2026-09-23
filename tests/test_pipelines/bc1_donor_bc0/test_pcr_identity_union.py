"""One physical PCR counts once, however many routes sequenced it (218 regression 1).

The same PCR sequenced 2x300 and 2x150 carries two identity strings when the two routes'
keyfiles declare different fields: the yL437/yL442 2x150 combined key declares no gDNA plate,
inner primers or polymerase, and the plate keys read by the 2x300 route do. Support is the
fewest physical PCRs consistent with the declarations, and an inferred merge is inexact.

Design: planning/218_magestic_production_corrections/reviews/merge_regressions_repair/DESIGN.md §2.
Tests that call ``data_loader.count_distinct_pcrs`` look it up inside the test, so at a parent
without the helper they fail on the missing name rather than on behaviour; every other test
goes through the merge and aggregation entry points that exist at the parent.
"""
import itertools
import random

import pandas as pd
import pytest

from magestic.pipelines.bc1_donor_bc0.core.merge import merge_data_sources
from magestic.pipelines.bc1_donor_bc0.pipeline.run_pipeline import aggregate_and_filter
from magestic.pipelines.bc1_donor_bc0.readers import data_loader
from magestic.pipelines.bc1_donor_bc0.readers.data_loader import (
    PCR_REPLICATE_ID_FIELDS, build_pcr_replicate_id, load_2x150_bc1_donor_bc0_data,
    load_2x300_data)

SPS = 'CGATAGACGACTGGACAGCA'
BC1, BC0 = 'G' * 20, 'A' * 10
DONOR = 'T' * 100 + 'ACGTACGTACGTACGTACGTACGTACGTAC'
FRAGMENT = DONOR[-30:] + SPS + BC0
DECLARED = dict(gDNA_plate_name='prep_A', inner_primers='inner_A', plate_polymerase='Q5')
WELL = dict(outer_primers='plate_1A', condition='plate wash', trafo_ID='trafo_3', sample_number=1)


def identity(**fields):
    return build_pcr_replicate_id({**WELL, **fields})


FULL_ID = identity(**DECLARED)
SHORT_ID = identity()  # what the yL437/yL442-shaped combined key declares for the same well


def aggregate(frame, collapse):
    return aggregate_and_filter(frame, min_reads=1, min_pcr_replicates=2,
                                collapse_sequences=collapse, show_progress=False)


def _one_pcr_on_two_routes(tmp_path, combined_key_extra):
    """One PCR (prep_A/inner_A/plate_1A/Q5, well 1, trafo_3): 2x300 on 20250101, 2x150 on 20250102."""
    long_data, long_keys = tmp_path / 'long_data', tmp_path / 'long_keys'
    short_data, short_keys = tmp_path / 'short_data', tmp_path / 'short_keys'
    for folder in (long_data, long_keys, short_data, short_keys):
        folder.mkdir()
    pd.DataFrame([dict(description='bc1-donor-bc0', sequencing_date='20250101',
                       outer_primers='plate_1A', inner_primers='inner_A',
                       gDNA_plate_name='prep_A', polymerase='Q5')]).to_csv(
        long_keys / 'plate_key.tsv', sep='\t', index=False)
    pd.DataFrame([dict(sample_number=1, trafo_ID='trafo_3', condition='plate wash')]).to_csv(
        long_keys / 'sample_key.tsv', sep='\t', index=False)
    pd.DataFrame([dict(counts=7, bc1=BC1, donor=DONOR, bc0=BC0, SPS=SPS,
                       donor_bc0_fragment=FRAGMENT)]).to_csv(
        long_data / '20250101_plate_1A_sample_1_bc1_donor_bc0.tsv', sep='\t', index=False)
    pd.DataFrame([dict(description='bc1-donor-bc0', sequencing_date='20250102',
                       outer_primers='plate_1A', sample_number=1, condition='plate wash',
                       trafo_ID='trafo_3', **combined_key_extra)]).to_csv(
        short_keys / 'bc1_donor_bc0_combined_key.tsv', sep='\t', index=False)
    pd.DataFrame([dict(counts=5, bc1=BC1, bc0=BC0, SPS=SPS, donor_bc0_fragment=FRAGMENT)]).to_csv(
        short_data / '20250102_plate_1A_1_bc1_donor_bc0_fragment_counts.tsv', sep='\t', index=False)
    df_2x300, _ = load_2x300_data(long_data, long_keys)
    df_2x150, _ = load_2x150_bc1_donor_bc0_data(short_data, short_keys)
    assert len(df_2x300) == len(df_2x150) == 1
    return merge_data_sources(df_2x150=df_2x150, df_2x300=df_2x300)


@pytest.mark.parametrize('collapse', [False, True])
def test_resequenced_pcr_counts_once_through_merge_and_aggregation(tmp_path, collapse):
    merged = _one_pcr_on_two_routes(tmp_path, {})
    assert len(merged) == 1
    row = merged.iloc[0]
    assert set(row['pcr_replicate_ids'].split(';')) == {FULL_ID, SHORT_ID}
    assert row['total_pcr_replicates'] == 1
    assert row['pcr_replicates_2x300'] == row['pcr_replicates_2x150'] == 1
    assert not row['pcr_support_exact']  # the 2x150 PCR could still be a different one
    result = aggregate(merged, collapse)
    assert len(result) == 1
    row = result.iloc[0]
    assert row['num_PCR_replicates'] == row['total_pcr_replicates'] == 1
    assert row['support_tier'] == 'provisional'
    assert not row['pcr_support_exact']


def test_fully_declared_routes_give_one_exact_pcr(tmp_path):
    merged = _one_pcr_on_two_routes(tmp_path, dict(gDNA_plate_name='prep_A',
                                                   inner_primers='inner_A', polymerase='Q5'))
    row = merged.iloc[0]
    assert row['pcr_replicate_ids'] == FULL_ID
    assert row['total_pcr_replicates'] == 1 and row['pcr_support_exact']
    row = aggregate(merged, True).iloc[0]
    assert row['num_PCR_replicates'] == 1 and row['pcr_support_exact']


def _full(count, pcr, donor=DONOR):
    return dict(counts=count, bc1=BC1, bc0=BC0, donor=donor, SPS=SPS,
                donor_bc0_fragment=donor[-30:] + SPS + BC0, pcr_replicate_id=pcr)


def _short(count, pcr, recovered=None):
    return dict(counts=count, bc1=BC1, bc0=BC0, SPS=SPS, donor_bc0_fragment=FRAGMENT,
                donor=recovered, pcr_replicate_id=pcr)


def test_recovered_donor_route_counts_the_pcr_once():
    """The second enabler: a recovered 2x150 donor shares the 2x300 link key directly."""
    merged = merge_data_sources(df_2x150=pd.DataFrame([_short(5, SHORT_ID, recovered=DONOR)]),
                                df_2x300=pd.DataFrame([_full(7, FULL_ID)]))
    assert len(merged) == 1
    row = merged.iloc[0]
    assert row['total_counts'] == 12 and row['total_pcr_replicates'] == 1
    assert not row['pcr_support_exact']


def test_two_preps_of_one_well_stay_two_pcrs():
    """prep_A and prep_B differ on a field both declare; the 2x150 PCR may be either one."""
    prep_b = identity(**dict(DECLARED, gDNA_plate_name='prep_B'))
    merged = merge_data_sources(df_2x150=pd.DataFrame([_short(5, SHORT_ID)]),
                                df_2x300=pd.DataFrame([_full(7, FULL_ID), _full(4, prep_b)]))
    assert len(merged) == 1
    row = merged.iloc[0]
    assert row['total_counts'] == 16 and row['total_pcr_replicates'] == 2
    assert not row['pcr_support_exact']


def test_aggregation_collapse_counts_a_shared_pcr_once():
    """A donor error variant from the other route joins at collapse; its PCR is the same one."""
    variant = DONOR[:50] + 'G' + DONOR[51:]
    merged = merge_data_sources(df_2x150=pd.DataFrame([_short(3, SHORT_ID, recovered=variant)]),
                                df_2x300=pd.DataFrame([_full(20, FULL_ID)]))
    assert len(merged) == 2
    row = aggregate(merged, True).iloc[0]
    assert row['counts'] == 23 and row['num_PCR_replicates'] == 1
    assert row['pcr_replicates_2x300'] == row['pcr_replicates_2x150'] == 1
    assert not row['pcr_support_exact']


def test_an_unknown_amplicon_level_is_undeclared():
    count = data_loader.count_distinct_pcrs
    no_trafo = build_pcr_replicate_id({k: v for k, v in {**WELL, **DECLARED}.items() if k != 'trafo_ID'})
    trafo, pooled = identity(**DECLARED), identity(**dict(DECLARED, trafo_ID='yL387'))
    assert count({no_trafo, trafo}) == (1, False)
    assert count({no_trafo, pooled}) == (1, False)
    assert count({trafo, pooled}) == (2, True)  # trafo-level and pooled_yL are both declared


def test_other_identity_encodings_stay_opaque():
    count = data_loader.count_distinct_pcrs
    legacy = 'prep_A_inner_A_plate_1A_Q5_plate wash_trafo_1'
    assert count({'p1', 'p2', legacy}) == (3, True)
    assert count({'p1', FULL_ID, SHORT_ID}) == (2, False)
    assert count(set()) == (0, True) and count({FULL_ID}) == (1, True)


def _trap():
    """First-fit in this order gives 3: u1, u2, v1 open two groups, v1 joins u1, v2 fits neither.

    The minimum is 2: {u1, v2} and {u2, v1}.
    """
    return [identity(gDNA_plate_name='A', plate_polymerase='P'),                 # u1
            identity(gDNA_plate_name='B', plate_polymerase='P'),                 # u2
            identity(inner_primers='X', plate_polymerase='P'),                   # v1
            identity(gDNA_plate_name='A', inner_primers='Y')]                    # v2


def test_minimum_is_found_where_first_fit_overcounts():
    for order in itertools.permutations(_trap()):
        assert data_loader.count_distinct_pcrs(order) == (2, False)


def test_exhausted_search_reports_a_proven_lower_bound(monkeypatch):
    monkeypatch.setattr(data_loader, 'PCR_COVER_SEARCH_LIMIT', 1)
    count, exact = data_loader.count_distinct_pcrs(_trap())
    assert 1 <= count <= 2 and not exact


def _minimum_cover_by_enumeration(declarations):
    """Every set partition of up to 7 declarations; the smallest whose blocks are consistent."""
    def partitions(items):
        if not items:
            yield []
            return
        first, rest = items[0], items[1:]
        for smaller in partitions(rest):
            for i in range(len(smaller)):
                yield smaller[:i] + [[first] + smaller[i]] + smaller[i + 1:]
            yield [[first]] + smaller

    def consistent(block):
        return all(len({d[k] for d in block if d[k] is not None}) <= 1 for k in range(len(block[0])))

    return min(len(p) for p in partitions(list(declarations)) if all(consistent(b) for b in p))


def test_count_is_the_minimum_and_never_exceeds_the_planted_pcrs():
    rng = random.Random(218)
    alphabet = dict(gDNA_plate_name='ab', inner_primers='xy', outer_primers='12',
                    plate_polymerase='QK', condition='cd', trafo_ID=('trafo_1', 'yL1'),
                    sample_number='34')
    for _ in range(300):
        planted = [{k: rng.choice(v) for k, v in alphabet.items()} for _ in range(rng.randint(1, 4))]
        observed, used = [], set()
        for _ in range(rng.randint(1, 7)):
            pcr = rng.randrange(len(planted))
            used.add(pcr)
            observed.append(build_pcr_replicate_id(
                {k: (None if rng.random() < 0.35 else v) for k, v in planted[pcr].items()}))
        declared = sorted({data_loader._declared_pcr_fields(i) for i in observed}, key=str)
        count, exact = data_loader.count_distinct_pcrs(observed)
        assert count == _minimum_cover_by_enumeration(declared)
        assert count <= len({tuple(sorted(planted[p].items())) for p in used})
        compatible_pair = any(all(x is None or y is None or x == y for x, y in zip(a, b))
                              for a, b in itertools.combinations(declared, 2))
        assert exact == (not compatible_pair)


def test_many_wells_on_two_routes_resolve_per_well():
    wells = range(1, 61)
    full = {identity(**dict(DECLARED, sample_number=w)) for w in wells}
    short = {identity(sample_number=w) for w in wells}
    assert data_loader.count_distinct_pcrs(full | short) == (60, False)
    assert data_loader.count_distinct_pcrs(full) == (60, True)
    assert len(PCR_REPLICATE_ID_FIELDS) == 7
