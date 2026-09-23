"""A BDB link is the plasmid linkage; where it was sampled is evidence (218 regression 2).

Kevin, 2026-09-23 (ledger §7 item 4): separate PCRs count toward a BDB link, because the concern
is PCR recombination. So a trafo-level plate and the pooled-yL outgrowth of that transformation
share a link, and so do the transformations of one construct library. The sampled population
(yeast_library, yeast_sublibrary, trafo_ID) is kept per link in ``population_counts`` and no
longer splits a link, its PCR support or its barcode's purity denominator. The construct scope
(library_ID, V_library, SPS) still separates links.

Design: planning/218_magestic_production_corrections/reviews/merge_regressions_repair/DESIGN.md §1.
Each test asserts the link grain first, so at a parent that splits by population it fails on
that behaviour before it reaches the new ``population_counts`` column.
"""
import json

import pandas as pd
import pytest

from magestic.pipelines.bc1_donor_bc0.core.merge import merge_data_sources
from magestic.pipelines.bc1_donor_bc0.pipeline.run_pipeline import aggregate_and_filter
from magestic.pipelines.bc1_donor_bc0.readers.data_loader import load_2x150_bc1_donor_bc0_data

SPS = 'CGTATGGAGATGACGTCAGG'
BC1, BC0 = 'G' * 20, 'A' * 10
DONOR = 'T' * 100 + 'ACGTACGTACGTACGTACGTACGTACGTAC'
FRAGMENT = DONOR[-30:] + SPS + BC0


def aggregate(frame, collapse=False):
    return aggregate_and_filter(frame, min_reads=1, min_pcr_replicates=2,
                                collapse_sequences=collapse, show_progress=False)


def populations(row):
    return {(item['trafo_ID'], item['counts']) for item in json.loads(row['population_counts'])}


def full(count, pcr, donor=DONOR, **scope):
    return dict(counts=count, bc1=BC1, bc0=BC0, donor=donor, SPS=SPS, pcr_replicate_id=pcr, **scope)


def test_trafo_plate_and_pooled_outgrowth_are_one_link(tmp_path):
    """yL406: plate trafo_2-1 is one of the three plates whose washes were pooled into yL387."""
    data, keys = tmp_path / 'data', tmp_path / 'keys'
    data.mkdir()
    keys.mkdir()
    pd.DataFrame([dict(description='bc1-donor-bc0', sequencing_date='20250101', outer_primers='plate_3F',
                       sample_number=well, condition='plate wash', trafo_ID=trafo,
                       gDNA_plate_name='prep_A', inner_primers='KR1967-KR1590', polymerase='Q5')
                  for well, trafo in ((2, 'trafo_2-1'), (14, 'yL387'))]).to_csv(
        keys / 'bc1_donor_bc0_combined_key.tsv', sep='\t', index=False)
    for well, count in ((2, 7), (14, 5)):
        pd.DataFrame([dict(counts=count, bc1=BC1, bc0=BC0, SPS=SPS, donor_bc0_fragment=FRAGMENT)]).to_csv(
            data / f'20250101_plate_3F_{well}_bc1_donor_bc0_fragment_counts.tsv', sep='\t', index=False)
    frame, _ = load_2x150_bc1_donor_bc0_data(data, keys)
    merged = merge_data_sources(df_2x150=frame)
    assert len(merged) == 1
    row = merged.iloc[0]
    assert row['total_counts'] == 12 and row['total_pcr_replicates'] == 2
    assert pd.isna(row['trafo_ID'])
    assert populations(row) == {('trafo_2-1', 7), ('yL387', 5)}
    row = aggregate(merged).iloc[0]
    assert row['num_PCR_replicates'] == 2 and row['support_tier'] == 'confident'
    assert row['bc1_purity'] == 1
    assert populations(row) == {('trafo_2-1', 7), ('yL387', 5)}


def test_transformations_of_one_construct_library_are_one_link():
    """20250912: trafo_3 (SpG) and trafo_4 (SpCas9) were transformed with the same cV586."""
    merged = merge_data_sources(df_2x300=pd.DataFrame([
        full(10, 'p1', trafo_ID='trafo_3'), full(20, 'p2', trafo_ID='trafo_4')]))
    assert len(merged) == 1
    row = merged.iloc[0]
    assert row['total_counts'] == 30 and row['total_pcr_replicates'] == 2
    assert populations(row) == {('trafo_3', 10), ('trafo_4', 20)}
    result = aggregate(merged)
    assert len(result) == 1
    assert result.iloc[0]['num_PCR_replicates'] == 2 and result.iloc[0]['bc1_purity'] == 1


def test_purity_denominator_spans_every_population_of_the_barcode():
    merged = merge_data_sources(df_2x300=pd.DataFrame([
        full(10, 'p1', trafo_ID='trafo_3'),
        full(90, 'p2', donor='C' * 130, trafo_ID='trafo_4')]))
    result = aggregate(merged)
    assert len(result) == 2
    assert result['total_bc1_counts'].tolist() == [100, 100]
    assert sorted(result['bc1_purity']) == pytest.approx([0.1, 0.9])


def test_a_population_column_keeps_a_value_only_when_all_observations_agree():
    same = merge_data_sources(df_2x300=pd.DataFrame([
        full(10, 'p1', trafo_ID='trafo_3'), full(5, 'p2', trafo_ID='trafo_3')])).iloc[0]
    mixed = merge_data_sources(df_2x300=pd.DataFrame([
        full(10, 'p1', trafo_ID='trafo_3'), full(5, 'p2', trafo_ID=None)]))
    assert len(mixed) == 1
    assert same['trafo_ID'] == 'trafo_3' and populations(same) == {('trafo_3', 15)}
    assert pd.isna(mixed.iloc[0]['trafo_ID'])
    assert populations(mixed.iloc[0]) == {('trafo_3', 10), (None, 5)}


def short(count, pcr, recovered=None, **scope):
    return dict(counts=count, bc1=BC1, bc0=BC0, SPS=SPS, donor_bc0_fragment=FRAGMENT,
                donor=recovered, pcr_replicate_id=pcr, **scope)


def test_a_fragment_joins_its_recovered_link_without_2x300_data():
    """20250912: yL442 has bc1-donor-bc0 samples but no donor-bc0 amplicon, so its fragments
    cannot recover a donor within their own population; the same bc1 recovered in trafo_1 can."""
    merged = merge_data_sources(df_2x150=pd.DataFrame([
        short(6, 'p1', recovered=DONOR, trafo_ID='trafo_1'), short(4, 'p2', trafo_ID='yL442')]))
    assert len(merged) == 1
    row = merged.iloc[0]
    assert row['donor'] == DONOR and row['total_counts'] == 10 and row['total_pcr_replicates'] == 2
    assert populations(row) == {('trafo_1', 6), ('yL442', 4)}


def test_a_fragment_matching_two_recovered_links_stays_separate():
    other = 'C' + DONOR[1:]  # same 30-nt tail, different donor
    merged = merge_data_sources(df_2x150=pd.DataFrame([
        short(6, 'p1', recovered=DONOR, trafo_ID='trafo_1'),
        short(5, 'p3', recovered=other, trafo_ID='trafo_3'),
        short(4, 'p2', trafo_ID='yL442')]))
    assert len(merged) == 3
    fragment_only = merged[merged['donor'].isna()].iloc[0]
    assert fragment_only['total_counts'] == 4 and fragment_only['short_link_candidates'] == 2


def test_gdb_annotation_ignores_the_sampled_population_but_not_the_construct():
    gdb = pd.DataFrame([dict(donor=DONOR, bc0=BC0, guide='T' * 20, library_ID='V546', trafo_ID='trafo_3',
                             counts=50)])
    seen = merge_data_sources(df_2x300=pd.DataFrame([full(10, 'p1', trafo_ID='trafo_4')]), df_gdb=gdb)
    assert seen.iloc[0]['has_gdb'] and seen.iloc[0]['guide'] == 'T' * 20
    other = merge_data_sources(df_2x300=pd.DataFrame([full(10, 'p1', library_ID='V545')]), df_gdb=gdb)
    assert not other.iloc[0]['has_gdb']


def test_construct_libraries_still_separate_links():
    merged = merge_data_sources(df_2x300=pd.DataFrame([
        full(10, 'p1', V_library='V545'), full(20, 'p2', V_library='V546')]))
    assert len(merged) == 2 and sorted(merged['total_pcr_replicates']) == [1, 1]
