"""A BDB link lives in one transformation (218 patch 45, superseding patch 41's link key).

Kevin, 2026-09-23: "For a given trafo, we require that a bc1 be matched up to a given donor by at
least two independent PCR wells. This guards against spurious associations due to PCR recombination.
When comparing two separate trafos, we simply want to make sure that the bc1 is observed in one trafo
with much higher read counts/purity than in the other."

The bc1 PCR fragment is joined to the cut guide-donor plasmid in yeast at each transformation, so a
linkage is one assembly event in one transformation. The loaders map each key row's sample label to
its transformation (through the screen's trafo key for a pooled outgrowth, and 'trafo_N-M' plates to
'trafo_N') and keep the label in trafo_label. Links, PCR support, the purity denominator and donor
recovery are keyed by transformation; the labels a link was seen under are evidence in
population_counts.

Design: planning/218_magestic_production_corrections/reviews/merge_regressions_repair/DESIGN.md §6.
Each test asserts the link grain first, so at a parent with patch 41's population-free key it fails
on that behaviour.
"""
import json
from urllib.parse import unquote

import pandas as pd
import pytest

from magestic.pipelines.bc1_donor_bc0.core.merge import merge_data_sources
from magestic.pipelines.bc1_donor_bc0.pipeline.run_pipeline import aggregate_and_filter
from magestic.pipelines.bc1_donor_bc0.readers.data_loader import (
    load_2x150_bc1_donor_bc0_data, load_donor_bc0_lookup, recover_full_donors)

SPS = 'CGTATGGAGATGACGTCAGG'
BC1, BC0 = 'G' * 20, 'A' * 10
DONOR = 'T' * 100 + 'ACGTACGTACGTACGTACGTACGTACGTAC'
FRAGMENT = DONOR[-30:] + SPS + BC0


def aggregate(frame, collapse=False):
    return aggregate_and_filter(frame, min_reads=1, min_pcr_replicates=2,
                                collapse_sequences=collapse, show_progress=False)


def labels(row):
    entries = json.loads(row['population_counts'])
    assert {item['trafo_ID'] for item in entries} == {row['trafo_ID']}
    return {(item['trafo_label'], item['counts']) for item in entries}


def _layout(tmp_path, wells, trafo_key=None, donor_wells=()):
    """A 2x150 run on plate_3F; `wells` is (sample, trafo label, reads); donor_wells add donor-bc0 rows.

    The donor-bc0 tables carry SPS, as the real complete_donor_bc0_counts tables do (20250912 plate_3A).
    """
    data, keys = tmp_path / 'data', tmp_path / 'keys'
    data.mkdir()
    keys.mkdir()
    rows = [dict(description='bc1-donor-bc0', sequencing_date='20250101', outer_primers='plate_3F',
                 sample_number=well, condition='plate wash', trafo_ID=label, gDNA_plate_name='prep_A',
                 inner_primers='KR1967-KR1590', polymerase='Q5') for well, label, _ in wells]
    rows += [dict(description='donor-bc0', sequencing_date='20250101', outer_primers='plate_3A',
                  sample_number=well, condition='plate wash', trafo_ID=label, gDNA_plate_name='prep_A',
                  inner_primers='KR1967-KR1590', polymerase='Q5') for well, label, _ in donor_wells]
    pd.DataFrame(rows).to_csv(keys / 'bc1_donor_bc0_combined_key.tsv', sep='\t', index=False)
    for well, _, count in wells:
        pd.DataFrame([dict(counts=count, bc1=BC1, bc0=BC0, SPS=SPS, donor_bc0_fragment=FRAGMENT)]).to_csv(
            data / f'20250101_plate_3F_{well}_bc1_donor_bc0_fragment_counts.tsv', sep='\t', index=False)
    for well, _, count in donor_wells:
        pd.DataFrame([dict(counts=count, donor=DONOR, donor_bc0_fragment=FRAGMENT, SPS=SPS, bc0=BC0)]).to_csv(
            data / f'20250101_plate_3A_{well}_complete_donor_bc0_counts.tsv', sep='\t', index=False)
    if trafo_key:
        pd.DataFrame([dict(yL=yl, trafo_ID=trafo, replicate='rep_1') for yl, trafo in trafo_key.items()]).to_csv(
            keys / 'bc1_donor_bc0_amplicons_trafo_key.tsv', sep='\t', index=False)
    return data, keys


YL406_LIKE = [(2, 'trafo_2-1', 7), (14, 'trafo_2-2', 3), (50, 'yL387', 5)]


def test_plates_and_the_pooled_outgrowth_of_one_transformation_are_one_link(tmp_path):
    """yL406: the washes of plates trafo_2-1 to 2-3 were pooled into the yL387 outgrowth."""
    data, keys = _layout(tmp_path, YL406_LIKE, trafo_key={'yL387': 'trafo_2'})
    frame, _ = load_2x150_bc1_donor_bc0_data(data, keys)
    merged = merge_data_sources(df_2x150=frame)
    assert len(merged) == 1
    row = merged.iloc[0]
    assert row['trafo_ID'] == 'trafo_2' and row['total_counts'] == 15 and row['total_pcr_replicates'] == 3
    assert labels(row) == {('trafo_2-1', 7), ('trafo_2-2', 3), ('yL387', 5)}
    row = aggregate(merged).iloc[0]
    assert row['num_PCR_replicates'] == 3 and row['support_tier'] == 'confident' and row['bc1_purity'] == 1


def test_without_a_trafo_key_a_pooled_outgrowth_is_its_own_scope(tmp_path):
    data, keys = _layout(tmp_path, YL406_LIKE)
    merged = merge_data_sources(df_2x150=load_2x150_bc1_donor_bc0_data(data, keys)[0])
    assert len(merged) == 2
    support = dict(zip(merged['trafo_ID'], merged['total_pcr_replicates']))
    assert support == {'trafo_2': 2, 'yL387': 1}


def test_a_pooled_outgrowth_keeps_its_amplicon_level_in_the_pcr_identity(tmp_path):
    data, keys = _layout(tmp_path, YL406_LIKE, trafo_key={'yL387': 'trafo_2'})
    frame, _ = load_2x150_bc1_donor_bc0_data(data, keys)
    levels = {label: json.loads(unquote(pcr[len('pcr:v1:'):]))[5]
              for label, pcr in zip(frame['trafo_label'], frame['pcr_replicate_id'])}
    assert levels == {'trafo_2-1': 'trafo', 'trafo_2-2': 'trafo', 'yL387': 'pooled_yL'}
    assert set(frame['trafo_ID']) == {'trafo_2'}


def test_conflicting_trafo_keys_are_refused(tmp_path):
    data, keys = _layout(tmp_path, YL406_LIKE, trafo_key={'yL387': 'trafo_2'})
    pd.DataFrame([dict(trafo_ID='trafo_3', yL='yL387')]).to_csv(keys / 'trafo_key.tsv', sep='\t', index=False)
    with pytest.raises(ValueError, match='yL387'):
        load_2x150_bc1_donor_bc0_data(data, keys)


def full(count, pcr, donor=DONOR, **scope):
    return dict(counts=count, bc1=BC1, bc0=BC0, donor=donor, SPS=SPS, pcr_replicate_id=pcr, **scope)


def test_different_transformations_of_one_construct_library_stay_separate():
    """20250912: trafo_3 (SpG) and trafo_4 (SpCas9) share cV586; each is its own assembly."""
    merged = merge_data_sources(df_2x300=pd.DataFrame([
        full(10, 'p1', trafo_ID='trafo_3'), full(20, 'p2', trafo_ID='trafo_4')]))
    assert len(merged) == 2
    assert dict(zip(merged['trafo_ID'], merged['total_pcr_replicates'])) == {'trafo_3': 1, 'trafo_4': 1}


def test_purity_is_per_transformation():
    merged = merge_data_sources(df_2x300=pd.DataFrame([
        full(10, 'p1', trafo_ID='trafo_3'), full(90, 'p2', donor='C' * 130, trafo_ID='trafo_4')]))
    result = aggregate(merged)
    assert sorted(result['total_bc1_counts']) == [10, 90]
    assert result['bc1_purity'].tolist() == [1, 1]


def test_an_empty_well_is_its_own_scope():
    """A read in a well with no sample cannot add a PCR to a real link."""
    merged = merge_data_sources(df_2x300=pd.DataFrame([
        full(10, 'p1', trafo_ID='trafo_1'), full(1, 'p2', trafo_ID='empty')]))
    assert len(merged) == 2
    real = merged[merged['trafo_ID'] == 'trafo_1'].iloc[0]
    assert real['total_pcr_replicates'] == 1
    assert aggregate(merged).set_index('trafo_ID').loc['trafo_1', 'support_tier'] == 'provisional'


def test_donor_recovery_crosses_the_labels_of_one_transformation(tmp_path):
    """The yL387 fragment recovers from the donor-bc0 amplicon of plate trafo_2-1 (same transformation)."""
    data, keys = _layout(tmp_path, [(50, 'yL387', 5)], trafo_key={'yL387': 'trafo_2'},
                         donor_wells=[(2, 'trafo_2-1', 11)])
    frame, _ = load_2x150_bc1_donor_bc0_data(data, keys)
    recovered = recover_full_donors(frame, load_donor_bc0_lookup(data, keys))
    assert recovered['donor_recovery_status'].tolist() == ['recovered_unique']
    assert recovered['donor'].tolist() == [DONOR]


def test_donor_recovery_does_not_cross_transformations(tmp_path):
    """Control: without the trafo key, yL387 is its own scope and trafo_2-1's amplicon does not apply."""
    data, keys = _layout(tmp_path, [(50, 'yL387', 5)], donor_wells=[(2, 'trafo_2-1', 11)])
    frame, _ = load_2x150_bc1_donor_bc0_data(data, keys)
    recovered = recover_full_donors(frame, load_donor_bc0_lookup(data, keys))
    assert recovered['donor_recovery_status'].tolist() == ['no_positive_candidate']


def test_no_short_link_attaches_after_a_2x150_only_merge():
    """A fragment spelled with a shorter tail stays its own observation (no stacked inference)."""
    short = DONOR[-25:] + SPS + BC0
    merged = merge_data_sources(df_2x150=pd.DataFrame([
        dict(counts=6, bc1=BC1, bc0=BC0, SPS=SPS, donor_bc0_fragment=FRAGMENT, donor=DONOR,
             pcr_replicate_id='p1', trafo_ID='trafo_1'),
        dict(counts=4, bc1=BC1, bc0=BC0, SPS=SPS, donor_bc0_fragment=short, donor=None,
             pcr_replicate_id='p2', trafo_ID='trafo_1')]))
    assert len(merged) == 2


def test_gdb_annotation_ignores_the_sample_scope_but_not_the_construct():
    """GDB is plasmid-level: a guide/donor/bc0 association exists before any transformation."""
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
