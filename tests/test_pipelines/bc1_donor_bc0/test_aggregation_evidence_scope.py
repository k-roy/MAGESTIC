"""Purity and PCR support describe all scoped observations, before row filters."""
import json
import pandas as pd
import pytest

from magestic.pipelines.bc1_donor_bc0.core.merge import merge_data_sources
from magestic.pipelines.bc1_donor_bc0.pipeline.run_pipeline import aggregate_and_filter


def aggregate(frame, **kwargs):
    return aggregate_and_filter(frame, min_reads=kwargs.pop('min_reads', 1),
                                collapse_sequences=kwargs.pop('collapse_sequences', False),
                                show_progress=False, **kwargs)


@pytest.mark.parametrize('min_reads', [1, 2])
def test_actual_merge_aggregate_does_not_filter_the_purity_denominator(min_reads):
    rows = pd.DataFrame(dict(bc1=['G'*20]*2, donor=['A'*129, 'C'*129], bc0=['A'*10, 'C'*10],
                             counts=[100, 1], pcr_replicate_id=['p1', 'p2']))
    result = aggregate(merge_data_sources(df_2x300=rows), min_reads=min_reads)
    top = result[result['donor'] == 'A'*129].iloc[0]
    assert top['counts'] == 100 and top['total_bc1_counts'] == 101
    assert top['bc1_purity'] == pytest.approx(100/101)
    assert top['support_tier'] == 'provisional'
    assert len(result) == (2 if min_reads == 1 else 1)


def test_unresolved_reads_stay_in_rows_and_in_denominator():
    frame = pd.DataFrame(dict(bc1=['G'*20]*2, donor=['A'*129, None], bc0=['A'*10]*2,
                              counts=[10, 90], pcr_replicate_id=['p1', 'p2']))
    result = aggregate(merge_data_sources(df_2x300=frame))
    assert len(result) == 2 and result['counts'].sum() == 100
    top = result[result['donor'].notna()].iloc[0]
    assert top['bc1_purity'] == pytest.approx(.1)
    assert top['unresolved_bc1_counts'] == 90
    assert top['num_PCR_replicates'] == 1


def test_transformation_scope_survives_aggregation():
    """Revised 2026-09-23 (218 patch 41; Kevin, ledger section 7 item 4).

    This asserted two rows with purity [1, 1] and one PCR each. One linkage in two
    transformations is now one row whose purity and PCR support span both, and the
    transformation scope survives as each population's reads in population_counts.
    """
    frame = pd.DataFrame(dict(bc1=['G'*20]*2, donor=['A'*129]*2, bc0=['A'*10]*2,
                              counts=[10, 20], trafo_ID=['trafo_3', 'trafo_4'],
                              pcr_replicate_id=['p1', 'p2']))
    result = aggregate(merge_data_sources(df_2x300=frame))
    assert len(result) == 1 and result['counts'].tolist() == [30]
    assert result['bc1_purity'].tolist() == [1]
    assert result['num_PCR_replicates'].tolist() == [2]
    assert {(p['trafo_ID'], p['counts']) for p in json.loads(result.iloc[0]['population_counts'])} == {
        ('trafo_3', 10), ('trafo_4', 20)}


def test_collapse_retains_guide_and_uses_pcr_union():
    frame = pd.DataFrame(dict(bc1=['G'*20]*2, donor=['A'*129, 'A'*128+'C'], bc0=['A'*10]*2,
                              counts=[20, 10], guide=['T'*20]*2, SPS=['C'*20]*2,
                              pcr_replicate_id=['p1', 'p1']))
    result = aggregate(frame, collapse_sequences=True)
    assert len(result) == 1
    assert result.iloc[0]['counts'] == 30
    assert result.iloc[0]['guide'] == 'T'*20 and result.iloc[0]['SPS'] == 'C'*20
    assert result.iloc[0]['num_PCR_replicates'] == 1


def test_anonymous_replicate_totals_are_not_summed_by_collapse():
    frame = pd.DataFrame(dict(bc1=['G'*20]*2, donor=['A'*129, 'A'*128+'C'], bc0=['A'*10]*2,
                              counts=[20, 10], num_PCR_replicates=[2, 2]))
    result = aggregate(frame, collapse_sequences=True).iloc[0]
    assert result['counts'] == 30 and result['num_PCR_replicates'] == 2
    assert not result['pcr_support_exact']


def test_json_pcr_identities_survive_exact_aggregation():
    frame = pd.DataFrame(dict(bc1=['G'*20]*2, donor=['A'*129]*2, bc0=['A'*10]*2,
                              counts=[20, 10], pcr_replicate_ids=['["p1", "p2"]', '["p2"]']))
    result = aggregate(frame).iloc[0]
    assert result['counts'] == 30 and result['num_PCR_replicates'] == 2
    assert result['pcr_replicate_ids'] == 'p1;p2'


def test_empty_merged_evidence_is_usable():
    result = aggregate(merge_data_sources())
    assert result.empty and {'counts', 'bc1_purity', 'support_tier'} <= set(result.columns)


def test_missing_barcode_has_no_defined_purity():
    frame = pd.DataFrame(dict(bc1=[None], donor=['A'*129], bc0=['A'*10], counts=[11]))
    result = aggregate(frame)
    assert result['counts'].tolist() == [11]
    assert result['bc1_purity'].isna().all()


def test_sps_does_not_split_a_transformation_barcode_denominator():
    frame = pd.DataFrame(dict(bc1=['G'*20]*2, donor=['A'*129, 'C'*129], bc0=['A'*10]*2,
                              counts=[10, 90], trafo_ID=['trafo_3']*2, SPS=['A'*20, 'C'*20],
                              pcr_replicate_id=['p1', 'p2']))
    result = aggregate(merge_data_sources(df_2x300=frame))
    assert result['total_bc1_counts'].tolist() == [100, 100]
    assert sorted(result['bc1_purity']) == pytest.approx([.1, .9])


@pytest.mark.parametrize('missing', ['donor', 'bc0'])
def test_unresolved_links_cannot_support_barcode_error_collapse(missing):
    frame = pd.DataFrame(dict(bc1=['G'*20, 'G'*19+'A'], donor=['A'*129]*2,
                              bc0=['A'*10]*2, counts=[100, 2], pcr_replicate_id=['p1', 'p2']))
    frame[missing] = None
    result = aggregate(frame, collapse_sequences=True)
    assert len(result) == 2
    assert dict(zip(result['bc1'], result['counts'])) == {'G'*20: 100, 'G'*19+'A': 2}
    assert result['unresolved_bc1_counts'].sum() == 102


def test_one_original_barcode_cannot_supply_two_purity_denominators():
    # The same observed child has two incompatible link-specific representatives.
    # It must not be counted into two supposedly resolved barcode denominators.
    child = 'G'*20
    frame = pd.DataFrame(dict(bc1=['A'+'G'*19, child, 'C'+'G'*19, child],
                              donor=['A'*129]*2+['C'*129]*2, bc0=['A'*10]*4,
                              counts=[100, 2, 80, 3], pcr_replicate_id=['p1','p2','p3','p4']))
    result = aggregate(frame, collapse_sequences=True)
    assert sorted(result['counts']) == [83, 102]
    assert result['counts'].sum() == 185
    assert result['barcode_identity_ambiguous'].all()
    assert result['total_bc1_counts'].isna().all() and result['bc1_purity'].isna().all()


def test_missing_bc0_and_low_read_filter_preserve_unresolved_denominator_mass():
    frame = pd.DataFrame(dict(bc1=['G'*20]*2, donor=['A'*129]*2,
                              bc0=['A'*10, None], counts=[100, 1], pcr_replicate_id=['p1','p2']))
    result = aggregate(frame, collapse_sequences=True, min_reads=2)
    assert len(result) == 1
    row = result.iloc[0]
    assert row['counts'] == 100 and row['total_bc1_counts'] == 101
    assert row['unresolved_bc1_counts'] == 1
    assert row['bc1_purity'] == pytest.approx(100/101)
    assert result.attrs['count_accounting']['reads_below_min_reads'] == 1


def test_actual_merge_collapse_keeps_every_pcr_support_column_coherent():
    frame = pd.DataFrame(dict(bc1=['G'*20]*2, donor=['A'*129, 'A'*128+'C'],
                              bc0=['A'*10]*2, counts=[20, 10], pcr_replicate_id=['p1','p2']))
    row = aggregate(merge_data_sources(df_2x300=frame), collapse_sequences=True).iloc[0]
    assert row['counts'] == row['total_counts'] == row['counts_2x300'] == 30
    assert row['pcr_replicate_ids'] == 'p1;p2'
    assert row['num_PCR_replicates'] == row['total_pcr_replicates'] == row['pcr_replicates_2x300'] == 2
    assert json.loads(row['pcr_ids_2x300']) == ['p1','p2']


def test_filtered_aggregate_cannot_silently_reset_denominator_when_reused():
    frame = pd.DataFrame(dict(bc1=['G'*20]*2, donor=['A'*129]*2,
                              bc0=['A'*10, None], counts=[100, 1], pcr_replicate_id=['p1','p2']))
    first = aggregate(frame, min_reads=2)
    assert first.iloc[0]['total_bc1_counts'] == 101
    try:
        second = aggregate(first)
    except ValueError as error:
        assert 'aggregat' in str(error).lower()
        return  # Refusing incomplete reaggregation is valid; declaring purity1 is not.
    assert second.iloc[0]['total_bc1_counts'] == 101
    assert second.iloc[0]['bc1_purity'] == pytest.approx(100/101)


@pytest.mark.parametrize('add_neighbor', [False, True])
def test_guide_candidates_survive_both_collapse_stages(add_neighbor):
    frame = pd.DataFrame(dict(bc1=['G'*20]*2, donor=['A'*129, 'A'*128+'C'],
                              bc0=['A'*10]*2, counts=[20, 10], pcr_replicate_id=['p1','p2']))
    gdb = pd.DataFrame(dict(donor=['A'*129, 'A'*128+'C'], bc0=['A'*10]*2,
                            guide=['T'*20, 'C'*20], counts=[200,100]))
    if add_neighbor:
        frame.loc[len(frame)] = ['G'*19+'A', 'A'*129, 'A'*10, 5, 'p3']
    row = aggregate(merge_data_sources(df_2x300=frame, df_gdb=gdb), collapse_sequences=True).iloc[0]
    assert pd.isna(row['guide']) and row['guide_evidence_conflict']
    assert set(json.loads(row['guide_candidates'])) == {'T'*20, 'C'*20}


def test_plasmid_pool_identity_is_not_a_separate_barcode_population():
    frame = pd.DataFrame(dict(bc1=['G'*20]*2, donor=['A'*129, 'C'*129],
                              bc0=['A'*10]*2, counts=[10,90], pcr_replicate_id=['p1','p2'],
                              yeast_library=['yL_fixture']*2, trafo_ID=['trafo_3']*2,
                              V_library=['V545','V546'], SPS=['A'*20,'C'*20]))
    result = aggregate(merge_data_sources(df_2x300=frame))
    assert len(result) == 2 and set(result['V_library']) == {'V545','V546'}
    assert result['total_bc1_counts'].tolist() == [100,100]
    assert sorted(result['bc1_purity']) == pytest.approx([.1,.9])


def test_source_specific_support_unions_do_not_double_count_shared_pcr():
    full = pd.DataFrame(dict(bc1=['G'*20]*2, donor=['A'*129, 'A'*128+'C'],
                             bc0=['A'*10]*2, counts=[20,10], pcr_replicate_id=['p1','p2']))
    recovered = pd.DataFrame(dict(bc1=['G'*20], donor=['A'*129], donor_bc0_fragment=['A'*40],
                                  bc0=['A'*10], counts=[5], pcr_replicate_id=['p1']))
    row = aggregate(merge_data_sources(df_2x300=full, df_2x150=recovered),
                    collapse_sequences=True).iloc[0]
    assert row['counts'] == 35 and row['counts_2x300'] == 30 and row['counts_2x150'] == 5
    assert row['num_PCR_replicates'] == row['total_pcr_replicates'] == 2
    assert row['pcr_replicates_2x300'] == 2 and row['pcr_replicates_2x150'] == 1
    assert json.loads(row['pcr_ids_2x300']) == ['p1','p2']
    assert json.loads(row['pcr_ids_2x150']) == ['p1']


@pytest.mark.parametrize('add_neighbor_donor', [False, True])
def test_ambiguous_gdb_candidate_evidence_survives_aggregation(add_neighbor_donor):
    frame = pd.DataFrame(dict(bc1=['G'*20], donor=['A'*129], bc0=['A'*10],
                              counts=[10], pcr_replicate_id=['p1']))
    gdb = pd.DataFrame(dict(donor=['A'*129]*2, bc0=['A'*10]*2,
                            guide=['T'*20,'C'*20], counts=[200,100]))
    expected = {'T'*20,'C'*20}
    if add_neighbor_donor:
        frame.loc[len(frame)] = ['G'*20,'A'*128+'C','A'*10,20,'p2']
        gdb.loc[len(gdb)] = ['A'*128+'C','A'*10,'G'*20,300]
        expected.add('G'*20)
    merged = merge_data_sources(df_2x300=frame, df_gdb=gdb)
    row = aggregate(merged, collapse_sequences=True).iloc[0]
    assert pd.isna(row['guide']) and row['guide_evidence_conflict']
    assert set(json.loads(row['guide_candidates'])) == expected
    assert {item['guide'] for item in json.loads(row['gdb_candidates'])} == expected


@pytest.mark.parametrize('collapse', [False, True])
def test_complete_native_annotation_batches_retain_all_candidate_records(collapse):
    # Complete, disjoint PCR observations annotated against separate GDB tables.
    # This is an accepted input composition, not the runner's single-table route.
    outputs = []
    for barcode, count, pcr, guide, library, scaffold, support in [
        ('G'*20, 20, 'p1', None, 'L1', 'S1', 1000),
        ('G'*19+'A' if collapse else 'G'*20, 10, 'p2', 'G'*20, 'L2', 'S2', 10),
    ]:
        observations = pd.DataFrame(dict(bc1=[barcode], donor=['A'*129], bc0=['A'*10],
                                         counts=[count], pcr_replicate_id=[pcr]))
        annotations = pd.DataFrame(dict(donor=['A'*129], bc0=['A'*10], guide=[guide],
                                        library_ID=[library], scaffold=[scaffold], counts=[support]))
        outputs.append(merge_data_sources(df_2x300=observations, df_gdb=annotations))
    complete = pd.concat(outputs, ignore_index=True)
    expected = {json.dumps(item, sort_keys=True) for value in complete['gdb_candidates']
                for item in json.loads(value)}
    result = aggregate(complete, collapse_sequences=collapse)
    assert len(result) == 1
    row = result.iloc[0]
    assert row['counts'] == 30 and row['num_PCR_replicates'] == 2
    assert row['guide'] == 'G'*20 and row['library_id'] == 'L2' and row['scaffold'] == 'S2'
    assert {json.dumps(item, sort_keys=True) for item in json.loads(row['gdb_candidates'])} == expected
