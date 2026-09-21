"""A reference-amplicon candidate must never multiply observed BDB read mass."""
import json
import pandas as pd
import pytest

from magestic.pipelines.bc1_donor_bc0.readers.data_loader import load_donor_bc0_lookup, recover_full_donors

FRAGMENT = 'C'*30 + 'G'*20 + 'A'*10
D1, D2 = 'A'*99+'C'*30, 'T'*99+'C'*30


def observation(**extra):
    return pd.DataFrame([dict(bc1='G'*20, bc0='A'*10, donor_bc0_fragment=FRAGMENT,
                               counts=20, pcr_replicate_id='p1', **extra)])


def lookup(donors, **extra):
    return pd.DataFrame(dict(donor_bc0_fragment=[FRAGMENT]*len(donors), donor=donors,
                              bc0=['A'*10]*len(donors), counts=[5]*len(donors), **extra))


def test_two_lookup_donors_do_not_duplicate_one_barcode_observation():
    result = recover_full_donors(observation(), lookup([D1,D2]))
    assert len(result) == 1 and result['counts'].sum() == 20
    assert pd.isna(result.iloc[0]['donor'])
    assert result.iloc[0]['donor_recovery_status'] == 'ambiguous'
    assert set(json.loads(result.iloc[0]['donor_recovery_candidates'])) == {D1,D2}


def test_repeated_same_donor_is_one_candidate_and_never_adds_reference_reads():
    result = recover_full_donors(observation(), lookup([D1,D1]))
    assert len(result) == 1 and result['counts'].sum() == 20
    assert result.iloc[0]['donor'] == D1


def test_unique_lookup_donor_control():
    result = recover_full_donors(observation(), lookup([D1]))
    assert len(result) == 1 and result.iloc[0]['donor'] == D1
    assert result['counts'].tolist() == [20]


def test_same_fragment_in_another_transformation_cannot_supply_a_donor():
    result = recover_full_donors(observation(trafo_ID='trafo_3'),
                                 lookup([D1], trafo_ID=['trafo_4']))
    assert len(result) == 1 and result['counts'].tolist() == [20]
    assert pd.isna(result.iloc[0]['donor'])


def test_known_donor_is_preserved_when_recovery_annotation_conflicts():
    result = recover_full_donors(observation(donor=D1), lookup([D2]))
    assert len(result) == 1 and result['counts'].tolist() == [20]
    assert result.iloc[0]['donor'] == D1
    assert result.iloc[0]['donor_recovery_status'] == 'conflict_with_observed'


def test_actual_donor_loader_preserves_competing_candidates(tmp_path):
    data, keys = tmp_path/'data', tmp_path/'keys'
    data.mkdir(); keys.mkdir()
    pd.DataFrame([dict(description='donor-bc0', sequencing_date='20250101',
                       outer_primers='plate_1A', sample_number=1)]).to_csv(
        keys/'bc1_donor_bc0_combined_key.tsv',sep='\t',index=False)
    lookup([D1,D2]).to_csv(data/'20250101_plate_1A_1_complete_donor_bc0_counts.tsv',sep='\t',index=False)
    candidates = load_donor_bc0_lookup(data,keys)
    assert len(candidates) == 2 and set(candidates['donor']) == {D1,D2}
    assert candidates['counts'].sum() == 10
    result = recover_full_donors(observation(), candidates)
    assert len(result) == 1 and result['counts'].sum() == 20 and result['donor'].isna().all()


def test_missing_fragment_is_not_a_matching_identity():
    obs, ref = observation(), lookup([D1])
    obs['donor_bc0_fragment'] = None
    ref['donor_bc0_fragment'] = None
    result = recover_full_donors(obs, ref)
    assert result['donor'].isna().all() and result['counts'].tolist() == [20]


def test_zero_reference_reads_do_not_establish_a_recovery_candidate():
    ref = lookup([D1,D2])
    ref['counts'] = [5,0]
    result = recover_full_donors(observation(), ref)
    assert result.iloc[0]['donor'] == D1
    assert result.iloc[0]['donor_recovery_support_reads'] == 5
    assert result['counts'].tolist() == [20]


def test_positive_unknown_reference_prevents_unique_recovery():
    result = recover_full_donors(observation(), lookup([D1,None]))
    row = result.iloc[0]
    assert pd.isna(row['donor']) and row['donor_recovery_status'] == 'unresolved_reference'
    assert row['donor_recovery_support_reads'] == 10 and row['donor_recovery_unresolved_reads'] == 5
    assert row['counts'] == 20


def test_unknown_scope_is_not_a_wildcard():
    result = recover_full_donors(observation(), lookup([D1], trafo_ID=['trafo_3']))
    assert result['donor'].isna().all() and result['counts'].tolist() == [20]


def test_empty_reference_and_empty_observation_have_usable_schemas():
    result = recover_full_donors(observation(), pd.DataFrame())
    assert result['donor'].isna().all() and result['counts'].tolist() == [20]
    empty = recover_full_donors(observation().iloc[:0], lookup([D1]))
    assert empty.empty and {'donor','donor_recovery_status','counts'} <= set(empty)


@pytest.mark.parametrize('count', [-1, .5, float('nan')])
def test_invalid_reference_support_is_refused(count):
    ref = lookup([D1]); ref['counts'] = count
    with pytest.raises(ValueError, match='count'):
        recover_full_donors(observation(), ref)


def test_duplicate_dataframe_indexes_do_not_change_row_mass():
    obs = pd.concat([observation(), observation()])
    result = recover_full_donors(obs, lookup([D1]))
    assert result.index.tolist() == [0,0] and result['counts'].tolist() == [20,20]
    assert result['donor'].tolist() == [D1,D1]


@pytest.mark.parametrize('conflicting', [False, True])
def test_one_physical_reference_file_cannot_count_twice_or_have_two_scopes(tmp_path, conflicting):
    data, keys = tmp_path/'data', tmp_path/'keys'
    data.mkdir(); keys.mkdir()
    row = dict(description='donor-bc0', sequencing_date='20250101', outer_primers='plate_1A',
               sample_number=1, trafo_ID='trafo_3')
    second = dict(row, trafo_ID='trafo_4') if conflicting else dict(row)
    pd.DataFrame([row,second]).to_csv(keys/'bc1_donor_bc0_combined_key.tsv',sep='\t',index=False)
    source = data/'20250101_plate_1A_1_complete_donor_bc0_counts.tsv'
    lookup([D1]).to_csv(source,sep='\t',index=False)
    if conflicting:
        with pytest.raises(ValueError, match='registration'):
            load_donor_bc0_lookup(data,keys)
    else:
        candidates = load_donor_bc0_lookup(data,keys)
        assert candidates['counts'].tolist() == [5]
        assert json.loads(candidates.iloc[0]['donor_recovery_source_files']) == [str(source.resolve())]
        result = recover_full_donors(observation(trafo_ID='trafo_3'),candidates)
        assert result['counts'].tolist() == [20]
        assert result.iloc[0]['donor_recovery_support_reads'] == 5
