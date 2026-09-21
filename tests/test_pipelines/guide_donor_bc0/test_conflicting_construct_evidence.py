"""A unique guide and a unique full donor must not name different constructs."""
import pandas as pd
import pytest

from magestic.pipelines.guide_donor_bc0.core import oligo_matching as matching

GUIDE_A, GUIDE_B = 'A'*20, 'C'*20
DONOR_A, DONOR_B = 'A'*30+'C'*60+'G'*39, 'C'*30+'G'*60+'T'*39


@pytest.fixture
def lookups():
    designs = pd.DataFrame([
        {'oligo_name':'A','guide':GUIDE_A,'donor':DONOR_A},
        {'oligo_name':'B','guide':GUIDE_B,'donor':DONOR_B},
    ])
    designs['middle_donor'] = designs.donor.str[30:90]
    return matching.build_lookup_dictionaries(designs)


@pytest.mark.parametrize('guide,donor', [(GUIDE_A,DONOR_B),(GUIDE_B,DONOR_A)])
def test_conflicting_unique_evidence_is_unassigned(lookups,tmp_path,guide,donor):
    result = matching.resolve_oligo_name(guide,donor,lookups)
    assert result.oligo_name is None
    assert result.confidence == 0
    assert result.match_type == 'conflicting'
    input_path, output_path = tmp_path/'counts.tsv',tmp_path/'matched.tsv'
    pd.DataFrame([dict(library_ID='V_test',total_counts=17,guide=guide,donor=donor,bc0='ACGTACGTAC')]).to_csv(input_path,sep='\t',index=False)
    stats = matching.process_counts_file(input_path,output_path,lookups,chunk_size=1)
    output = pd.read_csv(output_path,sep='\t',keep_default_na=False)
    assert len(output) == 1 and output['counts'].sum() == 17
    assert output.iloc[0]['oligo_name'] == 'NA'
    assert output.iloc[0]['match_type'] == 'conflicting'
    assert bool(output.iloc[0]['perfect_guide']) and bool(output.iloc[0]['perfect_donor'])
    assert stats == {'total':1,'matched':0,'perfect_match':0}


@pytest.mark.parametrize('guide,donor,name,kind',[
    (GUIDE_A,DONOR_A,'A','exact'),
    (GUIDE_A,'T'*129,'A','guide_only'),
    ('G'*20,DONOR_B,'B','donor_only'),
    ('G'*20,'T'+DONOR_A[1:],'A','partial'),
])
def test_compatible_and_single_evidence_routes_are_preserved(lookups,tmp_path,guide,donor,name,kind):
    result = matching.resolve_oligo_name(guide,donor,lookups)
    assert (result.oligo_name,result.match_type) == (name,kind)
    input_path, output_path = tmp_path/'counts.tsv',tmp_path/'matched.tsv'
    pd.DataFrame([dict(library_ID='V_test',total_counts=17,guide=guide,donor=donor,bc0='ACGTACGTAC')]).to_csv(input_path,sep='\t',index=False)
    stats = matching.process_counts_file(input_path,output_path,lookups,chunk_size=1)
    output = pd.read_csv(output_path,sep='\t',keep_default_na=False)
    assert output['oligo_name'].tolist() == [name]
    assert output['counts'].sum() == 17 and stats['matched'] == 1
