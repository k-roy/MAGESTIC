"""Explicit sample/assay/plate joins must not invent or multiply identities."""
import pandas as pd
import pytest
from magestic.pipelines.redi.readers.keyfiles import combine_keyfiles,load_redi_bc_keys

def write(tmp_path,name,rows):
    p=tmp_path/name;pd.DataFrame(rows).to_csv(p,sep='\t',index=False);return p

def keys(tmp_path,mutate=None):
    sample=[dict(gDNA_plate_name='dated_plate',sample_number='sample_1',REDI_plate_barcode='3500'),
            dict(gDNA_plate_name='dated_plate',sample_number='sample_2',REDI_plate_barcode='empty')]
    gdna=[dict(gDNA_plate_name='dated_plate',outer_primers='plate_4A',inner_primers='KR1967-KR2225',purpose='yT182 REDI'),
          dict(gDNA_plate_name='dated_plate',outer_primers='plate_4B',inner_primers='KR1967-KR2227',purpose='yT183 REDI')]
    plate=[dict(REDI_plate_barcode='3500',MAT_alpha_array_1536_barcode='6100',MAT_a_array_1536_barcode='6000')]
    if mutate:mutate(sample,gdna,plate)
    return [write(tmp_path,name,rows) for name,rows in zip(('sample.tsv','gdna.tsv','plate.tsv'),(sample,gdna,plate))]

def test_declared_multiple_assays_expand_but_empty_controls_do_not_join(tmp_path):
    out=combine_keyfiles(*keys(tmp_path))
    assert len(out)==4
    called=out[out.sample_number=='sample_1'];empty=out[out.sample_number=='sample_2']
    assert called.MAT_alpha_array_1536_barcode.tolist()==['6100','6100']
    assert empty.MAT_alpha_array_1536_barcode.isna().all()
    assert empty.registration_status.tolist()==['declared_empty','declared_empty']

@pytest.mark.parametrize('table',['sample','gdna','plate'])
def test_duplicate_registration_refused(tmp_path,table):
    def mutate(s,g,p):
        rows={'sample':s,'gdna':g,'plate':p}[table];rows.append(dict(rows[0]))
    with pytest.raises(ValueError,match='duplicate|unique'):
        combine_keyfiles(*keys(tmp_path,mutate))

@pytest.mark.parametrize('field,value',[('gDNA_plate_name','unknown'),('REDI_plate_barcode','typo'),('REDI_plate_barcode','')])
def test_unregistered_sample_metadata_refused(tmp_path,field,value):
    with pytest.raises(ValueError,match='missing|register|key'):
        combine_keyfiles(*keys(tmp_path,lambda s,g,p:s[0].update({field:value})))

def test_placeholder_strings_never_coerce_to_matching_nan(tmp_path):
    def mutate(s,g,p):
        s[0]['REDI_plate_barcode']='wrong';p[0]['REDI_plate_barcode']='different_wrong'
    with pytest.raises(ValueError,match='register|key'):
        combine_keyfiles(*keys(tmp_path,mutate))

def test_leading_zero_plate_identifier_preserved(tmp_path):
    def mutate(s,g,p):
        s[0]['REDI_plate_barcode']='03500';p[0]['REDI_plate_barcode']='03500'
    out=combine_keyfiles(*keys(tmp_path,mutate))
    assert set(out.loc[out.sample_number=='sample_1','REDI_plate_barcode'])=={'03500'}

def test_integral_spreadsheet_plate_normalization(tmp_path):
    out=combine_keyfiles(*keys(tmp_path,lambda s,g,p:s[0].update(REDI_plate_barcode='3500.0')))
    assert out.loc[out.sample_number=='sample_1','MAT_alpha_array_1536_barcode'].tolist()==['6100','6100']

def test_per_well_primer_assignment_prevents_cross_expansion(tmp_path):
    def mutate(s,g,p):
        for row in s:row['inner_primers']='KR1967-KR2225'
    out=combine_keyfiles(*keys(tmp_path,mutate))
    assert len(out)==2 and set(out.outer_primers)=={'plate_4A'}

def test_two_inner_primers_cannot_implicitly_share_one_output_sample(tmp_path):
    def mutate(s,g,p):g[1]['outer_primers']='plate_4A'
    with pytest.raises(ValueError,match='ambig|duplicate|unique'):
        combine_keyfiles(*keys(tmp_path,mutate))

def test_foreign_metadata_cannot_silently_get_suffixes(tmp_path):
    with pytest.raises(ValueError,match='overlap|collision'):
        combine_keyfiles(*keys(tmp_path,lambda s,g,p:p[0].update(purpose='false_purpose')))

def barcode_keys(tmp_path,rows=None):
    return write(tmp_path,'arbitrary_filename.tsv',rows or [dict(REDI_library='yT182',REDI_bc1_seq='A'*20,
        plate_384='RED_384',row_384='A',col_384=1)])

def test_barcode_registration_requires_unique_library_and_sequence(tmp_path):
    rows=[dict(REDI_library='yT182',REDI_bc1_seq='A'*20,plate_384='RED_384',row_384='A',col_384=c) for c in (1,2)]
    with pytest.raises(ValueError,match='duplicate|unique'):
        load_redi_bc_keys(yt182_yt183_key_path=barcode_keys(tmp_path,rows))

def test_library_namespace_keeps_same_barcode_in_two_libraries(tmp_path):
    rows=[dict(REDI_library=lib,REDI_bc1_seq='A'*20,plate_384='RED_384',row_384='A',col_384=c) for lib,c in [('yT182',1),('yT183',2)]]
    out=load_redi_bc_keys(yt182_yt183_key_path=barcode_keys(tmp_path,rows))
    assert len(out)==2 and set(out.REDI_library)=={'yT182','yT183'}

@pytest.mark.parametrize('column,value',[('REDI_bc1_seq','A'*19+'N'),('row_384','Q'),('col_384',25)])
def test_invalid_barcode_position_registration_refused(tmp_path,column,value):
    row=dict(REDI_library='yT182',REDI_bc1_seq='A'*20,plate_384='RED_384',row_384='A',col_384=1);row[column]=value
    with pytest.raises(ValueError):load_redi_bc_keys(yt182_yt183_key_path=barcode_keys(tmp_path,[row]))
