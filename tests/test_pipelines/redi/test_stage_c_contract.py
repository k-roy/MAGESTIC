"""Real TSV-to-TSV Stage C fixtures with independently declared molecules."""
import json
from pathlib import Path
import pandas as pd
import pytest
from magestic.pipelines.redi.core.config import REDIConfig,REDIPaths
from magestic.pipelines.redi.pipeline.step_c_process_bcs import run_step_c

BC='ACGTATGCCATTCAGGTAAC'; RB='ACCCATAGTCCATTCGGCTA'
def setup(tmp_path):
    paths=REDIPaths(tmp_path,tmp_path/'raw',tmp_path/'intermediate',tmp_path/'processed',tmp_path/'keys')
    paths.keyfiles.mkdir();paths.counts_tables.mkdir(parents=True)
    config=REDIConfig('synthetic','20250627',paths,redi_bc_keys=['arbitrary_keys.tsv'])
    sample='20250627_plate_4A_sample_1'
    meta=pd.DataFrame([dict(sample_number='sample_1',gDNA_plate_name='dated_plate',
        inner_primers='KR1967-KR2225',outer_primers='plate_4A',purpose='yT182 REDI',
        REDI_library='yT182',MAT_a_barcode_prefix='AAGTATGT',MAT_a_barcode_length=20,
        REDI_plate_barcode=3567,MAT_alpha_array_1536_barcode=6156,
        registration_status='registered')])
    combined=paths.keyfiles/'combined.tsv';meta.to_csv(combined,sep='\t',index=False)
    key=pd.DataFrame([dict(REDI_library='yT182',REDI_bc1_seq=RB,plate_384='RED_384',row_384='A',col_384=1)])
    key.to_csv(paths.keyfiles/'arbitrary_keys.tsv',sep='\t',index=False)
    rows=pd.DataFrame([dict(sample_name=sample,R1_seq='CAGGTCATGCTC'+BC+'CCCTAGGATACG',
                           R2_seq='AAGTATGT'+RB+'GATCATAACTTC',counts=7)])
    countpath=paths.counts_tables/f'{sample}_R1_R2_counts.tsv';rows.to_csv(countpath,sep='\t',index=False)
    return config,combined,meta,key,countpath,rows

def test_stage_c_loads_keys_by_schema_and_accounts_every_read(tmp_path):
    c,k,meta,key,p,rows=setup(tmp_path)
    out,_=run_step_c(c,k);result=pd.read_csv(out,sep='\t')
    assert len(result)==1 and result.bc1.iloc[0]==BC and result.counts.iloc[0]==7
    assert result.REDI_row_384.iloc[0]=='A'
    provenance=json.loads(Path(str(out)+'.provenance.json').read_text())
    assert provenance['counts']['input_reads']==7
    assert provenance['counts']['assigned_redi_reads']==7
    assert provenance['counts']['unassigned_reads']==0

def test_cross_library_same_barcode_does_not_multiply_rows(tmp_path):
    c,k,meta,key,p,rows=setup(tmp_path)
    other=key.copy();other.REDI_library='yT183';other.col_384=2
    pd.concat([key,other]).to_csv(c.paths.keyfiles/'arbitrary_keys.tsv',sep='\t',index=False)
    out,_=run_step_c(c,k);result=pd.read_csv(out,sep='\t')
    assert len(result)==1 and result.REDI_library.iloc[0]=='yT182' and result.REDI_col_384.iloc[0]==1

@pytest.mark.parametrize('kind',['duplicate_metadata','unregistered_sample','ambiguous_primer','missing_primer','negative_count','fractional_count','missing_sample_file'])
def test_stage_c_refuses_invalid_join_or_counts_before_publication(tmp_path,kind):
    c,k,meta,key,p,rows=setup(tmp_path)
    if kind=='duplicate_metadata':pd.concat([meta,meta]).to_csv(k,sep='\t',index=False)
    elif kind=='unregistered_sample':rows.sample_name='unknown';rows.to_csv(p,sep='\t',index=False)
    elif kind=='ambiguous_primer':c.primer_rules={'KR1967':'R1','KR2225':'R2'}
    elif kind=='missing_primer':c.primer_rules={'KR9999':'R1'}
    elif kind=='negative_count':rows.counts=-1;rows.to_csv(p,sep='\t',index=False)
    elif kind=='fractional_count':rows.counts=1.5;rows.to_csv(p,sep='\t',index=False)
    elif kind=='missing_sample_file':
        other=meta.copy();other.sample_number='sample_2';pd.concat([meta,other]).to_csv(k,sep='\t',index=False)
    with pytest.raises(ValueError):run_step_c(c,k)
    assert not (c.paths.combined_data/'combined_REDI_bc_MAGESTIC_bc_merged_df.tsv').exists()

def test_unmatched_molecular_evidence_is_counted_without_false_call(tmp_path):
    c,k,meta,key,p,rows=setup(tmp_path)
    bad=rows.copy();bad.R1_seq='T'*100;bad.counts=3
    unknown=rows.copy();unknown.R2_seq='AAGTATGT'+'T'*20+'GATCATAACTTC';unknown.counts=5
    pd.concat([rows,bad,unknown]).to_csv(p,sep='\t',index=False)
    out,_=run_step_c(c,k);result=pd.read_csv(out,sep='\t')
    assert len(result)==1 and result.counts.sum()==7
    provenance=json.loads(Path(str(out)+'.provenance.json').read_text())
    assert provenance['counts']['input_reads']==15
    assert provenance['counts']['assigned_redi_reads']==7
    assert provenance['counts']['unassigned_reads']==8

def test_declared_empty_control_cannot_become_a_physical_clone(tmp_path):
    c,k,meta,key,p,rows=setup(tmp_path)
    meta.registration_status='declared_empty';meta.REDI_plate_barcode='empty';meta.MAT_alpha_array_1536_barcode=''
    meta.to_csv(k,sep='\t',index=False)
    out,_=run_step_c(c,k)
    assert pd.read_csv(out,sep='\t').empty
    provenance=json.loads(Path(str(out)+'.provenance.json').read_text())
    assert provenance['counts']['input_reads']==7 and provenance['counts']['unassigned_reads']==7
