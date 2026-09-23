"""Stage D/E file publication preserves exact identities and count denominators."""
import json
from pathlib import Path
import pandas as pd
import pytest
from magestic.pipelines.redi.core.config import REDIConfig,REDIPaths
from magestic.pipelines.redi.core.annotate import annotate_with_bc1_reference
from magestic.pipelines.redi.pipeline.step_d_purity import run_step_d
from magestic.pipelines.redi.pipeline.step_e_annotate import run_step_e

BC='ACGTATGCCATTCAGGTAAC'
def setup(tmp_path):
 paths=REDIPaths(tmp_path,tmp_path/'raw',tmp_path/'int',tmp_path/'out',tmp_path/'keys')
 paths.combined_data.mkdir(parents=True)
 c=REDIConfig('fixture','20250627',paths)
 source=paths.combined_data/'counts.tsv'
 rows=pd.DataFrame([dict(sample_name='dated_001',REDI_bc='ACCT',bc1=BC,counts=9007199254740993,
   MAT_alpha_array_1536_barcode='00123',purpose='NA')])
 rows.to_csv(source,sep='\t',index=False)
 return c,source

def test_stage_d_attests_all_inputs_and_keeps_identifiers_exact(tmp_path):
 c,source=setup(tmp_path);out=run_step_d(c,source)
 frame=pd.read_csv(out,sep='\t',dtype=str,keep_default_na=False)
 assert frame.MAT_alpha_array_1536_barcode.iloc[0]=='00123'
 assert frame.purpose.iloc[0]=='NA'
 assert frame.counts.iloc[0]=='9007199254740993'
 meta=json.loads(Path(str(out)+'.provenance.json').read_text())
 assert meta['counts']['input_rows']==1 and meta['counts']['output_rows']==1
 assert meta['code']['redi_source_sha256'] and len(meta['outputs'])==2


def test_stage_e_preserves_identity_and_attests_both_input_sources(tmp_path):
 c,source=setup(tmp_path);top=run_step_d(c,source)
 ref=c.paths.project_root/'reference.tsv';pd.DataFrame([dict(bc1=BC,donor='NA',amino_acid_position=12)]).to_csv(ref,sep='\t',index=False)
 c.bc1_reference_table=ref;out=run_step_e(c,top)
 frame=pd.read_csv(out,sep='\t',dtype=str,keep_default_na=False)
 assert frame.MAT_alpha_array_1536_barcode.iloc[0]=='00123' and frame.donor.iloc[0]=='NA'
 assert frame.REDI_bc_bc1_counts.iloc[0]=='9007199254740993'
 meta=json.loads(Path(str(out)+'.provenance.json').read_text())
 assert len(meta['inputs'])==2 and meta['counts']['input_rows']==1 and meta['counts']['output_rows']==1

@pytest.mark.parametrize('column',['all_bc1_counts_json','observed_bc1_read_fraction','identity_ambiguous','colony_purity','MAT_alpha_array_1536_barcode','sample_name'])
def test_reference_cannot_replace_measured_distribution_or_physical_identity(column):
 top=pd.DataFrame([{'bc1':BC,column:'observed','counts':3,'total_counts':4}])
 ref=pd.DataFrame([{'bc1':BC,column:'reference'}])
 with pytest.raises(ValueError,match='observation|reserved|collision'):
  annotate_with_bc1_reference(top,ref)


def test_stage_b_completion_marker_attests_three_registration_sources(tmp_path):
 from magestic.pipelines.redi.pipeline.step_b_combine_keys import run_step_b
 c,_=setup(tmp_path);c.paths.keyfiles.mkdir()
 for name,rows in [
  (c.sample_key,[dict(gDNA_plate_name='dated',sample_number='sample_1',REDI_plate_barcode='03500')]),
  (c.gdna_plate_key,[dict(gDNA_plate_name='dated',outer_primers='plate_4A',inner_primers='KR1967',purpose='assay')]),
  (c.redi_plate_barcode_key,[dict(REDI_plate_barcode='03500',MAT_alpha_array_1536_barcode='06100')])]:
  pd.DataFrame(rows).to_csv(c.paths.keyfiles/name,sep='\t',index=False)
 out=run_step_b(c)
 assert pd.read_csv(out,sep='\t',dtype=str).MAT_alpha_array_1536_barcode.iloc[0]=='06100'
 meta=json.loads(Path(str(out)+'.provenance.json').read_text())
 assert len(meta['inputs'])==3 and meta['counts']['output_rows']==1
