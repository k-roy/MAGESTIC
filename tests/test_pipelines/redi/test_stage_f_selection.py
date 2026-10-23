"""Cross-array selection cannot silently absorb a third assay or era."""
import json
import pandas as pd
import pytest
from magestic.pipelines.redi.core.config import REDIArrayPair,REDIConfig,REDIPaths
from magestic.pipelines.redi.pipeline.step_f_cross_array import run_step_f

def setup(tmp_path):
    paths=REDIPaths(tmp_path,tmp_path/'raw',tmp_path/'intermediate',tmp_path/'processed',tmp_path/'keys')
    paths.keyfiles.mkdir();paths.combined_data.mkdir(parents=True)
    # These fixtures test assay selection, not prior physical-well validation.
    config=REDIConfig('fixture','20250627',paths,array_pairs=[REDIArrayPair('old','new','earlier',require_validation_a=False)])
    rows=[]
    for purpose,bc,n in [('earlier','A'*20,20),('later','A'*20,20),('unrelated_era','T'*20,100)]:
        rows.append(dict(sample_name=purpose,REDI_bc='C'*20,purpose=purpose,
            MAT_alpha_array_1536_barcode='6156',REDI_plate_color_384='RED_384',REDI_row_384='A',REDI_col_384=1,
            bc1=bc,REDI_bc_bc1_counts=n,REDI_bc_bc1_total_counts=n,all_bc1_counts_json=json.dumps({bc:n})))
    data=pd.DataFrame(rows);path=paths.combined_data/'annotated.tsv';data.to_csv(path,sep='\t',index=False)
    return config,path,data

def test_three_purposes_require_explicit_second_selector(tmp_path):
    config,path,_=setup(tmp_path)
    with pytest.raises(ValueError,match='purpose|selector|ambiguous'):run_step_f(config,path)

def test_explicit_second_selector_excludes_unrelated_era(tmp_path):
    config,path,_=setup(tmp_path);config.array_pairs[0].purpose_b='later'
    result=pd.read_csv(run_step_f(config,path)[0],sep='\t')
    assert len(result)==1 and bool(result.bc1_agreement.iloc[0])
    assert result.new_REDI_bc_bc1_counts.iloc[0]==20

def test_selector_cannot_validate_same_observation_against_itself(tmp_path):
    config,path,_=setup(tmp_path);config.array_pairs[0].purpose_b='earlier'
    with pytest.raises(ValueError,match='overlap|same|distinct'):run_step_f(config,path)

def test_geometry_duplicates_cannot_multiply_positions(tmp_path):
    config,path,data=setup(tmp_path);config.array_pairs[0].purpose_b='later'
    data.iloc[:2].to_csv(path,sep='\t',index=False)
    fmt=pd.DataFrame([dict(plate_384='RED_384',row_384='A',col_384=1,row_1536='A',col_1536=c) for c in (1,2)])
    fmt.to_csv(config.paths.keyfiles/'1536_PIXL_rearray_format.tsv',sep='\t',index=False)
    with pytest.raises(ValueError,match='duplicate|unique'):run_step_f(config,path)

def test_configured_barcode_length_reaches_comparison(tmp_path):
    config,path,data=setup(tmp_path);config.barcode.bc1_length=19
    data=data.iloc[:2].copy();data.bc1='A'*19;data.all_bc1_counts_json=json.dumps({'A'*19:20})
    data.to_csv(path,sep='\t',index=False)
    result=pd.read_csv(run_step_f(config,path)[0],sep='\t')
    assert result.old_bc1.iloc[0]=='A'*19 and bool(result.bc1_agreement.iloc[0])

def test_input_changed_after_read_cannot_receive_fresh_attestation(tmp_path,monkeypatch):
    from magestic.pipelines.redi.pipeline import step_f_cross_array as stage
    config,path,_=setup(tmp_path);config.array_pairs[0].purpose_b='later'
    original=stage._read_key
    def changing(p):
        frame=original(p)
        if p==path:path.write_text(path.read_text()+'\n')
        return frame
    monkeypatch.setattr(stage,'_read_key',changing)
    with pytest.raises(ValueError,match='input changed'):stage.run_step_f(config,path)
    assert not (config.paths.combined_data/'old_new_merged_REDI_df.tsv.provenance.json').exists()
