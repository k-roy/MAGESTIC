"""Primary NNS registration shape and lossless standalone-reader expectations.

These tests are a separate before-fix batch. They never infer a sample index
from row order; the positive geometry fixture supplies its index explicitly.
"""
import pandas as pd
import pytest
from magestic.pipelines.redi.readers.wgs_outcome import load_wgs_outcome,load_wgs_plate_key
from magestic.pipelines.redi.core.correspondence import load_type_i_wgs,join_correspondence


def actual_shape(tmp_path, explicit=True):
    key=pd.DataFrame({'WGS gDNA plate name':['20250227_plate_barcode_2461','20250227_plate_barcode_2462'],
        'Nextera outer primer set':['Tn5_plate_5C','Tn5_plate_6C'],
        'plate_96':['RED_96','BLUE_96'],'MAT_alpha_array_plate_384':['RED_384','RED_384'],
        'MAT_alpha_array_plate_384_barcode':['4994','4994'],
        'MAT_alpha_array_1536_barcode':['4996','4996'],'sequencing_date':['20250331','20250331']})
    fmt=pd.DataFrame({'row_384':['A','A'],'col_384':['1','2'],
                      'plate_96':['RED_96','BLUE_96'],'row_96':['A','A'],'col_96':['1','1']})
    if explicit:fmt['sample_number']=['sample_1','sample_1']
    kp=tmp_path/'key.tsv';fp=tmp_path/'geometry.tsv'
    key.to_csv(kp,sep='\t',index=False);fmt.to_csv(fp,sep='\t',index=False)
    return kp,fp


def test_actual_nns_schema_joins_by_96_quadrant_not_384_color(tmp_path):
    result=load_wgs_plate_key(*actual_shape(tmp_path))
    assert dict(zip(result.sample_name,result.col_384.astype(str)))=={
        '20250331_Tn5_plate_5C_sample_1':'1','20250331_Tn5_plate_6C_sample_1':'2'}
    assert result.MAT_alpha_array_plate_384.tolist()==['RED_384','RED_384']


def test_actual_nns_schema_without_explicit_sample_index_names_the_missing_artifact(tmp_path):
    with pytest.raises(ValueError,match='sample_number'):
        load_wgs_plate_key(*actual_shape(tmp_path,explicit=False))


def test_standalone_wgs_outcome_reader_preserves_zero_prefixed_sample_id(tmp_path):
    path=tmp_path/'wgs.tsv'
    path.write_text('sample_name\tbc1_from_wgs\n001\t'+'A'*20+'\n')
    assert load_wgs_outcome(path).sample_name.tolist()==['001']


def test_type_i_reader_preserves_zero_prefixed_sample_and_position_ids(tmp_path):
    path=tmp_path/'wgs.tsv'
    path.write_text('sample_name\tplate_and_position\tbc1_from_wgs\n001\t002\t'+'A'*20+'\n')
    result=load_type_i_wgs(str(path))
    assert result.sample_name.tolist()==['001'] and result.plate_and_position.tolist()==['002']


def test_type_i_reader_refuses_observation_rename_collision(tmp_path):
    path=tmp_path/'wgs.tsv'
    path.write_text('sample_name\tbc1\tbc1_from_wgs\n001\t'+'A'*20+'\t'+'C'*20+'\n')
    with pytest.raises(ValueError,match='duplicate|collision|reserved'):
        load_type_i_wgs(str(path))


def test_correspondence_refuses_reserved_internal_presence_column():
    left=pd.DataFrame({'plate_and_position':['p'],'bc1':['A'*20]})
    right=pd.DataFrame({'plate_and_position':['p'],'MAGESTIC_bc1':['A'*20],'_in_wgs':[False]})
    with pytest.raises(ValueError,match='reserved|collision'):
        join_correspondence(left,right)
