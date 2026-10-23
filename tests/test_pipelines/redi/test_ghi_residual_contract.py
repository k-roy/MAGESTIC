"""Further GHI counterexamples frozen separately after initial repairs."""
import importlib.util
from pathlib import Path

import pandas as pd
import pytest

from magestic.pipelines.redi.pipeline import step_h_plot as stage_h
from magestic.pipelines.redi.readers.plate_db import load_plate_db, resolve_plate
from test_ghi_followthrough_contract import config, cross_rows, plotting_stub


def test_plot_quality_must_be_supported_by_the_selected_callable_identity(tmp_path, monkeypatch):
    cfg=config(tmp_path);captured=plotting_stub(monkeypatch)
    frame=cross_rows(barcodes=('A'*20,))
    frame['old_REDI_bc_bc1_counts']=10
    # A is the callable identity but below the explicit20-read plot floor.
    # B's100 reads describe an uncalled/off-length winner, never support for A.
    frame['new_bc1']='';frame['new_REDI_bc_bc1_counts']=100
    frame.to_csv(cfg.paths.combined_data/'old_new_merged_REDI_df.tsv',sep='\t',index=False)
    outputs=stage_h.run_step_h(cfg)
    result=pd.read_csv(outputs['registered_position_counts'],sep='\t')
    assert len(result)==0 and captured==[]


def exporter():
    # The project-level snapshot exporter is intentionally outside src/.
    source=Path(__file__).resolve().parents[3]/'REDI/plate_database/export_plate_db.py'
    spec=importlib.util.spec_from_file_location('redi_plate_export_under_test',source)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


SHEETS=[('Source 1000-1999','rectangular plate barcode labels for PIXL-Spinnaker',1001),
        ('Target 2000-2999','target_plate_ID',2001),
        ('REDI_mating 3000-3999','barcode',3001),
        ('misc 4000-4999','PlusPlate barcode label',4001),
        ('Rearray_Target 5000-5999','384 YPD+G418 agar plate barcode',5001),
        ('Consolidated_384','PlusPlate barcode label',6001),
        ('misc_7000-7999','PlusPlate_barcode_label',7001)]


def workbook(path, sheets):
    with pd.ExcelWriter(path,engine='openpyxl') as writer:
        for name,header,identifier in sheets:
            pd.DataFrame({header:[identifier]}).to_excel(writer,sheet_name=name,index=False)


@pytest.mark.parametrize('reverse',[False,True])
def test_export_preserves_both_misc_eras_independent_of_sheet_order(tmp_path,reverse):
    source=tmp_path/'primary.xlsx';workbook(source,list(reversed(SHEETS)) if reverse else SHEETS)
    output=tmp_path/'snapshot';result=exporter().export(source,output)
    assert set(result)-{'__union_tsv'}=={'Source','Target','REDI_mating','Misc','Rearray_Target','Consolidated_384','Misc_7000'}
    snapshot=load_plate_db(output/'plates')
    assert resolve_plate(4001,snapshot)[0]=='Misc'
    assert resolve_plate(7001,snapshot)[0]=='Misc_7000'
    assert len(pd.read_csv(output/'plates.tsv',sep='\t'))==7


def test_export_refuses_ambiguous_primary_sheet_before_publishing(tmp_path):
    source=tmp_path/'ambiguous.xlsx'
    workbook(source,SHEETS+[('misc 4000 archived','PlusPlate barcode label',4002)])
    output=tmp_path/'snapshot'
    with pytest.raises(ValueError,match='ambiguous|multiple|unique'):
        exporter().export(source,output)
    assert not list(output.rglob('*.parquet'))


def test_export_refuses_missing_required_primary_sheet_before_publishing(tmp_path):
    source=tmp_path/'missing.xlsx';workbook(source,SHEETS[1:]);output=tmp_path/'snapshot'
    with pytest.raises(ValueError,match='missing|Source|required'):
        exporter().export(source,output)
    assert not list(output.rglob('*.parquet'))
