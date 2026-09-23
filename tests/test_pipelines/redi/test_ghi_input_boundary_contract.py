"""No physical observation row may disappear through a permissive TSV parse."""
import pandas as pd
import pytest

from magestic.pipelines.redi.pipeline import step_g_wgs_check as stage_g
from magestic.pipelines.redi.pipeline import step_h_plot as stage_h
from magestic.pipelines.redi.core.identity import POSITION_COLUMNS
from test_ghi_followthrough_contract import config, cross_rows, plotting_stub


def malformed_cross(cfg):
    frame=cross_rows()
    path=cfg.paths.combined_data/'old_new_merged_REDI_df.tsv'
    lines=frame.to_csv(sep='\t',index=False).splitlines()
    path.write_text('\n'.join(lines[:2]+['']+lines[2:])+'\n')
    return frame,path


def test_stage_h_refuses_blank_physical_row_in_observation_table(tmp_path,monkeypatch):
    cfg=config(tmp_path);plotting_stub(monkeypatch);malformed_cross(cfg)
    with pytest.raises(ValueError,match='blank|row'):
        stage_h.run_step_h(cfg)
    assert not list(cfg.paths.plots.glob('*.provenance.json'))


def test_stage_g_refuses_blank_physical_row_in_observation_table(tmp_path):
    cfg=config(tmp_path);frame,path=malformed_cross(cfg)
    key=frame[list(POSITION_COLUMNS)].assign(sample_name=['registered_1','registered_2'])
    key.to_csv(cfg.paths.keyfiles/'WGS_plate_key.tsv',sep='\t',index=False)
    cfg.wgs_outcome_table=tmp_path/'calls.tsv'
    pd.DataFrame({'sample':key.sample_name,'barcode':['A'*20]*2}).to_csv(cfg.wgs_outcome_table,sep='\t',index=False)
    with pytest.raises(ValueError,match='blank|row'):
        stage_g.run_step_g(cfg)
    assert not list(cfg.paths.combined_data.glob('*WGS*.provenance.json'))
