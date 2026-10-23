"""Independent G–I/workflow/plate interface expectations, frozen before repair.

Raw plate headers below are primary archived workbook A1 cells, not invented
normalized fixtures (bc311e534a3ce6d42d6da8d561efa6ee37e8c41cbcf4727fd72eddbf616822cc).
Plot stubs capture plotted numerical data; they do not claim visual rendering QA.
"""
from pathlib import Path
from types import ModuleType
import sys

import pandas as pd
import pytest

from magestic.pipelines.redi.core.config import REDIArrayPair, REDIConfig, REDIPaths
from magestic.pipelines.redi.pipeline import step_g_wgs_check as stage_g
from magestic.pipelines.redi.pipeline import step_h_plot as stage_h
from magestic.pipelines.redi.pipeline import step_i_select as stage_i
from magestic.pipelines.redi.pipeline import run_pipeline
from magestic.pipelines.redi.readers import plate_db


A='A'*20
C='C'*20


def config(tmp_path, pairs=None):
    paths=REDIPaths(tmp_path,tmp_path/'raw',tmp_path/'intermediate',tmp_path/'processed',tmp_path/'keys')
    paths.keyfiles.mkdir();paths.combined_data.mkdir(parents=True)
    return REDIConfig('independent','20250627',paths,
                      array_pairs=pairs or [REDIArrayPair('old','new','earlier',purpose_b='later')])


def cross_rows(a='old',b='new',barcodes=(A,A)):
    return pd.DataFrame([dict(MAT_alpha_array_1536_barcode='6156',REDI_plate_color_384='RED',
        REDI_row_384='A',REDI_col_384=str(i+1),num_independent_replicates=2,**{
            a+'_bc1':bc,b+'_bc1':bc,a+'_REDI_bc_bc1_counts':30,b+'_REDI_bc_bc1_counts':30,
            a+'_colony_purity':1.0,b+'_colony_purity':1.0}) for i,bc in enumerate(barcodes)])


def plotting_stub(monkeypatch):
    captured=[]
    matplotlib=ModuleType('matplotlib');pyplot=ModuleType('matplotlib.pyplot');sns=ModuleType('seaborn')
    matplotlib.use=lambda *a,**k:None
    matplotlib.pyplot=pyplot
    for name in ('figure','xlabel','ylabel','title','xticks','tight_layout','close','legend'):
        setattr(pyplot,name,lambda *a,**k:None)
    pyplot.savefig=lambda path,**kw:Path(path).write_bytes(b'plot-data-test-only')
    sns.barplot=lambda *,x,y,**kw:captured.append((list(x),list(y)))
    for name,module in [('matplotlib',matplotlib),('matplotlib.pyplot',pyplot),('seaborn',sns)]:
        monkeypatch.setitem(sys.modules,name,module)
    return captured


def test_stage_h_keeps_two_physical_clones_of_one_barcode(tmp_path,monkeypatch):
    cfg=config(tmp_path);captured=plotting_stub(monkeypatch)
    cross_rows().to_csv(cfg.paths.combined_data/'old_new_merged_REDI_df.tsv',sep='\t',index=False)
    cfg.bc1_reference_table=tmp_path/'reference.tsv'
    pd.DataFrame({'bc1':[A],'design':['one']}).to_csv(cfg.bc1_reference_table,sep='\t',index=False)
    stage_h.run_step_h(cfg)
    assert captured==[([2],[1])], 'two colonies with one barcode are two clones, not one barcode record'


def test_stage_h_ambiguous_reference_cannot_multiply_clone_count(tmp_path,monkeypatch):
    cfg=config(tmp_path);captured=plotting_stub(monkeypatch)
    cross_rows(barcodes=(A,)).to_csv(cfg.paths.combined_data/'old_new_merged_REDI_df.tsv',sep='\t',index=False)
    cfg.bc1_reference_table=tmp_path/'reference.tsv'
    pd.DataFrame({'bc1':[A,A],'design':['one','two']}).to_csv(cfg.bc1_reference_table,sep='\t',index=False)
    stage_h.run_step_h(cfg)
    assert captured==[([1],[1])], 'two design candidates do not create a second physical clone'


def test_stage_h_without_reference_still_has_observed_barcode(tmp_path,monkeypatch):
    cfg=config(tmp_path);captured=plotting_stub(monkeypatch)
    cross_rows().to_csv(cfg.paths.combined_data/'old_new_merged_REDI_df.tsv',sep='\t',index=False)
    stage_h.run_step_h(cfg)
    assert captured==[([2],[1])]


def test_stage_h_processes_every_configured_array_pair(tmp_path,monkeypatch):
    pairs=[REDIArrayPair('old','new','a',purpose_b='b'),REDIArrayPair('secondA','secondB','c',purpose_b='d')]
    cfg=config(tmp_path,pairs);captured=plotting_stub(monkeypatch)
    for pair in pairs:
        cross_rows(pair.name_a,pair.name_b).assign(bc1=A).to_csv(
            cfg.paths.combined_data/f'{pair.name_a}_{pair.name_b}_merged_REDI_df.tsv',sep='\t',index=False)
    outputs=stage_h.run_step_h(cfg)
    assert len(captured)==2
    assert len({Path(p) for p in outputs.values() if Path(p).suffix=='.png'})==2


def test_high_quality_filter_cannot_ignore_requested_replicate_gate():
    with pytest.raises(ValueError,match='replicate'):
        stage_h.filter_high_quality_clones(cross_rows().drop(columns='num_independent_replicates'),'old','new',min_replicates=2)


def test_high_quality_filter_accepts_literal_array_names():
    frame=cross_rows('old-era','repeat era')
    assert len(stage_h.filter_high_quality_clones(frame,'old-era','repeat era'))==2


def test_unique_identity_falls_back_from_empty_b_to_valid_a():
    frame=cross_rows(barcodes=(A,));frame['new_bc1']=''
    assert stage_h.unique_bc1_by_mean_purity(frame,'old','new').bc1.tolist()==[A]


def test_conflicting_valid_arrays_cannot_select_b_as_the_identity():
    frame=cross_rows(barcodes=(A,));frame['new_bc1']=C
    result=stage_h.unique_bc1_by_mean_purity(frame,'old','new')
    assert result.empty or result.bc1.isna().all(), 'conflicting observations remain unassessed'


@pytest.mark.parametrize('n',[-1,1.5,True])
def test_selection_rejects_invalid_n(n):
    with pytest.raises(ValueError,match='n|integer|positive'):
        stage_i.select_top_n_per_group(pd.DataFrame({'gene':['g']*3,'mean_purity':[.9,.8,.7]}),['gene'],n=n)


def test_selection_missing_group_is_not_silently_dropped():
    with pytest.raises(ValueError,match='group|missing'):
        stage_i.select_top_n_per_group(pd.DataFrame({'gene':['g',None],'mean_purity':[.9,1.0]}),['gene'])


def test_reversed_pipeline_range_cannot_succeed_without_work(tmp_path):
    with pytest.raises(ValueError,match='order|range|before|from'):
        run_pipeline.run(config(tmp_path),from_step='h',to_step='b')


def test_stage_g_input_mutation_cannot_attest_old_observations(tmp_path,monkeypatch):
    cfg=config(tmp_path)
    cross=cross_rows(barcodes=(A,))
    cross.to_csv(cfg.paths.combined_data/'old_new_merged_REDI_df.tsv',sep='\t',index=False)
    positions=['MAT_alpha_array_1536_barcode','REDI_plate_color_384','REDI_row_384','REDI_col_384']
    key=cross[positions].assign(sample_name='registered')
    key.to_csv(cfg.paths.keyfiles/'WGS_plate_key.tsv',sep='\t',index=False)
    cfg.wgs_outcome_table=tmp_path/'calls.tsv'
    pd.DataFrame({'sample':['registered'],'barcode':[A]}).to_csv(cfg.wgs_outcome_table,sep='\t',index=False)
    original=stage_g.load_wgs_bc1_calls
    def changed_after_read(path,**kwargs):
        result=original(path,**kwargs)
        Path(path).write_text(Path(path).read_text()+'\n')
        return result
    monkeypatch.setattr(stage_g,'load_wgs_bc1_calls',changed_after_read)
    with pytest.raises(ValueError,match='input changed'):
        stage_g.run_step_g(cfg)
    assert not list(cfg.paths.combined_data.glob('*WGS*.provenance.json'))


@pytest.mark.parametrize('stage,column,identifier',[
    ('Source','rectangular plate barcode labels for PIXL-Spinnaker',1000),
    ('Target','target_plate_ID',2000),
    ('REDI_mating','barcode',3000),
    ('Misc','PlusPlate barcode label',4000),
    ('Rearray_Target','384 YPD+G418 agar plate barcode',5000),
    ('Consolidated_384','PlusPlate barcode label',6000),
])
def test_plate_reader_accepts_actual_snapshot_primary_headers(stage,column,identifier):
    db={stage:pd.DataFrame({column:[str(identifier)+'.0']})}
    hit=plate_db.resolve_plate(identifier,db)
    assert hit is not None and hit[0]==stage


def test_plate_lookup_refuses_duplicate_identity_even_across_stages():
    db={'Target':pd.DataFrame({'plate_id':[2000]}),'Consolidated_384':pd.DataFrame({'plate_id':[2000]})}
    with pytest.raises(ValueError,match='duplicate|ambiguous'):
        plate_db.resolve_plate(2000,db)


def test_plate_lookup_refuses_duplicate_identity_in_one_sheet():
    with pytest.raises(ValueError,match='duplicate|ambiguous'):
        plate_db.resolve_plate(2000,{'Target':pd.DataFrame({'plate_id':[2000,2000]})})


def test_lineage_does_not_truncate_fractional_parent_identifier():
    db={'Target':pd.DataFrame({'plate_id':[2000],'source_plate_id':[1000.5]})}
    with pytest.raises(ValueError,match='integral|integer|plate'):
        plate_db.trace_lineage(2000,db)


def test_cyclic_lineage_is_not_reported_as_a_finished_chain():
    db={'Target':pd.DataFrame({'plate_id':[2000,2001],'source_plate_id':[2001,2000]})}
    with pytest.raises(ValueError,match='cycle'):
        plate_db.trace_lineage(2000,db)


def test_lineage_does_not_fabricate_a_missing_parent_record():
    db={'Target':pd.DataFrame({'plate_id':[2000],'source_plate_id':[1000]})}
    with pytest.raises(ValueError,match='missing|unresolved|found'):
        plate_db.trace_lineage(2000,db)
