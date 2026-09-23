"""Declared two-assay workflow oracle, independent of historical call tables.

The lightweight path supplies count tables with thirty molecules per assay.
The Snakefile path also supplies the full 150bp mates so its A-through-H fan-in
can be exercised with installed BBTools. No barcode is obtained from a caller.
"""
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest
import yaml

from magestic.pipelines.redi import _cli
from magestic.pipelines.redi.core.config import load_config


BC = 'ACGTATGCCATTCAGGTAAC'
OTHER = 'GATTCCGATACGTTCAGTAC'
RB = 'ACCCATAGTCCATTCGGCTA'


def declared_workflow(tmp_path):
    root = tmp_path/'screen with spaces'
    keys = root/'keys'; raw = root/'raw'
    counts = root/'processed/read_pair_counts_tables'
    for path in (keys,raw,counts):path.mkdir(parents=True)
    samples=[]; assays=[]; plates=[]; barcode_keys=[]
    rc=lambda value:value.translate(str.maketrans('ACGT','TGCA'))[::-1]
    molecule=('AGCTCAGGTCTAGTCGATCTAGGCTA'+'CAGGTCATGCTC'+BC+'CCCTAGGATACG'
              +'TGCATAACTTCGTATAATGTATGCTATACGAAGTTATGATC'+rc(RB)
              +'ACATACTTTTACTCCACGCACGGATCATCTGGCTA')
    for library,primer,plate in [('yT182','plate_4A','3567'),('yT183','plate_4B','3568')]:
        gdna='declared_'+library
        samples.append(dict(gDNA_plate_name=gdna,sample_number='sample_1',REDI_plate_barcode=plate))
        assays.append(dict(gDNA_plate_name=gdna,outer_primers=primer,inner_primers='KR1967-KR2225',
            purpose=library+' REDI',REDI_library=library,MAT_a_barcode_prefix='AAGTATGT',MAT_a_barcode_length='20'))
        plates.append(dict(REDI_plate_barcode=plate,MAT_alpha_array_1536_barcode='6156'))
        barcode_keys.append(dict(REDI_library=library,REDI_bc1_seq=RB,plate_384='RED_384',row_384='A',col_384='1'))
        name='20250627_'+primer+'_sample_1'
        for mate in (1,2):
            sequence=molecule[:150] if mate==1 else rc(molecule[-150:])
            with gzip.open(raw/f'{name}_R{mate}.fastq.gz','wt') as stream:
                for index in range(30):
                    stream.write(f'@declared_{index}/{mate}\n{sequence}\n+\n'+len(sequence)*'I'+'\n')
        pd.DataFrame([dict(sample_name=name,R1_seq='CAGGTCATGCTC'+BC+'CCCTAGGATACG',
            R2_seq='AAGTATGT'+RB+'GATCATAACTTC',counts=30)]).to_csv(
                counts/f'{name}_R1_R2_counts.tsv',sep='\t',index=False)
    for filename,rows in [('sample_key.tsv',samples),('gDNA_plate_key.tsv',assays),
                          ('REDI_plate_barcode_key.tsv',plates),('positions.tsv',barcode_keys)]:
        pd.DataFrame(rows).to_csv(keys/filename,sep='\t',index=False)
    pd.DataFrame([dict(bc1=BC,design='independently_declared_design')]).to_csv(root/'reference.tsv',sep='\t',index=False)
    pd.DataFrame([dict(sample='registered_wgs_001',barcode=BC)]).to_csv(root/'wgs.tsv',sep='\t',index=False)
    pd.DataFrame([dict(sample_name='registered_wgs_001',MAT_alpha_array_1536_barcode='6156',
        REDI_plate_color_384='RED_384',REDI_row_384='A',REDI_col_384='1')]).to_csv(keys/'WGS_plate_key.tsv',sep='\t',index=False)
    cfg=dict(screen_name='declared_workflow',sequencing_date='20250627',
        paths=dict(project_root='screen with spaces',raw_fastq='raw',intermediate='intermediate',processed='processed',keyfiles='keys'),
        redi_bc_keys=['positions.tsv'],bc1_reference_table='reference.tsv',wgs_outcome_table='wgs.tsv',
        array_pairs=[dict(name_a='old',name_b='new',purpose_a='yT182 REDI',purpose_b='yT183 REDI',require_validation_a=False)])
    path=tmp_path/'config with spaces.yaml';path.write_text(yaml.safe_dump(cfg))
    (root/'declared_oracle.json').write_text(json.dumps(dict(bc1=BC,REDI_bc=RB,assays=2,
        read_pairs_per_assay=30,total_read_pairs=60,registered_positions=1,WGS_positions=1,
        expected_wgs_agreements=1,independent_replicate_gate=False,molecule=molecule),indent=2)+'\n')
    return path,load_config(path)


def assert_workflow_outputs(cfg, *, expected_agreement=True):
    combined=cfg.paths.combined_data
    c=pd.read_csv(combined/'combined_REDI_bc_MAGESTIC_bc_merged_df.tsv',sep='\t')
    assert len(c)==2 and c.counts.sum()==60 and set(c.bc1)=={BC}
    f=pd.read_csv(combined/'old_new_merged_REDI_df.tsv',sep='\t')
    assert len(f)==1 and int(f.old_REDI_bc_bc1_counts.iloc[0])==30 and int(f.new_REDI_bc_bc1_counts.iloc[0])==30
    g=pd.read_csv(combined/'WGS_bc1_call_outer_merge_with_REDI_df_and_plate_key.tsv',sep='\t')
    assert len(g)==1
    for name in ('old','new'):
        assert bool(g[f'wgs_and_{name}_callable'].iloc[0])
        assert bool(g[f'wgs_and_{name}_bc1_agreement'].iloc[0]) is expected_agreement
    positions=cfg.paths.plots/'bc1_registered_position_counts.tsv'
    h=pd.read_csv(positions,sep='\t')
    assert h.to_dict(orient='records')==[dict(bc1=BC,registered_positions=1)]
    proof=json.loads(Path(str(positions)+'.provenance.json').read_text())
    assert proof['parameters']['replicate_policy']=='no independent-replicate gate'
    assert proof['counts']['registered_positions']==1 and proof['counts']['passing_positions']==1
    for artifact in proof['outputs']:
        path=Path(artifact['path'])
        assert path.is_file() and hashlib.md5(path.read_bytes()).hexdigest()==artifact['md5']
        if path.suffix=='.png':assert path.read_bytes().startswith(b'\x89PNG\r\n\x1a\n')


def test_real_cli_b_through_h_matches_declared_two_assay_oracle(tmp_path):
    path,cfg=declared_workflow(tmp_path)
    assert _cli.umbrella_main(['--config',str(path)])==0
    assert_workflow_outputs(cfg)


def test_cli_restart_g_h_recomputes_changed_wgs_without_changing_positions(tmp_path):
    path,cfg=declared_workflow(tmp_path)
    assert _cli.umbrella_main(['--config',str(path)])==0
    pd.DataFrame([dict(sample='registered_wgs_001',barcode=OTHER)]).to_csv(cfg.wgs_outcome_table,sep='\t',index=False)
    assert _cli.umbrella_main(['--config',str(path),'--from','g','--to','h'])==0
    assert_workflow_outputs(cfg,expected_agreement=False)


def _snake(tmp_path,config_path,*options):
    spec=importlib.util.find_spec('snakemake')
    if spec is None:pytest.skip('Snakemake is an optional workflow dependency; install it to run this parity acceptance test')
    import magestic.pipelines.redi as redi
    snakefile=Path(redi.__file__).parent/'snakemake/Snakefile'
    env=os.environ.copy()
    env['PYTHONPATH']=os.pathsep.join(str(p) for p in sys.path if p)
    env['XDG_CACHE_HOME']=str(tmp_path/'cache')
    env['MPLCONFIGDIR']=str(tmp_path/'mpl')
    env['REDI_TEST_SNAKE_CACHE']=str(tmp_path/'snake-cache')
    # Snakemake7's platformdirs ignores XDG_CACHE_HOME on macOS. Redirect only
    # its cache location; the workflow and scheduler remain the real module.
    bootstrap='''import os, runpy, sys
import snakemake.common
original = snakemake.common.get_appdirs()
class Cache:
    user_cache_dir = os.environ['REDI_TEST_SNAKE_CACHE']
    def __getattr__(self, name):
        return getattr(original, name)
snakemake.common.APPDIRS = Cache()
sys.argv[0] = 'snakemake'
runpy.run_module('snakemake', run_name='__main__', alter_sys=True)
'''
    command=[sys.executable,'-c',bootstrap,'--snakefile',str(snakefile),'--cores','1',
             '--config','config='+str(config_path),*options]
    result=subprocess.run(command,cwd=tmp_path,env=env,text=True,capture_output=True,timeout=300)
    with (tmp_path/'snakemake.log').open('a') as log:
        log.write(json.dumps(command)+'\n'+result.stdout+'\n'+result.stderr+'\n')
    assert result.returncode==0,result.stdout+'\n'+result.stderr
    return result


def test_snakefile_dag_includes_both_registered_assays_and_g_h(tmp_path):
    path,cfg=declared_workflow(tmp_path)
    result=_snake(tmp_path,path,'--dag')
    assert 'step_g_wgs_check' in result.stdout and 'step_h_plot' in result.stdout
    assert 'step_a_one_sample' in result.stdout
    assert '20250627_plate_4A_sample_1' in result.stdout and '20250627_plate_4B_sample_1' in result.stdout


def test_snakefile_outputs_match_same_declared_oracle(tmp_path):
    path,cfg=declared_workflow(tmp_path)
    _snake(tmp_path,path)
    assert_workflow_outputs(cfg)
    before={table:table.stat().st_mtime_ns for table in cfg.paths.counts_tables.glob('*.tsv')}
    assert len(before)==2
    _snake(tmp_path,path)
    assert {table:table.stat().st_mtime_ns for table in before}==before, 'unchanged native A work must remain reusable'
    pd.DataFrame([dict(sample='registered_wgs_001',barcode=OTHER)]).to_csv(cfg.wgs_outcome_table,sep='\t',index=False)
    _snake(tmp_path,path)
    assert_workflow_outputs(cfg,expected_agreement=False)
    assert {table:table.stat().st_mtime_ns for table in before}==before, 'WGS-only change must not repeat native A'
    image=cfg.paths.plots/'bc1_clones_distribution.png'
    image.unlink()
    _snake(tmp_path,path)
    assert_workflow_outputs(cfg,expected_agreement=False)
    assert {table:table.stat().st_mtime_ns for table in before}==before
    # A new completion guard is exercised against content tampering, not merely
    # an existence check. Restore the real PNG for subsequent visual inspection.
    from magestic.pipelines.redi.pipeline.run_pipeline import verify_completed_artifacts
    primary=cfg.paths.plots/'bc1_registered_position_counts.tsv'
    content=image.read_bytes();image.write_bytes(b'corrupted-output')
    try:
        with pytest.raises(ValueError,match='content changed'):
            verify_completed_artifacts([primary])
    finally:
        image.write_bytes(content)
