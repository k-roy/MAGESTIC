"""Stage A contracts with controlled merger output; actual BBTools smoke separate."""
import gzip,json
from pathlib import Path
import pandas as pd
import pytest
from magestic.pipelines.redi.core.config import REDIConfig,REDIPaths
from magestic.pipelines.redi.pipeline import step_a_trim_counts as stage

def setup(tmp_path,monkeypatch):
    paths=REDIPaths(tmp_path,tmp_path/'raw',tmp_path/'intermediate',tmp_path/'processed',tmp_path/'keys')
    paths.raw_fastq.mkdir();paths.counts_tables.mkdir(parents=True)
    cfg=REDIConfig('fixture','20250627',paths)
    prefix='20250627_plate_4A_sample_1';molecule='CAGGTCATGCTC'+'ACGTATGCCATTCAGGTAAC'+'CCCTAGGATACG'
    for mate in (1,2):
        with gzip.open(paths.raw_fastq/f'{prefix}_R{mate}.fastq.gz','wt') as f:
            for i in range(3):f.write(f'@read{i}/{mate}\n{molecule}\n+\n'+len(molecule)*'I'+'\n')
    calls=[]
    def trim(a,b,c,d,**kwargs):
        calls.append(('trim',kwargs))
        for source,target in ((a,c),(b,d)):
            with gzip.open(source,'rt') as f:Path(target).write_text(f.read())
    def merge(a,b,out,**kwargs):
        calls.append(('merge',kwargs));Path(out).write_text(Path(a).read_text())
    monkeypatch.setattr(stage,'_bbduk',trim);monkeypatch.setattr(stage,'_bbmerge',merge)
    return cfg,prefix,molecule,calls

def test_default_pipeline_has_native_counts_and_explicit_merged_representation(tmp_path,monkeypatch):
    cfg,prefix,molecule,calls=setup(tmp_path,monkeypatch)
    stage.run_step_a_one_sample(prefix,cfg)
    out=cfg.paths.counts_tables/f'{prefix}_R1_R2_counts.tsv';rows=pd.read_csv(out,sep='\t')
    assert len(rows)==1 and rows.counts.iloc[0]==3
    assert rows.R1_seq.iloc[0]==molecule
    assert rows.R2_seq.iloc[0]==molecule.translate(str.maketrans('ACGT','TGCA'))[::-1]
    meta=json.loads(Path(str(out)+'.provenance.json').read_text())
    assert meta['counts']['raw_read_pairs']==3 and meta['counts']['retained_merged_reads']==3
    assert meta['parameters']['sequence_representation']=='merged_consensus_and_reverse_complement'
    assert all(c[1]['threads']==1 and c[1]['memory_mb']==512 for c in calls)

def test_mate_name_mismatch_is_rejected_before_tools(tmp_path,monkeypatch):
    cfg,prefix,molecule,calls=setup(tmp_path,monkeypatch)
    with gzip.open(cfg.paths.raw_fastq/f'{prefix}_R2.fastq.gz','wt') as f:f.write('@wrong/2\nACGT\n+\nIIII\n')
    with pytest.raises(ValueError,match='mate|pair'):stage.run_step_a_one_sample(prefix,cfg)
    assert calls==[]

def test_unequal_mate_counts_cannot_be_silently_zipped(tmp_path,monkeypatch):
    cfg,prefix,molecule,calls=setup(tmp_path,monkeypatch)
    with gzip.open(cfg.paths.raw_fastq/f'{prefix}_R2.fastq.gz','wt') as f:
        f.write('@read0/2\nACGT\n+\nIIII\n')
    with pytest.raises(ValueError,match='mate|pair'):stage.run_step_a_one_sample(prefix,cfg)
    assert calls==[]

def test_existing_unattested_output_is_not_a_successful_cache(tmp_path,monkeypatch):
    cfg,prefix,molecule,calls=setup(tmp_path,monkeypatch)
    (cfg.paths.counts_tables/f'{prefix}_R1_R2_counts.tsv').write_text('header_only\n')
    with pytest.raises(FileExistsError):stage.run_step_a_one_sample(prefix,cfg,overwrite=False)

def test_minimum_count_filter_has_its_own_denominator(tmp_path,monkeypatch):
    cfg,prefix,molecule,calls=setup(tmp_path,monkeypatch)
    stage.run_step_a_one_sample(prefix,cfg,min_counts=4)
    out=cfg.paths.counts_tables/f'{prefix}_R1_R2_counts.tsv'
    assert pd.read_csv(out,sep='\t').empty
    meta=json.loads(Path(str(out)+'.provenance.json').read_text())
    assert meta['counts']['merged_reads']==3 and meta['counts']['below_minimum_count_reads']==3

def test_sample_path_traversal_is_refused(tmp_path,monkeypatch):
    cfg,prefix,molecule,calls=setup(tmp_path,monkeypatch)
    with pytest.raises(ValueError,match='sample|prefix'):stage.run_step_a_one_sample('../escape',cfg)
    assert calls==[]

@pytest.mark.parametrize('kind',['truncated','quality_length','bad_sequence','bad_plus'])
def test_malformed_fastq_refused_before_tools(tmp_path,monkeypatch,kind):
    cfg,prefix,molecule,calls=setup(tmp_path,monkeypatch)
    text={'truncated':'@read0/1\nACGT\n+\n',
          'quality_length':'@read0/1\nACGT\n+\nIII\n',
          'bad_sequence':'@read0/1\nACXT\n+\nIIII\n',
          'bad_plus':'@read0/1\nACGT\nminus\nIIII\n'}[kind]
    with gzip.open(cfg.paths.raw_fastq/f'{prefix}_R1.fastq.gz','wt') as f:f.write(text)
    with pytest.raises(ValueError,match='FASTQ'):stage.run_step_a_one_sample(prefix,cfg)
    assert calls==[]


def test_empty_paired_input_is_measured_zero_with_no_tool_invocation(tmp_path,monkeypatch):
    cfg,prefix,molecule,calls=setup(tmp_path,monkeypatch)
    for mate in (1,2):
        with gzip.open(cfg.paths.raw_fastq/f'{prefix}_R{mate}.fastq.gz','wt'):pass
    stage.run_step_a_one_sample(prefix,cfg)
    out=cfg.paths.counts_tables/f'{prefix}_R1_R2_counts.tsv'
    assert calls==[] and pd.read_csv(out,sep='\t').empty
    meta=json.loads(Path(str(out)+'.provenance.json').read_text())
    assert meta['counts']['raw_read_pairs']==meta['counts']['retained_merged_reads']==0


def test_altered_output_cannot_be_reused_as_attested_cache(tmp_path,monkeypatch):
    cfg,prefix,molecule,calls=setup(tmp_path,monkeypatch)
    stage.run_step_a_one_sample(prefix,cfg)
    out=cfg.paths.counts_tables/f'{prefix}_R1_R2_counts.tsv';out.write_text(out.read_text()+'\n')
    with pytest.raises(FileExistsError):stage.run_step_a_one_sample(prefix,cfg,overwrite=False)
    assert len(calls)==2
