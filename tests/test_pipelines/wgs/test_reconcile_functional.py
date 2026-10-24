"""Portable reconciliation molecular fixtures and exact expected outcomes.

All BAMs/FASTA inputs are generated here from independently specified molecules.
Requires the WGS dependencies and samtools for donor-cleaning cases; no external datasets or audit-workspace imports.
"""
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import pysam
import pytest

from trace_crispr.wgs import reconcile as r

from dataclasses import replace
import shutil
import yaml
from magestic.pipelines.wgs.core import reconcile_driver as driver
from magestic.pipelines.wgs.core.locus_config import build_tracks
from magestic.pipelines.wgs.pipeline import step_1_bc1_call as s1, step_4_reconcile as s4
from magestic.pipelines.wgs.pipeline._common import StepIncompleteError

BC = 'ACGTCAGTTCGATGACCTGA'
GLEFT, GRIGHT = 'ACCCTGTTAT', 'GAGCATGACC'
PLEFT, PRIGHT = 'GGTCATGCTC', 'ATAACAGGGT'

def make_read(header, name, contig, start, seq, cigar, flag=99, mate=None, mate_start=180, tlen=160, mapq=60, qual=35):
    r = pysam.AlignedSegment(header)
    r.query_name = name
    r.query_sequence = seq
    r.flag = flag
    r.reference_id = header.get_tid(contig)
    r.reference_start = start
    r.cigarstring = cigar
    r.mapping_quality = mapq
    r.next_reference_id = header.get_tid(mate or contig)
    r.next_reference_start = mate_start
    r.template_length = tlen
    r.query_qualities = pysam.qualitystring_to_array(chr(qual+33)*len(seq))
    return r

def write_bam(path, references, records):
    header_dict = pysam.AlignmentHeader.from_references(list(references),[len(v) if isinstance(v,str) else v for v in references.values()]).to_dict()
    header_dict['PG'] = [{'ID':'fixture_align', 'CL':f'fixture_align {path.parent / "ref.fa"}'}]
    header = pysam.AlignmentHeader.from_dict(header_dict)
    reads = [make_read(header, **d) for d in records]
    with pysam.AlignmentFile(str(path),'wb',header=header) as f:
        for r in sorted(reads,key=lambda r:(r.reference_id,r.reference_start)):
            f.write(r)
    pysam.index(str(path))
    return path

@pytest.fixture
def fixture_factory(tmp_path):
    def make(*, records=None, plasmid=True, bridge=True, window=(301,310), extra_refs=None):
        from magestic.pipelines.wgs.config import load_config
        root=tmp_path/('case'+str(len(list(tmp_path.iterdir()))));root.mkdir()
        refs={'chrII':'A'*90+GLEFT+'N'*20+GRIGHT+'A'*370}
        if plasmid:refs['plasmid_demo']='A'*90+PLEFT+'N'*20+PRIGHT+'A'*370
        refs.update(extra_refs or {})
        fa=root/'ref.fa';fa.write_text(''.join('>'+c+'\n'+s+'\n' for c,s in refs.items()));pysam.faidx(str(fa))
        if records is None:
            from trace_crispr.wgs.bc1_consensus import revcomp
            records=[dict(name='g1',contig='chrII',start=90,seq=GLEFT+revcomp(BC),cigar='10M20S'),dict(name='g2',contig='chrII',start=89,seq='A'+GLEFT+revcomp(BC),cigar='11M20S',tlen=161)]
            if plasmid:records += [dict(name='p1',contig='plasmid_demo',start=90,seq=PLEFT+BC,cigar='10M20S'),dict(name='p2',contig='plasmid_demo',start=89,seq='A'+PLEFT+BC,cigar='11M20S',tlen=161)]
        bam=write_bam(root/'different_from_sample.bam',refs,records)
        manifest=root/'manifest.tsv';manifest.write_text('sample\tbam\treference\nwell_A1\t'+str(bam)+'\t'+str(fa)+'\n')
        intended='[{"pos":301,"ref":"A","alt":"C"}]' if window[0] <= 301 <= window[1] else '[]'
        windows=root/'windows.tsv';windows.write_text('oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\tintended_variants\noligo_1\tchrII\t%d\t%d\t%s\n'%(*window,intended))
        bt=root/'bridge.tsv';bt.write_text('bc1\toligo_name\n'+BC+'\toligo_1\n')
        sidecar={'barcode_locus':'chrII:101-120','barcode_left_sequence':GLEFT,'barcode_right_sequence':GRIGHT,'plasmid_bc_locus':'plasmid_demo:101-120','plasmid_bc_left_sequence':PLEFT,'plasmid_bc_right_sequence':PRIGHT}
        flat=root/'locus.txt';flat.write_text(''.join(k+": '"+v+"'\n" for k,v in sidecar.items()))
        data={'cohort_name':'synthetic','sequencing_date':'20260919','paths':{'project_root':str(root),'manifest':str(manifest)},'locus':sidecar,'bc1':{'use_plasmid':plasmid,'bridge_table':str(bt) if bridge else None},'reconcile':{'oligo_windows':str(windows)}}
        config=root/'config.yaml';config.write_text(yaml.safe_dump(data))
        return dict(root=root,fa=fa,bam=bam,flat=flat,config=load_config(config),config_path=config,refs=refs,bridge=bt,windows=windows,records=records)
    return make

def read_rows(path):
    with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))


def read(name, seq, *, start=300, cigar=None, flag=0, qual=35):
    return dict(name=name, contig='chrII', start=start, seq=seq,
                cigar=cigar or f'{len(seq)}M', flag=flag, qual=qual, mapq=60)


def site(*alleles):
    return r.EditSite('chrII', 301, 320, 'oligo_1', tuple(r.IntendedAllele(*a) for a in alleles))


def score(f, es):
    with pysam.AlignmentFile(str(f['bam'])) as bam, pysam.FastaFile(str(f['fa'])) as fa:
        return r.score_edit_site(bam, fa, es)


def windows(f, alleles, *, name='oligo_1'):
    f['windows'].write_text('oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\tintended_variants\n'
        + name+'\tchrII\t301\t320\t'+json.dumps([dict(zip(('pos','ref','alt'), a)) for a in alleles])+'\n')


def test_assignment_and_counts_survive_actual_stage(fixture_factory):
    f=fixture_factory(bridge=False,plasmid=False)
    windows(f,[(301,'A','C')])
    write_bam(f['bam'],f['refs'],f['records']+[read(f'edit{i}','C'+'A'*19) for i in range(4)])
    assignments=f['root']/'assigned.tsv';assignments.write_text('sample\toligo_name\nwell_A1\toligo_1\n')
    config=replace(f['config'],reconcile=replace(f['config'].reconcile,oligo_assignments=assignments))
    s1.run_step_1(config)
    row=read_rows(s4.run_step_4(config))[0]
    assert (row['sample'],row['label'],row['n_full'],row['confidence'])==('well_A1','SINGLE_EDIT','1','high')
    assert json.loads(row['allele_counts'])==[{'oligo':'oligo_1','chrom':'chrII','pos':301,'ref':'A','alt':'C',
        'alt_count':4,'ref_count':0,'other_count':0,'depth':4,'low_mapq_reads':0,'unknown_mapq_reads':0}]


@pytest.mark.parametrize('seqs,label,status,depth,alts', [(['C'+'A'*19]*4,'SINGLE_EDIT','full',4,4),
    (['A'*20]*4,'NO_EDIT','wt',4,0),(['C'+'A'*19]*2+['A'*20]*2,'CONTAMINATION','partial',4,2),
    ([], 'LOW_CONF','nocov',0,0),(['AC'+'A'*18]*4,'NO_EDIT','wt',4,0),
    (['G'+'A'*19]*4,'NO_EDIT','other',4,0)])
def test_actual_stage_positive_negative_mixed_no_coverage_noise(fixture_factory,seqs,label,status,depth,alts):
    f=fixture_factory()
    windows(f,[(301,'A','C')])
    write_bam(f['bam'],f['refs'],f['records']+[read(f'edit{i}',seq) for i,seq in enumerate(seqs)])
    s1.run_step_1(f['config'])
    rows=read_rows(s4.run_step_4(f['config']))
    assert len(rows)==1
    row=rows[0]
    assert row['label']==label and row['edits'].endswith(':'+status)
    counts=json.loads(row['allele_counts'])
    assert (counts[0]['depth'],counts[0]['alt_count'])==(depth,alts)


def test_assignment_overrides_conflicting_bridge_without_hiding_conflict(fixture_factory):
    f=fixture_factory()
    windows(f,[(310,'A','G')],name='explicit')
    assignments=f['root']/'a.tsv';assignments.write_text('sample\toligo_name\nwell_A1\texplicit\n')
    config=replace(f['config'],reconcile=replace(f['config'].reconcile,oligo_assignments=assignments))
    write_bam(f['bam'],f['refs'],f['records']+[read(f'edit{i}','C'+'A'*19) for i in range(4)])
    s1.run_step_1(config)
    row=read_rows(s4.run_step_4(config))[0]
    assert row['label']=='NO_EDIT' and row['confidence']=='low'
    assert 'BRIDGE_ASSIGNMENT_CONFLICT' in row['flags']
    assert 'explicit@' in row['edits']


@pytest.mark.parametrize('rows',['other\toligo_1\n','well_A1\tmissing\n','well_A1\t\n',''])
def test_bad_assignment_fails_actual_stage(fixture_factory,rows):
    f=fixture_factory();windows(f,[(301,'A','C')])
    p=f['root']/'a.tsv';p.write_text('sample\toligo_name\n'+rows)
    config=replace(f['config'],reconcile=replace(f['config'].reconcile,oligo_assignments=p))
    s1.run_step_1(config)
    with pytest.raises(ValueError):s4.run_step_4(config)


def test_bridge_and_windows_load_once_for_cohort(fixture_factory,monkeypatch):
    f=fixture_factory();windows(f,[(301,'A','C')])
    s1.run_step_1(f['config'])
    calls={'windows':0,'bridge':0}
    oldw,oldb=s4.load_oligo_windows,s1.load_bridge
    def w(*a,**kw):calls['windows']+=1;return oldw(*a,**kw)
    def b(*a,**kw):calls['bridge']+=1;return oldb(*a,**kw)
    monkeypatch.setattr(s4,'load_oligo_windows',w);monkeypatch.setattr(s1,'load_bridge',b)
    s4.run_step_4(f['config'])
    assert calls=={'windows':1,'bridge':1}


def test_bridge_free_confidence_uses_caller_evidence():
    assert driver.bc_confidence([SimpleNamespace(caller_confidence='high',confidence='low')])=='high'
    assert driver.bc_confidence([])=='low'


def test_stage4_refreshes_stale_intended_assignment(fixture_factory):
    f=fixture_factory();windows(f,[(301,'A','C')])
    write_bam(f['bam'],f['refs'],f['records']+[read(f'edit{i}','C'+'A'*19) for i in range(4)])
    s1.run_step_1(f['config'])
    assert read_rows(s4.run_step_4(f['config']))[0]['label']=='SINGLE_EDIT'
    windows(f,[(310,'A','G')])
    row=read_rows(s4.run_step_4(f['config']))[0]
    assert row['label']=='NO_EDIT'
    assert json.loads(row['allele_counts'])[0]['pos']==310


def test_stage4_refreshes_stale_prerequisite_bam(fixture_factory):
    f=fixture_factory();windows(f,[(301,'A','C')])
    write_bam(f['bam'],f['refs'],f['records']+[read(f'edit{i}','C'+'A'*19) for i in range(4)])
    s1.run_step_1(f['config'])
    assert read_rows(s4.run_step_4(f['config']))[0]['label']=='SINGLE_EDIT'
    write_bam(f['bam'],f['refs'],[read(f'wt{i}','A'*20) for i in range(4)])
    row=read_rows(s4.run_step_4(f['config']))[0]
    assert (row['label'],row['n_genome_bc'],row['n_full'])==('LOW_CONF','0','0')


def test_repeat_equivalent_designs_in_different_windows_are_one_event(fixture_factory):
    f=fixture_factory(records=[read(f'r{i}','A'*19,cigar='8M1D11M') for i in range(4)])
    first=site((305,'AA','A'))
    second=r.EditSite('chrII',302,321,'second',(r.IntendedAllele(308,'AA','A'),))
    with pysam.AlignmentFile(str(f['bam'])) as bam,pysam.FastaFile(str(f['fa'])) as fa:
        calls=driver.score_sites(bam,fa,['a','b'],{'a':first,'b':second})
    assert len(calls)==1 and calls[0].status=='full'


def test_stage4_explicit_mapping_quality_sensitivity_override(fixture_factory):
    f=fixture_factory();windows(f,[(301,'A','C')])
    low=[read(f'r{i}','C'+'A'*19) for i in range(4)]
    for record in low:record['mapq']=0
    write_bam(f['bam'],f['refs'],f['records']+low)
    s1.run_step_1(f['config'])
    baseline=read_rows(s4.run_step_4(f['config']))[0]
    assert baseline['label']=='LOW_CONF'
    assert json.loads(baseline['allele_counts'])[0]['low_mapq_reads']==4
    cfg=replace(f['config'],reconcile=replace(f['config'].reconcile,min_mapq=0))
    sensitivity=read_rows(s4.run_step_4(cfg))[0]
    assert sensitivity['label']=='SINGLE_EDIT'
    assert json.loads(sensitivity['allele_counts'])[0]['alt_count']==4




def mixed_case(fixture_factory):
    f=fixture_factory()
    records=f['records']+[dict(name=f'template{i}',contig='chrII',start=300,seq='C'+'A'*9,
        cigar='10M',mate_start=300,tlen=10) for i in range(4)]+[
        dict(name=f'genome{i}',contig='chrII',start=295,seq='A'*20,cigar='20M',mate_start=350,tlen=75) for i in range(4)]
    write_bam(f['bam'],f['refs'],records)
    assignments=f['root']/'assigned.tsv';assignments.write_text('sample\toligo_name\nwell_A1\toligo_1\n')
    blacklist=f['root']/'blacklist.tsv';blacklist.write_text('plasmid_demo\t1\t500\n')
    cfg=replace(f['config'],reconcile=replace(f['config'].reconcile,oligo_assignments=assignments),
        sv_plasmid=replace(f['config'].sv_plasmid,blacklist=blacklist,blacklist_reference=f['fa'],samtools=shutil.which('samtools')))
    return f,cfg


def test_auto_uses_cleaned_alleles_and_retains_raw_barcode_evidence(fixture_factory):
    f,cfg=mixed_case(fixture_factory)
    row=read_rows(s4.run_step_4(cfg))[0]
    evidence=json.loads(row['allele_counts'])[0]
    assert (row['label'],evidence['alt_count'],evidence['ref_count'],evidence['depth'])==('NO_EDIT',0,4,4)
    assert (row['n_genome_bc'],row['n_plasmid_bc'])==('1','1')
    assert 'DONOR_CLEANED_READS' in row['flags']


def test_explicit_raw_policy_retains_mixture_and_labels_source(fixture_factory):
    f,cfg=mixed_case(fixture_factory)
    cfg=replace(cfg,reconcile=replace(cfg.reconcile,read_source='raw'))
    row=read_rows(s4.run_step_4(cfg))[0]
    evidence=json.loads(row['allele_counts'])[0]
    assert (row['label'],evidence['alt_count'],evidence['ref_count'],evidence['depth'])==('CONTAMINATION',4,4,8)
    assert 'RAW_DONOR_READS' in row['flags']
    assert row['confidence']=='low' and 'RAW_DONOR_ORIGIN_UNRESOLVED' in row['flags']


def test_explicit_cleaned_policy_cannot_silently_fall_back(fixture_factory):
    f=fixture_factory()
    cfg=replace(f['config'],reconcile=replace(f['config'].reconcile,read_source='donor_cleaned'))
    with pytest.raises(RuntimeError,match='cleaned|Half B'):
        s4.run_step_4(cfg)


def test_explicit_two_donor_assignments_reach_multi_edit(fixture_factory):
    f=fixture_factory(plasmid=False)
    windows(f,[(301,'A','C')])
    with f['windows'].open('a') as out:
        out.write('oligo_2\tchrII\t401\t410\t[{"pos":401,"ref":"A","alt":"G"}]\n')
    assignments=f['root']/'assigned.tsv'
    assignments.write_text('sample\toligo_name\nwell_A1\toligo_1\nwell_A1\toligo_2\n')
    cfg=replace(f['config'],reconcile=replace(f['config'].reconcile,oligo_assignments=assignments))
    write_bam(f['bam'],f['refs'],f['records']+[read(f'left{i}','C',start=300) for i in range(4)]+
              [read(f'right{i}','G',start=400) for i in range(4)])
    row=read_rows(s4.run_step_4(cfg))[0]
    assert (row['label'],row['n_full'],row['n_partial'])==('MULTI_EDIT','2','0')
    assert [(e['alt_count'],e['depth']) for e in json.loads(row['allele_counts'])]==[(4,4),(4,4)]


def test_unresolved_barcodes_are_not_an_implicit_design_expectation(fixture_factory):
    f=fixture_factory(bridge=False)
    write_bam(f['bam'],f['refs'],f['records']+[read(f'r{i}','C'+'A'*19) for i in range(4)])
    row=read_rows(s4.run_step_4(f['config']))[0]
    assert (row['label'],row['n_full'],row['allele_counts'])==('LOW_CONF','0','[]')
    assert 'UNRESOLVED_BC1' in row['flags']


def test_direct_driver_rejects_explicit_empty_assignment(fixture_factory):
    f=fixture_factory()
    gt,pt=build_tracks(f['flat'])
    with pysam.AlignmentFile(str(f['bam'])) as bam,pysam.FastaFile(str(f['fa'])) as fa:
        with pytest.raises(ValueError,match='explicit assignment is empty'):
            driver.reconcile_sample(bam,fa,'well_A1',genome_track=gt,plasmid_track=pt,
                assigned_oligos=[],oligo_windows={})


@pytest.mark.parametrize('text',['sample\tsample\toligo_name\nwell_A1\twell_A1\toligo_1\n',
    'sample\twrong\nwell_A1\toligo_1\n','sample\toligo_name\n   \toligo_1\n',
    'sample\toligo_name\nwell_A1\t   \n'])
def test_assignment_header_and_whitespace_errors_fail(fixture_factory,text):
    f=fixture_factory();p=f['root']/'assigned.tsv';p.write_text(text)
    cfg=replace(f['config'],reconcile=replace(f['config'].reconcile,oligo_assignments=p))
    with pytest.raises(ValueError):s4.run_step_4(cfg)
    assert not s4.shard_path(cfg,'well_A1').exists()


def test_single_sample_cache_reuses_verified_result_and_duplicate_assignment_once(fixture_factory,monkeypatch):
    f=fixture_factory(plasmid=False)
    p=f['root']/'assigned.tsv';p.write_text('sample\toligo_name\nwell_A1\toligo_1\nwell_A1\toligo_1\n')
    cfg=replace(f['config'],reconcile=replace(f['config'].reconcile,oligo_assignments=p))
    write_bam(f['bam'],f['refs'],f['records']+[read(f'r{i}','C'+'A'*19) for i in range(4)])
    output=s4.run_step_4(cfg,'well_A1')
    assert output==s4.shard_path(cfg,'well_A1')
    before=output.read_bytes()
    def unexpected(*args,**kwargs):raise AssertionError('verified cached output was recomputed')
    monkeypatch.setattr(s4,'reconcile_sample',unexpected)
    assert s4.run_step_4(cfg,'well_A1').read_bytes()==before
    assert read_rows(output)[0]['n_full']=='1'


def test_missing_prerequisite_after_failed_refresh_does_not_emit_result(fixture_factory,monkeypatch):
    f=fixture_factory()
    monkeypatch.setattr(s1,'run_step_1',lambda *args,**kwargs:None)
    with pytest.raises(StepIncompleteError,match='verified Stage-1'):
        s4.run_step_4(f['config'],'well_A1')
    assert not s4.shard_path(f['config'],'well_A1').exists()


def test_cleaned_zero_coverage_is_unassessable_and_keeps_raw_barcodes(fixture_factory):
    f,cfg=mixed_case(fixture_factory)
    with pysam.AlignmentFile(str(f['bam'])) as bam:
        records=[r for r in bam.fetch(until_eof=True) if not r.query_name.startswith('genome')]
        header=bam.header
    with pysam.AlignmentFile(str(f['bam']),'wb',header=header) as out:
        for record in records:out.write(record)
    pysam.index(str(f['bam']))
    row=read_rows(s4.run_step_4(cfg))[0]
    e=json.loads(row['allele_counts'])[0]
    assert (row['label'],row['n_nocov'],e['alt_count'],e['ref_count'],e['depth'])==('LOW_CONF','1',0,0,0)
    assert (row['n_genome_bc'],row['n_plasmid_bc'])==('1','1')


def test_switching_raw_to_cleaned_invalidates_result_and_preserves_barcode_source(fixture_factory):
    f,cfg=mixed_case(fixture_factory)
    raw=replace(cfg,reconcile=replace(cfg.reconcile,read_source='raw'))
    first=read_rows(s4.run_step_4(raw))[0]
    second=read_rows(s4.run_step_4(cfg))[0]
    assert [json.loads(x['allele_counts'])[0]['alt_count'] for x in (first,second)]==[4,0]
    assert [json.loads(x['allele_counts'])[0]['depth'] for x in (first,second)]==[8,4]
    assert (first['n_genome_bc'],first['n_plasmid_bc'])==(second['n_genome_bc'],second['n_plasmid_bc'])==('1','1')
def test_missing_windows_cannot_emit_unmeasured_no_edit_result(fixture_factory):
    f=fixture_factory()
    cfg=replace(f['config'],reconcile=replace(f['config'].reconcile,oligo_windows=None))
    with pytest.raises(StepIncompleteError,match='unmeasured sites are not no-edit'):
        s4.run_step_4(cfg)
    assert not s4.shard_path(cfg,'well_A1').exists()
    assert not s4.table_path(cfg).exists()
