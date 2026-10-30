"""Portable Stage 3 functional regression tests with independent synthetic BAMs.

Requires pytest, pysam and samtools on PATH. Fixtures create all data under pytest's
temporary directory; no private data, network, audit checkout or fixture download
is required. Expected read geometry, alleles and counts are fixed before calling
the implementation. Concurrent writer cases use two bounded threads.
"""
from dataclasses import replace
from pathlib import Path
import csv
import shutil
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor

import pysam
import pytest

from trace_crispr.wgs import sv_plasmid as sv

C = 'integrated_plasmid_oligo_1'
SAMTOOLS = shutil.which('samtools')


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

def read_rows(path):
    with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))


import yaml
from magestic.pipelines.wgs.pipeline import step_3_sv_plasmid as s3

BC = 'ACGTCAGTTCGATGACCTGA'
GLEFT, GRIGHT = 'ACCCTGTTAT', 'GAGCATGACC'
PLEFT, PRIGHT = 'GGTCATGCTC', 'ATAACAGGGT'
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

def alignment(name='r', start=95, cigar='20M', seq='A' * 20, **kwargs):
    return dict(name=name, contig=C, start=start, cigar=cigar, seq=seq, **kwargs)

def half_b_paths(f):
    blacklist = f['root'] / 'blacklist.tsv'; blacklist.write_text('plasmid_demo\t1\t500\n')
    return dict(in_bam=str(f['bam']), out_bam=str(f['root'] / 'out.bam'), chrom='chrII', donor_start=301,
                donor_end=310, blacklist_path=str(blacklist), samtools=SAMTOOLS,
                reference_path=str(f['fa']), blacklist_reference=str(f['fa']))

def test_full_stage_no_targets_reports_status_and_blank_counts(fixture_factory):
    f = fixture_factory(records=[])
    cfg = replace(f['config'], sv_plasmid=replace(f['config'].sv_plasmid, samtools=SAMTOOLS))
    paths = s3.run_step_3(cfg)
    row = read_rows(paths[0])[0]; summary = read_rows(paths[1])[0]
    assert row['half_a_status'] == summary['half_a_status'] == 'not_run_no_integrated_contigs'
    assert row['total_count'] == summary['total_count'] == ''

def test_full_half_b_real_samtools_drops_template_preserves_genome(fixture_factory):
    assert SAMTOOLS, 'real samtools is required for this acceptance test'
    reads = [dict(name='template', contig='chrII', start=300, seq='A' * 10, cigar='10M', mate_start=300, tlen=10), dict(name='genome', contig='chrII', start=295, seq='A' * 20, cigar='20M', mate_start=350, tlen=75)]
    f = fixture_factory(records=reads)
    assignment = f['root'] / 'assign.tsv'; assignment.write_text('sample\toligo_name\nwell_A1\toligo_1\n')
    blacklist = f['root'] / 'blacklist.tsv'; blacklist.write_text('plasmid_demo\t1\t500\n')
    cfg = replace(f['config'], sv_plasmid=replace(f['config'].sv_plasmid, samtools=SAMTOOLS, blacklist=blacklist, blacklist_reference=f['fa']), reconcile=replace(f['config'].reconcile, oligo_assignments=assignment))
    paths = s3.run_step_3(cfg); row = read_rows(paths[0])[0]
    assert [row[k] for k in ['read_count_mid_point', 'frags_genome', 'frags_plasmid', 'frags_chimeric']] == ['1', '1', '1', '0']
    with pysam.AlignmentFile(str(s3.cleaned_bam_path(cfg, 'well_A1'))) as bam:
        assert [r.query_name for r in bam] == ['genome']

@pytest.mark.parametrize('reads,expected', [([alignment(name='up'), alignment(name='down', start=295)], 2), ([alignment(start=100)], 0), ([alignment(start=85, cigar='5M30D15M')], 0)])
def test_full_stage_derives_and_validates_junctions(fixture_factory, reads, expected):
    contig = 'A' * 90 + 'C' * 10 + 'G' * 200 + 'C' * 10 + 'A' * 90
    f = fixture_factory(records=reads, extra_refs={C: contig})
    donor = f['root'] / 'donors.tsv'; donor.write_text('oligo_name\tdonor_sequence\noligo_1\tCCCCCCCCCC\n')
    cfg = replace(f['config'], sv_plasmid=replace(f['config'].sv_plasmid, donor_sequences=donor, bp_to_extend=90))
    paths = s3.run_step_3(cfg)
    assert read_rows(paths[0])[0]['total_count'] == str(expected)
    assert read_rows(paths[1])[0]['total_count'] == str(expected)
    assert read_rows(paths[0])[0]['half_a_status'] == 'run'

def test_derived_junction_checks_actual_donor_sequence(fixture_factory):
    f = fixture_factory(records=[], extra_refs={C: 'A' * 400})
    donor = f['root'] / 'donors.tsv'; donor.write_text('oligo_name\tdonor_sequence\noligo_1\tCCCCCCCCCC\n')
    cfg = replace(f['config'], sv_plasmid=replace(f['config'].sv_plasmid, donor_sequences=donor, bp_to_extend=90))
    with pytest.raises(s3.StepIncompleteError, match='reference sequence'):
        s3.run_step_3(cfg)
    assert not s3.shard_path(cfg, 'well_A1').exists()

def test_stage3_rejects_conflicting_sample_assignments(fixture_factory):
    f = fixture_factory(); path = f['root'] / 'assign.tsv'
    path.write_text('sample\toligo_name\nwell_A1\tone\nwell_A1\ttwo\n')
    cfg = replace(f['config'], reconcile=replace(f['config'].reconcile, oligo_assignments=path))
    with pytest.raises(s3.StepIncompleteError, match='conflicting'):
        s3._assigned_oligo(cfg, 'well_A1')

def test_stage3_assignment_cache_invalidates_content(fixture_factory):
    f = fixture_factory(); path = f['root'] / 'assign.tsv'
    path.write_text('sample\toligo_name\nwell_A1\tone\n')
    assert s3._assignment_index(path, 'oligo_name') == {'well_A1': 'one'}
    path.write_text('sample\toligo_name\nwell_A1\ttwo\n')
    assert s3._assignment_index(path, 'oligo_name') == {'well_A1': 'two'}

def test_stage3_cache_checks_cleaned_bam_membership(fixture_factory):
    f = fixture_factory(records=[]); args = half_b_paths(f)
    assignment = f['root'] / 'assign.tsv'; assignment.write_text('sample\toligo_name\nwell_A1\toligo_1\n')
    cfg = replace(f['config'], sv_plasmid=replace(f['config'].sv_plasmid, samtools=SAMTOOLS, blacklist=Path(args['blacklist_path']), blacklist_reference=f['fa']), reconcile=replace(f['config'].reconcile, oligo_assignments=assignment))
    s3.run_step_3(cfg)
    cleaned = s3.cleaned_bam_path(cfg, 'well_A1'); cleaned.unlink()
    s3.run_step_3(cfg)
    assert cleaned.is_file()
    with pysam.AlignmentFile(str(cleaned)) as bam:
        assert list(bam) == []

def test_stage3_recovers_truncated_shard_instead_of_trusting_exists(fixture_factory):
    f = fixture_factory(records=[])
    s3.run_step_3(f['config'])
    shard = s3.shard_path(f['config'], 'well_A1'); shard.write_text('')
    s3.run_step_3(f['config'])
    assert read_rows(shard)[0]['half_a_status'] == 'not_run_no_integrated_contigs'

def test_stage3_removes_obsolete_generated_cleaned_bam(fixture_factory):
    f = fixture_factory(records=[]); args = half_b_paths(f)
    assignment = f['root'] / 'assign.tsv'; assignment.write_text('sample\toligo_name\nwell_A1\toligo_1\n')
    cfg = replace(f['config'], sv_plasmid=replace(f['config'].sv_plasmid, samtools=SAMTOOLS, blacklist=Path(args['blacklist_path']), blacklist_reference=f['fa']), reconcile=replace(f['config'].reconcile, oligo_assignments=assignment))
    s3.run_step_3(cfg)
    cleaned = s3.cleaned_bam_path(cfg, 'well_A1'); assert cleaned.exists()
    cfg = replace(cfg, reconcile=replace(cfg.reconcile, oligo_assignments=None))
    s3.run_step_3(cfg)
    assert not cleaned.exists() and not Path(str(cleaned) + '.target_filtering.anno').exists()
    assert read_rows(s3.shard_path(cfg, 'well_A1'))[0]['half_b_status'] == 'not_run_no_oligo_assignment'

def test_stage3_rebuilds_missing_cleaned_index(fixture_factory):
    f = fixture_factory(records=[]); args = half_b_paths(f)
    assignment = f['root'] / 'assign.tsv'; assignment.write_text('sample\toligo_name\nwell_A1\toligo_1\n')
    cfg = replace(f['config'], sv_plasmid=replace(f['config'].sv_plasmid, samtools=SAMTOOLS, blacklist=Path(args['blacklist_path']), blacklist_reference=f['fa']), reconcile=replace(f['config'].reconcile, oligo_assignments=assignment))
    s3.run_step_3(cfg)
    index = s3.cleaned_bam_index_path(cfg, 'well_A1')
    assert index.is_file(); index.unlink()
    s3.run_step_3(cfg)
    assert index.is_file()
    with pysam.AlignmentFile(str(s3.cleaned_bam_path(cfg, 'well_A1'))) as bam:
        assert bam.check_index() and list(bam.fetch('chrII', 300, 310)) == []

def test_obsolete_cleaned_index_removed_with_bam(fixture_factory):
    f = fixture_factory(records=[]); args = half_b_paths(f)
    assignment = f['root'] / 'assign.tsv'; assignment.write_text('sample\toligo_name\nwell_A1\toligo_1\n')
    cfg = replace(f['config'], sv_plasmid=replace(f['config'].sv_plasmid, samtools=SAMTOOLS, blacklist=Path(args['blacklist_path']), blacklist_reference=f['fa']), reconcile=replace(f['config'].reconcile, oligo_assignments=assignment))
    s3.run_step_3(cfg)
    index = s3.cleaned_bam_index_path(cfg, 'well_A1'); assert index.exists()
    cfg = replace(cfg, reconcile=replace(cfg.reconcile, oligo_assignments=None))
    s3.run_step_3(cfg)
    assert not index.exists()

def test_same_sample_concurrent_stage_calls_compute_once(fixture_factory, monkeypatch):
    f = fixture_factory(records=[])
    first_inside = threading.Event(); second_inside = threading.Event(); release = threading.Event()
    state = {'calls': 0}
    def half_b(config, record):
        state['calls'] += 1
        (first_inside if state['calls'] == 1 else second_inside).set()
        assert release.wait(3), 'test did not release the held computation'
        return s3.HALF_B_NO_ASSIGNMENT, ['', '', '', ''], ['', '', '', '', '', '']
    monkeypatch.setattr(s3, '_run_half_b', half_b)
    monkeypatch.setattr(s3, 'cache_key', lambda *a, **kw: {'inputs': [], 'test_identity': 'same_sample'})
    with ThreadPoolExecutor(max_workers=2) as pool:
        one = pool.submit(s3.run_step_3, f['config'], 'well_A1')
        assert first_inside.wait(3)
        two = pool.submit(s3.run_step_3, f['config'], 'well_A1')
        overlapped = second_inside.wait(0.25)
        release.set(); one.result(timeout=3); two.result(timeout=3)
    assert not overlapped, 'two same-sample writers entered computation simultaneously'
    assert state['calls'] == 1

def test_failed_index_publication_never_sets_completion_and_retry_recovers(fixture_factory, monkeypatch):
    f = fixture_factory(records=[]); args = half_b_paths(f)
    assignment = f['root'] / 'assign.tsv'; assignment.write_text('sample\toligo_name\nwell_A1\toligo_1\n')
    cfg = replace(f['config'], sv_plasmid=replace(f['config'].sv_plasmid, samtools=SAMTOOLS, blacklist=Path(args['blacklist_path']), blacklist_reference=f['fa']), reconcile=replace(f['config'].reconcile, oligo_assignments=assignment))
    publish = s3.publish
    def fail_index(source, dest):
        if Path(dest).suffix in ('.bai', '.csi'):
            raise OSError('injected index publication failure')
        return publish(source, dest)
    monkeypatch.setattr(s3, 'publish', fail_index)
    with pytest.raises(OSError, match='publication failure'):
        s3.run_step_3(cfg, 'well_A1')
    shard = s3.shard_path(cfg, 'well_A1')
    assert not shard.exists() and not Path(str(shard) + '.cache.json').exists()
    monkeypatch.setattr(s3, 'publish', publish)
    s3.run_step_3(cfg, 'well_A1')
    with pysam.AlignmentFile(str(s3.cleaned_bam_path(cfg, 'well_A1'))) as bam:
        assert bam.check_index() and list(bam.fetch('chrII', 300, 310)) == []
    assert Path(str(shard) + '.cache.json').exists()

@pytest.mark.parametrize('text', [
    'oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\n',
    'oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\noligo_1\tchrII\tbad\t310\n',
    'oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\noligo_1\tchrII\t0\t310\n',
    'oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\noligo_1\tchrII\t301\t310\noligo_1\tchrII\t302\t310\n',
])
def test_invalid_stage3_windows_fail_before_output(fixture_factory, text):
    f = fixture_factory(records=[]); f['windows'].write_text(text)
    with pytest.raises(s3.StepIncompleteError):
        s3._donor_windows(f['config'])

def test_missing_assignment_file_is_not_a_biological_negative(fixture_factory):
    f = fixture_factory(records=[])
    cfg = replace(f['config'], reconcile=replace(f['config'].reconcile, oligo_assignments=f['root'] / 'missing.tsv'))
    with pytest.raises(s3.StepIncompleteError, match='assignments missing'):
        s3.run_step_3(cfg, 'well_A1')
    assert not s3.shard_path(cfg, 'well_A1').exists()

def test_stage3_missing_blacklist_provenance_fails_before_output(fixture_factory):
    f = fixture_factory(records=[]); args = half_b_paths(f)
    assignment = f['root'] / 'assign.tsv'; assignment.write_text('sample\toligo_name\nwell_A1\toligo_1\n')
    cfg = replace(f['config'], sv_plasmid=replace(f['config'].sv_plasmid, blacklist=Path(args['blacklist_path'])), reconcile=replace(f['config'].reconcile, oligo_assignments=assignment))
    with pytest.raises(s3.StepIncompleteError, match='blacklist_reference'):
        s3.run_step_3(cfg, 'well_A1')

def test_unknown_sample_rejected_before_lock_file_creation(fixture_factory):
    f = fixture_factory(records=[])
    with pytest.raises(KeyError, match='manifest'):
        s3.run_step_3(f['config'], '../../escaped_sample')
    assert not list(f['root'].rglob('*.lock'))

def test_stage_and_standalone_cleaner_share_final_output_lock(fixture_factory, monkeypatch):
    f = fixture_factory(records=[]); args = half_b_paths(f)
    assignments = f['root'] / 'assign.tsv'
    assignments.write_text('sample\toligo_name\nwell_A1\toligo_1\n')
    cfg = replace(f['config'], sv_plasmid=replace(f['config'].sv_plasmid,
        samtools=SAMTOOLS, blacklist=Path(args['blacklist_path']), blacklist_reference=f['fa']),
        reconcile=replace(f['config'].reconcile, oligo_assignments=assignments))
    args['out_bam'] = str(s3.cleaned_bam_path(cfg, 'well_A1'))
    original = sv.classify_donor_reads
    first_inside = threading.Event(); second_inside = threading.Event(); release = threading.Event()
    state = {'calls': 0}
    def classify(*a, **kw):
        state['calls'] += 1
        (first_inside if state['calls'] == 1 else second_inside).set()
        assert release.wait(5)
        return original(*a, **kw)
    monkeypatch.setattr(sv, 'classify_donor_reads', classify)
    with ThreadPoolExecutor(max_workers=2) as pool:
        one = pool.submit(s3.run_step_3, cfg, 'well_A1')
        assert first_inside.wait(5)
        two = pool.submit(sv.clean_donor_reads, **args)
        overlapped = second_inside.wait(0.25)
        release.set(); one.result(timeout=5); two.result(timeout=5)
    assert not overlapped, 'Stage 3 and standalone cleaner wrote the same destination concurrently'
    assert state['calls'] == 2
    with pysam.AlignmentFile(args['out_bam']) as bam:
        assert bam.check_index()


# planning/328 (Kevin, 2026-10-30): "count contained read pairs when the backbone is empty ... a per
# clone/per dataset option". The contained fragment is 301..310, inside the donor window 301..310.
CONTAINED = dict(name='template', contig='chrII', start=300, seq='A' * 10, cigar='10M', mate_start=300, tlen=10)
GENOMIC = dict(name='genome', contig='chrII', start=295, seq='A' * 20, cigar='20M', mate_start=350, tlen=75)


def _deep(n=3000):
    # nuclear depth away from the donor and the barcode locus, so the plasmid bound is informative
    return [dict(name=f'n{i}', contig='chrII', start=i % 60, seq='A' * 20, cigar='20M', mate_start=400, tlen=420)
            for i in range(n)]


def _policy_run(fixture_factory, records, **sv):
    f = fixture_factory(records=records)
    assignment = f['root'] / 'assign.tsv'; assignment.write_text('sample\toligo_name\nwell_A1\toligo_1\n')
    blacklist = f['root'] / 'blacklist.tsv'; blacklist.write_text('plasmid_demo\t1\t500\n')
    cfg = replace(f['config'], sv_plasmid=replace(
        f['config'].sv_plasmid, samtools=SAMTOOLS, blacklist=blacklist, blacklist_reference=f['fa'],
        donor_plasmid_contigs=['plasmid_demo'], plasmid_mask=['plasmid_demo:91-130'], nuclear_contigs=['chrII'], **sv),
        reconcile=replace(f['config'].reconcile, oligo_assignments=assignment))
    row = read_rows(s3.run_step_3(cfg)[0])[0]
    with pysam.AlignmentFile(str(s3.cleaned_bam_path(cfg, 'well_A1'))) as bam:
        kept = sorted({r.query_name for r in bam if not r.query_name.startswith('n')})
    return row, kept


def test_contained_default_is_the_span_rule(fixture_factory):
    row, kept = _policy_run(fixture_factory, [CONTAINED, GENOMIC] + _deep())
    assert row['contained_policy'] == 'discard' and row['frags_contained_kept'] == '0'
    assert row['frags_plasmid'] == '1' and kept == ['genome']


def test_contained_auto_keeps_when_the_backbone_is_provably_empty(fixture_factory):
    row, kept = _policy_run(fixture_factory, [CONTAINED, GENOMIC] + _deep(), contained_fragments='auto')
    assert row['contained_policy'] == 'keep_purged'
    assert row['plasmid_backbone_reads'] == '0' and row['donor_backbone_pairs'] == '0'
    assert float(row['plasmid_copies_upper95']) < 0.05
    assert row['frags_contained_kept'] == '1' and row['frags_plasmid'] == '0' and row['frags_genome'] == '1'
    assert kept == ['genome', 'template']


def test_contained_auto_discards_when_too_shallow_to_bound(fixture_factory):
    row, kept = _policy_run(fixture_factory, [CONTAINED, GENOMIC], contained_fragments='auto')
    assert row['contained_policy'] == 'discard_not_purged' and kept == ['genome']


def test_contained_keep_still_discards_a_sample_that_shows_plasmid(fixture_factory):
    pair = dict(name='donor_mate', contig='chrII', start=301, seq='A' * 8, cigar='8M', flag=97,
                mate='plasmid_demo', mate_start=300, tlen=0)
    row, kept = _policy_run(fixture_factory, [CONTAINED, GENOMIC, pair] + _deep(), contained_fragments='keep')
    assert row['contained_policy'] == 'discard_plasmid_seen' and row['donor_backbone_pairs'] == '1'
    assert 'template' not in kept


def test_contained_keep_is_honoured_without_plasmid(fixture_factory):
    row, kept = _policy_run(fixture_factory, [CONTAINED, GENOMIC], contained_fragments='keep')
    assert row['contained_policy'] == 'keep_declared' and kept == ['genome', 'template']


@pytest.mark.parametrize('sv,message', [
    ({'contained_fragments': 'sometimes'}, 'discard, auto or keep'),
    ({'contained_fragments': 'auto', 'plasmid_mask': []}, 'donor_plasmid_contigs'),
    ({'contained_fragments': 'auto', 'donor_plasmid_contigs': ['p']}, 'explicit plasmid_mask'),
    ({'contained_fragments': 'keep', 'donor_plasmid_contigs': ['p'], 'plasmid_mask': [], 'dup_mode': 'parity'}, 'standard-mode'),
    ({'contained_fragments': 'auto', 'donor_plasmid_contigs': ['p'], 'plasmid_mask': ['p:10']}, 'contig:start-end'),
    ({'purged_max_copies': 0}, 'positive number'),
])
def test_contained_policy_config_is_refused_when_incomplete(fixture_factory, sv, message):
    from magestic.pipelines.wgs.config import load_config
    f = fixture_factory(records=[])
    data = yaml.safe_load(f['config_path'].read_text()); data['sv_plasmid'] = sv
    path = f['root'] / 'bad.yaml'; path.write_text(yaml.safe_dump(data))
    with pytest.raises(ValueError, match=message):
        load_config(path)


# planning/328 v7 (Kevin, 2026-10-30): Stage 4 haplotype-level calling is opt-in and refuses an incomplete setup.
@pytest.mark.parametrize('recon,message', [
    ({'calling': 'sometimes'}, 'floors or v7'),
    ({'calling': 'v7'}, 'cds_table'),
    ({'calling': 'v7', 'cds_table': 'cds.tsv', 'allele_evidence': 'aligned'}, 'realignment'),
    ({'contamination_rate': 1.5}, r'\(0, 1\)'),
])
def test_v7_calling_config_is_refused_when_incomplete(fixture_factory, recon, message):
    from magestic.pipelines.wgs.config import load_config
    f = fixture_factory(records=[])
    data = yaml.safe_load(f['config_path'].read_text()); data['reconcile'].update(recon)
    path = f['root'] / 'bad_v7.yaml'; path.write_text(yaml.safe_dump(data))
    with pytest.raises(ValueError, match=message):
        load_config(path)


def test_v7_driver_refuses_a_window_without_a_target(fixture_factory):
    from magestic.pipelines.wgs.core import reconcile_driver as rd
    from trace_crispr.wgs.reconcile import EditSite, IntendedAllele
    f = fixture_factory(records=[])
    site_windows = {'oligo_1': None}
    class _TR:   # only build_edit_sites is reached before the refusal
        SITE_PARTIAL = 'partial'
        @staticmethod
        def build_edit_sites(oligos, windows, fasta=None):
            return [EditSite('chrII', 301, 310, 'oligo_1', (IntendedAllele(301, 'A', 'C'),))]
    import magestic.pipelines.wgs.core.reconcile_driver as mod
    real = mod._trace.reconcile
    mod._trace.reconcile = lambda: _TR
    try:
        with pysam.AlignmentFile(str(f['bam'])) as bam, pysam.FastaFile(str(f['fa'])) as fa:
            with pytest.raises(ValueError, match='target codon'):
                rd.score_sites(bam, fa, ['oligo_1'], site_windows, calling='v7', targets={})
    finally:
        mod._trace.reconcile = real
