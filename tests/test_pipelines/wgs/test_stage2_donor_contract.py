"""Portable Stage 2 integration contracts with real synthetic BAM/VCF fixtures.

Run with MAGESTIC's WGS dependencies and the corresponding TRACE checkout installed:
    python -m pytest tests/test_pipelines/wgs/test_stage2_donor_contract.py
All fixtures are built under tmp_path; no planning-directory artifact is read.
"""
from dataclasses import replace
import json
from pathlib import Path

import pysam
import pytest
import yaml

from trace_crispr.wgs import donor_outcome as d
from magestic.pipelines.wgs.config import load_config
from magestic.pipelines.wgs.core import samplesheet as sheets
from magestic.pipelines.wgs.core import vcf_adapter as adapter
from magestic.pipelines.wgs.pipeline import step_2_donor_precall as stage
from magestic.pipelines.wgs.pipeline._common import GateError, StepIncompleteError
from magestic.pipelines.wgs.readers import oligo_tables as ot

HEADER = '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n'
CHROM = 'integrated_plasmid_oligo_1'
WT = 'ACGTACGTAA'
REF = 'N' * 10 + WT + 'N' * 10


@pytest.mark.parametrize('strand,pos,ref,alt',[('+',11,'A','C'),('-',14,'T','G')])
def test_stage_palindromic_donor_strand_is_preserved_and_obeyed(stage_case,strand,pos,ref,alt):
    config,table,_,plate=stage_case
    record=stage.manifest_samples(config)[0]
    record.reference.write_text('>'+CHROM+'\n'+'N'*10+'ACGT'+'N'*16+'\n')
    pysam.faidx(str(record.reference))
    table.write_text('oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\twt_seq\tdesign_donor\tdonor_strand\n'+
                     f'oligo\t{CHROM}\t11\t14\tACGT\tCCGT\t{strand}\n')
    plate.write_text(HEADER+vcf_line(pos,ref,alt))
    row=stage.run_step_2(config).read_text().splitlines()[2].split('\t')
    assert row[:3] == ['s','designed_variant_installed',f'{CHROM}_{pos}_{ref}_to_{alt}']
    sheet=stage.samplesheet_path(config,'s')
    assert sheet.read_text().splitlines()[0].split('\t')[-1] == 'donor_strand'
    assert d.read_samplesheet(sheet)[0].donor_strand == strand
    # Keep observed variants fixed; changing only independent metadata must change
    # the answer and invalidate the cached outcome, not fit the design to the reads.
    opposite = '-' if strand == '+' else '+'
    table.write_text(table.read_text().replace('\t'+strand+'\n','\t'+opposite+'\n'))
    changed=stage.run_step_2(config).read_text().splitlines()[2].split('\t')
    assert changed[1] == 'unknown;spontaneous_mutations'
    assert d.read_samplesheet(sheet)[0].donor_strand == opposite


def test_stage_ambiguous_missing_strand_refused_without_outcome(stage_case):
    config,table,_,plate=stage_case
    record=stage.manifest_samples(config)[0]
    record.reference.write_text('>'+CHROM+'\n'+'N'*10+'ACGT'+'N'*16+'\n')
    pysam.faidx(str(record.reference))
    table.write_text('oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\twt_seq\tdesign_donor\n'+
                     f'oligo\t{CHROM}\t11\t14\tACGT\tCCGT\n')
    plate.write_text(HEADER+vcf_line(11,'A','C'))
    with pytest.raises(ValueError,match='donor_strand'):
        stage.run_step_2(config)
    assert not stage.shard_path(config,'s').exists()


def test_stage_contradictory_strand_refused_without_outcome(stage_case):
    config,table,_,_=stage_case
    lines=table.read_text().splitlines()
    table.write_text(lines[0]+'\tdonor_strand\n'+lines[1]+'\t-\n')
    with pytest.raises(ValueError,match='donor_strand'):
        stage.run_step_2(config)
    assert not stage.shard_path(config,'s').exists()


@pytest.mark.parametrize('gt', ['0/.', './.', '.'])
def test_plate_missing_genotype_is_preserved_as_unknown(stage_case, gt):
    config, _, _, plate = stage_case
    plate.write_text(HEADER+vcf_line(20, fmt='GT', value=gt))
    (config.paths.donor_precall_workdir/'gvcf/s.g.vcf').write_text(
        HEADER+vcf_line(11, alt='<NON_REF>', fmt='GT:MIN_DP', value='0:8', info='END=20'))
    row = stage.run_step_2(config).read_text().splitlines()[2].split('\t')
    assert row[1] == 'unknown;ambiguous_genotype'
    assert json.loads(row[5]) == [{'chrom': CHROM, 'pos': 20, 'GT': gt}]


def test_existing_sample_vcf_without_configured_plate(stage_case):
    config, _, _, plate = stage_case
    work = config.paths.donor_precall_workdir
    (work/'vcf').mkdir()
    (work/'vcf/s.fil.vcf').write_text(plate.read_text())
    config = replace(config, donor_precall=replace(config.donor_precall, plate_vcfs=[]))
    row = stage.run_step_2(config).read_text().splitlines()[2].split('\t')
    assert row[:3] == ['s', 'designed_variant_installed', f'{CHROM}_20_A_to_C']


@pytest.mark.parametrize('header', [HEADER.replace('\tFORMAT\t','\tOTHER\t'),
                                  HEADER.replace('\tPOS\t','\tPLACE\t')])
def test_plate_adapter_rejects_malformed_fixed_header_before_publication(tmp_path, header):
    plate = tmp_path/'bad.vcf'; plate.write_text(header+vcf_line(20))
    with pytest.raises(adapter.VcfAdapterError, match='header'):
        adapter.split_plate_vcf(plate, tmp_path/'work')
    assert not list((tmp_path/'work').rglob('*.fil.vcf'))


def vcf_line(pos, ref='A', alt='C', fmt='GT:AD', value='1:0,8', info='.', chrom=CHROM):
    return f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t{info}\t{fmt}\t{value}\n'


def call(tmp_path, vcf='', gvcf='', design='ACGTACGTAC', **kwargs):
    vpath, gpath = tmp_path/'s.vcf', tmp_path/'s.g.vcf'
    vpath.write_text(HEADER + vcf)
    gpath.write_text(HEADER + gvcf)
    spec = d.DonorSpec('s', CHROM, 11, 20, WT, design, design)
    return d.precall_sample(spec, vpath, gpath, REF, **kwargs)


def test_missing_inputs_fail_before_writing_any_output(tmp_path):
    sheet=tmp_path/'sheet.tsv'; ref=tmp_path/'ref.fa';out=tmp_path/'out.tsv'
    sheets.build_samplesheet([dict(sample='s',chrom=CHROM,donor_start_coord=11,donor_end_coord=20,wt_seq=WT,design_donor=WT[:-1]+'C')],sheet,synthetic_donor_source='design_donor')
    ref.write_text('>'+CHROM+'\n'+REF+'\n')
    with pytest.raises(FileNotFoundError,match='gvcf'):
        d.precall(sheet,tmp_path,ref,out)
    assert not out.exists()


@pytest.fixture
def stage_case(tmp_path):
    fa=tmp_path/'ref.fa'; fa.write_text('>'+CHROM+'\n'+REF+'\n');pysam.faidx(str(fa))
    bam=tmp_path/'s.bam'
    header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':CHROM,'LN':len(REF)}], 'PG':[{'ID':'fixture','PN':'bwa','CL':f'bwa mem {fa} fixture.fastq'}]}
    with pysam.AlignmentFile(str(bam),'wb',header=header): pass
    pysam.index(str(bam))
    manifest=tmp_path/'manifest.tsv';manifest.write_text(f'sample\tbam\treference\ns\t{bam}\t{fa}\n')
    table=tmp_path/'donors.tsv';table.write_text('oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\twt_seq\tdesign_donor\n'+f'oligo\t{CHROM}\t11\t20\t{WT}\t{WT[:-1]}C\n')
    assignment=tmp_path/'assign.tsv';assignment.write_text('sample\toligo_name\ns\toligo\n')
    plate=tmp_path/'plate.vcf';plate.write_text(HEADER+vcf_line(20,fmt='AD:GT:DP',value='0,9:1:9'))
    configpath=tmp_path/'config.yaml'
    configpath.write_text(yaml.safe_dump({'cohort_name':'donor_contract','sequencing_date':'20260920','paths':{'project_root':str(tmp_path),'manifest':str(manifest)},'donor_precall':{'donor_table':str(table),'plate_vcfs':[str(plate)]},'reconcile':{'oligo_assignments':str(assignment)}}))
    config=load_config(configpath)
    work=config.paths.donor_precall_workdir; (work/'gvcf').mkdir(parents=True)
    (work/'gvcf/s.g.vcf').write_text(HEADER)
    return config,table,assignment,plate


def test_real_stage2_entrypoint_with_plate_split_donor_join_axis_and_fanin(stage_case):
    config,*_=stage_case
    out=stage.run_step_2(config)
    lines=out.read_text().splitlines()
    assert lines[0] == '# synthetic_error_axis=disabled'
    assert len(lines) == 3
    # planning/268 merge: TRACE's six columns, then the 222 lane's appended outcome-reason
    # columns (step_2.OUTCOME_REASON_COLUMNS), which step_5 passes through by name.
    assert lines[1].split('\t')[6:] == list(stage.OUTCOME_REASON_COLUMNS)
    assert lines[2].split('\t')[:6] == ['s','designed_variant_installed',f'{CHROM}_20_A_to_C','','','[]']
    assert sheets.read_axis_stamp(out) == 'disabled'


def test_stage_missing_gvcf_reaches_real_gate_not_keyword_error(stage_case):
    config,*_=stage_case
    (config.paths.donor_precall_workdir/'gvcf/s.g.vcf').unlink()
    with pytest.raises(GateError,match='no gVCF'):
        stage.run_step_2(config,sample='s')
    assert not stage.shard_path(config,'s').exists()


def test_stage_chooses_sample_plate_and_refuses_multiple_sources(stage_case,tmp_path):
    config,_,_,plate=stage_case
    other=tmp_path/'other.vcf';other.write_text(HEADER.replace('\ts\n','\tother\n'))
    cfg=replace(config,donor_precall=replace(config.donor_precall,plate_vcfs=[other,plate]))
    assert stage.run_step_2(cfg,sample='s').is_file()
    stage.shard_path(cfg,'s').unlink()
    with pytest.raises(StepIncompleteError,match='multiple plate'):
        stage.run_step_2(replace(config,donor_precall=replace(config.donor_precall,plate_vcfs=[plate,plate])),sample='s')


@pytest.mark.parametrize('change', ['missing_assignment','duplicate_assignment','missing_oligo','duplicate_oligo','bad_coordinate','mismatched_reference'])
def test_donor_input_failures_are_explicit(stage_case,change):
    config,table,assignment,_=stage_case
    if change=='missing_assignment': assignment.write_text('sample\toligo_name\nother\toligo\n')
    elif change=='duplicate_assignment': assignment.write_text(assignment.read_text()+'s\toligo\n')
    elif change=='missing_oligo': assignment.write_text('sample\toligo_name\ns\tmissing\n')
    elif change=='duplicate_oligo': table.write_text(table.read_text()+table.read_text().splitlines()[1]+'\n')
    elif change=='bad_coordinate': table.write_text(table.read_text().replace('\t11\t20\t','\t0\t20\t'))
    else:
        lines=table.read_text().splitlines();table.write_text(lines[0]+'\treference\n'+lines[1]+'\t/other/ref.fa\n')
    with pytest.raises((ot.MissingOligoInputError,sheets.SamplesheetError)):
        stage.run_step_2(config,sample='s')
    assert not stage.shard_path(config,'s').exists()


def test_empirical_axis_requires_actual_nonempty_column(stage_case):
    config,*_=stage_case
    config=replace(config,donor_precall=replace(config.donor_precall,synthetic_donor_source='synthetic_donor'))
    with pytest.raises(ot.MissingOligoInputError,match='synthetic_donor'):
        stage.run_step_2(config,sample='s')


def test_explicit_junction_donor_reader(stage_case,tmp_path):
    config,*_=stage_case
    table=tmp_path/'seq.tsv';table.write_text('oligo_name\tdonor_sequence\na\tACGT\nb\tAAA\n')
    config=replace(config,sv_plasmid=replace(config.sv_plasmid,donor_sequences=table))
    assert ot.load_donor_sequences(config)=={'a':'ACGT','b':'AAA'}
    table.write_text(table.read_text()+'a\tTGCA\n')
    with pytest.raises(ot.MissingOligoInputError,match='duplicate'):
        ot.load_donor_sequences(config)


def test_adapter_uses_named_gt_with_ad_first(tmp_path):
    plate=tmp_path/'plate';plate.write_text(HEADER+vcf_line(20,fmt='AD:GT',value='0,8:1')+vcf_line(19,fmt='AD:GT',value='8,0:0'))
    result=adapter.split_plate_vcf(plate,tmp_path/'work',samples=['s'])
    assert result.n_records_in == 2 and result.n_records_out == {'s':1}
    assert '\t20\t' in result.vcf_paths['s'].read_text()
    assert not (tmp_path/'work/gvcf').exists()


def test_cohort_preparation_batches_and_cached_restart_avoids_scans(stage_case, monkeypatch):
    config, _, assignments, plate = stage_case
    record = stage.manifest_samples(config)[0]
    # Two independent outcomes share one source plate, donor table, and reference.
    manifest = config.paths.manifest
    manifest.write_text(manifest.read_text()+f't\t{record.bam}\t{record.reference}\n')
    assignments.write_text(assignments.read_text()+'t\toligo\n')
    plate.write_text(HEADER.replace('\ts\n','\ts\tt\n') +
                     vcf_line(20).rstrip('\n')+'\t0:8,0\n')
    gpath = config.paths.donor_precall_workdir/'gvcf/t.g.vcf'
    gpath.write_text(HEADER.replace('\ts\n','\tt\n')+
                    vcf_line(11,alt='<NON_REF>',fmt='GT:MIN_DP',value='0:8',info='END=20'))
    counts = dict(plate=0, donor=0, reference=0)
    split, read, load = stage.split_plate_vcf, ot.read_oligo_table, d.load_fasta
    def split_count(*args, **kwargs):
        counts['plate'] += 1
        return split(*args, **kwargs)
    def read_count(*args, **kwargs):
        counts['donor'] += 1
        return read(*args, **kwargs)
    def load_count(*args, **kwargs):
        counts['reference'] += 1
        return load(*args, **kwargs)
    monkeypatch.setattr(stage,'split_plate_vcf',split_count)
    monkeypatch.setattr(ot,'read_oligo_table',read_count)
    monkeypatch.setattr(d,'load_fasta',load_count)
    out = stage.run_step_2(config)
    first = out.read_bytes()
    # TRACE's six columns; the appended outcome-reason columns follow (planning/268).
    assert ['\t'.join(ln.split('\t')[:6]) for ln in first.decode().splitlines()[2:]] == [
        f's\tdesigned_variant_installed\t{CHROM}_20_A_to_C\t\t\t[]',
        't\tdesigned_variant_unedited\t\t\t\t[]',
    ]
    # The assignment table is read with the same utility, hence two table reads.
    assert counts == dict(plate=1, donor=2, reference=1)
    stage.run_step_2(config)
    assert out.read_bytes() == first
    assert counts == dict(plate=1, donor=2, reference=1)
    stage.shard_path(config,'t').unlink()
    stage.run_step_2(config)
    assert out.read_bytes() == first
    assert counts == dict(plate=2, donor=4, reference=2)


def test_adapter_duplicate_mapped_identity_refused(tmp_path):
    plate=tmp_path/'plate';plate.write_text(HEADER.replace('\ts\n','\ta\tb\n'))
    with pytest.raises(adapter.VcfAdapterError,match='many-to-one'):
        adapter.split_plate_vcf(plate,tmp_path/'work',sample_map={'a':'s','b':'s'})


def test_adapter_failure_keeps_existing_output(tmp_path):
    work=tmp_path/'work';(work/'vcf').mkdir(parents=True);out=work/'vcf/s.fil.vcf';out.write_text('previous completed output\n')
    plate=tmp_path/'plate';plate.write_text(HEADER+vcf_line(20)+'malformed\n')
    with pytest.raises(adapter.VcfAdapterError):adapter.split_plate_vcf(plate,work)
    assert out.read_text() == 'previous completed output\n'
    assert not list((work/'vcf').glob('*.tmp'))


def test_wrong_gvcf_sample_header_refused(stage_case):
    config,*_=stage_case
    path=config.paths.donor_precall_workdir/'gvcf/s.g.vcf'
    path.write_text(HEADER.replace('\ts\n','\tother_clone\n'))
    with pytest.raises(StepIncompleteError,match='sample header'):
        stage.run_step_2(config,sample='s')
    assert not stage.shard_path(config,'s').exists()


@pytest.mark.parametrize('change', ['duplicate_sample','short_row','zero_start','wrong_wt_length','empty_design','non_dna'])
def test_samplesheet_contract_failures(tmp_path,change):
    sheet=tmp_path/'sheet.tsv'
    row=dict(sample='s',chrom=CHROM,donor_start_coord=11,donor_end_coord=20,wt_seq=WT,design_donor=WT[:-1]+'C')
    rows=[row]
    if change=='duplicate_sample':rows=[row,row]
    elif change=='short_row':del row['chrom']
    elif change=='zero_start':row['donor_start_coord']=0
    elif change=='wrong_wt_length':row['wt_seq']=WT[:-1]
    elif change=='empty_design':row['design_donor']=''
    elif change=='non_dna':row['design_donor']='XXX'
    with pytest.raises(sheets.SamplesheetError):
        sheets.build_samplesheet(rows,sheet,synthetic_donor_source='design_donor')
    assert not sheet.exists()


def test_precall_late_failure_preserves_previous_complete_output(tmp_path):
    sheet=tmp_path/'sheet';ref=tmp_path/'ref';out=tmp_path/'out';out.write_text('previous\n')
    row=dict(sample='s',chrom=CHROM,donor_start_coord=11,donor_end_coord=20,wt_seq=WT,design_donor=WT[:-1]+'C')
    sheets.build_samplesheet([row,dict(row,sample='t')],sheet,synthetic_donor_source='design_donor')
    ref.write_text('>'+CHROM+'\n'+REF+'\n')
    (tmp_path/'vcf').mkdir();(tmp_path/'gvcf').mkdir()
    for sample in ('s','t'):
        (tmp_path/f'vcf/{sample}.fil.vcf').write_text(HEADER+vcf_line(20,'A' if sample=='s' else 'G','C'))
        (tmp_path/f'gvcf/{sample}.g.vcf').write_text(HEADER)
    with pytest.raises(ValueError,match='REF'):d.precall(sheet,tmp_path,ref,out)
    assert out.read_text() == 'previous\n'
    assert not list(tmp_path.glob('out.*.tmp'))


def test_adapter_validates_entire_gt_not_only_first_alt():
    with pytest.raises(adapter.VcfAdapterError):adapter.carries_alt('1/x')


def test_empirical_synthetic_axis_actual_stage(stage_case):
    config,table,_,plate=stage_case
    lines=table.read_text().splitlines()
    empirical=WT[:4]+'G'+WT[5:-1]+'C'
    table.write_text(lines[0]+'\tsynthetic_donor\n'+lines[1]+'\t'+empirical+'\n')
    plate.write_text(HEADER+vcf_line(15,'A','G')+vcf_line(20,'A','C'))
    config=replace(config,donor_precall=replace(config.donor_precall,synthetic_donor_source='synthetic_donor'))
    lines=stage.run_step_2(config).read_text().splitlines()
    assert lines[0]=='# synthetic_error_axis=enabled'
    assert lines[2].split('\t')[:6]==['s','designed_variant_installed;synthetic_error_installed',
                                  f'{CHROM}_20_A_to_C',f'{CHROM}_15_A_to_G','','[]']


def test_successful_stage_provenance_describes_actual_inputs(stage_case):
    config,*_=stage_case
    out=stage.run_step_2(config)
    provenance=json.loads(out.with_suffix(out.suffix+'.provenance.json').read_text())
    params=provenance['parameters']
    assert params['reference_depth_inputs']=='validated_per_sample_gvcf'
    assert params['synthetic_error_axis']=='disabled'
    assert 'gates' not in params
    assert 'diagnostic' in params['unedited_coverage_policy']
