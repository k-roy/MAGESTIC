"""Declared YAML settings and paths cannot silently change the assay."""
from pathlib import Path
import pytest,yaml
from magestic.pipelines.redi.core.config import load_config

def config(tmp_path,change=None):
 value=dict(screen_name='dated_assay',sequencing_date='20250707',paths=dict(project_root='project',raw_fastq='raw',intermediate='intermediate',processed='processed',keyfiles='keys'),bc1_reference_table='reference.tsv',wgs_outcome_table='wgs.tsv',wgs_sample_namespace='namespace.tsv')
 if change:change(value)
 p=tmp_path/'config.yaml';p.write_text(yaml.safe_dump(value));return p

def test_relative_paths_are_bound_to_config_then_project(tmp_path,monkeypatch):
 p=config(tmp_path);other=tmp_path/'elsewhere';other.mkdir();monkeypatch.chdir(other)
 actual=load_config(p);project=tmp_path/'project'
 assert actual.paths.project_root==project
 assert actual.paths.raw_fastq==project/'raw' and actual.paths.keyfiles==project/'keys'
 assert actual.bc1_reference_table==project/'reference.tsv'
 assert actual.wgs_sample_namespace==project/'namespace.tsv'

@pytest.mark.parametrize('name',['allow_partial_sample','reconcile','unexpected'])
def test_unknown_top_level_setting_is_not_ignored(tmp_path,name):
 p=config(tmp_path,lambda c:c.update({name:True}))
 with pytest.raises(ValueError,match='unknown|unsupported'):load_config(p)

@pytest.mark.parametrize('value',['false',0,1,[]])
def test_boolean_flags_are_actual_booleans(tmp_path,value):
 p=config(tmp_path,lambda c:c.update(allow_partial_samples=value))
 with pytest.raises(ValueError,match='boolean'):load_config(p)

@pytest.mark.parametrize('group,value',[('purity',{'min_counts':-1}),('purity',{'top_n_for_purity':0}),('barcode',{'bc1_length':0}),('barcode',{'bc1_length_tolerance':True})])
def test_invalid_thresholds_rejected_at_load(tmp_path,group,value):
 p=config(tmp_path,lambda c:c.update({group:value}))
 with pytest.raises(ValueError,match='integer|positive|nonnegative|length'):load_config(p)


def test_duplicate_yaml_keys_refused(tmp_path):
 p=config(tmp_path);p.write_text(p.read_text()+'sequencing_date: "20220101"\n')
 with pytest.raises(ValueError,match='duplicate'):load_config(p)


def test_same_array_names_cannot_generate_overlapping_columns(tmp_path):
 p=config(tmp_path,lambda c:c.update(array_pairs=[dict(name_a='same',name_b='same',purpose_a='a',purpose_b='b')]))
 with pytest.raises(ValueError,match='distinct|same|different'):load_config(p)


def test_zero_count_floor_is_refused_before_cross_array_execution(tmp_path):
 p=config(tmp_path,lambda c:c.update(purity={'min_counts':0}))
 with pytest.raises(ValueError,match='integer|positive'):load_config(p)


def test_distinct_pair_tuples_cannot_overwrite_same_output(tmp_path):
 p=config(tmp_path,lambda c:c.update(array_pairs=[dict(name_a='a_b',name_b='c',purpose_a='first',purpose_b='second'),dict(name_a='a',name_b='b_c',purpose_a='third',purpose_b='fourth')]))
 with pytest.raises(ValueError,match='duplicate|collision'):load_config(p)


def test_duplicate_second_purpose_is_refused_before_reading_data(tmp_path):
 p=config(tmp_path,lambda c:c.update(array_pairs=[dict(name_a='old',name_b='new',purpose_a='a',purpose_b=['b','b'])]))
 with pytest.raises(ValueError,match='duplicate|repeat'):load_config(p)


def test_standard_yaml_merge_can_override_one_array_name(tmp_path):
 p=config(tmp_path)
 p.write_text(p.read_text()+'''array_pairs:
  - &oldpair {name_a: old, name_b: first, purpose_a: earlier, purpose_b: later}
  - {<<: *oldpair, name_b: second}
''')
 actual=load_config(p)
 assert [(p.name_a,p.name_b) for p in actual.array_pairs]==[('old','first'),('old','second')]


def test_lowercase_anchors_preserve_supported_extraction_semantics(tmp_path):
 p=config(tmp_path,lambda c:c.update(barcode={'fwd_bc1_prefix':'acgt','fwd_bc1_suffix':'ttag'}))
 actual=load_config(p)
 assert actual.barcode.fwd_bc1_prefix.upper()=='ACGT'
 assert actual.barcode.fwd_bc1_suffix.upper()=='TTAG'


def test_explicit_prefix_only_assay_and_inactive_null_offset_remain_loadable(tmp_path):
 p=config(tmp_path,lambda c:c.update(barcode={'fwd_bc1_suffix':None,'fallback_prefix_length':None}))
 actual=load_config(p)
 assert actual.barcode.fwd_bc1_suffix is None and not actual.barcode.allow_positional_fallback


def test_active_positional_mode_cannot_use_missing_offset(tmp_path):
 p=config(tmp_path,lambda c:c.update(barcode={'allow_positional_fallback':True,'fallback_prefix_length':None}))
 with pytest.raises(ValueError,match='integer'):load_config(p)


@pytest.mark.parametrize('setting',[{'min_purity':float('nan')},{'min_purity':1.01},{'min_replicates':True},{'min_counts':-1}])
def test_plot_filters_refuse_invalid_values(tmp_path,setting):
 p=config(tmp_path,lambda c:c.update(plotting=setting))
 with pytest.raises(ValueError):load_config(p)


def test_standard_plot_does_not_invent_independent_replicates(tmp_path):
 actual=load_config(config(tmp_path))
 assert actual.plotting.min_replicates==0


@pytest.mark.parametrize('mapping',[
 'purity: {<<: {min_counts: 4, min_counts: 6}}\n',
 'purity: {<<: [{min_counts: 4, min_counts: 6}]}\n',
])
def test_duplicate_keys_cannot_hide_inside_inline_yaml_merge(tmp_path,mapping):
 p=config(tmp_path);p.write_text(p.read_text()+mapping)
 with pytest.raises(ValueError,match='duplicate'):load_config(p)
