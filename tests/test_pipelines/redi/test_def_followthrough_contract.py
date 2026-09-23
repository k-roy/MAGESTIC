"""Independent D→E→F physical-observation and exact-read contracts.

Fixtures spell out complete read distributions. No planning artifacts, barcode
lookup of observed winners, or generated production expectation is required.
"""
import json

import pandas as pd
import pytest

from magestic.pipelines.redi.core.annotate import annotate_with_bc1_reference
from magestic.pipelines.redi.core.config import REDIArrayPair, REDIConfig, REDIPaths
from magestic.pipelines.redi.core.cross_array import compare_arrays
from magestic.pipelines.redi.pipeline.step_d_purity import run_step_d
from magestic.pipelines.redi.pipeline.step_e_annotate import run_step_e
from magestic.pipelines.redi.pipeline.step_f_cross_array import run_step_f


A = 'A' * 20
C = 'C' * 20
BIG = 9007199254740993


def observation(sample, purpose, distribution, *, column='1', alias='G' * 20):
    maximum = max(distribution.values())
    winners = [bc for bc, count in distribution.items() if count == maximum and count]
    return dict(sample_name=sample, purpose=purpose, REDI_bc=alias,
        MAT_alpha_array_1536_barcode='00123', REDI_plate_color_384='RED_384',
        REDI_row_384='A', REDI_col_384=column,
        bc1=winners[0] if len(winners) == 1 else pd.NA,
        REDI_bc_bc1_counts=maximum, REDI_bc_bc1_total_counts=sum(distribution.values()),
        all_bc1_counts_json=json.dumps(distribution), identity_ambiguous=len(winners) > 1)


def compare(data, **kwargs):
    return compare_arrays(data, data.purpose.eq('earlier'), data.purpose.eq('later'),
                          'old', 'new', **kwargs)


def configuration(tmp_path, *, require_validation=False):
    paths = REDIPaths(tmp_path, tmp_path / 'raw', tmp_path / 'intermediate',
                      tmp_path / 'processed', tmp_path / 'keys')
    paths.combined_data.mkdir(parents=True)
    paths.keyfiles.mkdir()
    return REDIConfig('physical_contract', '20250627', paths,
        array_pairs=[REDIArrayPair('old', 'new', 'earlier', purpose_b='later',
                                  require_validation_a=require_validation)])


def test_outer_position_union_preserves_exact_large_counts_and_missingness():
    # Different physical wells are retained, never joined by their shared bc1.
    data = pd.DataFrame([observation('dated_old', 'earlier', {A: BIG}, column='1'),
                         observation('dated_new', 'later', {A: BIG}, column='2')])
    result = compare(data)
    assert len(result) == 2
    old = result.loc[result.REDI_col_384.eq('1')].iloc[0]
    new = result.loc[result.REDI_col_384.eq('2')].iloc[0]
    for row, side in ((old, 'old'), (new, 'new')):
        # int() exposes loss that float-vs-int comparison can hide by coercion.
        assert int(row[f'{side}_REDI_bc_bc1_counts']) == BIG
        assert int(row[f'{side}_REDI_bc_bc1_total_counts']) == BIG
        assert pd.isna(row['bc1_agreement'])
    assert pd.isna(old['new_REDI_bc_bc1_total_counts'])
    assert pd.isna(new['old_REDI_bc_bc1_total_counts'])


def test_one_sample_barcode_observation_cannot_reappear_under_another_purpose():
    data = pd.DataFrame([observation('same_dated_sample', 'earlier', {A: 20}),
                         observation('same_dated_sample', 'later', {A: 20})])
    with pytest.raises(ValueError, match='duplicate|observation|independent'):
        compare(data)


@pytest.mark.parametrize('column', [
    'validated_by_yT100_mating', 'confirmed_contaminated_by_yT100_mating'])
def test_design_reference_cannot_replace_physical_validation_flags(column):
    top = pd.DataFrame([{'bc1': A, 'counts': 20, 'total_counts': 20, column: False}])
    ref = pd.DataFrame([{'bc1': A, 'donor': 'design', column: True}])
    with pytest.raises(ValueError, match='reserved|observation|collision'):
        annotate_with_bc1_reference(top, ref)


@pytest.mark.parametrize('metadata', [{}, {'confirmed_contaminated_by_yT100_mating': False}])
def test_required_validation_cannot_be_inferred_from_missing_positive_flag(tmp_path, metadata):
    config = configuration(tmp_path, require_validation=True)
    data = pd.DataFrame([observation('dated_old', 'earlier', {A: 20}),
                         observation('dated_new', 'later', {A: 20})])
    for column, value in metadata.items():
        data[column] = value
    source = config.paths.combined_data / 'annotated.tsv'
    data.to_csv(source, sep='\t', index=False)
    with pytest.raises(ValueError, match='validation|validated'):
        run_step_f(config, source)
    assert not (config.paths.combined_data / 'old_new_merged_REDI_df.tsv.provenance.json').exists()


def test_explicit_no_validation_policy_accepts_missing_validation_metadata(tmp_path):
    config = configuration(tmp_path, require_validation=False)
    source = config.paths.combined_data / 'annotated.tsv'
    pd.DataFrame([observation('dated_old', 'earlier', {A: 20}),
                  observation('dated_new', 'later', {A: 20})]).to_csv(source, sep='\t', index=False)
    result = pd.read_csv(run_step_f(config, source)[0], sep='\t')
    assert len(result) == 1 and bool(result.bc1_agreement.iloc[0])


def test_top_only_legacy_observation_with_declared_tie_stays_uncalled():
    data = pd.DataFrame([observation('dated_old', 'earlier', {A: 10, C: 10}),
                         observation('dated_new', 'later', {A: 20})])
    data = data.drop(columns='all_bc1_counts_json')
    # Legacy exports may retain one tied spelling; the explicit ambiguity wins.
    data.loc[0, 'bc1'] = A
    result = compare(data).iloc[0]
    assert pd.isna(result.old_bc1)
    assert bool(result.old_identity_ambiguous)
    assert pd.isna(result.bc1_agreement)


@pytest.mark.parametrize('flag,ambiguous', [('True', True), ('False', False), ('1', True), ('0', False)])
def test_tsv_legacy_ambiguity_flags_do_not_use_string_truthiness(flag, ambiguous):
    data = pd.DataFrame([observation('dated_old', 'earlier', {A: 20}),
                         observation('dated_new', 'later', {A: 20})]).drop(columns='all_bc1_counts_json')
    data['identity_ambiguous'] = [flag, 'False']
    result = compare(data).iloc[0]
    assert bool(result.old_identity_ambiguous) == ambiguous
    if ambiguous:
        assert pd.isna(result.old_bc1) and pd.isna(result.bc1_agreement)
    else:
        assert result.old_bc1 == A and bool(result.bc1_agreement)


@pytest.mark.parametrize('validated,contaminated,callable_old', [
    ('True', 'False', True), ('False', 'False', False), ('True', 'True', False)])
def test_required_validation_uses_positive_and_contamination_evidence(
        tmp_path, validated, contaminated, callable_old):
    config = configuration(tmp_path, require_validation=True)
    data = pd.DataFrame([observation('dated_old', 'earlier', {A: 20}),
                         observation('dated_new', 'later', {A: 20})])
    data['validated_by_yT100_mating'] = [validated, 'False']
    data['confirmed_contaminated_by_yT100_mating'] = [contaminated, 'False']
    source = config.paths.combined_data / 'annotated.tsv'
    data.to_csv(source, sep='\t', index=False)
    result = pd.read_csv(run_step_f(config, source)[0], sep='\t').iloc[0]
    if callable_old:
        assert result.old_bc1 == A and bool(result.bc1_agreement)
    else:
        assert pd.isna(result.old_bc1) and pd.isna(result.bc1_agreement)
        assert result.new_bc1 == A


@pytest.mark.parametrize('first,second,expected_callable,expected_count,expected_total', [
    ({A: 6}, {A: 3}, False, 9, 9),
    ({A: 6}, {A: 4}, True, 10, 10),
    ({A: 10}, {C: 10}, False, 10, 20),
    ({A: 0}, {C: 0}, False, 0, 0),
])
def test_same_well_barcode_alias_pooling_uses_full_counts_before_threshold(
        first, second, expected_callable, expected_count, expected_total):
    data = pd.DataFrame([
        observation('dated_old', 'earlier', first, alias='G' * 20),
        observation('dated_old', 'earlier', second, alias='T' * 20),
        observation('dated_new', 'later', {A: 20}),
    ])
    result = compare(data, min_counts=10).iloc[0]
    assert bool(result.old_callable) == expected_callable
    assert result.old_REDI_bc_bc1_counts == expected_count
    assert result.old_REDI_bc_bc1_total_counts == expected_total
    assert result.old_observations == 2
    if expected_callable:
        assert result.old_bc1 == A and bool(result.bc1_agreement)
    else:
        assert pd.isna(result.old_bc1) and pd.isna(result.bc1_agreement)


@pytest.mark.parametrize('length', [19, 21])
def test_offlength_winner_is_uncalled_but_all_its_reads_stay_in_denominator(length):
    data = pd.DataFrame([observation('dated_old', 'earlier', {'T' * length: 12, A: 11}),
                         observation('dated_new', 'later', {A: 20})])
    result = compare(data, bc1_length=20).iloc[0]
    assert result.old_REDI_bc_bc1_counts == 12
    assert result.old_REDI_bc_bc1_total_counts == 23
    assert result.old_colony_purity == pytest.approx(12 / 23)
    assert pd.isna(result.old_bc1) and pd.isna(result.bc1_agreement)


def test_d_e_f_roundtrip_pools_aliases_and_preserves_ambiguous_design_candidates(tmp_path):
    config = configuration(tmp_path)
    source = config.paths.combined_data / 'registered_counts.tsv'
    distributions = [
        ('dated_old', 'earlier', 'G' * 20, {A: 6, C: 3}),
        ('dated_old', 'earlier', 'T' * 20, {A: 4, C: 5}),
        ('dated_new', 'later', 'G' * 20, {A: 10, C: 8}),
    ]
    raw = []
    for sample, purpose, alias, distribution in distributions:
        for barcode, count in distribution.items():
            raw.append(dict(sample_name=sample, purpose=purpose, REDI_bc=alias,
                MAT_alpha_array_1536_barcode='00123', REDI_plate_color_384='RED_384',
                REDI_row_384='A', REDI_col_384='1', bc1=barcode, counts=count))
    pd.DataFrame(raw).to_csv(source, sep='\t', index=False)
    config.bc1_reference_table = tmp_path / 'reference.tsv'
    pd.DataFrame([{'bc1': A, 'donor': 'one'}, {'bc1': A, 'donor': 'two'},
                  {'bc1': C, 'donor': 'other'}]).to_csv(config.bc1_reference_table, sep='\t', index=False)
    top_path = run_step_d(config, source)
    top = pd.read_csv(top_path, sep='\t', dtype=str, keep_default_na=False)
    assert len(raw) == 6 and len(top) == 3
    assert sum(int(value) for value in top.total_counts) == 36
    annotated_path = run_step_e(config, top_path)
    annotated = pd.read_csv(annotated_path, sep='\t', dtype=str, keep_default_na=False)
    assert len(annotated) == 3
    assert annotated.bc1_reference_ambiguous.eq('True').sum() == 2
    assert annotated.MAT_alpha_array_1536_barcode.eq('00123').all()
    result = pd.read_csv(run_step_f(config, annotated_path)[0], sep='\t').iloc[0]
    # Oracle: first assay has A=6+4=10 and C=3+5=8, including second-place reads.
    assert result.old_bc1 == A and result.new_bc1 == A
    assert result.old_REDI_bc_bc1_counts == result.new_REDI_bc_bc1_counts == 10
    assert result.old_REDI_bc_bc1_total_counts == result.new_REDI_bc_bc1_total_counts == 18
    assert result.old_colony_purity == pytest.approx(10 / 18)
    assert result.old_observations == 2 and result.new_observations == 1
    assert bool(result.bc1_agreement)
