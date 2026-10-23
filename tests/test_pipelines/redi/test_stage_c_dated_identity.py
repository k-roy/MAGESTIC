"""An explicit dated assay identity must agree with the configured run."""
import pandas as pd
import pytest

from magestic.pipelines.redi.core.config import REDIConfig, REDIPaths
from magestic.pipelines.redi.pipeline.step_c_process_bcs import BC_ONLY, run_step_c


def inputs(tmp_path, supplied_date):
    paths = REDIPaths(tmp_path, tmp_path/'raw', tmp_path/'intermediate',
                      tmp_path/'processed', tmp_path/'keys')
    paths.keyfiles.mkdir()
    paths.counts_tables.mkdir(parents=True)
    config = REDIConfig('date_fixture', '20250627', paths, redi_bc_keys=[])
    name = '20250627_plate_4A_sample_1'
    combined = paths.keyfiles/'combined.tsv'
    pd.DataFrame([dict(sample_number='sample_1', outer_primers='plate_4A',
        inner_primers='KR1967-KR2225', purpose=BC_ONLY,
        sequencing_date=supplied_date)]).to_csv(combined, sep='\t', index=False)
    pd.DataFrame([dict(sample_name=name, counts=2,
        R1_seq='CAGGTCATGCTC'+'ACGTATGCCATTCAGGTAAC'+'CCCTAGGATACG',
        R2_seq='ACGT')]).to_csv(paths.counts_tables/f'{name}_R1_R2_counts.tsv',
                              sep='\t', index=False)
    return config, combined


def test_explicit_registration_date_cannot_be_silently_replaced(tmp_path):
    config, combined = inputs(tmp_path, '20220101')
    with pytest.raises(ValueError, match='date|dated'):
        run_step_c(config, combined)
    assert not config.paths.combined_data.exists()


def test_matching_explicit_registration_date_preserves_valid_assay(tmp_path):
    config, combined = inputs(tmp_path, '20250627')
    _, output = run_step_c(config, combined)
    rows = pd.read_csv(output, sep='\t', dtype=str)
    assert rows.to_dict(orient='records') == [dict(
        sample_name='20250627_plate_4A_sample_1', bc1='ACGTATGCCATTCAGGTAAC', counts='2')]
