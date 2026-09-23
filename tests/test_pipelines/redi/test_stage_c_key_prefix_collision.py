"""A registered longer molecular key must not become a different shorter-key well."""
import pandas as pd
import pytest

from magestic.pipelines.redi.core.config import REDIConfig, REDIPaths
from magestic.pipelines.redi.pipeline.step_c_process_bcs import run_step_c

SHORT = 'ACCCATAGTCCATTCGGCTA'
LONG = SHORT + 'G'
CLONE = 'ACGTATGCCATTCAGGTAAC'


def inputs(tmp_path, long_library, molecular_key):
    paths = REDIPaths(tmp_path, tmp_path/'raw', tmp_path/'intermediate',
                      tmp_path/'processed', tmp_path/'keys')
    paths.keyfiles.mkdir()
    paths.counts_tables.mkdir(parents=True)
    config = REDIConfig('prefix_fixture', '20250627', paths, redi_bc_keys=['keys.tsv'])
    sample = '20250627_plate_4A_sample_1'
    metadata = paths.keyfiles/'combined.tsv'
    pd.DataFrame([dict(sample_number='sample_1', outer_primers='plate_4A',
        inner_primers='KR1967-KR2225', purpose='yT182 REDI', REDI_library='yT182',
        MAT_a_barcode_prefix='AAGTATGT', MAT_a_barcode_length='20',
        REDI_plate_barcode='3567', MAT_alpha_array_1536_barcode='6156',
        registration_status='registered')]).to_csv(metadata, sep='\t', index=False)
    # Two independently declared physical positions, never inferred from reads.
    pd.DataFrame([
        dict(REDI_library='yT182', REDI_bc1_seq=SHORT, plate_384='RED_384', row_384='A', col_384='1'),
        dict(REDI_library=long_library, REDI_bc1_seq=LONG, plate_384='RED_384', row_384='A', col_384='2'),
    ]).to_csv(paths.keyfiles/'keys.tsv', sep='\t', index=False)
    pd.DataFrame([dict(sample_name=sample, counts=7,
        R1_seq='CAGGTCATGCTC'+CLONE+'CCCTAGGATACG',
        R2_seq='AAGTATGT'+molecular_key+'GATCATAACTTC')]).to_csv(
            paths.counts_tables/f'{sample}_R1_R2_counts.tsv', sep='\t', index=False)
    return config, metadata


def test_longer_registered_key_cannot_be_sliced_into_another_registered_position(tmp_path):
    assert len(SHORT) == 20 and len(LONG) == 21 and LONG[:20] == SHORT
    config, metadata = inputs(tmp_path, 'yT182', LONG)
    # The independent input molecule carries the A2 key. A fixed20 extractor
    # cannot establish A1 identity when both keys belong to this same library.
    with pytest.raises(ValueError, match='ambig|prefix|length'):
        run_step_c(config, metadata)
    assert not config.paths.combined_data.exists()


def test_prefix_collision_in_another_declared_library_does_not_invalidate_assay(tmp_path):
    config, metadata = inputs(tmp_path, 'yT183', SHORT)
    output, _ = run_step_c(config, metadata)
    rows = pd.read_csv(output, sep='\t', dtype=str)
    assert len(rows) == 1
    assert rows.loc[0, 'REDI_library'] == 'yT182'
    assert rows.loc[0, 'REDI_col_384'] == '1'
    assert rows.loc[0, 'counts'] == '7'
