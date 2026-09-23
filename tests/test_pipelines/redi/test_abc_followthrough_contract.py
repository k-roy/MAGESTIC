"""Independent input-boundary contracts across REDI A–C.

Synthetic molecules and registrations are declared before invoking the stages;
none of the failure expectations comes from a production callset.
"""
import gzip
import json
from pathlib import Path

import pandas as pd
import pytest

from magestic.pipelines.redi.core.config import REDIConfig, REDIPaths
from magestic.pipelines.redi.pipeline import step_a_trim_counts as stage_a
from magestic.pipelines.redi.pipeline.step_c_process_bcs import run_step_c
from magestic.pipelines.redi.readers.keyfiles import _read_key, combine_keyfiles

BC = 'ACGTATGCCATTCAGGTAAC'
RB = 'ACCCATAGTCCATTCGGCTA'


def fixture(tmp_path):
    paths = REDIPaths(tmp_path, tmp_path/'raw', tmp_path/'intermediate',
                      tmp_path/'processed', tmp_path/'keys')
    paths.keyfiles.mkdir()
    paths.counts_tables.mkdir(parents=True)
    config = REDIConfig('declared_fixture', '20250627', paths, redi_bc_keys=['keys.tsv'])
    name = '20250627_plate_4A_sample_1'
    metadata = pd.DataFrame([dict(sample_number='sample_1', outer_primers='plate_4A',
        inner_primers='KR1967-KR2225', purpose='yT182 REDI', REDI_library='yT182',
        REDI_plate_barcode='3567', MAT_alpha_array_1536_barcode='6156',
        registration_status='registered', MAT_a_barcode_prefix='AAGTATGT',
        MAT_a_barcode_length='20')])
    combined = paths.keyfiles/'combined.tsv'
    metadata.to_csv(combined, sep='\t', index=False)
    pd.DataFrame([dict(REDI_library='yT182', REDI_bc1_seq=RB,
        plate_384='RED_384', row_384='A', col_384='1')]).to_csv(
            paths.keyfiles/'keys.tsv', sep='\t', index=False)
    reads = pd.DataFrame([dict(sample_name=name,
        R1_seq='CAGGTCATGCTC'+BC+'CCCTAGGATACG',
        R2_seq='AAGTATGT'+RB+'GATCATAACTTC', counts=7)])
    count_file = paths.counts_tables/f'{name}_R1_R2_counts.tsv'
    reads.to_csv(count_file, sep='\t', index=False)
    return config, combined, metadata, reads, count_file


@pytest.mark.parametrize('payload', [
    'a\tb\n1\t2\t3\n',
    'a\tb\n1\n',
    'a\tb\n1\t2\n\n3\t4\n',
])
def test_key_reader_refuses_ragged_or_blank_physical_rows(tmp_path, payload):
    source = tmp_path/'malformed.tsv'
    source.write_text(payload)
    with pytest.raises(ValueError, match='row|column|field|blank'):
        _read_key(source)


def test_key_reader_preserves_explicit_empty_cells_and_literal_identifiers(tmp_path):
    source = tmp_path/'valid.tsv'
    source.write_text('sample\toptional\tlabel\n0007\t\tNA\n')
    frame = _read_key(source)
    assert frame.to_dict(orient='records') == [dict(sample='0007', optional='', label='NA')]


@pytest.mark.parametrize('status,plate', [('registered', 'empty'), ('declared_empty', '3567')])
def test_stage_c_refuses_contradictory_control_registration(tmp_path, status, plate):
    config, combined, metadata, reads, count_file = fixture(tmp_path)
    metadata['registration_status'] = status
    metadata['REDI_plate_barcode'] = plate
    metadata.to_csv(combined, sep='\t', index=False)
    with pytest.raises(ValueError, match='empty|control|registration'):
        run_step_c(config, combined)
    assert not config.paths.combined_data.exists()


@pytest.mark.parametrize('changes', [
    {'REDI_library': 'yT182_typo'},
    {'MAT_a_barcode_prefix': '', 'MAT_a_barcode_length': '0'},
    {'MAT_a_barcode_length': '21'},
])
def test_enabled_assay_requires_registered_library_and_callable_key_length(tmp_path, changes):
    config, combined, metadata, reads, count_file = fixture(tmp_path)
    for column, value in changes.items():
        metadata[column] = value
    metadata.to_csv(combined, sep='\t', index=False)
    with pytest.raises(ValueError, match='library|length|anchor|schema|key'):
        run_step_c(config, combined)
    assert not config.paths.combined_data.exists()


@pytest.mark.parametrize('column,value', [('counts', '9000'), ('bc1', 'T'*20),
                                         ('REDI_row_384', 'P')])
def test_stage_c_refuses_reserved_metadata_even_without_assigned_reads(tmp_path, column, value):
    config, combined, metadata, reads, count_file = fixture(tmp_path)
    metadata[column] = value
    metadata.to_csv(combined, sep='\t', index=False)
    reads['R1_seq'] = 'T'*80
    reads.to_csv(count_file, sep='\t', index=False)
    with pytest.raises(ValueError, match='reserved|overlap|collision'):
        run_step_c(config, combined)
    assert not config.paths.combined_data.exists()


def test_stage_b_cannot_overwrite_a_source_control_classification(tmp_path):
    sample = tmp_path/'sample.tsv'
    sample.write_text('gDNA_plate_name\tsample_number\tREDI_plate_barcode\tregistration_status\n'
                      'plate\tsample_1\t3567\tdeclared_empty\n')
    gdna = tmp_path/'gdna.tsv'
    gdna.write_text('gDNA_plate_name\touter_primers\tinner_primers\nplate\tplate_4A\tKR1967-KR2225\n')
    plate = tmp_path/'plate.tsv'
    plate.write_text('REDI_plate_barcode\tMAT_alpha_array_1536_barcode\n3567\t6156\n')
    with pytest.raises(ValueError, match='registration|reserved|control'):
        combine_keyfiles(sample, gdna, plate)


@pytest.mark.parametrize('headers', [
    ('@same/2', '@same/1'),
    ('@same/1', '@same/1'),
    ('@same 2:N:0:ACGT', '@same 1:N:0:ACGT'),
    ('@same/1 2:N:0:ACGT', '@same/2 2:N:0:ACGT'),
])
def test_mate_identity_includes_declared_read_number(tmp_path, headers):
    paths = [tmp_path/'R1.fastq.gz', tmp_path/'R2.fastq.gz']
    for path, header in zip(paths, headers):
        with gzip.open(path, 'wt') as handle:
            handle.write(header+'\nACGT\n+\nIIII\n')
    with pytest.raises(ValueError, match='mate|read number|identity'):
        stage_a._validate_pairs(*paths)


@pytest.mark.parametrize('headers', [('@same/1', '@same/2'),
    ('@same 1:N:0:ACGT', '@same 2:N:0:ACGT'), ('@same', '@same')])
def test_unambiguous_or_unspecified_mate_labels_remain_supported(tmp_path, headers):
    paths = [tmp_path/'R1.fastq.gz', tmp_path/'R2.fastq.gz']
    for path, header in zip(paths, headers):
        with gzip.open(path, 'wt') as handle:
            handle.write(header+'\nACGT\n+\nIIII\n')
    assert stage_a._validate_pairs(*paths) == 1
