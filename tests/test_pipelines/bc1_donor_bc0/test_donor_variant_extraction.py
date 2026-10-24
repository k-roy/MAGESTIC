"""Step 06 donor variant extraction (register BDB-H16, BDB-H17).

A SNP must be found whether or not the aligner wrote an MD tag, and a soft clip is
unaligned sequence, not an insertion. Internal insertions and deletions are controls
whose VCF-style output is unchanged.
"""
import importlib.util
import sys
from pathlib import Path

import pytest

pysam = pytest.importorskip('pysam')

from magestic.pipelines.bc1_donor_bc0.readers import data_loader

REF = 'GATTACAGGCTTACCGATCGTTAGCCATGA'


def _step06():
    path = Path(data_loader.__file__).parents[1] / 'snakemake/06_annotate_unmatched_donors.py'
    spec = importlib.util.spec_from_file_location('test_variant_extraction_06', path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def reference(tmp_path):
    fasta = tmp_path / 'ref.fa'
    fasta.write_text('>chrT\n' + REF + '\n')
    pysam.faidx(str(fasta))
    return fasta


def _read(cigar, sequence, md=None):
    header = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': 'chrT', 'LN': len(REF)}]})
    read = pysam.AlignedSegment(header)
    read.query_name = 'donor'
    read.reference_id = 0
    read.reference_start = 0
    read.mapping_quality = 60
    read.cigarstring = cigar
    read.query_sequence = sequence
    if md is not None:
        read.set_tag('MD', md)
    return read


def test_snp_is_found_without_an_md_tag(reference):
    sequence = REF[:5] + 'T' + REF[6:12]
    variants = _step06().extract_variants_from_alignment(_read('12M', sequence), reference)
    assert variants == [('chrT', 6, 'C', 'T')]


def test_terminal_soft_clip_is_not_an_insertion(reference):
    variants = _step06().extract_variants_from_alignment(
        _read('10M5S', REF[:10] + 'TTTTT', md='10'), reference)
    assert variants == []


def test_internal_insertion_is_reported_anchored(reference):
    variants = _step06().extract_variants_from_alignment(
        _read('5M2I5M', REF[:5] + 'GG' + REF[5:10], md='10'), reference)
    assert variants == [('chrT', 5, 'A', 'AGG')]


def test_deletion_is_reported_anchored(reference):
    variants = _step06().extract_variants_from_alignment(
        _read('5M2D5M', REF[:5] + REF[7:12], md='5^CA5'), reference)
    assert variants == [('chrT', 5, 'ACA', 'A')]


def test_minimap2_is_asked_for_md_tags(tmp_path, monkeypatch):
    step06 = _step06()
    calls = []
    monkeypatch.setattr(step06.subprocess, 'run', lambda args, **kwargs: calls.append(args))
    assert step06.run_minimap2(tmp_path / 'donors.fa', tmp_path / 'ref.fa', tmp_path / 'out.sam')
    assert '--MD' in calls[0]
