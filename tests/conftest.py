"""
Shared pytest fixtures for MAGESTIC test suite.

Provides synthetic test data (no real Oak/Sherlock paths required) for unit tests.
Integration tests that require real data are marked with @pytest.mark.integration.
"""

import pytest
import pandas as pd
import numpy as np


@pytest.fixture
def tiny_oligo_df():
    """Minimal synthetic oligo design table for testing matching pipelines."""
    guide_seqs = [
        "ATCGATCGATCGATCGATCG",
        "GCTAGCTAGCTAGCTAGCTA",
        "TTTTAAAACCCCGGGGTTTT",
    ]
    donor_seqs = [
        "A" * 129,
        "T" * 129,
        "C" * 80 + "G" * 49,
    ]
    return pd.DataFrame({
        "oligo_name": ["oligo_001", "oligo_002", "oligo_003"],
        "guide_sequence": guide_seqs,
        "donor_sequence": donor_seqs,
        "full_sequence": [
            "N" * 20 + g + "N" * 11 + d + "N" * 20
            for g, d in zip(guide_seqs, donor_seqs)
        ],
    })


@pytest.fixture
def tmp_oligo_tsv(tmp_path, tiny_oligo_df):
    """Write oligo design table to a temp TSV file, return Path."""
    p = tmp_path / "oligos.tsv"
    tiny_oligo_df.to_csv(p, sep="\t", index=False)
    return p


@pytest.fixture
def tiny_gff_content():
    """Minimal synthetic GFF3 content for annotation testing."""
    return """\
##gff-version 3
chrI\tSGD\tgene\t1000\t2000\t.\t+\t.\tID=YAL001C;Name=TFC3
chrI\tSGD\tCDS\t1050\t1950\t.\t+\t.\tParent=YAL001C
chrI\tSGD\tgene\t3000\t4500\t.\t-\t.\tID=YAL002W;Name=VPS8
chrI\tSGD\tCDS\t3100\t4400\t.\t-\t.\tParent=YAL002W
"""


@pytest.fixture
def tmp_gff(tmp_path, tiny_gff_content):
    """Write synthetic GFF to a temp file, return Path."""
    p = tmp_path / "test.gff"
    p.write_text(tiny_gff_content)
    return p


@pytest.fixture
def tiny_fasta_content():
    """Minimal synthetic FASTA for genome testing (2 chroms, 5000 bp each)."""
    rng = np.random.default_rng(42)
    bases = list("ATCG")
    seq = "".join(rng.choice(bases, size=5000))
    return f">chrI\n{seq}\n"


@pytest.fixture
def tmp_fasta(tmp_path, tiny_fasta_content):
    """Write synthetic FASTA to a temp file, return Path."""
    p = tmp_path / "genome.fa"
    p.write_text(tiny_fasta_content)
    return p
