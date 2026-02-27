"""
Pytest fixtures for VEP pipeline unit tests.

Author: Kevin R. Roy
Date: 2026-02-23
"""

import pytest
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Tuple


# ============================================================================
# Test Data Fixtures
# ============================================================================

@pytest.fixture
def sample_dna_sequence():
    """Short DNA sequence for basic tests."""
    return "ATGCATGCATGCATGC"


@pytest.fixture
def sample_protein_sequence():
    """Short protein sequence for basic tests."""
    return "MHACMHACMHAC"


@pytest.fixture
def sample_cds_sequence():
    """CDS sequence that translates to a valid protein."""
    # ATG (M) + AAA (K) + TAT (Y) + TAA (stop) = MKYS*
    return "ATGAAATATTAA"


@pytest.fixture
def sample_cds_with_frameshift():
    """CDS sequence with a frameshift that needs extension."""
    # Original: ATG AAA TAT TAA (12 bp = 4 codons)
    # With insertion: ATG AAA TAT T (11 bp - not divisible by 3)
    return "ATGAAATATT"


@pytest.fixture
def sample_offset_table():
    """Sample offset table for coordinate transformation tests."""
    # Format: list of [s288c_pos, cumulative_offset]
    # This represents:
    #   - Positions 0-99: no offset
    #   - Position 100: 3bp insertion in strain (offset +3)
    #   - Position 200: 5bp deletion in strain (offset -2, net from +3)
    return [
        [0, 0],
        [100, 3],   # 3bp insertion
        [200, -2],  # 5bp deletion (cumulative: 3 - 5 = -2)
        [300, -2],  # No change
    ]


@pytest.fixture
def sample_variants():
    """Sample VCF-style variants for testing."""
    @dataclass
    class MockVariant:
        pos: int
        ref: str
        alt: str

    return [
        MockVariant(pos=50, ref="A", alt="G"),      # SNP
        MockVariant(pos=100, ref="AT", alt="A"),    # 1bp deletion
        MockVariant(pos=150, ref="C", alt="CGG"),   # 2bp insertion
        MockVariant(pos=200, ref="ATGC", alt="A"),  # 3bp deletion
    ]


@pytest.fixture
def sample_pairwise_variant():
    """Sample PairwiseVariant for testing."""
    from sequence_utils.pairwise_variants import PairwiseVariant
    return PairwiseVariant(
        background_strain="S288C",
        source_strain="RM11",
        bg_seq_pos=50,
        bg_allele="A",
        src_allele="G",
        s288c_pos=12350,
    )


# ============================================================================
# Biological Test Cases
# ============================================================================

@pytest.fixture
def ena1_test_case():
    """
    ENA1 complex indel case - 4 linked indels summing to 0bp net change.

    This is a regression test for linked variant detection.
    """
    return {
        "gene": "ENA1",
        "orf": "YDR040C",
        "chrom": "chrIV",
        "variants": [
            {"pos": 566100, "ref": "A", "alt": "AT"},      # +1bp
            {"pos": 566103, "ref": "TT", "alt": "T"},      # -1bp
            {"pos": 566110, "ref": "G", "alt": "GCC"},     # +2bp
            {"pos": 566115, "ref": "AAA", "alt": "A"},     # -2bp
        ],
        "expected_net_change": 0,  # Should sum to zero
        "expected_linked": True,   # Should be detected as linked
    }


@pytest.fixture
def multi_exon_gene_case():
    """
    Multi-exon gene test case for CDS extraction.

    RPS19A has 2 exons and is minus-strand.
    """
    return {
        "gene": "RPS19A",
        "orf": "YOL121C",
        "chrom": "chrXV",
        "strand": "-",
        "exons": [
            (206001, 206300),  # Exon 1
            (206400, 206700),  # Exon 2
        ],
        "expected_cds_length": 600,  # 300 + 300
    }


# ============================================================================
# Directory Fixtures
# ============================================================================

@pytest.fixture
def vep_base_dir():
    """Base directory for VEP pipeline."""
    return Path("/path/to/common/scripts/variant_effect_prediction")


@pytest.fixture
def test_data_dir(vep_base_dir):
    """Test data directory."""
    return vep_base_dir / "tests"


# ============================================================================
# Helper Fixtures
# ============================================================================

@pytest.fixture
def reverse_complement():
    """Reverse complement function for DNA sequences."""
    def _rc(seq: str) -> str:
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(b, 'N') for b in reversed(seq.upper()))
    return _rc


@pytest.fixture
def translate():
    """Translation function for DNA to protein."""
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }

    def _translate(dna: str) -> str:
        protein = []
        for i in range(0, len(dna) - 2, 3):
            codon = dna[i:i+3].upper()
            aa = codon_table.get(codon, 'X')
            if aa == '*':
                break
            protein.append(aa)
        return ''.join(protein)
    return _translate
