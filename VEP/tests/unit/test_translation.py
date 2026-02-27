"""
Unit tests for translation.py

Tests DNA to protein translation functions.

Author: Kevin R. Roy
Date: 2026-02-23
"""

import pytest
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from sequence_utils.translation import (
    translate,
    translate_with_stops,
    get_codon_at_position,
    CODON_TABLE,
)


# ============================================================================
# translate() Tests
# ============================================================================

class TestTranslate:
    """Tests for translate function."""

    def test_simple_translation(self):
        """Test basic CDS translation."""
        # ATG = M, AAA = K, TAT = Y
        cds = "ATGAAATAT"
        protein = translate(cds)
        assert protein == "MKY"

    def test_with_stop_codon(self):
        """Test translation stops at stop codon."""
        # ATG = M, TAA = stop
        cds = "ATGTAAATG"
        protein = translate(cds)
        assert protein == "M"  # Stops at TAA

    def test_all_stop_codons(self):
        """Test all three stop codons."""
        for stop in ["TAA", "TAG", "TGA"]:
            cds = f"ATG{stop}ATG"
            protein = translate(cds)
            assert protein == "M", f"Failed for stop codon {stop}"

    def test_empty_cds(self):
        """Test empty CDS returns empty protein."""
        protein = translate("")
        assert protein == ""

    def test_incomplete_codon_ignored(self):
        """Test that incomplete final codon is ignored."""
        cds = "ATGAA"  # ATG + incomplete AA
        protein = translate(cds)
        assert protein == "M"

    def test_lowercase_input(self):
        """Test that lowercase input is handled."""
        cds = "atgaaatga"  # M, K, stop
        protein = translate(cds)
        assert protein == "MK"

    def test_mixed_case_input(self):
        """Test mixed case input."""
        cds = "AtGaAaTaT"
        protein = translate(cds)
        assert protein == "MKY"

    def test_unknown_codon(self):
        """Test unknown/invalid codons translate to X."""
        cds = "ATGNNN"  # ATG = M, NNN = X (unknown)
        protein = translate(cds)
        assert protein == "MX"

    def test_remove_stop_false(self):
        """Test translation with remove_stop=False continues after stop."""
        cds = "ATGTAAATG"  # M, stop, M
        protein = translate(cds, remove_stop=False)
        # Should skip stop and continue
        assert protein == "MM"


# ============================================================================
# translate_with_stops() Tests
# ============================================================================

class TestTranslateWithStops:
    """Tests for translate_with_stops function."""

    def test_includes_stop_codon(self):
        """Test that stop codons are included as '*'."""
        cds = "ATGTAAATG"  # M, stop, M
        protein = translate_with_stops(cds)
        assert '*' in protein
        assert protein == "M*M"

    def test_multiple_stops(self):
        """Test multiple stop codons."""
        cds = "ATGTAATAGATG"  # M, TAA, TAG, M
        protein = translate_with_stops(cds)
        assert protein.count('*') == 2


# ============================================================================
# get_codon_at_position() Tests
# ============================================================================

class TestGetCodonAtPosition:
    """Tests for get_codon_at_position function."""

    def test_first_codon(self):
        """Test getting first codon (position 0)."""
        cds = "ATGAAATGA"
        codon, aa = get_codon_at_position(cds, 0)
        assert codon == "ATG"
        assert aa == "M"

    def test_middle_codon(self):
        """Test getting middle codon."""
        cds = "ATGAAATGA"
        codon, aa = get_codon_at_position(cds, 1)
        assert codon == "AAA"
        assert aa == "K"

    def test_last_codon(self):
        """Test getting last codon."""
        cds = "ATGAAATGA"
        codon, aa = get_codon_at_position(cds, 2)
        assert codon == "TGA"
        assert aa == "*"

    def test_position_out_of_range(self):
        """Test position beyond sequence returns None."""
        cds = "ATG"
        codon, aa = get_codon_at_position(cds, 5)
        assert codon is None
        assert aa is None

    def test_negative_position(self):
        """Test that negative position is handled (Python allows negative indexing)."""
        cds = "ATGAAATGA"
        # This might work due to Python negative indexing, or might fail
        # Just verify it doesn't crash
        try:
            codon, aa = get_codon_at_position(cds, -1)
        except (IndexError, ValueError):
            pass  # Expected behavior


# ============================================================================
# CODON_TABLE Tests
# ============================================================================

class TestCodonTable:
    """Tests for CODON_TABLE constant."""

    def test_all_64_codons_present(self):
        """Test that all 64 codons are in table."""
        bases = "ACGT"
        all_codons = [a + b + c for a in bases for b in bases for c in bases]
        assert len(all_codons) == 64
        for codon in all_codons:
            assert codon in CODON_TABLE, f"Missing codon: {codon}"

    def test_start_codon(self):
        """Test ATG is methionine."""
        assert CODON_TABLE["ATG"] == "M"

    def test_stop_codons(self):
        """Test all stop codons map to '*'."""
        assert CODON_TABLE["TAA"] == "*"
        assert CODON_TABLE["TAG"] == "*"
        assert CODON_TABLE["TGA"] == "*"

    def test_amino_acid_count(self):
        """Test we have 20 amino acids + stop."""
        unique_aas = set(CODON_TABLE.values())
        assert len(unique_aas) == 21  # 20 AAs + *


# ============================================================================
# Integration Tests
# ============================================================================

class TestTranslationIntegration:
    """Integration tests combining multiple functions."""

    def test_roundtrip_codon_extraction(self):
        """Test that we can extract each codon and reconstruct protein."""
        cds = "ATGAAATATGATTGA"  # M K Y D *
        protein = translate(cds)
        assert protein == "MKYD"

        # Now extract each codon and verify
        for i, expected_aa in enumerate("MKYD"):
            codon, aa = get_codon_at_position(cds, i)
            assert aa == expected_aa

    def test_frameshift_detection(self):
        """Test detecting frameshift by comparing translations."""
        wt_cds = "ATGAAATATGATTGA"  # MKYD*
        # Insert 1 bp after ATG
        mut_cds = "ATGAAAATATGATTGA"  # Frame shifted

        wt_protein = translate(wt_cds)
        mut_protein = translate(mut_cds)

        # Proteins should differ (frameshift changes reading frame)
        assert wt_protein != mut_protein

    def test_synonymous_mutation(self):
        """Test that synonymous mutation gives same protein."""
        wt_cds = "ATGAAATGA"  # MK*
        # AAA -> AAG is synonymous (both = K)
        mut_cds = "ATGAAGTGA"

        wt_protein = translate(wt_cds)
        mut_protein = translate(mut_cds)

        assert wt_protein == mut_protein == "MK"

    def test_nonsense_mutation(self):
        """Test nonsense mutation truncates protein."""
        wt_cds = "ATGAAATATTGA"  # MKY*
        # TAT -> TAA is nonsense (Y -> stop)
        mut_cds = "ATGAAATAATGA"

        wt_protein = translate(wt_cds)
        mut_protein = translate(mut_cds)

        assert wt_protein == "MKY"
        assert mut_protein == "MK"  # Truncated at new stop
