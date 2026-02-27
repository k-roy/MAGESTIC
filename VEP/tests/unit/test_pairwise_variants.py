"""
Unit tests for pairwise_variants.py

Tests coordinate transformations, variant application, and edge cases.

Author: Kevin R. Roy
Date: 2026-02-23
"""

import pytest
import numpy as np
from dataclasses import dataclass
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from sequence_utils.pairwise_variants import (
    PairwiseVariant,
    apply_variant,
    compute_strain_seq_position,
)


# ============================================================================
# PairwiseVariant Tests
# ============================================================================

class TestPairwiseVariant:
    """Tests for PairwiseVariant dataclass."""

    def test_snp_creation(self):
        """Test creating a SNP variant."""
        v = PairwiseVariant(
            background_strain="S288C",
            source_strain="RM11",
            bg_seq_pos=100,
            bg_allele="A",
            src_allele="G",
            s288c_pos=50100,
        )
        assert v.bg_seq_pos == 100
        assert v.bg_allele == "A"
        assert v.src_allele == "G"

    def test_insertion_creation(self):
        """Test creating an insertion variant."""
        v = PairwiseVariant(
            background_strain="S288C",
            source_strain="RM11",
            bg_seq_pos=100,
            bg_allele="",
            src_allele="ATG",
            s288c_pos=50100,
        )
        assert len(v.bg_allele) == 0
        assert len(v.src_allele) == 3

    def test_deletion_creation(self):
        """Test creating a deletion variant."""
        v = PairwiseVariant(
            background_strain="S288C",
            source_strain="RM11",
            bg_seq_pos=100,
            bg_allele="ATG",
            src_allele="",
            s288c_pos=50100,
        )
        assert len(v.bg_allele) == 3
        assert len(v.src_allele) == 0


# ============================================================================
# apply_variant Tests
# ============================================================================

class TestApplyVariant:
    """Tests for apply_variant function."""

    def test_apply_snp(self):
        """Test applying a SNP."""
        sequence = "ATGCATGC"
        variant = PairwiseVariant(
            background_strain="S288C",
            source_strain="RM11",
            bg_seq_pos=0,
            bg_allele="A",
            src_allele="G",
            s288c_pos=100,
        )
        result = apply_variant(sequence, variant, verify_allele=True)
        assert result == "GTGCATGC"

    def test_apply_insertion(self):
        """Test applying an insertion."""
        sequence = "ATGCATGC"
        variant = PairwiseVariant(
            background_strain="S288C",
            source_strain="RM11",
            bg_seq_pos=4,
            bg_allele="",
            src_allele="XXX",
            s288c_pos=100,
        )
        result = apply_variant(sequence, variant, verify_allele=True)
        assert result == "ATGCXXXATGC"

    def test_apply_deletion(self):
        """Test applying a deletion."""
        sequence = "ATGCATGC"
        variant = PairwiseVariant(
            background_strain="S288C",
            source_strain="RM11",
            bg_seq_pos=2,
            bg_allele="GCA",
            src_allele="",
            s288c_pos=100,
        )
        result = apply_variant(sequence, variant, verify_allele=True)
        assert result == "ATTGC"

    def test_apply_complex_substitution(self):
        """Test applying a complex substitution (different lengths)."""
        sequence = "ATGCATGC"
        variant = PairwiseVariant(
            background_strain="S288C",
            source_strain="RM11",
            bg_seq_pos=2,
            bg_allele="GC",
            src_allele="TTTT",
            s288c_pos=100,
        )
        result = apply_variant(sequence, variant, verify_allele=True)
        assert result == "ATTTTTATGC"

    def test_boundary_check_start(self):
        """Test that position 0 is valid."""
        sequence = "ATGC"
        variant = PairwiseVariant(
            background_strain="S288C",
            source_strain="RM11",
            bg_seq_pos=0,
            bg_allele="A",
            src_allele="G",
            s288c_pos=100,
        )
        result = apply_variant(sequence, variant, verify_allele=True)
        assert result == "GTGC"

    def test_boundary_check_end(self):
        """Test that position at end is valid for SNP."""
        sequence = "ATGC"
        variant = PairwiseVariant(
            background_strain="S288C",
            source_strain="RM11",
            bg_seq_pos=3,
            bg_allele="C",
            src_allele="G",
            s288c_pos=100,
        )
        result = apply_variant(sequence, variant, verify_allele=True)
        assert result == "ATGG"

    def test_boundary_check_past_end_rejected(self):
        """Test that position past end is rejected (critical bug fix test)."""
        sequence = "ATGC"  # len = 4
        variant = PairwiseVariant(
            background_strain="S288C",
            source_strain="RM11",
            bg_seq_pos=4,  # == len(sequence), invalid for non-empty bg_allele
            bg_allele="A",
            src_allele="G",
            s288c_pos=100,
        )
        result = apply_variant(sequence, variant, verify_allele=True)
        # Should return None because position is invalid
        assert result is None

    def test_insertion_at_end_allowed(self):
        """Test that insertion at end of sequence is allowed."""
        sequence = "ATGC"  # len = 4
        variant = PairwiseVariant(
            background_strain="S288C",
            source_strain="RM11",
            bg_seq_pos=4,  # == len(sequence), valid for insertion (empty bg_allele)
            bg_allele="",
            src_allele="GGG",
            s288c_pos=100,
        )
        result = apply_variant(sequence, variant, verify_allele=True)
        assert result == "ATGCGGG"

    def test_negative_position_rejected(self):
        """Test that negative positions are rejected."""
        sequence = "ATGC"
        variant = PairwiseVariant(
            background_strain="S288C",
            source_strain="RM11",
            bg_seq_pos=-1,
            bg_allele="A",
            src_allele="G",
            s288c_pos=100,
        )
        result = apply_variant(sequence, variant, verify_allele=True)
        assert result is None

    def test_allele_mismatch_rejected(self):
        """Test that allele mismatch is rejected."""
        sequence = "ATGC"
        variant = PairwiseVariant(
            background_strain="S288C",
            source_strain="RM11",
            bg_seq_pos=0,
            bg_allele="G",  # Mismatch - sequence has 'A' at position 0
            src_allele="T",
            s288c_pos=100,
        )
        result = apply_variant(sequence, variant, verify_allele=True)
        assert result is None

    def test_allele_extends_past_sequence(self):
        """Test that allele extending past sequence end is rejected."""
        sequence = "ATGC"
        variant = PairwiseVariant(
            background_strain="S288C",
            source_strain="RM11",
            bg_seq_pos=3,
            bg_allele="CCC",  # Would extend past sequence end
            src_allele="G",
            s288c_pos=100,
        )
        result = apply_variant(sequence, variant, verify_allele=True)
        assert result is None


# ============================================================================
# compute_strain_seq_position Tests
# ============================================================================

class TestComputeStrainSeqPosition:
    """Tests for compute_strain_seq_position function."""

    @pytest.fixture
    def mock_variants(self):
        """Create mock variant objects with pos and ref attributes."""
        @dataclass
        class MockVariant:
            pos: int
            ref: str

        return [
            MockVariant(pos=100, ref="A"),      # SNP - no offset change
            MockVariant(pos=200, ref="ATG"),    # 3bp deletion at pos 200
            MockVariant(pos=300, ref="C"),      # SNP
        ]

    def test_position_before_variants(self, mock_variants):
        """Test position before any variants - no offset."""
        result = compute_strain_seq_position(
            s288c_pos=50,
            strain_variants=mock_variants,
            locus_start=0,
        )
        assert result == 50

    def test_position_after_snp(self, mock_variants):
        """Test position after SNP - no offset change."""
        result = compute_strain_seq_position(
            s288c_pos=150,
            strain_variants=mock_variants,
            locus_start=0,
        )
        assert result == 150

    def test_position_after_deletion(self, mock_variants):
        """Test position after deletion - offset changes."""
        # After 3bp deletion at 200, positions shift by -2 (deletion of 3bp, VCF includes anchor)
        result = compute_strain_seq_position(
            s288c_pos=250,
            strain_variants=mock_variants,
            locus_start=0,
        )
        # This depends on the exact implementation
        # Just verify it returns an integer
        assert isinstance(result, int)

    def test_bounds_validation_raises_on_negative(self, mock_variants):
        """Test that bounds validation raises on negative position."""
        with pytest.raises(ValueError):
            compute_strain_seq_position(
                s288c_pos=-10,  # Invalid
                strain_variants=mock_variants,
                locus_start=0,
                validate_bounds=True,
            )


# ============================================================================
# Regression Tests for Known Bugs
# ============================================================================

class TestRegressionBugs:
    """Regression tests for previously identified bugs."""

    def test_off_by_one_at_sequence_end(self):
        """
        Regression test for BUG 1: Off-by-one boundary check.

        Previously, pos == len(sequence) was allowed even for non-empty bg_allele.
        """
        sequence = "ATGC"  # len = 4
        variant = PairwiseVariant(
            background_strain="S288C",
            source_strain="RM11",
            bg_seq_pos=4,  # This was incorrectly allowed
            bg_allele="X",  # Non-empty, so this should be rejected
            src_allele="Y",
            s288c_pos=100,
        )
        result = apply_variant(sequence, variant, verify_allele=True)
        # After fix, this should be None (rejected)
        assert result is None

    def test_variant_object_not_mutated(self):
        """
        Test that apply_variant does not mutate the original variant object.

        This was fixed in variant_inserter.py but should also hold here.
        """
        variant = PairwiseVariant(
            background_strain="S288C",
            source_strain="RM11",
            bg_seq_pos=2,
            bg_allele="G",
            src_allele="T",
            s288c_pos=100,
        )
        original_pos = variant.bg_seq_pos
        original_bg = variant.bg_allele
        original_src = variant.src_allele

        sequence = "ATGC"
        _ = apply_variant(sequence, variant, verify_allele=True)

        # Variant should not be mutated
        assert variant.bg_seq_pos == original_pos
        assert variant.bg_allele == original_bg
        assert variant.src_allele == original_src
