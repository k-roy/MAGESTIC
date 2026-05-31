"""
Shared variant annotation utilities.

This module provides common functions for classifying variant effects
(frameshift, missense, synonymous, etc.) used across multiple pipelines:
- bc1-donor-bc0 pipeline (unmatched donor annotation)
- variant_annotation_harmonization pipeline
- bc1_from_screen pipeline

Author: Kevin Roy / Claude
Date: 2026-02-02
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np


# =============================================================================
# Variant Type Classification
# =============================================================================

def classify_variant_type(
    ref_aas: Optional[str],
    alt_aas: Optional[str],
    aa_nums: Optional[str] = None,
    indel_length: Optional[int] = None,
    variant_pos: Optional[int] = None,
    cds_start: Optional[int] = None,
    cds_end: Optional[int] = None,
) -> str:
    """
    Classify the coding effect of a variant with granular details.

    This function uses a priority order to determine variant type:
    1. Handle indels (frameshift vs in-frame)
    2. Handle stop codons (stop_gain, stop_loss)
    3. Handle start codon loss
    4. Synonymous vs missense SNPs

    Args:
        ref_aas: Reference amino acid(s), comma-separated if multiple
        alt_aas: Alternate amino acid(s), comma-separated if multiple
        aa_nums: Amino acid position(s), comma-separated if multiple
        indel_length: Net indel length (positive=insertion, negative=deletion)
        variant_pos: Genomic position of variant
        cds_start: CDS start coordinate
        cds_end: CDS end coordinate

    Returns:
        Variant type: one of:
        - 'frameshift': indel not divisible by 3
        - 'in-frame_indel': indel divisible by 3
        - 'stop_gain': new stop codon introduced
        - 'stop_loss': stop codon removed
        - 'start_loss': start codon (M) removed at position 1
        - 'missense': amino acid change
        - 'synonymous': no amino acid change
        - 'non_coding': variant outside CDS or no AA info
        - 'complete_CDS_deletion': entire CDS deleted
    """
    # Handle indels FIRST (before checking variant_position)
    if pd.notna(indel_length) and indel_length != 0:
        # Check if we have CDS coordinates
        if pd.isna(cds_start) or pd.isna(cds_end):
            return "non_coding"

        # For indels, if we have ref_aas/alt_aas, it's coding
        if pd.notna(ref_aas) or pd.notna(alt_aas):
            return "frameshift" if abs(indel_length) % 3 != 0 else "in-frame_indel"
        else:
            return "non_coding"

    # Check if SNP variant is within CDS boundaries
    if pd.notna(variant_pos) and pd.notna(cds_start) and pd.notna(cds_end):
        if not (cds_start <= variant_pos <= cds_end):
            return "non_coding"

    # Handle SNPs and amino acid changes
    if pd.isna(ref_aas) or pd.isna(alt_aas):
        return "non_coding"

    ref_str = str(ref_aas)
    alt_str = str(alt_aas)

    # Check for stop gain (alt becomes stop codon)
    if "*" in alt_str and "*" not in ref_str:
        return "stop_gain"

    # Check for stop loss (ref has stop, alt doesn't)
    if "*" in ref_str and "*" not in alt_str:
        return "stop_loss"

    # Check for start loss (position 1, ref is M, alt is not M)
    if pd.notna(aa_nums):
        try:
            first_pos = int(str(aa_nums).split(",")[0])
            if first_pos == 1 and ref_str.startswith("M") and not alt_str.startswith("M"):
                return "start_loss"
        except (ValueError, IndexError):
            pass

    # Synonymous vs missense
    if ref_str == alt_str:
        return "synonymous"
    else:
        return "missense"


def classify_variant_type_simple(original_aa: str, new_aa: str) -> str:
    """
    Simple variant classification based on amino acid change only.

    This is a simplified version for cases where only AA info is available
    (no indel info or CDS coordinates).

    Args:
        original_aa: Original amino acid (single letter)
        new_aa: New amino acid (single letter or *)

    Returns:
        Variant type: 'nonsense', 'missense', 'synonymous', or 'unknown'
    """
    if not original_aa or not new_aa:
        return "unknown"

    if new_aa == "*":
        return "nonsense"  # Note: 'stop_gain' is the more specific term
    elif original_aa == new_aa:
        return "synonymous"
    else:
        return "missense"


# =============================================================================
# Indel Length Calculation
# =============================================================================

def compute_indel_length(individual_variants: str) -> Optional[int]:
    """
    Compute net indel length from individual_variants string.

    Format: position_ref_alt_INDEL or position_ref_alt (for SNPs)
    Example: "123_A_ATG_INDEL,456_GC_G_INDEL" -> +2 (net insertion of 2bp)

    Args:
        individual_variants: Comma-separated variant strings

    Returns:
        Net indel length (positive for insertions, negative for deletions),
        or None for SNPs only
    """
    if pd.isna(individual_variants):
        return None

    # Check if any variant is an indel
    if "INDEL" not in str(individual_variants):
        return None

    total_length = 0
    for var in str(individual_variants).split(","):
        if "_" not in var or not any(c.isdigit() for c in var):
            continue

        parts = var.split("_")
        if len(parts) < 3:
            continue

        try:
            ref_allele = parts[1]
            alt_allele = parts[2]

            # Calculate net length change
            length_change = len(alt_allele) - len(ref_allele)
            total_length += length_change
        except (ValueError, IndexError):
            continue

    return total_length if total_length != 0 else None


def compute_indel_length_within_cds(
    var_str: str,
    cds_start: int,
    cds_end: int
) -> int:
    """
    Sum the net length change of indels fully or partially contained in a CDS.

    For each variant in the comma-separated var_str, compute var_start and var_end
    (var_end = var_start + len(ref_allele) - 1). If positions are within or overlap
    the CDS boundaries, calculate the contribution to length change.

    Args:
        var_str: Comma-separated variant string
        cds_start: CDS start coordinate
        cds_end: CDS end coordinate

    Returns:
        Net indel length change within CDS
    """
    if pd.isna(var_str) or pd.isna(cds_start) or pd.isna(cds_end) or cds_start > cds_end:
        return 0

    total_delta = 0
    for var in str(var_str).split(","):
        if "_" not in var or not any(c.isdigit() for c in var):
            continue
        parts = var.split("_")
        if len(parts) < 3:
            continue

        try:
            var_start = int(parts[0])
        except ValueError:
            continue

        ref_allele = parts[1]
        alt_allele = parts[2]
        var_end = var_start + len(ref_allele) - 1

        # Fully contained within CDS
        if cds_start <= var_start and var_end <= cds_end:
            total_delta += len(alt_allele) - len(ref_allele)
            continue

        # Partially overlapping CDS boundaries
        overlap_start = max(var_start, cds_start)
        overlap_end = min(var_end, cds_end)
        if overlap_start <= overlap_end:
            ref_overlap = overlap_end - overlap_start + 1
            left_trim = max(cds_start - var_start, 0)
            right_trim = max(var_end - cds_end, 0)
            alt_overlap = max(len(alt_allele) - left_trim - right_trim, 0)
            total_delta += alt_overlap - ref_overlap

    return total_delta


# =============================================================================
# Variant Position Extraction
# =============================================================================

def extract_variant_positions(var_str: str) -> Tuple[Optional[int], Optional[int]]:
    """
    Extract start and end positions from variant string.

    Args:
        var_str: Comma-separated variant string (e.g., "123_A_G,456_GC_G_INDEL")

    Returns:
        Tuple of (min_position, max_position) or (None, None) if no valid positions
    """
    if pd.isna(var_str):
        return None, None
    try:
        positions = []
        for var in str(var_str).split(","):
            if "_" in var and any(c.isdigit() for c in var):
                parts = var.split("_")
                var_start = int(parts[0])
                positions.append(var_start)
                if parts[-1] == "INDEL":
                    positions.append(var_start + len(parts[1]) - 1)
        return (min(positions), max(positions)) if positions else (None, None)
    except (ValueError, IndexError):
        return None, None


# =============================================================================
# Amino Acid Change Formatting
# =============================================================================

def format_aa_change(
    ref_aas: Optional[str],
    aa_nums: Optional[str],
    alt_aas: Optional[str]
) -> Optional[str]:
    """
    Format amino acid change as human-readable string.

    Args:
        ref_aas: Reference amino acids (comma-separated)
        aa_nums: Amino acid positions (comma-separated)
        alt_aas: Alternate amino acids (comma-separated)

    Returns:
        Formatted string like "A1V_G2G_K3L" or None
    """
    if pd.isna(ref_aas) or pd.isna(aa_nums) or pd.isna(alt_aas):
        return None

    try:
        return "_".join(
            f"{r}{n}{a}"
            for r, n, a in zip(
                str(ref_aas).split(","),
                str(aa_nums).split(","),
                str(alt_aas).split(",")
            )
        )
    except Exception:
        return None


def extract_missense_changes(aa_change: str) -> Optional[str]:
    """
    Extract only missense (non-synonymous) changes from aa_change string.

    This extracts amino acid changes where ref != alt, excluding:
    - Synonymous changes (e.g., P28P where ref == alt)

    Example:
        "P28P_R29R_G30*" → "G30*"  (only G30* is missense/stop)
        "A1V_G2G_K3L" → "A1V_K3L"

    Args:
        aa_change: Underscore-separated amino acid change string

    Returns:
        Underscore-joined missense changes only, or None if no missense changes
    """
    if pd.isna(aa_change) or str(aa_change).strip() == "":
        return None

    missense_changes = []
    for change in str(aa_change).split("_"):
        change = change.strip()
        if not change or len(change) < 3:
            continue

        ref_aa = change[0]
        alt_aa = change[-1]

        # Include if ref != alt (missense, stop_gain, stop_loss, etc.)
        if ref_aa != alt_aa:
            missense_changes.append(change)

    if not missense_changes:
        return None

    return "_".join(missense_changes)


# =============================================================================
# CDS Overlap Checking
# =============================================================================

def variant_overlaps_cds(
    var_start: Optional[int],
    var_end: Optional[int],
    cds_start: int,
    cds_end: int
) -> bool:
    """
    Check if a variant overlaps with a CDS region.

    Args:
        var_start: Variant start position
        var_end: Variant end position
        cds_start: CDS start coordinate
        cds_end: CDS end coordinate

    Returns:
        True if variant overlaps CDS, False otherwise
    """
    if pd.isna(var_start) or pd.isna(var_end):
        return False

    # Overlap exists if variant doesn't end before CDS starts
    # and doesn't start after CDS ends
    return var_end >= cds_start and var_start <= cds_end


# =============================================================================
# Data Classes for Annotation Results
# =============================================================================

@dataclass
class VariantAnnotation:
    """
    Standardized variant annotation fields.

    This dataclass defines the common annotation fields used across
    pipelines for consistent output formatting.
    """
    # Core identifiers
    oligo_name: str = ""
    bc1: str = ""

    # Gene information
    systematic_gene_name: str = ""
    common_gene_name: str = ""
    chromosome: str = ""
    strand: str = ""

    # Variant description
    individual_variants: str = ""
    aa_change: str = ""
    shared_aa_missense_change: str = ""

    # Variant classification
    codon_variant_type: str = ""  # frameshift, missense, synonymous, etc.
    indel_length: Optional[int] = None
    variant_overlaps_CDS: bool = False

    # Sub-variant relationship
    is_sub_variant: bool = False
    is_incomplete_frameshift: bool = False

    def to_dict(self) -> Dict:
        """Convert to dictionary for DataFrame creation."""
        return {
            "oligo_name": self.oligo_name,
            "bc1": self.bc1,
            "systematic_gene_name": self.systematic_gene_name,
            "common_gene_name": self.common_gene_name,
            "chromosome": self.chromosome,
            "strand": self.strand,
            "individual_variants": self.individual_variants,
            "aa_change": self.aa_change,
            "shared_aa_missense_change": self.shared_aa_missense_change,
            "codon_variant_type": self.codon_variant_type,
            "indel_length": self.indel_length,
            "variant_overlaps_CDS": self.variant_overlaps_CDS,
            "is_sub_variant": self.is_sub_variant,
            "is_incomplete_frameshift": self.is_incomplete_frameshift,
        }


# =============================================================================
# Harmonized Column Names
# =============================================================================

# Standard column names for variant annotation outputs
# Use these consistently across all pipelines
VARIANT_ANNOTATION_COLUMNS = [
    # Identifiers
    "oligo_name",
    "bc1",

    # Gene info
    "systematic_gene_name_locus",
    "common_gene_name_locus",
    "chromosome",
    "strand",

    # Variant description
    "individual_variants",
    "aa_change",
    "shared_aa_missense_change",

    # Classification
    "codon_variant_type",  # frameshift, missense, synonymous, stop_gain, etc.
    "indel_length",
    "variant_overlaps_CDS",

    # Sub-variant info
    "is_sub_variant",
    "is_incomplete_frameshift",
]

# Valid values for codon_variant_type
CODON_VARIANT_TYPES = [
    "frameshift",
    "in-frame_indel",
    "stop_gain",
    "stop_loss",
    "start_loss",
    "missense",
    "synonymous",
    "non_coding",
    "complete_CDS_deletion",
]
