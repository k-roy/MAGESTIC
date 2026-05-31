#!/usr/bin/env python3
"""
Yeast Sublibrary Utilities
===========================

This module provides functions for working with yeast_sublibraries, which represent
nuclease-PAM combinations (transformation batches) in MAGESTIC screens.

Terminology (2026-02-14 update):
    - **yeast_sublibrary**: The nuclease-PAM combination at the transformation level.
      This is the trafo_ID (e.g., trafo_1, trafo_2) or individual yL (e.g., yL386, yL387).
      Examples: SpG_NGNG, SpG_NGNH, SpG_NGG, SpCas9_NGG

    - **guide_donor_oligo_subpool**: The oligo design pools with different SPS.
      These are the input oligo pools before transformation.

    - **step_1_plasmid_library**: The V guide-donor-bc0 libraries (e.g., V629, V631, V634).
      These are the plasmid libraries created from oligo subpools.

Key Concepts:
    - **PAM**: The DNA sequence the guide-donor oligo was DESIGNED for.
      Inferred from oligo_name prefix (e.g., SpG_Cas9 → NGNG/NGNH).

    - **Nuclease in Strain**: The actual nuclease protein in the yeast strain.
      Encoded in yKR ID/trafo_ID. CANNOT be inferred from oligo_name alone.

    - **yeast_sublibrary**: The nuclease-PAM combination that defines editing behavior.
      e.g., yL386 = SpG nuclease + NGNG PAM = SpG_NGNG
      e.g., yL391 = SpCas9 nuclease + NGG PAM = SpCas9_NGG

See PAM_NUCLEASE_REFERENCE.md for detailed documentation.

Author: Kevin R. Roy
"""

from dataclasses import dataclass
from typing import Dict, Optional, Tuple
import re
import logging

logger = logging.getLogger(__name__)


# ============================================================================
# Strain Nuclease Mappings
# ============================================================================

# yKR strain ID → nuclease protein mapping
# This is the CANONICAL mapping for strain nucleases
STRAIN_NUCLEASE_MAP = {
    # SpCas9 strains
    "yKR888": "WT SpCas9",
    "yKR1036": "WT SpCas9",

    # SpG strains (SpCas9-NG variant with broadened PAM)
    "yKR1033": "SpG",

    # SpG variants
    "yKR1034": "eSpCas9(1.1)-SpG",  # Enhanced specificity SpG
    "yKR1035": "Sniper2L SpG",      # Sniper2L variant

    # LbCas12a strains
    "yKR900": "WT LbCas12a",
    "yKR901": "impLbCas12a",        # Improved LbCas12a
}

# yL library ID → (yKR strain, PAM, notes) mapping
# Canonical yeast_sublibrary definitions for QTL screens
YEAST_SUBLIBRARY_MAP = {
    # SpCas9 NGG libraries
    "yL274": ("yKR888", "NGG", "SpCas9 library"),
    "yL275": ("yKR888", "NGG", "SpCas9 library"),
    "yL276": ("yKR888", "NGG", "SpCas9 library"),
    "yL391": ("yKR1036", "NGG", "NGG QTL library (rep_1)"),
    "yL397": ("yKR1036", "NGG", "NGG QTL library (rep_2)"),

    # SpG NGNG libraries
    "yL386": ("yKR1033", "NGNG", "NGNG QTL library"),
    "yL392": ("yKR1033", "NGNG", "NGNG QTL library (rep_2)"),

    # SpG NGNH libraries
    "yL387": ("yKR1033", "NGNH", "NGNH QTL library (rep_1)"),
    "yL393": ("yKR1033", "NGNH", "NGNH QTL library (rep_2)"),

    # LbCas12a TTTV libraries
    "yL278": ("yKR900", "TTTV", "WT LbCas12a TTTV library"),
    "yL279": ("yKR900", "TTTV", "WT LbCas12a TTTV library"),
    "yL398": ("yKR1033", "TTTV", "LbCas12a QTL library - non-editing control"),
    "yL401": ("yKR1036", "TTTV", "LbCas12a QTL library - non-editing control"),

    # impLbCas12a libraries
    "yL281": ("yKR901", "TTTV", "WT LbCas12a TTTV library, impLbCas12a nuclease"),
    "yL282": ("yKR901", "non-TTTV", "impLbCas12a library"),
    "yL283": ("yKR901", "TTTV", "WT LbCas12a TTTV library"),
    "yL284": ("yKR901", "non-TTTV", "impLbCas12a library"),

    # Non-editing controls (PAM/nuclease mismatch)
    "yL277": ("yKR900", "NGG", "SpCas9 library in LbCas12a strain - non-editing control"),
    "yL280": ("yKR901", "NGG", "SpCas9 library in impLbCas12a strain - non-editing control"),

    # Combined library (yL406)
    "yL406": ("yKR1033", "mixed", "Combined QTL library (NGNG:NGNH:NGG:TTTV:non-TTTV = 3:3:3:1:1)"),
}

# Nuclease categories for aggregation
NUCLEASE_CATEGORIES = {
    "SpCas9": ["WT SpCas9"],
    "SpG": ["SpG", "eSpCas9(1.1)-SpG", "Sniper2L SpG"],
    "LbCas12a": ["WT LbCas12a", "impLbCas12a"],
}

# Nuclease-PAM prefixes used in variant/oligo names
# These are the combinations that appear in oligo_name prefixes like "SpG_NGG_YOR202W_..."
# Used for:
#   1. Extracting nuclease type from variant IDs
#   2. Stripping prefixes to get variant identity across nuclease libraries
#   3. Regex pattern building for nuclease extraction
NUCLEASE_PAM_PREFIXES = [
    "SpG_NGNG",    # SpG nuclease with NGNG PAM
    "SpG_NGNH",    # SpG nuclease with NGNH PAM (lower editing efficiency)
    "SpG_NGG",     # SpG nuclease with NGG PAM
    "SpCas9_NGG",  # WT SpCas9 nuclease with NGG PAM
]

# Pre-built regex pattern for nuclease extraction (performance optimization)
NUCLEASE_PAM_PATTERN = r'(' + '|'.join(NUCLEASE_PAM_PREFIXES) + r')'


# ============================================================================
# Data Classes
# ============================================================================

@dataclass
class YeastSublibraryInfo:
    """Information about a yeast_sublibrary (nuclease-PAM combination).

    A yeast_sublibrary represents a transformation batch with a specific nuclease
    and PAM combination.

    Attributes:
        trafo_id: The yL library ID or trafo_ID (e.g., "yL386", "trafo_1")
        strain_id: The yKR strain ID (e.g., "yKR1033")
        nuclease: The nuclease protein in the strain (e.g., "SpG")
        pam: The PAM the oligo library was designed for (e.g., "NGNG")
        is_editing: Whether this is an editing combination (nuclease matches PAM)
        notes: Additional notes about this yeast_sublibrary
    """
    trafo_id: str
    strain_id: str
    nuclease: str
    pam: str
    is_editing: bool
    notes: str = ""

    @property
    def yeast_sublibrary_name(self) -> str:
        """Return a descriptive yeast_sublibrary name like 'SpG_NGNG'."""
        return f"{self.nuclease}_{self.pam}"

    @property
    def is_non_editing_control(self) -> bool:
        """Return True if this is a non-editing control (PAM/nuclease mismatch)."""
        return not self.is_editing


# ============================================================================
# PAM Inference (from oligo_name)
# ============================================================================

def infer_pam_type(oligo_name: str) -> str:
    """Infer PAM type from oligo name prefix.

    This function ONLY tells you what PAM the oligo was DESIGNED for.
    It does NOT tell you the nuclease in the strain.

    Args:
        oligo_name: The oligo name (e.g., "SpG_Cas9_NGNG_YDL195W_...")

    Returns:
        PAM type: "NGNG", "NGNH", "NGG", "TTTV", "non-TTTV", or "unknown"

    Examples:
        >>> infer_pam_type("SpG_Cas9_NGNG_YDL195W_A123T")
        "NGNG"
        >>> infer_pam_type("SpCas9_NGG_YBR121C_stop_30")
        "NGG"
        >>> infer_pam_type("LbCas12a_TTTV_YDL195W_A123T")
        "TTTV"
    """
    oligo_upper = oligo_name.upper()

    # SpG oligos - check NGNG vs NGNH
    if oligo_upper.startswith("SPG_CAS9_NGNG") or "_NGNG_" in oligo_upper:
        return "NGNG"
    if oligo_upper.startswith("SPG_CAS9_NGNH") or "_NGNH_" in oligo_upper:
        return "NGNH"
    if oligo_upper.startswith("SPG_CAS9") or oligo_upper.startswith("SPG_"):
        # Default SpG to NGNG if not specified
        return "NGNG"

    # SpCas9 oligos - always NGG
    if oligo_upper.startswith("SPCAS9_"):
        return "NGG"

    # LbCas12a oligos - TTTV
    if oligo_upper.startswith("LBCAS12A_"):
        return "TTTV"

    # impLbCas12a oligos - non-TTTV (relaxed PAM)
    if oligo_upper.startswith("IMPLBCAS12A_"):
        return "non-TTTV"

    # Fallback: try to parse from explicit PAM in name
    if "_NGG_" in oligo_upper:
        return "NGG"
    if "_NGNG_" in oligo_upper:
        return "NGNG"
    if "_NGNH_" in oligo_upper:
        return "NGNH"
    if "_TTTV_" in oligo_upper:
        return "TTTV"

    logger.warning(f"Could not infer PAM type from oligo name: {oligo_name}")
    return "unknown"


# ============================================================================
# Nuclease Lookup (from trafo_ID or strain ID)
# ============================================================================

def get_strain_nuclease(strain_or_trafo_id: str) -> Optional[str]:
    """Get the nuclease protein in a strain from its ID.

    This function looks up the nuclease from either:
    - yKR strain ID (e.g., "yKR1033" → "SpG")
    - yL library ID (e.g., "yL386" → "SpG" via yKR1033)

    Args:
        strain_or_trafo_id: Either a yKR strain ID or yL library ID

    Returns:
        Nuclease name (e.g., "SpG", "WT SpCas9", "WT LbCas12a")
        or None if not found in mapping.

    Examples:
        >>> get_strain_nuclease("yKR1033")
        "SpG"
        >>> get_strain_nuclease("yL386")
        "SpG"
        >>> get_strain_nuclease("yKR1036")
        "WT SpCas9"
    """
    # Normalize input
    id_upper = strain_or_trafo_id.strip().upper()

    # Try direct strain lookup (yKR...)
    for strain_id, nuclease in STRAIN_NUCLEASE_MAP.items():
        if strain_id.upper() == id_upper:
            return nuclease

    # Try sublibrary lookup (yL...)
    for lib_id, (strain_id, pam, notes) in YEAST_SUBLIBRARY_MAP.items():
        if lib_id.upper() == id_upper:
            return STRAIN_NUCLEASE_MAP.get(strain_id)

    logger.warning(f"Unknown strain/trafo ID: {strain_or_trafo_id}")
    return None


def get_nuclease_category(nuclease: str) -> Optional[str]:
    """Get the nuclease category (SpCas9, SpG, or LbCas12a).

    Groups nuclease variants into their major categories.

    Args:
        nuclease: Nuclease name (e.g., "SpG", "eSpCas9(1.1)-SpG")

    Returns:
        Category: "SpCas9", "SpG", or "LbCas12a", or None if unknown.
    """
    for category, variants in NUCLEASE_CATEGORIES.items():
        if nuclease in variants:
            return category
    return None


# ============================================================================
# Sublibrary Parsing
# ============================================================================

def parse_yeast_sublibrary_info(trafo_id: str) -> Optional[YeastSublibraryInfo]:
    """Parse full yeast_sublibrary information from a trafo_ID (yL library ID).

    Returns comprehensive information about the yeast_sublibrary including:
    - The strain nuclease
    - The PAM the oligos were designed for
    - Whether this is an editing or non-editing combination

    Args:
        trafo_id: The yL library ID (e.g., "yL386")

    Returns:
        YeastSublibraryInfo dataclass or None if not found.

    Examples:
        >>> info = parse_yeast_sublibrary_info("yL386")
        >>> info.nuclease
        "SpG"
        >>> info.pam
        "NGNG"
        >>> info.is_editing
        True
        >>> info.yeast_sublibrary_name
        "SpG_NGNG"
    """
    # Normalize
    trafo_upper = trafo_id.strip().upper()

    # Look up in sublibrary map
    for lib_id, (strain_id, pam, notes) in YEAST_SUBLIBRARY_MAP.items():
        if lib_id.upper() == trafo_upper:
            nuclease = STRAIN_NUCLEASE_MAP.get(strain_id)
            if nuclease is None:
                logger.warning(f"Strain {strain_id} not in nuclease map")
                return None

            # Determine if editing combination
            is_editing = _is_compatible_nuclease_pam(nuclease, pam)

            return YeastSublibraryInfo(
                trafo_id=lib_id,
                strain_id=strain_id,
                nuclease=nuclease,
                pam=pam,
                is_editing=is_editing,
                notes=notes,
            )

    logger.warning(f"Unknown trafo_ID: {trafo_id}")
    return None


def _is_compatible_nuclease_pam(nuclease: str, pam: str) -> bool:
    """Determine if a nuclease-PAM combination enables editing.

    Editing occurs when the nuclease can recognize the PAM:
    - SpCas9 → NGG only
    - SpG → NGG, NGNG, NGNH (broadened PAM)
    - LbCas12a → TTTV
    - impLbCas12a → TTTV, non-TTTV
    """
    nuclease_cat = get_nuclease_category(nuclease)

    if nuclease_cat == "SpCas9":
        return pam == "NGG"

    if nuclease_cat == "SpG":
        # SpG can edit NGG, NGNG, NGNH (all NG PAMs)
        return pam in ["NGG", "NGNG", "NGNH"]

    if nuclease_cat == "LbCas12a":
        if nuclease == "impLbCas12a":
            # Improved LbCas12a has relaxed PAM
            return pam in ["TTTV", "non-TTTV"]
        else:
            # WT LbCas12a requires TTTV
            return pam == "TTTV"

    # Unknown combination
    return False


# ============================================================================
# Sublibrary Extraction from oligo_name + trafo_ID
# ============================================================================

def get_yeast_sublibrary_from_context(
    oligo_name: str,
    trafo_id: Optional[str] = None,
) -> Tuple[str, str, bool]:
    """Get yeast_sublibrary information from oligo name and optional trafo_ID.

    This combines PAM inference from oligo_name with nuclease lookup from
    trafo_ID to determine the full yeast_sublibrary context.

    Args:
        oligo_name: The oligo name (for PAM inference)
        trafo_id: Optional yL library ID (for nuclease lookup)

    Returns:
        Tuple of (pam, nuclease, is_editing)
        - pam: The PAM type from oligo_name
        - nuclease: The nuclease from trafo_id (or "unknown")
        - is_editing: Whether this is an editing combination

    Examples:
        >>> get_yeast_sublibrary_from_context("SpG_Cas9_NGNG_YDL195W_A123T", "yL386")
        ("NGNG", "SpG", True)
        >>> get_yeast_sublibrary_from_context("LbCas12a_TTTV_YDL195W_A123T", "yL398")
        ("TTTV", "SpG", False)  # Non-editing control!
    """
    pam = infer_pam_type(oligo_name)

    if trafo_id:
        nuclease = get_strain_nuclease(trafo_id) or "unknown"
    else:
        nuclease = "unknown"

    if nuclease != "unknown":
        is_editing = _is_compatible_nuclease_pam(nuclease, pam)
    else:
        # Without trafo_id, assume PAM matches nuclease (standard case)
        is_editing = True

    return (pam, nuclease, is_editing)


# ============================================================================
# Utility Functions
# ============================================================================

def list_all_yeast_sublibraries() -> Dict[str, YeastSublibraryInfo]:
    """Return information for all known yeast_sublibraries.

    Returns:
        Dictionary mapping trafo_ID to YeastSublibraryInfo.
    """
    result = {}
    for lib_id in YEAST_SUBLIBRARY_MAP:
        info = parse_yeast_sublibrary_info(lib_id)
        if info:
            result[lib_id] = info
    return result


def get_editing_yeast_sublibraries() -> Dict[str, YeastSublibraryInfo]:
    """Return only editing yeast_sublibraries (PAM matches nuclease).

    Non-editing controls (e.g., LbCas12a oligos in SpG strain) are excluded.
    """
    return {k: v for k, v in list_all_yeast_sublibraries().items() if v.is_editing}


def get_non_editing_controls() -> Dict[str, YeastSublibraryInfo]:
    """Return only non-editing control yeast_sublibraries.

    These are intentional PAM/nuclease mismatches used as controls.
    """
    return {k: v for k, v in list_all_yeast_sublibraries().items() if not v.is_editing}


def validate_aggregation_group(trafo_ids: list) -> bool:
    """Validate that all trafo_IDs in a group can be aggregated together.

    Aggregation is only valid when ALL trafo_IDs are editing yeast_sublibraries
    OR ALL are non-editing controls. NEVER mix editing and non-editing.

    Args:
        trafo_ids: List of yL library IDs or trafo_IDs to check

    Returns:
        True if the group is valid for aggregation.

    Raises:
        ValueError: If the group mixes editing and non-editing yeast_sublibraries.
    """
    infos = [parse_yeast_sublibrary_info(tid) for tid in trafo_ids]
    infos = [i for i in infos if i is not None]

    if not infos:
        logger.warning("No valid sublibrary info found for aggregation check")
        return False

    editing_count = sum(1 for i in infos if i.is_editing)
    non_editing_count = len(infos) - editing_count

    if editing_count > 0 and non_editing_count > 0:
        editing_libs = [i.trafo_id for i in infos if i.is_editing]
        non_editing_libs = [i.trafo_id for i in infos if not i.is_editing]
        raise ValueError(
            f"Cannot aggregate editing libraries {editing_libs} with "
            f"non-editing controls {non_editing_libs}. "
            "Filter out non-editing controls before aggregation."
        )

    return True
