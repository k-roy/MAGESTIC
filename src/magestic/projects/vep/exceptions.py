"""
Custom exception hierarchy for VEP pipeline.

Provides specific exception classes for different error conditions,
enabling better error handling and debugging.

Author: Kevin R. Roy
Date: 2026-02-23

Exception Hierarchy:
    VEPError (base)
    ├── SequenceError
    │   ├── InvalidSequenceError
    │   ├── SequenceLengthError
    │   └── SequenceNotFoundError
    ├── VariantError
    │   ├── InvalidVariantError
    │   ├── VariantPositionError
    │   ├── AlleleMismatchError
    │   └── LinkedVariantError
    ├── CoordinateError
    │   ├── OffsetTableError
    │   ├── PositionOutOfBoundsError
    │   └── CoordinateConversionError
    ├── TranslationError
    │   ├── FrameshiftError
    │   ├── StartCodonError
    │   └── StopCodonError
    ├── ModelError
    │   ├── ModelNotLoadedError
    │   ├── ModelInferenceError
    │   └── ModelConfigError
    └── FileError
        ├── FastaParseError
        ├── GFFParseError
        └── VCFParseError
"""

from typing import Optional, Any


# ============================================================================
# Base Exception
# ============================================================================

class VEPError(Exception):
    """Base exception for all VEP pipeline errors."""

    def __init__(self, message: str, context: Optional[dict] = None):
        super().__init__(message)
        self.message = message
        self.context = context or {}

    def __str__(self):
        if self.context:
            context_str = ", ".join(f"{k}={v}" for k, v in self.context.items())
            return f"{self.message} [{context_str}]"
        return self.message


# ============================================================================
# Sequence Errors
# ============================================================================

class SequenceError(VEPError):
    """Base class for sequence-related errors."""
    pass


class InvalidSequenceError(SequenceError):
    """Raised when a sequence contains invalid characters."""

    def __init__(self, message: str, sequence: str = None, invalid_chars: str = None):
        context = {}
        if sequence:
            context["seq_len"] = len(sequence)
            context["seq_preview"] = sequence[:50] + "..." if len(sequence) > 50 else sequence
        if invalid_chars:
            context["invalid_chars"] = invalid_chars
        super().__init__(message, context)


class SequenceLengthError(SequenceError):
    """Raised when sequence length doesn't meet requirements."""

    def __init__(self, message: str, actual_length: int, expected_length: int = None,
                 min_length: int = None, max_length: int = None):
        context = {"actual_length": actual_length}
        if expected_length is not None:
            context["expected_length"] = expected_length
        if min_length is not None:
            context["min_length"] = min_length
        if max_length is not None:
            context["max_length"] = max_length
        super().__init__(message, context)


class SequenceNotFoundError(SequenceError):
    """Raised when a requested sequence is not found."""

    def __init__(self, message: str, sequence_id: str = None, source: str = None):
        context = {}
        if sequence_id:
            context["sequence_id"] = sequence_id
        if source:
            context["source"] = source
        super().__init__(message, context)


# ============================================================================
# Variant Errors
# ============================================================================

class VariantError(VEPError):
    """Base class for variant-related errors."""
    pass


class InvalidVariantError(VariantError):
    """Raised when a variant has invalid properties."""

    def __init__(self, message: str, variant: Any = None):
        context = {}
        if variant:
            context["variant"] = str(variant)
        super().__init__(message, context)


class VariantPositionError(VariantError):
    """Raised when a variant position is invalid."""

    def __init__(self, message: str, position: int, seq_length: int = None,
                 variant: Any = None):
        context = {"position": position}
        if seq_length is not None:
            context["seq_length"] = seq_length
        if variant:
            context["variant"] = str(variant)
        super().__init__(message, context)


class AlleleMismatchError(VariantError):
    """Raised when expected allele doesn't match sequence."""

    def __init__(self, message: str, expected_allele: str, actual_allele: str,
                 position: int = None):
        context = {
            "expected_allele": expected_allele,
            "actual_allele": actual_allele,
        }
        if position is not None:
            context["position"] = position
        super().__init__(message, context)


class LinkedVariantError(VariantError):
    """Raised for errors in linked variant detection/handling."""

    def __init__(self, message: str, variants: list = None, net_change: int = None):
        context = {}
        if variants:
            context["n_variants"] = len(variants)
        if net_change is not None:
            context["net_change"] = net_change
        super().__init__(message, context)


# ============================================================================
# Coordinate Errors
# ============================================================================

class CoordinateError(VEPError):
    """Base class for coordinate-related errors."""
    pass


class OffsetTableError(CoordinateError):
    """Raised for offset table loading/parsing errors."""

    def __init__(self, message: str, strain: str = None, chrom: str = None,
                 file_path: str = None):
        context = {}
        if strain:
            context["strain"] = strain
        if chrom:
            context["chrom"] = chrom
        if file_path:
            context["file_path"] = file_path
        super().__init__(message, context)


class PositionOutOfBoundsError(CoordinateError):
    """Raised when a position is outside valid bounds."""

    def __init__(self, message: str, position: int, bounds: tuple = None,
                 coordinate_system: str = None):
        context = {"position": position}
        if bounds:
            context["bounds"] = f"[{bounds[0]}, {bounds[1]}]"
        if coordinate_system:
            context["coordinate_system"] = coordinate_system
        super().__init__(message, context)


class CoordinateConversionError(CoordinateError):
    """Raised when coordinate conversion fails."""

    def __init__(self, message: str, source_pos: int = None, source_system: str = None,
                 target_system: str = None):
        context = {}
        if source_pos is not None:
            context["source_pos"] = source_pos
        if source_system:
            context["source_system"] = source_system
        if target_system:
            context["target_system"] = target_system
        super().__init__(message, context)


# ============================================================================
# Translation Errors
# ============================================================================

class TranslationError(VEPError):
    """Base class for translation-related errors."""
    pass


class FrameshiftError(TranslationError):
    """Raised for frameshift-related issues."""

    def __init__(self, message: str, cds_length: int = None, expected_mod: int = 0,
                 actual_mod: int = None, orf: str = None):
        context = {}
        if cds_length is not None:
            context["cds_length"] = cds_length
        if actual_mod is not None:
            context["length_mod_3"] = actual_mod
        if orf:
            context["orf"] = orf
        super().__init__(message, context)


class StartCodonError(TranslationError):
    """Raised when start codon is missing or mutated."""

    def __init__(self, message: str, actual_codon: str = None, orf: str = None):
        context = {}
        if actual_codon:
            context["actual_codon"] = actual_codon
        if orf:
            context["orf"] = orf
        super().__init__(message, context)


class StopCodonError(TranslationError):
    """Raised for stop codon-related issues."""

    def __init__(self, message: str, expected_stop: bool = True, orf: str = None):
        context = {"expected_stop": expected_stop}
        if orf:
            context["orf"] = orf
        super().__init__(message, context)


# ============================================================================
# Model Errors
# ============================================================================

class ModelError(VEPError):
    """Base class for model-related errors."""
    pass


class ModelNotLoadedError(ModelError):
    """Raised when trying to use a model that hasn't been loaded."""

    def __init__(self, message: str, model_name: str = None):
        context = {}
        if model_name:
            context["model_name"] = model_name
        super().__init__(message, context)


class ModelInferenceError(ModelError):
    """Raised when model inference fails."""

    def __init__(self, message: str, model_name: str = None, input_shape: tuple = None):
        context = {}
        if model_name:
            context["model_name"] = model_name
        if input_shape:
            context["input_shape"] = str(input_shape)
        super().__init__(message, context)


class ModelConfigError(ModelError):
    """Raised for model configuration errors."""

    def __init__(self, message: str, config_path: str = None, missing_key: str = None):
        context = {}
        if config_path:
            context["config_path"] = config_path
        if missing_key:
            context["missing_key"] = missing_key
        super().__init__(message, context)


# ============================================================================
# File Parsing Errors
# ============================================================================

class FileError(VEPError):
    """Base class for file parsing errors."""
    pass


class FastaParseError(FileError):
    """Raised for FASTA file parsing errors."""

    def __init__(self, message: str, file_path: str = None, line_number: int = None):
        context = {}
        if file_path:
            context["file_path"] = file_path
        if line_number:
            context["line_number"] = line_number
        super().__init__(message, context)


class GFFParseError(FileError):
    """Raised for GFF file parsing errors."""

    def __init__(self, message: str, file_path: str = None, line_number: int = None,
                 feature_type: str = None):
        context = {}
        if file_path:
            context["file_path"] = file_path
        if line_number:
            context["line_number"] = line_number
        if feature_type:
            context["feature_type"] = feature_type
        super().__init__(message, context)


class VCFParseError(FileError):
    """Raised for VCF file parsing errors."""

    def __init__(self, message: str, file_path: str = None, line_number: int = None,
                 sample: str = None):
        context = {}
        if file_path:
            context["file_path"] = file_path
        if line_number:
            context["line_number"] = line_number
        if sample:
            context["sample"] = sample
        super().__init__(message, context)
