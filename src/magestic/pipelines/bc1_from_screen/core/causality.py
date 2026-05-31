"""
Causality evidence framework for variant hit calling.

This module provides a structured framework for evaluating causality evidence
and making binary causality calls based on multiple lines of evidence:
- Statistical significance (empirical FDR-calibrated)
- Replication (BC1 count and concordance)
- Cross-validation (tier, time series, nuclease)
- Noise assessment (penalties for noisy conditions/loci)

Author: Kevin R. Roy
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np

from .sublibrary import NUCLEASE_PAM_PREFIXES, NUCLEASE_PAM_PATTERN


# ============================================================================
# Evidence Flags
# ============================================================================

class EvidenceFlag(Enum):
    """Evidence flag codes for causality assessment."""
    TS = "Time Series validated"
    CT = "Cross-Tier validated"
    KO = "Knockout variant"
    ML = "Multi-Library concordant"
    MN = "Multi-Nuclease concordant"
    HR = "High replication"  # ≥5 BC1s
    NC = "Noisy condition"
    NL = "Noisy locus"
    SV = "SV-prone region"


# ============================================================================
# Stringency Levels
# ============================================================================

@dataclass
class StringencyThresholds:
    """Thresholds for a given stringency level."""
    fdr: float
    syn_percentile: float
    min_bc1: int
    min_concordance: float


STRINGENCY_LEVELS = {
    "relaxed": StringencyThresholds(
        fdr=0.10,
        syn_percentile=95.0,
        min_bc1=2,
        min_concordance=0.67,
    ),
    "moderate": StringencyThresholds(
        fdr=0.05,
        syn_percentile=97.5,
        min_bc1=3,
        min_concordance=0.75,
    ),
    "stringent": StringencyThresholds(
        fdr=0.01,
        syn_percentile=99.0,
        min_bc1=5,
        min_concordance=0.85,
    ),
}


# ============================================================================
# Causality Evidence Dataclass
# ============================================================================

@dataclass
class CausalityEvidence:
    """
    Complete evidence record for a variant's causality assessment.

    This captures all lines of evidence used for making causality calls.
    """
    # Identifiers
    variant_id: str
    gene: str
    aa_change: Optional[str]
    condition: str
    library: str

    # Core statistics (from empirical calibration)
    lfc_centered: float
    z_score_empirical: float
    syn_percentile: float
    empirical_z_threshold: float

    # Replication evidence
    n_bc1s: int
    bc1_concordance: float
    n_oligos: int = 1

    # Cross-validation evidence
    cross_tier_validated: bool = False
    time_series_validated: bool = False
    time_series_r2: float = 0.0
    multi_nuclease_concordant: bool = False
    n_nuclease_libraries: int = 1

    # Noise flags
    noisy_condition: bool = False
    high_variance_locus: bool = False
    is_high_sv: bool = False
    noise_penalty: str = "none"
    confidence: str = "high"

    # Variant type
    codon_variant_type: str = ""
    is_ko_variant: bool = False  # stop_gain, frameshift, deletion

    # Evidence flags (computed)
    evidence_flags: List[str] = field(default_factory=list)

    def compute_evidence_flags(self) -> List[str]:
        """Compute list of evidence flags based on attributes."""
        flags = []

        if self.time_series_validated:
            flags.append("TS")
        if self.cross_tier_validated:
            flags.append("CT")
        if self.is_ko_variant:
            flags.append("KO")
        if self.multi_nuclease_concordant:
            flags.append("MN")
        if self.n_bc1s >= 5:
            flags.append("HR")
        if self.noisy_condition:
            flags.append("NC")
        if self.high_variance_locus:
            flags.append("NL")
        if self.is_high_sv:
            flags.append("SV")

        self.evidence_flags = flags
        return flags

    def is_causal(self, stringency: str = "moderate") -> bool:
        """
        Make binary causality call based on evidence.

        Args:
            stringency: 'relaxed', 'moderate', or 'stringent'

        Returns:
            True if variant passes causality threshold
        """
        thresholds = STRINGENCY_LEVELS.get(stringency, STRINGENCY_LEVELS["moderate"])

        # Must pass significance threshold
        passes_significance = abs(self.z_score_empirical) >= self.empirical_z_threshold

        # Must pass replication threshold
        passes_replication = (
            self.n_bc1s >= thresholds.min_bc1 and
            self.bc1_concordance >= thresholds.min_concordance
        )

        # Must pass percentile threshold
        passes_percentile = self.syn_percentile >= thresholds.syn_percentile

        return passes_significance and passes_replication and passes_percentile

    def to_dict(self) -> Dict:
        """Convert to dictionary for DataFrame construction."""
        self.compute_evidence_flags()
        return {
            "variant_id": self.variant_id,
            "gene": self.gene,
            "aa_change": self.aa_change,
            "condition": self.condition,
            "library": self.library,
            "lfc_centered": self.lfc_centered,
            "z_score_empirical": self.z_score_empirical,
            "syn_percentile": self.syn_percentile,
            "empirical_z_threshold": self.empirical_z_threshold,
            "n_bc1s": self.n_bc1s,
            "bc1_concordance": self.bc1_concordance,
            "n_oligos": self.n_oligos,
            "cross_tier_validated": self.cross_tier_validated,
            "time_series_validated": self.time_series_validated,
            "time_series_r2": self.time_series_r2,
            "multi_nuclease_concordant": self.multi_nuclease_concordant,
            "n_nuclease_libraries": self.n_nuclease_libraries,
            "noisy_condition": self.noisy_condition,
            "high_variance_locus": self.high_variance_locus,
            "is_high_sv": self.is_high_sv,
            "noise_penalty": self.noise_penalty,
            "confidence": self.confidence,
            "codon_variant_type": self.codon_variant_type,
            "is_ko_variant": self.is_ko_variant,
            "evidence_flags": ",".join(self.evidence_flags),
            "is_causal_relaxed": self.is_causal("relaxed"),
            "is_causal_moderate": self.is_causal("moderate"),
            "is_causal_stringent": self.is_causal("stringent"),
        }


# ============================================================================
# Causality Evaluator
# ============================================================================

class CausalityEvaluator:
    """
    Evaluator class for making causality calls on variant DataFrames.
    """

    def __init__(
        self,
        stringency: str = "moderate",
        ko_types: List[str] = None,
    ):
        """
        Initialize evaluator.

        Args:
            stringency: Default stringency level
            ko_types: Variant types considered knockout
        """
        self.stringency = stringency
        self.thresholds = STRINGENCY_LEVELS.get(stringency, STRINGENCY_LEVELS["moderate"])

        if ko_types is None:
            self.ko_types = ["stop_gain", "frameshift", "nonsense"]
        else:
            self.ko_types = ko_types

    def evaluate_variant(
        self,
        row: pd.Series,
        variant_id_col: str = "sublibrary_variant_combo",
        gene_col: str = "common_gene_name_locus",
        aa_change_col: str = "aa_change",
        condition_col: str = "comparison",
        library_col: str = "library",
        lfc_col: str = "lfc_centered",
        z_score_col: str = "z_score_empirical",
        syn_percentile_col: str = "syn_percentile",
        z_thresh_col: str = "empirical_z_threshold",
        n_bc1_col: str = "n_bc1",
        concordance_col: str = "concordance_fraction",
        variant_type_col: str = "codon_variant_type",
    ) -> CausalityEvidence:
        """
        Create CausalityEvidence from a DataFrame row.

        Args:
            row: Single row from variant DataFrame
            *_col: Column name mappings

        Returns:
            CausalityEvidence object
        """
        # Determine if KO variant
        variant_type = str(row.get(variant_type_col, "")).lower()
        is_ko = any(ko_type in variant_type for ko_type in self.ko_types)

        return CausalityEvidence(
            variant_id=str(row.get(variant_id_col, "")),
            gene=str(row.get(gene_col, "")),
            aa_change=row.get(aa_change_col),
            condition=str(row.get(condition_col, "")),
            library=str(row.get(library_col, "")),
            lfc_centered=float(row.get(lfc_col, 0)),
            z_score_empirical=float(row.get(z_score_col, 0)),
            syn_percentile=float(row.get(syn_percentile_col, 0)),
            empirical_z_threshold=float(row.get(z_thresh_col, 3.0)),
            n_bc1s=int(row.get(n_bc1_col, 0)),
            bc1_concordance=float(row.get(concordance_col, 0)),
            cross_tier_validated=bool(row.get("cross_tier_validated", False)),
            time_series_validated=bool(row.get("time_series_validated", False)),
            time_series_r2=float(row.get("time_series_r2", 0)),
            multi_nuclease_concordant=bool(row.get("multi_nuclease_concordant", False)),
            n_nuclease_libraries=int(row.get("n_nuclease_libraries", 1)),
            noisy_condition=bool(row.get("noisy_condition", False)),
            high_variance_locus=bool(row.get("high_variance_locus", False)),
            is_high_sv=bool(row.get("is_high_sv", False)),
            noise_penalty=str(row.get("noise_penalty", "none")),
            confidence=str(row.get("confidence", "high")),
            codon_variant_type=variant_type,
            is_ko_variant=is_ko,
        )

    def evaluate_dataframe(
        self,
        df: pd.DataFrame,
        **col_mappings,
    ) -> pd.DataFrame:
        """
        Evaluate causality for all variants in DataFrame.

        Args:
            df: DataFrame with variant-level results
            **col_mappings: Column name mappings

        Returns:
            DataFrame with causality calls
        """
        # Use itertuples for ~5-10x speedup over iterrows
        evidence_records = []

        for row in df.itertuples(index=False):
            # Convert namedtuple to dict-like access for evaluate_variant
            row_dict = row._asdict()
            evidence = self.evaluate_variant(row_dict, **col_mappings)
            evidence_records.append(evidence.to_dict())

        result_df = pd.DataFrame(evidence_records)

        # Set the final is_causal based on default stringency
        result_df["is_causal"] = result_df[f"is_causal_{self.stringency}"]

        return result_df


# ============================================================================
# Aggregation Functions
# ============================================================================

def aggregate_to_aa_level(
    df: pd.DataFrame,
    condition_col: str = "condition",
    aa_change_col: str = "aa_change",
    gene_col: str = "gene",
    is_causal_col: str = "is_causal",
    lfc_col: str = "lfc_centered",
    library_col: str = "library",
    min_concordance: float = 0.5,
) -> pd.DataFrame:
    """
    Aggregate variant-level calls to amino acid change level.

    Multiple DNA variants can produce the same AA change. This function
    combines evidence across DNA variants.

    Args:
        df: DataFrame with variant-level causality calls
        condition_col: Column with condition name
        aa_change_col: Column with AA change
        gene_col: Column with gene name
        is_causal_col: Column with causality call
        lfc_col: Column with LFC
        library_col: Column with library name
        min_concordance: Minimum fraction of variants agreeing

    Returns:
        DataFrame with AA-level causality calls
    """
    # Filter to variants with AA changes
    df_aa = df[df[aa_change_col].notna()].copy()

    if len(df_aa) == 0:
        return pd.DataFrame()

    # Group by (condition, gene, aa_change)
    agg_records = []

    for (condition, gene, aa_change), group in df_aa.groupby([condition_col, gene_col, aa_change_col]):
        n_variants = len(group)
        n_causal = group[is_causal_col].sum()
        concordance = n_causal / n_variants if n_variants > 0 else 0

        # LFC statistics
        mean_lfc = group[lfc_col].mean()
        std_lfc = group[lfc_col].std()

        # Direction check
        n_positive = (group[lfc_col] > 0).sum()
        n_negative = (group[lfc_col] < 0).sum()
        direction_concordance = max(n_positive, n_negative) / n_variants

        # Nuclease library check
        n_libraries = group[library_col].nunique() if library_col in group.columns else 1

        # Evidence flags
        evidence_flags = []
        if "evidence_flags" in group.columns:
            all_flags = ",".join(group["evidence_flags"].dropna()).split(",")
            evidence_flags = list(set(f for f in all_flags if f))

        # Confidence (take most conservative)
        if "confidence" in group.columns:
            confidence_order = {"high": 0, "medium": 1, "low": 2}
            confidences = group["confidence"].dropna()
            if len(confidences) > 0:
                confidence = max(confidences, key=lambda x: confidence_order.get(x, 1))
            else:
                confidence = "medium"
        else:
            confidence = "medium"

        # Make AA-level causality call
        is_causal_aa = concordance >= min_concordance and n_causal > 0

        agg_records.append({
            "condition": condition,
            "gene": gene,
            "aa_change": aa_change,
            "n_dna_variants": n_variants,
            "n_causal_variants": int(n_causal),
            "variant_concordance": concordance,
            "mean_lfc": mean_lfc,
            "std_lfc": std_lfc,
            "direction_concordance": direction_concordance,
            "n_nuclease_libraries": n_libraries,
            "evidence_flags": ",".join(evidence_flags),
            "confidence": confidence,
            "is_causal_aa": is_causal_aa,
        })

    return pd.DataFrame(agg_records)


def check_ngnh_exception(
    df: pd.DataFrame,
    variant_id_col: str = "variant_id",
    library_col: str = "library",
    is_causal_col: str = "is_causal",
) -> pd.DataFrame:
    """
    Check for SpG_NGNH exception cases.

    SpG_NGNH has lower editing efficiency. If a variant is causal in
    SpCas9_NGG, SpG_NGG, SpG_NGNG but NOT SpG_NGNH, we still consider
    it causal with a flag.

    Args:
        df: DataFrame with variant-level calls
        variant_id_col: Column with variant identifier
        library_col: Column with library/nuclease info
        is_causal_col: Column with causality call

    Returns:
        DataFrame with ngnh_exception column added
    """
    df = df.copy()
    df["ngnh_exception"] = False

    # Extract nuclease using vectorized str.extract (50x faster than apply)
    # Try library_col first, then variant_id_col as fallback
    # Uses NUCLEASE_PAM_PATTERN from sublibrary.py for consistency

    df["nuclease"] = df[library_col].astype(str).str.extract(NUCLEASE_PAM_PATTERN, expand=False)

    # Fill missing with variant_id_col extraction
    mask_missing = df["nuclease"].isna()
    if mask_missing.any():
        df.loc[mask_missing, "nuclease"] = (
            df.loc[mask_missing, variant_id_col]
            .astype(str)
            .str.extract(NUCLEASE_PAM_PATTERN, expand=False)
        )

    # Fill remaining with "unknown"
    df["nuclease"] = df["nuclease"].fillna("unknown")

    # Group by variant identity (without nuclease prefix)
    # Build pattern dynamically from NUCLEASE_PAM_PREFIXES
    prefix_pattern = r"^(" + "|".join(NUCLEASE_PAM_PREFIXES) + r")_"
    df["variant_identity"] = df[variant_id_col].str.replace(prefix_pattern, "", regex=True)

    # Check each variant identity
    for variant_identity, group in df.groupby("variant_identity"):
        nucleases = set(group["nuclease"])

        # Check if NGNH is the only one failing
        if "SpG_NGNH" in nucleases and len(nucleases) > 1:
            ngnh_mask = group["nuclease"] == "SpG_NGNH"
            others_mask = group["nuclease"] != "SpG_NGNH"

            ngnh_causal = group.loc[ngnh_mask, is_causal_col].any()
            others_causal = group.loc[others_mask, is_causal_col].any()

            if others_causal and not ngnh_causal:
                # Mark NGNH as exception
                df.loc[group.index[ngnh_mask], "ngnh_exception"] = True

    return df


# ============================================================================
# Output Generation
# ============================================================================

def generate_causal_variants_output(
    df: pd.DataFrame,
    output_stringency: str = "moderate",
) -> pd.DataFrame:
    """
    Generate final causal variants output file.

    Filters to causal variants and formats output columns.

    Args:
        df: DataFrame with causality evaluations
        output_stringency: Stringency level for output

    Returns:
        DataFrame ready for output
    """
    is_causal_col = f"is_causal_{output_stringency}"

    # Filter to causal
    causal_df = df[df[is_causal_col]].copy()

    # Select output columns
    output_cols = [
        "variant_id",
        "gene",
        "aa_change",
        "condition",
        "library",
        "is_causal",
        "confidence",
        "lfc_centered",
        "z_score_empirical",
        "syn_percentile",
        "n_bc1s",
        "bc1_concordance",
        "evidence_flags",
        "codon_variant_type",
        "is_ko_variant",
        "cross_tier_validated",
        "time_series_validated",
        "noisy_condition",
        "high_variance_locus",
        "is_high_sv",
    ]

    # Keep only columns that exist
    output_cols = [c for c in output_cols if c in causal_df.columns]

    return causal_df[output_cols].sort_values(
        ["gene", "condition", "syn_percentile"],
        ascending=[True, True, False],
    )


def generate_control_qc_output(
    calibration_report: pd.DataFrame,
    validation_stats: Dict,
    noise_qc_report: pd.DataFrame,
) -> pd.DataFrame:
    """
    Generate control QC output file.

    Combines calibration, validation, and noise QC metrics.

    Args:
        calibration_report: From empirical_thresholds
        validation_stats: From empirical_thresholds
        noise_qc_report: From noise_assessment

    Returns:
        DataFrame with QC metrics per condition
    """
    # Start with calibration report
    qc_df = calibration_report.copy()

    # Add KO enrichment (global)
    if "ko_enrichment" in validation_stats:
        ko_stats = validation_stats["ko_enrichment"]
        qc_df["ko_enrichment"] = ko_stats.get("ko_enrichment", 0)
        qc_df["ko_enrichment_status"] = ko_stats.get("status", "")

    # Merge noise QC if available
    if len(noise_qc_report) > 0:
        noise_cols = [
            "library", "condition",
            "is_noisy_condition", "n_high_variance_loci",
            "n_high_sv_variants", "pct_with_noise_flag",
            "n_hits_penalized", "pct_hits_penalized",
        ]
        noise_subset = noise_qc_report[[c for c in noise_cols if c in noise_qc_report.columns]]

        qc_df = qc_df.merge(
            noise_subset,
            on=["library", "condition"],
            how="left",
        )

    # Determine overall QC status using vectorized operations
    # Build issue flags as boolean columns
    has_high_fdr = qc_df.get("calibration_status", pd.Series(index=qc_df.index)) == "WARN_HIGH_FDR"
    has_low_n = qc_df.get("calibration_status", pd.Series(index=qc_df.index)) == "WARN_LOW_N"
    has_noisy = qc_df.get("is_noisy_condition", pd.Series(False, index=qc_df.index)).fillna(False)
    has_no_ko = qc_df.get("ko_enrichment", pd.Series(2.0, index=qc_df.index)).fillna(2.0) < 1.0

    # Count total issues per row
    issue_count = has_high_fdr.astype(int) + has_low_n.astype(int) + has_noisy.astype(int) + has_no_ko.astype(int)

    # Build status strings
    qc_df["qc_status"] = "PASS"

    # Single issue cases
    single_issue_mask = issue_count == 1
    qc_df.loc[single_issue_mask & has_high_fdr, "qc_status"] = "WARN:HIGH_FDR"
    qc_df.loc[single_issue_mask & has_low_n, "qc_status"] = "WARN:LOW_N"
    qc_df.loc[single_issue_mask & has_noisy, "qc_status"] = "WARN:NOISY_COND"
    qc_df.loc[single_issue_mask & has_no_ko, "qc_status"] = "WARN:NO_KO_ENRICH"

    # Multiple issues - use itertuples for complex string building (small DataFrame)
    multi_issue_mask = issue_count > 1
    if multi_issue_mask.any():
        for row in qc_df[multi_issue_mask].itertuples():
            issues = []
            if has_high_fdr.loc[row.Index]:
                issues.append("HIGH_FDR")
            if has_low_n.loc[row.Index]:
                issues.append("LOW_N")
            if has_noisy.loc[row.Index]:
                issues.append("NOISY_COND")
            if has_no_ko.loc[row.Index]:
                issues.append("NO_KO_ENRICH")
            qc_df.loc[row.Index, "qc_status"] = "FAIL:" + ",".join(issues)

    return qc_df
