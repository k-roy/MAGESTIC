"""
Count matrix functions for the BC1-From-Screen pipeline.

Creates DESeq2-ready count matrices from long-format BC1 counts.
"""

from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np

from ..config import LibraryConfig


def create_wide_count_matrix(
    long_table: pd.DataFrame,
    bc1_col: str = "bc1",
    sample_col: str = "sample_name",
    count_col: str = "sample_counts",
    fill_value: int = 0,
) -> pd.DataFrame:
    """
    Convert long-format BC1 counts to wide format for DESeq2.

    Args:
        long_table: Long-format BC1 counts
        bc1_col: Name of BC1 column
        sample_col: Name of sample column
        count_col: Name of count column
        fill_value: Value to use for missing counts

    Returns:
        Wide-format DataFrame with BC1s as rows, samples as columns
    """
    wide = long_table.pivot(
        index=bc1_col,
        columns=sample_col,
        values=count_col
    ).fillna(fill_value).astype(int)

    return wide


def consolidate_t0_samples(
    wide_matrix: pd.DataFrame,
    keyfile: pd.DataFrame,
    library_config: LibraryConfig,
    sample_col: str = "sample_name",
    condition_col: str = "condition",
    generation_col: str = "generations",
    screen_date_col: str = "screen_date",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Consolidate t0 samples into replicate groups.

    For example, 24 t0 samples -> 4 consolidated t0 replicates.

    Args:
        wide_matrix: Wide-format count matrix
        keyfile: Sample metadata
        library_config: Library configuration with consolidation settings
        sample_col: Name of sample column
        condition_col: Name of condition column
        generation_col: Name of generation column
        screen_date_col: Name of screen date column

    Returns:
        Tuple of (consolidated matrix, updated metadata)
    """
    # Identify t0 samples for this library
    # Handle generations as numeric for robust comparison
    generations_numeric = pd.to_numeric(keyfile[generation_col], errors="coerce").fillna(-1).astype(int)
    t0_mask = (
        (keyfile[condition_col] == "YPD") &
        (generations_numeric == 0) &
        (keyfile[screen_date_col].isin(library_config.screen_dates))
    )
    t0_samples = keyfile[t0_mask][sample_col].tolist()

    # Non-t0 samples (keep as-is)
    non_t0_samples = [col for col in wide_matrix.columns if col not in t0_samples]

    # Calculate samples per consolidated group
    n_t0 = len(t0_samples)
    n_groups = library_config.t0_consolidation_groups
    samples_per_group = n_t0 // n_groups

    if samples_per_group == 0:
        raise ValueError(
            f"Not enough t0 samples ({n_t0}) for {n_groups} consolidated groups"
        )

    # Warn if there are leftover samples that won't be included
    leftover_samples = n_t0 % n_groups
    if leftover_samples > 0:
        import logging
        logging.getLogger(__name__).warning(
            f"T0 consolidation: {leftover_samples} of {n_t0} samples will not be "
            f"included (using {samples_per_group} samples per group × {n_groups} groups)"
        )

    # Sort t0 samples for consistent grouping
    t0_samples_sorted = sorted(t0_samples)

    # Create consolidated t0 columns
    consolidated_t0_data = {}
    consolidated_metadata = []

    for i in range(n_groups):
        group_samples = t0_samples_sorted[i * samples_per_group:(i + 1) * samples_per_group]
        consolidated_name = f"YPD_0g_rep{i + 1}"

        # Sum counts across samples in group
        consolidated_t0_data[consolidated_name] = wide_matrix[group_samples].sum(axis=1)

        # Create metadata entry
        consolidated_metadata.append({
            sample_col: consolidated_name,
            condition_col: "YPD",
            generation_col: 0,
            "replicate": i + 1,
            "n_samples_consolidated": len(group_samples),
            "source_samples": ",".join(group_samples),
        })

    # Build consolidated matrix
    consolidated_matrix = pd.concat([
        wide_matrix[non_t0_samples],
        pd.DataFrame(consolidated_t0_data, index=wide_matrix.index)
    ], axis=1)

    # Build metadata DataFrame
    # Include non-t0 samples from original keyfile
    non_t0_metadata = keyfile[~t0_mask].to_dict('records')
    all_metadata = non_t0_metadata + consolidated_metadata
    metadata_df = pd.DataFrame(all_metadata)

    return consolidated_matrix, metadata_df


def split_by_comparison(
    wide_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    sample_col: str = "sample_name",
    condition_col: str = "condition",
    generation_col: str = "generations",
    reference_condition: str = "YPD",
    reference_generation: int = 0,
) -> Dict[str, Dict[str, pd.DataFrame]]:
    """
    Split count matrix by comparison (reference vs treatment).

    Args:
        wide_matrix: Wide-format count matrix
        metadata: Sample metadata
        sample_col: Name of sample column
        condition_col: Name of condition column
        generation_col: Name of generation column
        reference_condition: Reference condition name
        reference_generation: Reference generation (typically 0)

    Returns:
        Dictionary mapping comparison names to {reference, treatment} DataFrames
    """
    # Normalize generations column to numeric for robust comparison
    metadata = metadata.copy()
    metadata["_generations_numeric"] = pd.to_numeric(
        metadata[generation_col], errors="coerce"
    ).fillna(-1).astype(int)

    # Identify reference samples
    ref_mask = (
        (metadata[condition_col] == reference_condition) &
        (metadata["_generations_numeric"] == reference_generation)
    )
    ref_samples = metadata[ref_mask][sample_col].tolist()

    # Get unique treatment conditions (non-reference)
    treatment_conditions = metadata[~ref_mask].groupby(
        [condition_col, "_generations_numeric"]
    ).size().reset_index()[
        [condition_col, "_generations_numeric"]
    ].values.tolist()

    comparisons = {}

    for condition, generation in treatment_conditions:
        comparison_name = f"{condition}_{int(generation)}g"

        # Get treatment samples for this comparison
        treat_mask = (
            (metadata[condition_col] == condition) &
            (metadata["_generations_numeric"] == generation)
        )
        treat_samples = metadata[treat_mask][sample_col].tolist()

        comparisons[comparison_name] = {
            "reference": wide_matrix[ref_samples],
            "treatment": wide_matrix[treat_samples],
            "reference_samples": ref_samples,
            "treatment_samples": treat_samples,
        }

    return comparisons


def create_comparison_manifest(
    comparisons: Dict[str, Dict],
    library_name: str,
    output_path: Optional[str] = None,
) -> pd.DataFrame:
    """
    Create a manifest file listing all comparisons.

    Args:
        comparisons: Dictionary from split_by_comparison()
        library_name: Name of the library
        output_path: Optional path to save manifest

    Returns:
        DataFrame with comparison metadata
    """
    manifest_rows = []

    for comparison_name, data in comparisons.items():
        manifest_rows.append({
            "comparison": comparison_name,
            "library": library_name,
            "n_reference_samples": len(data["reference_samples"]),
            "n_treatment_samples": len(data["treatment_samples"]),
            "reference_samples": ",".join(data["reference_samples"]),
            "treatment_samples": ",".join(data["treatment_samples"]),
        })

    manifest = pd.DataFrame(manifest_rows)

    if output_path:
        manifest.to_csv(output_path, sep="\t", index=False)

    return manifest
