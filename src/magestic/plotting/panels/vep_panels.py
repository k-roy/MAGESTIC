"""
VEP (Variant Effect Prediction) panel drawing functions.

This module provides panel drawing functions for:
- Evo2 variant effect scores (colored by codon_variant_type)
- ESM1v missense effect predictions (colored by codon_variant_type)
- Shorkie per-gene expression effects (colored by gene_type)
- Yorzoi per-gene expression effects (colored by gene_type)

Styling conventions:
- Genomic/protein models (Evo2, ESM1v): Color by codon_variant_type
- RNA expression models (Yorzoi, Shorkie): Color by gene_type
- Only complete variants are shown (sub-variants filtered out)

Author: Kevin R. Roy
Date: 2026-02-10
"""

from typing import Dict, Optional, List
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

# Import from package config
from ..config import (
    CODON_VARIANT_COLORS,
    GENE_TYPE_COLORS,
    GENE_TYPE_MARKERS,
)


def extract_position_from_variants(variants_str: str) -> Optional[int]:
    """Extract genomic position from individual_variants string."""
    if pd.isna(variants_str) or not variants_str:
        return None
    try:
        first_variant = str(variants_str).split(',')[0].strip()
        parts = first_variant.split('_')
        if parts and parts[0].isdigit():
            return int(parts[0])
    except (ValueError, IndexError):
        pass
    return None


def _filter_complete_variants(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter to only complete variants (exclude sub-variants).

    Complete variants have is_sub_variant=False or missing.
    Also excludes variants with is_filtered=True.
    """
    if df.empty:
        return df

    df = df.copy()

    # Filter out sub-variants
    if 'is_sub_variant' in df.columns:
        df = df[~df['is_sub_variant'].fillna(False)]

    # Filter out filtered variants
    if 'is_filtered' in df.columns:
        df = df[~df['is_filtered'].fillna(False)]

    return df


def _add_codon_type_legend(ax, variant_types_shown: List[str]):
    """Add legend for codon variant types."""
    handles = []
    for vtype in sorted(set(variant_types_shown)):
        if vtype in CODON_VARIANT_COLORS:
            color = CODON_VARIANT_COLORS[vtype]
            label = vtype.replace('_', ' ').title()
            handles.append(Line2D([0], [0], marker='o', color='w',
                                  markerfacecolor=color, markersize=6, label=label))

    if handles:
        ax.legend(handles=handles, loc='upper right', fontsize=5,
                  framealpha=0.8, ncol=min(len(handles), 3),
                  columnspacing=0.5, handletextpad=0.3)


def _add_gene_type_legend(ax, gene_labels: Dict[str, str]):
    """Add legend for gene types (for per-gene expression panels)."""
    handles = []
    for gene_type in ['target', 'upstream_1', 'upstream_2', 'downstream_1', 'downstream_2']:
        if gene_type in gene_labels:
            color = GENE_TYPE_COLORS.get(gene_type, '#95a5a6')
            marker = GENE_TYPE_MARKERS.get(gene_type, 'o')
            label = gene_labels[gene_type]
            handles.append(Line2D([0], [0], marker=marker, color='w',
                                  markerfacecolor=color, markersize=6, label=label))

    if handles:
        ax.legend(handles=handles, loc='upper right', fontsize=5,
                  framealpha=0.8, ncol=min(len(handles), 3),
                  columnspacing=0.5, handletextpad=0.3)


def draw_evo2_panel(
    ax,
    df_evo2: pd.DataFrame,
    region_start: int,
    region_end: int,
    show_legend: bool = True,
):
    """
    Draw Evo2 variant effect scores panel.

    Colors variants by codon_variant_type. Only shows complete variants
    (sub-variants and filtered variants are excluded).

    Args:
        ax: Matplotlib axes to draw on
        df_evo2: DataFrame with evo2 scores
        region_start: Start of genomic region
        region_end: End of genomic region
        show_legend: Whether to show variant type legend
    """
    ax.set_xlim(region_start, region_end)

    if df_evo2.empty:
        ax.text(0.5, 0.5, 'No Evo2 data', ha='center', va='center',
                transform=ax.transAxes, fontsize=8, color='gray')
        ax.set_ylabel('Evo2', fontsize=7)
        return

    # Filter to complete variants only
    df_evo2 = _filter_complete_variants(df_evo2)

    if df_evo2.empty:
        ax.text(0.5, 0.5, 'No complete variants', ha='center', va='center',
                transform=ax.transAxes, fontsize=8, color='gray')
        ax.set_ylabel('Evo2', fontsize=7)
        return

    # Get score column
    score_col = None
    for col in ['variant_score', 'raw_likelihood_score', 'evo2_score']:
        if col in df_evo2.columns:
            score_col = col
            break

    if score_col is None:
        ax.text(0.5, 0.5, 'No score column', ha='center', va='center',
                transform=ax.transAxes, fontsize=8, color='gray')
        return

    # Get position column
    pos_col = 'variant_coord' if 'variant_coord' in df_evo2.columns else 'POS'

    # Track variant types for legend
    variant_types_shown = []

    for idx, row in df_evo2.iterrows():
        x = row.get(pos_col)
        y = row.get(score_col)

        if pd.isna(x) or pd.isna(y):
            continue

        # Color by codon variant type
        vtype = str(row.get('codon_variant_type', 'unknown')).lower()
        color = CODON_VARIANT_COLORS.get(vtype, CODON_VARIANT_COLORS.get('unknown', '#95a5a6'))
        variant_types_shown.append(vtype)

        ax.scatter(x, y, c=color, s=25, alpha=0.8,
                   edgecolors='black', linewidth=0.3)

    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_ylabel('Evo2\nScore', fontsize=7)
    ax.grid(True, alpha=0.3)

    # Add legend
    if show_legend and variant_types_shown:
        _add_codon_type_legend(ax, variant_types_shown)


def draw_esm1v_panel(
    ax,
    df_esm1v: pd.DataFrame,
    region_start: int,
    region_end: int,
    show_legend: bool = True,
):
    """
    Draw ESM1v missense effect predictions panel.

    Colors variants by codon_variant_type. Only shows complete variants
    (sub-variants and filtered variants are excluded).

    Args:
        ax: Matplotlib axes to draw on
        df_esm1v: DataFrame with ESM1v scores
        region_start: Start of genomic region
        region_end: End of genomic region
        show_legend: Whether to show variant type legend
    """
    ax.set_xlim(region_start, region_end)

    if df_esm1v.empty:
        ax.text(0.5, 0.5, 'No ESM1v data', ha='center', va='center',
                transform=ax.transAxes, fontsize=8, color='gray')
        ax.set_ylabel('ESM1v', fontsize=7)
        return

    # Filter to complete variants only
    df_esm1v = _filter_complete_variants(df_esm1v)

    if df_esm1v.empty:
        ax.text(0.5, 0.5, 'No complete variants', ha='center', va='center',
                transform=ax.transAxes, fontsize=8, color='gray')
        ax.set_ylabel('ESM1v', fontsize=7)
        return

    # Find score column
    score_col = 'ESM1v_Score' if 'ESM1v_Score' in df_esm1v.columns else 'esm1v_score'
    if score_col not in df_esm1v.columns:
        ax.text(0.5, 0.5, 'No score column', ha='center', va='center',
                transform=ax.transAxes, fontsize=8, color='gray')
        return

    # Find position column
    df_esm1v = df_esm1v.copy()
    if 'genomic_coords' in df_esm1v.columns:
        df_esm1v['pos'] = pd.to_numeric(df_esm1v['genomic_coords'], errors='coerce')
    elif 'individual_variants' in df_esm1v.columns:
        df_esm1v['pos'] = df_esm1v['individual_variants'].apply(extract_position_from_variants)
    elif 'POS' in df_esm1v.columns:
        df_esm1v['pos'] = df_esm1v['POS']
    elif 'variant_coord' in df_esm1v.columns:
        df_esm1v['pos'] = df_esm1v['variant_coord']
    else:
        ax.text(0.5, 0.5, 'No position column', ha='center', va='center',
                transform=ax.transAxes, fontsize=8, color='gray')
        return

    # Track variant types for legend
    variant_types_shown = []

    # Plot each variant
    for idx, row in df_esm1v.iterrows():
        x = row.get('pos')
        y = row.get(score_col)

        if pd.isna(x) or pd.isna(y):
            continue

        # Color by codon variant type
        vtype = str(row.get('codon_variant_type', 'missense')).lower()
        color = CODON_VARIANT_COLORS.get(vtype, CODON_VARIANT_COLORS.get('missense', '#3498db'))
        variant_types_shown.append(vtype)

        ax.scatter(x, y, c=color, s=25, alpha=0.8,
                   edgecolors='black', linewidth=0.3)

    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_ylabel('ESM1v\nScore', fontsize=7)
    ax.grid(True, alpha=0.3)

    # Add legend
    if show_legend and variant_types_shown:
        _add_codon_type_legend(ax, variant_types_shown)


def draw_shorkie_pergene_panel(
    ax,
    df_pergene: pd.DataFrame,
    region_start: int,
    region_end: int,
    gene_name: str,
    genes_in_region: Optional[Dict[str, str]] = None,
    show_legend: bool = True,
):
    """
    Draw Shorkie per-gene expression effects panel.

    Shows logSED (log2 scaled expression difference) for genes in the zoomed region,
    with each gene represented by a different color and marker shape.

    Colors by gene_type (target=red, upstream_1=blue, etc.).

    Args:
        ax: Matplotlib axes to draw on
        df_pergene: DataFrame with columns: variant_description, gene_systematic,
                    gene_common, gene_type, shorkie_logSED
        region_start: Start of region in plot coordinates
        region_end: End of region in plot coordinates
        gene_name: Target gene name for labeling
        genes_in_region: Dict mapping gene systematic name -> common name for genes to include
        show_legend: Whether to show gene type legend
    """
    ax.set_xlim(region_start, region_end)

    if df_pergene is None or df_pergene.empty:
        ax.text(0.5, 0.5, 'No Shorkie per-gene data', ha='center', va='center',
                transform=ax.transAxes, fontsize=8, color='gray')
        ax.set_ylabel('Shorkie\nlogSED', fontsize=7)
        return

    # Extract variant position from variant_description
    def extract_pos(desc):
        try:
            return int(str(desc).split('_')[0])
        except (ValueError, IndexError):
            return None

    df_pergene = df_pergene.copy()
    df_pergene['pos'] = df_pergene['variant_description'].apply(extract_pos)
    df_pergene = df_pergene.dropna(subset=['pos', 'shorkie_logSED'])

    if df_pergene.empty:
        ax.text(0.5, 0.5, 'No valid Shorkie data', ha='center', va='center',
                transform=ax.transAxes, fontsize=8, color='gray')
        ax.set_ylabel('Shorkie\nlogSED', fontsize=7)
        return

    # Filter to genes in the zoomed region if provided
    if genes_in_region and 'gene_systematic' in df_pergene.columns:
        region_genes = set(genes_in_region.keys())
        mask = df_pergene['gene_systematic'].isin(region_genes)
        df_pergene = df_pergene[mask]

        if df_pergene.empty:
            ax.text(0.5, 0.5, 'No genes in region', ha='center', va='center',
                    transform=ax.transAxes, fontsize=8, color='gray')
            ax.set_ylabel('Shorkie\nlogSED', fontsize=7)
            return

    # Track gene labels for legend
    gene_labels = {}

    # Plot each gene type with its own color and marker
    gene_types_in_data = df_pergene['gene_type'].unique()

    for gene_type in ['target', 'upstream_1', 'upstream_2', 'downstream_1', 'downstream_2']:
        if gene_type not in gene_types_in_data:
            continue

        subset = df_pergene[df_pergene['gene_type'] == gene_type]
        if subset.empty:
            continue

        color = GENE_TYPE_COLORS.get(gene_type, '#95a5a6')
        marker = GENE_TYPE_MARKERS.get(gene_type, 'o')

        # Get gene name for legend
        gene_label = subset['gene_common'].iloc[0] if not subset['gene_common'].isna().all() else gene_type
        gene_labels[gene_type] = gene_label

        ax.scatter(
            subset['pos'],
            subset['shorkie_logSED'],
            c=color,
            marker=marker,
            s=25 if gene_type == 'target' else 18,
            alpha=0.85 if gene_type == 'target' else 0.7,
            edgecolors='black',
            linewidth=0.3,
            zorder=10 if gene_type == 'target' else 5
        )

    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_ylabel('Shorkie\nlogSED', fontsize=7)
    ax.grid(True, alpha=0.3)

    # Add legend
    if show_legend and gene_labels:
        _add_gene_type_legend(ax, gene_labels)


def draw_yorzoi_pergene_panel(
    ax,
    df_pergene: pd.DataFrame,
    region_start: int,
    region_end: int,
    gene_name: str,
    genes_in_region: Optional[Dict[str, str]] = None,
    show_legend: bool = True,
):
    """
    Draw Yorzoi per-gene expression effects panel.

    Shows expression ratio (log2) for genes in the zoomed region,
    with each gene represented by a different color and marker shape.

    Colors by gene_type (target=red, upstream_1=blue, etc.).

    Args:
        ax: Matplotlib axes to draw on
        df_pergene: DataFrame with columns: variant_description, gene_systematic,
                    gene_common, gene_type, yorzoi_ratio
        region_start: Start of region in plot coordinates
        region_end: End of region in plot coordinates
        gene_name: Target gene name for labeling
        genes_in_region: Dict mapping gene systematic name -> common name for genes to include
        show_legend: Whether to show gene type legend
    """
    ax.set_xlim(region_start, region_end)

    if df_pergene is None or df_pergene.empty:
        ax.text(0.5, 0.5, 'No Yorzoi per-gene data', ha='center', va='center',
                transform=ax.transAxes, fontsize=8, color='gray')
        ax.set_ylabel('Yorzoi\nlog₂(ratio)', fontsize=7)
        return

    # Extract variant position from variant_description
    def extract_pos(desc):
        try:
            return int(str(desc).split('_')[0])
        except (ValueError, IndexError):
            return None

    df_pergene = df_pergene.copy()
    df_pergene['pos'] = df_pergene['variant_description'].apply(extract_pos)

    # Calculate log2 ratio
    ratio_col = 'yorzoi_ratio'
    if ratio_col not in df_pergene.columns:
        ax.text(0.5, 0.5, 'No Yorzoi ratio column', ha='center', va='center',
                transform=ax.transAxes, fontsize=8, color='gray')
        ax.set_ylabel('Yorzoi\nlog₂(ratio)', fontsize=7)
        return

    df_pergene['log2_ratio'] = np.log2(df_pergene[ratio_col].clip(lower=0.01))
    df_pergene = df_pergene.dropna(subset=['pos', 'log2_ratio'])

    if df_pergene.empty:
        ax.text(0.5, 0.5, 'No valid Yorzoi data', ha='center', va='center',
                transform=ax.transAxes, fontsize=8, color='gray')
        ax.set_ylabel('Yorzoi\nlog₂(ratio)', fontsize=7)
        return

    # Filter to genes in the zoomed region if provided
    if genes_in_region and 'gene_systematic' in df_pergene.columns:
        region_genes = set(genes_in_region.keys())
        mask = df_pergene['gene_systematic'].isin(region_genes)
        df_pergene = df_pergene[mask]

        if df_pergene.empty:
            ax.text(0.5, 0.5, 'No genes in region', ha='center', va='center',
                    transform=ax.transAxes, fontsize=8, color='gray')
            ax.set_ylabel('Yorzoi\nlog₂(ratio)', fontsize=7)
            return

    # Track gene labels for legend
    gene_labels = {}

    # Plot each gene type with its own color and marker
    gene_types_in_data = df_pergene['gene_type'].unique()

    for gene_type in ['target', 'upstream_1', 'upstream_2', 'downstream_1', 'downstream_2']:
        if gene_type not in gene_types_in_data:
            continue

        subset = df_pergene[df_pergene['gene_type'] == gene_type]
        if subset.empty:
            continue

        color = GENE_TYPE_COLORS.get(gene_type, '#95a5a6')
        marker = GENE_TYPE_MARKERS.get(gene_type, 'o')

        # Get gene name for legend
        gene_label = subset['gene_common'].iloc[0] if not subset['gene_common'].isna().all() else gene_type
        gene_labels[gene_type] = gene_label

        ax.scatter(
            subset['pos'],
            subset['log2_ratio'],
            c=color,
            marker=marker,
            s=25 if gene_type == 'target' else 18,
            alpha=0.85 if gene_type == 'target' else 0.7,
            edgecolors='black',
            linewidth=0.3,
            zorder=10 if gene_type == 'target' else 5
        )

    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_ylabel('Yorzoi\nlog₂(ratio)', fontsize=7)
    ax.grid(True, alpha=0.3)

    # Add legend
    if show_legend and gene_labels:
        _add_gene_type_legend(ax, gene_labels)


def draw_vep_panel_styled(
    ax,
    df_vep: pd.DataFrame,
    score_col: str,
    ylabel: str,
    region_start: int,
    region_end: int,
    use_log2: bool = False,
    show_legend: bool = True,
):
    """
    Draw a generic VEP panel with consistent styling.

    Colors by codon_variant_type. Only shows complete variants
    (sub-variants and filtered variants are excluded).

    Args:
        ax: Matplotlib axes to draw on
        df_vep: DataFrame with VEP scores
        score_col: Name of the column containing scores
        ylabel: Y-axis label
        region_start: Start of genomic region
        region_end: End of genomic region
        use_log2: Whether to apply log2 transform to scores
        show_legend: Whether to show variant type legend
    """
    ax.set_xlim(region_start, region_end)

    if df_vep.empty or score_col not in df_vep.columns:
        ax.text(0.5, 0.5, f'No {ylabel} data', ha='center', va='center',
                transform=ax.transAxes, fontsize=8, color='gray')
        ax.set_ylabel(ylabel, fontsize=7)
        return

    # Filter to complete variants only
    df_vep = _filter_complete_variants(df_vep)

    if df_vep.empty:
        ax.text(0.5, 0.5, 'No complete variants', ha='center', va='center',
                transform=ax.transAxes, fontsize=8, color='gray')
        ax.set_ylabel(ylabel, fontsize=7)
        return

    # Find position column
    pos_col = None
    for col in ['POS', 'variant_coord', 'pos', 'position']:
        if col in df_vep.columns:
            pos_col = col
            break

    if pos_col is None:
        # Try to extract from variant_description
        if 'variant_description' in df_vep.columns:
            df_vep = df_vep.copy()
            df_vep['pos'] = df_vep['variant_description'].apply(extract_position_from_variants)
            pos_col = 'pos'
        else:
            ax.text(0.5, 0.5, 'No position column', ha='center', va='center',
                    transform=ax.transAxes, fontsize=8, color='gray')
            return

    # Track variant types for legend
    variant_types_shown = []

    # Plot each variant
    for idx, row in df_vep.iterrows():
        x = row.get(pos_col)
        y = row.get(score_col)

        if pd.isna(x) or pd.isna(y):
            continue

        # Transform y if needed
        if use_log2 and y > 0:
            y = np.log2(y)
        elif use_log2:
            continue  # Skip non-positive values for log transform

        # Color by codon variant type
        vtype = str(row.get('codon_variant_type', 'unknown')).lower()
        color = CODON_VARIANT_COLORS.get(vtype, CODON_VARIANT_COLORS.get('unknown', '#95a5a6'))
        variant_types_shown.append(vtype)

        ax.scatter(x, y, c=color, s=25, alpha=0.8,
                   edgecolors='black', linewidth=0.3)

    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_ylabel(ylabel, fontsize=7)
    ax.grid(True, alpha=0.3)

    # Add legend
    if show_legend and variant_types_shown:
        _add_codon_type_legend(ax, variant_types_shown)
