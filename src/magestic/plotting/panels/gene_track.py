"""
Gene track panel drawing functions.

This module provides functions for drawing gene structure tracks with
box-arrow representations of CDS features.

Features:
- Pentagon/arrow shapes indicating strand direction
- Vertical staggering to avoid overlaps
- Consistent coloring with per-gene expression panels
- Support for common and systematic gene names

Author: Kevin R. Roy
Date: 2026-02-09
"""

from typing import Dict, List, Optional
from matplotlib.patches import Polygon

# Import from package config
from ..config import GENE_TYPE_COLORS


def assign_feature_levels(features: List[Dict], min_gap: int = 500) -> List[Dict]:
    """
    Assign y-levels to features to avoid overlaps.

    Args:
        features: List of feature dicts with 'start' and 'end'
        min_gap: Minimum gap between features to avoid overlap

    Returns:
        List of feature dicts with added 'level' key
    """
    if not features:
        return features

    # Sort by start position
    sorted_features = sorted(features, key=lambda x: x['start'])

    # Track end position of last feature at each level
    level_ends = []

    for feature in sorted_features:
        # Find first level where this feature fits
        placed = False
        for level_idx, level_end in enumerate(level_ends):
            if feature['start'] > level_end + min_gap:
                # Feature fits at this level
                feature['level'] = level_idx
                level_ends[level_idx] = feature['end']
                placed = True
                break

        if not placed:
            # Need a new level
            feature['level'] = len(level_ends)
            level_ends.append(feature['end'])

    return sorted_features


def draw_gene_track(
    ax,
    gene_name: str,
    gff_features: Dict,
    region_start: int,
    region_end: int,
    highlight_gene: bool = True,
    show_all_features: bool = True,
    show_labels: bool = True,
    use_common_names: bool = True,
    neighbor_genes: Optional[Dict[str, str]] = None,
):
    """
    Draw gene structure track with all CDS annotations in region.

    Features are staggered vertically to avoid overlaps.
    The target gene is highlighted, and neighbor genes (if specified) are
    colored to match the per-gene expression panels.

    Args:
        ax: Matplotlib axes to draw on
        gene_name: Target gene name to highlight
        gff_features: Dict of chromosome -> list of GeneFeature objects
        region_start: Start of genomic region
        region_end: End of genomic region
        highlight_gene: If True, highlight the target gene
        show_all_features: If True, show all CDS in region (not just target)
        show_labels: If False, omit gene name labels (for crowded broad views)
        use_common_names: If True, use common gene names instead of systematic
        neighbor_genes: Optional dict mapping gene_type to gene_name, e.g.,
                       {'upstream_1': 'RSM10', 'downstream_1': 'ENA2'}
                       Used to color genes consistently with expression panels.
    """
    ax.set_xlim(region_start, region_end)

    # Get the gene name lookup table (systematic -> common)
    gene_name_lookup = gff_features.get('_gene_name_lookup', {})

    # Find the chromosome for this gene first
    target_chrom = None
    for chrom, features in gff_features.items():
        if chrom.startswith('_'):  # Skip special keys like _gene_name_lookup
            continue
        for f in features:
            # Check both the feature's common_name and the lookup table
            feature_common = f.common_name or gene_name_lookup.get(f.gene_id, '')
            if feature_common == gene_name or f.gene_id == gene_name:
                target_chrom = chrom
                break
        if target_chrom:
            break

    if target_chrom is None:
        ax.text(0.5, 0.5, f'{gene_name}\n(no GFF)', ha='center', va='center',
                transform=ax.transAxes, fontsize=8, color='gray')
        ax.set_yticks([])
        ax.set_ylim(0, 1)
        return

    # Collect all CDS features in the region
    region_features = []
    for f in gff_features.get(target_chrom, []):
        if f.feature_type == 'CDS':
            # Check if feature overlaps with region
            if f.end >= region_start and f.start <= region_end:
                # Clip to region
                feat_start = max(f.start, region_start)
                feat_end = min(f.end, region_end)

                # Determine names
                systematic_name = f.gene_id or 'unknown'
                # Look up common name from the parent gene feature
                common_name = gene_name_lookup.get(systematic_name, None)
                if not common_name and '_' in systematic_name:
                    # Try without the _CDS suffix
                    parent_name = systematic_name.rsplit('_', 1)[0]
                    common_name = gene_name_lookup.get(parent_name, None)
                # Extract just the systematic name
                if systematic_name and '_' in systematic_name:
                    systematic_name = systematic_name.split('_')[0]
                display_name = common_name if common_name else systematic_name
                is_target = (display_name == gene_name) or (f.gene_id == gene_name) or (common_name == gene_name)

                region_features.append({
                    'start': feat_start,
                    'end': feat_end,
                    'strand': f.strand,
                    'name': display_name,
                    'common_name': common_name,
                    'systematic_name': systematic_name,
                    'is_target': is_target,
                })

    if not region_features:
        ax.text(0.5, 0.5, 'No CDS in region', ha='center', va='center',
                transform=ax.transAxes, fontsize=8, color='gray')
        ax.set_yticks([])
        ax.set_ylim(0, 1)
        return

    # Assign levels to avoid overlaps
    region_features = assign_feature_levels(region_features, min_gap=300)
    max_level = max(f['level'] for f in region_features)

    # Set y limits based on number of levels
    ax.set_ylim(-0.5, max_level + 1.5)

    # Build reverse lookup for neighbor genes (gene_name -> gene_type)
    gene_type_lookup = {}
    if neighbor_genes:
        for gene_type, g_name in neighbor_genes.items():
            if g_name:
                gene_type_lookup[g_name] = gene_type
                gene_type_lookup[g_name.upper()] = gene_type

    # Draw each feature as a pentagon/arrow
    for feature in region_features:
        start = feature['start']
        end = feature['end']
        width = end - start
        y_level = feature['level']

        # Determine gene type for coloring
        feature_common = feature.get('common_name', '')
        feature_systematic = feature.get('systematic_name', '')
        gene_type = None

        if feature['is_target']:
            gene_type = 'target'
        elif neighbor_genes:
            gene_type = gene_type_lookup.get(feature_common) or gene_type_lookup.get(feature_systematic)

        # Colors based on gene type (for consistency with expression panels)
        if gene_type and gene_type in GENE_TYPE_COLORS:
            color = GENE_TYPE_COLORS[gene_type]
            alpha = 0.9 if gene_type == 'target' else 0.8
            linewidth = 1.5 if gene_type == 'target' else 1.0
        elif feature['is_target']:
            color = '#e74c3c'  # Red
            alpha = 0.95
            linewidth = 1.5
        else:
            color = '#95a5a6'  # Gray
            alpha = 0.5
            linewidth = 0.8

        # Arrow tip size
        arrow_tip = min(width * 0.15, 400)

        # Y coordinates for pentagon
        y_bot = y_level - 0.4
        y_mid = y_level
        y_top = y_level + 0.4

        # Create pentagonal arrow shape
        if feature['strand'] == '+':
            vertices = [
                (start, y_bot),
                (start, y_top),
                (end - arrow_tip, y_top),
                (end, y_mid),
                (end - arrow_tip, y_bot),
            ]
        else:
            vertices = [
                (start + arrow_tip, y_top),
                (end, y_top),
                (end, y_bot),
                (start + arrow_tip, y_bot),
                (start, y_mid),
            ]

        arrow_patch = Polygon(vertices, closed=True, facecolor=color,
                              edgecolor='black', linewidth=linewidth, alpha=alpha)
        ax.add_patch(arrow_patch)

        # Label (skip if show_labels is False)
        if show_labels:
            label_x = (start + end) / 2
            label_y = y_level

            # Determine label text based on use_common_names
            if use_common_names:
                label_text = feature.get('common_name') or feature.get('systematic_name', feature['name'])
            else:
                label_text = feature.get('systematic_name', feature['name'])

            # Shorter labels for non-target genes
            if feature['is_target']:
                fontsize = 8
                fontweight = 'bold'
                fontcolor = 'white'
            else:
                label_text = label_text[:10] if len(label_text) > 10 else label_text
                fontsize = 6
                fontweight = 'normal'
                fontcolor = 'white'

            ax.text(label_x, label_y, label_text, ha='center', va='center',
                    fontsize=fontsize, fontweight=fontweight, color=fontcolor, alpha=0.95)

    ax.set_yticks([])
    ax.set_ylabel('Genes', fontsize=8)


def get_genes_in_region(
    gff_features: Dict,
    chrom: str,
    region_start: int,
    region_end: int,
) -> Dict[str, str]:
    """
    Get all genes (CDS features) that overlap a genomic region.

    Args:
        gff_features: Dict of chromosome -> list of GeneFeature objects
        chrom: Chromosome name
        region_start: Start of region
        region_end: End of region

    Returns:
        Dict mapping systematic name -> common name for genes in region
    """
    gene_name_lookup = gff_features.get('_gene_name_lookup', {})
    genes_in_region = {}

    for f in gff_features.get(chrom, []):
        if f.feature_type == 'CDS':
            if f.end >= region_start and f.start <= region_end:
                systematic_name = f.gene_id or 'unknown'
                if '_' in systematic_name:
                    systematic_name = systematic_name.split('_')[0]

                common_name = gene_name_lookup.get(systematic_name, None)
                if not common_name and f.gene_id and '_' in f.gene_id:
                    parent_name = f.gene_id.rsplit('_', 1)[0]
                    common_name = gene_name_lookup.get(parent_name, None)

                genes_in_region[systematic_name] = common_name or systematic_name

    return genes_in_region
