"""
Annotation module for BC1-From-Screen pipeline.

Uses pre-parsed harmonized QTL oligo table for variant annotations,
eliminating fragile runtime oligo name parsing.

The harmonized table contains all variant information already parsed:
- 70,656 oligos from 2021 and 2024 Twist pools
- Complete variant annotations (gene, AA change, coordinates, etc.)
- Additional metadata (SCORE-SV, linkage relationships, etc.)

For backward compatibility, the VariantAnnotation dataclass and function
signatures are preserved, but all parsing now uses table lookups.
"""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

# Import from the unified oligo_annotations module
from magestic.utils.oligo_annotations import (
    load_oligo_annotations,
    annotate_dataframe as _annotate_dataframe,
    get_annotation_columns_for_bc1_from_screen,
    COLUMN_MAPPING_BC1_FROM_SCREEN,
)

logger = logging.getLogger(__name__)

# Cache for annotation lookups
_ANNOTATION_CACHE: Optional[Dict[str, dict]] = None


@dataclass
class VariantAnnotation:
    """Parsed variant annotation from oligo name.

    Kept for backward compatibility. Fields are now populated from
    the pre-parsed harmonized table rather than runtime parsing.
    """
    oligo_name: str
    cas_type: str = ""
    gene_systematic: str = ""
    gene_common: str = ""
    codon_position: int = 0
    original_aa: str = ""
    new_aa: str = ""
    original_codon: str = ""
    new_codon: str = ""
    upstream_changes: int = 0
    downstream_changes: int = 0
    strand: str = ""
    chromosome: str = ""
    guide_strand: str = ""
    guide_position: int = 0
    pam_sequence: str = ""
    edit_start: int = 0
    edit_end: int = 0
    variant_type: str = ""  # missense, nonsense, synonymous, control
    is_control: bool = False

    def to_dict(self) -> Dict:
        """Convert to dictionary for DataFrame creation."""
        return {
            "oligo_name": self.oligo_name,
            "cas_type": self.cas_type,
            "gene_systematic": self.gene_systematic,
            "gene_common": self.gene_common,
            "codon_position": self.codon_position,
            "original_aa": self.original_aa,
            "new_aa": self.new_aa,
            "original_codon": self.original_codon,
            "new_codon": self.new_codon,
            "upstream_changes": self.upstream_changes,
            "downstream_changes": self.downstream_changes,
            "strand": self.strand,
            "chromosome": self.chromosome,
            "guide_strand": self.guide_strand,
            "guide_position": self.guide_position,
            "pam_sequence": self.pam_sequence,
            "edit_start": self.edit_start,
            "edit_end": self.edit_end,
            "variant_type": self.variant_type,
            "is_control": self.is_control,
        }


def _get_annotation_lookup() -> Dict[str, dict]:
    """
    Load and cache the annotation lookup dictionary.

    Returns:
        Dictionary mapping oligo_name -> annotation dict
    """
    global _ANNOTATION_CACHE

    if _ANNOTATION_CACHE is not None:
        return _ANNOTATION_CACHE

    logger.info("Loading annotation lookup from harmonized table...")

    # Load with legacy column names for backward compatibility
    df = load_oligo_annotations(
        columns=get_annotation_columns_for_bc1_from_screen(include_extras=False),
        rename_to_legacy=True,
        add_is_control=True,
    )

    # Build lookup dictionary using itertuples for ~5x speedup over iterrows
    _ANNOTATION_CACHE = {}
    for row in df.itertuples(index=False):
        oligo_name = getattr(row, "oligo_name", "")
        if oligo_name and pd.notna(oligo_name):
            _ANNOTATION_CACHE[oligo_name] = row._asdict()

    logger.info(f"  Loaded {len(_ANNOTATION_CACHE):,} oligo annotations")

    return _ANNOTATION_CACHE


def classify_variant_type(original_aa: str, new_aa: str) -> str:
    """
    Classify the variant type based on amino acid change.

    Args:
        original_aa: Original amino acid (single letter)
        new_aa: New amino acid (single letter or *)

    Returns:
        Variant type: 'nonsense', 'missense', 'synonymous', or 'unknown'

    Note:
        This function is kept for backward compatibility. The harmonized
        table already contains the variant_type column (codon_variant_type).
    """
    if not original_aa or not new_aa:
        return "unknown"

    if new_aa == "*":
        return "nonsense"
    elif original_aa == new_aa:
        return "synonymous"
    else:
        return "missense"


def parse_oligo_name(oligo_name: str) -> VariantAnnotation:
    """
    Get variant annotation for an oligo by lookup.

    This function now uses the pre-parsed harmonized table instead of
    runtime parsing. Kept for backward compatibility.

    Args:
        oligo_name: Full oligo name string

    Returns:
        VariantAnnotation with fields from harmonized table
    """
    annot = VariantAnnotation(oligo_name=oligo_name)

    if not oligo_name or pd.isna(oligo_name):
        return annot

    # Look up in harmonized table
    lookup = _get_annotation_lookup()

    if oligo_name in lookup:
        row = lookup[oligo_name]

        # Populate from pre-parsed data
        annot.cas_type = str(row.get("cas_type", "")) if pd.notna(row.get("cas_type")) else ""
        annot.gene_systematic = str(row.get("gene_systematic", "")) if pd.notna(row.get("gene_systematic")) else ""
        annot.gene_common = str(row.get("gene_common", "")) if pd.notna(row.get("gene_common")) else ""

        # Numeric fields with safe conversion
        annot.codon_position = int(row.get("codon_position", 0)) if pd.notna(row.get("codon_position")) else 0
        annot.upstream_changes = int(row.get("upstream_changes", 0)) if pd.notna(row.get("upstream_changes")) else 0
        annot.downstream_changes = int(row.get("downstream_changes", 0)) if pd.notna(row.get("downstream_changes")) else 0
        annot.guide_position = int(row.get("guide_position", 0)) if pd.notna(row.get("guide_position")) else 0
        annot.edit_start = int(row.get("edit_start", 0)) if pd.notna(row.get("edit_start")) else 0
        annot.edit_end = int(row.get("edit_end", 0)) if pd.notna(row.get("edit_end")) else 0

        # String fields
        annot.original_aa = str(row.get("original_aa", "")) if pd.notna(row.get("original_aa")) else ""
        annot.new_aa = str(row.get("new_aa", "")) if pd.notna(row.get("new_aa")) else ""
        annot.original_codon = str(row.get("original_codon", "")) if pd.notna(row.get("original_codon")) else ""
        annot.new_codon = str(row.get("new_codon", "")) if pd.notna(row.get("new_codon")) else ""
        annot.strand = str(row.get("strand", "")) if pd.notna(row.get("strand")) else ""
        annot.chromosome = str(row.get("chromosome", "")) if pd.notna(row.get("chromosome")) else ""
        annot.guide_strand = str(row.get("guide_strand", "")) if pd.notna(row.get("guide_strand")) else ""
        annot.pam_sequence = str(row.get("pam_sequence", "")) if pd.notna(row.get("pam_sequence")) else ""
        annot.variant_type = str(row.get("variant_type", "")) if pd.notna(row.get("variant_type")) else ""

        # Boolean
        annot.is_control = bool(row.get("is_control", False))

    else:
        # Oligo not found in harmonized table - check if it's a control
        if "control" in oligo_name.lower() or "background" in oligo_name.lower():
            annot.is_control = True
            annot.variant_type = "control"
            logger.debug(f"Control oligo not in harmonized table: {oligo_name}")
        else:
            logger.warning(f"Oligo not found in harmonized table: {oligo_name}")

    return annot


def parse_oligo_names_to_df(oligo_names: List[str]) -> pd.DataFrame:
    """
    Parse a list of oligo names and return a DataFrame with annotations.

    Uses the pre-parsed harmonized table for efficient lookup.

    Args:
        oligo_names: List of oligo name strings

    Returns:
        DataFrame with annotation columns
    """
    # For large lists, direct merge is more efficient than iterating
    if len(oligo_names) > 100:
        logger.info(f"Annotating {len(oligo_names):,} oligos from harmonized table...")

        # Create DataFrame with oligo names
        df = pd.DataFrame({"oligo_name": oligo_names})

        # Merge with harmonized table
        annot_df = load_oligo_annotations(
            columns=get_annotation_columns_for_bc1_from_screen(include_extras=False),
            rename_to_legacy=True,
            add_is_control=True,
        )

        result = df.merge(annot_df, on="oligo_name", how="left")

        # Report match rate
        n_matched = result["gene_systematic"].notna().sum()
        logger.info(f"  Matched {n_matched:,}/{len(oligo_names):,} oligos")

        return result

    # For small lists, use the parse function
    annotations = [parse_oligo_name(name).to_dict() for name in oligo_names]
    return pd.DataFrame(annotations)


def create_annotation_lookup(oligo_design_files: List[Path] = None) -> pd.DataFrame:
    """
    Create an annotation lookup table from harmonized oligo table.

    This function now uses the pre-parsed harmonized table instead of
    loading raw oligo design files and parsing oligo names.

    Args:
        oligo_design_files: Ignored - kept for backward compatibility.
            The harmonized table already contains all oligo designs.

    Returns:
        DataFrame with oligo_name and annotation columns
    """
    if oligo_design_files:
        logger.info(
            "Note: oligo_design_files argument is deprecated. "
            "Using pre-parsed harmonized table."
        )

    logger.info("Loading annotation lookup from harmonized table...")

    annotations = load_oligo_annotations(
        columns=get_annotation_columns_for_bc1_from_screen(include_extras=True),
        rename_to_legacy=True,
        add_is_control=True,
    )

    logger.info(f"  Loaded {len(annotations):,} oligo annotations")

    return annotations


def annotate_deseq2_results(
    deseq2_df: pd.DataFrame,
    annotation_df: pd.DataFrame = None,
    bc1_ref_df: Optional[pd.DataFrame] = None,
    oligo_col: str = "oligo_name",
) -> pd.DataFrame:
    """
    Annotate DESeq2 results with variant information.

    Args:
        deseq2_df: DESeq2 results DataFrame
        annotation_df: Optional pre-loaded annotation DataFrame.
            If None, loads from harmonized table.
        bc1_ref_df: Optional BC1 reference table for additional info
        oligo_col: Column name containing oligo names

    Returns:
        Annotated DataFrame
    """
    logger.info(f"Annotating DESeq2 results ({len(deseq2_df)} rows)...")

    # Load annotations if not provided
    if annotation_df is None:
        annotation_df = create_annotation_lookup()

    # Merge with annotations
    annotated = deseq2_df.merge(
        annotation_df,
        on=oligo_col,
        how='left'
    )

    # Add BC1 reference info if provided
    if bc1_ref_df is not None and 'bc1' in deseq2_df.columns:
        bc1_cols = ['bc1', 'library_ID', 'yeast_sublibrary', 'bc1_purity', 'num_PCR_replicates']
        available_cols = [c for c in bc1_cols if c in bc1_ref_df.columns]
        if available_cols:
            annotated = annotated.merge(
                bc1_ref_df[available_cols].drop_duplicates(subset=['bc1']),
                on='bc1',
                how='left'
            )

    # Count annotations
    n_annotated = annotated['gene_systematic'].notna().sum()
    pct_annotated = 100 * n_annotated / len(annotated) if len(annotated) > 0 else 0
    logger.info(f"  Annotated {n_annotated}/{len(annotated)} rows ({pct_annotated:.1f}%)")

    # Summary by variant type
    if 'variant_type' in annotated.columns:
        type_counts = annotated['variant_type'].value_counts()
        logger.info("  Variant types:")
        for vtype, count in type_counts.items():
            logger.info(f"    {vtype}: {count}")

    return annotated


def annotate_hits(
    hits_df: pd.DataFrame,
    annotation_df: pd.DataFrame = None,
    oligo_col: str = "oligo_name",
) -> pd.DataFrame:
    """
    Annotate hit calling results with variant information.

    Args:
        hits_df: Hit calling results DataFrame
        annotation_df: Optional pre-loaded annotation DataFrame.
            If None, loads from harmonized table.
        oligo_col: Column name containing oligo names

    Returns:
        Annotated hits DataFrame
    """
    logger.info(f"Annotating hits ({len(hits_df)} rows)...")

    # Load annotations if not provided
    if annotation_df is None:
        annotation_df = create_annotation_lookup()

    annotated = hits_df.merge(
        annotation_df,
        on=oligo_col,
        how='left'
    )

    # Add gene summary for multi-BC1 oligos
    if 'gene_systematic' in annotated.columns:
        # Count unique genes per hit (should be 1, but check)
        n_genes = annotated.groupby(oligo_col)['gene_systematic'].nunique()
        multi_gene = n_genes[n_genes > 1]
        if len(multi_gene) > 0:
            logger.warning(f"  {len(multi_gene)} oligos map to multiple genes (unexpected)")

    return annotated


def create_gene_summary(
    annotated_df: pd.DataFrame,
    group_cols: List[str] = ['gene_systematic', 'gene_common'],
) -> pd.DataFrame:
    """
    Create a gene-level summary from annotated results.

    Args:
        annotated_df: Annotated DataFrame with variant info
        group_cols: Columns to group by for gene summary

    Returns:
        Gene-level summary DataFrame
    """
    # Filter to valid genes
    has_gene = annotated_df['gene_systematic'].notna()
    gene_df = annotated_df[has_gene].copy()

    if gene_df.empty:
        return pd.DataFrame()

    # Aggregate by gene
    agg_dict = {}

    if 'log2FoldChange' in gene_df.columns:
        agg_dict['log2FoldChange'] = ['mean', 'std', 'count']
    if 'padj' in gene_df.columns:
        agg_dict['padj'] = 'min'
    if 'is_hit' in gene_df.columns:
        agg_dict['is_hit'] = ['sum', 'count']

    if not agg_dict:
        # Just count variants per gene
        agg_dict['oligo_name'] = 'count'

    summary = gene_df.groupby(group_cols).agg(agg_dict)
    summary.columns = ['_'.join(col).strip('_') for col in summary.columns]
    summary = summary.reset_index()

    # Add variant type breakdown
    for vtype in ['missense', 'nonsense', 'synonymous']:
        vtype_count = gene_df[gene_df['variant_type'] == vtype].groupby(group_cols).size()
        summary[f'n_{vtype}'] = summary.set_index(group_cols).index.map(vtype_count).fillna(0).astype(int).values

    return summary


def clear_annotation_cache():
    """Clear the annotation lookup cache.

    Useful for testing or when the harmonized table has been updated.
    """
    global _ANNOTATION_CACHE
    _ANNOTATION_CACHE = None
    logger.info("Annotation cache cleared")
