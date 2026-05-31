"""Harmonize per-nuclease oligo info files into a single annotated TSV.

Ported into the magestic package on 2026-05-30 from the off-tree project
script at projects/QTL/20260116_variant_annotation_harmonization/scripts/
combine_and_annotate_designed_variants.py. CLI entry: magestic-harmonize-oligos
(see pyproject.toml). Auto-gsheet-sync (copy_to_google_drive) is now opt-in
via --gsheet-sync.

Changes from the original project script:
  - Wrapped the original top-level main flow in a main(args=None) function.
  - Added load_v547_sv_prone_oligo_info() + V547 ingestion (677 SV-prone oligos
    joined from the 2024 archived subpool subset + parent Bloom-QTL info file).
  - Added SPS->V547 mapping ("CTGATGCACACGACGTCTTC") to SPS_TO_LIBRARY and
    SPS->oligo_pool ("20240411_SpG_QTL_SV_prone") to SPS_TO_OLIGO_POOL.
  - Added V_aliases column tagging 2025 reclones (V629/V631/V634) as aliases
    of their 2024 originals (V545/V546/V429); harmonized table still keys by
    the 2024 original V.
  - Added argparse-based _cli() entry; copy_to_google_drive() is now opt-in.
"""


from pathlib import Path
import subprocess
import pandas as pd
import numpy as np
from dataclasses import dataclass
from typing import Optional, Tuple

# Configuration
PROJECT = "QTL"
SCREEN_NAME = "20210422_and_20240411_Bloom_et_al_16_strains_QTL"

# Directory setup - using environment variables per Sherlock best practices
# $OAK - archival storage (use for final outputs)
# $SCRATCH - high-bandwidth scratch (use for job I/O if running as SLURM job)
import os
OAK_BASE = Path(os.environ.get("OAK", "/path/to"))
BASE_DIR = OAK_BASE / "Users/kevinroy"
PROJECT_DIR = BASE_DIR / "projects/QTL/20260116_variant_annotation_harmonization"
SCRIPTS_DIR = PROJECT_DIR / "scripts"
KEYFILES_DIR = PROJECT_DIR / "keyfiles"
OUTPUT_DIR = BASE_DIR / "common/annotation_files"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


@dataclass
class FilePaths:
    """File paths for input and output files."""

    spg_oligo_info: str
    spcas9_oligo_info: str
    camp_ras_pka_oligo_info: str
    variants_to_aa: str
    orf_info: str
    score_sv_prediction: str
    output_file: Path
    sample_output_file: Path
    stratified_sample_file: Path
    debug_output_file: Path


# File paths - using BASE_DIR for portability
PATHS = FilePaths(
    spg_oligo_info=str(BASE_DIR / "common/oligo_pools/Bloom_et_al_16_strains_QTL_genes_SNPs_INDELs_with_3_bp_linked_20_bp_from_NGNN_PAM_guide_donor_200_bp_oligos_info.txt"),
    spcas9_oligo_info=str(BASE_DIR / "projects/QTL/20210822_yL274-yL284_QTL_screen/keyfiles/SpCas9_LbCas12a_impLbCas12a_QTL_oligos_info.txt"),
    camp_ras_pka_oligo_info=str(BASE_DIR / "common/oligo_designs/20240514_cAMP_RAS_PKA_saturation_editing_combined_nat_var_codon_change_PTC_syn_control_oligo_info.tsv"),
    variants_to_aa=str(BASE_DIR / "common/published_data/bloom_et_al_2019/Bloom_et_al_16_strains_QTL_genes_5bp_linked_SNPs_INDELs_codon_changes_to_make.tsv"),
    orf_info=str(BASE_DIR / "common/annotation_files/MAGESTIC_background_strain_ORF_info.txt"),
    score_sv_prediction=str(BASE_DIR / "common/published_data/li_et_al_2025/s288c_NGG_PAM_sites_editing_SCORE.tsv"),
    output_file=OUTPUT_DIR / f"{SCREEN_NAME}_harmonized_designed_variant_oligos.tsv",
    sample_output_file=OUTPUT_DIR / f"{SCREEN_NAME}_harmonized_designed_variant_oligos_subsample.tsv",
    stratified_sample_file=OUTPUT_DIR / f"{SCREEN_NAME}_harmonized_designed_variant_oligos_stratified_sample.tsv",
    debug_output_file=OUTPUT_DIR / f"{SCREEN_NAME}_ENA1_oligos.tsv",
)


# ============================================================================
# Section: cAMP/Ras/PKA saturation oligo loader (20240514 design, V574-V597)
# ============================================================================
# The cAMP/Ras/PKA file is a combined 20240514 design covering 24 genes
# (yL410-yL433 trafos -> V574-V597 plasmid libraries). Its schema differs
# from the Bloom-QTL SpG file in three ways:
#   1. Amino-acid columns are suffixed with `_str` (aa_nums_str, ref_aas_str,
#      alt_aas_str, ref_codons_str, alt_codons_str).
#   2. `common_ORF_name` is NaN for nat_var rows; the authoritative gene
#      name is in `ORF_locus` (100% populated for all 4 oligo_types).
#   3. Subpool plumbing columns (subpool_idx, subpool_priming_seq,
#      oligo_type, FWD_PRIMING_SITE etc.) are cAMP-only.
#
# yL410-yL433 -> V574-V597 mapping is from
# library_lineage/LINEAGE_MASTER.tsv (descending oligo count per gene).
CAMP_RAS_PKA_GENE_TO_V = {
    "IRA2": "V574",  "IRA1": "V575",  "RIM15": "V576", "CYR1": "V577",
    "GPB2": "V578",  "GPR1": "V579",  "CDC25": "V580", "GPB1": "V581",
    "GLK1": "V582",  "HXK1": "V583",  "PLC1": "V584",  "RAS1": "V585",
    "SRV2": "V586",  "GPA2": "V587",  "PDE2": "V588",  "RAS2": "V589",
    "RGS2": "V590",  "BCY1": "V591",  "TPK1": "V592",  "TPK3": "V593",
    "PDE1": "V594",  "TFS1": "V595",  "HXK2": "V596",  "TPK2": "V597",
}
CAMP_RAS_PKA_24_GENES = set(CAMP_RAS_PKA_GENE_TO_V.keys())


def load_camp_ras_pka_oligo_info(filepath: str) -> pd.DataFrame:
    """Load the 20240514 cAMP/Ras/PKA saturation oligo info TSV and reshape it
    to the SpG-input schema expected by harmonize_oligo_table().

    Returns a DataFrame with the same columns as spg_oligo_info_df (plus a few
    extras used downstream for per-pool tagging), restricted to the 24
    cAMP/Ras/PKA pathway genes via ORF_locus."""
    df = pd.read_csv(filepath, sep="\t", low_memory=False)

    # Filter to the 24 pathway genes (ORF_locus is 100% populated; common_ORF_name
    # is NaN for all 53,479 nat_var rows so it cannot be used as the discriminator).
    df = df[df["ORF_locus"].isin(CAMP_RAS_PKA_24_GENES)].copy()

    # Rename the _str-suffixed AA columns so downstream compute_aa_change_cols
    # finds them.
    df = df.rename(columns={
        "aa_nums_str": "aa_nums",
        "ref_aas_str": "ref_aas",
        "alt_aas_str": "alt_aas",
        "ref_codons_str": "ref_codons",
        "alt_codons_str": "alt_codons",
    })

    # ensure_oligo_name's fallback infers ORF from oligo_name.split("_")[1].
    # The cAMP oligo_name layout does not match that scheme, so populate ORF/NAME
    # explicitly here. Use ORF_locus (always populated) for both, since
    # common_ORF_name is missing for nat_var rows; downstream tagging will use
    # ORF_locus directly for V-library lookup.
    df["ORF"] = df["ORF_locus"]
    df["NAME"] = df["ORF_locus"]

    return df


# ============================================================================
# Section: V547 SV-prone loader (2024 Twist SpG_SV_prone subpool, 677 oligos)
# ============================================================================
# The V547 SV-prone subpool is a re-cloned SCORE-SV-flagged subset of the
# parent V545/V546 Bloom-QTL design: 677 oligos whose guide+donor content is
# identical to 677 specific oligos in the V545/V546 design, but with a
# V547-specific SPS (CTGATGCACACGACGTCTTC) appended for separate physical
# cloning into V547 (pF988 + cV324).
#
# Source files (Sherlock):
#   subset_filename: the 2-col [oligo_name, oligo_sequence] TSV with V547 SPS
#     lab_shared/.../oligo_array_orders/20240411_Twist_200mer/archived_subpools/
#     Bloom_..._oligos_info_SV_prone_sites_YEF1_oligos_removed.txt
#   parent_info_filename: the rich Bloom-QTL info table that carries
#     guide, donor, PAM, ORF, NAME, individual_variants, etc. (the same one
#     loaded via PATHS.spg_oligo_info).
#
# Strategy: load the 2-col subset, strip the "_SV_prone_site" suffix from each
# oligo_name, INNER-JOIN to the parent info file on oligo_name -> recover all
# annotation columns. Then OVERWRITE oligo_sequence with the subset's
# V547-SPS-tagged version (so downstream SPS extraction yields V547 SPS), and
# re-suffix oligo_name with "_SV_prone_site" so the V547 rows remain
# distinguishable from the V545/V546 rows in the output.
def load_v547_sv_prone_oligo_info(subset_filename: str, parent_info_filename: str) -> pd.DataFrame:
    """Load V547 SV-prone oligos: 2-col subset joined to the rich parent info file.

    Returns a DataFrame in the same schema as the parent SpG info file (so it
    can flow through harmonize_oligo_table() unchanged), with oligo_sequence
    replaced by the V547-SPS version and oligo_name re-suffixed for traceability.
    """
    subset = pd.read_csv(subset_filename, sep="\t")
    if "oligo_sequence" not in subset.columns:
        raise ValueError(
            f"V547 subset file missing oligo_sequence column: {subset_filename}"
        )
    subset["oligo_name_stripped"] = subset["oligo_name"].str.replace(
        r"_SV_prone_site$", "", regex=True
    )

    parent = pd.read_csv(parent_info_filename, sep="\t", low_memory=False)
    n_parent = len(parent)
    n_subset = len(subset)

    merged = parent.merge(
        subset[["oligo_name_stripped", "oligo_sequence"]],
        left_on="oligo_name",
        right_on="oligo_name_stripped",
        how="inner",
        suffixes=("_parent", "_subset"),
    )
    # The subset's oligo_sequence carries the V547 SPS. The parent's carries V545/V546.
    # Use the subset's (V547 SPS), and drop the parent variant.
    if "oligo_sequence_parent" in merged.columns:
        merged = merged.drop(columns=["oligo_sequence_parent"])
    merged = merged.rename(columns={"oligo_sequence_subset": "oligo_sequence"})
    merged = merged.drop(columns=["oligo_name_stripped"])

    # Re-suffix so V547 rows are visually distinct from V545/V546.
    merged["oligo_name"] = merged["oligo_name"] + "_SV_prone_site"

    print(
        f"  V547 SV-prone load: subset={n_subset:,} oligos, parent_info={n_parent:,} oligos, "
        f"merged={len(merged):,} oligos (expect 677)"
    )
    return merged



# ============================================================================
# Section: SCORE-SV Assignment and Imputation
# ============================================================================


def assign_sv_scores(
    variants_df: pd.DataFrame, score_sv_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Merge SV prediction scores from SCORE-SV results with variant data.

    Args:
        variants_df: DataFrame with variant information including 'guide' column
        score_sv_df: DataFrame with SV scores including 'guide' and 'SCORE.SV.minmax' columns

    Returns:
        DataFrame with added 'SV_prediction_score' column
    """
    # Merge on guide sequence
    merged = pd.merge(
        variants_df,
        score_sv_df[["guide", "SCORE.SV.minmax"]],
        on="guide",
        how="left",
    )

    # Rename SCORE column
    merged = merged.rename(columns={"SCORE.SV.minmax": "SV_prediction_score"})

    return merged


def impute_missing_sv_scores(
    variants_df: pd.DataFrame, score_sv_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Impute missing SV scores for non-NGG guides by averaging the two nearest NGG guides.

    Args:
        variants_df: DataFrame with 'chrom', 'PAM_coord', 'SV_prediction_score'
        score_sv_df: DataFrame with SV scores for NGG guides including 'chr', 'pos', 'SCORE.SV.minmax'

    Returns:
        DataFrame with imputed SV_prediction_score values
    """
    print("\nImputing missing SV scores...")

    # Get guides with missing scores
    missing_mask = variants_df["SV_prediction_score"].isna()
    missing_guides = variants_df[missing_mask].copy()

    if len(missing_guides) == 0:
        print("  No missing scores to impute")
        return variants_df

    print(f"  Imputing scores for {len(missing_guides)} guides")

    # Prepare score_sv reference with chromosome and position
    score_sv_ref = score_sv_df[["chr", "pos", "SCORE.SV.minmax"]].copy()
    score_sv_ref = score_sv_ref.dropna(subset=["SCORE.SV.minmax"])

    # Standardize chromosome names: variants use chrI, chrII, etc.; SCORE-SV uses I, II, etc.
    # Add 'chr' prefix to SCORE-SV chromosome names if not already present
    def standardize_chrom(chrom):
        if pd.isna(chrom):
            return chrom
        chrom_str = str(chrom)
        if not chrom_str.startswith("chr"):
            return f"chr{chrom_str}"
        return chrom_str

    score_sv_ref["chr"] = score_sv_ref["chr"].apply(standardize_chrom)

    imputed_scores = []

    # Group by chromosome for efficiency
    for chrom in missing_guides["chrom"].unique():
        chrom_missing = missing_guides[missing_guides["chrom"] == chrom]
        chrom_ref = score_sv_ref[score_sv_ref["chr"] == chrom].sort_values("pos")

        if len(chrom_ref) == 0:
            # No reference guides on this chromosome, use NaN
            imputed_scores.extend([np.nan] * len(chrom_missing))
            continue

        ref_positions = chrom_ref["pos"].values
        ref_scores = chrom_ref["SCORE.SV.minmax"].values

        for _, guide in chrom_missing.iterrows():
            target_pos = guide["PAM_coord"]

            if pd.isna(target_pos):
                imputed_scores.append(np.nan)
                continue

            # Find two nearest positions
            distances = np.abs(ref_positions - target_pos)

            if len(distances) == 1:
                # Only one reference guide on chromosome
                imputed_scores.append(ref_scores[0])
            else:
                # Get two nearest
                nearest_indices = np.argpartition(
                    distances, min(2, len(distances) - 1)
                )[:2]
                nearest_scores = ref_scores[nearest_indices]
                imputed_scores.append(np.mean(nearest_scores))

    # Assign imputed scores
    variants_df.loc[missing_mask, "SV_prediction_score"] = imputed_scores

    # Report imputation statistics
    still_missing = variants_df["SV_prediction_score"].isna().sum()
    successfully_imputed = len(missing_guides) - still_missing

    print(f"  Successfully imputed: {successfully_imputed}")
    print(f"  Still missing (no reference on chromosome): {still_missing}")

    return variants_df


# ============================================================================
# Section: Variant Parsing and Classification
# ============================================================================


def compute_indel_length(individual_variants: str) -> Optional[int]:
    """
    Compute net indel length from individual_variants string.

    Format: position_ref_alt_INDEL or position_ref_alt (for SNPs)

    Args:
        individual_variants: Comma-separated variant strings

    Returns:
        Net indel length (positive for insertions, negative for deletions), or None for SNPs
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


# ============================================================================
# Section: Data Harmonization Functions
# ============================================================================


def ensure_oligo_name(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure oligo_name column exists by using alternative columns if needed."""
    if "oligo_name" not in df.columns:
        for alt_col in ["bdb_oligo_name", "gdb_oligo_name"]:
            if alt_col in df.columns:
                df["oligo_name"] = df[alt_col]
                return df
        raise ValueError(
            "No oligo_name, bdb_oligo_name, or gdb_oligo_name column found."
        )
    return df


def add_orf_and_name(df: pd.DataFrame) -> pd.DataFrame:
    """Infer ORF from oligo_name and populate NAME column."""
    if "ORF" not in df.columns:
        df["ORF"] = None

    inferred_orf = df["oligo_name"].astype(str).str.split("_").str[1]
    df["ORF"] = df["ORF"].where(df["ORF"].notna(), inferred_orf)

    if "NAME" not in df.columns:
        df["NAME"] = None
    df["NAME"] = df["NAME"].where(df["NAME"].notna(), df["ORF"])
    return df


def compute_aa_change_cols(df: pd.DataFrame) -> pd.DataFrame:
    """
    Create aa_change, indel_length, and codon_variant_type columns.

    This function computes:
    1. indel_length: Net length change from individual_variants
    2. aa_change: Human-readable amino acid change string
    3. codon_variant_type: Detailed variant classification
    """

    # First, compute indel_length from individual_variants
    df["indel_length"] = df["individual_variants"].apply(compute_indel_length)

    def mk_aa_change(row):
        ref, nums, alt = row.get("ref_aas"), row.get("aa_nums"), row.get("alt_aas")
        if pd.notna(ref) and pd.notna(nums) and pd.notna(alt):
            try:
                return "_".join(
                    f"{r}{n}{a}"
                    for r, n, a in zip(
                        str(ref).split(","), str(nums).split(","), str(alt).split(",")
                    )
                )
            except Exception:
                return None
        return None

    df["aa_change"] = df.apply(mk_aa_change, axis=1)

    def classify_variant_type(row):
        """
        Classify coding effect of variant with granular details.

        Priority order:
        1. Check if variant is within CDS (if variant_position available)
        2. Handle indels (frameshift vs in-frame)
        3. Handle stop codons (stop_gain, stop_loss)
        4. Handle start codon loss
        5. Synonymous vs missense SNPs
        """
        ref_aas = row.get("ref_aas")
        alt_aas = row.get("alt_aas")
        aa_nums = row.get("aa_nums")
        variant_pos = row.get("variant_position")
        cds_start, cds_end = row.get("start"), row.get("end")
        indel_length = row.get("indel_length")

        # Handle indels FIRST (before checking variant_position)
        # Indels are identified by indel_length being non-null and non-zero
        if pd.notna(indel_length) and indel_length != 0:
            # Check if we have CDS coordinates
            if pd.isna(cds_start) or pd.isna(cds_end):
                return "non_coding"

            # For indels, we need to check if variant overlaps CDS
            # This is better handled by variant_overlaps_CDS later, but for now
            # we can use a simple heuristic: if we have ref_aas/alt_aas, it's coding
            if pd.notna(ref_aas) or pd.notna(alt_aas):
                # Indel affects coding sequence
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
                if (
                    first_pos == 1
                    and ref_str.startswith("M")
                    and not alt_str.startswith("M")
                ):
                    return "start_loss"
            except (ValueError, IndexError):
                pass

        # Synonymous vs missense
        if ref_str == alt_str:
            return "synonymous"
        else:
            return "missense"

    df["codon_variant_type"] = df.apply(classify_variant_type, axis=1)
    return df


def merge_variants_info(
    oligo_df: pd.DataFrame, variants_df: pd.DataFrame
) -> pd.DataFrame:
    """Merge variant annotations into oligo table, avoiding duplicate columns."""
    common_cols = set(oligo_df.columns) & set(variants_df.columns)
    cols_to_drop = [c for c in common_cols if c != "individual_variants"]
    variants_trimmed = variants_df.drop(columns=cols_to_drop)

    return pd.merge(oligo_df, variants_trimmed, how="left", on="individual_variants")


def harmonize_oligo_table(
    oligo_df: pd.DataFrame, variants_df: pd.DataFrame
) -> pd.DataFrame:
    """Apply all harmonization steps to an oligo table."""
    df = ensure_oligo_name(oligo_df)
    df = merge_variants_info(df, variants_df)
    df = add_orf_and_name(df)
    df = compute_aa_change_cols(df)
    return df


def extract_variant_positions(var_str: str) -> tuple:
    """Extract start and end positions from variant string."""
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


def compute_indel_length_change_within_cds(
    var_str: str, cds_start: int, cds_end: int
) -> int:
    """
    Sum the net length change of indels fully contained in a CDS.

    For each variant in the comma-separated var_str, compute var_start and var_end
    (var_end = var_start + len(ref_allele) - 1). If both positions are within the
    CDS boundaries, add len(alt_allele) - len(ref_allele) to the running total.
    For variants overlapping the CDS boundaries, trim to the CDS span before
    computing the delta. Introns are not handled (would require a GFF).
    """
    if (
        pd.isna(var_str)
        or pd.isna(cds_start)
        or pd.isna(cds_end)
        or cds_start > cds_end
    ):
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


# ============================================================================
# Section: Sub-Variant and Amino Acid Level Pre-Processing
# ============================================================================


def extract_shared_aa_missense_change(aa_change: str) -> Optional[str]:
    """
    Extract only missense (non-synonymous, non-silent) changes from aa_change string.

    This extracts amino acid changes where ref != alt, excluding:
    - Synonymous changes (e.g., P28P where ref == alt)
    - Empty or missing values

    Examples:
        "P28P_R29R_G30*" → "G30*"  (only G30* is missense/stop)
        "R29R_G30*" → "G30*"
        "A1V_G2G_K3L" → "A1V_K3L"
        "S10S" → None (all synonymous)

    This allows easy matching of variants that produce the same missense effect
    regardless of linked synonymous changes.

    Args:
        aa_change: Underscore-separated amino acid change string (e.g., "A1V_G2G_K3L")

    Returns:
        Underscore-joined missense changes only, or None if no missense changes
    """
    if pd.isna(aa_change) or str(aa_change).strip() == "":
        return None

    missense_changes = []
    for change in str(aa_change).split("_"):
        change = change.strip()
        if not change:
            continue

        # Parse format: RefAA + Position + AltAA (e.g., "G30*", "A1V")
        # Handle multi-character positions
        if len(change) < 3:
            continue

        ref_aa = change[0]
        alt_aa = change[-1]

        # Include if ref != alt (missense, stop_gain, stop_loss, etc.)
        if ref_aa != alt_aa:
            missense_changes.append(change)

    if not missense_changes:
        return None

    return "_".join(missense_changes)


def identify_sub_variant_relationships(df: pd.DataFrame) -> pd.DataFrame:
    """
    Identify variants that are subsets of other variants at the same locus.

    A sub-variant contains a subset of the mutations in a complete linkage block.
    This is important because:
    - Sub-variants may not represent the full natural variant
    - Frameshift sub-variants of non-frameshift parents are artifacts

    For each gene, variants are compared by their individual_variants sets.
    A variant is marked as a sub-variant if its individual_variants set is a
    proper subset of another variant's set at the same locus.

    Important: Don't link frameshift sub-variants to non-frameshift complete variants.

    Args:
        df: DataFrame with columns: individual_variants, common_gene_name_locus,
            codon_variant_type, variant_overlaps_CDS, oligo_name

    Returns:
        DataFrame with added columns:
        - is_sub_variant: True if this variant's mutations are a subset of another
        - sub_variant_of: List of oligo_names of parent variants (or None)
    """
    print("\nIdentifying sub-variant relationships...")

    df = df.copy()

    # Parse individual_variants into frozensets of mutations
    def parse_variants(var_str):
        if pd.isna(var_str):
            return frozenset()
        return frozenset(str(var_str).split(","))

    df["_variant_set"] = df["individual_variants"].apply(parse_variants)
    df["_variant_set_size"] = df["_variant_set"].apply(len)

    # Initialize columns
    df["is_sub_variant"] = False

    # Store sub_variant_of relationships in a dict to avoid pandas assignment issues
    sub_variant_of_dict = {}

    # Process by gene to limit comparisons
    genes_processed = 0
    sub_variants_found = 0

    for gene, gene_df in df.groupby("common_gene_name_locus"):
        if pd.isna(gene):
            continue

        # Get variants with non-empty sets, sorted by size descending
        gene_variants = gene_df[gene_df["_variant_set_size"] > 0].sort_values(
            "_variant_set_size", ascending=False
        )

        if len(gene_variants) < 2:
            continue

        genes_processed += 1

        # Build lookup of larger variant sets for quick subset checking
        # (idx, oligo_name, variant_set, is_frameshift, overlaps_CDS)
        larger_variants = []
        for idx, row in gene_variants.iterrows():
            is_fs = row.get("codon_variant_type") == "frameshift"
            overlaps_cds = row.get("variant_overlaps_CDS", True)
            larger_variants.append(
                (
                    idx,
                    row["oligo_name"],
                    row["_variant_set"],
                    is_fs,
                    overlaps_cds,
                )
            )

        # For each variant, check if it's a subset of any larger variant
        for i, (idx, oligo, var_set, is_fs, overlaps_cds) in enumerate(larger_variants):
            if len(var_set) == 0:
                continue

            parent_oligos = []
            # Only compare with variants that have more mutations (earlier in sorted list)
            for j in range(i):
                idx2, oligo2, var_set2, is_fs2, _ = larger_variants[j]
                if len(var_set2) <= len(var_set):
                    continue

                if var_set.issubset(var_set2):
                    # Include ALL subset relationships, even frameshift-incompatible ones
                    # The incompatible ones will be flagged as incomplete frameshifts later
                    parent_oligos.append(oligo2)

            if parent_oligos:
                df.at[idx, "is_sub_variant"] = True
                sub_variant_of_dict[idx] = parent_oligos
                sub_variants_found += 1

    # Assign sub_variant_of relationships using map
    df["sub_variant_of"] = df.index.map(lambda idx: sub_variant_of_dict.get(idx, None))

    print(f"  Genes processed: {genes_processed:,}")
    print(f"  Sub-variants identified: {sub_variants_found:,}")

    # Clean up temporary columns
    df = df.drop(columns=["_variant_set", "_variant_set_size"])

    return df


def flag_incomplete_frameshifts(df: pd.DataFrame) -> pd.DataFrame:
    """
    Flag sub-variants that cause frameshifts within a CDS.

    These are typically artifacts of how VCF calls split complex natural variants
    into multiple indels. A sub-variant that causes a frameshift when the full
    variant doesn't represents an "incomplete" or artificial variant.

    Example:
        Full variant: 3bp deletion + 1bp insertion (net -2bp = frameshift)
        Sub-variant 1: 3bp deletion alone (net -3bp = in-frame)
        Sub-variant 2: 1bp insertion alone (net +1bp = frameshift) ← flag this

    Args:
        df: DataFrame with is_sub_variant, sub_variant_of, codon_variant_type,
            variant_overlaps_CDS columns

    Returns:
        DataFrame with added is_incomplete_frameshift column
    """
    print("\nFlagging incomplete frameshift sub-variants...")

    df = df.copy()

    # A sub-variant is an incomplete frameshift if:
    # 1. It is a sub-variant
    # 2. It overlaps CDS
    # 3. It causes a frameshift
    # 4. At least one of its parent variants is NOT a frameshift
    df["is_incomplete_frameshift"] = False

    sub_variant_mask = df["is_sub_variant"] == True
    if not sub_variant_mask.any():
        print("  No sub-variants found, skipping frameshift flagging")
        return df

    # Build lookup of codon_variant_type by oligo_name
    oligo_to_vartype = df.set_index("oligo_name")["codon_variant_type"].to_dict()

    incomplete_count = 0
    for idx in df[sub_variant_mask].index:
        row = df.loc[idx]

        # Skip if not a frameshift or doesn't overlap CDS
        if row["codon_variant_type"] != "frameshift":
            continue
        if not row.get("variant_overlaps_CDS", True):
            continue

        # Check if any parent is NOT a frameshift
        parent_oligos = row.get("sub_variant_of", [])
        if parent_oligos is None:
            continue

        for parent_oligo in parent_oligos:
            parent_vartype = oligo_to_vartype.get(parent_oligo)
            if parent_vartype is not None and parent_vartype != "frameshift":
                df.at[idx, "is_incomplete_frameshift"] = True
                incomplete_count += 1
                break

    print(f"  Incomplete frameshifts flagged: {incomplete_count:,}")

    return df



def main(args=None):
    """Run the harmonization pipeline.

    args is an argparse Namespace (or None for default behavior); BASE_DIR
    and PATHS resolve from environment / module-top defaults.
    """
    """Main execution block."""
    # Load input files
    print("Loading input files...")
    spg_oligo_info_df = pd.read_csv(PATHS.spg_oligo_info, sep="\t")
    spcas9_lbcas12a_oligo_info_df = pd.read_csv(PATHS.spcas9_oligo_info, sep="\t")
    variants_to_aa_info_df = pd.read_csv(PATHS.variants_to_aa, sep="\t")
    orf_info_df = pd.read_csv(PATHS.orf_info, sep="\t")
    # Rename columns to match expected names in the script
    orf_info_df = orf_info_df.rename(columns={
        "systematic_gene_name": "systematic_ORF_name",
        "common_gene_name": "common_ORF_name",
        "CDS_strand": "ORF_strand",
        "CDS_aa_length": "ORF_aa_length"
    })
    score_sv_df = pd.read_csv(PATHS.score_sv_prediction, sep="\t")

    print(f"Loaded {len(spg_oligo_info_df)} SpG oligos")
    print(f"Loaded {len(spcas9_lbcas12a_oligo_info_df)} SpCas9/LbCas12a/impLbCas12a oligos")
    # Report breakdown by nuclease
    nuclease_counts = spcas9_lbcas12a_oligo_info_df['nuclease_PAM'].value_counts()
    for nuclease, count in nuclease_counts.items():
        print(f"  - {nuclease}: {count:,}")
    print(f"Loaded {len(variants_to_aa_info_df)} variant annotations")
    print(f"Loaded {len(score_sv_df)} SCORE-SV predictions")

    # Merge ORF coordinates
    orf_coords = orf_info_df[["systematic_ORF_name", "start", "end"]].drop_duplicates()
    variants_to_aa_info_df = variants_to_aa_info_df.merge(
        orf_coords, on="systematic_ORF_name", how="left"
    )
    variants_to_aa_info_df = variants_to_aa_info_df.drop_duplicates(
        subset=["individual_variants"]
    )

    # Load the cAMP/Ras/PKA saturation oligo info (20240514 SpG design; V574-V597).
    # Filtered to the 24 pathway genes via ORF_locus inside the loader.
    print("\nLoading cAMP/Ras/PKA saturation oligos (24 genes, V574-V597)...")
    camp_oligo_info_df = load_camp_ras_pka_oligo_info(PATHS.camp_ras_pka_oligo_info)
    print(f"Loaded {len(camp_oligo_info_df)} cAMP/Ras/PKA oligos (24 pathway genes)")
    camp_oligo_type_counts = camp_oligo_info_df["oligo_type"].value_counts()
    for ot, n in camp_oligo_type_counts.items():
        print(f"  - {ot}: {n:,}")

    # Load V547 SV-prone oligos (20240411 SpG_SV_prone subpool, 677 oligos).
    # Uses PATHS.spg_oligo_info as the lookup for guide/donor/PAM/ORF/NAME
    # annotations (same parent design, just V547 SPS).
    v547_subset_filename = str(
        BASE_DIR / "lab_shared/SGTC_genome_editing_group/oligo_array_orders/20240411_Twist_200mer/"
        "archived_subpools/Bloom_et_al_16_strains_QTL_genes_SNPs_INDELs_with_3_bp_linked_20_bp_from_"
        "NGNN_PAM_guide_donor_200_bp_oligos_info_SV_prone_sites_YEF1_oligos_removed.txt"
    )
    v547_oligo_info_df = load_v547_sv_prone_oligo_info(v547_subset_filename, PATHS.spg_oligo_info)
    print(f"Loaded {len(v547_oligo_info_df)} V547 SV-prone oligos")


    # Harmonize all three libraries
    print("\nHarmonizing SpG oligos...")
    spg_harmonized = harmonize_oligo_table(spg_oligo_info_df, variants_to_aa_info_df)
    print("\nHarmonizing SpCas9/LbCas12a/impLbCas12a oligos...")
    spcas9_lbcas12a_harmonized = harmonize_oligo_table(
        spcas9_lbcas12a_oligo_info_df, variants_to_aa_info_df
    )
    print("\nHarmonizing cAMP/Ras/PKA oligos...")
    # cAMP's aa_nums/ref_aas/alt_aas are pre-populated from load_camp_ras_pka_oligo_info,
    # so merge_variants_info's left-join won't overwrite them (it drops the variants_df
    # columns that conflict with oligo_df). The Bloom variants_to_aa won't have matching
    # individual_variants strings for cAMP rows, so the merge just adds NaN for
    # variants-only cols (variant_position, start, end) -- those get used in
    # classify_variant_type to flag "non_coding" when missing, which is acceptable for
    # the cAMP path; the lookup table is the gating consumer.
    camp_harmonized = harmonize_oligo_table(camp_oligo_info_df, variants_to_aa_info_df)
    print("\nHarmonizing V547 SV-prone oligos...")
    v547_harmonized = harmonize_oligo_table(v547_oligo_info_df, variants_to_aa_info_df)

    # Align columns before concatenation (all three libraries)
    all_cols = sorted(
        set(spg_harmonized.columns)
        | set(spcas9_lbcas12a_harmonized.columns)
        | set(camp_harmonized.columns)
        | set(v547_harmonized.columns)
    )
    spg_harmonized = spg_harmonized.reindex(columns=all_cols)
    spcas9_lbcas12a_harmonized = spcas9_lbcas12a_harmonized.reindex(columns=all_cols)
    camp_harmonized = camp_harmonized.reindex(columns=all_cols)
    v547_harmonized = v547_harmonized.reindex(columns=all_cols)

    spg_harmonized["category"] = "natural_variants"
    spg_harmonized["nuclease"] = "SpG"

    # The spcas9_lbcas12a_oligo_info_df has a column "nuclease_PAM" but the spg_oligo_info_df does not.
    # Use the "nuclease_PAM" to get the correct nuclease value. Easy to get from nuclease_PAM, as it is always in the format "Nuclease_PAMsequence", e.g. "SpCas9_NGG", "LbCas12a_TTTV", etc.


    def infer_nuclease_from_pam(row):
        nuclease_pam = row.get("nuclease_PAM")
        if pd.notna(nuclease_pam):
            return nuclease_pam.split("_")[0]
        return None


    spcas9_lbcas12a_harmonized["nuclease"] = spcas9_lbcas12a_harmonized.apply(
        infer_nuclease_from_pam, axis=1
    )

    # Per-pool tagging for cAMP/Ras/PKA rows: only the fields the downstream SPS-based
    # block does NOT recompute. SPS / guide_donor_bc0_plasmid_library / oligo_pool are
    # filled later (after combined["SPS"] is extracted from oligo_sequence), so don't
    # set them here -- they would be clobbered.
    camp_harmonized["nuclease"] = "SpG"
    camp_harmonized["category"] = camp_harmonized["oligo_type"]
    v547_harmonized["nuclease"] = "SpG"
    v547_harmonized["category"] = "natural_variants"
    # guide_nuclease_PAM is computed downstream by infer_guide_nuclease_pam, which keys
    # off (nuclease, PAM) and correctly returns SpG_NGNG vs SpG_NGNH for these rows.

    # Concatenate (cAMP appended after the two Bloom-QTL libraries)
    print("\nConcatenating all three libraries...")
    combined = pd.concat(
        [spg_harmonized, spcas9_lbcas12a_harmonized, camp_harmonized, v547_harmonized], ignore_index=True
    )
    combined = combined.reindex(
        columns=list(spcas9_lbcas12a_oligo_info_df.columns)
        + [c for c in combined.columns if c not in spcas9_lbcas12a_oligo_info_df.columns]
    )

    # Assign SCORE-SV predictions
    print("\nAssigning SCORE-SV predictions...")
    print(f"Before SV assignment: {len(combined)} oligos")
    combined = assign_sv_scores(combined, score_sv_df)
    print(f"  NGG guides with SV scores: {(~combined['SV_prediction_score'].isna()).sum()}")
    print(f"  Missing scores: {combined['SV_prediction_score'].isna().sum()}")

    # Impute missing SV scores for non-NGG guides
    combined = impute_missing_sv_scores(combined, score_sv_df)
    print(f"\nAfter imputation:")
    print(f"  Total with SV scores: {(~combined['SV_prediction_score'].isna()).sum()}")
    print(f"  Still missing: {combined['SV_prediction_score'].isna().sum()}")

    # Rename variant_type to natural_variant_type to clarify it applies only to natural variants
    # (has values of LINKED, SNP, INDEL for natural_variants category; NA for other categories)
    if "variant_type" in combined.columns:
        combined = combined.rename(columns={"variant_type": "natural_variant_type"})

    # Map gene names
    gene_to_sysorf = dict(
        zip(orf_info_df["common_ORF_name"], orf_info_df["systematic_ORF_name"])
    )
    combined["ORF"] = combined.apply(
        lambda row: (
            gene_to_sysorf.get(row["ORF"], row["ORF"])
            if row["category"] == "natural_variants" and row["nuclease"] == "SpCas9"
            else row["ORF"]
        ),
        axis=1,
    )

    combined = combined.rename(columns={"ORF": "systematic_gene_name_locus"})
    sysorf_to_gene = dict(
        zip(orf_info_df["systematic_ORF_name"], orf_info_df["common_ORF_name"])
    )
    combined["common_gene_name_locus"] = combined["systematic_gene_name_locus"].map(
        sysorf_to_gene
    )

    # Assess variant overlap with CDS
    combined["variant_overlaps_CDS"] = False
    combined["overlapping_systematic_gene_name"] = None
    combined["overlapping_common_gene_name"] = None

    # Fast path for controls
    control_mask = combined["individual_variants"].isna()
    combined.loc[control_mask, "overlapping_systematic_gene_name"] = combined.loc[
        control_mask, "systematic_gene_name_locus"
    ]
    combined.loc[control_mask, "overlapping_common_gene_name"] = combined.loc[
        control_mask, "common_gene_name_locus"
    ]
    combined.loc[control_mask, "variant_overlaps_CDS"] = True
    # Vectorized CDS overlap checks for variants
    variant_mask = ~control_mask
    if variant_mask.any():
        combined[["_start", "_end"]] = combined.loc[
            variant_mask, "individual_variants"
        ].apply(lambda x: pd.Series(extract_variant_positions(x)))

        for _, orf in orf_info_df.iterrows():
            mask = (
                variant_mask
                & combined["_start"].notna()
                & (combined["chrom"] == orf["chrom"])
                & (combined["_end"] >= orf["start"])
                & (combined["_start"] <= orf["end"])
            )
            combined.loc[mask, "variant_overlaps_CDS"] = True
            combined.loc[mask, "overlapping_systematic_gene_name"] = orf[
                "systematic_ORF_name"
            ]
            combined.loc[mask, "overlapping_common_gene_name"] = orf["common_ORF_name"]
            # For the portion of overlapping INDEL variants that reside entirely within CDS,
            # determine if they are frameshift or in-frame indels.
            indel_mask = (
                mask
                & combined["indel_length"].notna()
                & (combined["indel_length"] != 0)
                & (combined["_start"] >= orf["start"])
                & (combined["_end"] <= orf["end"])
            )
            combined.loc[indel_mask, "codon_variant_type"] = combined.loc[indel_mask].apply(
                lambda row: (
                    "frameshift" if abs(row["indel_length"]) % 3 != 0 else "in-frame_indel"
                ),
                axis=1,
            )

        combined = combined.drop(columns=["_start", "_end"])

    combined["category"].unique()
    #
    # When category is one of these: PTC_del_syn_controls, PDR5_controls, DHY214_controls,
    # the variant indeed overlaps the CDS of the gene being targeted.
    # use the systematic_gene_name_locus and common_gene_name_locus columns for these.
    control_categories = [
        "PTC_del_syn_controls",
        "PDR5_controls",
        "DHY214_controls",
    ]
    control_cat_mask = combined["category"].isin(control_categories)
    combined.loc[control_cat_mask, "variant_overlaps_CDS"] = True
    combined.loc[control_cat_mask, "overlapping_systematic_gene_name"] = combined.loc[
        control_cat_mask, "systematic_gene_name_locus"
    ]
    combined.loc[control_cat_mask, "overlapping_common_gene_name"] = combined.loc[
        control_cat_mask, "common_gene_name_locus"
    ]

    # Classify complete CDS deletions in PTC_del_syn_controls
    # These are oligos that delete the entire ORF (negative indel_length covering whole CDS)
    # Identify by: category is PTC_del_syn_controls, has indel_length (is a deletion)
    print("\nDebug: Checking for complete CDS deletions...")
    print(
        f"  Total PTC_del_syn_controls: {(combined['category'] == 'PTC_del_syn_controls').sum():,}"
    )
    if (combined["category"] == "PTC_del_syn_controls").any():
        ptc_rows = combined[combined["category"] == "PTC_del_syn_controls"]
        print(
            f"  PTC rows with individual_variants: {ptc_rows['individual_variants'].notna().sum():,}"
        )
        print(f"  individual_variants value counts:")
        for val, count in ptc_rows["individual_variants"].value_counts().head(5).items():
            print(f"    '{val}': {count:,}")

    # Identify complete CDS deletions by the 'deletion' label
    ptc_del_mask = (combined["category"] == "PTC_del_syn_controls") & (
        combined["individual_variants"] == "deletion"
    )
    combined.loc[ptc_del_mask, "codon_variant_type"] = "complete_CDS_deletion"
    print(
        f"\nClassified {ptc_del_mask.sum():,} complete CDS deletions in PTC_del_syn_controls"
    )

    # Note: natural_variant_type is kept for natural_variants category (values: LINKED, SNP, INDEL)
    # It will be NA for other categories (codon_changes, DHY214_controls, etc.)

    # Fill missing values
    combined["first_aa_num"] = combined["first_aa_num"].where(
        combined["first_aa_num"].notna(),
        combined["aa_nums"].apply(lambda x: str(x).split(",")[0] if pd.notna(x) else None),
    )

    # Rename nuclease_PAM to guide_nuclease_PAM to avoid confusion with strain background nuclease
    # This column refers to the nuclease + PAM compatibility of the designed guide
    combined = combined.rename(columns={"nuclease_PAM": "guide_nuclease_PAM"})


    def infer_guide_nuclease_pam(row):
        """
        Infer guide_nuclease_PAM from nuclease and PAM columns.

        Returns the nuclease_PAM format used in oligo design:
        - SpG_NGNG, SpG_NGNH for SpG nuclease
        - SpCas9_NGG for SpCas9 nuclease
        - LbCas12a_TTTV for LbCas12a nuclease
        - impLbCas12a_nonTTTV for impLbCas12a nuclease
        """
        nuclease = row.get("nuclease")
        pam = row.get("PAM")

        if nuclease == "SpG":
            if pd.notna(pam) and len(str(pam)) == 4:
                pam_str = str(pam)
                if pam_str[1] == "G" and pam_str[3] == "G":
                    return "SpG_NGNG"
                else:
                    return "SpG_NGNH"
            return "SpG_NGNH"  # Default for SpG if PAM not parseable
        elif nuclease == "SpCas9":
            return "SpCas9_NGG"
        elif nuclease == "LbCas12a":
            return "LbCas12a_TTTV"
        elif nuclease == "impLbCas12a":
            return "impLbCas12a_nonTTTV"
        else:
            return None


    combined["guide_nuclease_PAM"] = combined["guide_nuclease_PAM"].where(
        combined["guide_nuclease_PAM"].notna(),
        combined.apply(infer_guide_nuclease_pam, axis=1),
    )

    # ============================================================================
    # SPS, guide_donor_bc0_plasmid_library, and oligo_pool columns
    # ============================================================================

    # Extract SPS (Sublibrary Primer Sequence) directly from oligo_sequence
    # SPS is the last 20 bases of the oligo, uppercased
    print("\nExtracting SPS from oligo sequences...")
    combined["SPS"] = combined["oligo_sequence"].str[-20:].str.upper()

    # SPS to guide_donor_bc0_plasmid_library mapping (V format, not cV)
    # Based on actual oligo pool designs:
    # - 2024 SpG pool: V520 (NGNG), V521 (NGNH)
    # - 2021 SpCas9/LbCas12a pool: V416 (SpCas9), V419 (LbCas12a, impLbCas12a)
    SPS_TO_LIBRARY = {
        # SPS -> step-1 V library (the V that holds the oligos), per PROJECT.md
        # "MAGESTIC Library Construction Reference / Step 1: Plasmid Gibson Assembly".
        # Construction chain is: oligo pool -> step-1 V (pF amplicon + cV backbone via Gibson)
        # -> step-2 V (step-1 V cut + step-2 cV insert via ligation) -> step-2 cV (step-2 V + I-SceI)
        # + bc1 pF -> yeast trafo. The legacy bc1_reference_table.tsv labels with the
        # step-1 V (e.g. V545/V546/V429/V454), so harmonized should match.
        # Prior mapping used V520/V521/V416/V419 which are unrelated plasmids (V520
        # in V.tsv is a 2023 impLbCas12a step-2 build; V416 is a 2021 backbone) --
        # corrected 2026-05-30.
        "CTTGACTCGACAGTGAGAGC": "V545",   # SpG_NGNG step-1 V (2024 Twist, 10,439 oligos)
        "CGTATGGAGATGACGTCAGG": "V546",   # SpG_NGNH step-1 V (2024 Twist, 20,026 oligos)
        "CGATAGACGACTGGACAGCA": "V429",   # SpCas9_NGG step-1 V (2021 Twist, 17,782 oligos)
        "GAACTACGTGTCTCGTCAGG": "V454",   # LbCas12a_TTTV step-1 V (2021 Twist, WT LbCas12a scaffold)
        "GTTATGCTGGTCCTAGGTCG": "V455",   # impLbCas12a_nonTTTV step-1 V (2021 Twist, modified scaffold)
        "CTGATGCACACGACGTCTTC": "V547",   # SpG_SV_prone step-1 V (2024 Twist, 677 oligos; SCORE-SV-flagged subset of V545/V546 parent design)
        # NOTE: 2025 reclones V629/V631/V634 share SPS with V545/V546/V429 respectively
        # (same oligo content, different backbone pS1439). Disambiguation by SPS alone
        # is not possible -- downstream consumers using the screen-level keyfile know
        # which physical library was transformed.
    }

    # SPS to oligo_pool mapping (year + description of the Twist oligo pool)
    SPS_TO_OLIGO_POOL = {
        "CTTGACTCGACAGTGAGAGC": "20240411_SpG_QTL",
        "CGTATGGAGATGACGTCAGG": "20240411_SpG_QTL",
        "CTGATGCACACGACGTCTTC": "20240411_SpG_QTL_SV_prone",
        "CGATAGACGACTGGACAGCA": "20210422_SpCas9_LbCas12a_QTL",
        "GAACTACGTGTCTCGTCAGG": "20210422_SpCas9_LbCas12a_QTL",
        "GTTATGCTGGTCCTAGGTCG": "20210422_SpCas9_LbCas12a_QTL",
    }

    # Map SPS to guide_donor_bc0_plasmid_library
    combined["guide_donor_bc0_plasmid_library"] = combined["SPS"].map(SPS_TO_LIBRARY)

    # Map SPS to oligo_pool
    combined["oligo_pool"] = combined["SPS"].map(SPS_TO_OLIGO_POOL)

    # Fill cAMP/Ras/PKA rows' guide_donor_bc0_plasmid_library and oligo_pool from
    # ORF_locus -> V (each gene has its own per-subpool V; SPS varies by gene so the
    # SPS_TO_* dicts above would need 24 entries -- using ORF_locus is the same data
    # in one line). ORF_locus is dropped by the cleanup block below, so this MUST
    # run before that block.
    camp_mask = combined["ORF_locus"].isin(CAMP_RAS_PKA_24_GENES) if "ORF_locus" in combined.columns else pd.Series(False, index=combined.index)
    if camp_mask.any():
        combined.loc[camp_mask, "guide_donor_bc0_plasmid_library"] = (
            combined.loc[camp_mask, "ORF_locus"].map(CAMP_RAS_PKA_GENE_TO_V)
        )
        combined.loc[camp_mask, "oligo_pool"] = "20240514_cAMP_Ras_PKA_saturation"
        print(f"  Tagged {int(camp_mask.sum())} cAMP/Ras/PKA rows with V574-V597 + oligo_pool '20240514_cAMP_Ras_PKA_saturation'")

    # ============================================================================
    # V_aliases column: list 2025 reclones that share SPS with 2024 originals.
    # Strategy A (selected 2026-05-30): keep SPS_TO_LIBRARY pointing at 2024
    # originals (V545/V546/V429); record 2025 reclones (V629/V631/V634) as
    # pipe-delimited aliases so downstream consumers can remap if they need to
    # distinguish backbone (pS1380 vs pS1439). Screens that physically used a
    # 2025 reclone do the remap in their own keyfile join.
    # ============================================================================
    V_LIBRARY_ALIASES = {
        "V545": "V545|V629",  # SpG_NGNG: 2024 V545 (pS1380) / 2025 reclone V629 (pS1439)
        "V546": "V546|V631",  # SpG_NGNH: 2024 V546 (pS1380) / 2025 reclone V631 (pS1439)
        "V429": "V429|V634",  # SpCas9_NGG: 2024 V429 (pS1380) / 2025 reclone V634 (pS1439)
    }
    combined["V_aliases"] = combined["guide_donor_bc0_plasmid_library"].map(V_LIBRARY_ALIASES).fillna(
        combined["guide_donor_bc0_plasmid_library"]
    )
    n_with_aliases = (combined["V_aliases"] != combined["guide_donor_bc0_plasmid_library"]).sum()
    print(f"  V_aliases column: {n_with_aliases:,} rows tagged with 2025-reclone aliases (V629/V631/V634)")


    # Report SPS mapping results
    print(f"  SPS values: {combined['SPS'].value_counts().to_dict()}")
    print(f"  guide_donor_bc0_plasmid_library: {combined['guide_donor_bc0_plasmid_library'].value_counts().to_dict()}")
    print(f"  oligo_pool: {combined['oligo_pool'].value_counts().to_dict()}")

    # Check for unmapped SPS values (excludes cAMP rows we just patched)
    unmapped_sps = combined[combined["guide_donor_bc0_plasmid_library"].isna()]["SPS"].unique()
    if len(unmapped_sps) > 0:
        print(f"  WARNING: {len(unmapped_sps)} unmapped SPS values: {unmapped_sps[:5]}")

    # Clean up by removing these cols:
    cols_to_remove = [
        "samples.1",
        "common_ORF_name",
        "systematic_ORF_name",
        "ORF_strand",
        "ORF_locus",
        "NAME_locus",
        "ORF_aa_length",
        "in_frame_ATGs",
        "genomic_target",
        "NAME",
        "end",
        "gene_name",
        "start",
    ]
    combined = combined.drop(columns=[c for c in cols_to_remove if c in combined.columns])

    # ============================================================================
    # Sub-Variant and Amino Acid Level Pre-Processing
    # ============================================================================

    # Compute shared_aa_missense_change (only missense changes from aa_change)
    print("\nComputing shared amino acid missense changes...")
    combined["shared_aa_missense_change"] = combined["aa_change"].apply(
        extract_shared_aa_missense_change
    )
    n_with_shared = combined["shared_aa_missense_change"].notna().sum()
    print(f"  Variants with shared_aa_missense_change: {n_with_shared:,}")

    # Identify sub-variant relationships (variants that are subsets of linkage blocks)
    combined = identify_sub_variant_relationships(combined)

    # Flag incomplete frameshift sub-variants (artifacts of VCF splitting)
    combined = flag_incomplete_frameshifts(combined)

    # Write outputs
    # Write outputs
    print("\nWriting output files...")
    combined.to_csv(PATHS.output_file, sep="\t", index=False)
    combined.sample(n=min(100, len(combined)), random_state=42).to_csv(
        PATHS.sample_output_file, sep="\t", index=False
    )

    # Create stratified sample: up to 5 rows per (guide_nuclease_PAM, category, codon_variant_type, natural_variant_type) combination
    print("\nCreating stratified sample...")
    strat_cols = ["guide_nuclease_PAM", "category", "codon_variant_type", "natural_variant_type"]
    # Fill NA values for grouping
    combined_for_strat = combined.copy()
    for col in strat_cols:
        if col in combined_for_strat.columns:
            # Simply select the column - if there are duplicates, pandas will return a DataFrame
            temp = combined_for_strat[col]
            if isinstance(temp, pd.DataFrame):
                # Multiple columns with same name - take first
                temp_series = temp.iloc[:, 0]
            else:
                # Single column - it's already a Series
                temp_series = temp
            # Convert to string and replace 'nan' and 'None' strings with '_NA_'
            temp_series = temp_series.astype(str).replace(["nan", "None", ""], "_NA_")
            combined_for_strat[f"_strat_{col}"] = temp_series
        else:
            combined_for_strat[f"_strat_{col}"] = "_missing_col_"

    strat_group_cols = [f"_strat_{col}" for col in strat_cols]
    stratified_sample = (
        combined_for_strat.groupby(strat_group_cols, dropna=False)
        .apply(lambda x: x.head(5), include_groups=False)
        .reset_index(drop=True)
    )
    # Drop the temporary stratification columns
    stratified_sample = stratified_sample.drop(
        columns=[c for c in stratified_sample.columns if c.startswith("_strat_")]
    )
    stratified_sample.to_csv(PATHS.stratified_sample_file, sep="\t", index=False)

    # Report stratified sample statistics
    n_combos = combined_for_strat.groupby(strat_group_cols, dropna=False).ngroups
    print(
        f"  Unique (guide_nuclease_PAM, category, codon_variant_type, natural_variant_type) combinations: {n_combos}"
    )
    print(f"  Stratified sample rows: {len(stratified_sample):,}")
    print(f"  Stratified sample file: {PATHS.stratified_sample_file}")

    # Write out all lines with "ENA1" in oligo_name for debugging
    combined[combined["oligo_name"].str.contains("ENA1", na=False)].to_csv(
        PATHS.debug_output_file, sep="\t", index=False
    )

    # Print comprehensive summary
    print(f"\n{'='*80}")
    print("HARMONIZATION SUMMARY")
    print(f"{'='*80}")
    print(f"Output file: {PATHS.output_file}")
    print(f"Total oligos: {combined.shape[0]:,}")
    print(f"Total columns: {combined.shape[1]}")
    if "oligo_name" in combined.columns:
        print(f"Unique oligo names: {combined['oligo_name'].nunique():,}")

    # Variant type distribution
    if "codon_variant_type" in combined.columns:
        print(f"\nVariant type distribution:")
        for vtype, count in combined["codon_variant_type"].value_counts().head(10).items():
            print(f"  {vtype:20s}: {count:7,} ({100*count/len(combined):5.1f}%)")

        # Count frameshifts specifically
        frameshift_count = (combined["codon_variant_type"] == "frameshift").sum()
        print(f"\n  Total frameshifts: {frameshift_count:,}")

    # Indel statistics
    if "indel_length" in combined.columns:
        indels = combined[combined["indel_length"].notna()]
        if len(indels) > 0:
            print(f"\nIndel statistics:")
            print(f"  Total indels: {len(indels):,}")
            print(f"  Mean length: {indels['indel_length'].mean():.2f} bp")
            print(f"  Frameshift indels: {((indels['indel_length'] % 3) != 0).sum():,}")
            print(f"  In-frame indels: {((indels['indel_length'] % 3) == 0).sum():,}")

    # SV score statistics
    if "SV_prediction_score" in combined.columns:
        print(f"\nSCORE-SV statistics:")
        print(f"  Total with scores: {combined['SV_prediction_score'].notna().sum():,}")
        print(f"  Missing scores: {combined['SV_prediction_score'].isna().sum():,}")
        if combined["SV_prediction_score"].notna().any():
            print(f"  Mean score: {combined['SV_prediction_score'].mean():.4f}")
            print(f"  High SV (>0.15): {(combined['SV_prediction_score'] > 0.15).sum():,}")

    if {"ref_aas", "aa_nums", "alt_aas"}.issubset(combined.columns):
        print(f"\nAmino acid annotations:")
        print(f"  Oligos with aa_change: {combined['aa_change'].notna().sum():,}")

    # Sub-variant statistics
    if "is_sub_variant" in combined.columns:
        print(f"\nSub-variant analysis:")
        n_sub = combined["is_sub_variant"].sum()
        n_incomplete_fs = combined.get("is_incomplete_frameshift", pd.Series([False])).sum()
        n_shared_aa = combined["shared_aa_missense_change"].notna().sum()
        print(f"  Total sub-variants: {n_sub:,}")
        print(f"  Incomplete frameshifts (artifacts): {n_incomplete_fs:,}")
        print(f"  Variants with shared_aa_missense_change: {n_shared_aa:,}")

        # Show some examples of sub-variant relationships
        if n_sub > 0:
            print(f"\n  Example sub-variant relationships (first 5):")
            sub_examples = combined[combined["is_sub_variant"]].head(5)
            for _, row in sub_examples.iterrows():
                gene = row.get("common_gene_name_locus", "?")
                oligo = row.get("oligo_name", "?")
                parent = row.get("sub_variant_of", [])
                print(f"    {gene}: {oligo} is sub-variant of {parent}")

    print(f"{'='*80}\n")




def copy_to_google_drive():
    """Sync directory to Google Drive using rclone."""
    folder = "projects/QTL/20260116_variant_annotation_harmonization/"
    drive_dir = f"stanford_google_drive:sequencing_data/{folder}"
    remote_dir = f"/path/to/{folder}"

    rclone_cmd = [
        "rclone",
        "copy",
        remote_dir,
        drive_dir,
        "--progress",
        "--checkers",
        "24",
        "--transfers",
        "24",
        "--update",
    ]

    try:
        subprocess.run(rclone_cmd, check=True)
        print(f"Successfully synced {remote_dir} to {drive_dir}")
    except subprocess.CalledProcessError as e:
        print(f"Error running rclone: {e}")




def _cli():
    """Argparse entry point for the magestic-harmonize-oligos CLI."""
    import argparse
    parser = argparse.ArgumentParser(
        description=(
            "Harmonize per-nuclease oligo info files into a single annotated TSV. "
            "Inputs/outputs are resolved from BASE_DIR (env var OAK or --oak), see "
            "module-top PATHS definition. Run with --gsheet-sync to push the result "
            "directory to the lab Google Drive via rclone."
        )
    )
    parser.add_argument(
        "--gsheet-sync",
        action="store_true",
        help=(
            "After successful harmonization, sync the project dir to "
            "stanford_google_drive: via rclone. Default off (was on in the original "
            "off-tree project script)."
        ),
    )
    args = parser.parse_args()
    main(args)
    if args.gsheet_sync:
        copy_to_google_drive()


if __name__ == "__main__":
    _cli()
