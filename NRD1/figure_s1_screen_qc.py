#!/usr/bin/env python3
"""
figure_s1_screen_qc.py

Supplementary Figure S1: Screen QC for NRD1 saturation mutagenesis manuscript.

Panels:
  A. Barcode coverage per amino acid position, by mutation type (missense / synonymous / PTC)
  B. Replicate correlation: YPD+IAA R1 vs R2 at generation 0 (log10 CPM+1)
  C. Replicate correlation: YPD R1 vs R2 at generation 0 (log10 CPM+1)
  D. Distribution of barcodes per position (histograms by mutation type)
  E. Coverage summary table: % designed positions with ≥1/3/5/10 barcodes

Design context:
  - 574 positions with ≥1 designed missense oligo
  - 575 positions with ≥1 designed synonymous oligo
  - 5 sentinel PTC positions (aa 1, 50, 100, 200, 300) — negative controls only

Author: Kevin R. Roy
Email: kevinrjroy@gmail.com
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats

# ============================================================================
# PATHS
# ============================================================================
BASE_DIR = Path("/path/to")
PROJECT_DIR = BASE_DIR / "projects/NNS/20250813_NRD1_MAGESTIC_screen_reps_1-3"

# Full counts table (373 MB) — used only for barcode metadata (usecols); available on Zenodo.
# For local runs you can point this at the filtered table; coverage numbers will reflect
# only barcodes with mean G0 ≥ 5 CPM (slightly stricter than all-barcode coverage).
COUNTS_TABLE         = PROJECT_DIR / "processed_data/consolidated_normalized_counts.csv"
# Filtered counts table (9 MB, GitHub) — YPD + YPD+IAA columns,
# barcodes with mean YPD G0 ≥ 5 CPM (14,967 barcodes).
# Threshold applied to YPD G0 only — avoids discarding lethal/fitness-defective
# barcodes that deplete rapidly at later timepoints.
# Used for replicate correlation panels (B & C).
FILTERED_COUNTS_TABLE = PROJECT_DIR / "processed_data/NRD1_bc1_counts_YPD_filtered.csv"
DESIGN_TABLE  = BASE_DIR / "software/MAGESTIC/NRD1/NRD1_SpCas9_NGG_guide_donor_200bp_info.tsv"
MANUSCRIPT_DIR = BASE_DIR / "manuscripts/NRD1_saturation_mutagenesis/v3_2026"
OUTPUT_PDF = MANUSCRIPT_DIR / "Figure S1_screen_QC.pdf"
OUTPUT_PNG = MANUSCRIPT_DIR / "Figure S1_screen_QC.png"

# ============================================================================
# CONSTANTS
# ============================================================================
NRD1_LENGTH = 575
GENE_NAME   = "NRD1"

DOMAINS = [
    ("CID",   1,   168, "#1f77b4"),
    ("RS/RE", 169, 301, "#7f7f7f"),
    ("RRM",   302, 464, "#2ca02c"),
    ("PrD",   465, 575, "#ff7f0e"),
]

MUT_ORDER  = ["missense", "synonymous", "PTC"]
MUT_COLORS = {"missense": "#d62728", "synonymous": "#1f77b4", "PTC": "#9467bd"}
MUT_LABELS = {"missense": "Missense", "synonymous": "Synonymous", "PTC": "PTC (sentinel)"}

COVERAGE_THRESHOLDS = [1, 3, 5, 10]


# ============================================================================
# HELPERS
# ============================================================================

def load_designed_positions(design_path: Path) -> dict:
    """
    Return {mutation_type -> frozenset of designed aa positions}.
    mutation_type: 'missense', 'synonymous', 'PTC'
    """
    d = pd.read_csv(design_path, sep="\t", usecols=["aa_nums_str", "ref_aas_str", "alt_aas_str"])
    d = d.dropna(subset=["aa_nums_str", "ref_aas_str", "alt_aas_str"])

    def classify(row):
        if row["alt_aas_str"] == "*":
            return "PTC"
        if row["ref_aas_str"] == row["alt_aas_str"]:
            return "synonymous"
        return "missense"

    d["mutation_type"] = d.apply(classify, axis=1)

    designed = {}
    for mt in MUT_ORDER:
        positions = set(d.loc[d["mutation_type"] == mt, "aa_nums_str"].astype(int))
        designed[mt] = positions
    return designed


def assign_mutation_type(df: pd.DataFrame) -> pd.DataFrame:
    """Derive mutation_type from ref_amino_acid / mutant_amino_acid."""
    df = df.copy()
    ref = df["ref_amino_acid"].astype(str).str.strip()
    alt = df["mutant_amino_acid"].astype(str).str.strip()

    conditions = [alt == "*", ref == alt]
    choices    = ["PTC", "synonymous"]
    df["mutation_type"] = np.select(conditions, choices, default="missense")
    return df


def count_barcodes_per_position(annotated: pd.DataFrame) -> pd.DataFrame:
    """Return DataFrame: amino_acid_position × mutation_type → n_barcodes."""
    return (
        annotated
        .groupby(["amino_acid_position", "mutation_type"])
        .size()
        .reset_index(name="n_barcodes")
    )


def coverage_summary(grp: pd.DataFrame, mut_type: str, designed_positions: set) -> dict:
    """
    For a single mutation type, compute % of designed positions covered at each threshold.
    Positions not in designed_positions are excluded (not shown as '0 barcodes').
    """
    sub = grp[grp["mutation_type"] == mut_type]
    bc_per_pos = sub.set_index("amino_acid_position")["n_barcodes"]

    # Build series over designed positions only
    base = pd.Series(0, index=sorted(designed_positions))
    bc_per_pos = base.add(bc_per_pos.reindex(base.index, fill_value=0), fill_value=0)

    result = {"n_designed": len(designed_positions)}
    for t in COVERAGE_THRESHOLDS:
        result[f">={t}"] = (bc_per_pos >= t).mean() * 100
    return result


def get_rep_corr_cols(df: pd.DataFrame, treatment: str, generation: int):
    r1 = f"{treatment}_R1_Normalized_Count_{generation}"
    r2 = f"{treatment}_R2_Normalized_Count_{generation}"
    if r1 not in df.columns or r2 not in df.columns:
        raise KeyError(
            f"Columns not found: {r1}, {r2}. "
            f"Available: {[c for c in df.columns if 'Normalized' in c][:10]}"
        )
    return r1, r2


# ============================================================================
# PLOTTING
# ============================================================================

def plot_coverage_per_position(ax, grp: pd.DataFrame):
    """Panel A: stacked bar of n_barcodes per position by mutation type."""
    pivot = grp.pivot_table(
        index="amino_acid_position", columns="mutation_type",
        values="n_barcodes", fill_value=0
    )
    for m in MUT_ORDER:
        if m not in pivot.columns:
            pivot[m] = 0
    pivot = pivot[MUT_ORDER]

    x      = pivot.index.values
    bottom = np.zeros(len(x))
    for mut in MUT_ORDER:
        vals = pivot[mut].values.astype(float)
        ax.bar(x, vals, bottom=bottom, width=1.0,
               color=MUT_COLORS[mut], alpha=0.85, label=MUT_LABELS[mut])
        bottom += vals

    for name, start, end, color in DOMAINS:
        ax.axvspan(start - 0.5, end + 0.5, color=color, alpha=0.08, zorder=0)

    ax.set_xlim(0.5, NRD1_LENGTH + 0.5)
    ax.set_xlabel("Amino acid position (NRD1)", fontsize=9)
    ax.set_ylabel("Barcodes (n)", fontsize=9)
    ax.set_title("A  Barcode coverage per position", loc="left", fontsize=10, fontweight="bold")
    ax.legend(fontsize=8, loc="upper right")
    ax.tick_params(labelsize=8)

    ymax = ax.get_ylim()[1]
    for name, start, end, color in DOMAINS:
        mid = (start + end) / 2
        ax.text(mid, ymax * 0.93, name, ha="center", va="top",
                fontsize=7, color=color, fontweight="bold")


def plot_replicate_corr(ax, df: pd.DataFrame, treatment: str, generation: int,
                        panel_label: str):
    """Panels B & C: R1 vs R2 scatter at generation 0 for NRD1 barcodes."""
    r1_col, r2_col = get_rep_corr_cols(df, treatment, generation)

    sub  = df[df["gene_name"] == GENE_NAME].copy()
    x    = np.log10(sub[r1_col] + 1)
    y    = np.log10(sub[r2_col] + 1)
    mask = (x > 0) | (y > 0)
    x, y = x[mask], y[mask]

    r, _ = stats.pearsonr(x, y)

    ax.scatter(x, y, s=4, alpha=0.3, linewidths=0, color=MUT_COLORS["missense"], rasterized=True)

    lim = max(x.max(), y.max()) * 1.05
    ax.plot([0, lim], [0, lim], "k--", linewidth=0.8, alpha=0.6)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("R1 log₁₀(CPM+1)", fontsize=9)
    ax.set_ylabel("R2 log₁₀(CPM+1)", fontsize=9)
    ax.set_title(
        f"{panel_label}  {treatment} G{generation} replicate correlation",
        loc="left", fontsize=10, fontweight="bold"
    )
    ax.text(0.05, 0.92, f"Pearson r = {r:.3f}\nn = {len(x):,}",
            transform=ax.transAxes, fontsize=8, va="top")
    ax.tick_params(labelsize=8)
    ax.set_aspect("equal", adjustable="box")


def plot_bc_distribution(ax, grp: pd.DataFrame, designed: dict):
    """Panel D: histogram of barcodes per designed position."""
    bins = np.arange(0, 31) - 0.5

    for mut in MUT_ORDER:
        # Fill designed positions with 0 if absent
        base = pd.Series(0, index=sorted(designed[mut]))
        sub  = grp[grp["mutation_type"] == mut].set_index("amino_acid_position")["n_barcodes"]
        bc   = base.add(sub.reindex(base.index, fill_value=0), fill_value=0)

        ax.hist(bc.clip(upper=30), bins=bins,
                color=MUT_COLORS[mut], alpha=0.6, label=MUT_LABELS[mut],
                histtype="stepfilled", edgecolor="none")

    ax.set_xlabel("Barcodes per designed position", fontsize=9)
    ax.set_ylabel("Positions (n)", fontsize=9)
    ax.set_title("D  Barcode distribution per designed position",
                 loc="left", fontsize=10, fontweight="bold")
    ax.legend(fontsize=8)
    ax.tick_params(labelsize=8)

    for t in COVERAGE_THRESHOLDS:
        ax.axvline(t, color="gray", linewidth=0.8, linestyle="--", alpha=0.7)


def plot_coverage_table(ax, grp: pd.DataFrame, designed: dict):
    """Panel E: text summary of coverage against designed positions."""
    ax.axis("off")
    ax.set_title("E  Coverage of designed positions", loc="left",
                 fontsize=10, fontweight="bold")

    header = ["Mutation type", "Designed\npositions"] + [f"≥{t} bc" for t in COVERAGE_THRESHOLDS]
    rows   = [header]

    for mut in MUT_ORDER:
        summ = coverage_summary(grp, mut, designed[mut])
        row  = [MUT_LABELS[mut], str(summ["n_designed"])]
        row += [f"{summ[f'>={t}']:.1f}%" for t in COVERAGE_THRESHOLDS]
        rows.append(row)

    table = ax.table(
        cellText=rows[1:],
        colLabels=rows[0],
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.2, 1.7)

    for j in range(len(header)):
        table[(0, j)].set_facecolor("#e8e8e8")


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("Loading design table...")
    designed = load_designed_positions(DESIGN_TABLE)
    for mt, pos in designed.items():
        print(f"  {mt}: {len(pos)} designed positions")

    # Load only metadata columns from full table (fast; avoids reading 373 MB of count data)
    META_COLS = ["bc1", "gene_name", "amino_acid_position", "ref_amino_acid",
                 "mutant_amino_acid", "multiple_donor_bc0_fragments_passing_filter"]
    print("\nLoading barcode metadata from full counts table...")
    df_meta = pd.read_csv(COUNTS_TABLE, usecols=META_COLS, low_memory=False)
    print(f"  {len(df_meta):,} barcodes loaded (metadata only)")

    # Load filtered counts table for replicate correlation panels (YPD + YPD+IAA, mean G0 ≥ 5)
    print("Loading filtered counts table for replicate correlation...")
    df_corr = pd.read_csv(FILTERED_COUNTS_TABLE, low_memory=False)
    print(f"  {len(df_corr):,} barcodes × {len(df_corr.columns)} columns")

    nrd1 = df_meta[df_meta["gene_name"] == GENE_NAME].copy()
    print(f"  {len(nrd1):,} NRD1 barcodes")

    nrd1 = assign_mutation_type(nrd1)
    print("  Mutation type counts:", nrd1["mutation_type"].value_counts().to_dict())

    col      = nrd1["multiple_donor_bc0_fragments_passing_filter"]
    nrd1_hq  = nrd1[col.isin([True, "TRUE", "True", 1])].copy()
    print(f"  {len(nrd1_hq):,} high-quality (single donor) barcodes")

    grp = count_barcodes_per_position(nrd1_hq)
    print(f"  {grp['amino_acid_position'].nunique()} positions with ≥1 HQ barcode")

    print("\nCoverage of designed positions (HQ barcodes):")
    for mut in MUT_ORDER:
        summ = coverage_summary(grp, mut, designed[mut])
        line = f"  {MUT_LABELS[mut]} ({summ['n_designed']} positions): "
        line += ", ".join(f"≥{t}: {summ[f'>={t}']:.1f}%" for t in COVERAGE_THRESHOLDS)
        print(line)

    # ------------------------------------------------------------------
    # Figure layout
    # ------------------------------------------------------------------
    fig = plt.figure(figsize=(18, 14))
    gs  = gridspec.GridSpec(3, 2, height_ratios=[1.2, 1.0, 1.0], hspace=0.45, wspace=0.3)

    ax_cov   = fig.add_subplot(gs[0, :])
    ax_corrB = fig.add_subplot(gs[1, 0])
    ax_corrC = fig.add_subplot(gs[1, 1])
    ax_hist  = fig.add_subplot(gs[2, 0])
    ax_tbl   = fig.add_subplot(gs[2, 1])

    plot_coverage_per_position(ax_cov, grp)

    for ax, treatment, label in [
        (ax_corrB, "YPD+IAA", "B"),
        (ax_corrC, "YPD",     "C"),
    ]:
        try:
            plot_replicate_corr(ax, df_corr, treatment, 0, label)
        except KeyError as e:
            ax.text(0.5, 0.5, str(e), ha="center", va="center",
                    transform=ax.transAxes, fontsize=7, wrap=True)
            ax.set_title(f"{label}  {treatment} G0 replicate correlation",
                         loc="left", fontsize=10, fontweight="bold")

    plot_bc_distribution(ax_hist, grp, designed)
    plot_coverage_table(ax_tbl, grp, designed)

    fig.suptitle(
        "Figure S1 — NRD1 MAGESTIC screen quality control\n"
        "(MAGESTIC 3.0, biological replicates 1–3; "
        "PTC at 5 sentinel positions only: aa 1, 50, 100, 200, 300)",
        fontsize=12, y=0.99,
    )

    MANUSCRIPT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PDF, dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT_PNG, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {OUTPUT_PDF}")
    print(f"Saved: {OUTPUT_PNG}")


if __name__ == "__main__":
    main()
