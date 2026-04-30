#!/usr/bin/env python3
"""
figure1_fitness_heatmap.py

Publication-quality fitness heatmap for Figure 1 of the NRD1 saturation
mutagenesis manuscript (MAGESTIC 3.0).

Shows median-centered log2FC slope (YPD+IAA condition) for each amino acid
position × mutation type combination, with domain annotations.

Author: Kevin R. Roy
Email: kevinrjroy@gmail.com
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm

# ============================================================================
# PATHS
# ============================================================================
BASE_DIR = Path("/path/to")
PROJECT_DIR = BASE_DIR / "projects/NNS/20250813_NRD1_MAGESTIC_screen_reps_1-3"

SLOPE_TABLE = PROJECT_DIR / "processed_data/fitness_results/slope_table_log2FC.csv"
HITS_TABLE  = PROJECT_DIR / "processed_data/fitness_results/hits_summary.csv"
MANUSCRIPT_DIR = BASE_DIR / "manuscripts/NRD1_saturation_mutagenesis/v3_2026"
OUTPUT_PDF = MANUSCRIPT_DIR / "Figure 1_fitness_heatmap.pdf"

# ============================================================================
# DOMAIN BOUNDARIES  (from NRD1_domains.tsv, Mar 2026)
# ============================================================================
NRD1_LENGTH = 575
DOMAINS = [
    ("CID",   1,   168, "#1f77b4"),   # blue
    ("RS/RE", 169, 301, "#7f7f7f"),   # gray
    ("RRM",   302, 464, "#2ca02c"),   # green
    ("PrD",   465, 575, "#ff7f0e"),   # orange
]

TREAT_YPD_IAA = "YPD+IAA"

# Colormap: white at 0, red negative (depleted), blue positive
CMAP = "RdBu"
VMIN, VMAX = -3.0, 1.0          # slope_adj range from screen
VCENTER = 0.0

MUT_ORDER  = ["missense", "synonymous", "PTC"]
MUT_LABELS = ["Missense", "Synonymous", "PTC"]
MUT_COLORS = ["#d62728", "#1f77b4", "#9467bd"]   # row label colors


def build_matrix(slope_df: pd.DataFrame) -> dict[str, np.ndarray]:
    """
    Build a dict of {mutation_type: 1×NRD1_LENGTH array} for the YPD+IAA treatment.
    NaN where no data.
    """
    ypd = slope_df[slope_df["treatment"] == TREAT_YPD_IAA].copy()
    matrices = {}
    for mut_type in MUT_ORDER:
        arr = np.full(NRD1_LENGTH, np.nan)
        sub = ypd[ypd["mutation_type"] == mut_type]
        for _, row in sub.iterrows():
            pos = int(row["amino_acid_position"]) - 1   # 0-indexed
            if 0 <= pos < NRD1_LENGTH:
                arr[pos] = row["slope_adj"]
        matrices[mut_type] = arr
    return matrices


def main():
    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    slope_df = pd.read_csv(SLOPE_TABLE)
    hits_df  = pd.read_csv(HITS_TABLE)

    matrices = build_matrix(slope_df)

    # Hit positions (missense + PTC, for tick marks on heatmap)
    hit_positions = sorted(hits_df["amino_acid_position"].dropna().astype(int).unique())

    # ------------------------------------------------------------------
    # Figure layout
    # ------------------------------------------------------------------
    fig = plt.figure(figsize=(18, 5))

    # Rows: [heatmap_missense, heatmap_synonymous, heatmap_PTC, domain_bar]
    # Heights proportional
    gs = fig.add_gridspec(
        4, 2,
        height_ratios=[1, 1, 0.4, 0.3],
        width_ratios=[20, 1],
        hspace=0.04,
        wspace=0.02,
    )

    ax_miss = fig.add_subplot(gs[0, 0])
    ax_syn  = fig.add_subplot(gs[1, 0])
    ax_ptc  = fig.add_subplot(gs[2, 0])
    ax_dom  = fig.add_subplot(gs[3, 0])
    ax_cbar = fig.add_subplot(gs[0:3, 1])

    axes_heatmap = [ax_miss, ax_syn, ax_ptc]
    norm = TwoSlopeNorm(vmin=VMIN, vcenter=VCENTER, vmax=VMAX)

    # ------------------------------------------------------------------
    # Heatmap rows
    # ------------------------------------------------------------------
    x = np.arange(1, NRD1_LENGTH + 1)
    im = None

    for ax, mut_type, label, lcolor in zip(
        axes_heatmap, MUT_ORDER, MUT_LABELS, MUT_COLORS
    ):
        arr = matrices[mut_type].reshape(1, -1)
        im = ax.imshow(
            arr,
            aspect="auto",
            cmap=CMAP,
            norm=norm,
            extent=[0.5, NRD1_LENGTH + 0.5, -0.5, 0.5],
            interpolation="nearest",
        )

        ax.set_yticks([0])
        ax.set_yticklabels([label], fontsize=9, color=lcolor, fontweight="bold")
        ax.tick_params(axis="y", left=False, pad=2)
        ax.set_xlim(0.5, NRD1_LENGTH + 0.5)
        ax.set_xticks([])

    # Mark hit positions with ticks on the bottom of each row
    for ax in axes_heatmap:
        for pos in hit_positions:
            ax.axvline(pos, color="black", linewidth=0.6, alpha=0.8, ymin=0, ymax=0.12)

    # ------------------------------------------------------------------
    # Colorbar
    # ------------------------------------------------------------------
    cbar = fig.colorbar(im, cax=ax_cbar, orientation="vertical")
    cbar.set_label("Adjusted slope\n(log₂FC / generation)", fontsize=8)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_ticks([VMIN, -2, -1, VCENTER, 1, VMAX])
    cbar.set_ticklabels([f"{VMIN:.0f}", "−2", "−1", "0", "+1", f"+{VMAX:.0f}"])

    # ------------------------------------------------------------------
    # Domain annotation bar
    # ------------------------------------------------------------------
    ax_dom.set_xlim(0.5, NRD1_LENGTH + 0.5)
    ax_dom.set_ylim(0, 1)
    ax_dom.set_yticks([])

    for name, start, end, color in DOMAINS:
        ax_dom.axvspan(start - 0.5, end + 0.5, color=color, alpha=0.85)
        mid = (start + end) / 2
        ax_dom.text(
            mid, 0.5, name,
            ha="center", va="center",
            fontsize=9, fontweight="bold", color="white",
        )

    ax_dom.set_xlabel("Amino acid position (NRD1)", fontsize=11)

    # X-axis ticks only on domain bar
    xtick_locs = list(range(50, NRD1_LENGTH + 1, 50))
    ax_dom.set_xticks(xtick_locs)
    ax_dom.set_xticklabels([str(x) for x in xtick_locs], fontsize=8)
    ax_dom.tick_params(axis="x", length=3)

    # ------------------------------------------------------------------
    # Title and legend for hit markers
    # ------------------------------------------------------------------
    ax_miss.set_title(
        f"NRD1 fitness landscape (MAGESTIC screen, YPD+IAA, n={len(slope_df[slope_df['treatment']==TREAT_YPD_IAA]):,} variant–position pairs)",
        fontsize=11, pad=6,
    )

    # Hit tick legend
    hit_patch = mpatches.Patch(facecolor="none", edgecolor="black", linewidth=0.8,
                               label=f"Hit position (n={len(hit_positions)})")
    ax_miss.legend(
        handles=[hit_patch],
        loc="upper right", fontsize=8,
        framealpha=0.9, edgecolor="gray",
    )

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    MANUSCRIPT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PDF, dpi=300, bbox_inches="tight")
    fig.savefig(str(OUTPUT_PDF).replace(".pdf", ".png"), dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {OUTPUT_PDF}")
    print(f"Saved: {str(OUTPUT_PDF).replace('.pdf', '.png')}")

    # ------------------------------------------------------------------
    # Print hit summary
    # ------------------------------------------------------------------
    print(f"\nHit positions ({len(hit_positions)}): {hit_positions}")
    ypd_hits = hits_df.sort_values("slope_adj_YPD+IAA")
    print(ypd_hits[["amino_acid_position", "mutation_type",
                     "slope_adj_YPD+IAA", "n_barcodes", "pct_support"]].to_string(index=False))


if __name__ == "__main__":
    main()
