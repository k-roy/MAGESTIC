"""
Step H — plot REDI-isolated clone counts + per-aa-position coverage.

Distilled from legacy `step_h_plot_REDI_isolated_clones.py`. The plotting
helpers are split out from the project-specific NRD1/NAB3/SEN1 ORF-length
constants, which now belong to per-screen drivers under
`REDI/screens/<screen>/`.

This module only writes:
    plots/bc1_clones_distribution.png
    plots/sublibrary_ID_counts_barchart.png

Per-gene amino-acid-position coverage plots are driven by the screen's
`driver.py` since they depend on the per-screen sublibrary → gene mapping
(`sublib_to_gene`) and ORF-length table.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..core.config import REDIConfig


def filter_high_quality_clones(
    annotated_cross_array_df: pd.DataFrame,
    array_a_name: str,
    array_b_name: str,
    min_counts: int = 20,
    min_purity: float = 0.95,
    min_replicates: int = 2,
) -> pd.DataFrame:
    """Filter cross-array merged dataframe to high-confidence clones.

    Mirrors legacy step_h's filter:
        ((A_counts >= 20 and A_purity >= 0.95) or
         (B_counts >= 20 and B_purity >= 0.95))
        and num_independent_replicates > 1
    """
    a_pass = (
        f"({array_a_name}_REDI_bc_bc1_counts >= {min_counts} and "
        f"{array_a_name}_colony_purity >= {min_purity})"
    )
    b_pass = (
        f"({array_b_name}_REDI_bc_bc1_counts >= {min_counts} and "
        f"{array_b_name}_colony_purity >= {min_purity})"
    )
    out = annotated_cross_array_df.query(f"{a_pass} or {b_pass}")
    if "num_independent_replicates" in out.columns:
        out = out.query(f"num_independent_replicates > {min_replicates - 1}")
    return out


def unique_bc1_by_mean_purity(
    df: pd.DataFrame,
    array_a_name: str,
    array_b_name: str,
) -> pd.DataFrame:
    """Add mean_purity, fill bc1 from B then A, dedupe by bc1."""
    df = df.copy()
    df["mean_purity"] = df[
        [f"{array_a_name}_colony_purity", f"{array_b_name}_colony_purity"]
    ].mean(axis=1)
    df["bc1"] = df[f"{array_b_name}_bc1"].fillna(df[f"{array_a_name}_bc1"])
    return df.sort_values(by="mean_purity", ascending=False).drop_duplicates(
        subset=["bc1"]
    )


def run_step_h(
    config: REDIConfig,
    annotated_cross_array_path: Path | None = None,
) -> dict[str, Path]:
    """Run step_h. Returns dict of plot/output paths."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "step_h requires matplotlib + seaborn; install magestic[full] "
            "or `pip install matplotlib seaborn`"
        ) from e

    combined_data = config.paths.combined_data
    plots = config.paths.plots
    plots.mkdir(parents=True, exist_ok=True)

    if not config.array_pairs:
        raise ValueError("step_h needs at least one configured array pair")
    pair = config.array_pairs[0]

    if annotated_cross_array_path is None:
        annotated_cross_array_path = (
            combined_data / f"{pair.name_a}_{pair.name_b}_merged_REDI_df.tsv"
        )

    df = pd.read_csv(annotated_cross_array_path, sep="\t")

    # Join with bc1 reference if available
    if config.bc1_reference_table is not None and Path(
        config.bc1_reference_table
    ).exists():
        bc1_ref = pd.read_csv(config.bc1_reference_table, sep="\t")
        df = unique_bc1_by_mean_purity(df, pair.name_a, pair.name_b)
        df = df.merge(bc1_ref, on="bc1", how="left")

    filtered = filter_high_quality_clones(df, pair.name_a, pair.name_b)

    clones_per_bc1 = filtered["bc1"].value_counts().reset_index()
    clones_per_bc1.columns = ["bc1", "num_clones"]

    dist = clones_per_bc1["num_clones"].value_counts().sort_index()
    plt.figure(figsize=(10, 6))
    sns.barplot(x=dist.index, y=dist.values)
    plt.xlabel("Number of Clones per bc1")
    plt.ylabel("Count")
    plt.title("Distribution of Number of Clones per bc1")
    plt.xticks(rotation=45)
    plt.tight_layout()
    dist_path = plots / "bc1_clones_distribution.png"
    plt.savefig(dist_path)
    plt.close()

    # Sublibrary counts barchart (only if sublibrary_ID column is present)
    out: dict[str, Path] = {"bc1_clones_distribution": dist_path}
    if "sublibrary_ID" in filtered.columns:
        sublib_path = plots / "sublibrary_ID_counts_barchart.png"
        _plot_sublibrary_counts(df, filtered, sublib_path)
        out["sublibrary_ID_counts_barchart"] = sublib_path
    return out


def _plot_sublibrary_counts(
    all_df: pd.DataFrame,
    filtered_df: pd.DataFrame,
    out_path: Path,
) -> None:
    import matplotlib.pyplot as plt

    counts_df = pd.DataFrame(
        {
            "all_picked_colonies": all_df["sublibrary_ID"].value_counts().sort_index(),
            "all_picked_colonies_passing_filter": filtered_df[
                "sublibrary_ID"
            ]
            .value_counts()
            .sort_index(),
        }
    ).fillna(0).astype(int)
    colors = ["#1f77b4", "#ff7f0e"]
    counts_df.plot(kind="bar", color=colors[: counts_df.shape[1]])
    plt.xlabel("sublibrary_ID")
    plt.ylabel("Count")
    plt.title("Counts by sublibrary_ID")
    plt.tight_layout()
    plt.legend(title="Source", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.savefig(out_path, bbox_inches="tight")
    plt.close()
