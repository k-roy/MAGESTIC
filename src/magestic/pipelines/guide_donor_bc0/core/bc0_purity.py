"""
Step 04 core logic: Calculate BC0 purity for guide-donor-BC0 combinations.

For each unique (guide, BC0) pair, purity is defined as the fraction of reads
that support the single most common donor sequence.  High-purity pairs represent
clean, single-donor BC0 assignments; low-purity pairs indicate contamination
from multiple donor sequences sharing the same BC0.

Author: Kevin R. Roy
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from magestic.utils.provenance import StepProvenance


INPUT_SUFFIX = "_matched_guide_donor_bc0_counts.tsv"
OUTPUT_SUFFIX = "_guide_donor_bc0_purity_filtered.tsv"


@dataclass
class PurityConfig:
    """Configuration for BC0 purity calculation."""
    min_counts: int = 1
    """Minimum read count for a row to be included in the analysis."""
    library_prefix: str = "V"
    """Library ID prefix used when searching for input files (e.g. ``"V"`` for V629)."""


def find_input_files(
    input_dir: Path,
    suffix: str = INPUT_SUFFIX,
    prefix: str = "",
) -> List[Path]:
    """
    Return sorted paths matching ``{prefix}*{suffix}`` in *input_dir*.

    Parameters
    ----------
    input_dir : Path
        Directory to search.
    suffix : str
        File suffix pattern.
    prefix : str
        Optional library-ID prefix filter.

    Returns
    -------
    list[Path]
        Sorted matching file paths.
    """
    pattern = f"{prefix}*{suffix}" if prefix else f"*{suffix}"
    return sorted(input_dir.glob(pattern))


def calculate_bc0_purity(
    input_file: Path,
    output_file: Path,
    config: Optional[PurityConfig] = None,
) -> Dict[str, object]:
    """
    Calculate purity for every guide-BC0 pair in one matched-counts file.

    For each (guide, bc0) group:

    1. Sum counts across all donor variants.
    2. Sum counts for identical observed donors, then select the dominant donor.
    3. ``purity = dominant_donor_counts / total_count``.

    ``counts`` and other sequence-row fields come from the highest-count row
    supporting the selected donor. ``dominant_donor_counts`` includes all of
    that donor's rows, including different leaders. Equal-support donor ties
    retain first appearance, as before; purity alone does not assert uniqueness.
    Missing donor reads remain in the denominator but cannot establish a donor
    sequence. Groups without an observed donor retain their representative row
    and denominator, with ``donor_observed=False`` and undefined purity.

    Parameters
    ----------
    input_file : Path
        Matched guide-donor-BC0 counts file from Step 03.
        Must have columns: ``guide``, ``bc0``, ``donor``, ``counts``.
    output_file : Path
        Destination TSV file (created or overwritten).
    config : PurityConfig, optional
        Purity configuration.  Defaults to ``PurityConfig()``.

    Returns
    -------
    dict
        Statistics: total_rows, unique_guide_bc0, mean_purity, median_purity,
        high_purity_count (>= 0.9), low_purity_count (< 0.5).
    """
    if config is None:
        config = PurityConfig()

    prov = StepProvenance("guide_donor_bc0", "04_purity", output_file,
                          description="Dominant donor support per guide/bc0, aggregated across sequence rows")
    prov.add_input(input_file, role="matched_guide_donor_bc0_counts")
    prov.set_parameters(min_counts=config.min_counts, high_purity_threshold=0.9,
                        low_purity_threshold=0.5, donor_key=["guide", "bc0", "donor"],
                        representative_row="highest count within dominant donor; first appearance on ties",
                        missing_group_keys="rows with missing guide or bc0 excluded and counted",
                        missing_donors="denominator only; groups without observed donor have undefined purity",
                        min_counts_scope="individual sequence rows before donor aggregation")
    df = pd.read_csv(input_file, sep="\t")
    required = {"guide", "bc0", "donor", "counts"}
    if required - set(df.columns):
        raise ValueError(f"Missing purity input columns: {sorted(required - set(df.columns))}")
    rows_in, reads_in = len(df), int(df["counts"].sum())

    if config.min_counts > 1:
        df = df[df["counts"] >= config.min_counts]

    filtered_rows, filtered_reads = rows_in - len(df), reads_in - int(df["counts"].sum())
    missing_keys = df[["guide", "bc0"]].isna().any(axis=1)
    missing_key_rows = int(missing_keys.sum())
    missing_key_reads = int(df.loc[missing_keys, "counts"].sum())
    eligible = df.loc[~missing_keys]
    missing_donors = eligible["donor"].isna()

    # Total counts per guide-bc0 pair
    summed = (
        eligible.groupby(["guide", "bc0"])["counts"]
        .sum()
        .reset_index()
        .rename(columns={"counts": "grouped_guide_bc0_counts"})
    )

    # Select a donor using all its supporting rows, while retaining a real
    # representative row for leader/assignment flags. Missing identities cannot
    # be combined into evidence for one identical donor sequence.
    if eligible.empty:
        top = eligible.copy()
        top["dominant_donor_counts"] = pd.Series(dtype="int64")
        top["donor_observed"] = pd.Series(dtype="bool")
    else:
        donors = eligible.loc[~missing_donors].groupby(["guide", "bc0", "donor"], sort=False)["counts"].agg(
            dominant_donor_counts="sum", representative_index="idxmax").reset_index()
        winners = donors.loc[donors.groupby(["guide", "bc0"], sort=False)["dominant_donor_counts"].idxmax()]
        top = df.loc[winners["representative_index"]].copy()
        top["dominant_donor_counts"] = winners["dominant_donor_counts"].to_numpy()
        top["donor_observed"] = True
        has_observed = eligible.groupby(["guide", "bc0"])["donor"].transform("count") > 0
        unknown = eligible.loc[~has_observed]
        if not unknown.empty:
            unassessable = unknown.loc[unknown.groupby(["guide", "bc0"], sort=False)["counts"].idxmax()].copy()
            unassessable["dominant_donor_counts"] = 0
            unassessable["donor_observed"] = False
            top = pd.concat([top, unassessable], ignore_index=True)

    merged = top.merge(summed, on=["guide", "bc0"])
    merged["guide_bc0_purity"] = merged["dominant_donor_counts"] / merged["grouped_guide_bc0_counts"]
    merged.loc[~merged["donor_observed"], "guide_bc0_purity"] = float("nan")

    merged.to_csv(output_file, sep="\t", index=False)

    stats = {
        "total_rows": len(df),
        "unique_guide_bc0": len(merged),
        "mean_purity": float(merged["guide_bc0_purity"].mean()) if merged["guide_bc0_purity"].notna().any() else None,
        "median_purity": float(merged["guide_bc0_purity"].median()) if merged["guide_bc0_purity"].notna().any() else None,
        "unassessable_groups": int((~merged["donor_observed"]).sum()),
        "high_purity_count": int((merged["guide_bc0_purity"] >= 0.9).sum()),
        "low_purity_count": int((merged["guide_bc0_purity"] < 0.5).sum()),
    }
    prov.set_counts(rows_in=rows_in, reads_in=reads_in, rows_out=len(merged),
                    eligible_group_reads=int(merged["grouped_guide_bc0_counts"].sum()),
                    filtered_rows=filtered_rows, filtered_reads=filtered_reads,
                    missing_key_rows=missing_key_rows, missing_key_reads=missing_key_reads,
                    missing_donor_rows=int(missing_donors.sum()),
                    missing_donor_reads=int(eligible.loc[missing_donors, "counts"].sum()),
                    **stats)
    prov.note("This is plasmid donor-sequence purity, not genomic editing efficiency. Representative sequence counts and aggregate dominant-donor counts are distinct.")
    prov.add_output(output_file, rows=len(merged))
    prov.write()
    return stats


def write_purity_summary(stats: List[Dict[str, object]], input_files: List[Path],
                         output_file: Path, config: PurityConfig) -> None:
    """Write per-library summary statistics and their producing-input identity."""
    prov = StepProvenance("guide_donor_bc0", "04_purity_summary", output_file)
    for path in input_files:
        prov.add_input(path, role="matched_guide_donor_bc0_counts")
    prov.set_parameters(min_counts=config.min_counts, high_purity_threshold=0.9,
                        low_purity_threshold=0.5,
                        aggregation="per library; purity uses aggregated exact donor support")
    pd.DataFrame(stats).to_csv(output_file, sep="\t", index=False)
    prov.set_counts(input_libraries=len(input_files), rows_out=len(stats),
                    input_rows_after_count_filter=sum(int(s["total_rows"]) for s in stats),
                    guide_bc0_groups=sum(int(s["unique_guide_bc0"]) for s in stats))
    prov.add_output(output_file, rows=len(stats))
    prov.write()
