#!/usr/bin/env python3
"""
Step 10: Split Count Tables by Comparison

Splits the wide-format count matrix into per-comparison reference/treatment
pairs for DESeq2 analysis.

IMPORTANT: Uses screen_date-based t0 consolidation.
All gen=0 samples from the same screen_date are from the same glycerol aliquot
(technical replicates), so they are consolidated as reference for all
comparisons on that screen_date.

Creates:
- comparison_manifest.tsv: List of all comparisons with metadata
- reference/{comparison}_reference.tsv: Reference (t0) counts
- treatment/{comparison}_treatment.tsv: Treatment counts

Usage:
    python 10_split_comparisons.py \
        --count-matrix count_matrix.tsv \
        --sample-key sample_key.tsv \
        --library yL406 \
        --output-dir counts_for_deseq/yL406

Author: Kevin R. Roy
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

# Shared rules (sibling modules free of magestic imports, like config.py), so steps
# 08, 10 and 11 cannot drift apart.
try:
    from magestic.pipelines.bc1_from_screen.library_scope import LIBRARY_COL_CANDIDATES
    from magestic.pipelines.bc1_from_screen.transformation_scope import (
        SAMPLE_TRANSFORMATION_COLUMNS, is_design_pool_column, parse_generation)
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from library_scope import LIBRARY_COL_CANDIDATES
    from transformation_scope import (
        SAMPLE_TRANSFORMATION_COLUMNS, is_design_pool_column, parse_generation)

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Split count matrix into per-comparison tables"
    )
    parser.add_argument(
        "--count-matrix", type=Path, required=True,
        help="Input wide-format count matrix TSV"
    )
    parser.add_argument(
        "--sample-key", type=Path, required=True,
        help="Sample key TSV with condition/generation/screen_date info"
    )
    parser.add_argument(
        "--library", type=str, required=True,
        help="Library name"
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True,
        help="Output directory for split tables"
    )
    parser.add_argument(
        "--screen-date-column", type=str, default="screen_date",
        help="Column name for screen date (default: screen_date)"
    )
    parser.add_argument(
        "--condition-column", type=str, default="condition",
        help="Column name for condition (default: condition)"
    )
    parser.add_argument(
        "--generation-column", type=str, default="generations",
        help="Column name for generation (default: generations)"
    )
    parser.add_argument(
        "--sample-name-column", type=str, default="sample_name",
        help="Column name for sample name (default: sample_name)"
    )
    parser.add_argument(
        "--consolidate-by-screen-date", action="store_true", default=True,
        help="Consolidate all gen=0 samples from same screen_date as reference (default: True)"
    )
    parser.add_argument(
        "--no-consolidate", dest="consolidate_by_screen_date", action="store_false",
        help="Pool every generation-0 column into one reference and pool treatments across "
             "screen dates (explicit pooled-dates estimand)"
    )
    # Exactly one transformation scope must be declared; pooling is not a default.
    scope = parser.add_mutually_exclusive_group()
    scope.add_argument(
        "--split-by-yeast-sublibrary", action="store_true", default=False,
        help="Split comparisons by a construction-derived bc1 transformation column "
             "(requires --bc1-reference-table and --yeast-sublibrary-column)"
    )
    scope.add_argument(
        "--pool-transformations", action="store_true", default=False,
        help="Analyse every barcode together across transformations (explicit pooling)"
    )
    parser.add_argument(
        "--bc1-reference-table", type=Path,
        help="BC1 reference table carrying the transformation column (required if --split-by-yeast-sublibrary)"
    )
    parser.add_argument(
        "--yeast-sublibrary-column", type=str, default=None,
        help="Construction-derived bc1 -> transformation column in the reference table, e.g. trafo_ID"
    )
    parser.add_argument(
        "--exclude-samples-file", type=Path, default=None,
        help="Count-matrix sample columns to drop before analysis, one per line (# comments allowed)"
    )
    return parser.parse_args()


def _declared(value):
    """A key cell as stripped text, or None when blank or missing."""
    if value is None or pd.isna(value) or str(value).strip() == "":
        return None
    return value


def build_sample_info(sample_key: pd.DataFrame,
                      sample_columns: List[str],
                      sample_name_col: str,
                      condition_col: str,
                      generation_col: str,
                      screen_date_col: str,
                      *,
                      require_screen_date: bool = True,
                      scope_columns=()) -> Dict:
    """
    Build sample info dictionary from sample key.

    Refuses (ValueError naming the columns or samples) instead of defaulting: a missing
    required column; a count column with no key row; key rows for one sample that
    disagree on condition, generation, date or scope columns (identical rows are
    accepted); a missing or non-integer generation ("6.0" is 6; "none", text and 0.5
    are refused); a missing condition for a non-baseline sample; and, when
    ``require_screen_date``, a missing screen date.

    Returns:
        {sample_name: {
            "condition": str,
            "generation": int,
            "screen_date": str
        }}
    """
    gen_col = generation_col
    if gen_col not in sample_key.columns and "generation" in sample_key.columns:
        gen_col = "generation"
    required = [sample_name_col, condition_col, gen_col]
    if require_screen_date:
        required.append(screen_date_col)
    missing = [c for c in required if c not in sample_key.columns]
    if missing:
        raise ValueError(
            f"Sample key lacks required column(s): {missing}. Available: {list(sample_key.columns)}. "
            f"Map them with --condition-column / --generation-column / --screen-date-column, or pass "
            f"--no-consolidate to pool every date explicitly (pooled-dates estimand).")

    wanted = {str(c) for c in sample_columns}
    names = sample_key[sample_name_col]
    rows = sample_key[names.notna() & names.astype(str).isin(wanted)]
    uncovered = sorted(wanted - set(rows[sample_name_col].astype(str)))
    if uncovered:
        raise ValueError(
            f"Count-matrix column(s) have no row in the sample key: {uncovered[:20]}"
            f"{' ...' if len(uncovered) > 20 else ''}. Step 08 consolidated or pseudo-replicate "
            f"columns (t0_<date>_consolidated, t0_pseudorep_<i>) need their metadata; "
            f"pseudo-replicates are technical splits of one aliquot.")

    scope_columns = [c for c in scope_columns if c in sample_key.columns]
    has_date = screen_date_col in sample_key.columns
    parsed: Dict[str, set] = {}
    bad_generation = []
    for _, row in rows.iterrows():
        sample = str(row[sample_name_col])
        try:
            generation = parse_generation(row[gen_col])
        except ValueError:
            bad_generation.append(f"{sample}={row[gen_col]!r}")
            continue
        condition = _declared(row[condition_col])
        date = _declared(row[screen_date_col]) if has_date else None
        scope = tuple(None if _declared(row[c]) is None else str(row[c]).strip() for c in scope_columns)
        parsed.setdefault(sample, set()).add(
            (None if condition is None else condition, generation,
             None if date is None else str(date), scope))
    if bad_generation:
        raise ValueError(
            f"Generation must be an integer for every analysis sample; refused: "
            f"{', '.join(sorted(bad_generation)[:20])} (\"6.0\" is accepted as 6; missing, \"none\", "
            f"text and fractional values are not).")

    conflicts = {s: sorted(map(str, v)) for s, v in parsed.items() if len(v) > 1}
    if conflicts:
        shown = "; ".join(f"{s}: {' vs '.join(v)}" for s, v in sorted(conflicts.items())[:10])
        raise ValueError(f"Conflicting key rows for sample(s): {shown}. Resolve the key; step 10 will not pick one.")

    sample_info = {}
    no_condition, no_date = [], []
    for col in sample_columns:
        condition, generation, date, _ = next(iter(parsed[str(col)]))
        if condition is None and generation != 0:
            no_condition.append(str(col))
        if require_screen_date and date is None:
            no_date.append(str(col))
        sample_info[col] = {
            "condition": condition if condition is not None else "unknown",
            "generation": generation,
            "screen_date": date if date is not None else "unknown",
        }
    if no_condition:
        raise ValueError(f"Column '{condition_col}' is empty for non-baseline sample(s): {sorted(no_condition)[:20]}.")
    if no_date:
        raise ValueError(
            f"Column '{screen_date_col}' is empty for sample(s): {sorted(no_date)[:20]}. "
            f"Pass --no-consolidate to pool every date explicitly.")
    return sample_info


def extract_yeast_sublibrary_mapping(bc1_reference: pd.DataFrame, yeast_sublibrary_col: str) -> Dict[str, str]:
    """
    bc1 -> transformation from a construction-derived column of the reference table.

    Raises if the column is absent or names a design pool (oligo_name, SPS, V/cV
    library IDs...), or if a bc1 maps to more than one value; step 10 never chooses.
    Values are used verbatim after stripping. A bc1 with no value belongs to no
    stratum and is absent from the mapping.

    Returns:
        {bc1: transformation}
    """
    if yeast_sublibrary_col not in bc1_reference.columns:
        raise ValueError(
            f"--yeast-sublibrary-column '{yeast_sublibrary_col}' is not a column of the bc1 "
            f"reference table; available: {list(bc1_reference.columns)}.")
    if is_design_pool_column(yeast_sublibrary_col):
        raise ValueError(
            f"--yeast-sublibrary-column '{yeast_sublibrary_col}' identifies a design pool, not a "
            f"transformation; name a construction-derived column such as trafo_ID.")

    df = bc1_reference[["bc1", yeast_sublibrary_col]].dropna(subset=["bc1"]).copy()
    df["value"] = df[yeast_sublibrary_col].map(
        lambda v: None if _declared(v) is None else str(v).strip())
    df = df.dropna(subset=["value"])
    n_values = df.groupby("bc1")["value"].nunique()
    multi = list(n_values[n_values > 1].index)
    if multi:
        shown = ", ".join(
            f"{b}=" + "|".join(sorted(df.loc[df["bc1"] == b, "value"].unique())) for b in multi[:5])
        raise ValueError(
            f"{len(multi)} bc1(s) map to more than one transformation in '{yeast_sublibrary_col}': "
            f"{shown}. Resolve the reference table; step 10 will not choose one.")
    return dict(zip(df["bc1"], df["value"]))


def check_sample_scope(sample_key: pd.DataFrame, analysis_samples, sample_name_col: str,
                       pool_transformations: bool) -> None:
    """Refuse analysis samples that span transformations (unless pooling) or library captures.

    A row split cannot separate columns, so transformation labels on the sample axis can
    only refuse. Library captures are checked for uniqueness only; --library is an
    output label here, not a filter (step 08 scopes the matrix).
    """
    wanted = {str(s) for s in analysis_samples}
    rows = sample_key[sample_key[sample_name_col].astype(str).isin(wanted)]

    def per_value(column, missing_is_value):
        counts: Dict[str, set] = {}
        for sample, value in zip(rows[sample_name_col].astype(str), rows[column]):
            label = _declared(value)
            if label is None:
                if not missing_is_value:
                    continue
                label = "<missing>"
            counts.setdefault(str(label).strip(), set()).add(sample)
        return {k: len(v) for k, v in sorted(counts.items())}

    if not pool_transformations:
        for column in SAMPLE_TRANSFORMATION_COLUMNS:
            if column in rows.columns:
                values = per_value(column, missing_is_value=False)
                if len(values) > 1:
                    shown = ", ".join(f"{v} ({n} samples)" for v, n in values.items())
                    raise ValueError(
                        f"Analysis samples span {len(values)} transformations in key column '{column}': "
                        f"{shown}. A row split cannot separate columns; remove columns with "
                        f"--exclude-samples-file, or pass --pool-transformations to pool them explicitly.")
    for column in LIBRARY_COL_CANDIDATES:
        if column in rows.columns:
            values = per_value(column, missing_is_value=True)
            if len(values) > 1:
                shown = ", ".join(f"{v} ({n})" for v, n in values.items())
                raise ValueError(
                    f"Analysis samples span {len(values)} library captures in key column '{column}': "
                    f"{shown}. Step 10 compares within one capture; scope the matrix with step 08 "
                    f"--library or --exclude-samples-file.")


def check_baselines(sample_info: Dict, consolidate: bool) -> None:
    """Refuse treatment groups without a generation-0 baseline in scope, and runs with no treatment.

    Groups exactly as the comparison builders below do; nothing is skipped silently.
    """
    missing = []
    if consolidate:
        by_date = group_samples_by_screen_date(sample_info)
        n_treatment = sum(len(v) for data in by_date.values() for v in data["treat"].values())
        for date, data in sorted(by_date.items(), key=lambda item: str(item[0])):
            if not data["ref"]:
                missing += [f"{date}/{key} ({len(v)} samples)" for key, v in sorted(data["treat"].items())]
    else:
        references = [s for s, info in sample_info.items() if info["generation"] == 0]
        groups: Dict[str, list] = {}
        for sample, info in sample_info.items():
            if info["generation"] > 0:
                groups.setdefault(f"{info['condition']}_{info['generation']}", []).append(sample)
        n_treatment = sum(len(v) for v in groups.values())
        if not references:
            missing += [f"{key} ({len(v)} samples)" for key, v in sorted(groups.items())]
    if n_treatment == 0:
        raise ValueError("No treatment samples (generation other than 0) among the analysis samples; "
                         "nothing to compare.")
    if missing:
        raise ValueError(f"Treatment group(s) without a generation-0 baseline in scope: "
                         f"{', '.join(missing)}. Refusing instead of skipping them.")


def group_samples_by_screen_date(sample_info: Dict) -> Dict[str, Dict]:
    """
    Group samples by screen_date for t0 consolidation.

    All gen=0 samples from the same screen_date are from the same
    glycerol aliquot (technical replicates), so we consolidate them.

    Returns:
        {screen_date: {
            "ref": [list of gen=0 sample names],
            "treat": {comparison_id: [list of treatment sample names]}
        }}
    """
    by_screen_date = {}

    for sample, info in sample_info.items():
        screen_date = info["screen_date"]

        if screen_date not in by_screen_date:
            by_screen_date[screen_date] = {"ref": [], "treat": {}}

        if info["generation"] == 0:
            by_screen_date[screen_date]["ref"].append(sample)
        else:
            condition = info["condition"]
            gen = info["generation"]
            key = f"{condition}_{gen}"

            if key not in by_screen_date[screen_date]["treat"]:
                by_screen_date[screen_date]["treat"][key] = []
            by_screen_date[screen_date]["treat"][key].append(sample)

    return by_screen_date


def create_comparisons_by_screen_date(
    count_matrix_filtered: pd.DataFrame,
    sample_info: Dict,
    args,
    reference_dir: Path,
    treatment_dir: Path,
    sublib_suffix: str,
    yeast_sublib
) -> List[Dict]:
    """
    Create comparisons using screen_date-based t0 consolidation.

    Returns list of manifest rows.
    """
    manifest_rows = []
    bc1_col = "bc1"

    by_screen_date = group_samples_by_screen_date(sample_info)

    logger.info("Screen-date summary:")
    for sd, data in sorted(by_screen_date.items()):
        n_ref = len(data["ref"])
        n_treat_groups = len(data["treat"])
        logger.info(f"  {sd}: {n_ref} reference samples, {n_treat_groups} treatment groups")

    # Create comparisons grouped by screen_date
    for screen_date, data in by_screen_date.items():
        reference_samples = data["ref"]

        if not reference_samples:
            logger.warning(f"  No reference (gen=0) samples for screen_date {screen_date}, skipping")
            continue

        for comparison_id, treatment_samples in data["treat"].items():
            if not treatment_samples:
                continue

            # Sanitize for filename - include screen_date to make unique
            safe_id = comparison_id.replace("/", "_").replace(" ", "_").replace("(", "").replace(")", "")
            full_id = f"{screen_date}_{safe_id}"  # Include screen_date for uniqueness

            logger.info(f"  Processing: {comparison_id} (screen_date={screen_date})")
            logger.info(f"    Reference: {len(reference_samples)} samples (all gen=0 from {screen_date})")
            logger.info(f"    Treatment: {len(treatment_samples)} samples")

            # Reference counts - ALL gen=0 from this screen_date
            ref_cols = [bc1_col] + reference_samples
            ref_df = count_matrix_filtered[ref_cols].copy()
            ref_df = ref_df.set_index(bc1_col)

            ref_file = reference_dir / f"{args.library}{sublib_suffix}_{full_id}_reference.tsv"
            ref_df.to_csv(ref_file, sep="\t")

            # Treatment counts
            treat_cols = [bc1_col] + treatment_samples
            treat_df = count_matrix_filtered[treat_cols].copy()
            treat_df = treat_df.set_index(bc1_col)

            treat_file = treatment_dir / f"{args.library}{sublib_suffix}_{full_id}_treatment.tsv"
            treat_df.to_csv(treat_file, sep="\t")

            # Parse condition and generation
            parts = comparison_id.rsplit("_", 1)
            condition = parts[0] if len(parts) > 1 else comparison_id
            generation = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 0

            manifest_rows.append({
                "library": args.library,
                "yeast_sublibrary": yeast_sublib if yeast_sublib else "all",
                "comparison_id": f"{args.library}{sublib_suffix}_{full_id}",
                "condition": condition,
                "generation": generation,
                "screen_date": screen_date,
                "reference_file": ref_file.name,
                "treatment_file": treat_file.name,
                "n_reference_samples": len(reference_samples),
                "n_treatment_samples": len(treatment_samples),
                "n_bc1s": len(ref_df),
            })

    return manifest_rows


def create_comparisons_no_consolidation(
    count_matrix_filtered: pd.DataFrame,
    sample_info: Dict,
    args,
    reference_dir: Path,
    treatment_dir: Path,
    sublib_suffix: str,
    yeast_sublib
) -> List[Dict]:
    """
    Create comparisons using ALL gen=0 samples as reference (no screen_date grouping).

    Returns list of manifest rows.
    """
    manifest_rows = []
    bc1_col = "bc1"

    logger.info("Using all gen=0 samples as pooled reference (no screen_date consolidation)")

    # Identify all reference and treatment samples
    reference_samples = [s for s, info in sample_info.items() if info["generation"] == 0]
    treatment_samples = [s for s, info in sample_info.items() if info["generation"] > 0]

    logger.info(f"  Reference samples (all gen=0): {len(reference_samples)}")
    logger.info(f"  Treatment samples: {len(treatment_samples)}")

    if not reference_samples:
        logger.warning("  No reference (gen=0) samples found!")
        return manifest_rows

    # Group treatment by condition and generation
    comparisons = {}
    for sample in treatment_samples:
        info = sample_info[sample]
        key = f"{info['condition']}_{info['generation']}"
        if key not in comparisons:
            comparisons[key] = []
        comparisons[key].append(sample)

    logger.info(f"  Creating {len(comparisons)} comparisons")

    for comparison_id, treatment_list in comparisons.items():
        safe_id = comparison_id.replace("/", "_").replace(" ", "_").replace("(", "").replace(")", "")

        # Reference counts (ALL gen=0)
        ref_cols = [bc1_col] + reference_samples
        ref_df = count_matrix_filtered[ref_cols].copy()
        ref_df = ref_df.set_index(bc1_col)

        ref_file = reference_dir / f"{args.library}{sublib_suffix}_{safe_id}_reference.tsv"
        ref_df.to_csv(ref_file, sep="\t")

        # Treatment counts
        treat_cols = [bc1_col] + treatment_list
        treat_df = count_matrix_filtered[treat_cols].copy()
        treat_df = treat_df.set_index(bc1_col)

        treat_file = treatment_dir / f"{args.library}{sublib_suffix}_{safe_id}_treatment.tsv"
        treat_df.to_csv(treat_file, sep="\t")

        parts = comparison_id.rsplit("_", 1)
        condition = parts[0] if len(parts) > 1 else comparison_id
        generation = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 0

        manifest_rows.append({
            "library": args.library,
            "yeast_sublibrary": yeast_sublib if yeast_sublib else "all",
            "comparison_id": f"{args.library}{sublib_suffix}_{safe_id}",
            "condition": condition,
            "generation": generation,
            "screen_date": "pooled",
            "reference_file": ref_file.name,
            "treatment_file": treat_file.name,
            "n_reference_samples": len(reference_samples),
            "n_treatment_samples": len(treatment_list),
            "n_bc1s": len(ref_df),
        })

        logger.info(f"  {comparison_id}: ref={len(reference_samples)}, treat={len(treatment_list)}")

    return manifest_rows


def main():
    args = parse_args()

    logger.info(f"Step 10: Splitting comparisons for library {args.library}")
    logger.info(f"  Count matrix: {args.count_matrix}")
    logger.info(f"  Sample key: {args.sample_key}")
    logger.info(f"  Output directory: {args.output_dir}")
    logger.info(f"  Screen-date consolidation: {args.consolidate_by_screen_date}")
    logger.info(f"  Screen-date column: {args.screen_date_column}")

    # Load data
    logger.info("Loading count matrix...")
    count_matrix = pd.read_csv(args.count_matrix, sep="\t")
    logger.info(f"  Loaded {len(count_matrix):,} BC1s x {len(count_matrix.columns)} columns")

    logger.info("Loading sample key...")
    sample_key = pd.read_csv(args.sample_key, sep="\t")
    logger.info(f"  Loaded {len(sample_key):,} sample entries")
    logger.info(f"  Sample key columns: {list(sample_key.columns)}")

    # Identify sample columns (excluding metadata)
    meta_cols = ["bc1", "tier", "oligo_name"]
    sample_columns = [c for c in count_matrix.columns if c not in meta_cols]

    # Every check runs before any file is written; a refusal leaves earlier outputs untouched.
    yeast_sublibrary_bc1_mapping = None
    yeast_sublibraries_to_process = [None]  # pooled: all BC1s together
    try:
        if not (args.split_by_yeast_sublibrary or args.pool_transformations):
            raise ValueError(
                "Step 10 needs an explicit transformation scope: pass --split-by-yeast-sublibrary "
                "--bc1-reference-table <table> --yeast-sublibrary-column <construction-derived column, "
                "e.g. trafo_ID>, or --pool-transformations to analyse every barcode together. "
                "Pooling is no longer the default.")

        if args.exclude_samples_file:
            listed = [line.split("#", 1)[0].strip()
                      for line in Path(args.exclude_samples_file).read_text().splitlines()]
            excluded = {name for name in listed if name}
            unknown = sorted(excluded - {str(c) for c in sample_columns})
            if unknown:
                raise ValueError(f"--exclude-samples-file lists name(s) that are not count-matrix columns: {unknown}.")
            sample_columns = [c for c in sample_columns if str(c) not in excluded]
            count_matrix = count_matrix.drop(columns=[c for c in count_matrix.columns if str(c) in excluded])
            logger.info(f"  Excluded {len(excluded)} sample column(s) listed in {args.exclude_samples_file}")
        logger.info(f"  Sample columns: {len(sample_columns)}")

        scope_columns = [c for c in (*SAMPLE_TRANSFORMATION_COLUMNS, *LIBRARY_COL_CANDIDATES)
                         if c in sample_key.columns]
        sample_info = build_sample_info(
            sample_key,
            sample_columns,
            args.sample_name_column,
            args.condition_column,
            args.generation_column,
            args.screen_date_column,
            require_screen_date=args.consolidate_by_screen_date,
            scope_columns=scope_columns,
        )
        check_sample_scope(sample_key, sample_columns, args.sample_name_column, args.pool_transformations)

        if args.split_by_yeast_sublibrary:
            if not args.bc1_reference_table:
                raise ValueError("--split-by-yeast-sublibrary needs --bc1-reference-table.")
            if not args.yeast_sublibrary_column:
                raise ValueError(
                    "--split-by-yeast-sublibrary needs --yeast-sublibrary-column naming a construction-derived "
                    "transformation column (e.g. trafo_ID). The former default 'oligo_name' is not accepted: "
                    "oligo-name prefixes, SPS and V/cV library IDs identify design pools, and trafo_3 and "
                    "trafo_4 share V634.")
            logger.info("Loading BC1 reference table for transformation annotations...")
            bc1_reference = pd.read_csv(args.bc1_reference_table, sep="\t")
            logger.info(f"  Loaded {len(bc1_reference):,} BC1s")
            yeast_sublibrary_bc1_mapping = extract_yeast_sublibrary_mapping(
                bc1_reference, args.yeast_sublibrary_column)
            yeast_sublibraries_to_process = sorted(set(yeast_sublibrary_bc1_mapping.values()))
            logger.info(f"  Transformations to process: {yeast_sublibraries_to_process}")
            from collections import Counter
            for sublib, n in sorted(Counter(yeast_sublibrary_bc1_mapping.values()).items()):
                logger.info(f"    {sublib}: {n:,} BC1s")
            unassigned = ~count_matrix["bc1"].isin(set(yeast_sublibrary_bc1_mapping))
            if unassigned.any():
                reads = int(count_matrix.loc[unassigned, sample_columns].to_numpy().sum())
                logger.warning(
                    f"  {int(unassigned.sum()):,} count-matrix BC1(s) ({reads:,} reads) have no value in "
                    f"'{args.yeast_sublibrary_column}' and are excluded from every stratum")

        check_baselines(sample_info, args.consolidate_by_screen_date)
    except ValueError as exc:
        logger.error(str(exc))
        return 1

    screen_dates = set(info["screen_date"] for info in sample_info.values())
    logger.info(f"  Unique screen_dates: {screen_dates}")
    scope_label = ("pooled" if args.pool_transformations
                   else f"split:{args.yeast_sublibrary_column}")

    # Create output directories
    reference_dir = args.output_dir / "reference"
    treatment_dir = args.output_dir / "treatment"
    reference_dir.mkdir(parents=True, exist_ok=True)
    treatment_dir.mkdir(parents=True, exist_ok=True)

    # Manifest rows
    manifest_rows = []
    bc1_col = "bc1"

    # Loop over sublibraries (or process all BC1s together if not splitting)
    for yeast_sublib in yeast_sublibraries_to_process:
        # Filter count matrix by yeast_sublibrary if splitting
        if yeast_sublib is None:
            # Process all BC1s together
            count_matrix_filtered = count_matrix.copy()
            sublib_suffix = ""
            logger.info(f"\nProcessing all BC1s together...")
        else:
            # Filter to BC1s from this yeast_sublibrary
            bc1s_in_yeast_sublibrary = [bc1 for bc1, sl in yeast_sublibrary_bc1_mapping.items() if sl == yeast_sublib]
            count_matrix_filtered = count_matrix[count_matrix[bc1_col].isin(bc1s_in_yeast_sublibrary)].copy()
            sublib_suffix = f"_{yeast_sublib}"
            logger.info(f"\nProcessing yeast_sublibrary: {yeast_sublib}")
            logger.info(f"  Filtered to {len(count_matrix_filtered):,} BC1s")

            if len(count_matrix_filtered) == 0:
                logger.warning(f"  No BC1s found for yeast_sublibrary {yeast_sublib}, skipping")
                continue

        # Create comparisons using appropriate method
        if args.consolidate_by_screen_date:
            rows = create_comparisons_by_screen_date(
                count_matrix_filtered, sample_info, args,
                reference_dir, treatment_dir, sublib_suffix, yeast_sublib
            )
        else:
            rows = create_comparisons_no_consolidation(
                count_matrix_filtered, sample_info, args,
                reference_dir, treatment_dir, sublib_suffix, yeast_sublib
            )

        for row in rows:
            row["transformation_scope"] = scope_label  # appended after the 11 positional columns
        manifest_rows.extend(rows)

    if not manifest_rows:
        logger.error("No comparisons were produced; refusing to write an empty manifest")
        return 1

    # Save manifest
    manifest_df = pd.DataFrame(manifest_rows)
    manifest_file = args.output_dir / "comparison_manifest.tsv"
    manifest_df.to_csv(manifest_file, sep="\t", index=False)
    logger.info(f"\nSaved manifest with {len(manifest_df)} comparisons to {manifest_file}")

    # Summary
    if len(manifest_df) > 0:
        ref_counts = manifest_df["n_reference_samples"]
        logger.info(f"Reference sample counts: min={ref_counts.min()}, max={ref_counts.max()}, "
                   f"median={ref_counts.median():.0f}")

    logger.info("Step 10: Split complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
