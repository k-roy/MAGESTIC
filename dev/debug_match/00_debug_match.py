#!/usr/bin/env python3
"""
Local M1 debug script for the bc1_donor_bc0 matching pipeline.

What it does (run end-to-end):
  1. Load the harmonized oligo table (V545/V546/V547/V429/V454/V455 + cAMP V574-V597).
  2. Build the OligoLookupTables (perfect-donor + sliding-window + guide).
  3. Load 1-4 parsed bc1_donor_bc0 TSVs from yL406 (small subset).
  4. Apply collapse + min_reads + min_pcr_replicates gates (matching what STEP 3 does in run_pipeline).
  5. Run match_dataframe and report:
     - match_method distribution
     - num_candidates distribution
     - donor length distribution
     - Sample 5 unmatched rows + sample 5 ambiguous-window rows for manual inspection.
  6. Optional experiments:
     - --pre-trim N: strip first/last N bp from observed donor before matching
     - --window-size W: use a non-default sliding window (default 60)
     - --use-guide: extract guide from donor (positions 20-39 if donor has known structure) and pass to matcher

Usage:
  cd ~/MAGESTIC_debug
  python3 scripts/00_debug_match.py                          # baseline
  python3 scripts/00_debug_match.py --pre-trim 36            # strip 36bp flanking from each end
  python3 scripts/00_debug_match.py --window-size 30         # tighter window
  python3 scripts/00_debug_match.py --pre-trim 36 --window-size 30   # combine
"""
import argparse
import sys
from pathlib import Path

# Allow the script to find the magestic package
REPO = Path(__file__).resolve().parent.parent
SRC = REPO / "src" / "MAGESTIC_src"
sys.path.insert(0, str(SRC))

import pandas as pd
import numpy as np

from magestic.pipelines.bc1_donor_bc0.core.matching import (
    build_lookup_tables,
    match_dataframe,
    match_sequence,
)
from magestic.pipelines.bc1_donor_bc0.core.merge import merge_data_sources

DATA = REPO / "data"
HARMONIZED = DATA / "20210422_and_20240411_Bloom_et_al_16_strains_QTL_harmonized_designed_variant_oligos.tsv"


def load_oligos(filepath: Path) -> pd.DataFrame:
    """Load harmonized oligos with the columns matching.py needs (+ variant).

    individual_variants is the canonical per-edit identity the matcher collapses
    by (VARIANT_COLUMN in matching.py); it MUST be loaded or the matcher falls
    back to donor-sequence identity for every oligo. donor_changes kept only for
    the legacy/diagnostic comparison.
    """
    cols = ["oligo_name", "guide", "donor", "individual_variants", "donor_changes"]
    df = pd.read_csv(filepath, sep="\t", usecols=cols, low_memory=False)
    df = df.dropna(subset=["guide", "donor"]).reset_index(drop=True)
    print(f"  Loaded harmonized table: {len(df):,} oligos (rows with guide+donor)")
    print(f"    designed donor length: median={df['donor'].str.len().median():.0f} bp")
    return df


def load_parsed_tsvs(data_dir: Path, n_samples: int) -> pd.DataFrame:
    """Concatenate 1-N parsed bc1_donor_bc0 TSVs and tag each row with sample_number."""
    tsvs = sorted(data_dir.glob("*_bc1_donor_bc0.tsv"))[:n_samples]
    if not tsvs:
        raise FileNotFoundError(f"No parsed TSVs in {data_dir}")
    dfs = []
    for i, p in enumerate(tsvs):
        d = pd.read_csv(p, sep="\t", low_memory=False)
        d["sample_number"] = p.stem.split("_sample_")[-1].split("_")[0]
        dfs.append(d)
    df = pd.concat(dfs, ignore_index=True)
    print(f"  Loaded {len(tsvs)} sample TSVs → {len(df):,} rows; cols={list(df.columns)}")
    print(f"    observed donor length: median={df['donor'].str.len().median():.0f} bp")
    return df


def pre_trim_donor(df: pd.DataFrame, n_trim: int) -> pd.DataFrame:
    """Strip first n_trim and last n_trim bp from each observed donor."""
    if n_trim <= 0:
        return df
    df = df.copy()
    df["donor"] = df["donor"].str[n_trim:-n_trim]
    print(f"  Pre-trimmed observed donor by {n_trim} bp each side → median length now {df['donor'].str.len().median():.0f} bp")
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-samples", type=int, default=1,
                        help="number of parsed TSVs to load (1-4)")
    parser.add_argument("--pre-trim", type=int, default=0,
                        help="legacy: strip first/last N bp from observed donor "
                             "before matching (disables the matcher's anchor-trim)")
    parser.add_argument("--no-anchor-trim", action="store_true",
                        help="disable the matcher's internal anchor-based trim "
                             "(reproduces the pre-fix baseline)")
    parser.add_argument("--window-size", type=int, default=60,
                        help="sliding window size (default 60)")
    parser.add_argument("--top-n", type=int, default=200_000,
                        help="cap dataframe size for speed (default 200K rows)")
    args = parser.parse_args()

    print("=" * 70)
    print(f"DEBUG RUN: n_samples={args.n_samples} pre_trim={args.pre_trim} window={args.window_size}")
    print("=" * 70)

    print("\n--- Step 1: load harmonized oligos ---")
    oligos = load_oligos(HARMONIZED)

    print("\n--- Step 2: build lookup tables ---")
    tables = build_lookup_tables(oligos, window_size=args.window_size)

    print("\n--- Step 3: load parsed TSVs ---")
    obs = load_parsed_tsvs(DATA, n_samples=args.n_samples)
    if len(obs) > args.top_n:
        obs = obs.head(args.top_n)
        print(f"  Capped to first {args.top_n:,} rows for speed")

    print("\n--- Step 4: optional fixed pre-trim ---")
    # The matcher applies an anchor-based trim internally (pre_trim=True). A
    # fixed-bp --pre-trim is a separate legacy experiment; if requested, apply
    # it here and turn OFF the anchor trim to avoid double-trimming.
    anchor_trim = not args.no_anchor_trim
    if args.pre_trim:
        obs = pre_trim_donor(obs, args.pre_trim)
        anchor_trim = False
        print("  (anchor-trim disabled because fixed --pre-trim was given)")
    else:
        print(f"  anchor-trim in matcher: {anchor_trim}")

    print("\n--- Step 5: match ---")
    matched = match_dataframe(
        obs, tables, donor_col="donor", guide_col="guide",
        prefix="", n_jobs=1, chunk_size=5000, pre_trim=anchor_trim,
    )

    print("\n--- Step 6: report ---")
    print(f"Total rows: {len(matched):,}")
    print()
    print("match_method value counts:")
    print(matched["match_method"].value_counts(dropna=False).to_string())
    print()
    print("num_candidates distribution:")
    print(matched["num_candidates"].describe().to_string())
    print()
    # DENOMINATOR FIX: a NaN/empty observed donor is an UPSTREAM PARSE FAILURE
    # (no sequence to match), not a matcher miss. Report it separately and judge
    # the matcher on the MATCHEABLE population (non-null, non-empty donor).
    n_total = len(matched)
    donor_str = matched["donor"].astype("string")
    is_parse_fail = donor_str.isna() | (donor_str.str.len() == 0)
    n_parse_fail = int(is_parse_fail.sum())
    m = matched[~is_parse_fail]
    n = len(m)
    print(f"parse failures (donor NaN/empty, EXCLUDED from matcher denom): "
          f"{n_parse_fail:,} ({n_parse_fail/n_total:.1%} of {n_total:,})")
    print(f"matcheable population (denominator): {n:,}")
    print()
    n_perfect = int(m["perfect_donor"].sum())
    n_middle = int(m["perfect_middle_region"].sum())
    n_unambig = int(m["is_unambiguous"].sum())
    n_ambig_sw = int((m["match_method"] == "sliding_window_ambiguous").sum())
    n_none = int((m["match_method"] == "none").sum())
    print(f"perfect_donor=True: {n_perfect:,} ({n_perfect/n:.1%})")
    print(f"perfect_middle_region=True (sliding-or-perfect): {n_middle:,}")
    print(f"is_unambiguous=True (VARIANT level, individual_variants): {n_unambig:,} ({n_unambig/n:.1%})")
    print(f"sliding_window_ambiguous: {n_ambig_sw:,} ({n_ambig_sw/n:.1%})")
    print(f"unmatched-with-data (none): {n_none:,} ({n_none/n:.1%})")
    print()
    print("--- SUCCESS CRITERIA (on matcheable population) ---")
    print(f"  perfect_donor >= 50%:           {n_perfect/n:.1%}  {'PASS' if n_perfect/n>=0.50 else 'MISS'}")
    print(f"  sliding_window_ambiguous <= 20%:{n_ambig_sw/n:.1%}  {'PASS' if n_ambig_sw/n<=0.20 else 'MISS'}")
    print(f"  is_unambiguous >= 80%:          {n_unambig/n:.1%}  {'PASS' if n_unambig/n>=0.80 else 'MISS'}")
    print(f"  unmatched <= 10%:               {n_none/n:.1%}  {'PASS' if n_none/n<=0.10 else 'MISS'}")
    # num_candidate_variants: how concentrated is the variant collapse?
    if "num_candidate_variants" in m.columns:
        matched_rows = m[m["match_method"] != "none"]
        print()
        print("num_candidate_variants among matched (1 = unambiguous edit):")
        print(matched_rows["num_candidate_variants"].describe().to_string())

    # --- Residual characterization: aberrant-edit vs parse/sequencing artifact ---
    # An unmatched read with the anchor present but a post-anchor donor not in
    # the design = candidate aberrant edit. No anchor present = parse/seq artifact.
    ANCHOR = "GTTTGAAGAGC"
    unmatched = matched[matched["match_method"] == "none"].copy()
    print(f"\n--- Residual (unmatched) characterization: n={len(unmatched):,} ---")
    if len(unmatched) > 0:
        # RANDOM sample, not head() -- NaN-donor parse failures cluster at the
        # file top and bias a head() sample badly (seed fixed for reproducibility).
        sample = unmatched.sample(n=min(500, len(unmatched)), random_state=0)
        has_anchor = 0; missing_donor = 0; no_anchor = 0
        for d in sample["donor"]:
            if not isinstance(d, str) or not d:
                missing_donor += 1
            elif ANCHOR in d:
                has_anchor += 1
            else:
                no_anchor += 1
        n_s = len(sample)
        print(f"  sampled {n_s} unmatched rows:")
        print(f"    anchor present but donor not in design (ABERRANT-EDIT candidate): "
              f"{has_anchor} ({has_anchor/n_s:.1%})")
        print(f"    donor present but no anchor (parse/seq artifact): "
              f"{no_anchor} ({no_anchor/n_s:.1%})")
        print(f"    donor missing/NaN (parse failure artifact): "
              f"{missing_donor} ({missing_donor/n_s:.1%})")
        print("  examples:")
        for _, row in sample.head(3).iterrows():
            dval = row["donor"]
            dstr = dval if isinstance(dval, str) else f"<{type(dval).__name__}:{dval}>"
            print(f"    bc1={row['bc1']}  donor[0:80]={dstr[:80]}  len={len(dstr) if isinstance(dval,str) else 0}")


if __name__ == "__main__":
    main()
