#!/usr/bin/env python3
"""
Map bc1-donor-bc0 data to designed oligos via the CANONICAL rich matcher
(magestic.pipelines.bc1_donor_bc0.core.matching.match_dataframe).

REWIRED 2026-06-02 (user decision: core/matching is the canonical matcher).
Previously this step did a simple exact pd.merge of a fixed-length RE-site-trimmed
donor against the raw Twist oligo pools. It now delegates to the rich, vectorized
match_dataframe (variant-level unambiguity + tolerant Hamming anchoring + dual-anchor
SpCas9/SpG + LbCas12a trim + middle-region sliding-window rescue), matched against the
harmonized oligo table. ~12x faster than the old per-row matcher AND semantically
richer than the old exact merge (which lacked all of the above).

Output preserves the prior contract columns consumed by assign_05:
    oligo_name, oligo_guide, oligo_pool, perfect_donor_match
plus the rich extras: variant, is_unambiguous, perfect_middle_region, match_method,
num_candidates, num_candidate_variants, guide_match.

Usage:
    python 04_map_to_oligos.py --input <linked_tsv> [...] --output <out_tsv> \\
        [--donor-col donor] [--v-libraries V629 V631 ...] [--nuclease SpG_NGNG ...]

Author: Kevin R. Roy
"""
import argparse
import sys
import pandas as pd
from pathlib import Path

# Canonical package matcher. Works whether magestic is installed (workflow PYTHON)
# or this file is run directly by path (fallback adds .../src to sys.path).
try:
    from magestic.pipelines.bc1_donor_bc0.core.matching import load_and_build_lookup, match_dataframe
    from magestic.utils import path_utils
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # .../src
    from magestic.pipelines.bc1_donor_bc0.core.matching import load_and_build_lookup, match_dataframe
    from magestic.utils import path_utils


def _build_oligo_attr_maps():
    """oligo_name -> guide, oligo_name -> library, from the harmonized oligo table,
    for the oligo_guide / oligo_pool contract columns consumed downstream."""
    hp = path_utils.HARMONIZED_VARIANTS_FILE
    cols = pd.read_csv(hp, sep="\t", nrows=0).columns
    guide_col = "guide" if "guide" in cols else None
    lib_col = next((c for c in ("guide_donor_bc0_plasmid_library", "V_library", "twist_pool") if c in cols), None)
    use = ["oligo_name"] + [c for c in (guide_col, lib_col) if c]
    h = pd.read_csv(hp, sep="\t", usecols=use).drop_duplicates(subset=["oligo_name"])
    o2g = dict(zip(h["oligo_name"], h[guide_col])) if guide_col else {}
    o2l = dict(zip(h["oligo_name"], h[lib_col])) if lib_col else {}
    return o2g, o2l


def main():
    parser = argparse.ArgumentParser(description="Map bc1-donor-bc0 to designed oligos (canonical rich matcher)")
    parser.add_argument("--input", nargs="+", required=True, help="Input TSV file(s) from Step 03 / combine")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--donor-col", default=None, help="Donor column (default: auto donor/top_donor)")
    parser.add_argument("--guide-col", default="guide", help="Guide column if present")
    parser.add_argument("--v-libraries", nargs="*", default=None, help="Restrict harmonized lookup to these V libraries")
    parser.add_argument("--nuclease", nargs="*", default=None, help="Restrict harmonized lookup to these nuclease/PAM types")
    parser.add_argument("--no-pre-trim", action="store_true", help="Disable anchor pre-trim (already-trimmed donors)")
    # Back-compat (accepted + IGNORED): the old exact-merge args, so existing rule
    # shells keep working without a rule edit. The harmonized table is used instead.
    parser.add_argument("--oligo-pools", nargs="*", default=None, help="(ignored; harmonized table used)")
    parser.add_argument("--guide-start", type=int, default=None, help="(ignored)")
    parser.add_argument("--guide-end", type=int, default=None, help="(ignored)")
    parser.add_argument("--donor-start", type=int, default=None, help="(ignored)")
    parser.add_argument("--donor-end", type=int, default=None, help="(ignored)")
    parser.add_argument("--lbcas12a-donor-start", type=int, default=None, help="(ignored)")
    args = parser.parse_args()

    print("Loading bc1-donor-bc0 data...")
    bc1_data = pd.concat([pd.read_csv(f, sep="\t") for f in args.input], ignore_index=True)
    nbc1 = bc1_data["bc1"].nunique() if "bc1" in bc1_data.columns else "NA"
    print(f"Loaded {len(bc1_data)} rows; unique bc1s: {nbc1}")

    donor_col = args.donor_col
    if donor_col is None:
        donor_col = "donor" if "donor" in bc1_data.columns else ("top_donor" if "top_donor" in bc1_data.columns else None)
    if donor_col is None:
        print("Warning: no donor column found; writing input unchanged.")
        bc1_data.to_csv(args.output, sep="\t", index=False)
        return

    print("Building harmonized oligo lookup (canonical matcher)...")
    tables = load_and_build_lookup(
        use_harmonized_table=True,
        v_libraries=args.v_libraries,
        nuclease_types=args.nuclease,
    )

    print(f"Matching donors (col '{donor_col}') via core.matching.match_dataframe...")
    out = match_dataframe(bc1_data, tables, donor_col=donor_col, guide_col=args.guide_col,
                          pre_trim=not args.no_pre_trim)

    # Preserve the prior output contract (assign_05 consumes these).
    out["perfect_donor_match"] = out["perfect_donor"]
    o2g, o2l = _build_oligo_attr_maps()
    if o2g:
        out["oligo_guide"] = out["oligo_name"].map(o2g)
    if o2l:
        out["oligo_pool"] = out["oligo_name"].map(o2l)

    n = len(out)
    m = int(out["oligo_name"].notna().sum())
    unamb = int(out["is_unambiguous"].sum()) if "is_unambiguous" in out.columns else -1
    print(f"Matched {m}/{n} ({100*m/max(n,1):.1f}%) to designed oligos; unambiguous={unamb}")
    out.to_csv(args.output, sep="\t", index=False)
    print(f"Saved mapped data to: {args.output}")


if __name__ == "__main__":
    main()
