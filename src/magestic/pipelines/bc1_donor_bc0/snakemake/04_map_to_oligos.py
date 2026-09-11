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
    from magestic.pipelines.bc1_donor_bc0.core.matching import (
        load_and_build_lookup, match_dataframe, repair_gdb01_shifted_guides)
    from magestic.utils import path_utils
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # .../src
    from magestic.pipelines.bc1_donor_bc0.core.matching import (
        load_and_build_lookup, match_dataframe, repair_gdb01_shifted_guides)
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
    parser.add_argument("--no-gdb01-repair", action="store_true",
                        help="Do not restore Cas12a guides that gdb Step 02 extracted 1 nt early (planning 30_ GDB-01)")
    parser.add_argument("--sps-mapping", default=None,
                        help="SPS -> library_id TSV (cols: sps, library_id). When given, each row is matched ONLY "
                             "against the subpool its own SPS names (planning 34_/36_). Without it, behaviour is "
                             "unchanged: every row is matched against the whole --v-libraries scope.")
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

    # GDB-01: gdb Step 02 shifted every G-ending Cas12a guide 1 nt early ('T' + guide[:19]); bdb 03 copied the shifted
    # string. Restore it before matching so guide narrowing works, keep the raw string, and label the row.
    if args.guide_col in bc1_data.columns and not args.no_gdb01_repair:
        fixed, flag = repair_gdb01_shifted_guides(bc1_data[args.guide_col], tables)
        bc1_data["guide_gdb_raw"] = bc1_data[args.guide_col].where(flag)
        bc1_data[args.guide_col] = fixed
        bc1_data["guide_repaired_gdb01"] = flag
        nr = int(flag.sum()); ng = int(bc1_data[args.guide_col].notna().sum())
        print(f"GDB-01 repair: restored {nr:,} of {ng:,} guide-carrying rows ({100*nr/max(ng,1):.2f}%) to a designed Cas12a guide")

    # --- SPS scoping (planning 36_): a bc1-donor-bc0 read carries the subpool priming site between donor and bc0,
    # and it names which design subpool the molecule was amplified from. Matching a row against all subpools at once
    # lets the donor tie be broken by alphabetical order instead: 313,570 named rows (7.0%) ended up named from a
    # subpool their own read did not come from. Post-filtering the candidate list cannot fix that, because an
    # UNAMBIGUOUS row carries no candidate list -- the right design was never a candidate. So scope the search itself.
    sps2v = {}
    if args.sps_mapping:
        m = pd.read_csv(args.sps_mapping, sep="	", dtype=str)
        sps2v = {str(k).upper(): v for k, v in zip(m["sps"], m["library_id"]) if isinstance(v, str) and v}
        if args.v_libraries:
            sps2v = {k: v for k, v in sps2v.items() if v in args.v_libraries}
        print(f"SPS scoping: {len(sps2v)} SPS -> subpool ({sorted(set(sps2v.values()))})")

    if sps2v and "SPS" in bc1_data.columns:
        scope = bc1_data["SPS"].astype(str).str.upper().map(sps2v)
        bc1_data["sps_scope_V"] = scope
        bc1_data["sps_scoped"] = scope.notna()
        print(f"  {int(scope.notna().sum()):,} of {len(bc1_data):,} rows scoped by their own SPS; "
              f"{int(scope.isna().sum()):,} fall back to the full scope")
        print("  rows per subpool: " + str(scope.value_counts(dropna=False).to_dict()))
        parts = []
        for v, sub in bc1_data.groupby(scope.fillna("__UNSCOPED__"), sort=False):
            if v == "__UNSCOPED__":
                t = tables
                print(f"  matching {len(sub):,} unscoped rows against the full scope...")
            else:
                print(f"  building lookup for {v} and matching {len(sub):,} rows...")
                t = load_and_build_lookup(use_harmonized_table=True, v_libraries=[v],
                                          nuclease_types=args.nuclease)
            parts.append(match_dataframe(sub, t, donor_col=donor_col, guide_col=args.guide_col,
                                         pre_trim=not args.no_pre_trim))
            if v != "__UNSCOPED__":
                del t
        out = pd.concat(parts).sort_index()
        del parts
    else:
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
