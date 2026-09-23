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
        load_and_build_lookup, match_dataframe, match_dataframe_scoped, repair_gdb01_shifted_guides,
        assign_sps_scope, load_sps_mapping)
    from magestic.utils import path_utils
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # .../src
    from magestic.pipelines.bc1_donor_bc0.core.matching import (
        load_and_build_lookup, match_dataframe, match_dataframe_scoped, repair_gdb01_shifted_guides,
        assign_sps_scope, load_sps_mapping)
    from magestic.utils import path_utils
from magestic.utils.provenance import StepProvenance


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


def build_parser():
    """The CLI contract, factored out so a test can assert the defaults (notably that the GDB-01
    bridge is OFF) without running the step."""
    parser = argparse.ArgumentParser(description="Map bc1-donor-bc0 to designed oligos (canonical rich matcher)")
    parser.add_argument("--input", nargs="+", required=True, help="Input TSV file(s) from Step 03 / combine")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--donor-col", default=None, help="Donor column (default: auto donor/top_donor)")
    parser.add_argument("--guide-col", default="guide", help="Guide column if present")
    parser.add_argument("--v-libraries", nargs="*", default=None, help="Restrict harmonized lookup to these V libraries")
    parser.add_argument("--nuclease", nargs="*", default=None, help="Restrict harmonized lookup to these nuclease/PAM types")
    parser.add_argument("--no-pre-trim", action="store_true", help="Disable anchor pre-trim (already-trimmed donors)")
    # GDB-01 bridge. The root fix landed in gdb Step 02 on 2026-09-11 (sequence_anchors.find_re_site),
    # so this consumer-side repair is now OFF by default and exists only so tables built BEFORE that fix
    # stay readable. See PIPELINE_CORRECTNESS_POLICY.md §1: a downstream repair is a bridge, never the fix.
    parser.add_argument("--gdb01-repair", action="store_true",
                        help="BRIDGE, off by default. Restore Cas12a guides that gdb Step 02 extracted 1 nt early "
                             "(planning 30_ GDB-01). Only for guide-donor-bc0 tables built before the 2026-09-11 "
                             "root fix; against corrected tables it should restore ~0 rows and will warn loudly.")
    parser.add_argument("--no-gdb01-repair", action="store_true",
                        help="(deprecated no-op; the repair is off by default since 2026-09-11)")
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
    return parser


def main():
    args = build_parser().parse_args()

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

    # GDB-01 BRIDGE (off unless --gdb01-repair). gdb Step 02 used to shift every G-ending Cas12a guide 1 nt
    # early ('T' + guide[:19]) and bdb 03 copied the shifted string. That is fixed at the source now, so this
    # only runs against pre-2026-09-11 tables; it restores the guide, keeps the raw string, and labels the row.
    if args.no_gdb01_repair:
        print("NOTE: --no-gdb01-repair is a deprecated no-op; the GDB-01 bridge is off by default since 2026-09-11.")
    if args.guide_col in bc1_data.columns and args.gdb01_repair:
        fixed, flag = repair_gdb01_shifted_guides(bc1_data[args.guide_col], tables)
        bc1_data["guide_gdb_raw"] = bc1_data[args.guide_col].where(flag)
        bc1_data[args.guide_col] = fixed
        bc1_data["guide_repaired_gdb01"] = flag
        nr = int(flag.sum()); ng = int(bc1_data[args.guide_col].notna().sum())
        pct = 100 * nr / max(ng, 1)
        print(f"GDB-01 BRIDGE: restored {nr:,} of {ng:,} guide-carrying rows ({pct:.2f}%) to a designed Cas12a guide")
        if nr == 0:
            print("GDB-01 BRIDGE: 0 rows restored -- these tables already carry corrected guides. Drop --gdb01-repair.")
        else:
            print("=" * 100)
            print(f"WARNING: the GDB-01 bridge fired on {nr:,} rows ({pct:.2f}%). That means this input was built by")
            print("         gdb Step 02 BEFORE the 2026-09-11 root fix (sequence_anchors.find_re_site). The bridge")
            print("         repairs only guides it can prove are shifted; it cannot repair the donor or the counts,")
            print("         and it cannot help any other consumer of those tables. The correct action is to re-run")
            print("         gdb Step 02 on the fix and rebuild from there. See PIPELINE_CORRECTNESS_POLICY.md §1.")
            print("=" * 100)

    # --- SPS scoping (core.matching; planning 34_/36_) -------------------------------------------
    # Every read carries the subpool priming site between donor and bc0, naming which design subpool the
    # molecule was amplified from. Matching a row against all subpools at once let the donor tie be broken
    # alphabetically across subpools instead: 313,570 of 4,458,480 named rows with a readable SPS (7.0%)
    # took a name from a subpool their own read did not come from. A post-hoc candidate filter cannot fix
    # an unambiguous row -- the right design was never a candidate -- so the SEARCH is scoped.
    if args.sps_mapping:
        sps2v = load_sps_mapping(args.sps_mapping, args.v_libraries)
        print(f"SPS scoping: {len(sps2v)} SPS -> subpool ({sorted(set(sps2v.values()))})")
        scope = assign_sps_scope(bc1_data, sps2v)
        bc1_data["sps_scope_V"] = scope
        bc1_data["sps_scoped"] = scope.notna()
        n_scoped = int(scope.notna().sum())
        print(f"  {n_scoped:,} of {len(bc1_data):,} rows scoped by their own SPS; "
              f"{len(bc1_data) - n_scoped:,} fall back to the full scope")
        print("  rows per subpool: " + str(scope.value_counts(dropna=False).to_dict()))

        def _build(v):
            if v is None:
                print(f"  matching unscoped rows against the full scope...")
                return tables
            print(f"  building lookup for {v}...")
            return load_and_build_lookup(use_harmonized_table=True, v_libraries=[v],
                                         nuclease_types=args.nuclease)

        out = match_dataframe_scoped(bc1_data, scope, _build, donor_col=donor_col,
                                     guide_col=args.guide_col, pre_trim=not args.no_pre_trim)
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

    # --- step sidecar (PIPELINE_CORRECTNESS_POLICY.md Sec. 2) ---
    prov = StepProvenance(pipeline="bc1_donor_bc0", step="04_map_to_oligos", output=args.output,
                          description="each bc1-donor-bc0 link matched to a designed oligo of the harmonized table "
                                      "(canonical rich matcher); with --sps-mapping the SEARCH is scoped to the "
                                      "subpool the row's own SPS names (planning 34_/36_)")
    for f in args.input:
        prov.add_input(f, role="bc1_linked_guide_donor_bc0")
    prov.add_input(path_utils.HARMONIZED_VARIANTS_FILE, role="harmonized_oligo_table")
    if args.sps_mapping:
        prov.add_input(args.sps_mapping, role="sps_mapping")
    prov.set_parameters(donor_col=donor_col, guide_col=args.guide_col, v_libraries=args.v_libraries,
                        nuclease=args.nuclease, pre_trim=not args.no_pre_trim,
                        gdb01_repair_bridge=bool(args.gdb01_repair), sps_scoping=bool(args.sps_mapping))
    counts = {"rows": int(n), "unique_bc1": (int(nbc1) if nbc1 != "NA" else None),
              "oligo_named": int(m), "unambiguous": int(unamb)}
    if "sps_scoped" in out.columns:
        counts["sps_scoped_rows"] = int(out["sps_scoped"].sum())
        counts["rows_per_sps_scope"] = {str(k): int(v) for k, v in
                                       out["sps_scope_V"].value_counts(dropna=False).items()}
    if "guide_repaired_gdb01" in out.columns:
        counts["gdb01_bridge_repaired_rows"] = int(out["guide_repaired_gdb01"].sum())
    if "guide_match" in out.columns:
        counts["guide_match"] = {str(k): int(v) for k, v in out["guide_match"].value_counts(dropna=False).items()}
    if "edit_evidence" in out.columns:
        # planning 218 item 10: the before/after measure of how many rows the edit certification changes.
        counts["edit_evidence"] = {str(k): int(v) for k, v in out["edit_evidence"].value_counts().items()}
    prov.set_counts(**counts)
    prov.add_output(args.output, rows=int(n))
    print(f"sidecar: {prov.write()}")


if __name__ == "__main__":
    main()
