#!/usr/bin/env python3
"""
Combine and collapse bc1-donor-bc0 data across samples -- and across sequencing runs.

Rewritten 2026-09-10 for the 20210822 rebuild (planning ``28_``; Kevin's decisions in
``27_`` Sec. 7). What changed, and why each change is load-bearing:

* **LABEL, never filter.** Every bc1-donor-bc0 link is retained. ``support_tier`` =
  ``confident`` (>= ``--min-pcr-replicates`` distinct PCR replicates, default 2) or
  ``provisional``. The old ``num_PCR_replicates >= 2`` row filter deleted the singleton
  tier Kevin asked to see (measured on yL406 at 91.47% of links, ``21_``).
* **Replicate support is an exact SET.** Identity is the approved composite
  ``build_pcr_replicate_id`` (gDNA plate + inner + outer primers + polymerase +
  condition + amplicon level + well; sequencing date EXCLUDED, so re-sequencing a
  library is not new evidence). Sets are unioned through every collapse -- the old
  code summed ``num_PCR_replicates`` across merged rows.
* **One invocation spans every run** of a screen via ``--manifest`` (explicit parsed
  path + replicate fields per sample). The legacy per-``--run-date`` mode remains, but
  any missing input or an empty load is now an ERROR -- the old code printed
  "No data files found!" and exited 0.
* **Per-tier read counts** (``n_<tier>``, ``23_`` option b) are carried as additive
  columns: ``counts == sum(n_<tier>)`` on every row.
* **Cross-read-length reconciliation**: a 2x150 donor (~88 nt) is attached to the
  2x300 donor (238 nt) it is a suffix of (``core.link_collapse.attach_short_donors``).
* **SPS and donor_bc0_fragment survive** (from the top-count row); the old groupby
  dropped them, leaving V_library assignable only for guide-linked bc1s.
* Rows with no donor or bc0 are **set aside and counted**, not dropped by a groupby.

Usage (manifest mode -- the rebuild path):
    python 02_combine_and_collapse.py --manifest manifest.tsv \\
        --output-collapsed combined/bc1_donor_bc0_collapsed.tsv \\
        --output-purity combined/bc1_purity.tsv [--n-jobs 32]

Manifest columns: required ``sample_id parsed_path sequencing_date sample_number
outer_primers``; replicate-identity fields ``gDNA_plate_name inner_primers
plate_polymerase condition amplicon_level`` (absent, blank or "NA" -> missing, warned;
a missing amplicon_level is derived from trafo_ID as the loaders do); optional
passthrough ``trafo_ID V_library nuclease_PAM yeast_sublibrary``.

Legacy mode (one run): --parsed-dir --sample-key --plate-key --run-date [--trafo-key].

Author: Kevin Roy Lab
"""

import argparse
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from rapidfuzz.distance import Levenshtein as _lev

from magestic.utils.bc1_extraction import TIER_COUNT_COLUMNS, BC1_ORIENTATION
from magestic.pipelines.bc1_donor_bc0.readers.data_loader import (
    PCR_REPLICATE_ID_FIELDS, _cross_merge_keys, build_pcr_replicate_id, declared_value,
    derive_amplicon_level)
from magestic.pipelines.bc1_donor_bc0.core.link_collapse import (
    aggregate_links, collapse_links, attach_short_donors, merge_orphan_bc1,
    mask_columns, masks_for_codes, popcount_rows, decode_masks)
from magestic.utils.provenance import StepProvenance, md5_file

TIER_COLS = list(TIER_COUNT_COLUMNS)
SUM_COLS = ["counts", *TIER_COLS]
LINK_KEY = ["bc1", "donor", "bc0"]
SEQ_PRESERVE = ["SPS", "donor_bc0_fragment"]
TRAFO_COLS = ["trafo_ID", "V_library", "nuclease_PAM", "yeast_sublibrary"]
PARSED_REQUIRED = ["counts", "bc1", "donor", "bc0", "SPS", "donor_bc0_fragment", *TIER_COLS]
MANIFEST_REQUIRED = ["sample_id", "parsed_path", "sequencing_date", "sample_number", "outer_primers"]
ACCUMULATE_EVERY = 24   # default samples per incremental fold (sums/ORs are associative; preserved SPS/fragment come
                        # from the top row of each fold, which is batch-order-dependent only on exact read-count ties)
LIBRARY_COL = "yeast_library"   # manifest column -> per-link yeast_libraries_observed (WGS lane ask D.3)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr, flush=True)
    sys.exit(2)


# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------

def build_manifest(args) -> pd.DataFrame:
    """The explicit per-sample input list. Fails loudly on anything missing."""
    if args.manifest:
        man = pd.read_csv(args.manifest, sep="\t", dtype=str, keep_default_na=False)
        missing = [c for c in MANIFEST_REQUIRED if c not in man.columns]
        if missing:
            die(f"manifest lacks required columns {missing}")
    else:
        if not all([args.parsed_dir, args.sample_key, args.plate_key, args.run_date]):
            die("give --manifest, or all of --parsed-dir --sample-key --plate-key --run-date")
        plate_key = pd.read_csv(args.plate_key, sep="\t", dtype=str, keep_default_na=False)
        sample_key = pd.read_csv(args.sample_key, sep="\t", dtype=str, keep_default_na=False)
        plate_key = plate_key[plate_key["description"] == "bc1-donor-bc0"]
        # Same explicit renames and shared-column resolution as the 2x300 loader, so both
        # producers read one pair of keyfiles into the same declarations.
        plate_key = plate_key.rename(columns={"polymerase": "plate_polymerase",
                                              "plate_category": "plate_category_from_plate_key"})
        sample_key = sample_key.rename(columns={"polymerase": "sample_listed_polymerase",
                                                "plate_category": "plate_category_from_sample_key"})
        man = _cross_merge_keys(plate_key, sample_key)
        man["sequencing_date"] = args.run_date
        man["sample_id"] = (args.run_date + "_" + man["outer_primers"].astype(str)
                            + "_" + man["sample_number"].astype(str))
        man["parsed_path"] = [str(Path(args.parsed_dir) / f"{sid}_parsed.tsv") for sid in man["sample_id"]]
    if len(man) == 0:
        die("zero samples in the input definition")
    # normalize: surrounding whitespace/CR, 'nan'/'None' spellings, zero-padded wells -- two rows of one library
    # must produce the SAME replicate id (review REP-07)
    for c in man.columns:
        man[c] = man[c].astype(str).str.strip().replace({"nan": "", "None": "", "NaN": ""})
    man["sample_number"] = man["sample_number"].str.replace(r"^sample_", "", regex=True).map(
        lambda v: str(int(v)) if v.isdigit() else v)
    if man["sample_id"].duplicated().any():
        die(f"duplicate sample_id in manifest: {man.loc[man.sample_id.duplicated(), 'sample_id'].tolist()[:10]}")
    if man["parsed_path"].duplicated().any():
        # the same file twice = a fake second replicate + doubled reads (review F5)
        die(f"the same parsed file is listed twice: {man.loc[man.parsed_path.duplicated(), 'parsed_path'].tolist()[:5]}")

    # optional supplementary trafo key (sample_number -> trafo fields)
    if args.trafo_key:
        tk = pd.read_csv(args.trafo_key, sep="\t", dtype=str, keep_default_na=False)
        if "sample_number" not in tk.columns:
            die(f"--trafo-key lacks sample_number (has {list(tk.columns)})")
        if tk["sample_number"].duplicated().any():
            die("--trafo-key has duplicate sample_number rows")
        cols = [c for c in TRAFO_COLS if c in tk.columns]
        # OVERRIDE per sample where the trafo key has a value; keep the manifest's value
        # elsewhere (review F4: a plain replace silently blanked unmatched samples).
        sup = man[["sample_number"]].merge(tk[["sample_number", *cols]], on="sample_number", how="left")
        for c in cols:
            new = sup[c].fillna("").to_numpy()
            old = man[c].to_numpy() if c in man.columns else [""] * len(man)
            man[c] = [n if n != "" else o for n, o in zip(new, old)]

    absent = [f for f in PCR_REPLICATE_ID_FIELDS if f not in man.columns]
    for f in absent:
        man[f] = ""
    if absent:
        log(f"WARNING: replicate-identity fields absent from the input definition, treated as missing: {absent}")
    # A missing amplicon_level is derived from trafo_ID, as both loaders derive it.
    trafo = man["trafo_ID"] if "trafo_ID" in man.columns else pd.Series("", index=man.index)
    man["amplicon_level"] = [
        level if declared_value(level) is not None else derive_amplicon_level(declared_value(t))
        for level, t in zip(man["amplicon_level"], trafo)]
    no_prep = man["gDNA_plate_name"].map(declared_value).isna()
    if no_prep.any():
        # review F6: partial blanks were silent; two different preps sharing primers + well
        # would then collapse to ONE replicate id.
        log(f"WARNING: gDNA_plate_name unavailable for {int(no_prep.sum())}/{len(man)} samples "
            f"(e.g. {man.loc[no_prep, 'sample_id'].tolist()[:5]}) -> the re-prep guard is inactive "
            f"for them; an independent re-prep that reused the same primers would count as ONE replicate.")
    # The shared builder normalizes blanks and NA spellings to missing, so this producer and
    # the loaders encode the same declared PCR identically.
    man["pcr_replicate_id"] = man.apply(build_pcr_replicate_id, axis=1)

    missing_files = [p for p in man["parsed_path"] if not Path(p).is_file()]
    if missing_files:
        if not args.allow_missing:
            die(f"{len(missing_files)} of {len(man)} parsed files missing (first: {missing_files[:5]}). "
                f"Pass --allow-missing to proceed without them (they are then recorded as missing).")
        log(f"WARNING: {len(missing_files)} parsed files missing -- proceeding (--allow-missing): "
            f"{missing_files[:5]}")
    man["present"] = [Path(p).is_file() for p in man["parsed_path"]]
    if "_key_conflicts" in man.columns:
        conflicted = man["present"] & (man["_key_conflicts"] != "")
        if conflicted.any():
            die(f"plate_key and sample_key declare conflicting values for present samples: "
                f"{man.loc[conflicted, ['sample_id', '_key_conflicts']].head(5).values.tolist()}")
        man = man.drop(columns="_key_conflicts")
    return man


# -----------------------------------------------------------------------------
# Per-sample: load, set aside unusable rows, collapse within sample (Steps 1-3)
# -----------------------------------------------------------------------------

def process_sample(task):
    (sample_id, path, rep_code, n_ids, run_code, n_runs, lib_code, n_libs, trafo, thr) = task
    df = pd.read_csv(path, sep="\t", dtype={"bc1": str, "donor": str, "bc0": str, "SPS": str,
                                             "donor_bc0_fragment": str},
                     keep_default_na=False, na_values=[""])
    missing = [c for c in PARSED_REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: parsed file lacks {missing} -- produced by a pre-2026-09-10 "
                         f"step 01? Re-parse; mixing schemas would silently lose the tier counts.")
    aside = {"sample_id": sample_id, "rows_in": len(df), "reads_in": int(df["counts"].sum())}
    for reason, mask in (("no_bc1", df["bc1"].isna()),
                         ("no_donor", df["bc1"].notna() & df["donor"].isna()),
                         ("no_bc0", df["bc1"].notna() & df["donor"].notna() & df["bc0"].isna())):
        aside[f"rows_{reason}"] = int(mask.sum())
        aside[f"reads_{reason}"] = int(df.loc[mask, "counts"].sum())
    df = df[df["bc1"].notna() & df["donor"].notna() & df["bc0"].notna()]
    for c in SEQ_PRESERVE:
        df[c] = df[c].fillna("")
    if len(df) == 0:
        aside["rows_out"] = 0
        return None, aside

    df = aggregate_links(df, LINK_KEY, "counts", SUM_COLS, [], SEQ_PRESERVE)
    # Step 1: similar donors within (bc1, bc0); Step 2: similar bc1 within (donor, bc0);
    # Step 3: similar bc0 within (bc1, donor) -- all within this one sample.
    df = collapse_links(df, ["bc1", "bc0"], "donor", "counts", SUM_COLS, [], SEQ_PRESERVE,
                        threshold=thr["donor"])
    df = collapse_links(df, ["donor", "bc0"], "bc1", "counts", SUM_COLS, [], SEQ_PRESERVE,
                        threshold=thr["bc1"])
    df = collapse_links(df, ["bc1", "donor"], "bc0", "counts", SUM_COLS, [], SEQ_PRESERVE,
                        threshold=thr["bc0"])
    rm = masks_for_codes(np.full(len(df), rep_code), n_ids)
    for k, c in enumerate(mask_columns(n_ids)):
        df[c] = rm[:, k]
    um = masks_for_codes(np.full(len(df), run_code), n_runs)
    for k, c in enumerate(run_mask_columns(n_runs)):
        df[c] = um[:, k]
    lm = masks_for_codes(np.full(len(df), lib_code), n_libs)
    for k, c in enumerate(lib_mask_columns(n_libs)):
        df[c] = lm[:, k]
    for c, v in trafo.items():
        df[c] = v
    aside["rows_out"] = len(df)
    return df, aside


def run_mask_columns(n_runs):
    return [c.replace("_repmask_", "_runmask_") for c in mask_columns(n_runs)]


def lib_mask_columns(n_libs):
    return [c.replace("_repmask_", "_libmask_") for c in mask_columns(n_libs)]


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Combine and collapse bc1-donor-bc0 data",
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=__doc__)
    parser.add_argument("--manifest", help="Per-sample input definition TSV (rebuild path)")
    parser.add_argument("--parsed-dir", help="[legacy] Directory with parsed TSV files")
    parser.add_argument("--sample-key", help="[legacy] Sample key TSV file")
    parser.add_argument("--plate-key", help="[legacy] Plate key TSV file")
    parser.add_argument("--run-date", help="[legacy] Sequencing run date")
    parser.add_argument("--trafo-key", help="Optional trafo mapping TSV (sample_number -> trafo_ID ...)")
    parser.add_argument("--output-collapsed", required=True, help="Output collapsed TSV")
    parser.add_argument("--output-purity", required=True, help="Output purity summary TSV")
    parser.add_argument("--donor-edit-distance", type=int, default=5)
    parser.add_argument("--bc1-edit-distance", type=int, default=1)
    parser.add_argument("--bc0-edit-distance", type=int, default=1)
    parser.add_argument("--min-pcr-replicates", type=int, default=2,
                        help="LABEL threshold for support_tier=confident (no row is dropped)")
    parser.add_argument("--allow-missing", action="store_true",
                        help="Proceed when some parsed files are absent (recorded, not silent)")
    parser.add_argument("--n-jobs", type=int, default=24)
    parser.add_argument("--no-orphan-merge", action="store_true",
                        help="Skip the orphan-bc1 merge (Kevin 2026-09-10, decision C). By default a bc1 at "
                             "Hamming 1 from another whose bc0 is within one edit and whose donor is compatible "
                             "is merged into it: the three exact-key collapses above cannot catch these because "
                             "a 238-nt donor read rarely repeats exactly (the X1 defect).")
    parser.add_argument("--no-orphan-star-guard", action="store_true",
                        help="Merge whole connected components instead of requiring direct evidence against the "
                             "component root. NOT recommended: accepted pairs form a graph, and transitive "
                             "merging produced components spanning Hamming distance 8.")
    parser.add_argument("--attach-dominance", type=float, default=3.0,
                        help="attach a short donor to the top of several long candidates when it has >= this "
                             "multiple of the runner-up's reads (inf = never; 29_)")
    parser.add_argument("--accumulate-every", type=int, default=ACCUMULATE_EVERY,
                        help="Samples per incremental fold (exact at any value; small values force the "
                             "fold branch in a smoke test)")
    args = parser.parse_args()
    t0 = time.time()

    man = build_manifest(args)
    present = man[man["present"]].reset_index(drop=True)
    rep_ids = sorted(present["pcr_replicate_id"].unique())
    rep_code = {r: i for i, r in enumerate(rep_ids)}
    runs = sorted(present["sequencing_date"].unique())
    run_code = {r: i for i, r in enumerate(runs)}
    if LIBRARY_COL not in present.columns:
        present[LIBRARY_COL] = "NA"
    present[LIBRARY_COL] = present[LIBRARY_COL].replace("", "NA")
    libs = sorted(present[LIBRARY_COL].unique())
    lib_code = {r: i for i, r in enumerate(libs)}
    mcols, ucols, lcols = mask_columns(len(rep_ids)), run_mask_columns(len(runs)), lib_mask_columns(len(libs))
    log(f"{len(present)} samples present ({len(man) - len(present)} missing); "
        f"{len(rep_ids)} distinct PCR replicates; {len(runs)} sequencing runs: {runs}")
    trafo_present = [c for c in TRAFO_COLS if c in present.columns]
    thr = {"donor": args.donor_edit_distance, "bc1": args.bc1_edit_distance, "bc0": args.bc0_edit_distance}

    tasks = [(r.sample_id, r.parsed_path, rep_code[r.pcr_replicate_id], len(rep_ids),
              run_code[r.sequencing_date], len(runs), lib_code[r[LIBRARY_COL]], len(libs),
              {c: r[c] for c in trafo_present}, thr) for _, r in present.iterrows()]

    # --- per-sample collapse (parallel) + incremental cross-sample aggregation ---
    acc, pending, asides = None, [], []
    preserve_x = SEQ_PRESERVE + trafo_present
    or_cols = mcols + ucols + lcols

    def fold(acc, frames):
        stack = frames if acc is None else [acc, *frames]
        return aggregate_links(pd.concat(stack, ignore_index=True), LINK_KEY, "counts",
                               SUM_COLS, or_cols, preserve_x)

    done = 0
    with ProcessPoolExecutor(max_workers=max(1, args.n_jobs)) as ex:
        for df, aside in ex.map(process_sample, tasks, chunksize=1):
            asides.append(aside)
            done += 1
            if df is not None:
                pending.append(df)
            if len(pending) >= args.accumulate_every:
                acc, pending = fold(acc, pending), []
                log(f"  {done}/{len(tasks)} samples folded; {len(acc):,} links so far")
    if pending:
        acc = fold(acc, pending)
    if acc is None or len(acc) == 0:
        die("no usable bc1-donor-bc0 rows in any sample")
    zero = [a["sample_id"] for a in asides if a.get("rows_out", 0) == 0]
    if zero:
        log(f"WARNING: {len(zero)} samples contributed ZERO usable links (e.g. {zero[:5]}); see .setaside.tsv")
    n_links_agg = len(acc)
    log(f"cross-sample aggregation: {n_links_agg:,} links, {acc['bc1'].nunique():,} bc1")

    # --- within-length collapses FIRST, so each long parent is unique before attach
    #     (review F2/F3 2026-09-10: attaching first let a sequencing-error variant of the
    #     238-nt parent -- or of its bc0 -- make the 88-nt donor "ambiguous", leaving the
    #     link split across eras). Short and long donors cannot merge here (length gap). ---
    acc = collapse_links(acc, ["bc1", "bc0"], "donor", "counts", SUM_COLS, or_cols, preserve_x,
                         threshold=args.donor_edit_distance, n_jobs=args.n_jobs)
    acc = collapse_links(acc, ["bc1", "donor"], "bc0", "counts", SUM_COLS, or_cols, preserve_x,
                         threshold=args.bc0_edit_distance, n_jobs=args.n_jobs)

    # --- cross-read-length reconciliation (88-nt 2x150 donors -> their 238-nt parents) ---
    acc, attach_stats = attach_short_donors(acc, "donor", "counts", SUM_COLS, or_cols,
                                            max_mismatch=args.donor_edit_distance,
                                            bc0_max_edit=args.bc0_edit_distance,
                                            dominance=args.attach_dominance)
    log(f"suffix-attach: {attach_stats}; {len(acc):,} links remain")

    # --- Step 4: collapse across all samples, replicate sets UNIONED ---
    px = preserve_x + ["donor_suffix_ambiguous", "donor_lengths_observed"]
    acc = collapse_links(acc, ["bc1", "bc0"], "donor", "counts", SUM_COLS, or_cols, px,
                         threshold=args.donor_edit_distance, n_jobs=args.n_jobs)
    acc = collapse_links(acc, ["donor", "bc0"], "bc1", "counts", SUM_COLS, or_cols, px,
                         threshold=args.bc1_edit_distance, n_jobs=args.n_jobs)
    acc = collapse_links(acc, ["bc1", "donor"], "bc0", "counts", SUM_COLS, or_cols, px,
                         threshold=args.bc0_edit_distance, n_jobs=args.n_jobs)
    log(f"after Step-4 collapses: {len(acc):,} links, {acc['bc1'].nunique():,} bc1")

    # --- Step 5: orphan-bc1 merge (Kevin 2026-09-10, decision C) ---
    # The collapses above group by an EXACT (donor, bc0) key, so a bc1 read error whose donor read also
    # carries an error is never grouped with its parent and survives as its own barcode. Merge it here,
    # beside the other collapses -- not in a post-hoc script against the shipped table.
    orphan_stats, bc1_to_root = {"skipped": True}, pd.DataFrame(columns=["bc1", "bc1_root"])
    if not args.no_orphan_merge:
        acc, orphan_stats, bc1_to_root = merge_orphan_bc1(
            acc, "counts", SUM_COLS, or_cols, px,
            donor_max_mismatch=args.donor_edit_distance,
            star_guard=not args.no_orphan_star_guard)
        log(f"orphan merge: {orphan_stats.get('accepted_pairs', 0):,} accepted pairs removed "
            f"{orphan_stats.get('bc1_removed', 0):,} bc1 and moved {orphan_stats.get('reads_moved', 0):,} reads; "
            f"{len(acc):,} links, {acc['bc1'].nunique():,} bc1 remain")

    # --- replicate support: exact, labeled ---
    masks = acc[mcols].to_numpy(dtype=np.uint64)
    acc["num_PCR_replicates"] = popcount_rows(masks)
    acc["pcr_replicate_ids"] = decode_masks(masks, rep_ids)
    acc["runs_observed"] = decode_masks(acc[ucols].to_numpy(dtype=np.uint64), runs)
    acc["yeast_libraries_observed"] = decode_masks(acc[lcols].to_numpy(dtype=np.uint64), libs)
    acc = acc.drop(columns=or_cols)
    acc["support_tier"] = np.where(acc["num_PCR_replicates"] >= args.min_pcr_replicates,
                                   "confident", "provisional")
    if not (acc["counts"] == acc[TIER_COLS].sum(axis=1)).all():   # explicit: survives python -O
        die("tier invariant broken: counts != sum(n_<tier>) on "
            f"{int((acc['counts'] != acc[TIER_COLS].sum(axis=1)).sum())} links")

    # --- top link per bc1, purity (over ALL links; nothing filtered) ---
    acc = acc.sort_values(["counts", "bc1", "donor", "bc0"], ascending=[False, True, True, True],
                          kind="mergesort").reset_index(drop=True)
    top = acc.drop_duplicates(subset=["bc1"], keep="first")[
        ["bc1", "donor", "bc0", "counts", "num_PCR_replicates"]].rename(columns={
            "donor": "top_donor", "bc0": "top_bc0", "counts": "top_counts",
            "num_PCR_replicates": "top_num_PCR_replicates"})
    acc = acc.merge(top, on="bc1", how="left")
    acc["total_bc1_counts"] = acc.groupby("bc1")["counts"].transform("sum")
    acc["bc1_purity"] = acc["counts"] / acc["total_bc1_counts"]
    acc["num_distinct_donors"] = acc.groupby("bc1")["donor"].transform("nunique")
    acc["edit_distance_to_top"] = [_lev.distance(a, b) for a, b in zip(acc["donor"], acc["top_donor"])]
    acc["bc1_orientation"] = BC1_ORIENTATION
    acc["donor_suffix_ambiguous"] = acc["donor_suffix_ambiguous"].astype(bool)

    lead = ["bc1", "donor", "bc0", "SPS", "donor_bc0_fragment", "counts", *TIER_COLS,
            "num_PCR_replicates", "support_tier", "pcr_replicate_ids", "runs_observed",
            "yeast_libraries_observed",
            "total_bc1_counts", "bc1_purity", "num_distinct_donors", "edit_distance_to_top",
            "top_donor", "top_bc0", "top_counts", "top_num_PCR_replicates",
            "donor_suffix_ambiguous", "donor_lengths_observed", "bc1_orientation", *trafo_present]
    acc = acc[lead]
    os.makedirs(os.path.dirname(os.path.abspath(args.output_collapsed)), exist_ok=True)
    acc.to_csv(args.output_collapsed, sep="\t", index=False)
    log(f"saved {len(acc):,} links -> {args.output_collapsed}")

    purity = acc.groupby("bc1").agg(total_counts=("counts", "sum"),
                                    max_PCR_replicates=("num_PCR_replicates", "max"),
                                    max_purity=("bc1_purity", "max"),
                                    num_donors=("num_distinct_donors", "first")).reset_index()
    purity.to_csv(args.output_purity, sep="\t", index=False)

    # --- accounting / provenance sidecars ---
    aside_df = pd.DataFrame(asides)
    aside_df.to_csv(f"{args.output_collapsed}.setaside.tsv", sep="\t", index=False)
    rep_legend = present.groupby("pcr_replicate_id").agg(
        n_samples=("sample_id", "size"), sample_ids=("sample_id", lambda s: ";".join(sorted(s))),
        sequencing_dates=("sequencing_date", lambda s: ";".join(sorted(set(s))))).reset_index()
    rep_legend.to_csv(f"{args.output_collapsed}.replicate_ids.tsv", sep="\t", index=False)
    # The bc1 -> root map is a first-class OUTPUT, not a diagnostic: every downstream table keyed on bc1
    # (the screen counts matrix above all) must be re-keyed with exactly this map rather than re-deriving
    # the merge and possibly disagreeing with it.
    bc1_to_root.to_csv(f"{args.output_collapsed}.bc1_to_root.tsv", sep="\t", index=False)
    tier_counts = acc["support_tier"].value_counts().to_dict()
    stats = {
        "argv": sys.argv, "manifest": args.manifest,
        "manifest_md5": hashlib.md5(Path(args.manifest).read_bytes()).hexdigest() if args.manifest else None,
        "samples_defined": int(len(man)), "samples_present": int(len(present)),
        "samples_missing": man.loc[~man["present"], "sample_id"].tolist(),
        "pcr_replicate_ids": len(rep_ids), "sequencing_runs": runs,
        "support_threshold": args.min_pcr_replicates, "support_tier_links": tier_counts,
        "support_tier_bc1": acc.drop_duplicates("bc1")["support_tier"].value_counts().to_dict(),
        "links_after_aggregation": n_links_agg, "suffix_attach": attach_stats,
        "links_out": int(len(acc)), "bc1_out": int(acc["bc1"].nunique()),
        "reads_in": int(aside_df["reads_in"].sum()),
        "reads_set_aside": {k: int(aside_df[k].sum()) for k in aside_df.columns if k.startswith("reads_no_")},
        "reads_out": int(acc["counts"].sum()),
        "bc1_orientation": BC1_ORIENTATION, "tier_columns": TIER_COLS,
        "orphan_merge": orphan_stats,
        "bc1_to_root_map": f"{args.output_collapsed}.bc1_to_root.tsv",
        "edit_distances": thr, "attach_dominance": args.attach_dominance,
        "zero_link_samples": [a["sample_id"] for a in asides if a.get("rows_out", 0) == 0],
        "elapsed_s": round(time.time() - t0, 1),
    }
    with open(f"{args.output_collapsed}.stats.json", "w") as fh:
        json.dump(stats, fh, indent=1, default=str)
    set_aside = sum(stats["reads_set_aside"].values())
    if stats["reads_in"] != stats["reads_out"] + set_aside:
        die(f"read accounting broken: in {stats['reads_in']} != out {stats['reads_out']} + "
            f"set aside {set_aside}")
    log(f"support_tier (links): {tier_counts}; reads in {stats['reads_in']:,} = out "
        f"{stats['reads_out']:,} + set aside {set_aside:,}")

    # --- step sidecar (PIPELINE_CORRECTNESS_POLICY.md Sec. 2): every input with its md5, the code sha,
    # the exact argv, every parameter that shapes the result, and the accounting. The parsed inputs are
    # hashed in parallel; a path is not a version.
    prov = StepProvenance(pipeline="bc1_donor_bc0", step="02_combine_and_collapse", output=args.output_collapsed,
                          description="per-sample parsed bc1-donor-bc0 links combined across every sequencing "
                                      "run, collapsed (donor/bc1/bc0 within one edit), short donors attached, "
                                      "orphan bc1 merged, replicate support LABELED (never filtered)")
    if args.manifest:
        prov.add_input(args.manifest, role="manifest", md5=stats["manifest_md5"])
    parsed_paths = present["parsed_path"].tolist()
    with ProcessPoolExecutor(max_workers=max(1, min(args.n_jobs, 16))) as ex:
        parsed_md5s = list(ex.map(md5_file, parsed_paths))
    for (_, r), m in zip(present.iterrows(), parsed_md5s):
        prov.add_input(r["parsed_path"], role="parsed_bc1_donor_bc0", md5=m, sample_id=r["sample_id"],
                       sequencing_date=str(r["sequencing_date"]), pcr_replicate_id=r["pcr_replicate_id"])
    prov.set_parameters(edit_distances=thr, min_pcr_replicates=args.min_pcr_replicates,
                        attach_dominance=args.attach_dominance, accumulate_every=args.accumulate_every,
                        orphan_merge=not args.no_orphan_merge, orphan_star_guard=not args.no_orphan_star_guard,
                        n_jobs=args.n_jobs, bc1_orientation=BC1_ORIENTATION, tier_columns=TIER_COLS)
    prov.set_counts(**{k: stats[k] for k in (
        "samples_defined", "samples_present", "pcr_replicate_ids", "sequencing_runs",
        "links_after_aggregation", "suffix_attach", "orphan_merge", "links_out", "bc1_out",
        "reads_in", "reads_set_aside", "reads_out", "support_tier_links", "support_tier_bc1")})
    prov.add_output(args.output_collapsed, rows=int(len(acc)))
    for extra in (args.output_purity, f"{args.output_collapsed}.stats.json",
                  f"{args.output_collapsed}.replicate_ids.tsv", f"{args.output_collapsed}.bc1_to_root.tsv",
                  f"{args.output_collapsed}.setaside.tsv"):
        if Path(extra).is_file():
            prov.add_output(extra)
    log(f"sidecar: {prov.write()}")
    log("done")


if __name__ == "__main__":
    main()
