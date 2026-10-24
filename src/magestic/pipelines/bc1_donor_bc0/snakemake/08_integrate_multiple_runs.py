#!/usr/bin/env python3
"""
Integrate bc1-donor-bc0 data from multiple sequencing runs.

This script combines data from 2x150 bp and 2x300 bp runs to maximize bc1 recovery:
- 2x300 bp runs capture full donor in merged reads
- 2x150 bp runs capture partial donor; full donor recovered via donor-bc0 linkage

Integration strategy:
1. Load collapsed data from each run
2. For 2x150 data, link bc1 to full donor via donor_bc0_fragment
3. Combine datasets
4. Collapse similar sequences across all runs
5. Apply quality filters

Usage:
    python integrate_multiple_runs.py \\
        --input-files <file1.tsv> <file2.tsv> ... \\
        --output <integrated.tsv> \\
        [options]

Author: Kevin Roy Lab
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from rapidfuzz.distance import Levenshtein as _rapidfuzz_lev



def _peel_off(seqs, cvs, threshold):
    """Greedy top-by-count peel-off on plain lists.

    Behavior-preserving reimplementation of the original pandas process_group inner
    loop: repeatedly pick the highest-total-count surviving sequence, collapse every
    surviving sequence within `threshold` edit distance of it, aggregate counts, repeat.
    Uses rapidfuzz `score_cutoff` for early-exit distance. No pandas churn.

    Args:
        seqs: list[str | non-str(NaN)]   sequences in the group (original within-group order)
        cvs:  list[tuple[int, ...]]      per-row count vector aligned to count cols
        threshold: int                    max edit distance to collapse
    Returns:
        list of (rep_local_idx, summed_count_tuple) — one per output cluster
    """
    n = len(seqs)
    alive = [True] * n
    totals = [sum(cv) for cv in cvs]
    ncc = len(cvs[0]) if n else 0
    out = []
    remaining = n
    while remaining:
        # top = FIRST surviving index with the max total (matches pandas idxmax tie-break)
        top = -1
        best = None
        for i in range(n):
            if alive[i] and (best is None or totals[i] > best):
                best = totals[i]
                top = i
        ts = seqs[top]
        agg = list(cvs[top])
        alive[top] = False
        remaining -= 1
        ts_ok = isinstance(ts, str)  # NaN top (non-str) collapses only itself, as in the original
        if ts_ok:
            for i in range(n):
                if alive[i]:
                    si = seqs[i]
                    if isinstance(si, str) and (
                        si == ts or _rapidfuzz_lev.distance(si, ts, score_cutoff=threshold) <= threshold
                    ):
                        cv = cvs[i]
                        for k in range(ncc):
                            agg[k] += cv[k]
                        alive[i] = False
                        remaining -= 1
        out.append((top, tuple(agg)))
    return out


def _collapse_payload(payload, group_cols, seq_col, cc, pcols, threshold):
    """Worker: collapse one group given plain-python payload (no pandas/DataFrames).

    payload = (gkey_values, seqs_list, count_vectors_list, preserve_values_list|None)
    Returns a list of result-row dicts identical in shape to the original output.
    """
    gkey, seqs, cvs, pres = payload
    base = {group_cols[k]: gkey[k] for k in range(len(group_cols))}
    rows = []
    for rep, agg in _peel_off(seqs, cvs, threshold):
        row = dict(base)
        row[seq_col] = seqs[rep]
        for k in range(len(cc)):
            row[cc[k]] = agg[k]
        if pcols and pres is not None:
            pr = pres[rep]
            for k in range(len(pcols)):
                row[pcols[k]] = pr[k]
        rows.append(row)
    return rows


def collapse_sequences(df, group_cols, seq_col, count_cols, edit_distance_threshold,
                       show_progress=True, n_jobs=1, preserve_cols=None):
    """
    Collapse similar sequences within each group based on edit distance threshold.

    Args:
        df: Input DataFrame
        group_cols: Columns to group by
        seq_col: Sequence column to collapse
        count_cols: Count column(s) to aggregate
        edit_distance_threshold: Maximum edit distance for collapsing
        show_progress: Show progress bar
        n_jobs: Number of parallel jobs
        preserve_cols: Additional columns to preserve (take value from top-count row)

    Returns:
        DataFrame with collapsed sequences

    OPTIMIZED 2026-06-02 (behavior-preserving vs the original greedy top-by-count
    peel-off; validated row/count-equal on real data):
      - SINGLETON FAST-PATH: size-1 groups (the large majority of bc1/donor groups)
        skip ALL edit-distance work and dispatch.
      - PURE-LIST peel-off (_peel_off): no per-iteration pandas .apply/.loc/.drop;
        rapidfuzz Levenshtein with score_cutoff for early-exit.
      - BATCHED PLAIN-TUPLE DISPATCH: workers receive small (lists) payloads, never
        pickled per-group DataFrames (the old defect that pinned the parent on IPC
        and starved workers, so --n-jobs gave ~no speedup). Workers are numpy/MKL-free
        (AVX-512 safe); numpy is used only in the parent for factorize/argsort.
    """
    if preserve_cols is None:
        preserve_cols = []
    cc = count_cols if isinstance(count_cols, list) else [count_cols]
    pcols = [c for c in preserve_cols if c in df.columns]
    out_cols = list(group_cols) + [seq_col] + list(cc) + pcols

    n = len(df)
    if n == 0:
        return pd.DataFrame(columns=out_cols)

    # Factorize the group key, then stable-sort rows so each group is a contiguous slice
    # (avoids materializing/pickling one pandas sub-DataFrame per group).
    if len(group_cols) == 1:
        codes, _ = pd.factorize(df[group_cols[0]], sort=False)
    else:
        codes, _ = pd.factorize(
            pd.Series(list(zip(*[df[c].tolist() for c in group_cols])), index=df.index),
            sort=False,
        )
    order = np.argsort(codes, kind="stable")
    codes_s = codes[order]
    seqs_all = df[seq_col].to_numpy(dtype=object)[order]
    counts_all = df[list(cc)].to_numpy()[order]
    gvals_all = df[list(group_cols)].to_numpy(dtype=object)[order]
    pres_all = df[pcols].to_numpy(dtype=object)[order] if pcols else None

    bounds = np.flatnonzero(np.diff(codes_s)) + 1
    starts = np.concatenate(([0], bounds)).tolist()
    ends = np.concatenate((bounds, [n])).tolist()

    ngc = len(group_cols)
    ncc = len(cc)
    npc = len(pcols)
    singleton_rows = []
    payloads = []
    for a, b in zip(starts, ends):
        if b - a == 1:
            row = {group_cols[k]: gvals_all[a][k] for k in range(ngc)}
            row[seq_col] = seqs_all[a]
            crow = counts_all[a]
            for k in range(ncc):
                v = crow[k]
                row[cc[k]] = v.item() if hasattr(v, "item") else v
            if npc:
                pr = pres_all[a]
                for k in range(npc):
                    row[pcols[k]] = pr[k]
            singleton_rows.append(row)
        else:
            gkey = [gvals_all[a][k] for k in range(ngc)]
            seqs = list(seqs_all[a:b])
            cvs = [tuple(int(counts_all[i][k]) for k in range(ncc)) for i in range(a, b)]
            pres = ([[pres_all[i][k] for k in range(npc)] for i in range(a, b)] if npc else None)
            payloads.append((gkey, seqs, cvs, pres))

    worker = partial(_collapse_payload, group_cols=list(group_cols), seq_col=seq_col,
                     cc=list(cc), pcols=pcols, threshold=edit_distance_threshold)
    multi_rows = []
    if payloads:
        if n_jobs and n_jobs > 1 and len(payloads) > 1:
            chunksize = max(1, len(payloads) // (n_jobs * 8))
            with ProcessPoolExecutor(max_workers=n_jobs) as ex:
                it = ex.map(worker, payloads, chunksize=chunksize)
                if show_progress:
                    it = tqdm(it, total=len(payloads), desc=f"Collapsing {seq_col}")
                for rows in it:
                    multi_rows.extend(rows)
        else:
            it = tqdm(payloads, desc=f"Collapsing {seq_col}") if show_progress else payloads
            for p in it:
                multi_rows.extend(worker(p))

    result = pd.DataFrame(singleton_rows + multi_rows)
    return result[out_cols] if len(result) else pd.DataFrame(columns=out_cols)




def main():
    parser = argparse.ArgumentParser(description="Integrate bc1 data from multiple runs")
    parser.add_argument("--input-files", nargs="+", required=True, help="Input collapsed TSV files")
    parser.add_argument("--output", required=True, help="Output integrated TSV")
    parser.add_argument("--donor-edit-distance", type=int, default=5)
    parser.add_argument("--bc1-edit-distance", type=int, default=1)
    parser.add_argument("--min-pcr-replicates", type=int, default=2)
    parser.add_argument("--n-jobs", type=int, default=24)

    args = parser.parse_args()

    print("=" * 60)
    print("BC1-DONOR-BC0 INTEGRATION PIPELINE")
    print("=" * 60)

    # Load all input files
    all_dfs = []
    for f in args.input_files:
        print(f"Loading: {f}")
        df = pd.read_csv(f, sep="\t")
        df["source_file"] = Path(f).name
        all_dfs.append(df)
        print(f"  {len(df)} rows, {df['bc1'].nunique()} unique bc1s")

    combined = pd.concat(all_dfs, ignore_index=True)
    print(f"\nCombined: {len(combined)} total rows, {combined['bc1'].nunique()} unique bc1s")

    # Standardize columns
    if "donor" not in combined.columns and "top_donor" in combined.columns:
        combined["donor"] = combined["top_donor"]
    if "bc0" not in combined.columns and "top_bc0" in combined.columns:
        combined["bc0"] = combined["top_bc0"]

    # Filter for entries with full donor
    combined = combined[combined["donor"].notna()].copy()
    print(f"After filtering for full donor: {combined['bc1'].nunique()} unique bc1s")

    # 2026-09-10 (planning 28_ e6): replicate support is a SET unioned across runs --
    # summing num_PCR_replicates across runs counted a re-sequenced library twice, and
    # the >= min_pcr_replicates FILTER deleted the singleton tier. Now: exact union when
    # every input carries pcr_replicate_ids (step 02 >= 2026-09-10), a loudly-flagged
    # upper bound otherwise; and a support_tier LABEL instead of a filter.
    from magestic.pipelines.bc1_donor_bc0.core.link_collapse import (
        aggregate_links, collapse_links, mask_columns, masks_for_codes, popcount_rows, decode_masks)
    try:
        from magestic.utils.bc1_extraction import TIER_COUNT_COLUMNS
    except ImportError:
        TIER_COUNT_COLUMNS = ()
    tier_cols = [c for c in TIER_COUNT_COLUMNS if c in combined.columns]
    if tier_cols and combined[tier_cols].isna().any().any():
        # mixed-schema inputs (some pre-2026-09-10): NaN would be cast to INT64_MIN (review F3)
        print("WARNING: tier columns missing for some inputs -- dropping n_<tier> from the output")
        combined = combined.drop(columns=tier_cols)
        tier_cols = []
    sum_cols = ["counts", *tier_cols]
    exact = "pcr_replicate_ids" in combined.columns and combined["pcr_replicate_ids"].notna().all()
    if exact:
        combined = combined.reset_index(drop=True)
        ex = combined["pcr_replicate_ids"].astype(str).str.split(";").explode()
        ex = ex[ex != ""]
        # A pcr:v1: identity and a legacy underscore-joined identity for the same PCR are
        # different strings, so unioning them would count one PCR twice.
        v1 = ex.str.startswith("pcr:v1:")
        if v1.any() and not v1.all():
            raise ValueError(
                "inputs mix pcr:v1: and legacy PCR replicate-id encodings; rebuild every "
                "step-02 output with the same code version before integrating runs")
        universe = sorted(set(ex))
        code = {r: k for k, r in enumerate(universe)}
        or_cols = mask_columns(len(universe))
        m = np.zeros((len(combined), len(or_cols)), dtype=np.uint64)
        np.bitwise_or.at(m, ex.index.to_numpy(), masks_for_codes(ex.map(code).to_numpy(), len(universe)))
        for k, c in enumerate(or_cols):
            combined[c] = m[:, k]
        print(f"Replicate support: EXACT union over {len(universe)} replicate ids")
    else:
        or_cols = []
        sum_cols.append("num_PCR_replicates")
        print("WARNING: inputs lack pcr_replicate_ids (pre-2026-09-10 step 02) -- num_PCR_replicates "
              "is SUMMED across runs and is an UPPER BOUND (a re-sequenced library counts twice).")
    combined["source_file"] = combined["source_file"].astype(str)
    agg_df = aggregate_links(combined, ["bc1", "donor", "bc0"], "counts", sum_cols, or_cols, [])
    src = combined.groupby(["bc1", "donor", "bc0"])["source_file"].agg(lambda x: ",".join(sorted(set(x))))
    agg_df = agg_df.merge(src.rename("source_file").reset_index(), on=["bc1", "donor", "bc0"], how="left")

    print(f"\nAfter aggregation: {len(agg_df)} rows")

    print("\nCollapsing similar donors / bc1s / bc0s (replicate sets unioned)...")
    collapsed_bc0 = collapse_links(agg_df, ["bc1", "bc0"], "donor", "counts", sum_cols, or_cols,
                                   ["source_file"], threshold=args.donor_edit_distance, n_jobs=args.n_jobs)
    collapsed_bc0 = collapse_links(collapsed_bc0, ["bc0", "donor"], "bc1", "counts", sum_cols, or_cols,
                                   ["source_file"], threshold=args.bc1_edit_distance, n_jobs=args.n_jobs)
    collapsed_bc0 = collapse_links(collapsed_bc0, ["bc1", "donor"], "bc0", "counts", sum_cols, or_cols,
                                   ["source_file"], threshold=1, n_jobs=args.n_jobs)
    if exact:
        mm = collapsed_bc0[or_cols].to_numpy(dtype=np.uint64)
        collapsed_bc0["num_PCR_replicates"] = popcount_rows(mm)
        collapsed_bc0["pcr_replicate_ids"] = decode_masks(mm, universe)
        collapsed_bc0 = collapsed_bc0.drop(columns=or_cols)
    collapsed_bc0["num_PCR_replicates_is_upper_bound"] = not exact

    # Sort and get top per bc1 (deterministic on ties)
    collapsed_bc0 = collapsed_bc0.sort_values(["counts", "bc1", "donor", "bc0"],
                                              ascending=[False, True, True, True], kind="mergesort")

    top_per_bc1 = collapsed_bc0.drop_duplicates(subset=["bc1"], keep="first").copy()
    top_per_bc1 = top_per_bc1.rename(columns={
        "donor": "top_donor",
        "bc0": "top_bc0",
        "counts": "top_counts",
        "num_PCR_replicates": "top_num_PCR_replicates"
    })

    result = collapsed_bc0.merge(
        top_per_bc1[["bc1", "top_donor", "top_bc0", "top_counts", "top_num_PCR_replicates"]],
        on="bc1",
        how="left"
    )

    # LABEL, never filter (Kevin 2026-09-10: N=2 as a label; every link retained)
    filtered = result.copy()
    filtered["support_tier"] = np.where(filtered["num_PCR_replicates"] >= args.min_pcr_replicates,
                                        "confident", "provisional")
    print(f"\nsupport_tier at >= {args.min_pcr_replicates} PCR replicates (no row dropped): "
          f"{filtered.drop_duplicates('bc1')['support_tier'].value_counts().to_dict()} bc1s")

    # Add metrics
    filtered["total_bc1_counts"] = filtered.groupby("bc1")["counts"].transform("sum")
    filtered["bc1_purity"] = filtered["counts"] / filtered["total_bc1_counts"]
    filtered["num_distinct_donors"] = filtered.groupby("bc1")["donor"].transform("nunique")

    # Save
    filtered.to_csv(args.output, sep="\t", index=False)
    print(f"\nSaved integrated data to: {args.output}")

    # Summary
    print("\n" + "=" * 60)
    print("INTEGRATION SUMMARY")
    print("=" * 60)
    print(f"Total unique bc1s: {filtered['bc1'].nunique()}")
    print(f"Total rows: {len(filtered)}")
    print(f"Total counts: {filtered['counts'].sum()}")


if __name__ == "__main__":
    main()
