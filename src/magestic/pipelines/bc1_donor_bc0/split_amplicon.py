"""
Split-amplicon bc1 reference table: the route that built the NRD1 2025 screen table.

LAYOUT. The 2025 NRD1 library (yL435) was linked with TWO amplicons sequenced 2x150 on AVITI, because a
2x150 read pair cannot span bc1 -> full donor -> bc0 (the middle ~41 nt of the donor is never read):

* **bc1-donor-bc0 fragment** (KR1967-KR1590, plate 5B, unmerged R1/R2): R1 carries bc1 at fixed offsets
  after ``CAGGTCATGCTC``; R2 carries msd | bc0 | SPS | the bc0-proximal donor end. The join key it yields is
  ``donor_bc0_fragment`` = the last 30 nt of the donor + SPS + bc0 (60 nt).
* **donor-bc0** (KR1882-KR1883, plate 1C, merged + collapsed): the FULL donor + SPS + bc0, so each
  ``donor_bc0_fragment`` is resolved to a full donor and then to a designed oligo.
* **guide-donor-bc0** plasmid map (the GDB pipeline's ``V625_guide_donor_bc0_purity_filtered.tsv``) adds the
  guide side through the same 60-nt key.

STAGES (the legacy step letters of ``projects/NNS/20250813_NRD1_MAGESTIC_bc0_donor_bc1/scripts/`` are kept, so
a reader can put the two side by side):

    b  parse_bc1_fragment_pair      per-sample R1/R2 -> fragment counts
    e  parse_donor_bc0_collapsed    per-sample collapsed merged reads -> donor-bc0 records
    c  combine_bc1_fragments        all samples -> edit-distance collapse (<=3) -> purity -> top fragment per bc1
    f  combine_donor_bc0            all samples -> summed donor-bc0 counts
    g  map_donor_bc0_to_oligos      donor -> designed oligo (perfect donor, else the central 60 nt)
    h  link_reference_table         c + g + GDB + oligo annotations -> the reference table

TIE ORDER, AND WHY THE LEGACY TABLE CANNOT BE REBUILT BIT FOR BIT (measured, planning/378). Stages c and g
choose a representative among TIED rows by the order an earlier sort left them in. The legacy scripts used
pandas' default sort (numpy quicksort), whose order among ties depends on the CPU: numpy >= 1.25 sorts with
AVX-512 kernels where the CPU has them. The same legacy scripts on the same inputs gave tables that differ in
1,758 bc1s at stage c and in 5,339 donor_bc0_fragments at stage g between an AMD Milan and an Intel Skylake
node. Stage c also gathered results with ``Pool.imap_unordered``, so its per-bc1 top fragment among fragments
tied on counts depended on worker scheduling (117 bc1s between two runs on one node).

* ``tie_order="designed"`` (the runner's default; Kevin 2026-11-02): as ``stable``, and in addition, whenever the
  collapse or the final top-fragment choice has to pick among EQUALLY supported fragments, it prefers one that
  the donor-bc0 map (stage g) or the GDB plasmid map resolves to a designed donor: rank 0 a perfect full donor,
  rank 1 a perfect central window, rank 2 neither. Rows that are not tied are never affected. Without it, a bc1
  whose reads are two one-read fragments three edits apart was named by whichever sorted first, often the one
  carrying the sequencing errors, and so lost its designed edit (planning/378 §2c).
* ``tie_order="stable"``: every sort is stable and results are gathered in group order, so ties go to the
  earlier row of a fully specified order and the table is a function of the inputs alone, on any CPU.
* ``tie_order="legacy"``: the legacy sorts, so the route equals the legacy scripts run on the same CPU class,
  except stage c's scheduling ties (measured: 105 bc1s, all listed in the ties file).

Every bc1 whose top fragment was a tie is listed in ``<out>/data_tables/top_fragment_ties.tsv``; the table
schema is unchanged.

Behaviours kept for parity and documented rather than changed (changing any of them changes the table):

* step b's fixed-offset fallbacks (bc1 at 16:36, SaCas9/SceI fragment at 36) when an anchor is absent;
* step g's central-window index stops one window short (``range(len(donor) - 60)``), and takes the FIRST oligo
  in oligo-order file order when a donor matches several (flagged ``unambiguous_oligo_name = False``);
* step h's ``multiple_donor_bc0_fragments_passing_filter`` is True when the bc1 is UNIQUE in the table (its name
  says the opposite). The delivered table also carries ``bc1_unique``, equal to it on all 36,628 rows; the
  surviving script does not write it, so it is restored here from the delivered header;
* TSV round trips between stages use pandas defaults, because the delivered floats carry the default parser's
  last-digit rounding (e.g. bc1_purity ...803 in stage c becomes ...804 in the table).

Author: Kevin R. Roy
"""

from __future__ import annotations

import argparse
import json
import logging
import multiprocessing
import os
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

logger = logging.getLogger(__name__)

# --------------------------------------------------------------------------------------------- constants
BC1_LENGTH = 20
FWD_BC1_PREFIX = "CAGGTCATGCTC"
SACAS9_SCEI_FRAGMENT = "CCCTAGGATACGTGC"  # just downstream of bc1
RE_SITE = "GTTTGAAGAGC"  # SpCas9 cloning site, upstream of the donor
EDIT_DISTANCE_THRESHOLD = 3
MIDDLE_DONOR_START_IDX = 30
MIDDLE_DONOR_END_IDX = 90
MIDDLE_DONOR_LENGTH = MIDDLE_DONOR_END_IDX - MIDDLE_DONOR_START_IDX
MIN_BC1_PURITY = 0.05
MIN_PCR_REPLICATES = 2

FRAGMENT_SUFFIX = "_bc1_donor_bc0_fragment_counts.tsv"
DONOR_BC0_SUFFIX = "_complete_donor_bc0_counts.tsv"
FRAGMENT_HEADER = ["counts", "PS1", "bc1", "SaCas9_SceI_fragment", "SaCas9_SceI_fragment_found",
                   "bc1_proximal_MAGESTIC_3_homology", "donor_fragment", "SPS", "bc0", "msd", "donor_bc0_fragment"]
DONOR_BC0_HEADER = ["counts", "RE_site_found", "RE_site", "donor", "SPS", "bc0", "donor_bc0_fragment"]

_COMP = str.maketrans("ACGTacgtRYMKrymkVBHDvbhdNn ", "TGCAtgcaYRKMyrkmBVDHbvdhNn ")


def rev_comp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def _levenshtein():
    from rapidfuzz.distance import Levenshtein  # same metric as python-Levenshtein's `distance`
    return Levenshtein.distance


# =========================================================================================== stage b
def _split_r1(seq: str) -> list:
    if FWD_BC1_PREFIX in seq:
        bc1_start = seq.index(FWD_BC1_PREFIX) + len(FWD_BC1_PREFIX)
        bc1 = seq[bc1_start:bc1_start + BC1_LENGTH]
        ps1 = seq[:bc1_start]
    else:
        bc1 = seq[16:16 + BC1_LENGTH]
        ps1 = seq[:16]
    idx = seq.find(SACAS9_SCEI_FRAGMENT)
    found = True
    if idx == -1:
        idx, found = 36, False
    return [ps1, bc1, seq[idx:idx + len(SACAS9_SCEI_FRAGMENT)], found, seq[idx:]]


def _split_r2(seq: str) -> list:
    return [rev_comp(e) for e in (seq[42:], seq[22:42], seq[12:22], seq[:12], seq[12:72])]


def parse_bc1_fragment_pair(r1_path: str, r2_path: str, out_path: str) -> int:
    """Stage b. Count every distinct (R1 parts + R2 parts) tuple; write by count, descending, ties in first-seen
    order (Python's sort is stable). Returns the number of read pairs."""
    import pysam
    counts: Dict[tuple, int] = {}
    n = 0
    with pysam.FastxFile(r1_path) as r1, pysam.FastxFile(r2_path) as r2:
        for a, b in zip(r1, r2):
            n += 1
            key = tuple(_split_r1(a.sequence) + _split_r2(b.sequence))
            counts[key] = counts.get(key, 0) + 1
    with open(out_path, "w") as fh:
        fh.write("\t".join(FRAGMENT_HEADER) + "\n")
        for seqs, c in sorted(counts.items(), key=lambda x: x[1], reverse=True):
            fh.write(f"{c}\t" + "\t".join(str(e) for e in seqs) + "\n")
    return n


# =========================================================================================== stage e
def _split_merged(seq: str) -> list:
    seq = rev_comp(seq)
    idx = seq.find(RE_SITE)
    found = True
    if idx == -1:
        idx, found = 3, False
    return [found, seq[idx:idx + len(RE_SITE)], seq[idx + len(RE_SITE):-39], seq[-39:-19], seq[-19:-9],
            seq[-69:-9]]


def parse_donor_bc0_collapsed(collapsed_fastq: str, out_path: str) -> int:
    """Stage e. One output row per collapsed record (the read NAME is its count), in file order."""
    import pysam
    n = 0
    with open(out_path, "w") as fh:
        fh.write("\t".join(DONOR_BC0_HEADER) + "\n")
        with pysam.FastxFile(collapsed_fastq) as fq:
            for rec in fq:
                n += 1
                fh.write(f"{int(rec.name)}\t" + "\t".join(str(e) for e in _split_merged(rec.sequence)) + "\n")
    return n


# =========================================================================================== stage c
def _collapse_group(args):
    """Greedy collapse inside one group: take the top row (by count, or by the row sum of several count
    columns), absorb every sequence within the edit-distance threshold, repeat. Verbatim legacy logic."""
    group, group_cols, seq_col, count_col, thr = args
    distance = _levenshtein()
    group = group.copy()
    rows = []
    ranked = "_rank" in group.columns
    while not group.empty:
        score = group[count_col].sum(axis=1) if isinstance(count_col, list) else group[count_col]
        if ranked:  # among rows tied on the score, the best designed rank, then group order
            top_idx = group.loc[score == score.max(), "_rank"].idxmin()
        else:
            top_idx = score.idxmax()
        top_seq = group.loc[top_idx, seq_col]
        group["edit_distance"] = group[seq_col].apply(lambda x: distance(x, top_seq))
        mask = group["edit_distance"] <= thr
        cols = count_col if isinstance(count_col, list) else [count_col]
        summed = {c: group.loc[mask, c].sum() for c in cols}
        rows.append({**{c: group.iloc[0][c] for c in group_cols}, seq_col: top_seq, **summed})
        group = group.loc[~mask].drop(columns=["edit_distance"])
    return rows


def collapse_sequences(df: pd.DataFrame, group_cols: List[str], seq_col: str, count_col, thr: int,
                       n_jobs: int = 1) -> pd.DataFrame:
    """Collapse within each group. Results are gathered IN GROUP ORDER (``imap``), never ``imap_unordered``,
    so the output order -- and every tie broken by it downstream -- is a function of the input alone."""
    args = [(g, group_cols, seq_col, count_col, thr) for _, g in df.groupby(group_cols)]
    rows: List[dict] = []
    if n_jobs > 1 and len(args) > 1:
        with multiprocessing.get_context("fork").Pool(n_jobs) as pool:
            for r in pool.imap(_collapse_group, args, chunksize=max(1, min(512, len(args) // (n_jobs * 8) or 1))):
                rows.extend(r)
    else:
        for a in args:
            rows.extend(_collapse_group(a))
    return pd.DataFrame(rows)


TIE_ORDERS = ("designed", "stable", "legacy")


def _kind(tie_order: str) -> str:
    if tie_order not in TIE_ORDERS:
        raise ValueError(f"tie_order must be one of {TIE_ORDERS}, not {tie_order!r}")
    return "quicksort" if tie_order == "legacy" else "stable"


def designed_fragment_rank(donor_mapping: str, gdb_table: Optional[str] = None) -> Dict[str, int]:
    """``donor_bc0_fragment -> 0`` (resolves to a perfect designed donor) or ``1`` (a perfect central window),
    from the stage g map and, when given, the GDB plasmid map. Fragments absent from the dict rank 2."""
    rank: Dict[str, int] = {}
    srcs = [(donor_mapping, "db_perfect_donor", "db_perfect_middle_donor_region")]
    if gdb_table:
        srcs.append((gdb_table, "perfect_donor", "perfect_middle_donor_region"))
    for path, pcol, mcol in srcs:
        d = pd.read_csv(path, sep="\t", usecols=["donor_bc0_fragment", pcol, mcol])
        d = d[d[pcol].astype(str).eq("True") | d[mcol].astype(str).eq("True")]
        r = (~d[pcol].astype(str).eq("True")).astype(int)
        for f, v in zip(d["donor_bc0_fragment"], r):
            if v < rank.get(f, 2):
                rank[f] = int(v)
    return rank


def combine_bc1_fragments(key: pd.DataFrame, fragment_dir: str, out_dir: str, n_jobs: int = 1,
                          tie_order: str = "stable", rank: Optional[Dict[str, int]] = None) -> Dict[str, int]:
    """Stage c. Writes ``top_bc1_donor_bc0_fragments.tsv`` and ``top_fragment_ties.tsv`` in ``out_dir``.
    ``tie_order="designed"`` needs ``rank`` (``designed_fragment_rank``)."""
    _kind(tie_order)
    designed = tie_order == "designed"
    if designed and rank is None:
        raise ValueError("tie_order='designed' needs the designed-fragment rank (run stage g first)")

    def with_rank(df):
        if designed:
            df["_rank"] = df["donor_bc0_fragment"].map(rank).fillna(2).astype(int)
        return df

    key = key[key["description"] == "bc1-donor-bc0"]
    dfs = []
    for _, row in key.iterrows():
        p = os.path.join(fragment_dir, f"{row['sequencing_date']}_{row['outer_primers']}_{row['sample_number']}"
                         f"{FRAGMENT_SUFFIX}")
        if os.path.exists(p):
            d = pd.read_csv(p, sep="\t")
            d["sample_number"] = row["sample_number"]
            dfs.append(d)
    if not dfs:
        raise FileNotFoundError(f"no bc1-donor-bc0 fragment tables for the key's samples in {fragment_dir}")
    combined = pd.concat(dfs, ignore_index=True).merge(key, how="left")
    combined = (combined.groupby(["sample_number", "polymerase", "bc1", "donor_bc0_fragment"])["counts"]
                .sum().reset_index())

    per_sample = collapse_sequences(with_rank(combined), ["sample_number", "bc1"], "donor_bc0_fragment", "counts",
                                    EDIT_DISTANCE_THRESHOLD, n_jobs)
    per_sample = per_sample.merge(key, how="left")

    allsamp = per_sample.groupby(["bc1", "donor_bc0_fragment"])["counts"].sum().reset_index()
    allsamp["num_PCR_replicates"] = (per_sample.groupby(["bc1", "donor_bc0_fragment"])["sample_number"]
                                     .nunique().reset_index(drop=True))
    allsamp.sort_values(by=["num_PCR_replicates"], ascending=False, inplace=True, kind=_kind(tie_order))

    frag = collapse_sequences(with_rank(allsamp), ["bc1"], "donor_bc0_fragment", ["counts", "num_PCR_replicates"],
                              EDIT_DISTANCE_THRESHOLD, n_jobs)
    if designed:  # equal counts: the designed fragment first; then imap order (stable)
        frag["designed_rank"] = frag["donor_bc0_fragment"].map(rank).fillna(2).astype(int)
        frag = frag.sort_values(by=["counts", "designed_rank"], ascending=[False, True], kind="stable")
    else:
        frag.sort_values(by=["counts"], ascending=False, inplace=True, kind="stable")  # imap order: no scheduling ties
    frag.rename(columns={"counts": "total_bc1_donor_bc0_counts"}, inplace=True)

    # Ties the legacy code broke by worker scheduling: report them.
    mx = frag.groupby("bc1")["total_bc1_donor_bc0_counts"].transform("max")
    at_max = frag[frag["total_bc1_donor_bc0_counts"] == mx]
    tied_bc1 = at_max["bc1"].value_counts()
    tied_bc1 = tied_bc1[tied_bc1 > 1].index
    ties = at_max[at_max["bc1"].isin(tied_bc1)]
    frag = frag.drop(columns=["designed_rank"], errors="ignore")  # reported in the ties file, not the table

    top = frag.drop_duplicates(subset=["bc1"], keep="first").reset_index(drop=True)
    top = top.rename(columns={"donor_bc0_fragment": "top_bc1_donor_fragment",
                              "total_bc1_donor_bc0_counts": "top_total_bc1_donor_bc0_counts",
                              "num_PCR_replicates": "top_num_PCR_replicates"})
    frag = frag.merge(top, how="left")
    ties = ties.merge(top[["bc1", "top_bc1_donor_fragment"]], how="left")

    out = frag[frag["num_PCR_replicates"] >= MIN_PCR_REPLICATES].copy()
    distance = _levenshtein()
    out["edit_distance"] = [distance(a, b) for a, b in zip(out["donor_bc0_fragment"], out["top_bc1_donor_fragment"])]
    out["num_distinct_donor_fragments"] = out.groupby("bc1")["donor_bc0_fragment"].transform("nunique")
    out["total_bc1_counts"] = out.groupby("bc1")["total_bc1_donor_bc0_counts"].transform("sum")
    out["bc1_purity"] = out["total_bc1_donor_bc0_counts"] / out["total_bc1_counts"]
    os.makedirs(out_dir, exist_ok=True)
    out.to_csv(os.path.join(out_dir, "top_bc1_donor_bc0_fragments.tsv"), sep="\t", index=False)
    ties.to_csv(os.path.join(out_dir, "top_fragment_ties.tsv"), sep="\t", index=False)
    return {"samples": len(dfs), "bc1_all": int(frag["bc1"].nunique()), "rows_out": int(len(out)),
            "bc1_out": int(out["bc1"].nunique()), "bc1_top_tied": int(len(tied_bc1)),
            "bc1_top_tied_in_output": int(out[out["bc1"].isin(tied_bc1)]["bc1"].nunique())}


# =========================================================================================== stage f
def combine_donor_bc0(key: pd.DataFrame, donor_bc0_dir: str, out_dir: str) -> Dict[str, int]:
    """Stage f. Sum donor-bc0 records over every donor-bc0 sample of the key."""
    key = key[key["description"] == "donor-bc0"]
    dfs = []
    for _, row in key.iterrows():
        p = os.path.join(donor_bc0_dir, f"{row['sequencing_date']}_{row['outer_primers']}_{row['sample_number']}"
                         f"{DONOR_BC0_SUFFIX}")
        if os.path.exists(p):
            d = pd.read_csv(p, sep="\t")
            d["sample_number"] = row["sample_number"]
            dfs.append(d)
    if not dfs:
        raise FileNotFoundError(f"no donor-bc0 tables for the key's samples in {donor_bc0_dir}")
    combined = pd.concat(dfs, ignore_index=True).merge(key, how="left")
    combined = combined.drop(columns=["sequencing_date", "inner_product_size_bp", "outer_product_size_bp"],
                             errors="ignore")
    grouped = combined.groupby(["RE_site_found", "RE_site", "donor", "SPS", "bc0", "donor_bc0_fragment"],
                               as_index=False).sum(numeric_only=True)
    os.makedirs(out_dir, exist_ok=True)
    grouped.to_csv(os.path.join(out_dir, "combined_donor_bc0_counts.tsv"), sep="\t", index=False)
    return {"samples": len(dfs), "rows_out": int(len(grouped))}


# =========================================================================================== stage g
def _oligo_lookups(oligo_pool: pd.DataFrame):
    pool = oligo_pool.copy()
    pool["donor"] = pool["oligo_seq"].str.slice(51, 180).str.upper()
    perfect: Dict[str, List[str]] = {}
    partial: Dict[str, List[str]] = {}
    for donor, name in zip(pool["donor"], pool["oligo_name"]):
        perfect.setdefault(donor, []).append(name)
        for i in range(0, len(donor) - MIDDLE_DONOR_LENGTH):  # legacy bound: the last window is never indexed
            partial.setdefault(donor[i:i + MIDDLE_DONOR_LENGTH], []).append(name)
    return perfect, partial


def map_donor_bc0_to_oligos(combined_donor_bc0: str, oligo_pool_path: str, out_dir: str,
                            tie_order: str = "stable") -> Dict[str, int]:
    """Stage g. One row per donor_bc0_fragment: its highest-count full donor (ties resolved toward a perfect
    donor, then a perfect central window) and the designed oligo that donor names."""
    df = pd.read_csv(combined_donor_bc0, sep="\t")
    perfect, partial = _oligo_lookups(pd.read_csv(oligo_pool_path, sep="\t"))
    df.sort_values("counts", ascending=False, inplace=True, kind=_kind(tie_order))
    keep = df.groupby("donor_bc0_fragment")["counts"].transform("max") == df["counts"]
    df = df[keep].reset_index(drop=True)
    donors = df["donor"].values
    mids = [d[MIDDLE_DONOR_START_IDX:MIDDLE_DONOR_END_IDX] for d in donors]
    is_perfect = [d in perfect for d in donors]
    is_mid = [m in partial for m in mids]
    names, unamb = [], []
    for d, m, p, q in zip(donors, mids, is_perfect, is_mid):
        name, u = "NA", False
        if p:
            cands = perfect[d]
        elif q:
            cands = partial[m]
        else:
            cands = None
        if cands is not None:
            u = len(cands) == 1
            name = cands[0]
        names.append(name)
        unamb.append(u)
    res = pd.DataFrame({"donor_bc0_fragment": df["donor_bc0_fragment"].values, "db_counts": df["counts"].values,
                        "db_perfect_donor": is_perfect, "db_perfect_middle_donor_region": is_mid,
                        "db_donor": donors, "db_oligo_name": names, "unambiguous_oligo_name": unamb})
    res = (res.sort_values(by=["donor_bc0_fragment", "db_counts", "db_perfect_donor",
                               "db_perfect_middle_donor_region"], ascending=[True, False, False, False])
           .drop_duplicates(subset=["donor_bc0_fragment"]).reset_index(drop=True))
    res = res.sort_values(by=["db_counts"], ascending=False, kind=_kind(tie_order))
    os.makedirs(out_dir, exist_ok=True)
    res.to_csv(os.path.join(out_dir, "donor_bc0_to_oligo_mapping.tsv"), sep="\t", index=False)
    return {"rows_out": int(len(res)), "ambiguous_name": int((~res["unambiguous_oligo_name"]
                                                               & (res["db_oligo_name"] != "NA")).sum())}


# =========================================================================================== stage h
def link_reference_table(top_fragments: str, donor_mapping: str, gdb_table: str, oligo_info: str,
                         out_path: str, min_purity: float = MIN_BC1_PURITY) -> Dict[str, int]:
    """Stage h. bc1 fragment rows with purity > ``min_purity``, joined on ``donor_bc0_fragment`` to the full
    donor (stage g) and to the plasmid guide-donor-bc0 map (columns prefixed ``gdb_``), then annotated with
    the designed oligo's information by ``db_oligo_name``. A fragment that several plasmid rows share yields
    one table row per plasmid row (36,628 rows for 34,500 bc1 in the 2025 table)."""
    bdb = pd.read_csv(top_fragments, sep="\t")
    db = pd.read_csv(donor_mapping, sep="\t")
    gdb = pd.read_csv(gdb_table, sep="\t")
    info = pd.read_csv(oligo_info, sep="\t")
    bdb = bdb[bdb["bc1_purity"] > min_purity].copy()
    gdb = gdb.rename(columns=lambda x: "gdb_" + x if x != "donor_bc0_fragment" else x)
    res = pd.merge(bdb, db, how="left")
    res = pd.merge(res, gdb, how="left", on="donor_bc0_fragment")
    res["donor_bc0_df_matched"] = res["db_donor"].notnull()
    res["guide_donor_bc0_df_matched"] = res["gdb_donor"].notnull()
    res["db_oligo_name"] = res["db_oligo_name"].fillna(res["gdb_oligo_name"])
    res["gdb_oligo_name"] = res["gdb_oligo_name"].fillna(res["db_oligo_name"])
    res = pd.merge(res, info, how="left", left_on="db_oligo_name", right_on="oligo_name")
    unique = ~res["bc1"].duplicated(keep=False)
    res["bc1_unique"] = unique
    res["multiple_donor_bc0_fragments_passing_filter"] = unique  # legacy name; True means the bc1 is UNIQUE
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    res.to_csv(out_path, sep="\t", index=False)
    return {"rows_out": int(len(res)), "bc1_out": int(res["bc1"].nunique()),
            "with_donor_match": int(res["donor_bc0_df_matched"].sum()),
            "with_gdb_match": int(res["guide_donor_bc0_df_matched"].sum())}


# =========================================================================================== runner
def _samples(directory: str, suffix: str) -> List[str]:
    return sorted(f[: -len(suffix)] for f in os.listdir(directory) if f.endswith(suffix))


def _run_b(args):
    r1, r2, out = args
    return out, parse_bc1_fragment_pair(r1, r2, out)


def _run_e(args):
    fq, out = args
    return out, parse_donor_bc0_collapsed(fq, out)


def run(cfg: dict, stages: Sequence[str] = ("b", "e", "f", "g", "c", "h"), n_jobs: int = 1) -> dict:
    """Run the requested stages. ``cfg`` keys: ``out_dir``, ``key``, ``fragment_fastq_dir`` (trimmed R1/R2),
    ``donor_bc0_fastq_dir`` (collapsed), ``gdb_table``, ``oligo_pool``, ``oligo_info``, ``table_name``; optional
    ``fragment_counts_dir`` / ``donor_bc0_counts_dir`` to start from existing stage b / e outputs, and
    ``tie_order`` (``designed``, the default, ``stable`` or ``legacy``; see the module docstring). Stages run in
    the order b, e, f, g, c, h whatever order they are named in, because ``designed`` reads stage g in stage c."""
    from magestic.utils.provenance import StepProvenance, write_run_provenance

    out = Path(cfg["out_dir"])
    frag_dir = cfg.get("fragment_counts_dir") or str(out / "bc1_donor_bc0_fragment_counts")
    dbc_dir = cfg.get("donor_bc0_counts_dir") or str(out / "complete_donor_bc0_counts")
    tables = str(out / "data_tables")
    tie_order = cfg.get("tie_order", "designed")
    _kind(tie_order)
    table = os.path.join(tables, cfg.get("table_name", "bc1_reference_table.tsv"))
    key = pd.read_csv(cfg["key"], sep="\t")
    sidecars, summary = [], {}

    def step(name, output, desc):
        p = StepProvenance("bc1_reference_split_amplicon", name, output, description=desc)
        p.add_input(cfg["key"], role="sample_key")
        return p

    if "b" in stages:
        os.makedirs(frag_dir, exist_ok=True)
        d = cfg["fragment_fastq_dir"]
        jobs = [(f"{d}/{p}_trimmed_R1.fastq.gz", f"{d}/{p}_trimmed_R2.fastq.gz", f"{frag_dir}/{p}{FRAGMENT_SUFFIX}")
                for p in _samples(d, "_trimmed_R1.fastq.gz")]
        with multiprocessing.get_context("fork").Pool(max(1, n_jobs)) as pool:
            pairs = dict(pool.map(_run_b, jobs))
        p = step("b_parse_bc1_fragment_pairs", frag_dir, "per-sample R1/R2 -> bc1-donor-bc0 fragment counts")
        p.set_counts(samples=len(pairs), read_pairs=int(sum(pairs.values()))).write(f"{frag_dir}/.b.provenance.json")
        sidecars.append(f"{frag_dir}/.b.provenance.json")
        summary["b"] = {"samples": len(pairs), "read_pairs": int(sum(pairs.values()))}
    if "e" in stages:
        os.makedirs(dbc_dir, exist_ok=True)
        d = cfg["donor_bc0_fastq_dir"]
        jobs = [(f"{d}/{p}_collapsed.fastq", f"{dbc_dir}/{p}{DONOR_BC0_SUFFIX}") for p in _samples(d, "_collapsed.fastq")]
        with multiprocessing.get_context("fork").Pool(max(1, n_jobs)) as pool:
            recs = dict(pool.map(_run_e, jobs))
        p = step("e_parse_donor_bc0", dbc_dir, "per-sample collapsed merged reads -> donor-bc0 records")
        p.set_counts(samples=len(recs), records=int(sum(recs.values()))).write(f"{dbc_dir}/.e.provenance.json")
        sidecars.append(f"{dbc_dir}/.e.provenance.json")
        summary["e"] = {"samples": len(recs), "records": int(sum(recs.values()))}
    if "f" in stages:
        o = os.path.join(tables, "combined_donor_bc0_counts.tsv")
        p = step("f_combine_donor_bc0", o, "donor-bc0 counts summed over samples")
        summary["f"] = combine_donor_bc0(key, dbc_dir, tables)
        p.set_counts(**summary["f"]).add_output(o).write()
        sidecars.append(p.sidecar_path())
    if "g" in stages:
        o = os.path.join(tables, "donor_bc0_to_oligo_mapping.tsv")
        p = step("g_map_donor_bc0_to_oligos", o, "full donor -> designed oligo")
        p.add_input(cfg["oligo_pool"], role="oligo_pool")
        p.set_parameters(tie_order=tie_order)
        summary["g"] = map_donor_bc0_to_oligos(os.path.join(tables, "combined_donor_bc0_counts.tsv"),
                                               cfg["oligo_pool"], tables, tie_order)
        p.set_counts(**summary["g"]).add_output(o).write()
        sidecars.append(p.sidecar_path())
    if "c" in stages:
        o = os.path.join(tables, "top_bc1_donor_bc0_fragments.tsv")
        p = step("c_combine_bc1_fragments", o, "edit-distance collapse, PCR-replicate support, bc1 purity")
        p.set_parameters(edit_distance_threshold=EDIT_DISTANCE_THRESHOLD, min_pcr_replicates=MIN_PCR_REPLICATES,
                         tie_order=tie_order)
        rank = None
        if tie_order == "designed":
            mapping = os.path.join(tables, "donor_bc0_to_oligo_mapping.tsv")
            p.add_input(mapping, role="donor_bc0_map").add_input(cfg["gdb_table"], role="gdb_plasmid_map")
            rank = designed_fragment_rank(mapping, cfg["gdb_table"])
        summary["c"] = combine_bc1_fragments(key, frag_dir, tables, n_jobs, tie_order, rank)
        p.set_counts(**summary["c"]).add_output(o).write()
        sidecars.append(p.sidecar_path())
    if "h" in stages:
        p = step("h_link_reference_table", table, "bc1 reference table")
        p.add_input(cfg["gdb_table"], role="gdb_plasmid_map").add_input(cfg["oligo_info"], role="oligo_info")
        p.set_parameters(min_bc1_purity=cfg.get("min_purity", MIN_BC1_PURITY))
        summary["h"] = link_reference_table(os.path.join(tables, "top_bc1_donor_bc0_fragments.tsv"),
                                            os.path.join(tables, "donor_bc0_to_oligo_mapping.tsv"),
                                            cfg["gdb_table"], cfg["oligo_info"], table,
                                            cfg.get("min_purity", MIN_BC1_PURITY))
        p.set_counts(**summary["h"]).add_output(table, rows=summary["h"]["rows_out"]).write()
        sidecars.append(p.sidecar_path())
    if sidecars:
        write_run_provenance(str(out), pipeline="bc1_reference_split_amplicon", run_name=out.name,
                             step_sidecars=sidecars, description="split-amplicon bc1 reference table")
    return summary


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(prog="magestic-bc1-reference",
                                 description="Build a bc1 reference table from a split-amplicon (2x150) layout: "
                                             "bc1-donor-bc0 fragment reads + donor-bc0 reads + a GDB plasmid map.")
    ap.add_argument("--config", required=True, help="JSON or YAML with the keys documented in run()")
    ap.add_argument("--stages", default="b,e,f,g,c,h")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    text = Path(a.config).read_text()
    if a.config.endswith((".yaml", ".yml")):
        import yaml
        cfg = yaml.safe_load(text)
    else:
        cfg = json.loads(text)
    t = time.time()
    summary = run(cfg, [s.strip() for s in a.stages.split(",") if s.strip()], a.jobs)
    print(json.dumps({"summary": summary, "seconds": round(time.time() - t)}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
