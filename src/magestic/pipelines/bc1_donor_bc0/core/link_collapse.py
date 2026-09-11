"""Link aggregation and collapsing for bc1-donor-bc0 tables, with EXACT replicate sets.

Written 2026-09-10 for the 20210822 rebuild (planning ``28_``). Three rules this module
exists to enforce, each of which the previous code broke:

1. **Replicate support is a SET, unioned -- never a sum.** A link seen in replicates
   {p1, p2} merged with a sequencing-error variant seen in {p1} has support {p1, p2}
   = 2, not 3. Summing ``num_PCR_replicates`` across collapsed rows (the old step-02
   Step 4 and step 08) inflates exactly the singleton tier Kevin asked to see.
   Sets are carried as fixed-width uint64 BITMASKS (one bit per replicate id), so
   union is a vectorized bitwise OR and support is a popcount.

2. **The collapse representative is chosen by read count alone.** The old greedy
   peel-off ranked by the sum of *every* count column; add per-tier or bitmask
   columns and the ranking silently changes. ``rank_col`` is explicit here.

3. **A shorter donor is the same link as the longer donor it is a suffix of.**
   2x150 reads yield ~88 nt of donor, 2x300 reads 238 nt; measured on 20210822, the
   2021 donor is an exact suffix of the AVITI donor for 98.3% of dual-era bc1s. An
   edit-distance collapse can never merge them, so without :func:`attach_short_donors`
   every such link is split in two across eras (1,002,360 bc1s in the 2026-03-22 table).
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from functools import partial
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from rapidfuzz.distance import Levenshtein as _lev

# ---------------------------------------------------------------------------
# Replicate bitmasks
# ---------------------------------------------------------------------------

def mask_columns(n_ids: int) -> List[str]:
    """Column names for an ``n_ids``-bit replicate mask (64 ids per uint64 word)."""
    return [f"_repmask_{w}" for w in range((max(n_ids, 1) + 63) // 64)]


def masks_for_codes(codes: np.ndarray, n_ids: int) -> np.ndarray:
    """One-hot bitmask rows (shape ``[len(codes), n_words]``, uint64) for replicate codes."""
    n_words = len(mask_columns(n_ids))
    codes = np.asarray(codes, dtype=np.int64)
    if (codes < 0).any() or (codes >= n_ids).any():
        raise ValueError("replicate code out of range")
    out = np.zeros((len(codes), n_words), dtype=np.uint64)
    word = codes // 64
    bit = (codes % 64).astype(np.uint64)
    out[np.arange(len(codes)), word] = np.left_shift(np.uint64(1), bit)
    return out


_POP8 = np.array([bin(i).count("1") for i in range(256)], dtype=np.uint8)


def popcount_rows(masks: np.ndarray) -> np.ndarray:
    """Number of set bits per row of a uint64 mask array -> int64."""
    if masks.size == 0:
        return np.zeros(masks.shape[0], dtype=np.int64)
    as_bytes = np.ascontiguousarray(masks).view(np.uint8).reshape(masks.shape[0], -1)
    return _POP8[as_bytes].sum(axis=1).astype(np.int64)


def decode_masks(masks: np.ndarray, id_names: Sequence[str], sep: str = ";") -> List[str]:
    """Sorted, ``sep``-joined replicate-id names per row (for the output table)."""
    names = list(id_names)
    out = []
    for row in masks:
        ids = []
        for w, word in enumerate(row):
            word = int(word)
            while word:
                low = word & -word
                ids.append(names[w * 64 + low.bit_length() - 1])
                word ^= low
        out.append(sep.join(sorted(ids)))
    return out


# ---------------------------------------------------------------------------
# Vectorized aggregation by key (associative: safe to apply batch-wise)
# ---------------------------------------------------------------------------

def aggregate_links(df: pd.DataFrame, key_cols: Sequence[str], rank_col: str,
                    sum_cols: Sequence[str], or_cols: Sequence[str],
                    preserve_cols: Sequence[str] = ()) -> pd.DataFrame:
    """Group rows by ``key_cols``: SUM ``sum_cols`` (``rank_col`` must be among them),
    bitwise-OR ``or_cols``, and take ``preserve_cols`` from the row with the highest
    ``rank_col`` in each group (ties -> the lexicographically smallest preserve tuple,
    so the answer does not depend on how rows were batched).

    Associative and commutative in the aggregated columns, so it can be applied to
    batches and then to the concatenation of batch results with the same answer.
    Rows with a NaN in any key column are NOT accepted -- callers must set them
    aside explicitly (pandas would otherwise drop them silently in the groupby).
    """
    if len(df) == 0:
        return df.copy()
    key_cols = list(key_cols)
    if df[key_cols].isna().any().any():
        raise ValueError("aggregate_links: NaN in key columns; set those rows aside first")
    if len(sum_cols) and df[list(sum_cols)].isna().any().any():
        # NaN cast to int64 is INT64_MIN -- a silent wrong number (review F4/F5, 2026-09-10)
        raise ValueError(f"aggregate_links: NaN in count columns {[c for c in sum_cols if df[c].isna().any()]}")
    codes, _ = pd.factorize(pd.MultiIndex.from_frame(df[key_cols]), sort=False)
    order = np.argsort(codes, kind="stable")
    codes_s = codes[order]
    starts = np.flatnonzero(np.r_[True, codes_s[1:] != codes_s[:-1]])

    out = pd.DataFrame(index=np.arange(len(starts)))
    # key values from each group's first row (version-independent; factorize on a
    # MultiIndex does not preserve level names in every pandas release)
    first_rows = order[starts]
    for c in key_cols:
        out[c] = df[c].to_numpy(dtype=object)[first_rows]
    for c in sum_cols:
        vals = df[c].to_numpy()[order]
        out[c] = np.add.reduceat(vals.astype(np.int64), starts)
    for c in or_cols:
        vals = df[c].to_numpy(dtype=np.uint64)[order]
        out[c] = np.bitwise_or.reduceat(vals, starts)
    if preserve_cols:
        # representative row = max rank within group (deterministic tie-break below)
        rank = df[rank_col].to_numpy()[order]
        # position (within sorted order) of the max-rank row in each group
        seg_max = np.maximum.reduceat(rank, starts)
        seg_id = np.repeat(np.arange(len(starts)), np.diff(np.r_[starts, len(order)]))
        idx = np.flatnonzero(rank == seg_max[seg_id])
        # Among max-rank rows, take the lexicographically smallest preserve tuple --
        # NOT "first in input order", which depends on how rows were batched and would
        # make batch folding non-associative for these columns (review F4, 2026-09-10).
        cand = pd.DataFrame({"_seg": seg_id[idx], "_row": order[idx]})
        for c in preserve_cols:
            cand[c] = df[c].to_numpy(dtype=object)[order[idx]]
            cand[c] = cand[c].map(lambda v: "" if v is None or (isinstance(v, float) and v != v) else str(v))
        cand = cand.sort_values(["_seg", *preserve_cols], kind="mergesort").drop_duplicates("_seg")
        src = cand["_row"].to_numpy()
        for c in preserve_cols:
            out[c] = df[c].to_numpy(dtype=object)[src]
    return out[[*key_cols, *sum_cols, *or_cols, *preserve_cols]]


# ---------------------------------------------------------------------------
# Greedy edit-distance collapse (behavior of the 2026-06-02 _peel_off, made exact)
# ---------------------------------------------------------------------------

def _peel(seqs, ranks, sums, ors, threshold):
    """Greedy top-by-RANK peel-off. Returns [(rep_idx, summed_vec, or_vec), ...]."""
    n = len(seqs)
    alive = [True] * n
    out = []
    remaining = n
    ns, no = (len(sums[0]) if n else 0), (len(ors[0]) if n and ors else 0)
    while remaining:
        top, best = -1, None
        for i in range(n):
            if alive[i] and (best is None or ranks[i] > best):
                best, top = ranks[i], i
        agg = list(sums[top])
        oagg = list(ors[top]) if no else []
        alive[top] = False
        remaining -= 1
        ts = seqs[top]
        if isinstance(ts, str):
            for i in range(n):
                if not alive[i]:
                    continue
                si = seqs[i]
                if isinstance(si, str) and (si == ts or
                                            _lev.distance(si, ts, score_cutoff=threshold) <= threshold):
                    sv = sums[i]
                    for k in range(ns):
                        agg[k] += sv[k]
                    if no:
                        ov = ors[i]
                        for k in range(no):
                            oagg[k] |= ov[k]
                    alive[i] = False
                    remaining -= 1
        out.append((top, agg, oagg))
    return out


def _collapse_payload(payload, group_cols, seq_col, sum_cols, or_cols, pcols, threshold):
    gkey, seqs, ranks, sums, ors, pres = payload
    base = {group_cols[k]: gkey[k] for k in range(len(group_cols))}
    rows = []
    for rep, agg, oagg in _peel(seqs, ranks, sums, ors, threshold):
        row = dict(base)
        row[seq_col] = seqs[rep]
        for k, c in enumerate(sum_cols):
            row[c] = agg[k]
        for k, c in enumerate(or_cols):
            row[c] = oagg[k]
        if pcols:
            for k, c in enumerate(pcols):
                row[c] = pres[rep][k]
        rows.append(row)
    return rows


def collapse_links(df: pd.DataFrame, group_cols: Sequence[str], seq_col: str,
                   rank_col: str, sum_cols: Sequence[str], or_cols: Sequence[str] = (),
                   preserve_cols: Sequence[str] = (), threshold: int = 1,
                   n_jobs: int = 1) -> pd.DataFrame:
    """Within each ``group_cols`` group, greedily merge ``seq_col`` values within
    ``threshold`` edits of the highest-``rank_col`` survivor: SUM ``sum_cols``,
    OR ``or_cols``, keep ``preserve_cols`` (and ``seq_col``) from the representative.

    Row cardinality only ever DECREASES by merging; nothing is filtered.
    """
    group_cols, sum_cols, or_cols = list(group_cols), list(sum_cols), list(or_cols)
    pcols = [c for c in preserve_cols if c in df.columns]
    out_cols = [*group_cols, seq_col, *sum_cols, *or_cols, *pcols]
    if rank_col not in sum_cols:
        raise ValueError("rank_col must be one of sum_cols")
    if len(df) and df[sum_cols].isna().any().any():
        raise ValueError(f"collapse_links: NaN in count columns {[c for c in sum_cols if df[c].isna().any()]}")
    if len(df) == 0:
        return pd.DataFrame(columns=out_cols)
    ri = sum_cols.index(rank_col)

    codes, _ = pd.factorize(pd.MultiIndex.from_frame(df[group_cols]), sort=False)
    order = np.argsort(codes, kind="stable")
    codes_s = codes[order]
    n = len(df)
    bounds = np.flatnonzero(np.diff(codes_s)) + 1
    starts = np.r_[0, bounds].tolist()
    ends = np.r_[bounds, n].tolist()

    gv = df[group_cols].to_numpy(dtype=object)[order]
    sq = df[seq_col].to_numpy(dtype=object)[order]
    sv = df[sum_cols].to_numpy(dtype=np.int64)[order]
    ov = df[or_cols].to_numpy(dtype=np.uint64)[order] if or_cols else None
    pv = df[pcols].to_numpy(dtype=object)[order] if pcols else None

    singles_idx = []
    payloads = []
    for a, b in zip(starts, ends):
        if b - a == 1:
            singles_idx.append(a)
        else:
            payloads.append((
                list(gv[a]), list(sq[a:b]),
                [int(sv[i][ri]) for i in range(a, b)],
                [[int(x) for x in sv[i]] for i in range(a, b)],
                [[int(x) for x in ov[i]] for i in range(a, b)] if or_cols else [[] for _ in range(a, b)],
                [list(pv[i]) for i in range(a, b)] if pcols else None,
            ))

    # singletons: vectorized passthrough
    si = np.asarray(singles_idx, dtype=np.int64)
    single = pd.DataFrame({c: gv[si, k] for k, c in enumerate(group_cols)})
    single[seq_col] = sq[si]
    for k, c in enumerate(sum_cols):
        single[c] = sv[si, k]
    for k, c in enumerate(or_cols):
        single[c] = ov[si, k]
    for k, c in enumerate(pcols):
        single[c] = pv[si, k]

    worker = partial(_collapse_payload, group_cols=group_cols, seq_col=seq_col,
                     sum_cols=sum_cols, or_cols=or_cols, pcols=pcols, threshold=threshold)
    multi = []
    if payloads:
        if n_jobs and n_jobs > 1 and len(payloads) > 1:
            chunksize = max(1, len(payloads) // (n_jobs * 8))
            with ProcessPoolExecutor(max_workers=n_jobs) as ex:
                for rows in ex.map(worker, payloads, chunksize=chunksize):
                    multi.extend(rows)
        else:
            for p in payloads:
                multi.extend(worker(p))
    if not multi:
        return single[out_cols]
    # Build the merged rows' numeric columns as EXPLICIT arrays. Letting pandas infer
    # them from Python ints makes a mask word int64 when every value is < 2**63; the
    # concat with the uint64 singletons then promotes to float64, which silently drops
    # low bits whenever a word also carries a bit >= 2**53 (reproduced 2026-09-10 on
    # pandas 1.5.1 and 2.2.2: a link seen in replicates {0, 60} reported support 1).
    mdf = pd.DataFrame({c: [r[c] for r in multi] for c in [*group_cols, seq_col, *pcols]})
    for c in sum_cols:
        mdf[c] = np.array([r[c] for r in multi], dtype=np.int64)
    for c in or_cols:
        mdf[c] = np.array([r[c] for r in multi], dtype=np.uint64)
    res = pd.concat([single[out_cols], mdf[out_cols]], ignore_index=True)
    for c in or_cols:
        if res[c].dtype != np.uint64:
            raise AssertionError(f"mask column {c} lost its uint64 dtype ({res[c].dtype})")
    return res


# ---------------------------------------------------------------------------
# Cross-read-length reconciliation
# ---------------------------------------------------------------------------

def _hamming_le(a: str, b: str, k: int) -> bool:
    if len(a) != len(b):
        return False
    d = 0
    for x, y in zip(a, b):
        if x != y:
            d += 1
            if d > k:
                return False
    return True


def attach_short_donors(df: pd.DataFrame, seq_col: str, rank_col: str,
                        sum_cols: Sequence[str], or_cols: Sequence[str],
                        max_mismatch: int = 5, bc0_max_edit: int = 1,
                        min_attach_len: int = 30,
                        dominance: float = 3.0) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """Within each bc1, merge a shorter donor into the UNIQUE longer donor of which it
    is a suffix (<= ``max_mismatch`` Hamming mismatches over the aligned suffix) and
    whose bc0 is within ``bc0_max_edit`` edits. Donors are processed longest-first, so
    an 88-mer whose 150-mer parent was itself absorbed into a 238-mer sees one
    candidate, not two.

    On attach: ``sum_cols`` summed, ``or_cols`` OR-ed, the long donor kept.
    0 candidates -> the short link stays as its own row. >1 -> attach to the top
    candidate IF it carries >= ``dominance`` x the runner-up's ``rank_col``
    (posterior >= d/(d+1) that the short read came from it, under abundance-
    proportional origin), else it stays, flagged ``donor_suffix_ambiguous``.
    Measured on the 20210822 smoke (29_): a strict any->1-is-ambiguous rule left
    18.2% of reads unreconciled, the runner-up was single-replicate 82% of the
    time, and the top had >=3x the runner-up's reads 57% of the time.
    ``dominance=float('inf')`` restores the strict rule. Donors shorter than
    ``min_attach_len`` never attach.
    Adds ``donor_suffix_ambiguous`` (bool) and ``donor_lengths_observed`` (str).
    Returns (frame, stats).
    """
    sum_cols, or_cols = list(sum_cols), list(or_cols)
    stats = {"groups_examined": 0, "attached": 0, "attached_by_dominance": 0,
             "ambiguous": 0, "no_parent": 0}
    df = df.reset_index(drop=True).copy()
    lens = df[seq_col].map(lambda s: len(s) if isinstance(s, str) else 0).to_numpy()
    df["donor_suffix_ambiguous"] = False
    df["donor_lengths_observed"] = [str(x) for x in lens]
    if len(df) == 0:
        return df, stats

    codes, _ = pd.factorize(df["bc1"], sort=False)
    order = np.argsort(codes, kind="stable")
    cs = codes[order]
    bounds = np.flatnonzero(np.diff(cs)) + 1
    drop = np.zeros(len(df), dtype=bool)
    seqs = df[seq_col].to_numpy(dtype=object)
    bc0s = df["bc0"].to_numpy(dtype=object)
    ranks = df[rank_col].to_numpy()
    sums = df[sum_cols].to_numpy(dtype=np.int64)
    ors = df[or_cols].to_numpy(dtype=np.uint64) if or_cols else None
    amb = df["donor_suffix_ambiguous"].to_numpy()
    lobs = df["donor_lengths_observed"].to_numpy(dtype=object)

    for a, b in zip(np.r_[0, bounds], np.r_[bounds, len(df)]):
        if b - a < 2:
            continue
        idx = order[a:b]
        glens = lens[idx]
        if glens.max() == glens.min():
            continue
        stats["groups_examined"] += 1
        # longest first, then most reads
        idx = sorted(idx, key=lambda i: (-lens[i], -ranks[i]))
        survivors = []
        for i in idx:
            s = seqs[i]
            L = lens[i]
            cands = []
            if L >= min_attach_len:
                for j in survivors:
                    if lens[j] > L and _hamming_le(s, seqs[j][-L:], max_mismatch) and \
                            isinstance(bc0s[i], str) and isinstance(bc0s[j], str) and \
                            _lev.distance(bc0s[i], bc0s[j], score_cutoff=bc0_max_edit) <= bc0_max_edit:
                        cands.append(j)
            chosen = None
            if len(cands) == 1:
                chosen = cands[0]
            elif len(cands) > 1:
                top2 = sorted(cands, key=lambda c: -ranks[c])[:2]
                if ranks[top2[0]] >= dominance * max(ranks[top2[1]], 1):
                    chosen = top2[0]
                    stats["attached_by_dominance"] += 1
            if chosen is not None:
                j = chosen
                sums[j] += sums[i]
                if or_cols:
                    ors[j] |= ors[i]
                lobs[j] = ";".join(sorted(set(lobs[j].split(";")) | set(lobs[i].split(";")),
                                          key=int, reverse=True))
                drop[i] = True
                stats["attached"] += 1
            else:
                if len(cands) > 1:
                    amb[i] = True
                    stats["ambiguous"] += 1
                elif lens[i] < glens.max():
                    stats["no_parent"] += 1
                survivors.append(i)

    df[sum_cols] = sums
    if or_cols:
        for k, c in enumerate(or_cols):
            df[c] = ors[:, k]
    df["donor_suffix_ambiguous"] = amb
    df["donor_lengths_observed"] = lobs
    return df.loc[~drop].reset_index(drop=True), stats


# ---------------------------------------------------------------------------
# Orphan bc1 merge (the X1 defect)
# ---------------------------------------------------------------------------
#
# WHY THIS EXISTS. The three collapses above merge a bc1 sequencing error into its parent only when the
# other two fields match EXACTLY (they group by (donor, bc0) and peel similar bc1 within the group). A
# 238-nt donor read rarely repeats exactly, so a bc1 error whose donor read also carries an error is
# never grouped with its parent and survives as its own "orphan" barcode. Evidence that the minor member
# really is a read error of the major: planning 31b_x1_mechanism.py.
#
# Kevin, 2026-09-10: "can we merge them back into their parents for count purposes?" -- yes, and this is
# where it belongs, beside the other collapses, rather than in a post-hoc script run against the shipped
# reference table (PIPELINE_CORRECTNESS_POLICY.md §1: a repair applied after the fact protects one
# consumer and leaves the split barcodes on disk for everyone else).
#
# 🔴 THE STAR-ONLY GUARD. Accepted pairs form a GRAPH: a barcode can be the minor of one pair and the
# major of another, so naive connected components merge transitively and produce components whose
# members are up to Hamming distance 8 apart -- far beyond anything a single read error explains. Every
# merged member must therefore have DIRECT pairwise evidence against the component root (a "star", not a
# chain). Measured on the 2026-09-10 reference table, that keeps 91.5% of merges and removes the chains.
#
# WHAT IS AND IS NOT CHANGED. Counts and the per-tier count columns are SUMMED; replicate/run/library
# membership is bitwise-ORed, i.e. a set UNION, never a sum -- re-sequencing the same PCR replicate is
# not new evidence, and support_tier is recomputed from the union afterwards. Measured effect on the
# 2026-09-10 table: 1,188,500 accepted pairs remove 1,075,901 barcodes (20.2% of them) but move only
# 1,308,605 reads (0.34%). It is a barcode-list correction, not a count correction.


def _bc1_hamming1_pairs(bc1: np.ndarray) -> np.ndarray:
    """Every pair of indices whose 20-nt bc1 differ by exactly one SUBSTITUTION.

    Wildcard-key method: a 20-mer has 20 one-position-masked keys, and two sequences are Hamming-1 iff
    they share a masked key. Sorting 20N keys finds every pair in O(N log N) instead of O(N^2).
    Requires clean 20-nt ACGT barcodes; callers filter first.
    """
    n = len(bc1)
    if n < 2:
        return np.empty((0, 2), dtype=np.int64)
    lut = np.zeros(256, dtype=np.uint64)
    lut[[ord(c) for c in "ACGT"]] = np.arange(4, dtype=np.uint64)
    A = lut[np.frombuffer("".join(bc1).encode(), dtype=np.uint8).reshape(n, 20)]
    code = np.zeros(n, dtype=np.uint64)
    for j in range(20):
        code |= A[:, j] << np.uint64(2 * j)
    keys = np.concatenate([(code & ~np.uint64(3 << (2 * j))) | np.uint64(j << 42) for j in range(20)])
    idx = np.tile(np.arange(n, dtype=np.int64), 20)
    order = np.argsort(keys, kind="stable")
    keys, idx = keys[order], idx[order]
    brk = np.flatnonzero(np.diff(keys)) + 1
    out = []
    for s, e in zip(np.r_[0, brk], np.r_[brk, len(keys)]):
        if e - s < 2:
            continue
        m = idx[s:e]
        for u in range(len(m)):
            for v in range(u + 1, len(m)):
                a, b = m[u], m[v]
                out.append((a, b) if a < b else (b, a))
    if not out:
        return np.empty((0, 2), dtype=np.int64)
    return np.unique(np.array(out, dtype=np.int64), axis=0)


def _edit_le1(a: str, b: str) -> bool:
    """True when a and b are within one substitution, insertion or deletion."""
    if a == b:
        return True
    la, lb = len(a), len(b)
    if la == lb:
        d = 0
        for x, y in zip(a, b):
            if x != y:
                d += 1
                if d > 1:
                    return False
        return d == 1
    if abs(la - lb) != 1:
        return False
    if la > lb:
        a, b = b, a
    i = 0
    while i < len(a) and a[i] == b[i]:
        i += 1
    return a[i:] == b[i + 1:]


def _donor_compatible(x, y, max_mismatch: int) -> bool:
    """Two donor reads can be the same molecule: identical, one a suffix of the other (a 2x150 read of a
    2x300 parent), or the same length within ``max_mismatch`` substitutions."""
    if not (isinstance(x, str) and isinstance(y, str)) or not x or not y:
        return False
    if x == y:
        return True
    s, l = (x, y) if len(x) <= len(y) else (y, x)
    if len(s) != len(l):
        return l.endswith(s)
    d = 0
    for p, q in zip(s, l):
        if p != q:
            d += 1
            if d > max_mismatch:
                return False
    return True


def merge_orphan_bc1(df: pd.DataFrame, rank_col: str, sum_cols: Sequence[str],
                     or_cols: Sequence[str], preserve_cols: Sequence[str] = (),
                     bc1_col: str = "bc1", donor_col: str = "donor", bc0_col: str = "bc0",
                     donor_max_mismatch: int = 5, star_guard: bool = True):
    """Merge orphan bc1 into their parents, then re-aggregate. Returns ``(df, stats, bc1_to_root)``.

    A candidate pair is two bc1 at Hamming 1. It is ACCEPTED when **some** link of one and **some** link
    of the other have bc0 within one edit and compatible donors -- every link of each barcode is tried,
    not just its most abundant one, because an orphan's dominant link is not always the one carrying the
    evidence. The accepted root link (highest-count matching row of the root) is remembered per pair.

    The component root is the bc1 with the highest total ``rank_col``; ties break lexicographically so
    reruns are reproducible. With ``star_guard`` (the default) a member joins only if the (member, root)
    pair was itself accepted, so every merge rests on direct evidence rather than a transitive chain --
    naive connected components produced components spanning Hamming distance 8.

    Re-keying the bc1 alone does not merge the link: the premise of an orphan is that its donor AND bc0
    reads also carry errors, which is exactly why the exact-key collapses missed it, so after re-keying
    the rows still differ and ``aggregate_links`` keeps them apart. Each row of a re-keyed bc1 is
    therefore snapped onto the accepted root link when it is individually compatible with it, by the same
    two tests that accepted the pair. Any other link of that barcode keeps its own donor and bc0.

    ``bc1_to_root`` is returned so downstream tables keyed on bc1 -- the screen counts matrix above all --
    can be re-keyed identically instead of re-deriving the merge and possibly disagreeing.
    """
    stats: Dict[str, object] = {"rows_in": int(len(df)), "bc1_in": int(df[bc1_col].nunique()),
                                "reads_in": int(df[rank_col].sum()), "star_guard": bool(star_guard),
                                "donor_max_mismatch": int(donor_max_mismatch)}
    empty_map = pd.DataFrame({bc1_col: pd.Series(dtype=str), "bc1_root": pd.Series(dtype=str)})

    def _unchanged(extra=None):
        stats.update(accepted_pairs=0, components_merged=0, bc1_removed=0, reads_moved=0,
                     rows_snapped_to_root_link=0, rows_out=int(len(df)),
                     bc1_out=int(df[bc1_col].nunique()) if len(df) else 0,
                     reads_out=int(df[rank_col].sum()) if len(df) else 0)
        if extra:
            stats.update(extra)
        return df, stats, empty_map

    if df.empty:
        return _unchanged()

    clean_mask = df[bc1_col].astype(str).str.fullmatch("[ACGT]{20}").fillna(False).to_numpy()
    sub = df[clean_mask]
    bc1v = np.array(sorted(sub[bc1_col].unique()))
    stats["bc1_clean_20nt"] = int(len(bc1v))
    if len(bc1v) < 2:
        return _unchanged()

    pos = {b: i for i, b in enumerate(bc1v)}
    # per-bc1 list of (donor, bc0, counts), most abundant first
    links: List[List[Tuple[object, object, int]]] = [[] for _ in bc1v]
    srt = sub.sort_values(rank_col, ascending=False)
    for b, d, c0, ct in zip(srt[bc1_col].to_numpy(), srt[donor_col].to_numpy(),
                            srt[bc0_col].to_numpy(), srt[rank_col].to_numpy()):
        links[pos[b]].append((d, c0, int(ct)))
    totv = np.array([sum(x[2] for x in L) for L in links], dtype=np.int64)

    pairs = _bc1_hamming1_pairs(bc1v)
    stats["hamming1_bc1_pairs"] = int(len(pairs))

    accepted: Dict[Tuple[int, int], None] = {}
    best_link: Dict[Tuple[int, int], Tuple[object, object]] = {}   # (member, root) -> root's matched link
    for a, b in pairs:
        a, b = int(a), int(b)
        hit = None
        for da, ba, _ca in links[a]:
            for db, bb, _cb in links[b]:
                if (isinstance(ba, str) and isinstance(bb, str) and _edit_le1(ba, bb)
                        and _donor_compatible(da, db, donor_max_mismatch)):
                    hit = ((da, ba), (db, bb))
                    break
            if hit:
                break
        if hit is None:
            continue
        accepted[(a, b)] = None
        best_link[(a, b)] = hit[1]      # b's link, used when b is the root
        best_link[(b, a)] = hit[0]      # a's link, used when a is the root
    stats["accepted_pairs"] = len(accepted)
    if not accepted:
        return _unchanged()

    parent = np.arange(len(bc1v), dtype=np.int64)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in accepted:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    comp = np.array([find(i) for i in range(len(bc1v))], dtype=np.int64)

    order = np.lexsort((bc1v, -totv))
    root_of_comp: Dict[int, int] = {}
    for i in order:
        c = int(comp[i])
        if c not in root_of_comp:
            root_of_comp[c] = int(i)

    remap: Dict[str, str] = {}
    root_link: Dict[str, Tuple[object, object]] = {}
    dropped_by_guard = 0
    for i in range(len(bc1v)):
        r = root_of_comp[int(comp[i])]
        if i == r:
            continue
        key = (i, r) if (i, r) in best_link else ((r, i) if (r, i) in best_link else None)
        direct = ((min(i, r), max(i, r)) in accepted)
        if star_guard and not direct:
            dropped_by_guard += 1
            continue
        remap[bc1v[i]] = bc1v[r]
        if direct:
            root_link[bc1v[i]] = best_link[(i, r)]
        else:                       # star_guard off: no direct evidence, so do not snap the link
            root_link[bc1v[i]] = (None, None)
    stats["members_dropped_by_star_guard"] = int(dropped_by_guard)
    stats["bc1_remapped"] = int(len(remap))

    bc1_to_root = pd.DataFrame({bc1_col: list(remap), "bc1_root": [remap[b] for b in remap]})
    if bc1_to_root.empty:
        return _unchanged()

    totals = pd.Series(totv, index=bc1v)
    reads_moved = int(totals.reindex(bc1_to_root[bc1_col]).fillna(0).sum())

    out = df.copy()
    bc1_arr = out[bc1_col].to_numpy()
    aff = np.array([b in remap for b in bc1_arr])
    snapped = 0
    if aff.any():
        d_new = out[donor_col].to_numpy(copy=True)
        b_new = out[bc0_col].to_numpy(copy=True)
        for i in np.flatnonzero(aff):
            rd, rb = root_link[bc1_arr[i]]
            if rd is None:
                continue
            if (_donor_compatible(d_new[i], rd, donor_max_mismatch)
                    and isinstance(b_new[i], str) and isinstance(rb, str) and _edit_le1(b_new[i], rb)):
                d_new[i], b_new[i] = rd, rb
                snapped += 1
        out[donor_col] = d_new
        out[bc0_col] = b_new
    stats["rows_snapped_to_root_link"] = int(snapped)
    out[bc1_col] = out[bc1_col].map(lambda x: remap.get(x, x))
    out = aggregate_links(out, [bc1_col, donor_col, bc0_col], rank_col, sum_cols, or_cols, preserve_cols)

    stats.update(
        components_merged=int(len({root_of_comp[int(comp[pos[b]])] for b in remap})),
        bc1_removed=int(len(remap)),
        reads_moved=reads_moved,
        rows_out=int(len(out)), bc1_out=int(out[bc1_col].nunique()),
        reads_out=int(out[rank_col].sum()),
    )
    if stats["reads_out"] != stats["reads_in"]:   # explicit: survives python -O
        raise RuntimeError(f"orphan merge lost reads: in {stats['reads_in']} != out {stats['reads_out']}")
    return out, stats, bc1_to_root
