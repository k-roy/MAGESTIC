"""
Pooled-slope fitness: the method behind the NRD1 fitness screen (Supplementary Table 1, Figure 1D).

A Python implementation of ``NRD1/fitness_slope_heatmap.R`` (Umberto Aiello, 2025-08-25) plus the two steps
around it that the published numbers also depend on: the consolidation of per-set normalized counts
(``planning/289_consolidate_umberto.py``: pseudocount, fill, the unique-donor filter) and criterion 4 of
``planning/333c_build_supp_table1.py`` (a hit needs at least one barcode that meets criteria 1 to 3 alone).
It sits beside the per-barcode regression of step 12 (``core.slope_fitness``); it does not replace it.

THE METHOD, per edit (a group of barcodes sharing ``key_cols``) and treatment:

1. A barcode is kept when its mean normalized count over every generation-0 sample is at least ``min_g0``;
   generation-0 values are replaced by that mean, and each value becomes log2(max(x, eps) / max(g0, eps)).
2. The edit's series is the mean log2 fold change across its barcodes and replicates at each generation,
   weighted by max(g0, 1).
3. The series is truncated at the first run of ``trim_k`` consecutive post-generation-0 points below
   log2(``trim_fc``) (inclusive), and the raw slope is the least-squares slope over the kept points
   (at least ``min_points``).
3b. ``index="endpoint"`` replaces the slope with the pooled log2 fold change at ONE generation (the series'
   last finite post-zero generation, or ``endpoint_generation``), so an edit measured at generation 0 and one
   later point gets an index where the slope route gets nothing. Truncation is off by default on this route,
   because a series that dies IS the depletion an endpoint measure reports. ``min_points`` does not apply.
   Steps 4 to 7 are shared verbatim between the two routes, so they cannot drift apart.
4. The fitness index is the slope minus the median slope of all edits in the same treatment.
5. Depletion: the series' last finite point is below log2(``depletion_fc``).
6. Hit (criteria 1-3): fitness index in ``treatment`` below its ``tail_prob`` quantile over all edits,
   depletion in ``treatment``, and index(treatment) / ``ratio_factor`` < index(``control``).
7. Barcode support: each barcode's own slope over its edit's window (one value per generation: the FIRST
   replicate in input order, as R's ``distinct()`` keeps), centered on the median barcode slope of the
   treatment, its own last-point depletion, and the same three tests. Criterion 4: a hit needs at least
   ``min_support`` supporting barcodes.

PARITY (planning/380): on the shipped NRD1 input this reproduces the R script's six tables and the 23 hits.
Behaviours kept from the R script on purpose, because the published numbers rest on them:
* the hit quantile and the medians are taken over ALL edits of the treatment, including rows whose mutation
  type is missing; only the hit table is restricted to ``mutation_types``;
* per-barcode slopes use the first replicate at each generation (input column order), not their mean;
* an edit with no generation-0 point is not trimmed (R), but a dataset with no generation 0 at all is
  REFUSED here, since every fold change is relative to generation 0.

NOT DECIDED HERE (planning/380 §4, Kevin): the edit key for QTL libraries, a hit rule for conditions without
an auxin-style control, enrichment (beneficial) hits, and the reference generation of screens whose key
uses negative generations. ``control=None`` therefore returns fitness indices with NO hit calls.

Author: Kevin R. Roy
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

NORM_TAG = "_Normalized_Count_"


@dataclass
class PooledSlopeConfig:
    treatment: str = "YPD"
    control: Optional[str] = "YPD_noIAA"
    key_cols: Tuple[str, ...] = ("amino_acid_position", "mutation_type")
    min_g0: float = 10.0
    eps: float = 1e-8
    weighted: bool = True
    trim_fc: float = 0.01
    trim_k: int = 3
    min_points: int = 3
    #: the R script's defaults; the NRD1 Supplementary Table 1 used 0.05 and 0.1 (planning/373)
    tail_prob: float = 0.02
    depletion_fc: float = 1e-4
    ratio_factor: float = 3.0
    min_support: int = 1
    mutation_types: Tuple[str, ...] = ("missense", "synonymous", "PTC", "Deletion")
    #: "depleted" (the R script) or "both": also call the enriched tail by the mirror of every test --
    #: index above the (1 - tail_prob) quantile, final fold change above 1 / depletion_fc, and
    #: index(treatment) / ratio_factor > index(control). PROVISIONAL for QTL screens (planning/384 §0).
    direction: str = "depleted"
    #: "slope" (the R script: OLS over the truncated series, needs >= min_points generations) or
    #: "endpoint" (the pooled log2 fold change at ONE generation relative to generation 0, needs only that
    #: generation). Kevin asked for the endpoint route on 2026-11-03 because two QTL screens yield no slope
    #: at all: the December 2021 screen has none of its 32 conditions fittable and 42 of 48 conditions in
    #: yL437/yL442 are endpoint only (planning/384 §3). Every later step is shared between the two routes.
    index: str = "slope"
    #: the generation the endpoint index is read at; None takes each series' last finite generation, so
    #: series of unequal length are each read at their own last point and ``endpoint_gen`` records it.
    endpoint_generation: Optional[int] = None
    #: apply the ``trim_fc``/``trim_k`` truncation before reading the endpoint. Default False: a series that
    #: dies IS the depletion signal an endpoint measure is meant to report, and truncating would discard it.
    endpoint_trim: bool = False

    def __post_init__(self) -> None:
        if self.index not in ("slope", "endpoint"):
            raise ValueError(f"index must be 'slope' or 'endpoint', got {self.index!r}")
        if self.direction not in ("depleted", "both"):
            raise ValueError(f"direction must be 'depleted' or 'both', got {self.direction!r}")
        if self.endpoint_generation is not None and self.endpoint_generation == 0:
            raise ValueError("endpoint_generation 0 is generation 0 itself: its log2 fold change is 0 by "
                             "construction and carries no fitness information")

    def to_dict(self) -> dict:
        d = asdict(self)
        d["key_cols"] = list(self.key_cols)
        d["mutation_types"] = list(self.mutation_types)
        return d


# ===================================================================================== consolidation
def consolidate_normalized_counts(sets: Sequence[Tuple[str, str, pd.DataFrame]], reference: pd.DataFrame,
                                  annot: Mapping[str, str], pseudocount: float = 0.0, fill: str = "zero",
                                  mdf: str = "keep", join: str = "inner",
                                  mdf_col: str = "multiple_donor_bc0_fragments_passing_filter") -> pd.DataFrame:
    """Wide table of per-set normalized counts joined to the bc1 reference (``planning/289`` exactly).

    ``sets``: ``(treatment, replicate, table)`` in the order the columns should appear; each table has
    ``Barcode`` and ``Normalized_Count_<g>`` columns, renamed to ``<treatment>_set_<replicate>_Normalized_Count_<g>``.
    ``reference``: the bc1 reference table read with ``dtype=str``; ``annot`` maps its columns to the names the
    fitness step expects (gene_name, amino_acid_position, ref_amino_acid, mutant_amino_acid, the mdf flag).
    """
    if fill not in ("zero", "na") or mdf not in ("keep", "drop") or join not in ("inner", "left"):
        raise ValueError("fill must be zero|na, mdf keep|drop, join inner|left")
    wide = None
    for t, r, d in sets:
        norm = [c for c in d.columns if c.startswith("Normalized_Count_")]
        d = d[["Barcode"] + norm].rename(columns={c: f"{t}_set_{r}_{c}" for c in norm})
        wide = d if wide is None else wide.merge(d, on="Barcode", how="outer")
    if wide is None:
        raise ValueError("no count sets given")
    if fill == "zero":
        wide = wide.fillna(0)
    if pseudocount:
        cnt = [c for c in wide.columns if c != "Barcode"]
        wide[cnt] = wide[cnt] + pseudocount
    ref = reference.copy()
    ref["bc1"] = ref["bc1"].str.strip()
    cols = {k: v for k, v in annot.items() if k in ref.columns}
    if mdf == "drop":
        cols.pop(mdf_col, None)
    sub = ref[["bc1"] + list(cols)].rename(columns=cols)
    w = wide.rename(columns={"Barcode": "bc1"})
    out = sub.merge(w, on="bc1", how="inner") if join == "inner" else w.merge(sub, on="bc1", how="left")
    lead = ["bc1"] + list(cols.values())
    return out[lead + [c for c in out.columns if c not in lead]]


# ===================================================================================== long format
def _truthy(s: pd.Series) -> pd.Series:
    return s.isin([True, "TRUE", "True", "true"])


def mutation_type_of(ref: pd.Series, mut: pd.Series) -> pd.Series:
    """PTC when the mutant is ``*``, synonymous when ref == mut, else missense; missing when either is."""
    r = ref.where(ref.isna(), ref.astype(str).str.strip())
    m = mut.where(mut.isna(), mut.astype(str).str.strip())
    out = np.where(m == "*", "PTC", np.where(r == m, "synonymous", "missense")).astype(object)
    out[(r.isna() | m.isna()).to_numpy()] = None
    return pd.Series(out, index=ref.index, dtype=object)


def long_from_wide(wide: pd.DataFrame, mdf_col: str = "multiple_donor_bc0_fragments_passing_filter") -> pd.DataFrame:
    """The R script's long format from a consolidated wide table: one row per (barcode, sample), ordered
    row by row and, within a row, in column order (``pivot_longer``), which fixes "the first replicate"."""
    w = wide.copy()
    w["bc1"] = w["bc1"].astype(str)
    if mdf_col in w.columns:
        w = w[_truthy(w[mdf_col])]
    pos = pd.to_numeric(w["amino_acid_position"], errors="coerce")
    w["amino_acid_position"] = np.trunc(pos).astype("Int64")
    w["mutation_type"] = mutation_type_of(w.get("ref_amino_acid", pd.Series(np.nan, index=w.index)),
                                          w.get("mutant_amino_acid", pd.Series(np.nan, index=w.index)))
    norm = [c for c in w.columns if NORM_TAG in c]
    if not norm:
        raise ValueError(f"no '*{NORM_TAG}<generation>' columns")
    w = w.reset_index(drop=True)
    w["_row"] = np.arange(len(w))
    ids = ["_row", "bc1", "gene_name", "amino_acid_position", "mutation_type"]
    ids = [c for c in ids if c in w.columns]
    long = w[ids + norm].melt(id_vars=ids, value_vars=norm, var_name="sample", value_name="value")
    long["_col"] = long["sample"].map({c: i for i, c in enumerate(norm)})
    long = long.sort_values(["_row", "_col"], kind="stable").reset_index(drop=True)
    raw = long["sample"].str.replace(r"_Normalized_Count_.*", "", regex=True)
    long["replicate"] = raw.str.extract(r"_R(\d+)$", expand=False).astype("Int64")
    long["treatment"] = (raw.str.replace(r"_R\d+$", "", regex=True).str.replace(r"_set$", "", regex=True)
                         .str.replace(r"\s+", "_", regex=True))
    long["generation"] = long["sample"].str.extract(r"Normalized_Count_(\d+)", expand=False).astype(int)
    long["value"] = pd.to_numeric(long["value"], errors="coerce").astype(float)
    long = long[long["amino_acid_position"].notna()]
    return long.drop(columns=["_row", "_col"]).reset_index(drop=True)


def long_from_count_matrix(counts: pd.DataFrame, samples: pd.DataFrame, bc1_keys: pd.DataFrame,
                           pseudocount: float = 0.0, normalize: str = "cpm") -> pd.DataFrame:
    """Long table from a raw count matrix (bc1 x sample, e.g. step 08's), an EXPLICIT sample table and a bc1 ->
    edit-key table. Nothing is inferred from sample names.

    ``samples``: one row per count-matrix column to use, columns ``sample``, ``treatment``, ``replicate``,
    ``generation`` (integer); the row order sets "the first replicate" for per-barcode slopes.
    ``bc1_keys``: ``bc1`` plus the edit-key columns; only barcodes listed are kept (one row per bc1 required).
    ``normalize``: ``cpm`` (counts per million of the sample's total over ALL barcodes of the matrix, as the
    NRD1 per-set counts were), ``median_depth`` (rescaled to the median sample depth, so ``min_g0`` stays a read
    floor in libraries of millions of barcodes, where a CPM floor of 10 would drop nearly all of them) or ``none``. ``pseudocount`` is added after normalizing, as in the NRD1 run.
    """
    need = {"sample", "treatment", "replicate", "generation"}
    if not need <= set(samples.columns):
        raise ValueError(f"samples needs columns {sorted(need)}")
    missing = [s for s in samples["sample"] if s not in counts.columns]
    if missing:
        raise ValueError(f"samples not in the count matrix: {missing[:10]}")
    if bc1_keys["bc1"].duplicated().any():
        raise ValueError("bc1_keys has duplicate bc1: resolve the edit of each barcode before fitting")
    gens = pd.to_numeric(samples["generation"], errors="coerce")
    if gens.isna().any() or (gens != np.round(gens)).any():
        raise ValueError("every sample needs an integer generation")
    m = counts[list(samples["sample"])].astype(float)
    if normalize == "cpm":
        m = m / m.sum(axis=0) * 1e6
    elif normalize == "median_depth":  # counts rescaled to the median sample depth: a floor in reads, not per million
        depth = m.sum(axis=0)
        m = m / depth * float(np.median(depth))
    elif normalize != "none":
        raise ValueError("normalize must be cpm, median_depth or none")
    m = m + pseudocount
    m = m.loc[m.index.intersection(pd.Index(bc1_keys["bc1"]))]
    meta = samples.assign(generation=gens.astype(int)).set_index("sample")
    long = m.rename_axis("bc1").reset_index().melt(id_vars="bc1", var_name="sample", value_name="value")
    order = {s: i for i, s in enumerate(samples["sample"])}
    long["_c"] = long["sample"].map(order)
    long = long.sort_values(["bc1", "_c"], kind="stable").drop(columns="_c")
    long = long.join(meta[["treatment", "replicate", "generation"]], on="sample")
    return long.merge(bc1_keys, on="bc1", how="inner").reset_index(drop=True)


# ===================================================================================== the method
def _ols_slope(x: np.ndarray, y: np.ndarray, min_points: int) -> float:
    if len(x) < min_points:
        return np.nan
    xm, ym = x.mean(), y.mean()
    sxx = ((x - xm) ** 2).sum()
    return float(((x - xm) * (y - ym)).sum() / sxx) if sxx > 0 else np.nan


def _first_per_x_slope(x: np.ndarray, y: np.ndarray, min_points: int) -> float:
    """R ``safe_lm_slope``: finite pairs, the FIRST row at each x, then OLS."""
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) == 0:
        return np.nan
    _, first = np.unique(x, return_index=True)
    first.sort()
    x, y = x[first], y[first]
    o = np.argsort(x, kind="stable")
    return _ols_slope(x[o], y[o], min_points)


def _cut_generation(g: np.ndarray, y: np.ndarray, trim_fc: float, k: int) -> float:
    """R ``trim_series`` reduced to the generation it keeps up to (inclusive); -inf when nothing is finite."""
    o = np.argsort(g, kind="stable")
    g, y = g[o], y[o]
    ok = np.isfinite(y) & np.isfinite(g)
    g, y = g[ok], y[ok]
    if len(g) == 0:
        return -np.inf
    if not np.any(g == 0):
        return float(g.max())
    post = g > 0
    gp, yp = g[post], y[post]
    run = 0
    for gi, below in zip(gp, yp < np.log2(trim_fc)):
        run = run + 1 if below else 0
        if run >= k:
            return float(gi)
    return float(g.max())


def _wide_by_treatment(d: pd.DataFrame, index: List[str], values: List[str], treatments: List[str],
                       name=lambda v, t: f"{v}_{t}") -> pd.DataFrame:
    """``pivot_wider`` without pandas' pivot_table semantics: one row per index value present in ``d``
    (missing keys kept), one column per (value, treatment), in ``values`` x ``treatments`` order."""
    base = d[index].drop_duplicates().reset_index(drop=True)
    for v in values:
        for t in treatments:
            part = d.loc[d["treatment"] == t, index + [v]].drop_duplicates(index).rename(columns={v: name(v, t)})
            base = base.merge(part, on=index, how="left")
    return base


def _quantile7(v: np.ndarray, p: float) -> float:
    v = v[np.isfinite(v)]
    return float(np.quantile(v, p)) if len(v) else np.nan


def _endpoint(x: np.ndarray, y: np.ndarray, at: Optional[int]) -> Tuple[float, float]:
    """The log2 fold change at one generation, and the generation it was read at.

    ``at`` None takes the last generation carrying a finite value, excluding generation 0, whose fold change
    is 0 by construction. A requested generation that this series does not carry gives (nan, nan) rather
    than the nearest one, so a missing timepoint is never silently substituted."""
    m = np.isfinite(y) & (x > 0)
    if not m.any():
        return np.nan, np.nan
    xs, ys = x[m], y[m]
    if at is None:
        i = int(np.argmax(xs))
        return float(ys[i]), float(xs[i])
    hit = np.flatnonzero(xs == float(at))
    if not len(hit):
        return np.nan, np.nan
    return float(ys[hit[0]]), float(at)


def pooled_slope_fitness(long: pd.DataFrame, cfg: PooledSlopeConfig) -> Dict[str, pd.DataFrame]:
    """Run the method on a long table (``bc1``, ``key_cols``, ``treatment``, ``generation``, ``value``; row
    order = input column order within each barcode). Returns the R script's six tables under its file
    stems plus ``hits_with_support`` (criteria 1-4)."""
    K = list(cfg.key_cols)
    KT = K + ["treatment"]
    if not (long["generation"] == 0).any():
        raise ValueError("no generation-0 samples: every fold change is relative to generation 0, refusing")

    g0 = long[long["generation"] == 0].groupby("bc1")["value"].mean().rename("g0")
    keep = g0[np.isfinite(g0) & (g0 >= cfg.min_g0)].index
    df = long[long["bc1"].isin(keep)].merge(g0, left_on="bc1", right_index=True, how="left")
    df["value"] = np.where(df["generation"] == 0, df["g0"], df["value"])
    df["log2fc"] = np.log2(np.maximum(df["value"], cfg.eps) / np.maximum(df["g0"], cfg.eps))
    df.loc[df["value"].isna(), "log2fc"] = np.nan
    df["w"] = np.maximum(df["g0"], 1.0) if cfg.weighted else 1.0

    # 2. the pooled series
    fin = df["log2fc"].notna()
    tmp = df.assign(wx=np.where(fin, df["w"] * df["log2fc"], 0.0), ww=np.where(fin, df["w"], 0.0))
    agg = tmp.groupby(KT + ["generation"], dropna=False)[["wx", "ww"]].sum().reset_index()
    with np.errstate(invalid="ignore", divide="ignore"):
        agg["mean_log2fc"] = agg["wx"] / agg["ww"]
    mean_tab = agg[KT + ["generation", "mean_log2fc"]].copy()

    # 3. the truncation point of each series, then its raw index (OLS slope, or the endpoint fold change)
    endpoint = cfg.index == "endpoint"
    if endpoint and not (long["generation"] > 0).any():
        raise ValueError("no generation above 0: the endpoint index is a fold change against generation 0, refusing")
    cuts, slopes = [], []
    for key, g in mean_tab.groupby(KT, dropna=False, sort=True):
        x, y = g["generation"].to_numpy(float), g["mean_log2fc"].to_numpy(float)
        cut = _cut_generation(x, y, cfg.trim_fc, cfg.trim_k)
        cuts.append((*key, cut))
        if endpoint:
            xw, yw = (x, y) if not cfg.endpoint_trim else (x[x <= cut], y[x <= cut])
            v, at = _endpoint(xw, yw, cfg.endpoint_generation)
            if np.isfinite(v):   # a series with no usable post-zero point gets no index row
                slopes.append((*key, v, at))
        else:
            m = np.isfinite(y) & (x <= cut)
            if m.any():  # R summarises only the rows that survive the filter: an empty series has no slope row
                slopes.append((*key, _first_per_x_slope(x[m], y[m], cfg.min_points), np.nan))
    trim = pd.DataFrame(cuts, columns=KT + ["cut_gen"])
    st = pd.DataFrame(slopes, columns=KT + ["slope", "endpoint_gen"])
    if not endpoint:
        st = st.drop(columns=["endpoint_gen"])
    st["median_slope"] = st.groupby("treatment")["slope"].transform(lambda s: np.nanmedian(s) if s.notna().any() else np.nan)
    st["slope_adj"] = st["slope"] - st["median_slope"]

    # 5. depletion
    dep_rows = []
    for key, g in mean_tab.groupby(KT, dropna=False, sort=True):
        f = g[np.isfinite(g["mean_log2fc"])]
        if f.empty:
            continue
        last = f.loc[f["generation"] == f["generation"].max(), "mean_log2fc"].iloc[0]
        z = g.loc[g["generation"] == 0, "mean_log2fc"]
        dep_rows.append((*key, z.iloc[0] if len(z) else np.nan, last, bool(np.isfinite(last) and last < np.log2(cfg.depletion_fc))))
    dep = pd.DataFrame(dep_rows, columns=KT + ["log2fc_gen0", "log2fc_last", "depletion"])

    out: Dict[str, pd.DataFrame] = {"mean_log2FC_by_group": mean_tab, "depletion_log2FC": dep}
    T, C = cfg.treatment, cfg.control
    thr = _quantile7(st.loc[st["treatment"] == T, "slope_adj"].to_numpy(float), cfg.tail_prob)

    # 6. criteria 1-3
    hits = None
    if C is not None:
        sw = st[st["treatment"].isin([T, C]) & st["mutation_type"].isin(cfg.mutation_types)] \
            if "mutation_type" in K else st[st["treatment"].isin([T, C])]
        sw = _wide_by_treatment(sw, K, ["slope", "median_slope", "slope_adj"], [T, C])
        sw = sw.merge(dep.loc[dep["treatment"] == T, K + ["depletion"]], on=K, how="left")
        a, b = sw[f"slope_adj_{T}"], sw[f"slope_adj_{C}"]
        sel = a.notna() & b.notna() & sw["depletion"].eq(True) & (a < thr) & (a / cfg.ratio_factor < b)
        hits = sw[sel].copy()
        hits["hit"] = True

    # 7. barcode support
    d2 = df.merge(trim, on=KT, how="inner")
    keep_window = np.isfinite(d2["log2fc"])
    if not endpoint or cfg.endpoint_trim:
        keep_window &= d2["generation"] <= d2["cut_gen"]
    d2 = d2[keep_window]
    rows = []
    for key, g in d2.groupby(KT + ["bc1"], dropna=False, sort=False):
        gx, gy = g["generation"].to_numpy(float), g["log2fc"].to_numpy(float)
        if endpoint:
            # the same first-replicate-per-generation rule the slope route inherits from R's distinct()
            _, first = np.unique(gx, return_index=True)
            v = _endpoint(gx[first], gy[first], cfg.endpoint_generation)[0]
        else:
            v = _first_per_x_slope(gx, gy, cfg.min_points)
        rows.append((*key, v, g["w"].iloc[0]))
    pb = pd.DataFrame(rows, columns=KT + ["bc1", "slope", "w"])
    pb["slope_adj"] = pb["slope"] - pb.groupby("treatment")["slope"].transform(
        lambda s: np.nanmedian(s) if s.notna().any() else np.nan)

    dt = df[df["treatment"] == T].copy()
    dt["_o"] = np.arange(len(dt))
    dt = dt.sort_values(K + ["bc1", "generation", "_o"], kind="stable", na_position="last")
    brow = []
    for key, g in dt.groupby(K + ["bc1"], dropna=False, sort=False):
        f = g[np.isfinite(g["log2fc"])]
        last = f["log2fc"].iloc[int(np.argmax(f["generation"].to_numpy()))] if len(f) else np.nan
        brow.append((*key, last, bool(np.isfinite(last) and last < np.log2(cfg.depletion_fc))))
    bdep = pd.DataFrame(brow, columns=K + ["bc1", "log2fc_last", "depletion_bc"])

    keep_t = [T] + ([C] if C is not None else [])
    bw = _wide_by_treatment(pb[pb["treatment"].isin(keep_t)], K + ["bc1"], ["slope_adj"], keep_t,
                            name=lambda v, t: t)
    bw = bw.merge(bdep, on=K + ["bc1"], how="left")
    ta = bw[T]
    dep_bc = bw["depletion_bc"].astype("boolean")
    if C is not None:
        cb = bw[C]
        bw["supports_hit"] = (ta.notna() & cb.notna()).astype("boolean") & dep_bc & (ta < thr).astype("boolean") \
            & (ta / cfg.ratio_factor < cb).astype("boolean")
        bw["agrees_direction"] = ta.notna() & cb.notna() & (ta < cb) & (ta < 0)
    else:
        bw["supports_hit"] = ta.notna().astype("boolean") & dep_bc & (ta < thr).astype("boolean")
        bw["agrees_direction"] = ta.notna() & (ta < 0)
    sup = bw.groupby(K, dropna=False).agg(n_barcodes=("bc1", "size"),
                                         n_support=("supports_hit", lambda s: int(s.fillna(False).sum())),
                                         n_agree=("agrees_direction", "sum")).reset_index()
    sup["pct_support"] = sup["n_support"] / sup["n_barcodes"]
    sup["pct_agree"] = sup["n_agree"] / sup["n_barcodes"]

    st = st.merge(sup, on=K, how="left")
    if hits is not None:
        hits = hits.merge(sup, on=K, how="left")
        hk = hits[K].assign(hit=True)
        st = st.merge(hk, on=K, how="left")
        st["hit"] = st["hit"].eq(True)
        hrows = st[st["hit"]].drop_duplicates(K)[K].assign(treatment="HIT", slope=np.nan, median_slope=np.nan,
                                                            slope_adj=np.nan, hit=True)
        st_out = pd.concat([st, hrows], ignore_index=True)
        final = hits[hits["n_support"].fillna(0) >= cfg.min_support].copy()
    else:
        st_out, final = st, None
    out.update({"slope_table_log2FC": st_out, "barcode_support_flags": bw, "hit_barcode_support_summary": sup})
    if hits is not None:
        out["hits_summary"] = hits
        out["hits_with_support"] = final
    thr_e = np.nan
    if cfg.direction == "both":
        thr_e = _quantile7(st.loc[st["treatment"] == T, "slope_adj"].to_numpy(float), 1 - cfg.tail_prob)
        up = -np.log2(cfg.depletion_fc)
        bw_e = bw[K + ["bc1", T] + ([C] if C is not None else []) + ["log2fc_last"]].copy()
        ta = bw_e[T]
        enr_bc = (bw_e["log2fc_last"] > up).astype("boolean").mask(bw_e["log2fc_last"].isna())
        if C is not None:
            cb = bw_e[C]
            bw_e["supports_hit"] = (ta.notna() & cb.notna()).astype("boolean") & enr_bc & (ta > thr_e).astype("boolean") \
                & (ta / cfg.ratio_factor > cb).astype("boolean")
            bw_e["agrees_direction"] = ta.notna() & cb.notna() & (ta > cb) & (ta > 0)
        else:
            bw_e["supports_hit"] = ta.notna().astype("boolean") & enr_bc & (ta > thr_e).astype("boolean")
            bw_e["agrees_direction"] = ta.notna() & (ta > 0)
        sup_e = bw_e.groupby(K, dropna=False).agg(n_barcodes=("bc1", "size"),
                                                 n_support=("supports_hit", lambda s: int(s.fillna(False).sum())),
                                                 n_agree=("agrees_direction", "sum")).reset_index()
        sup_e["pct_support"] = sup_e["n_support"] / sup_e["n_barcodes"]
        sup_e["pct_agree"] = sup_e["n_agree"] / sup_e["n_barcodes"]
        out.update({"barcode_support_flags_enriched": bw_e, "hit_barcode_support_summary_enriched": sup_e})
        if C is not None:
            sw = st[st["treatment"].isin([T, C]) & st["mutation_type"].isin(cfg.mutation_types)] \
                if "mutation_type" in K else st[st["treatment"].isin([T, C])]
            sw = _wide_by_treatment(sw, K, ["slope", "median_slope", "slope_adj"], [T, C])
            last = dep.loc[dep["treatment"] == T, K + ["log2fc_last"]]
            sw = sw.merge(last, on=K, how="left")
            sw["enrichment"] = sw["log2fc_last"] > up
            a, b = sw[f"slope_adj_{T}"], sw[f"slope_adj_{C}"]
            sel = a.notna() & b.notna() & sw["enrichment"] & (a > thr_e) & (a / cfg.ratio_factor > b)
            he = sw[sel].copy()
            he["hit"] = True
            he = he.merge(sup_e, on=K, how="left")
            out["hits_summary_enriched"] = he
            out["hits_with_support_enriched"] = he[he["n_support"].fillna(0) >= cfg.min_support].copy()
    out["_parameters"] = pd.DataFrame([{**cfg.to_dict(), "hit_threshold": thr, "hit_threshold_enriched": thr_e}])
    return out


def write_outputs(res: Dict[str, pd.DataFrame], out_dir: str) -> List[str]:
    import os
    os.makedirs(out_dir, exist_ok=True)
    paths = []
    for name, d in res.items():
        p = os.path.join(out_dir, f"{name.lstrip('_')}.csv")
        d.to_csv(p, index=False)
        paths.append(p)
    return paths
