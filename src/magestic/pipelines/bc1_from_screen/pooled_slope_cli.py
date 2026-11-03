"""``magestic-pooled-slope``: the NRD1 pooled-slope fitness method (core.pooled_slope_fitness) from the shell.

Three ways in, all explicit (nothing is inferred from file or sample names beyond the documented column grammar):

  --consolidated CSV                          a wide table of ``<treatment>_set_<R>_Normalized_Count_<g>``
                                              columns, as the NRD1 fitness screen used
  --sets JSON --reference TSV                 build that table first (planning/289's consolidation): JSON =
                                              {"sets": [[treatment, replicate, csv], ...], "annot": {...},
                                               "pseudocount": 1, "fill": "zero", "mdf": "keep", "join": "inner"}
  --count-matrix TSV --samples TSV --bc1-keys TSV
                                              a raw bc1 x sample matrix (step 08), an explicit sample table
                                              (sample, treatment, replicate, generation) and a bc1 -> edit-key
                                              table; CPM-normalized, pseudocount from --config

``--config`` (JSON) sets PooledSlopeConfig fields; anything not named keeps the R script's default, and the
resolved parameters are written to ``parameters.csv`` and the provenance sidecar. The NRD1 Supplementary Table 1
settings are in ``reproduce/20250813_NRD1_fitness/config.json``.

The ENDPOINT route needs no new flag: ``{"index": "endpoint"}`` in that JSON scores each edit by its pooled log2
fold change against generation 0 instead of a slope, so a screen measured at generation 0 and one later point is
callable (``{"index": "endpoint", "endpoint_generation": 12}`` pins the generation). Everything after the index is
shared with the slope route, and ``parameters.csv`` records which route ran.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import List, Optional

import pandas as pd

from .core import pooled_slope_fitness as psf


def _config(path: Optional[str]) -> dict:
    return json.load(open(path)) if path else {}


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(prog="magestic-pooled-slope", description=__doc__.split("\n\n")[0])
    ap.add_argument("--consolidated")
    ap.add_argument("--sets")
    ap.add_argument("--reference")
    ap.add_argument("--count-matrix")
    ap.add_argument("--samples")
    ap.add_argument("--bc1-keys")
    ap.add_argument("--config")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args(argv)

    from magestic.utils.provenance import StepProvenance

    raw = _config(a.config)
    fields = {k: v for k, v in raw.items() if k in psf.PooledSlopeConfig.__dataclass_fields__}
    for k in ("key_cols", "mutation_types"):
        if k in fields:
            fields[k] = tuple(fields[k])
    cfg = psf.PooledSlopeConfig(**fields)
    os.makedirs(a.out_dir, exist_ok=True)
    prov = StepProvenance("bc1_from_screen", "pooled_slope_fitness", os.path.join(a.out_dir, "hits_with_support.csv"),
                          description="NRD1 pooled-slope fitness: pooled weighted series, truncation, median-centered "
                                      "slope, criteria 1-3 and barcode support (criterion 4)")
    if a.config:
        prov.add_input(a.config, role="config")

    if a.consolidated:
        wide = pd.read_csv(a.consolidated)
        prov.add_input(a.consolidated, role="consolidated_counts")
        long = psf.long_from_wide(wide)
    elif a.sets:
        if not a.reference:
            ap.error("--sets needs --reference")
        spec = json.load(open(a.sets))
        sets = [(t, r, pd.read_csv(p)) for t, r, p in spec["sets"]]
        for _, _, p in spec["sets"]:
            prov.add_input(p, role="per_set_normalized_counts")
        ref = pd.read_csv(a.reference, sep="\t", low_memory=False, dtype=str)
        prov.add_input(a.reference, role="bc1_reference")
        wide = psf.consolidate_normalized_counts(sets, ref, spec["annot"], spec.get("pseudocount", 0.0),
                                                 spec.get("fill", "zero"), spec.get("mdf", "keep"),
                                                 spec.get("join", "inner"))
        wide.to_csv(os.path.join(a.out_dir, "consolidated_normalized_counts.csv"), index=False)
        long = psf.long_from_wide(wide)
    elif a.count_matrix:
        if not (a.samples and a.bc1_keys):
            ap.error("--count-matrix needs --samples and --bc1-keys")
        counts = pd.read_csv(a.count_matrix, sep="\t", index_col=0)
        samples = pd.read_csv(a.samples, sep="\t")
        keys = pd.read_csv(a.bc1_keys, sep="\t", dtype={"bc1": str})
        for p, r in ((a.count_matrix, "count_matrix"), (a.samples, "samples"), (a.bc1_keys, "bc1_keys")):
            prov.add_input(p, role=r)
        long = psf.long_from_count_matrix(counts, samples, keys, raw.get("pseudocount", 0.0),
                                          raw.get("normalize", "cpm"))
    else:
        ap.error("give --consolidated, --sets with --reference, or --count-matrix with --samples and --bc1-keys")

    res = psf.pooled_slope_fitness(long, cfg)
    paths = psf.write_outputs(res, a.out_dir)
    prov.set_parameters(**cfg.to_dict())
    counts = {k: int(len(v)) for k, v in res.items() if not k.startswith("_")}
    prov.set_counts(**counts)
    for p in paths:
        prov.add_output(p)
    prov.write()
    print(json.dumps({"out_dir": a.out_dir, "rows": counts,
                      "hit_threshold": float(res["_parameters"]["hit_threshold"].iloc[0])}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
