"""Steps 10-18 of the screen route, run the way the shipped rules run them, with every output checked.

This is what ``magestic-screen-analysis --route-config <yaml>`` executes (planning/272). It is a thin
orchestrator: each stage is the numbered script under ``snakemake/``, called with the argv its rule in
``Snakefile.rules`` renders, so a stage behaves the same from Snakemake and from here. It adds two things
the rules cannot give a caller:

* **No silent success** (docs/ARCHITECTURE_POLICY.md section 7). A stage passes only when its process exits 0
  AND every output it declares exists with at least one data row. A script that exits 0 without writing is
  a failure here, not a result.
* **A run manifest** (``<output_dir>/run_manifest.json``): per stage the argv, exit code, wall time and each
  output's bytes and row count, plus the code identity, written after every stage so a failed run still
  says how far it got.

The route config is explicit. Nothing is inferred from a project layout, and there are no screen presets:
``--project-dir`` silently loading the SpG preset is exactly the false-default this module avoids.

Route config (YAML)::

    screen: 20240806_yL406_SpG_QTL_screen
    library: yL406
    output_dir: /oak/.../tmp/272_armD/runs/yL406
    inputs:
      count_matrix: .../yL406_bc1_counts_wide_combined_tiers.tsv     # step 08/09 output
      sample_key: .../sample_key.tsv
      bc1_reference_table: .../bc1_reference_table.tsv
      tier_mapping: .../bc1_tier_mapping_yL406.tsv
      bc1_counts_long: .../bc1_counts_long_format.tsv                 # step 18 only
      common_dir: .../common            # MAGESTIC_COMMON_DIR for step 13's registry (optional)
    split: {transformation_scope: pooled}       # or: split + transformation_column; pooled_dates: true
    comparisons: [yL406_Azauracil_12, ...]      # optional: run 11/13/14/17 on these only (default: all)
    deseq2: {min_basemean: 0, size_factors_fit_type: ratio}   # ratio | poscounts | iterative
    slope: {conditions: [Azauracil]}
    hit_calling: {tier_1_z: 3.0, tier_1_concordance: 0.75, tier_2_z: 3.75, tier_2_concordance: 0.85,
                  use_synonymous_null: true}
    rarefaction: {samples: {<sample_name>: <parsed _bc1_counts.tsv>}, data_type: bc1_timepoints}
    deseq_rarefaction: {comparisons: [...], fractions: [0.5, 1.0], replicates: 2}
    composition: {consolidated_library: yL406, yeast_sublibraries: [yL386, ...],
                  expected_proportions: [3, 3, 3, 1, 1], filter_condition: ..., min_counts: 5}

Author: Kevin R. Roy
"""

from __future__ import annotations

import json
import logging
import os
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import pandas as pd

logger = logging.getLogger(__name__)

SCRIPTS = Path(__file__).resolve().parent.parent / "snakemake"

#: Stage name -> numbered script. Order is execution order.
STAGES = {
    "split": "10_split_comparisons.py",
    "deseq2": "11_deseq2_single.py",
    "slope": "12_slope_fitness.py",
    "annotate": "13_annotate_results.py",
    "hitcall": "14_hit_calling.py",
    "aggregate": "15_aggregate_hits.py",
    "rarefaction": "16_rarefaction.py",
    "deseq_rarefaction": "17_deseq_rarefaction.py",
    "composition": "18_validate_sublibrary_composition.py",
}


class StageFailed(RuntimeError):
    """A stage exited non-zero, or exited 0 without writing every output it declares."""


@dataclass
class RouteConfig:
    screen: str
    library: str
    output_dir: Path
    inputs: Dict[str, Path]
    split: Dict[str, Any] = field(default_factory=dict)
    comparisons: Optional[List[str]] = None
    deseq2: Dict[str, Any] = field(default_factory=dict)
    slope: Dict[str, Any] = field(default_factory=dict)
    hit_calling: Dict[str, Any] = field(default_factory=dict)
    rarefaction: Dict[str, Any] = field(default_factory=dict)
    deseq_rarefaction: Dict[str, Any] = field(default_factory=dict)
    composition: Dict[str, Any] = field(default_factory=dict)
    python: str = sys.executable

    @classmethod
    def from_yaml(cls, path: Path) -> "RouteConfig":
        import yaml
        with open(path) as fh:
            raw = yaml.safe_load(fh)
        missing = [k for k in ("screen", "library", "output_dir", "inputs") if k not in (raw or {})]
        if missing:
            raise ValueError(f"route config {path} lacks {missing}")
        raw["output_dir"] = Path(raw["output_dir"])
        raw["inputs"] = {k: Path(v) for k, v in raw["inputs"].items()}
        known = set(cls.__dataclass_fields__)
        unknown = sorted(set(raw) - known)
        if unknown:
            raise ValueError(f"route config {path} has unknown keys {unknown}; known: {sorted(known)}")
        return cls(**raw)

    # Output layout, mirroring Snakefile.rules (BC1_COUNTS_DIR / DESEQ2_DIR / HIT_CALLING_DIR).
    @property
    def split_dir(self) -> Path:
        return self.output_dir / "counts_for_deseq" / self.library

    @property
    def deseq2_dir(self) -> Path:
        return self.output_dir / "DESeq2" / self.library

    @property
    def annotated_dir(self) -> Path:
        return self.output_dir / "DESeq2" / "per_comparison" / self.library

    @property
    def hits_dir(self) -> Path:
        return self.output_dir / "hit_calling" / self.library

    def need(self, *names: str) -> List[Path]:
        absent = [n for n in names if n not in self.inputs]
        if absent:
            raise StageFailed(f"route config inputs lack {absent}")
        paths = [self.inputs[n] for n in names]
        missing = [str(p) for p in paths if not p.exists()]
        if missing:
            raise StageFailed(f"input file(s) do not exist: {missing}")
        return paths


def _data_rows(path: Path) -> int:
    """Data rows of a TSV (header excluded); -1 for a non-table output."""
    if path.suffix not in (".tsv", ".txt", ".csv"):
        return -1
    with open(path, "rb") as fh:
        return max(sum(1 for _ in fh) - 1, 0)


class Route:
    """Runs stages in order and keeps the manifest."""

    def __init__(self, config: RouteConfig, *, git_state: Optional[dict] = None):
        self.config = config
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        self.logs = self.config.output_dir / "logs"
        self.logs.mkdir(exist_ok=True)
        self.manifest: Dict[str, Any] = {
            "screen": config.screen, "library": config.library, "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "code": git_state, "python": config.python,
            "inputs": {k: str(v) for k, v in config.inputs.items()}, "stages": [],
        }
        self._write_manifest()

    def _write_manifest(self) -> None:
        target = self.config.output_dir / "run_manifest.json"
        tmp = target.with_name(target.name + ".tmp")
        tmp.write_text(json.dumps(self.manifest, indent=1, default=str) + "\n")
        os.replace(tmp, target)

    def _env(self) -> Dict[str, str]:
        env = dict(os.environ)
        if "common_dir" in self.config.inputs:
            env["MAGESTIC_COMMON_DIR"] = str(self.config.inputs["common_dir"])
        env.setdefault("MPLBACKEND", "Agg")
        return env

    def call(self, stage: str, label: str, args: Sequence[Any], outputs: Sequence[Path],
             *, min_rows: int = 1) -> List[Dict[str, Any]]:
        """Run one invocation of a stage's script; refuse unless every declared output has data."""
        argv = [self.config.python, str(SCRIPTS / STAGES[stage]), *map(str, args)]
        log = self.logs / f"{stage}__{label}.log"
        record: Dict[str, Any] = {"stage": stage, "label": label, "argv": argv, "log": str(log)}
        self.manifest["stages"].append(record)
        start = time.time()
        with open(log, "w") as fh:
            rc = subprocess.run(argv, stdout=fh, stderr=subprocess.STDOUT, env=self._env()).returncode
        record.update(rc=rc, seconds=round(time.time() - start, 1))
        outs = []
        for path in outputs:
            path = Path(path)
            exists = path.exists()
            outs.append({"path": str(path), "bytes": path.stat().st_size if exists else None,
                         "rows": _data_rows(path) if exists else None})
        record["outputs"] = outs
        self._write_manifest()
        if rc != 0:
            tail = "".join(open(log).readlines()[-15:])
            raise StageFailed(f"{stage} [{label}] exited {rc}; log {log}:\n{tail}")
        empty = [o["path"] for o in outs if o["bytes"] in (None, 0) or (o["rows"] is not None and 0 <= o["rows"] < min_rows)]
        if empty:
            raise StageFailed(f"{stage} [{label}] exited 0 but wrote no data to {empty}; log {log}")
        logger.info(f"  {stage} [{label}]: ok in {record['seconds']} s; "
                    + ", ".join(f"{Path(o['path']).name}={o['rows']} rows" for o in outs))
        return outs

    # -------------------------------------------------------------------------------------------
    def comparisons(self) -> pd.DataFrame:
        manifest = self.config.split_dir / "comparison_manifest.tsv"
        if not manifest.exists():
            raise StageFailed(f"no step-10 manifest at {manifest}; run the 'split' stage first")
        table = pd.read_csv(manifest, sep="\t")
        wanted = self.config.comparisons
        if wanted:
            unknown = sorted(set(wanted) - set(table["comparison_id"]))
            if unknown:
                raise StageFailed(f"comparisons {unknown} are not in the step-10 manifest {manifest}")
            table = table[table["comparison_id"].isin(wanted)]
        if table.empty:
            raise StageFailed(f"no comparisons to run from {manifest}")
        return table

    def run_split(self) -> None:
        c = self.config
        count_matrix, sample_key = c.need("count_matrix", "sample_key")
        scope = c.split.get("transformation_scope")
        if scope == "pooled":
            scope_args = ["--pool-transformations"]
        elif scope == "split":
            (reference,) = c.need("bc1_reference_table")
            scope_args = ["--split-by-yeast-sublibrary", "--bc1-reference-table", reference,
                          "--yeast-sublibrary-column", c.split["transformation_column"]]
        else:
            raise StageFailed("split.transformation_scope must be 'pooled' or 'split' (step 10 has no default)")
        args = ["--count-matrix", count_matrix, "--sample-key", sample_key, "--library", c.library,
                "--output-dir", c.split_dir, *scope_args]
        if c.split.get("pooled_dates"):
            args.append("--no-consolidate")
        if c.split.get("exclude_samples_file"):
            args += ["--exclude-samples-file", c.split["exclude_samples_file"]]
        self.call("split", c.library, args, [c.split_dir / "comparison_manifest.tsv"])

    def run_deseq2(self) -> None:
        c = self.config
        for row in self.comparisons().itertuples():
            out = c.deseq2_dir / f"{row.comparison_id}_DESeq2_results.tsv"
            self.call("deseq2", row.comparison_id, [
                "--reference", c.split_dir / "reference" / row.reference_file,
                "--treatment", c.split_dir / "treatment" / row.treatment_file,
                "--output", out, "--comparison", row.comparison_id, "--library", c.library,
                "--min-basemean", c.deseq2.get("min_basemean", 0),
                "--size-factors-fit-type", c.deseq2.get("size_factors_fit_type", "ratio")], [out])

    def run_slope(self) -> None:
        c = self.config
        count_matrix, sample_key, reference = c.need("count_matrix", "sample_key", "bc1_reference_table")
        conditions = c.slope.get("conditions")
        if not conditions:
            raise StageFailed("slope.conditions is empty: name the condition(s) to fit (step 12 fits one each)")
        out_dir = c.deseq2_dir / "slope_fitness"
        for condition in conditions:
            stem = f"{c.screen}_{condition}".replace(" ", "_")
            outs = [out_dir / f"{stem}_bc1_slope_fitness.tsv", out_dir / f"{stem}_variant_slope_fitness.tsv",
                    out_dir / f"{stem}_slope_fitness_summary.tsv"]
            self.call("slope", stem, [
                "--count-matrix", count_matrix, "--sample-key", sample_key, "--bc1-reference", reference,
                "--output-bc1", outs[0], "--output-variant", outs[1], "--output-summary", outs[2],
                "--library", c.library, "--condition", condition], outs)

    def run_annotate(self) -> None:
        c = self.config
        reference, tiers = c.need("bc1_reference_table", "tier_mapping")
        registry = c.inputs.get("annotation_registry", reference)
        for row in self.comparisons().itertuples():
            out = c.annotated_dir / f"{row.comparison_id}_annotated.tsv"
            self.call("annotate", row.comparison_id, [
                "--deseq2-results", c.deseq2_dir / f"{row.comparison_id}_DESeq2_results.tsv",
                "--ref-table", reference, "--tier-mapping", tiers, "--oligo-design", registry,
                "--output", out, "--library", c.library], [out])

    def run_hitcall(self) -> None:
        c, h = self.config, self.config.hit_calling
        for row in self.comparisons().itertuples():
            hits, stats = c.hits_dir / f"hits_{row.comparison_id}.tsv", c.hits_dir / f"{row.comparison_id}_stats.tsv"
            self.call("hitcall", row.comparison_id, [
                "--input", c.annotated_dir / f"{row.comparison_id}_annotated.tsv",
                "--output-hits", hits, "--output-stats", stats, "--library", c.library,
                "--comparison", row.comparison_id,
                "--tier-1-z", h.get("tier_1_z", 3.0), "--tier-1-concordance", h.get("tier_1_concordance", 0.75),
                "--tier-2-z", h.get("tier_2_z", 3.75), "--tier-2-concordance", h.get("tier_2_concordance", 0.85),
                "--use-synonymous-null", h.get("use_synonymous_null", True)], [hits, stats])

    def run_aggregate(self) -> None:
        c = self.config
        out_dir = c.hits_dir / "aggregated"
        self.call("aggregate", c.library, ["--input-dir", c.hits_dir, "--output-dir", out_dir,
                                           "--screen-name", c.library],
                  [out_dir / "hits_summary_by_tier.tsv"])

    def run_rarefaction(self) -> None:
        c = self.config
        samples = c.rarefaction.get("samples") or {}
        if not samples:
            raise StageFailed("rarefaction.samples is empty: map sample names to their parsed _bc1_counts.tsv")
        out_dir = c.output_dir / "qc" / "rarefaction"
        for sample, table in samples.items():
            if not Path(table).exists():
                raise StageFailed(f"rarefaction input for {sample} does not exist: {table}")
            out = out_dir / f"{sample}_rarefaction.tsv"
            self.call("rarefaction", sample, [
                "--input", table, "--output", out, "--data-type", c.rarefaction.get("data_type", "bc1_timepoints"),
                "--sample-name", sample, "--barcode-column", "bc1", "--count-column", "counts"],
                [out, out.with_suffix(".summary.tsv")])

    def run_deseq_rarefaction(self) -> None:
        c, r = self.config, self.config.deseq_rarefaction
        table = self.comparisons()
        wanted = r.get("comparisons")
        if wanted:
            table = table[table["comparison_id"].isin(wanted)]
            if table.empty:
                raise StageFailed(f"deseq_rarefaction.comparisons {wanted} match no comparison being run")
        for row in table.itertuples():
            out_dir = c.output_dir / "DESeq2_rarefaction" / c.library / row.comparison_id
            self.call("deseq_rarefaction", row.comparison_id, [
                "--reference", c.split_dir / "reference" / row.reference_file,
                "--treatment", c.split_dir / "treatment" / row.treatment_file,
                "--output-dir", out_dir, "--comparison", row.comparison_id, "--library", c.library,
                "--fractions", *r.get("fractions", [0.1, 0.25, 0.5, 0.75, 1.0]),
                "--replicates", r.get("replicates", 5), "--seed", r.get("seed", 42),
                "--size-factors-fit-type", r.get("size_factors_fit_type",
                                                 c.deseq2.get("size_factors_fit_type", "ratio"))],
                [out_dir / f"{row.comparison_id}_rarefaction_summary.tsv",
                 out_dir / f"{row.comparison_id}_rarefaction_stats.tsv"])

    def run_composition(self) -> None:
        c, k = self.config, self.config.composition
        long_table, sample_key, reference = c.need("bc1_counts_long", "sample_key", "bc1_reference_table")
        for key in ("consolidated_library", "yeast_sublibraries", "expected_proportions"):
            if not k.get(key):
                raise StageFailed(f"composition.{key} is required (step 18 has no default)")
        out_dir = c.output_dir / "composition" / k["consolidated_library"]
        args = ["--bc1-counts-long", long_table, "--sample-key", sample_key, "--bc1-reference", reference,
                "--consolidated-library", k["consolidated_library"],
                "--yeast-sublibraries", ",".join(k["yeast_sublibraries"]),
                "--expected-proportions", ",".join(str(x) for x in k["expected_proportions"]),
                "--output-dir", out_dir, "--min-counts", k.get("min_counts", 5)]
        if k.get("filter_condition"):
            args += ["--filter-condition", k["filter_condition"]]
        if k.get("filter_gDNA_plate"):
            args += ["--filter-gDNA-plate", k["filter_gDNA_plate"]]
        self.call("composition", k["consolidated_library"], args, [out_dir / "yeast_sublibrary_proportions.tsv"])


def run_downstream(config: RouteConfig, steps: Sequence[str], *, git_state: Optional[dict] = None,
                   keep_going: bool = False) -> Path:
    """Run `steps` (names from STAGES, executed in STAGES order) and return the run manifest path.

    Raises StageFailed on the first stage that fails or writes nothing; the manifest records how far it got.
    With keep_going, later stages still run (a stage whose inputs a failed stage should have written fails
    on its own), and StageFailed is raised at the end naming every failed stage -- the run is still a failure.
    """
    unknown = [s for s in steps if s not in STAGES]
    if unknown:
        raise ValueError(f"unknown downstream step(s) {unknown}; available: {list(STAGES)}")
    route = Route(config, git_state=git_state)
    failed: Dict[str, str] = {}
    try:
        for name in STAGES:
            if name not in steps:
                continue
            logger.info(f"Stage {STAGES[name][:2]} ({name})")
            try:
                getattr(route, f"run_{name}")()
            except StageFailed as exc:
                failed[name] = str(exc)
                route.manifest["failed_stages"] = failed
                if not keep_going:
                    raise
                logger.error(f"Stage {STAGES[name][:2]} ({name}) FAILED; continuing (--keep-going): {exc}")
        if failed:
            raise StageFailed(f"{len(failed)} stage(s) failed: {', '.join(failed)}")
        route.manifest["status"] = "completed"
    except StageFailed as exc:
        route.manifest["status"] = "failed"
        route.manifest["error"] = str(exc)
        raise
    finally:
        route.manifest["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        route._write_manifest()
    return config.output_dir / "run_manifest.json"
