"""Route A harness: run the real bc1_donor_bc0 stage CLIs over a synthetic library and
score the output against planted truth.

Stages are run as SUBPROCESSES with MAGESTIC_COMMON_DIR set, because
`path_utils.HARMONIZED_VARIANTS_FILE` is bound at import and `oligo_annotations._load_cached`
is lru-cached -- in-process patching is fragile.
"""

from __future__ import annotations

import os
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

REPO = Path(__file__).resolve().parents[3]
SNK = REPO / "src" / "magestic" / "pipelines" / "bc1_donor_bc0" / "snakemake"
PY = sys.executable

import synth_bc1_library as S  # noqa: E402  (same directory)


def _env(common: Path) -> dict:
    e = dict(os.environ)
    e["MAGESTIC_COMMON_DIR"] = str(common)
    e["MAGESTIC_BASE_DIR"] = str(common.parent)
    e["PYTHONPATH"] = f"{REPO / 'src'}:{REPO}"
    e["PYTHONDONTWRITEBYTECODE"] = "1"
    return e


def run(cmd: List[str], common: Path, log: Path, check: bool = True) -> subprocess.CompletedProcess:
    p = subprocess.run([str(c) for c in cmd], env=_env(common), capture_output=True, text=True)
    log.write_text(f"$ {' '.join(str(c) for c in cmd)}\n--- stdout ---\n{p.stdout}\n"
                   f"--- stderr ---\n{p.stderr}\n--- rc {p.returncode} ---\n")
    if check and p.returncode != 0:
        raise RuntimeError(f"stage failed rc={p.returncode}: {' '.join(str(c) for c in cmd)}\n"
                           f"{p.stdout[-3000:]}\n{p.stderr[-3000:]}")
    return p


# ---------------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------------

CLASSES = ["exact_oligo", "same_variant_other_guide", "wrong_variant",
           "ambiguous_incl_truth", "ambiguous_excl_truth", "missing"]


@dataclass
class Score:
    scope: str
    denominator: int
    counts: Dict[str, int] = field(default_factory=dict)
    spurious_bc1: int = 0
    spurious_examples: List[str] = field(default_factory=list)
    detail: Dict[str, str] = field(default_factory=dict)

    @property
    def recovered(self) -> int:
        """Any non-missing outcome: the bc1 came out with SOME oligo call."""
        return self.denominator - self.counts.get("missing", 0)

    @property
    def recovery(self) -> float:
        return self.recovered / max(self.denominator, 1)

    @property
    def precision_oligo(self) -> float:
        """Of RECOVERED links, the fraction that are the planted oligo exactly."""
        return self.counts.get("exact_oligo", 0) / max(self.recovered, 1)

    @property
    def precision_variant(self) -> float:
        """Of RECOVERED links, the fraction correct at EDIT identity."""
        ok = self.counts.get("exact_oligo", 0) + self.counts.get("same_variant_other_guide", 0)
        return ok / max(self.recovered, 1)

    def line(self) -> str:
        c = self.counts
        return (f"[{self.scope}] n={self.denominator} "
                f"recovery={100*self.recovery:.1f}% ({self.recovered}/{self.denominator}) "
                f"prec_oligo={100*self.precision_oligo:.1f}% "
                f"prec_variant={100*self.precision_variant:.1f}% | "
                + " ".join(f"{k}={c.get(k,0)}" for k in CLASSES)
                + f" spurious_bc1={self.spurious_bc1}")


def score_table(df: pd.DataFrame, design: S.Design, scope: str,
                oligo_col: str = "oligo_name") -> Score:
    """Classify every PLANTED bc1 against the stage output. Denominator = planted bc1s."""
    sc = Score(scope=scope, denominator=len(design.links))
    sc.counts = {k: 0 for k in CLASSES}
    if df is None or len(df) == 0 or "bc1" not in df.columns:
        sc.counts["missing"] = sc.denominator
        return sc

    has_oligo = oligo_col in df.columns
    cand_col = "candidate_oligos" if "candidate_oligos" in df.columns else None
    unamb_col = "is_unambiguous" if "is_unambiguous" in df.columns else None

    # one verdict per bc1: the row with the most reads
    cnt = "counts" if "counts" in df.columns else None
    d = df.sort_values(cnt, ascending=False, kind="mergesort") if cnt else df
    best = {}
    for _, r in d.iterrows():
        b = str(r["bc1"])
        if b not in best:
            best[b] = r

    for bc1, planted_names in design.links.items():
        r = best.get(bc1)
        if r is None:
            sc.counts["missing"] += 1
            continue
        got = r.get(oligo_col) if has_oligo else None
        if got is None or pd.isna(got) or got == "":
            sc.counts["missing"] += 1
            sc.detail[bc1] = "present_but_unnamed"
            continue
        got = str(got)
        planted_variants = {design.variant_of(n) for n in planted_names}
        unamb = bool(r.get(unamb_col)) if unamb_col else True
        if not unamb:
            cands = str(r.get(cand_col) or got)
            if any(n in cands for n in planted_names):
                sc.counts["ambiguous_incl_truth"] += 1
            else:
                sc.counts["ambiguous_excl_truth"] += 1
            continue
        if got in planted_names:
            sc.counts["exact_oligo"] += 1
        elif design.variant_of(got) in planted_variants:
            sc.counts["same_variant_other_guide"] += 1
            sc.detail[bc1] = f"same_variant: planted={planted_names} got={got}"
        else:
            sc.counts["wrong_variant"] += 1
            sc.detail[bc1] = f"WRONG: planted={planted_names} got={got}"

    planted = set(design.links)
    extra = [b for b in df["bc1"].astype(str).unique() if b not in planted]
    sc.spurious_bc1 = len(extra)
    sc.spurious_examples = extra[:5]
    return sc


def survival(df: Optional[pd.DataFrame], design: S.Design) -> Tuple[int, int, int]:
    """(planted bc1s present, spurious bc1s present, total rows)."""
    if df is None or len(df) == 0 or "bc1" not in df.columns:
        return 0, 0, 0 if df is None else len(df)
    seen = set(df["bc1"].astype(str))
    planted = set(design.links)
    return len(seen & planted), len(seen - planted), len(df)


# ---------------------------------------------------------------------------------
# Route A
# ---------------------------------------------------------------------------------

@dataclass
class RouteAResult:
    parsed: List[Path]
    collapsed: Path
    linked: Path
    mapped: Path
    reference: Optional[Path]
    stage_frames: Dict[str, pd.DataFrame]
    survival: Dict[str, Tuple[int, int, int]]
    logs: Dict[str, Path]


def run_route_a(design: S.Design, mols: List[S.Molecule], workdir: Path, *,
                n_samples: int = 2, seed: int = 11, gzipped: bool = False,
                header_style: str = "plain",
                gdb_kwargs: Optional[dict] = None,
                run_05: bool = True,
                extra_04: Optional[List[str]] = None,
                tails: Optional[List[str]] = None,
                sps_mapping: bool = False) -> RouteAResult:
    """Split `mols` across `n_samples` PCR replicates, then 01 -> 02 -> 03 -> 04 -> 05."""
    workdir.mkdir(parents=True, exist_ok=True)
    common = workdir / "common"
    S.write_common_dir(design, common)
    logs: Dict[str, Path] = {}
    frames: Dict[str, pd.DataFrame] = {}
    surv: Dict[str, Tuple[int, int, int]] = {}

    # --- 01, once per sample. Every molecule is emitted into EVERY sample so each link has
    #     n_samples PCR replicates (02's support_tier threshold is 2).
    parsed: List[Path] = []
    samples = []
    pdir = workdir / "parsed"
    pdir.mkdir(exist_ok=True)
    for i in range(n_samples):
        tail = (tails[i] if tails and i < len(tails) else S.TAIL_C)
        smols = [S.Molecule(**{**m.__dict__, "tail": tail}) for m in mols]
        reads = S.emit_reads(smols, seed=seed + i)
        ext = ".fastq.gz" if gzipped else ".fastq"
        r1 = workdir / f"s{i}_R1{ext}"
        r2 = workdir / f"s{i}_R2{ext}"
        S.write_fastq_pair(reads, r1, r2, gzipped=gzipped, header_style=header_style)
        out = pdir / f"synth_sample_{i}_bc1_donor_bc0.tsv"
        logs[f"01_s{i}"] = workdir / f"log_01_s{i}.txt"
        run([PY, SNK / "01_parse_bc1_donor_bc0.py", r1, r2, out], common, logs[f"01_s{i}"])
        parsed.append(out)
        samples.append({
            "sample_id": f"synth_sample_{i}", "include": "True", "raw_r1": str(r1),
            "raw_r2": str(r2), "r1_bytes": r1.stat().st_size, "parsed_path": str(out),
            "sequencing_date": "20261001", "run_id": "SYNTH", "sample_number": str(i + 1),
            "outer_primers": "KR1967-KR1590", "gDNA_plate_name": f"synth_plate_{i}",
            "inner_primers": "KR1970-KR1965", "plate_polymerase": "Q5",
            "condition": "plate_wash", "amplicon_level": "pooled_yL",
            "yeast_library": "yLsynth", "source_rule": "synthetic",
        })
    frames["01"] = pd.concat([pd.read_csv(p, sep="\t", dtype=str) for p in parsed],
                             ignore_index=True)
    surv["01_parse"] = survival(frames["01"], design)

    # --- 02
    man = S.write_manifest02(samples, workdir / "manifest_02.tsv")
    collapsed = workdir / "collapsed.tsv"
    purity = workdir / "purity.tsv"
    logs["02"] = workdir / "log_02.txt"
    run([PY, SNK / "02_combine_and_collapse.py", "--manifest", man,
         "--output-collapsed", collapsed, "--output-purity", purity,
         "--n-jobs", "1", "--accumulate-every", "1"], common, logs["02"])
    frames["02"] = pd.read_csv(collapsed, sep="\t", low_memory=False)
    surv["02_collapse"] = survival(frames["02"], design)

    # --- 03
    gdb = workdir / "gdb"
    S.write_gdb_dir(design, gdb, **(gdb_kwargs or {}))
    linked = workdir / "linked.tsv"
    logs["03"] = workdir / "log_03.txt"
    run([PY, SNK / "03_link_guide_donor_bc0.py", "--input", collapsed,
         "--output", linked, "--guide-dir", gdb], common, logs["03"])
    frames["03"] = pd.read_csv(linked, sep="\t", low_memory=False)
    surv["03_link_guide"] = survival(frames["03"], design)

    # --- 04
    mapped = workdir / "mapped.tsv"
    logs["04"] = workdir / "log_04.txt"
    cmd = [PY, SNK / "04_map_to_oligos.py", "--input", linked, "--output", mapped]
    if sps_mapping:
        cmd += ["--sps-mapping", S.write_sps_mapping(design, workdir / "sps_map.tsv")]
    cmd += extra_04 or []
    run(cmd, common, logs["04"])
    frames["04"] = pd.read_csv(mapped, sep="\t", low_memory=False)
    surv["04_map_oligos"] = survival(frames["04"], design)

    # --- 05
    reference = None
    if run_05:
        reference = workdir / "bc1_reference_table.tsv"
        logs["05"] = workdir / "log_05.txt"
        p = run([PY, SNK / "05_assign_sublibrary.py", "--input", mapped,
                 "--sps-mapping", S.write_sps_mapping(design, workdir / "sps_map.tsv"),
                 "--output-reference", reference,
                 "--output-summary", workdir / "summary_05.tsv"],
                common, logs["05"], check=False)
        if p.returncode == 0 and reference.exists():
            frames["05"] = pd.read_csv(reference, sep="\t", low_memory=False)
            surv["05_assign"] = survival(frames["05"], design)
        else:
            reference = None

    return RouteAResult(parsed=parsed, collapsed=collapsed, linked=linked, mapped=mapped,
                        reference=reference, stage_frames=frames, survival=surv, logs=logs)
