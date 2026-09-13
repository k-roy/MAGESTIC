# HANDOFF — Wire the `bc1_from_screen` CLI (R1) + reconcile hit-calling (R3)

**For:** dedicated implementation agent(s), monitored by Kevin + lead session.
**Repo:** `/oak/stanford/groups/larsms/Users/kevinroy/software/MAGESTIC` (branch
`gdb-bdb-harmonized-pipeline`; origin = public `github.com/k-roy/MAGESTIC`).
**Staged read-only copy (no ssh needed for reading):** `/tmp/magestic_review/bc1_from_screen/`
(re-stage: `rsync -az sherlock:/oak/.../software/MAGESTIC/src/magestic/pipelines/bc1_from_screen/ /tmp/magestic_review/bc1_from_screen/`).
**Companion:** `planning/R2_FIX_RESULT_20260611.md` (R2/R4 already done + committed this session).

---

## 0. Goal & the one thing to confirm first (Decision D0)

Make `magestic-screen-analysis` actually run the screen analysis instead of being a
**silent no-op** (today every step fn just logs and `return True` — it reports
"PIPELINE COMPLETED SUCCESSFULLY" having produced nothing). The CLI must either do the
work or fail loudly; it must never again return success without producing outputs.

**D0 — orchestration architecture (confirm before building):**

| Option | What the CLI does | Pros | Cons | Effort |
|---|---|---|---|---|
| **2b-scripts (RECOMMENDED)** | `run_pipeline` orchestrates the numbered `snakemake/01..18` scripts in dependency order, as subprocesses, with paths derived from `BC1FromScreenConfig`. | Reuses the **validated real executor** (the numbered scripts); no duplicate logic; no Snakemake project-file needed; easiest to validate (same scripts the paper used). | Subprocess orchestration; must map config→each script's argparse; loses Snakemake's free DAG/restart (acceptable — it's a linear-ish DAG). | Medium |
| 2a-core | Each step fn loads inputs, calls `core/` functions directly, writes outputs (pure Python, no subprocess). | Pure-Python, no Snakemake; cleanest for pip users. | **Reimplements** the per-step I/O the numbered scripts already do → highest divergence risk; must re-validate every step against the scripts. | High |
| 2c-snakemake | CLI shells out to `snakemake -s Snakefile.rules` with a generated project Snakefile. | Uses Snakemake's DAG/parallelism/restart. | Needs a project Snakefile + path vars; only 11/18 steps have rules (see §3); adds a Snakemake runtime dep. | Medium-High |

**Lead's recommendation: 2b-scripts.** Rationale: the numbered scripts are the only
*validated* implementation; orchestrating them is the lowest-divergence path to a working,
paper-faithful CLI, and it's the cheapest to validate (outputs must match an existing
Snakemake-produced screen run). Treat `core/` as the library the scripts already use, not a
second execution path. **Do not start building until D0 is confirmed.**

---

## 1. Current state (precise)

- **Entry point** (`pyproject.toml:86`): `magestic-screen-analysis = "magestic.pipelines.bc1_from_screen._cli:main"`
  → `_cli.py:main` → `pipeline/run_pipeline.py:main` (`run_pipeline.py:324`).
- **7 stub step fns** in `run_pipeline.py` (all log-only, `return True`, never call the
  functions they import):
  | CLI step | fn (file:line) | imports (uncalled) | working dir | maps to numbered script |
  |---|---|---|---|---|
  | `tier` | `run_tier_classification` (62-82) | `core.tier_classification`: `calculate_t0_counts, classify_bc1_tiers, create_tier_summary` | `get_working_dir(...)` | `07_tier_classification.py` |
  | `matrix` | `run_count_matrix` (85-101) | `core.count_matrix`: `create_wide_count_matrix, consolidate_t0_samples` | `bc1_counts` | `08_create_count_matrix.py` |
  | `split` | `run_split_comparisons` (104-117) | `core.count_matrix`: `split_by_comparison, create_comparison_manifest` | `counts_for_deseq` | `10_split_comparisons.py` |
  | `deseq2` | `run_deseq2` (120-129) | **none** (comment only) | `DESeq2_consolidated_t0` | `11_deseq2_single.py` (pydeseq2) |
  | `annotate` | `run_annotation` (132-137) | **none** (empty) | — | `13_annotate_results.py` |
  | `hitcall` | `run_hit_calling` (140-159) | `core.hit_calling`: `call_hits_single_comparison, aggregate_hit_calls` | `hit_calling_tiered` | `14_hit_calling.py` + `15_aggregate_hits.py` |
  | `plots` | `run_plots` (162-167) | **none** (empty) | — | (no script; see §6) |
- **Config** (`config.py`): `BC1FromScreenConfig` exposes `get_working_dir(subdir)`,
  `processed_data_dir`, `bc1_counts_dir`, `deseq2_dir`, `hit_calling_dir`, `keyfiles_dir`,
  `libraries` (`LibraryConfig`: name, screen_dates, timepoints, t0_consolidation_groups),
  `hit_calling_config` (per-tier z/concordance), `tier_config`. Pre-defined factory
  `create_spg_qtl_screen_config()` (config.py:364) targets `20250912_SpG_QTL_screen`
  (libs yL437, yL442). `use_scratch` toggles working dirs to `$SCRATCH/magestic/<screen>`.
- **`run_pipeline` driver** (250-321) already has: config.validate(), scratch staging
  (`stage_inputs_to_scratch`/`stage_outputs_to_oak`), the `STEP_FUNCTIONS` loop with
  per-step try/except. **This scaffolding is reusable** — only the step bodies are empty.

---

## 2. The numbered scripts = source of truth

`snakemake/01..18` are the real implementations (each a standalone argparse script). The
AGENTS' first discovery task per step: read the script's argparse + its `core` imports +
its input/output file contract, then wire `run_pipeline`'s step fn to invoke it (2b) or
replicate its core calls (2a). Inventory:

01 validate_collapsed · 02 parse_bc1_reads · 03 sample_qc · 04 combine_counts ·
05 error_detection · 06 absorb_error_counts · 07 tier_classification · 08 create_count_matrix ·
09 combine_technical_replicates · 10 split_comparisons · 11 deseq2_single · 12 slope_fitness ·
13 annotate_results · 14 hit_calling · 15 aggregate_hits · 16 rarefaction ·
17 deseq_rarefaction · 18 validate_sublibrary_composition

**The CLI's 7 steps only cover 07,08,10,11,13,14/15 — the count-building steps 02-06 and
09 are NOT in the CLI's step list.** A real end-to-end CLI must include them (esp. **05/06
error-detection + absorption, which change the counts** → affect every downstream number).
Decide explicitly: does `magestic-screen-analysis` start from raw collapsed FASTQ (02) or
from an existing count matrix? Recommended: support both via `--from-step`, default to the
full chain.

---

## 3. Snakefile.rules coverage (for 2c, and as the dependency DAG reference)

`snakemake/Snakefile.rules` has rules for **02,03,07,08,10,11,12,13,14,15,16** (11/18).
**No rule:** 01,04,05,06,09,17,18. (So 2c would still need these 7 wired.) Even under 2b,
use the rules' input/output declarations as the authoritative **dependency order**.

---

## 4. Decision D1 (BLOCKS the `hitcall` step) — hit-calling reconciliation (R3)

The `hitcall` step cannot be finalized until Kevin resolves R3. Summary for the agent:
- **Live** `snakemake/14_hit_calling.py` gates significance on the **per-tier, mean-centered
  synonymous-null z-score ONLY** (`is_significant = abs(z) >= z_threshold`, ~line 233) — **no
  padj/FDR gate**, although DESeq2 (`11_deseq2_single.py`) does compute `padj`.
- **Dead** `core/hit_calling.py:call_hits_single_comparison` (~line 202) IS padj-gated
  (`padj<0.1 & abs(z)>=z_thr & baseMean>=thr`) but **never called**.
- **IHW is broken:** `14` imports `apply_ihw_to_dataframe` from `core/ihw.py`, which **does
  not define it** (it defines `run_ihw`, `apply_ihw_to_robust_wls_results`); the import sits
  in a `try/except` so `--use-ihw` is a **silently-swallowed no-op**.
- **`core/variance_calibration.py`** is bypassed — `14` imports only the
  `VarianceCalibrationResult` dataclass and reimplements null estimation inline (~40 lines).

**Pending Kevin decision (the significance criterion):** (A) z AND padj/FDR; (B) z-only
(per-tier synonymous null IS the FDR control, padj informational); (C) adopt the padj-gated
core caller as canonical. **Regardless of A/B/C, the agent will:** fix the IHW import (define
`apply_ihw_to_dataframe` wrapping `run_ihw`, or repoint), and **unify to ONE hit-caller** (no
dead second caller). Hold the gate logic until D1 is set; everything else in `hitcall` can
proceed.

---

## 5. Validation plan (the gate for "done" — non-negotiable)

A wired CLI that "runs" is not done; it must **reproduce a known-good reference**.
1. **Establish the reference.** Find the existing Snakemake-produced outputs for a real
   screen (start with `20250912_SpG_QTL_screen`: `projects/QTL/20250912_SpG_QTL_screen/
   processed_data/bc1_from_screen/`). If a clean full run doesn't exist, the lead session
   will produce one via the numbered scripts to serve as the oracle.
2. **Per-step parity.** For each wired step, run the CLI step and **diff its output against
   the reference** (counts/tier/matrix/DESeq2 LFC+padj/hits). Tolerances: exact for
   counts/tiers/matrix; DESeq2 LFC/padj within float tolerance (pydeseq2 determinism);
   hits exact once D1 is fixed.
3. **End-to-end.** Full `magestic-screen-analysis --use-spg-config` (or a small fixture
   screen) reproduces the reference hit list.
4. **No-silent-success test.** Force a missing input → assert the CLI **exits non-zero**
   and does NOT print "COMPLETED SUCCESSFULLY".
5. **Idempotency / `--from-step` / `--steps` subsetting** behave (re-run = same outputs).

Heavy steps run via `sbatch` (sherlock-sbatch skill; account `larsms`, AVX-512 constraint,
`set -u` after `source ~/.bashrc`). pydeseq2/rpy2 env: use
`/home/groups/larsms/users/kevinroy/anaconda3/bin/python` (confirm pydeseq2 present; IHW
needs rpy2+R — gate behind availability check, never hard-fail the whole run).

---

## 6. Known traps / notes

- **pydeseq2** (deseq2 step) and **rpy2/R** (IHW) are optional-heavy deps — detect-and-skip
  with a loud warning, never a silent pass.
- **`plots` step** has no numbered script and no core module — either implement minimal QC
  plots (run the `figure-qa` skill on any figure) or drop `plots` from `AVAILABLE_STEPS`
  rather than ship an empty stub. Decide explicitly.
- **`scaffold`/path hardcodes**: `config.base_dir` default is the scrubbed `/path/to`
  placeholder (public-repo safe). Real runs pass `--project-dir` or a YAML. Don't reintroduce
  real paths (origin is public; the repo was history-scrubbed — see memory
  `public_repos_scrubbed`).
- **Sherlock ControlMaster discipline:** serial ssh only, never close the master, heavy work
  → sbatch (see project CLAUDE.md). Background long jobs; watch via one completion-loop.
- **~5 orphan core modules** (`causality, empirical_thresholds, library_calibration,
  noise_assessment, sublibrary`, ~2,100 LOC) are invoked by nothing — out of scope here;
  wire or delete in a separate pass (note, don't touch).

---

## 7. Suggested agent decomposition (after D0 confirmed)

- **Agent A — count-building chain (02-09):** wire/validate the steps that build the count
  matrix (incl. 05/06 error absorption). Highest correctness leverage.
- **Agent B — stats chain (10,11,13):** split → DESeq2 → annotate; validate LFC/padj parity.
- **Agent C — hitcall + R3 (14/15):** unify the hit-caller, fix IHW, implement D1's gate
  (gated on Kevin's decision); validate hit-list parity.
- **Agent D — CLI shell + tests:** `run_pipeline` orchestration (D0 architecture),
  `--from-step`/`--steps`, the no-silent-success test, end-to-end harness.
Agents A/B can run in parallel; C waits on D1; D integrates last. Lead reviews each before merge.

## 8. Done criteria

- [ ] D0 + D1 confirmed by Kevin.
- [ ] `magestic-screen-analysis` reproduces the reference screen's hit list end-to-end.
- [ ] Per-step parity diffs green (§5).
- [ ] No-silent-success test passes (fails loudly on bad input).
- [ ] One hit-caller; IHW import fixed (no swallowed no-op); `--use-ihw` either works or is removed.
- [ ] `plots` either implemented (figure-qa clean) or removed from steps.
- [ ] Commit(s) on `gdb-bdb-harmonized-pipeline`, path-clean (no `/oak`,`/scratch`,`larsms`,
      `kevinroy`), HELD for Kevin's push.

---

_Render date 2026-06-12. Source: lead session R1/R3 mapping (Explore passes) + direct reads
of run_pipeline.py / config.py / pyproject.toml. R2 + R4 already committed (0ea8110, a5043c4)._
