# HANDOFF — MAGESTIC dual-path program (policy conformance across all modules)

**For:** dedicated implementation agent(s), monitored by Kevin + lead session.
**Policy contract:** `docs/ARCHITECTURE.md` (ratified 2026-06-12) — every module must be
dual-path (direct CLI **and** Snakemake), shared-core, chunked, HPC-array-idempotent, with
**byte-identical** CLI↔Snakemake parity, and no silent success. Read it first.
**Repo:** `software/MAGESTIC` branch `gdb-bdb-harmonized-pipeline` (origin = public k-roy/MAGESTIC).
**Already done this session (committed, unpushed):** R2 replicate fix + multi-plate_key loader
(`0ea8110`), R4 vectorized matcher (`a5043c4`). See `planning/R2_FIX_RESULT_20260611.md`.

---

## 0. Decisions locked (2026-06-12)

- **Scope:** package-wide dual-path program (this doc), not R1-only.
- **wgs:** a **new first-class magestic module** (`src/magestic/pipelines/wgs/`).
- **wgs ↔ TRACE:** `magestic.wgs` **depends on TRACE** (one-directional, pinned). TRACE =
  source-of-truth WGS algorithms (stays standalone); `wgs` = orchestration/integration over it;
  `redi/wgs_check` refactors to **consume** the shared module (one algorithm, one orchestrator).
  Supersedes the prior "not-magestic-depends-on-trace" stance. (`docs/ARCHITECTURE.md §11`.)
- **Parity contract:** **byte-identical** CLI↔Snakemake outputs per module (`§4`).

---

## 1. Gap matrix (current vs policy)

| Module | (A) direct CLI | (B) Snakemake | chunking | HPC array | parity test | Primary gaps |
|---|---|---|---|---|---|---|
| **gdb** `guide_donor_bc0` | ✓ `magestic-gdb-pipeline` | **✗ no Snakefile** | ✓ (vectorized) | partial | ✗ | **add Snakefile**; array-ize; parity test |
| **bdb** `bc1_donor_bc0` | ✓ `magestic-bc1-link` | ✓ Snakefile+13 | ✓ | partial | ✗ | array uniformity; parity test; 2x150 gDNA-guard limitation (documented) |
| **bc1** `bc1_from_screen` | **stub** (R1) | ✓ Snakefile.rules+19 | ✓ | ✓ deseq2 array + `run_full_pipeline.sh` | ✗ | **wire direct CLI (R1)**; **R3 hit-caller**; parity test; reconcile w/ `run_full_pipeline.sh` |
| **redi** | ✓ umbrella+per-step | ✓ Snakefile | ? audit | ✓ array CLIs | ✗ | refactor `wgs_check` to consume `magestic.wgs`; chunking audit; parity test |
| **wgs** *(new)* | ✗ build | ✗ build | ✗ build | ✗ build | ✗ build | **create module** depending on TRACE; full policy build-out |
| pool_edit_characterization | ✗ skeleton | ✗ | ✗ | ✗ | ✗ | out of scope this round (note) |

---

## 2. Per-module work items

### WI-1 · gdb — add the Snakemake path (close the dual-path gap)
- Author `snakemake/Snakefile` over the existing numbered scripts (01-08); the direct
  `run_pipeline` already works. Ensure both call identical `core/` (the just-vectorized
  `process_counts_file` etc.).
- Array-ize the per-library matching as idempotent array tasks (skip-if-matched-output-exists).
- Add `tests/test_dualpath_parity.py`: byte-identical CLI vs Snakemake on a small fixture pool+counts.

### WI-2 · bc1 — wire the direct CLI (R1) + reconcile hit-calling (R3)
- **R1:** replace the 7 stub step fns in `pipeline/run_pipeline.py` so the CLI actually runs.
  Architecture = **orchestrate the numbered `01..18` scripts** (the validated executor) in
  dependency order, reconciling with the existing `run_full_pipeline.sh` (already does
  start/end-step + SLURM/local) — fold that bash orchestrator's behavior into `run_pipeline`
  (or have it shell out), don't reinvent. Include the count-building steps **02–06, 09**
  (esp. 05/06 error-absorption — they change counts), not just the 7 currently-listed steps.
  **No silent success** (§7): fail loud on missing inputs.
- **R3 (BLOCKED on Kevin's gate decision):** the live `14_hit_calling.py` gates on the per-tier
  synonymous-null z ONLY (no padj/FDR); the padj-gated `core/hit_calling.call_hits_single_comparison`
  is never called; `--use-ihw` is a swallowed no-op (`apply_ihw_to_dataframe` undefined in
  `core/ihw.py`). Regardless of the gate choice: **fix the IHW import**, **unify to ONE
  hit-caller**, wire/justify `variance_calibration`. Gate logic waits on D1 (see §4).
- Parity test (byte-identical CLI vs Snakemake) + validation against an existing
  Snakemake-produced screen run as the oracle (`20250912_SpG_QTL_screen`).
- Full detail: `planning/HANDOFF_R1_screen_cli_wiring_20260612.md` (now a sub-doc of this program;
  its "D0 pick-one" framing is SUPERSEDED — the answer is dual-path parity, orchestrate the scripts).

### WI-3 · wgs — create the new module (depends on TRACE)
- Scaffold `src/magestic/pipelines/wgs/` per the §1 module layout. `core/` is **thin over
  TRACE**: import TRACE's WGS editing-outcome functions (donor_outcome 4-Levenshtein,
  wgs-backend); declare TRACE as a **pinned one-directional dependency** in `pyproject.toml`
  (`magestic[wgs]`). TRACE never imports magestic.
- Dual-path: `_cli.py` (`magestic-wgs[-step]`) + `pipeline/run_pipeline.py` + `snakemake/`.
- magestic-side glue = the orchestration + integration not in TRACE: per-build **reference
  design** (T2T-per-build, yKR genealogy → offsets; see memory `wgs_t2t_reference_and_sv_outgrowths`,
  `wgs_per_dataset_reference_design`), lineage ledgers, array workflows over BAMs, REDI hand-off.
- **Scope boundary:** TRACE = algorithms; `wgs` = orchestration/refs/lineage over TRACE.
- Parity + no-silent-success + array tests. **First sub-task: confirm TRACE's installable
  API surface** (what functions `wgs` will call) before building — TRACE is on branch
  `wgs-backend`; pin a commit.

### WI-4 · redi — consume the shared wgs orchestrator
- Refactor `redi/core/wgs_check.py` + `redi/readers/wgs_outcome.py` to **call** `magestic.wgs`
  / TRACE rather than keep an independent WGS implementation (one algorithm home).
- Audit redi steps for chunking; add the parity test.

### WI-5 · cross-cutting — chunking + array uniformity + parity-test harness
- A shared `tests/` helper for the byte-identical dual-path parity check (reused by every module).
- A shared array-sbatch template (encodes the `sherlock-sbatch` traps) + an idempotent
  skip-check convention used by all long steps.
- Audit every module's heavy steps for the chunking + array-idempotency requirements (§5, §6).

---

## 3. Agent decomposition (after policy review)

Parallelizable; one module per agent (worktree isolation where they touch the same files):
- **Agent G** → WI-1 (gdb Snakefile + array + parity).
- **Agent B** → WI-2 R1 (bc1 direct CLI orchestration + no-silent-success + parity). *(R3 gate
  held; do everything except the final gate logic.)*
- **Agent W** → WI-3 (wgs module scaffold + TRACE dependency + dual-path). *(Highest unknowns;
  start with the TRACE API-surface scouting sub-task and report before building.)*
- **Agent R** → WI-4 (redi consume + audit) — sequence after W's wgs API stabilizes.
- **Agent X** → WI-5 (shared parity/array harness) — do FIRST; G/B/W depend on the parity helper.
- Lead session reviews each module's parity + no-silent-success tests before merge; nothing
  pushed to public origin without Kevin (carryover item C).

**Suggested order:** X (harness) → G + B + W in parallel → R → integration. R3 gate + any
push gated on Kevin.

---

## 4. Open decisions still needed from Kevin

- **D1 (R3 hit gate):** per-tier synonymous-null z ONLY, **or** z AND padj/FDR, **or** adopt the
  padj-gated core caller. Blocks only the final `hitcall` gate logic (WI-2). (Options + file:lines
  in `planning/HANDOFF_R1_screen_cli_wiring_20260612.md §4`.)
- **TRACE pin:** which TRACE commit/branch (`wgs-backend`) `magestic[wgs]` pins to (WI-3).
- **pool_edit_characterization:** in or out of this program (currently a skeleton; proposed: out,
  noted as tech-debt).

## 5. Done criteria (program-level)

- [ ] `docs/ARCHITECTURE.md` committed; every module's README links it.
- [ ] Every in-scope module passes: dual-path **byte-identical** parity, no-silent-success,
      all-branch unit, one real-data regression.
- [ ] gdb has a Snakefile; bc1 has a working direct CLI; wgs exists as a TRACE-backed module;
      redi consumes the shared wgs; chunking + array-idempotency uniform.
- [ ] Path-clean diffs (no `/oak|/scratch|larsms|kevinroy`); commits HELD for Kevin's push.

---

_Render date 2026-06-12. Built on the package scan + the four ratified decisions (scope,
wgs-as-module, wgs→TRACE dependency, byte-identical parity). Supersedes the standalone R1
handoff's orchestration framing._
