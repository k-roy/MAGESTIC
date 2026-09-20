# MAGESTIC — Package Architecture & Contribution Policy

Canonical policy for every pipeline module in `magestic`. New code and refactors MUST
conform; deviations require an explicit, documented exception. (Status: ratified
2026-06-12. Supersedes the previously-implicit convention.)

---

## 1. Module model

Every pipeline lives under `src/magestic/pipelines/<module>/` with a fixed layout:

```
<module>/
  _cli.py            # console_script entry point(s) -> pipeline/run_pipeline
  config.py          # dataclass config (+ from_yaml, factory presets); paths via config, never hardcoded
  core/              # ALL analysis logic, as pure library functions (no I/O orchestration)
  readers/           # input loaders / output writers
  pipeline/
    run_pipeline.py  # DIRECT path: orchestrates core/steps in dependency order
  snakemake/
    NN_*.py          # numbered step scripts (thin CLI wrappers over core/)
    Snakefile[.rules]# SNAKEMAKE path: rules over the numbered scripts
  qc/                # QC + provenance
  tests/             # unit (all-branch) + dual-path parity + no-silent-success
```

Current peer modules: **gdb** (`guide_donor_bc0`), **bdb** (`bc1_donor_bc0`),
**bc1** (`bc1_from_screen`), **redi**, **wgs** (new — see program plan).

---

## 2. Dual-path requirement (the core contract)

Every module MUST expose **both**:
- **(A) a direct CLI** — a pip-installed `console_script` (`magestic-<module>[-<step>]`) driving
  `pipeline/run_pipeline.py`; runnable with no Snakemake install.
- **(B) a Snakemake workflow** — a `Snakefile` over the numbered step scripts.

Neither path may contain analysis logic the other lacks. A module with only one path is
**incomplete**, not "done." (Anti-patterns today: `bc1` direct path is a stub; `gdb` has no
Snakefile.)

## 3. Shared-core principle

All analysis logic lives in `core/`. Both orchestrators (the direct `run_pipeline` and the
Snakemake numbered scripts) are **thin wiring** over identical `core/` functions — they sequence
steps and move files, nothing more. A numbered script is a CLI wrapper that parses args, calls
the same `core/` function the direct path calls, and writes the same output. **No logic is
duplicated** between a script and `core/`, or between the two orchestrators.

## 4. Byte-identical parity contract

For the same inputs + config, path (A) and path (B) MUST produce **byte-identical output
files**. Enforced per module by `tests/test_dualpath_parity.py` (run both paths on a small
fixture, `diff` every output). This is achievable because both paths call identical `core/`
functions (§3). Sole allowed carve-out: an external nondeterministic dependency (e.g.
`pydeseq2`) must be **seeded/pinned** so output is still byte-identical across paths; if a
column is genuinely unpinnable, it is excluded from the diff with an explicit, reviewed
comment naming the source of nondeterminism. (Precedent: R2/R4 were validated byte-identical.)

## 5. Chunking

Any step whose memory scales with read/row count MUST process in **bounded chunks**
(`chunk_size` parameter; stream via `pd.read_csv(..., chunksize=)` or equivalent) — never an
unbounded load-all-into-memory. Each step documents its peak-memory target. (Precedent:
`process_counts_file` chunked + vectorized.)

## 6. HPC array-idempotency

Any step that can exceed ~20 min MUST be decomposable into **idempotent SLURM-array tasks**:
- each array task **skip-checks its own output at start** (re-run safe; no work if output exists),
- target 20–60 min per task (preemptible `owners`-style queues reclaim nodes; monolithic jobs lose all WIP),
- ship an array `sbatch` template (per the `sherlock-sbatch` skill traps: account, `set -u`
  after `source ~/.bashrc`, AVX-512 constraint, env propagation),
- both orchestrators support `--start-step/--end-step` (or `--from-step`) and re-run idempotently.

## 7. No silent success

Orchestrators MUST exit **non-zero** on a missing input or failed step, and MUST NOT log
"completed successfully" without having produced outputs. (Anti-pattern: the bc1 stub returns
`True` from every step.) A `tests/` case asserts loud failure on bad input.

## 8. Path & provenance hygiene (public repo)

Origin is the public `github.com/k-roy/MAGESTIC` (history-scrubbed). No real paths, usernames,
or accounts in code: paths come from config/env/CLI args; defaults are placeholders
(`/path/to`, `account=mylab`). Every data output carries a `PROVENANCE.json`; every figure a
provenance footer (global figure policy). CI/pre-commit greps the diff for internal
cluster paths, usernames, and account names before any push.

## 9. CLI conventions

`console_scripts` named `magestic-<module>` (umbrella) and `magestic-<module>-<step>` (per-step,
the REDI pattern). Standard flags: `--config <yaml>`, `--project-dir`, `--steps a,b,c`,
`--from-step`/`--start-step`/`--end-step`, `--use-scratch`, `--dry-run`, `--debug`.

## 10. Testing (per module, required before "done")

- **Unit:** all-branch coverage of `core/` (every match/decision branch).
- **Dual-path parity:** §4 (byte-identical A vs B on a fixture).
- **No-silent-success:** §7 (loud failure on bad input).
- **Regression:** at least one real (or realistic) dataset reproduces a known-good reference.

## 11. Cross-package dependencies (sibling published packages)

A module MAY depend on a sibling published package (e.g. **TRACE**) when that package is the
source-of-truth for an algorithm, rather than duplicating it:
- the dependency is **one-directional and pinned** (declared in `pyproject.toml`, e.g.
  `magestic[wgs]` → a pinned TRACE); the sibling **never imports magestic** (no cycles),
- the sibling stays **independently installable and usable** on its own,
- the magestic module is the **orchestration/integration layer** (dual-path CLI + Snakemake +
  array + magestic config/glue); the sibling owns the **algorithms** (shared-core, §3),
- **parity (§4) still holds** because both the CLI and Snakemake paths call the *same pinned*
  sibling functions → byte-identical output,
- exactly **one** home per algorithm: do not let a third location keep its own copy — internal
  consumers (e.g. `redi/wgs_check`) call the shared module, they don't re-implement it.

**Canonical instance:** `magestic.pipelines.wgs` → depends on **TRACE** (TRACE = WGS
editing-outcome algorithms / source of truth; `wgs` = orchestration over it; `redi/wgs_check`
is refactored to consume it). This consciously supersedes the earlier
"not-magestic-depends-on-trace" stance, now that `wgs` is a first-class magestic module.

---

_A module is "policy-complete" only when all of §2–§11 hold and the three test classes are green.
The per-module gap matrix + closure plan live in `planning/HANDOFF_dualpath_program_20260612.md`._
