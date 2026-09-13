# HANDOFF — MAGESTIC gdb / bdb pipeline (harmonized matcher + bridge) — 2026-06-08

Scope: get all gdb (guide-donor-bc0) data mapped through the package CLI
(`magestic-gdb-pipeline`), cross-era validate (SpCas9/LbCas12a/SpG), iron out bugs,
wire the bdb (bc1-donor-bc0) `06` annotate-unmatched into the CLI DAG, and test the
bdb→gdb bridge. Companion findings doc (numbers + tables):
`/oak/.../projects/QTL/gdb_pipeline_findings_20260608.md` (mirrored in Drive
`planning/`). Auto-memories: `gdb_reclone_scope_policy`, `bdb_bridge_validation`.

All pipeline edits are **[uncommitted]** working-tree changes in the git repo
`/oak/stanford/groups/larsms/Users/kevinroy/software/MAGESTIC` (editable-installed, so
edits are live). Each has a dated `.bak_*` sibling for rollback. Nothing committed yet —
a fresh agent should review the 5 modified files and commit (or revert) deliberately.

## 1. What was done

- **gdb harmonized matcher** — `oligo_matching.py` `load_oligo_pool_harmonized()` +
  `config.py` (`use_harmonized_table`, `harmonized_v_libraries`,
  `harmonized_nuclease_types`, `library_scope`) + `run_pipeline.py`
  `run_oligo_matching` per-library scoped path + CLI flags (`--use-harmonized-table`,
  `--nuclease`, `--v-libraries`, `--library-scope-json`). Fixes LbCas12a/impLbCas12a
  extraction (per-nuclease guide/donor from the harmonized table vs the broken fixed
  200mer slice) and lets each library be scoped to its own design sub-pool.
- **gdb match-step glob fix** — excludes its own `*_matched_*` outputs (was re-ingesting
  them → `_matched_matched_` junk + doubled runtime).
- **gdb represent/plots harmonized mode** — `_canonical_twist_designs()` routes Step 05/06
  to the canonical Twist designs (`20210422_Twist_200mer.tsv` + `20240411_...`) which
  carry `oligo_seq`→SPS; gracefully skips missing-oligo fill when no `library_sps_map`.
- **bdb STEP 6 wired** — `--annotate-unmatched` runs `06_annotate_unmatched_donors.py`
  (+ `--integrate-annotations` → `07`) after STEP 5, using the **MAGESTIC_cs design
  reference** (`common/reference_genomes/MAGESTIC/MAGESTIC_background_strain.fasta` +
  `common/annotation_files/MAGESTIC_background_strain_annotations.gff`), NOT the yKR1033
  WGS assembly. Fixed 2 stale template paths (REF `/MAGESTIC/` subdir; SpCas9 design
  `20210422_Twist_200mer.tsv`).
- **bdb bridge guide-preservation fix** (the core bridge bug) — `aggregate_and_filter`
  groupby was dropping the `guide` the gdb bridge sets (`merge.py record.guide`), so
  STEP 4 matched donor-only and `guide_match` was False for every bc1. Now preserves
  `guide`/`donor_from_gdb` with `'first'`.
- **bdb merge NaN guards (6 sites)** in `merge.py` — `build_donor_bc0_key` `len(float)`;
  5 truthy-on-NaN donor-field checks; `donor_portion + self.bc0` concat. NA-donor/bc0
  rows (present in matched gdb tables + 2x300) crashed the merge.
- **All gdb pools mapped** (scratch fresh dirs): 20240404 reclones (V629/631/634),
  2021 pool (11 libs), yL406 (V545/546/547 originals + reused V429/V454).
- **bdb→gdb bridge + STEP 6 validated** on yL406 (debug 4-sample).

## 2. What's verified

- **gdb cross-era vs March oracle** (read-level from matched count tables):
  SpCas9 V429 75.4% perfect = oracle 75.4% (**exact parity**); LbCas12a V454
  **0.10% → 54.5% perfect**; impLbCas12a never matched in March → now 54.5%; SpG V629
  perfect 69.3% = oracle 69.3%; V634 (SpCas9 reclone) 24.9%(wrong-pool) → 80.4%.
- **Per-library scope isolation**: V634 lookup = 16,878 oligos (SpCas9 only), no SpG leak
  (lookup sizes logged per library).
- **represent harmonized mode**: 2021 pool Step 05 COMPLETED, 11 libs, 140,887
  oligo-library rows (`reptest_*` job).
- **bdb full pipeline + STEP 6**: debug run exits clean (~2 min); `06` triple-aligner →
  `unmatched_donor_annotations.tsv` → `07` → `final/bc1_reference_table_with_unmatched_annotations.tsv`.
- **bridge correctness at EDIT identity** (chr_pos_ref_alt, debug 4-sample vs oracle):
  perfect 91.3%, sliding_window 95.6%, guide_only 100%, ALL 82.8% (ambiguous categories
  lower, expected). guide_match 0 → 4,499 True after the guide-preservation fix.
- **NOT VERIFIED:** full-scale (all-sample, collapse-on) bdb oracle number — the full
  STEP 1-5 OOMs at 128GB (61M 2x300 + 44M gdb rows). Mechanism validated on debug only.
- **NOT VERIFIED:** cross-library/cross-era contamination (Task #4) — not started.

## 3. Open items

- **Full-scale bridge oracle number** — deferred. Why: full STEP 1-5 OOMs at 128GB; the
  *mechanism* is already validated at edit level on the debug subset, so this is
  completeness-only. To do it: dedup the gdb bridge reference to unique
  `(donor_bc0_fragment → oligo_name, guide)` before the merge (collapses ~44M rows to
  ~design size — `load_gdb_data` in `bc1_donor_bc0/readers/data_loader.py` is where the
  44M rows enter) and/or run on bigmem (`--partition bigmem --mem 240G`).
- **Task #4 — contamination review** — not started. Note: the original framing
  ("V629/631/634 contamination in the 20210422 dir") is GARBLED — the 20210422 dir is a
  V455-only stub; V629/631/634 live in 20240404. Proper version: match each library
  against *other* libraries' scopes and quantify cross-bleed. Why deferred: fresh
  investigation, better in a clean context.
- **Task #6 — bc0-conditioned synthesis-error QC report** — core question answered
  (original V545 = 77.8% perfect > reclone V629 69.3%; the "91.62%" was the
  high-abundance figure; gap = low-abundance synthesis-error tail). The per-bc0
  perfect-vs-1-edit-vs-N-edit breakdown report is the remaining deliverable; overlaps
  with the bdb `06` alignment annotation now wired. Why deferred: nice-to-have.
- **Commit the 5 pipeline edits** — all [uncommitted] in the MAGESTIC repo. Why deferred:
  left for deliberate review (this session validated behavior, not a code review).
- **Scratch outputs purge** — everything under `/scratch/users/kevinroy/gdb_pilot/` is
  reproducible (re-run the sbatch scripts there) but will purge ~90 days untouched.

## 4. Resume command

**Resume:** nothing is running (`ssh sherlock 'squeue -u kevinroy'` → only the
code-server). First decide on the 5 uncommitted edits: `ssh sherlock 'cd
/oak/.../software/MAGESTIC && git status -s && git diff --stat'`. If the diff looks
right, commit (gdb + bdb fixes are independent — consider two commits). Then pick ONE
open item:
- For the **full-scale bridge number**: edit `load_gdb_data`
  (`bc1_donor_bc0/readers/data_loader.py`) to `drop_duplicates` the gdb frame on
  `donor_bc0_fragment` (keep `oligo_name`,`guide`) before returning, then resubmit
  `/scratch/users/kevinroy/gdb_pilot/run_bridge_full.sbatch` with `--mem 240G`; when
  COMPLETED, re-run `/scratch/users/kevinroy/gdb_pilot/cut_edit_identity.py` (compares at
  edit identity — DO NOT use raw donor/oligo_name equality, it gives a false ~8%).
- For **Task #4 contamination**: take a matched gdb table (e.g. `run_yl406_fresh` V545)
  and re-match its reads against a different scope (e.g. `[V429]`); a meaningful match
  fraction = contamination.
All diagnostic/run scripts are in `/scratch/users/kevinroy/gdb_pilot/` (named
`run_*.sbatch`, `cut_*.py`, `diag_*.py`, `compare_*.py`, `patch_*.py`).

## 5. Files touched

In `/oak/.../software/MAGESTIC/src/magestic/pipelines/` — all **[uncommitted]**, each
with a dated `.bak_*` sibling:

- `guide_donor_bc0/core/oligo_matching.py` — `load_oligo_pool_harmonized()`
  (`.bak_20260607_gdb_harmonized`)
- `guide_donor_bc0/config.py` — harmonized/scope config fields (`.bak_20260607_gdb_harmonized`)
- `guide_donor_bc0/pipeline/run_pipeline.py` — harmonized scoped matching + CLI flags;
  glob self-ingestion fix; represent/plots harmonized fallback
  (`.bak_20260607_gdb_harmonized`, `_globfix`, `.bak_20260608_representfix`)
- `bc1_donor_bc0/pipeline/run_pipeline.py` — STEP 6 (06+07) wiring + path resolution;
  guide-preservation in `aggregate_and_filter` (`.bak_20260607_step6`,
  `.bak_20260608_preserveguide`)
- `bc1_donor_bc0/core/merge.py` — 6 NaN guards (`.bak_20260608_nanguard`, `_nanguard2`,
  `_bc0guard`)

Not git (durable outputs):
- `/oak/.../projects/QTL/gdb_pipeline_findings_20260608.md` — findings doc (+ Drive mirror)
- `/scratch/users/kevinroy/gdb_pilot/` — fresh-dir run outputs (reproducible) + all
  sbatch/diagnostic scripts
