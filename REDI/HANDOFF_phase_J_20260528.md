# Phase J — REDI software submodule migration — HANDOFF (2026-05-28)

**Location:** `/path/to/software/MAGESTIC/REDI/HANDOFF_phase_J_20260528.md`
**Author:** Opus agent (resumption of stalled 2026-05-28 agent)
**Parent plan:** `/oak/.../QTL/MASTER_PLAN.md` §9d

---

## 1. State on disk at session start

Prior agent (stalled mid-write) left **scaffolding directories only** (all empty except one `__init__.py`):

```
software/MAGESTIC/REDI/
  configs/  data/  docs/  plate_database/  reference/
  screens/  snakemake/  tests/                  (all empty)

software/MAGESTIC/src/magestic/pipelines/redi/
  __init__.py    (architecture docstring + __all__ = ["config","core","readers","pipeline","qc"])
  core/  readers/  pipeline/  qc/  scripts/  slurm/  snakemake/   (all empty)
```

The `__init__.py` docstring is the only authoritative signal from prior intent:
> Operates on paired-end reads from REDI mating plates and yields:
> - Per-(REDI_bc, MAGESTIC_bc1) read-count tables
> - Colony purity scores
> - Cross-array (SHA435 vs yT182/yT183) agreement metrics
> - Annotation against bc1-donor-bc0 reference table
> - Type-(i) bc1↔edit (WGS) and type-(ii) REDI-bc↔MAGESTIC-bc1 (re-REDI) via unified `core.correspondence`

Architecture intent matches MASTER_PLAN §9d. Building on that.

---

## 2. Verified facts

| Item | Status |
|---|---|
| Canonical pipeline source `20250312_NNS_complex_satmut/scripts/step_a–step_k` | EXISTS — 17 step files (a through k); sizes 1–15 KB |
| `20240305_RBYG_REDI/scripts/` | EXISTS but EMPTY (`raw_data/`, `keyfiles/`, `processed_data/` populated); per OQ-22, scripts never written there |
| **`20240305_RBYG_REDI/keyfiles/` is EMPTY** | **BLOCKER for J.4 step 2** — can't re-run without keyfiles. See §5. |
| `20250312_NNS/keyfiles/` | EXISTS — `REDI_plate_barcode_key.tsv`, `20230513_updated_REDI_bc_V5_key.tsv`, `sample_key.tsv`, `gDNA_plate_key.tsv`, `WGS_plate_key.tsv`, `20250615_yT182_yT183_keyfile.tsv`, `1536_PIXL_rearray_format.tsv`, `MAT_alpha_array_1536_to_384_key.tsv`, `384_PIXL_rearray_format.tsv`, `remaining_NRD1_hits.tsv`, `agar_phenotyping/` |
| Shengdi WGS table | EXISTS at `/oak/.../collaborator/projects/WGS_MAGESTIC_QTL_bc_calling/output3/summaries/wgs_outcome_full_table.v0.tsv` (also `v20250701.tsv`) |
| NRD1 README pattern | Captured; per-stage step_a–step_i naming, `.sh` + `.py` pairs, `data/` for processed artifacts |
| `bc1_donor_bc0/` engine pattern | `core/`, `readers/`, `pipeline/`, `qc/`, `scripts/`, `snakemake/`, `templates/`, `_cli.py`, `config.py`, `__init__.py` — REDI engine should mirror this |
| `pyproject.toml` mtime | May 8 — **no concurrent edits from WGS B.2 agent visible**. Safe to edit (re-check mtime immediately before edit). |
| Plate DB master xlsx | EXISTS at `/oak/.../lab_shared/SGTC_genome_editing_group/Rectangular_plate_barcode_database_for_PIXL_Spinnaker.xlsx` (7 sheets per spec) |
| Git state | `src/`, `tests/`, `pyproject.toml`, `docs/`, `recipes/`, `reproduce/`, `dist/`, `.github/`, `CHANGELOG.md` ALL UNTRACKED at top level — package repo in unusual state. **Do not `git add -A` / `git add .` — add specific paths only.** |

---

## 3. Scope for THIS session (explicitly bounded)

The prior agent crashed trying to write the entire engine in one go. This session writes **a subset**, commits incrementally (~every 10 file writes), and defers the rest with clear instructions for the next agent.

**In scope:**
- J.1 engine scaffolding: `core/{config.py, correspondence.py, purity.py, annotate.py, cross_array.py, wgs_check.py}`, `readers/{fastq.py, keyfiles.py, plate_db.py, wgs_outcome.py}`, `pipeline/{step_a..step_i}_*.py` ported as importable modules, `qc/`, `__init__.py` files
- J.2 plate DB snapshot tool: `REDI/plate_database/export_plate_db.py` (xlsx → TSV + parquet); **does NOT modify the live xlsx**
- CLI registration: `_cli.py` + new entries in `pyproject.toml [project.scripts]` (re-check mtime before edit)
- `REDI/README.md` NRD1-style
- `REDI/configs/20250312_NNS_complex_satmut.yaml` config (reproduces the canonical run)
- 1–2 sanity tests under `tests/redi/` (smoke import, plate_db load)
- This HANDOFF.md (updated at session end with completed/deferred status)

**Deferred to next session(s):**
- J.4 re-runs (20250312 NNS bit-identical verification; 20240305_RBYG re-run — blocked by missing keyfiles anyway)
- J.5 v2/v3/v_noncoding script-dedup audit (requires reading historical `process_REDI_bc_counts_table*.py` variants across screens)
- step_j/step_k per-screen rearray drivers — these are **project-specific, not engine** (NRD1/NAB3/SEN1/auxin/plate_2668 collections); they go to `REDI/screens/<screen>/driver.py`, not into `pipeline/`
- 20240305_RBYG re-run (blocked: missing `keyfiles/`)
- 20230616_yL310-yL315, 20220511_PIXL_REDI, 20250627_NNS_validation migrations
- `analyze_agar_phenotyping.R` port decision (OQ-26)

**Bit-identical caveat:** step_c/step_d/step_e logic is dict-iteration + collections.Counter heavy. Even pure-Python with Python 3.7+ dict-insertion-order this is reproducible across reruns of the SAME script, but a refactored engine may emit columns in different order or use sorted vs insertion-order keys. **Documented diff strategy:** verify equivalence by joined-on-keys row count + per-cell value equality on a sorted projection — not byte-identical TSV. Update README + exit criterion 2 to reflect this.

---

## 4. Architecture (concrete file plan)

**Engine: `src/magestic/pipelines/redi/`**
```
__init__.py              (exists)
_cli.py                  Umbrella + per-step CLIs
config.py                YAML config schema (pydantic-light dataclass)
core/
  __init__.py
  config.py              REDIConfig dataclass; load_config(yaml_path)
  purity.py              top-bc1 / top-4 scoring (from step_d)
  correspondence.py      UNIFIED API for type-(i) WGS + type-(ii) re-REDI
  annotate.py            join with bc1_donor_bc0 reference table (from step_e)
  cross_array.py         SHA435 vs yT182/yT183 agreement (from step_f, generalized)
  wgs_check.py           REDI_bc1 vs WGS bc1 (from step_g)
readers/
  __init__.py
  fastq.py               trimmed/merged/collapsed-counts FASTQ reader (consumes step_a output)
  keyfiles.py            combined keyfile loader (from step_b)
  plate_db.py            load_plate_db / resolve_plate / trace_lineage (snapshot consumer)
  wgs_outcome.py         shengdi wgs_outcome_full_table.v0.tsv wrapper
pipeline/
  __init__.py
  step_a_trim_counts.py        (wraps step_a_trim_merge_collapse_counts_tables.sh logic in Python orchestrator)
  step_b_combine_keyfiles.py
  step_c_process_bcs.py        REDI_bc + MAGESTIC_bc extraction & counting
  step_d_purity.py             purity + bc1 edit distance
  step_e_annotate.py           annotate with bc1 reference table
  step_f_cross_array.py        SHA435/yT182/yT183 agreement
  step_g_wgs_check.py
  step_h_plot.py               plot REDI isolated clones (renamed from step_h_plot_REDI_isolated_clones.py)
  step_i_select.py             generic clone selection (parametrized; collection-specific in screens/)
  run_pipeline.py              orchestrator
qc/
  __init__.py
  pipeline_qc.py
scripts/
  __init__.py
slurm/
  template.sbatch              uses sherlock-sbatch conventions
snakemake/
  Snakefile.rules
```

**Project dir: `software/MAGESTIC/REDI/`**
```
README.md                          NRD1-style
configs/
  20250312_NNS_complex_satmut.yaml
  20240305_RBYG_REDI.yaml          (skeleton; keyfiles missing - see §5)
data/
  bc1_to_edit_from_WGS.tsv         (J.3 type-i output; placeholder, run pending)
  REDI_bc_to_MAGESTIC_bc1.tsv      (J.3 type-ii output; placeholder)
  REDI_bc1_correspondence_master.tsv  (joined; placeholder)
plate_database/
  export_plate_db.py               Reads xlsx, writes TSV + parquet
  plates.tsv                       (generated snapshot — not committed initially)
  plates/{Source,Target,REDI_mating,Misc,Rearray_Target,Consolidated_384}.parquet
reference/
  20230513_updated_REDI_bc_V5_key.tsv  (copy from 20250312 keyfiles)
docs/
  pipeline.md                      step_a–step_i prose
  config_schema.md
screens/
  20250312_NNS_complex_satmut/
    driver.py                      per-screen orchestrator + step_j/step_k clone-selection logic
tests/
  test_smoke.py
  test_plate_db.py
snakemake/
  Snakefile
```

**CLI entry points added to `pyproject.toml`:**
```
magestic-redi             = magestic.pipelines.redi._cli:umbrella_main
magestic-redi-trim        = magestic.pipelines.redi._cli:trim_main
magestic-redi-combine-keys = magestic.pipelines.redi._cli:combine_keys_main
magestic-redi-extract-bc  = magestic.pipelines.redi._cli:extract_bc_main
magestic-redi-purity      = magestic.pipelines.redi._cli:purity_main
magestic-redi-annotate    = magestic.pipelines.redi._cli:annotate_main
magestic-redi-cross-array = magestic.pipelines.redi._cli:cross_array_main
magestic-redi-wgs-check   = magestic.pipelines.redi._cli:wgs_check_main
magestic-redi-rearray     = magestic.pipelines.redi._cli:rearray_main
```

`[project.optional-dependencies]` adds `redi = ["openpyxl>=3.1", "pyarrow>=14"]`.

---

## 5. Open Questions / Blockers for next agent

1. **`20240305_RBYG_REDI/keyfiles/` is EMPTY.** Cannot re-run J.4 step 2 (20240305_RBYG) without the REDI_bc_V?_key, sample_key, plate keys. **Required from Kevin or notebook reconstruction.** Likely OQ-22 + OQ-23 root cause: the 20240305 run may have used V5 (`20230513_updated_REDI_bc_V5_key.tsv`) shared with 20250312, but `sample_key.tsv`, `WGS_plate_key.tsv`, etc. are screen-specific. Surface as the first thing the next agent asks Kevin.
2. **OQ-22 confirmed:** scripts dir empty AND keyfiles dir empty → 20240305_RBYG REDI was likely run elsewhere (different repo / informal scripts) or stalled at raw FASTQ. Raw data present though — `raw_data/` is populated.
3. **Bit-identical exit criterion (J exit-2):** likely impossible due to engine refactor changing column ordering / dict iteration; document equivalence-on-sorted-projection in README and revise exit criterion. (Already flagged in §3.)
4. **OQ-23 REDI_bc keyfile versioning:** 20250312 uses V5 (`20230513_updated_REDI_bc_V5_key.tsv`). Need V6/V7 inventory before final `reference/` content lock — see other REDI screen dirs.
5. **OQ-25 cross-array scope:** step_f hardcodes SHA435 vs yT182. Engine generalizes via `cross_array.compare(array_a, array_b)`; configs declare the pairs. Defer decision on which pairs to surface in `data/`.
6. **OQ-26 R port:** `analyze_agar_phenotyping.R` left as a sibling `.R` file (NRD1 precedent: `fitness_slope_heatmap.R`). Decision: minimal-disruption — keep as `.R`, document as a sibling, don't port.
7. **OQ-27 plate DB write-back:** Snapshot reads xlsx only. Writes to TSV/parquet. **No write-back to xlsx in this engine.** Sets the answer.
8. **OQ-29:** Snapshot, not live xlsx — answered same.

---

## 6. Resume instructions for the next agent

1. Read this HANDOFF in full + `__init__.py` docstring.
2. Verify `pyproject.toml` mtime — if newer than this session's last commit timestamp, check `git status` for WGS B.2 agent's `pool_edit_characterization` entries and rebase your edits around them.
3. Read each remaining step file you need to port: `cat /oak/.../REDI/20250312_NNS_complex_saturation_mutagenesis/scripts/step_<x>_*.{py,sh}`.
4. Continue from §3 "Deferred" list. Highest priority deferred: ask Kevin for 20240305 keyfiles (§5 item 1), then run J.4 step 1 (20250312 NNS re-run for equivalence check).
5. Use `git add` with **explicit paths** — top-level repo is in an unusual untracked state.
6. Use `sherlock-sbatch` skill for any SLURM job.
7. Commit every ~10 file writes.

---

## 7. Files touched / created this session (will be appended at end)

(populated at session end)

