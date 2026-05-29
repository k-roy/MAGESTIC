# REDI — barcode-isolated clone array engine for MAGESTIC 3.0

Code + processed data for **REDI (Restriction Enzyme Diagnostic Indexing)**, the arrayed-clone validation platform for MAGESTIC 3.0 yeast saturation editing screens. Operates on paired-end Illumina reads from REDI mating plates and yields:

1.  Per-(REDI_bc, MAGESTIC_bc1) read-count tables
2.  Colony purity scores (top-bc1 / top-4-sum within each clone)
3.  Cross-array agreement metrics (e.g., SHA345 vs yT182/yT183 reference arrays)
4.  Annotation against the bc1-donor-bc0 reference table (output of `magestic.pipelines.bc1_donor_bc0`)
5.  **Two complementary correspondence types** (per `MASTER_PLAN.md §9d`):
    - **Type (i)** `bc1 ↔ edit` from WGS (low-throughput, ground truth)
    - **Type (ii)** `REDI-bc ↔ MAGESTIC-bc1` from re-REDI (high-throughput, same array)

Mirrors the NRD1 precedent: project-level dir at the package root holds configs + screen drivers; reusable engine lives in `src/magestic/pipelines/redi/`.

---

## Repository structure

```
REDI/
├── README.md                  (this file)
├── HANDOFF_phase_J_20260528.md  Migration session state
│
├── configs/                   per-screen YAML configs
│   ├── 20250312_NNS_complex_satmut.yaml
│   └── 20240305_RBYG_REDI.yaml      (skeleton; keyfiles pending — see HANDOFF §5)
│
├── data/                      published output artifacts
│   ├── bc1_to_edit_from_WGS.tsv         Type-(i) correspondence
│   ├── REDI_bc_to_MAGESTIC_bc1.tsv      Type-(ii) correspondence
│   └── REDI_bc1_correspondence_master.tsv   Joined; published
│
├── plate_database/            agar plate barcode DB snapshot
│   ├── export_plate_db.py     reads lab_shared xlsx, writes parquet + TSV
│   ├── plates.tsv             flat snapshot (regenerated)
│   └── plates/{Source,Target,REDI_mating,Misc,Rearray_Target,Consolidated_384}.parquet
│
├── reference/                 versioned REDI_bc keyfiles
│   └── 20230513_updated_REDI_bc_V5_key.tsv
│
├── docs/
│   ├── pipeline.md            step_a–step_i prose with primer table
│   └── config_schema.md       YAML schema reference
│
├── screens/                   per-screen drivers (step_j/step_k logic)
│   └── 20250312_NNS_complex_satmut/
│       └── driver.py
│
├── snakemake/                 (reserved — workflow lives in src/.../snakemake/)
└── tests/
    ├── test_smoke.py
    └── test_plate_db.py
```

The engine itself is `src/magestic/pipelines/redi/`:

```
src/magestic/pipelines/redi/
├── __init__.py
├── _cli.py                CLI entry points (one per step + umbrella)
├── core/
│   ├── config.py          REDIConfig + load_config(yaml)
│   ├── purity.py          step_d
│   ├── annotate.py        step_e
│   ├── cross_array.py     step_f (generalized beyond SHA345/yT182)
│   ├── wgs_check.py       step_g
│   └── correspondence.py  Unified two-type API (J.3)
├── readers/
│   ├── fastq.py           bc1/REDI_bc extraction (step_c primitives)
│   ├── keyfiles.py        step_b combine_keyfiles + V5/yT182 consolidation
│   ├── plate_db.py        snapshot consumer (load / resolve / lineage)
│   └── wgs_outcome.py     Shengdi WGS table wrapper
├── pipeline/
│   ├── step_a_trim_counts.py    (per-sample shell orchestrator)
│   ├── step_b_combine_keys.py
│   ├── step_c_process_bcs.py
│   ├── step_d_purity.py
│   ├── step_e_annotate.py
│   ├── step_f_cross_array.py
│   ├── step_g_wgs_check.py
│   ├── step_h_plot.py
│   ├── step_i_select.py
│   └── run_pipeline.py     end-to-end orchestrator
├── qc/
│   └── pipeline_qc.py
├── snakemake/Snakefile     per-plate fan-out
└── slurm/template.sbatch   per-sample step_a SLURM template
```

---

## Pipeline (step_a → step_i)

| Step | CLI                          | Input                               | Output                                                  |
|------|------------------------------|-------------------------------------|---------------------------------------------------------|
| a    | `magestic-redi-trim`         | raw fastq.gz pair                   | `<sample>_R1_R2_counts.tsv` (collapsed counts)          |
| b    | `magestic-redi-combine-keys` | sample_key + gDNA + REDI_plate_bc   | `combined_key.tsv`                                      |
| c    | `magestic-redi-extract-bc`   | counts tables + combined_key        | `combined_REDI_bc_MAGESTIC_bc_merged_df.tsv` + `MAGESTIC_bc_only_grouped_df.tsv` |
| d    | `magestic-redi-purity`       | step_c output                        | `top_REDI_bc_MAGESTIC_bc1_df.tsv`                        |
| e    | `magestic-redi-annotate`     | step_d + bc1_donor_bc0 ref           | `annotated_top_REDI_bc_MAGESTIC_bc1_df.tsv`             |
| f    | `magestic-redi-cross-array`  | step_e + 1536/384 plate format       | `<A>_<B>_merged_REDI_df.tsv`                            |
| g    | `magestic-redi-wgs-check`    | step_f + WGS outcome + WGS plate key | `WGS_bc1_call_outer_merge_with_REDI_df_and_plate_key.tsv` |
| h    | `magestic-redi-rearray`      | step_e/f                             | `plots/bc1_clones_distribution.png` + per-screen rearray |
| i    | (per-screen driver)          | step_h filtered df                   | per-screen clone-selection TSVs                          |

The umbrella `magestic-redi --config <yaml>` runs step_b through step_h; step_a is per-sample and parallelized via the Snakemake workflow or SLURM array.

### Primer table (default `primer_rules`)

| Inner primer set | bc1 on read | REDI_bc on read |
|------------------|-------------|------------------|
| KR1965           | R2          | R1                |
| KR1967           | R1          | R2                |

Override per-screen in YAML if the screen used a different primer set.

---

## Data outputs (`data/`)

| Artifact                              | Source                       | Notes                                            |
|---------------------------------------|------------------------------|--------------------------------------------------|
| `bc1_to_edit_from_WGS.tsv`            | step_g + WGS outcome table   | Type-(i) ground truth                            |
| `REDI_bc_to_MAGESTIC_bc1.tsv`         | step_f                       | Type-(ii) high-throughput                        |
| `REDI_bc1_correspondence_master.tsv`  | Join of the above             | `concordance` flag: `WGS_only` / `reREDI_only` / `both_agree` / `both_disagree` |

---

## Reproducibility & equivalence

Exit criterion 2 in `MASTER_PLAN.md §9d` calls for bit-identical reproduction of the legacy `combined_data_tables/` outputs. Due to engine refactors that change column ordering and dict-iteration order, **bit-identical equality is not the success criterion** — equivalence is verified on a sorted projection:

```python
left.sort_values(KEY_COLS).reset_index(drop=True).equals(
    right.sort_values(KEY_COLS).reset_index(drop=True)
)
```

where `KEY_COLS = ['sample_name', 'REDI_bc', 'bc1']` for the merged table.

---

## Status (2026-05-28)

- Engine scaffold + step_b/c/d/e/f/g/h/i modules: complete.
- Plate-DB snapshot tool: complete.
- 20250312_NNS bit-identical / equivalence run: **pending** (see HANDOFF §3).
- 20240305_RBYG re-run: **blocked** on missing `keyfiles/` (see HANDOFF §5 item 1).
- Snakemake workflow: skeleton only.
- Per-screen step_j/step_k drivers: not migrated (project-specific; live under `screens/`).

See `HANDOFF_phase_J_20260528.md` for the full session log and next-agent instructions.
