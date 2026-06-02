# MAGESTIC unified workflow (Phase 1b prototype)

**One packaged workflow, parameterized by a per-screen config — not one Snakefile per screen.**

The per-screen-Snakefile pattern is what produced the fork divergence the
consolidation is undoing: each screen carried its own copy of the rules, so a
single central change (e.g. consolidating the parsers into the auto-detecting
`01_parse_bc1_donor_bc0.py`) silently broke every per-screen copy. Here the
workflow logic lives once; a screen is just a config file.

## Layout

```
software/MAGESTIC/
  workflow/
    Snakefile                      # the single entry point (config-driven)
    config/
      config.example.yaml          # template + worked example (20250912 yL437/yL442)
      screens/                      # one <screen>.yaml per dataset
    profiles/sherlock/             # (to add) cluster profile w/ the mylab+AVX-512+set-u traps
  src/magestic/pipelines/bc1_donor_bc0/snakemake/
    Snakefile.rules                # the canonical rule library (parse_01 .. assign_05)
    01_parse_bc1_donor_bc0.py ...  # the step scripts
```

`workflow/Snakefile` resolves `Snakefile.rules` **relative to itself** (repo-internal),
so it does not depend on `magestic` being importable inside the snakemake conda env.

## Run a screen

```bash
cd software/MAGESTIC/workflow
snakemake --snakefile Snakefile \
          --configfile config/screens/<screen>.yaml \
          --profile profiles/sherlock          # or: -j N for a local/sdev run
```

## Add a screen — no Snakefile edits

```bash
cp config/config.example.yaml config/screens/20240806_yL406.yaml
# edit: project_dir, sequencing_runs, active_runs, libraries, keyfiles, resources
snakemake --snakefile Snakefile --configfile config/screens/20240806_yL406.yaml ...
```

## Re-run every dataset (the consolidation goal)

```bash
for c in config/screens/*.yaml; do
  snakemake --snakefile Snakefile --configfile "$c" --profile profiles/sherlock
done
```

## What this prototype changes vs the legacy per-screen Snakefiles

- **Rules: one canonical copy** (`Snakefile.rules`) instead of N drifting copies.
- **Parsing: raw R1/R2, auto-detected.** `01_parse_bc1_donor_bc0.py` auto-detects
  read length / orientation / bc1 prefix, so the legacy `step_00` (bbduk+bbmerge+
  collapse) and most of the `sequence_parsing` config block are gone. The config
  is correspondingly leaner.
- **A screen = a config.** Migrating a screen means deleting its
  `workflow/bc1_donor_bc0/Snakefile` and keeping only its config here.

## Prototype caveats (before production cutover)

1. **bbmerge-vs-raw-R1/R2 parity (2x300).** The legacy 2x300 path merged reads
   before parsing; the unified path parses raw R1/R2. Parity-test the unified
   parser output against the legacy merge-then-parse output on a 2x300 sample
   before cutting a screen over (logged in HANDOFF_MAGESTIC_CONSOLIDATION §7).
2. **`profiles/sherlock/`** is not yet committed — add a Snakemake profile that
   encodes the `mylab,owners` account, the AVX-512 `--constraint`, and the
   `set -u`-after-`source ~/.bashrc` trap before running on the cluster.
3. **Other pipelines** (guide_donor_bc0, bc1_from_screen, editing_outcome) get
   the same treatment incrementally — each gets a `Snakefile` here that includes
   its own `src/.../snakemake/Snakefile.rules`.
