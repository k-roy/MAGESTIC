"""magestic.pipelines.wgs.snakemake — Snakemake rules for per-sample fan-out over WGS BAMs.

The Snakefile lives as a flat file alongside this module so the workflow can be
referenced without importing `magestic` inside the snakemake environment:

    snakemake --snakefile $(python -c "import magestic.pipelines.wgs.snakemake; print(magestic.pipelines.wgs.snakemake.__path__[0])")/Snakefile \
              --config config=/path/to/<cohort>_wgs.yaml

`workflow/Snakefile.wgs` is the repo-level entry point; it resolves this same
Snakefile relative to itself and adds nothing.
"""
