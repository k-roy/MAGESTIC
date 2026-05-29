"""magestic.pipelines.redi.snakemake — Snakemake rules for per-plate fan-out.

Snakefile and rules live as flat files alongside this module so the workflow
can be referenced via `snakemake --snakefile $(python -c "import magestic.pipelines.redi.snakemake; print(magestic.pipelines.redi.snakemake.__path__[0])")/Snakefile`.
"""
