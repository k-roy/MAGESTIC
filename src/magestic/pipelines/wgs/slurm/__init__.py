"""magestic.pipelines.wgs.slurm — SLURM array templates for the per-sample fan-out.

`template.sbatch` is the §6 idempotent array task: one sample per array index, read from
a manifest, with a skip-check on its own output at task start so a preempted and requeued
task costs one sample rather than the whole job.
"""
