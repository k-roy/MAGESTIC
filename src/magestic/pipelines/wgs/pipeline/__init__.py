"""magestic.pipelines.wgs.pipeline — step orchestrators (step_1..step_5, DIRECT path A).

The step modules are imported lazily by `run_pipeline.run()` rather than here,
because each pulls in TRACE and pysam — the `magestic[wgs]` extra — and
`magestic-wgs --help` must not require them to be installed. (`--dry-run` does
import them: it resolves each step's output path from the step module that owns
it, so "would this run?" answers honestly when the extra is missing.)
"""

from .run_pipeline import STEPS, main, run

__all__ = ["STEPS", "run", "main"]
