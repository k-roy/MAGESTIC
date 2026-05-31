"""
magestic.pipelines.bc1_donor_bc0._cli — CLI entry point.

Entry point registered in pyproject.toml:
    magestic-bc1-link -> main()

Author: Kevin R. Roy
"""

from .pipeline.run_pipeline import main

__all__ = ["main"]
