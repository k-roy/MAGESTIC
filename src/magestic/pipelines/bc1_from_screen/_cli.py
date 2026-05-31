"""
magestic.pipelines.bc1_from_screen._cli — CLI entry point.

Entry point registered in pyproject.toml:
    magestic-screen-analysis -> main()

Author: Kevin R. Roy
"""

from .pipeline.run_pipeline import main

__all__ = ["main"]
