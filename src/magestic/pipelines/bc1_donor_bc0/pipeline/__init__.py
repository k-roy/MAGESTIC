"""
Pipeline step modules for the BC1-Donor-BC0 pipeline.

The main entry point is run_pipeline.py which orchestrates:
1. Data loading from available sources
2. Early merge of data sources
3. Unified sliding window oligo matching
4. Output of bc1 reference table
"""

from .run_pipeline import run_pipeline, aggregate_and_filter

__all__ = [
    'run_pipeline',
    'aggregate_and_filter',
]
