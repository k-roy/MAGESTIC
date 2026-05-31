"""
Thread Limiting Utility for SLURM Job Arrays

This module MUST be imported BEFORE numpy, scipy, sklearn, pydeseq2, or any other
numerical package that uses threading or joblib.

Usage:
    # At the VERY TOP of your script, before other imports:
    from magestic.utils.thread_limiting import limit_threads
    limit_threads()

    # Now safe to import numerical packages
    import numpy as np
    import pandas as pd
    from pydeseq2.dds import DeseqDataSet

The limit_threads() function:
1. Reads SLURM_CPUS_PER_TASK from environment (defaults to 1 if not set)
2. Sets all thread-limiting environment variables
3. Configures joblib.parallel_config if joblib is available

This prevents the resource usage mismatch warnings from Sherlock when running
job arrays, where each task spawns more threads than its allocated cores.

Key insight: joblib uses the 'loky' backend which ignores JOBLIB_WORKERS.
You MUST set LOKY_MAX_CPU_COUNT to limit loky worker processes.
"""

import os
from typing import Optional


def limit_threads(n_threads: Optional[int] = None) -> int:
    """
    Set all thread-limiting environment variables for numerical packages.

    MUST be called BEFORE importing numpy, scipy, sklearn, pydeseq2, etc.

    Args:
        n_threads: Number of threads to allow. If None, reads from
                   SLURM_CPUS_PER_TASK environment variable (defaults to 1).

    Returns:
        The number of threads that was set.

    Example:
        # At the very top of your script:
        from magestic.utils.thread_limiting import limit_threads
        limit_threads()

        # Now import numerical packages
        import numpy as np
        from pydeseq2.dds import DeseqDataSet
    """
    if n_threads is None:
        n_threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))

    n_str = str(n_threads)

    # OpenMP threads (used by numpy/scipy BLAS implementations)
    os.environ["OMP_NUM_THREADS"] = n_str

    # OpenBLAS threads
    os.environ["OPENBLAS_NUM_THREADS"] = n_str

    # Intel MKL threads
    os.environ["MKL_NUM_THREADS"] = n_str

    # macOS Accelerate framework
    os.environ["VECLIB_MAXIMUM_THREADS"] = n_str

    # NumExpr threads
    os.environ["NUMEXPR_NUM_THREADS"] = n_str

    # CRITICAL: loky (joblib's backend) uses LOKY_MAX_CPU_COUNT, NOT JOBLIB_WORKERS
    # This is the key variable that actually limits joblib worker processes
    os.environ["LOKY_MAX_CPU_COUNT"] = n_str

    # Also configure joblib if it's available
    try:
        import joblib
        joblib.parallel_config(n_jobs=n_threads)
    except ImportError:
        pass
    except Exception:
        # joblib.parallel_config might not be available in older versions
        pass

    return n_threads


def get_slurm_cpus() -> int:
    """
    Get the number of CPUs allocated to this SLURM task.

    Returns 1 if not running under SLURM.
    """
    return int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))


def is_slurm_array_task() -> bool:
    """Check if we're running as part of a SLURM job array."""
    return "SLURM_ARRAY_TASK_ID" in os.environ


# Auto-configure when this module is imported
# This provides a safety net, but scripts should still call limit_threads()
# explicitly for clarity
_auto_configured = False

def _auto_configure():
    """
    Automatically configure thread limits when this module is imported.

    This is a safety net - scripts should still call limit_threads() explicitly.
    """
    global _auto_configured
    if not _auto_configured:
        # Only auto-configure if we're in a SLURM environment
        if "SLURM_JOB_ID" in os.environ:
            limit_threads()
        _auto_configured = True


# Run auto-configuration on import
_auto_configure()
