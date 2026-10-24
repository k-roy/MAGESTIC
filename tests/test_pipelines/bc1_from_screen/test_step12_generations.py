"""Step 12 parses generations like steps 10 and 11, and counts a timepoint once (218 patch 49).

Ledger item 9 (SCR-D13) and §3 item 5: step 12's sample-key parser turned 0.5 into generation 0 (a
baseline), raised on "6.0", silently skipped "none" and let a later key row overwrite an earlier one;
the time-series check counted replicate rows of one generation as extra timepoints. The estimand
changes in item 9 (zeros, cultures, strata, shared t0) wait on Kevin:
planning/218_magestic_production_corrections/reviews/step12_slope_trajectory/DESIGN.md.
"""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from magestic.pipelines.bc1_from_screen.core.time_series import TimeSeriesConfig, analyze_time_series_for_variant

pytest.importorskip("statsmodels", reason="12_slope_fitness.py imports core/slope_fitness.py, which needs statsmodels")
STEP12 = (Path(__file__).resolve().parents[3] / "src/magestic/pipelines/bc1_from_screen/snakemake/12_slope_fitness.py")


def step12():
    spec = importlib.util.spec_from_file_location("screen_step12", STEP12)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def key(tmp_path, generations, names=None):
    names = names or [f"s{i}" for i in range(len(generations))]
    path = tmp_path / "key.tsv"
    pd.DataFrame({"sample_name": names, "condition": "Caffeine", "library": "yL437",
                  "generations": generations}).to_csv(path, sep="\t", index=False)
    return path


def test_a_fractional_generation_is_refused_not_made_a_baseline(tmp_path):
    with pytest.raises(ValueError, match="s1=0.5"):
        step12().parse_sample_key(key(tmp_path, [0, 0.5, 6, 12]), "Caffeine", "yL437")


def test_a_float_formatted_generation_is_read(tmp_path):
    mapping = step12().parse_sample_key(key(tmp_path, ["0", "6.0", "12"]), "Caffeine", "yL437")
    assert mapping == {"s0": 0, "s1": 6, "s2": 12}


def test_a_missing_generation_is_refused_not_skipped(tmp_path):
    with pytest.raises(ValueError, match="s1="):
        step12().parse_sample_key(key(tmp_path, ["0", "none", "12"]), "Caffeine", "yL437")


def test_conflicting_key_rows_for_one_sample_are_refused(tmp_path):
    with pytest.raises(ValueError, match="s1"):
        step12().parse_sample_key(key(tmp_path, [0, 6, 12], names=["s0", "s1", "s1"]), "Caffeine", "yL437")


def test_identical_repeated_key_rows_are_accepted(tmp_path):
    mapping = step12().parse_sample_key(key(tmp_path, [0, 6, 6], names=["s0", "s1", "s1"]), "Caffeine", "yL437")
    assert mapping == {"s0": 0, "s1": 6}


def variant(generations):
    return pd.DataFrame({"comparison": [f"20251014_Caffeine_{g}" for g in generations],
                         "mean_log2fc": np.linspace(-0.5, -2.0, len(generations))})


def test_replicate_rows_of_one_generation_are_not_extra_timepoints():
    """Two generations with two replicates each: four rows, two timepoints, below a minimum of three."""
    result = analyze_time_series_for_variant(variant([6, 6, 12, 12]), TimeSeriesConfig(min_timepoints=3))
    assert (result["n_timepoints"], result["n_observations"]) == (2, 4)
    assert np.isnan(result["fitness_slope"]) and not result["time_series_validated"]


def test_three_generations_with_replicates_are_fitted():
    result = analyze_time_series_for_variant(variant([6, 6, 12, 12, 21, 21]), TimeSeriesConfig(min_timepoints=3))
    assert (result["n_timepoints"], result["n_observations"]) == (3, 6)
    assert result["fitness_slope"] < 0
