"""Step 11 names its DESeq2 size-factor estimator (planning/272).

pydeseq2's default "ratio" silently becomes the iterative estimator when every feature has a zero in some
sample. A screen matrix with a consolidated t0 always does (yL406: 0 of 64,695 bc1s zero-free), and the
iterative fit ran for hours on one comparison. The estimator is now an explicit option, and the default is
unchanged; the run log says which estimator was actually used.
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("pydeseq2")
SCRIPT = (Path(__file__).resolve().parents[3] / "src/magestic/pipelines/bc1_from_screen/snakemake"
          / "11_deseq2_single.py")


def _tables(tmp_path):
    rng = np.random.default_rng(0)
    n = 150
    bc1 = [f"BC{i:04d}" for i in range(n)]
    counts = rng.poisson(30, size=(n, 6)) + 1
    counts[np.arange(n), np.arange(n) % 6] = 0          # one zero per bc1, element-wise (not iloc's outer grid)
    ref = pd.DataFrame(counts, index=bc1, columns=[f"r{i}" for i in range(6)])
    treat = pd.DataFrame(rng.poisson(30, size=(n, 3)) + 1, index=bc1, columns=[f"t{i}" for i in range(3)])
    ref.index.name = treat.index.name = "bc1"
    paths = tmp_path / "ref.tsv", tmp_path / "treat.tsv"
    ref.to_csv(paths[0], sep="\t"); treat.to_csv(paths[1], sep="\t")
    return paths


@pytest.mark.parametrize("fit_type", ["ratio", "poscounts"])
def test_step11_runs_and_reports_its_size_factor_estimator(tmp_path, fit_type):
    ref, treat = _tables(tmp_path)
    out = tmp_path / f"{fit_type}.tsv"
    r = subprocess.run([sys.executable, str(SCRIPT), "--reference", ref, "--treatment", treat, "--output", out,
                        "--comparison", "c1", "--library", "yLX", "--size-factors-fit-type", fit_type],
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr
    log = r.stdout + r.stderr
    if fit_type == "ratio":
        assert "0 of 150 BC1s are zero-free" in log and "ITERATIVE" in log
    else:
        assert "size factors: poscounts" in log
    res = pd.read_csv(out, sep="\t")
    assert len(res) == 150 and res["log2FoldChange"].notna().all()


def test_the_default_estimator_is_unchanged(tmp_path):
    """Without the flag, the log reports pydeseq2's default (ratio), so no existing run changes."""
    ref, treat = _tables(tmp_path)
    r = subprocess.run([sys.executable, str(SCRIPT), "--reference", ref, "--treatment", treat,
                        "--output", tmp_path / "d.tsv", "--comparison", "c1", "--library", "yLX"],
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "ratio requested" in r.stdout + r.stderr
