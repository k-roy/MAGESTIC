"""Regression tests for the BDB order and hash-seed dependences arm G found at master (planning/288 §2 and
the STEP08 ruling), fixed in planning/378."""
import importlib.util
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("rapidfuzz")

SNAKE = Path(__file__).resolve().parents[3] / "src/magestic/pipelines/bc1_donor_bc0/snakemake"


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, SNAKE / fname)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _collapse(frame):
    m = _load("bdb08", "08_integrate_multiple_runs.py")
    out = m.collapse_sequences(frame, ["bc1"], "donor", ["counts"], 1, show_progress=False, n_jobs=1,
                               preserve_cols=["oligo_name"])
    return out.sort_values(["bc1", "donor"]).reset_index(drop=True)


def test_step08_peel_off_does_not_depend_on_row_order():
    # X and Y tie at 5 reads, 1 edit apart; Z (1 read) is 1 edit from Y, 2 from X. Top = X gives {X+Y},{Z};
    # top = Y gives {Y+X+Z}. The tie must be broken by the rows' content, not by which arrived first.
    rows = [{"bc1": "B1", "donor": "ACGTACGTAC", "counts": 5, "oligo_name": "oX"},
            {"bc1": "B1", "donor": "ACGTACGTAA", "counts": 5, "oligo_name": "oY"},
            {"bc1": "B1", "donor": "ACGTACGTTA", "counts": 1, "oligo_name": "oZ"}]
    a, b = _collapse(pd.DataFrame(rows)), _collapse(pd.DataFrame(rows[::-1]))
    assert a.to_dict("records") == b.to_dict("records")
    assert a.to_dict("records") == [{"bc1": "B1", "donor": "ACGTACGTAA", "counts": 11, "oligo_name": "oY"}]


def test_step08_untied_groups_are_unchanged():
    rows = [{"bc1": "B1", "donor": "ACGTACGTAC", "counts": 9, "oligo_name": "o1"},
            {"bc1": "B1", "donor": "ACGTACGTAA", "counts": 2, "oligo_name": "o2"}]
    out = _collapse(pd.DataFrame(rows[::-1])).to_dict("records")
    assert out == [{"bc1": "B1", "donor": "ACGTACGTAC", "counts": 11, "oligo_name": "o1"}]


def test_step06_closest_design_does_not_depend_on_the_hash_seed(tmp_path):
    src = Path(__file__).resolve().parents[3] / "src"
    probe = tmp_path / "probe.py"
    probe.write_text(textwrap.dedent(f"""
        import importlib.util, sys, random
        sys.path.insert(0, {str(src)!r})
        spec = importlib.util.spec_from_file_location("step06", {str(SNAKE / '06_annotate_unmatched_donors.py')!r})
        m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
        rnd = random.Random(5)
        donor = "".join(rnd.choice("ACGT") for _ in range(m.DONOR_LENGTH))
        lookup = {{}}
        for i, name in enumerate(["CDS_DEL_guide_D", "CDS_DEL_guide_B", "CDS_DEL_guide_A", "CDS_DEL_guide_C"]):
            lookup[name] = m.OligoInfo(oligo_name=name, guide=f"G{{i}}" * 5, donor=donor,
                                       middle_donor=donor[m.MIDDLE_DONOR_START:m.MIDDLE_DONOR_END], oligo_seq=donor)
        observed = donor[:60] + ("A" if donor[60] != "A" else "C") + donor[61:]
        print(m.find_closest_designed_donor(observed, lookup, m.build_kmer_index(lookup))[0])
    """))
    winners = set()
    for seed in ("0", "1", "2", "3", "4", "5", "6", "7"):
        out = subprocess.run([sys.executable, str(probe)], env=dict(os.environ, PYTHONHASHSEED=seed),
                             capture_output=True, text=True)
        assert out.returncode == 0, out.stderr[-1500:]
        winners.add(out.stdout.strip().splitlines()[-1])
    assert winners == {"CDS_DEL_guide_A"}
