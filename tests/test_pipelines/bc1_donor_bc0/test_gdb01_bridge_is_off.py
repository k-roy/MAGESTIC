"""The GDB-01 consumer-side repair in bdb 04 is a BRIDGE and must default OFF (policy §1).

The root fix landed in gdb Step 02 on 2026-09-11 (sequence_anchors.find_re_site). Leaving a downstream
repair enabled by default would silently paper over any future recurrence and would keep only this one
consumer correct while wrong guides stayed on disk for every other reader. It stays reachable so tables
built before the fix remain readable.
"""
import importlib.util
from pathlib import Path

import pandas as pd
import pytest

import magestic
from magestic.pipelines.bc1_donor_bc0.core.matching import (
    build_lookup_tables, repair_gdb01_shifted_guides)

_STEP04 = (Path(magestic.__file__).resolve().parent / "pipelines" / "bc1_donor_bc0" / "snakemake"
           / "04_map_to_oligos.py")


def _load_step04():
    spec = importlib.util.spec_from_file_location("step04_map_to_oligos", _STEP04)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_bridge_is_off_by_default_and_still_reachable():
    p = _load_step04().build_parser()
    base = ["--input", "x.tsv", "--output", "y.tsv"]
    assert p.parse_args(base).gdb01_repair is False, "the GDB-01 bridge must be OFF by default"
    assert p.parse_args(base + ["--gdb01-repair"]).gdb01_repair is True


def test_the_old_flag_is_still_accepted_so_existing_rule_shells_do_not_die():
    """Snakefile.rules and the screen's sbatch scripts may still pass --no-gdb01-repair; argparse would
    SystemExit on an unknown flag and kill the rule."""
    p = _load_step04().build_parser()
    args = p.parse_args(["--input", "x.tsv", "--output", "y.tsv", "--no-gdb01-repair"])
    assert args.no_gdb01_repair is True and args.gdb01_repair is False


def test_the_bridge_restores_nothing_on_corrected_guides():
    """Against tables built after the root fix the bridge is a no-op -- which is what makes it safe to
    leave reachable, and what its 'drop --gdb01-repair' message keys on."""
    oligos = pd.DataFrame({
        "oligo_name": ["LbCas12a_a", "impLbCas12a_b", "SpCas9_c"],
        "guide": ["ACGTACGTACGTACGTACGG", "TTGGTTGGTTGGTTGGTTGG", "CCAACCAACCAACCAACCAA"],
        "donor": ["A" * 129, "C" * 129, "G" * 129],
        "individual_variants": ["10_A_C_SNP", "20_G_T_SNP", "30_C_G_SNP"],
    })
    tables = build_lookup_tables(oligos)
    corrected = pd.Series(oligos["guide"].tolist())
    out, flag = repair_gdb01_shifted_guides(corrected, tables)
    assert not flag.any()
    assert out.tolist() == corrected.tolist()
    # and it still repairs a genuinely shifted pre-fix guide
    shifted = pd.Series(["T" + oligos["guide"][0][:19]])
    out2, flag2 = repair_gdb01_shifted_guides(shifted, tables)
    assert flag2.tolist() == [True] and out2.tolist() == [oligos["guide"][0]]
