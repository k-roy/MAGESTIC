"""Smoke tests — imports + config round-trip + tiny in-memory pipeline pieces."""

from __future__ import annotations

import textwrap
from pathlib import Path

import pandas as pd
import pytest


def test_import_engine():
    import magestic.pipelines.redi as r

    assert r.__name__ == "magestic.pipelines.redi"
    from magestic.pipelines.redi import core, readers, pipeline, qc  # noqa: F401


def test_config_round_trip(tmp_path: Path):
    from magestic.pipelines.redi.core.config import load_config

    yaml_text = textwrap.dedent(
        f"""
        screen_name: test_screen
        sequencing_date: "20250101"
        paths:
          project_root: {tmp_path}/proj/
          raw_fastq:    {tmp_path}/proj/fastq/
          intermediate: {tmp_path}/proj/intermediate/
          processed:    {tmp_path}/proj/processed/
          keyfiles:     {tmp_path}/proj/keyfiles/
        barcode:
          bc1_length: 20
          fwd_bc1_prefix: CAGGTCATGCTC
          fwd_bc1_suffix: CCCTAGGATACG
          fallback_prefix_length: 16
        purity:
          top_n_for_purity: 4
          edit_distance_drop: [1, 2]
          min_counts: 10
        array_pairs:
          - name_a: A
            name_b: B
            purpose_a: "A REDI"
        """
    )
    yaml_path = tmp_path / "test.yaml"
    yaml_path.write_text(yaml_text)
    cfg = load_config(yaml_path)
    assert cfg.screen_name == "test_screen"
    assert cfg.sequencing_date == "20250101"
    assert cfg.barcode.bc1_length == 20
    assert cfg.purity.top_n_for_purity == 4
    assert cfg.array_pairs[0].name_a == "A"


def test_extract_bc1_from_seq():
    from magestic.pipelines.redi.readers.fastq import extract_bc1_from_seq, _rc

    # 20 bp prefix-anchored bc1 (default args, no suffix)
    bc1 = "A" * 20
    seq = "GGGG" + "CAGGTCATGCTC" + bc1 + "CCCTAGGATACG" + "TTTT"
    assert extract_bc1_from_seq(seq) == bc1
    # With suffix passed, exact suffix match at expected length
    assert extract_bc1_from_seq(seq, fwd_bc1_suffix="CCCTAGGATACG") == bc1
    # 19 bp bc1 with length_tolerance=1
    bc1_19 = "C" * 19
    seq_19 = "GGGG" + "CAGGTCATGCTC" + bc1_19 + "CCCTAGGATACG" + "TTTT"
    assert extract_bc1_from_seq(seq_19, fwd_bc1_suffix="CCCTAGGATACG") == bc1_19
    # 21 bp bc1 with length_tolerance=1
    bc1_21 = "G" * 21
    seq_21 = "GGGG" + "CAGGTCATGCTC" + bc1_21 + "CCCTAGGATACG" + "TTTT"
    assert extract_bc1_from_seq(seq_21, fwd_bc1_suffix="CCCTAGGATACG") == bc1_21
    # Reverse orientation: same construct read from the other strand
    bc1_fwd = "ACGTACGTACGTACGTACGT"
    fwd = "GGGG" + "CAGGTCATGCTC" + bc1_fwd + "CCCTAGGATACG" + "TTTT"
    assert extract_bc1_from_seq(_rc(fwd), fwd_bc1_suffix="CCCTAGGATACG") == bc1_fwd
    # Missing prefix → fallback slice
    short_seq = "X" * 16 + "Y" * 20
    assert extract_bc1_from_seq(short_seq) == "Y" * 20
    # None / NaN passthrough
    assert extract_bc1_from_seq(None) is None


def test_purity_basic():
    """End-to-end purity calc on a 6-row synthetic dataset."""
    from magestic.pipelines.redi.core.purity import compute_purity

    # 1 sample, 1 REDI_bc, 3 distinct bc1's (all edit-distance >= 5 from each
    # other), counts {100, 30, 10, 5}. top_n=4 → total_counts = 145.
    # colony_purity for top = 100/145 ≈ 0.6897
    df = pd.DataFrame(
        {
            "bc1": ["AAAAAAAAAAAAAAAAAAAA", "GGGGGGGGGGGGGGGGGGGG", "TTTTTTTTTTTTTTTTTTTT", "CCCCCCCCCCCCCCCCCCCC"],
            "REDI_bc": ["R1", "R1", "R1", "R1"],
            "sample_name": ["s1", "s1", "s1", "s1"],
            "counts": [100, 30, 10, 5],
            "REDI_plate_color_384": ["pink"] * 4,  # required by purity filter
        }
    )
    out = compute_purity(df)
    assert len(out) == 1
    row = out.iloc[0]
    assert row["top_bc1"] == "AAAAAAAAAAAAAAAAAAAA"
    assert row["total_counts"] == 145
    assert abs(row["colony_purity"] - 100 / 145) < 1e-9
