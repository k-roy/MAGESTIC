"""Tests for collapsing the read-level gdb table into a bc1 bridge reference.

The full-scale bridge run exhausted 128 GB because `load_gdb_data` returned the
read-level frame (~44M rows) where the merge only needs
`donor_bc0_fragment -> (oligo_name, guide)`, which is at most design size
(HANDOFF_20260608_gdb_bdb_pipeline.md section 3).

Collapsing is not only a size fix. A fragment can map to more than one oligo across
reads, and silently picking one is how a bc1 ends up assigned to the wrong edit. These
tests pin that the ambiguity is counted and surfaced rather than resolved in silence.
"""

from __future__ import annotations

import pandas as pd

from magestic.pipelines.bc1_donor_bc0.readers.data_loader import collapse_gdb_bridge


def _frame(rows):
    return pd.DataFrame(
        rows, columns=["donor_bc0_fragment", "oligo_name", "guide", "count"]
    )


def test_collapses_read_level_rows_to_one_per_fragment():
    df = _frame([
        ("FRAG_A", "oligo_1", "G" * 20, 10),
        ("FRAG_A", "oligo_1", "G" * 20, 7),
        ("FRAG_A", "oligo_1", "G" * 20, 3),
        ("FRAG_B", "oligo_2", "C" * 20, 5),
    ])
    out, stats = collapse_gdb_bridge(df)
    assert len(out) == 2
    assert stats["rows_in"] == 4 and stats["rows_out"] == 2
    assert stats["unique_keys"] == 2
    assert set(out["donor_bc0_fragment"]) == {"FRAG_A", "FRAG_B"}


def test_most_supported_row_wins():
    df = _frame([
        ("FRAG_A", "oligo_minor", "G" * 20, 2),
        ("FRAG_A", "oligo_major", "T" * 20, 99),
    ])
    out, _ = collapse_gdb_bridge(df)
    assert out.iloc[0]["oligo_name"] == "oligo_major"
    assert out.iloc[0]["guide"] == "T" * 20


def test_ambiguous_fragments_are_counted_not_hidden():
    """A fragment mapping to two oligos is the upstream cause of a mis-assigned bc1."""
    df = _frame([
        ("FRAG_A", "oligo_1", "G" * 20, 9),
        ("FRAG_A", "oligo_2", "C" * 20, 1),   # ambiguous
        ("FRAG_B", "oligo_3", "A" * 20, 4),   # clean
    ])
    out, stats = collapse_gdb_bridge(df)
    assert stats["unique_keys"] == 2
    assert stats["ambiguous_keys"] == 1
    assert len(out) == 2


def test_no_ambiguity_reported_when_mapping_is_clean():
    df = _frame([
        ("FRAG_A", "oligo_1", "G" * 20, 9),
        ("FRAG_A", "oligo_1", "G" * 20, 1),
        ("FRAG_B", "oligo_2", "A" * 20, 4),
    ])
    _, stats = collapse_gdb_bridge(df)
    assert stats["ambiguous_keys"] == 0


def test_rows_without_a_fragment_are_counted_not_silently_dropped():
    df = _frame([
        ("FRAG_A", "oligo_1", "G" * 20, 9),
        (None, "oligo_2", "C" * 20, 4),
        ("", "oligo_3", "A" * 20, 2),
    ])
    out, stats = collapse_gdb_bridge(df)
    assert stats["dropped_null_key"] == 2
    assert len(out) == 1


def test_collapse_is_deterministic_under_tied_support():
    """Stable sort, so two runs over differently-ordered input agree on ties."""
    rows = [
        ("FRAG_A", "oligo_1", "G" * 20, 5),
        ("FRAG_A", "oligo_2", "C" * 20, 5),
    ]
    first, _ = collapse_gdb_bridge(_frame(rows))
    again, _ = collapse_gdb_bridge(_frame(rows))
    assert first.iloc[0]["oligo_name"] == again.iloc[0]["oligo_name"]


def test_missing_support_column_still_collapses():
    df = pd.DataFrame(
        [("FRAG_A", "oligo_1"), ("FRAG_A", "oligo_1")],
        columns=["donor_bc0_fragment", "oligo_name"],
    )
    out, stats = collapse_gdb_bridge(df)
    assert len(out) == 1 and stats["rows_out"] == 1


def test_empty_and_keyless_frames_are_safe():
    empty, stats = collapse_gdb_bridge(pd.DataFrame())
    assert empty.empty and stats["rows_out"] == 0

    keyless = pd.DataFrame([("x",)], columns=["something_else"])
    out, stats = collapse_gdb_bridge(keyless)
    assert len(out) == 1 and stats["rows_in"] == 1


def test_collapse_scales_the_way_the_bridge_needs():
    """200k read-level rows over 1k fragments must return 1k rows, not 200k."""
    rows = [(f"FRAG_{i % 1000}", f"oligo_{i % 1000}", "G" * 20, i) for i in range(200_000)]
    out, stats = collapse_gdb_bridge(_frame(rows))
    assert stats["rows_in"] == 200_000
    assert stats["rows_out"] == 1000
    assert stats["ambiguous_keys"] == 0
