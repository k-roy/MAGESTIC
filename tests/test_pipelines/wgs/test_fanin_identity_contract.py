"""Fan-in checks each shard's sample identity, not only the row count (audit W04).

At 9d654b8 `collect_shards` checked only the aggregate row count, so a shard carrying another
sample's row, or a duplicated sample, passed as long as the total matched the manifest.
"""
from pathlib import Path

import pytest

from magestic.pipelines.wgs.core.manifest import WGSSample
from magestic.pipelines.wgs.pipeline._common import StepIncompleteError, collect_shards


def _records(tmp_path, names):
    return [WGSSample(sample=n, bam=tmp_path / f"{n}.bam", reference=tmp_path / "ref.fa") for n in names]


def _shards(tmp_path, bodies):
    for name, body in bodies.items():
        (tmp_path / f"{name}.tsv").write_text("sample\tvalue\n" + body)
    return lambda record: tmp_path / f"{record.sample}.tsv"


def test_fanin_returns_rows_in_manifest_order(tmp_path):
    path = _shards(tmp_path, {"s2": "s2\t2\n", "s1": "s1\t1\n"})
    header, rows = collect_shards(_records(tmp_path, ["s1", "s2"]), path, label="t")
    assert header == "sample\tvalue" and rows == ["s1\t1", "s2\t2"]


@pytest.mark.parametrize("bodies", [
    {"s1": "s1\t1\n", "s2": "s1\t2\n"},        # s2's shard carries s1: a duplicate identity
    {"s1": "s1\t1\n", "s2": "other\t2\n"},     # a foreign sample in the right place
    {"s1": "s1\t1\n", "s2": "s2\t2\ns2\t3\n"}, # two rows where one is owed
], ids=["duplicate_identity", "foreign_identity", "two_rows"])
def test_fanin_refuses_a_shard_whose_identity_or_shape_is_wrong(tmp_path, bodies):
    path = _shards(tmp_path, bodies)
    with pytest.raises(StepIncompleteError):
        collect_shards(_records(tmp_path, ["s1", "s2"]), path, label="t")
