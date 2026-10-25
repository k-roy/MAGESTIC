"""The WGS cache sees a same-size edit made inside one second (planning/268 §4d; 284 O-W03).

Oak stores whole-second mtime and ctime (arm F measured a sub-second part of 0 and one
distinct mtime over five writes 50 ms apart). `_content_hash` is memoized on the stat
tuple, so a same-size rewrite of a small input inside that second, checked in the same
process, returned the OLD hash and the cache was reused (2/2 in arm F's round 1). The
filesystem here keeps nanoseconds, so the test reproduces Oak by reporting whole-second
times for the one file under test.
"""
import time
from pathlib import Path
from types import SimpleNamespace

from magestic.pipelines.wgs.pipeline import _common as cm


def _whole_second_stat(monkeypatch, target: Path, second_ns: int):
    real = Path.stat

    def fake(self, *args, **kwargs):
        st = real(self, *args, **kwargs)
        if Path(self) == target:
            return SimpleNamespace(st_mode=st.st_mode, st_dev=st.st_dev, st_ino=st.st_ino,
                                   st_size=st.st_size, st_mtime_ns=second_ns,
                                   st_ctime_ns=second_ns)
        return st

    monkeypatch.setattr(Path, "stat", fake)


def test_a_same_size_edit_inside_one_second_changes_the_identity(tmp_path, monkeypatch):
    small = (tmp_path / "keyfile.tsv").resolve()
    small.write_text("sample\toligo\nA1\tO1\n")
    _whole_second_stat(monkeypatch, small, time.time_ns() // 10**9 * 10**9)
    before = cm._file_identity(small)
    small.write_text("sample\toligo\nA1\tO2\n")  # same size, same reported second
    after = cm._file_identity(small)
    assert before["stat"] == after["stat"]  # the Oak condition: stat cannot tell them apart
    assert before["sha256"] != after["sha256"]


def test_an_old_file_still_uses_the_memo(tmp_path, monkeypatch):
    """The memo is kept for files outside the racy window, so resumes do not re-read
    every small input; positive control that the stat key still works there."""
    small = (tmp_path / "old.tsv").resolve()
    small.write_text("x\n")
    old = (time.time_ns() - 10 * cm.RACY_WINDOW_NS) // 10**9 * 10**9
    _whole_second_stat(monkeypatch, small, old)
    calls = []
    real_hash = cm._hash_file
    monkeypatch.setattr(cm, "_hash_file", lambda p: calls.append(p) or real_hash(p))
    cm._content_hash.cache_clear()
    first = cm._file_identity(small)
    second = cm._file_identity(small)
    assert first == second and len(calls) == 1
