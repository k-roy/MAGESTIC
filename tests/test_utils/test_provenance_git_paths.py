"""Git's fixed-width status prefix must retain its significant leading space."""
import subprocess

import pytest

from magestic.utils.provenance import git_state


def git(repo, *args):
    return subprocess.run(["git", *args], cwd=repo, capture_output=True,
                          text=True, check=True).stdout


@pytest.mark.parametrize("state", ["clean", "unstaged", "staged", "two_unstaged"])
def test_actual_git_state_preserves_changed_paths(tmp_path, state):
    git(tmp_path, "init", "-q")
    for name in ("source.py", "second.py"):
        (tmp_path / name).write_text("VALUE = 1\n")
    git(tmp_path, "add", "--", "source.py", "second.py")
    git(tmp_path, "-c", "user.name=Audit fixture", "-c", "user.email=audit@example.invalid",
        "-c", "commit.gpgsign=false", "commit", "-qm", "Synthetic provenance fixture")
    expected_paths = []
    if state != "clean":
        (tmp_path / "source.py").write_text("VALUE = 2\n")
        expected_paths.append("source.py")
    if state == "staged":
        git(tmp_path, "add", "--", "source.py")
    if state == "two_unstaged":
        (tmp_path / "second.py").write_text("VALUE = 2\n")
        expected_paths.append("second.py")
    observed = git_state(tmp_path)
    assert observed["commit"] == git(tmp_path, "rev-parse", "HEAD").strip()
    assert observed["dirty"] == bool(expected_paths)
    assert sorted(observed["dirty_files"]) == sorted(expected_paths)
