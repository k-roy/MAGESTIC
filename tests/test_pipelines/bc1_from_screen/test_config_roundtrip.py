"""Configuration is reproducible only when every public setting survives saving."""
from dataclasses import asdict
import json
from pathlib import Path

import pytest
import yaml

from magestic.pipelines.bc1_from_screen.config import (
    BC1FromScreenConfig, HitCallingConfig, LibraryConfig, TierConfig,
)


@pytest.mark.parametrize("configured", [False, True])
def test_all_public_settings_survive_yaml(tmp_path, configured):
    config = BC1FromScreenConfig("project", "screen", tmp_path)
    if configured:
        config.libraries = [LibraryConfig("lib", ["20210102"], True, [0, 4, 8], 3, "test")]
        config.tier_config = TierConfig(111, 23, 2)
        config.hit_calling_config = HitCallingConfig(
            3.1, 0.76, 3.8, 0.86, 4, False, 17, False, 6, 0.2, 0.95,
        )
        config.sequencing_dates = ["20210103"]
        config.base_dir = tmp_path / "base"
        config.min_total_reads_per_sample = 123
        config.min_replicates_per_condition = 3
        config.use_scratch = True
        config.scratch_dir = tmp_path / "scratch"
    output = tmp_path / "config.yaml"
    config.to_yaml(output)
    restored = BC1FromScreenConfig.from_yaml(output)
    assert asdict(restored) == asdict(config)
    assert restored.get_working_dir("counts") == config.get_working_dir("counts")
    assert restored.library_mapping == ({"20210102": "lib"} if configured else {})


def test_legacy_yaml_defaults_and_optional_scratch(tmp_path):
    output = tmp_path / "legacy.yaml"
    data = dict(project_name="project", screen_name="screen", project_dir=str(tmp_path))
    output.write_text(yaml.safe_dump(data))
    restored = BC1FromScreenConfig.from_yaml(output)
    assert not restored.use_scratch and restored.scratch_dir is None
    assert restored.hit_calling_config == HitCallingConfig()
    data.update(use_scratch=True, scratch_dir=str(tmp_path / "scratch"))
    output.write_text(yaml.safe_dump(data))
    restored = BC1FromScreenConfig.from_yaml(output)
    assert restored.get_working_dir("counts") == tmp_path / "scratch" / "counts"


def test_config_writer_records_effective_configuration(tmp_path):
    output = tmp_path / "config.yaml"
    config = BC1FromScreenConfig("project", "screen", tmp_path)
    config.to_yaml(output)
    sidecar = json.loads(Path(str(output) + ".provenance.json").read_text())
    assert sidecar["parameters"]["configuration"] == yaml.safe_load(output.read_text())
    assert sidecar["counts"]["libraries"] == 0


def test_config_writer_does_not_suppress_sidecar_failure(tmp_path, monkeypatch):
    from magestic.utils.provenance import StepProvenance

    def fail_write(self):
        raise OSError("test sidecar write failure")

    monkeypatch.setattr(StepProvenance, "write", fail_write)
    with pytest.raises(OSError, match="test sidecar write failure"):
        BC1FromScreenConfig("project", "screen", tmp_path).to_yaml(tmp_path / "config.yaml")
