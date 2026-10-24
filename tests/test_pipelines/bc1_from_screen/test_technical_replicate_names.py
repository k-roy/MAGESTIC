"""Technical aggregation must never overwrite a distinct sample or metadata."""
import importlib.util
from pathlib import Path

import magestic
import pandas as pd
import pytest


def module():
    source = Path(magestic.__file__).parent / "pipelines/bc1_from_screen/snakemake/09_combine_technical_replicates.py"
    spec = importlib.util.spec_from_file_location("technical_replicate_names_test", source)
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


def fixture(group):
    counts = pd.DataFrame([dict(bc1="A" * 20, tier="tier_1", A=10, B=20, pooled_AB=99)])
    key = pd.DataFrame(dict(sample_name=["A", "B"], group=[group, group]))
    return counts, key


@pytest.mark.parametrize("name", ["pooled_AB", "bc1", "tier", "oligo_name"])
def test_output_namespace_collision_is_refused(name):
    counts, key = fixture(name)
    original = counts.copy(deep=True)
    with pytest.raises(ValueError, match="Output column collision"):
        module().combine_technical_replicates(counts, key, "group")
    pd.testing.assert_frame_equal(counts, original)


@pytest.mark.parametrize("name", ["merged_AB", "A"])
def test_distinct_output_or_group_named_after_own_member_conserves_reads(name):
    counts, key = fixture(name)
    combined, correlations, outliers = module().combine_technical_replicates(counts, key, "group")
    assert combined.loc[0, name] == 30
    assert combined.loc[0, "pooled_AB"] == 99
    assert combined[[name, "pooled_AB"]].sum(axis=1).iloc[0] == 129
    assert combined.loc[0, "bc1"] == "A" * 20
    assert combined.loc[0, "tier"] == "tier_1"
    assert len(correlations) == len(outliers) == 1
    assert outliers.iloc[0]["reason"] == "insufficient_data"
