"""support_tier is a LABEL, never a filter, in the two non-DAG code paths (planning 21_/28_, 2026-09-10)."""
import pandas as pd
import pytest


def test_run_pipeline_and_merge_label_not_filter():
    from magestic.pipelines.bc1_donor_bc0.pipeline.run_pipeline import aggregate_and_filter
    from magestic.pipelines.bc1_donor_bc0.core.merge import filter_by_quality
    df = pd.DataFrame({"bc1": ["A" * 20, "C" * 20, "G" * 20], "donor": ["T" * 40, "A" * 40, "C" * 40],
                       "bc0": ["G" * 10] * 3, "counts": [10, 5, 3], "total_pcr_replicates": [1, 2, 3]})
    a = aggregate_and_filter(df, min_reads=1, min_pcr_replicates=2, collapse_sequences=False,
                             show_progress=False)
    assert len(a) == 3
    assert a.set_index("bc1").support_tier.to_dict() == {"A" * 20: "provisional", "C" * 20: "confident",
                                                         "G" * 20: "confident"}
    m = df.rename(columns={"counts": "total_counts"})
    assert len(filter_by_quality(m, min_reads=1, min_pcr_replicates=2)) == 3
    assert len(filter_by_quality(m, min_reads=1, min_pcr_replicates=2, mode="filter")) == 2
    with pytest.raises(ValueError):
        filter_by_quality(m, mode="bogus")
