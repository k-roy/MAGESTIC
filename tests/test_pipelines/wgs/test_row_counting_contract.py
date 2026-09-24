"""Data-row counting for written TSVs: one rule for `count_data_rows` and `require_nonempty_output`.

A written table may carry a `# key=value` stamp, then a header that is either commented
(TRACE writes `#sample_name\\t...`) or plain (magestic's own tables), then data rows. Both
helpers used to find the header by LINE POSITION and were wrong in opposite directions:
`count_data_rows` counted a plain header as data when a stamp preceded it (audit W14b), and
`require_nonempty_output` subtracted a header that had been commented out, so a one-row
TRACE table read as "ZERO data rows" (audit W14a). The header is now recognised by content.
"""
import pytest

from magestic.pipelines.wgs.core.emit import count_data_rows
from magestic.pipelines.wgs.qc.pipeline_qc import QCFailure, require_nonempty_output

STAMP = "# synthetic_error_axis=disabled\n"
PLAIN = "sample\tbarcode\n"
COMMENTED = "#sample\tbarcode\n"
ROWS = "s1\tACGT\ns2\tTTGA\n"

CASES = [
    # (id, text, expected data rows)
    ("plain_header", PLAIN + ROWS, 2),
    ("commented_header", COMMENTED + ROWS, 2),
    ("stamp_then_plain_header", STAMP + PLAIN + ROWS, 2),
    ("stamp_then_commented_header", STAMP + COMMENTED + ROWS, 2),
    ("commented_header_one_row", COMMENTED + "s1\tACGT\n", 1),
    ("stamp_plain_header_one_row", STAMP + PLAIN + "s1\tACGT\n", 1),
    ("blank_lines_ignored", PLAIN + "\n" + ROWS + "\n", 2),
    ("plain_header_only", PLAIN, 0),
    ("commented_header_only", COMMENTED, 0),
    ("stamp_and_plain_header_only", STAMP + PLAIN, 0),
    ("stamp_and_commented_header_only", STAMP + COMMENTED, 0),
]


@pytest.mark.parametrize("case_id,text,expected", CASES, ids=[c[0] for c in CASES])
def test_count_data_rows_finds_the_header_by_content(tmp_path, case_id, text, expected):
    p = tmp_path / f"{case_id}.tsv"
    p.write_text(text)
    assert count_data_rows(p) == expected


@pytest.mark.parametrize("case_id,text,expected", [c for c in CASES if c[2] > 0],
                         ids=[c[0] for c in CASES if c[2] > 0])
def test_require_nonempty_output_counts_what_count_data_rows_counts(tmp_path, case_id, text, expected):
    p = tmp_path / f"{case_id}.tsv"
    p.write_text(text)
    assert require_nonempty_output(p, label=case_id) == expected == count_data_rows(p)


@pytest.mark.parametrize("case_id,text", [(c[0], c[1]) for c in CASES if c[2] == 0],
                         ids=[c[0] for c in CASES if c[2] == 0])
def test_header_only_table_is_still_a_loud_failure(tmp_path, case_id, text):
    p = tmp_path / f"{case_id}.tsv"
    p.write_text(text)
    with pytest.raises(QCFailure):
        require_nonempty_output(p, label=case_id)


def test_missing_table_is_a_loud_failure(tmp_path):
    with pytest.raises(QCFailure):
        require_nonempty_output(tmp_path / "absent.tsv", label="absent")
