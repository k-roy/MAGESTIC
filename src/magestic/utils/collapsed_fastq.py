"""
Canonical reader for MAGESTIC collapsed FASTQ — ONE definition of "how many
reads is this record worth".

A collapsed FASTQ carries the number of identical merged reads that produced
each record IN THE READ NAME. Two producers exist in this codebase and they
write DIFFERENT names:

  * ``collapse_identical_merged_reads_min_counts.py`` (2019, repo root)
        read_name = '@' + str(counts)          ->  ``@12``
    Used for the 2021, 2022, 2022-10 REDI, 2024 cAMP and 2024 REDI-manual
    collapsed trees -- 15,500 of the 19,765 in-scope collapsed FASTQ.

  * ``bc1_from_screen/snakemake/01_trim_merge_collapse.sh`` (current, inline)
        f"@seq_{seq_num}_count={count}"        ->  ``@seq_1_count=208``
    Used for yL406 (2024) and yL437/yL442 (2025).

Both eras coexist INSIDE one directory (yL406's ``bc1_from_screen/collapsed_fastq``
holds 3,528 files of both kinds), so the format is resolved PER RECORD, never
per run, per directory or per era.

Why this module exists (planning/210_bc1_recall_code_gate.md):
    ``02_parse_bc1_reads.py`` recognized only the ``count=`` form and silently
    fell back to ``count = 1`` for everything else. On the legacy ``@N`` form
    that discards the collapse counts entirely. Measured on
    ``20210820_plate_2B_sample_25_collapsed.fastq``: true depth 79,474 reads,
    reported depth 57,306 (the record count), and the top barcode's abundance
    read 4 instead of 12 -- the top-10 lists shared no members. The distortion
    SATURATES rather than scaling, so it is worst on the abundant lineages.

    It also carried an ``elif "-" in header`` heuristic that would turn any
    future dash-bearing read name into a plausible integer. That is removed.

The contract here is the opposite of the old one: an unrecognized header is an
ERROR that names the offending line. A caller that genuinely wants records
weighted 1 must say so with ``allow_unknown_as_one=True``, and the returned
format counter still records how many took that path.
"""

import gzip
import re
from collections import Counter
from pathlib import Path
from typing import Iterator, Tuple, Union

#: Recognized read-name formats, IN PRIORITY ORDER. ``count=`` is tried first
#: because a ``@seq_1_count=208`` name also contains digits. Extend this tuple
#: -- do not add another ``elif`` at a call site.
COLLAPSED_HEADER_FORMATS = (
    # @seq_1_count=208  (01_trim_merge_collapse.sh). Anchored at the END so a
    # producer that ever appends a field cannot have a LATER number silently
    # read as the count -- it raises instead.
    ("count_kv", re.compile(r"count=(\d+)\s*$")),
    # @12              (collapse_identical_merged_reads_min_counts.py)
    ("bare_int", re.compile(r"^@?(\d+)\s*$")),
)

UNKNOWN_FORMAT = "unknown"


class CollapsedHeaderError(ValueError):
    """A collapsed-FASTQ read name matched no known count format."""


def parse_collapsed_header_count(header: str) -> Tuple[int, str]:
    """
    Return ``(count, format_name)`` for one collapsed-FASTQ read name.

    Raises:
        CollapsedHeaderError: naming the offending header, if no format matches.
    """
    h = header.strip()
    for name, pattern in COLLAPSED_HEADER_FORMATS:
        m = pattern.search(h)
        if m:
            count = int(m.group(1))
            # A zero-count record is not a thing either producer writes (both
            # filter on ``count >= min_counts``, min 1). If one appears, the
            # read name does not mean what we think it means -- and a
            # 0-weighted record would contribute nothing while still counting
            # toward n_records, depressing the collapse factor. That is the
            # same silent-miscount class this module exists to prevent.
            if count < 1:
                raise CollapsedHeaderError(
                    "collapsed-FASTQ read name %r parsed as count=%d via the "
                    "%r format. A collapsed record represents at least one "
                    "read, so this read name does not mean what the parser "
                    "assumes." % (header, count, name)
                )
            return count, name
    raise CollapsedHeaderError(
        "unrecognized collapsed-FASTQ read name: %r. Known formats: %s. "
        "A collapsed read name MUST carry the number of reads collapsed into "
        "it; defaulting it to 1 silently discards abundance "
        "(planning/210_bc1_recall_code_gate.md)."
        % (header, ", ".join(n for n, _ in COLLAPSED_HEADER_FORMATS))
    )


def open_collapsed(path: Union[str, Path]):
    """Open a collapsed FASTQ, transparently gzipped."""
    return (gzip.open if str(path).endswith(".gz") else open)(path, "rt")


def iter_collapsed_fastq(
    path: Union[str, Path],
    limit: int = None,
    allow_unknown_as_one: bool = False,
    format_counter: Counter = None,
) -> Iterator[Tuple[str, int, str, str]]:
    """
    Yield ``(sequence, count, quality, format_name)`` per collapsed record.

    Args:
        path: collapsed FASTQ, plain or ``.gz``.
        limit: stop after this many RECORDS (not reads). ``None`` = whole file.
        allow_unknown_as_one: weight an unrecognized read name as 1 instead of
            raising. Opt-in only, and every such record is counted as
            ``UNKNOWN_FORMAT`` in ``format_counter``.
        format_counter: a ``Counter`` the caller supplies to receive the
            per-record format tally, for the provenance sidecar.

    Raises:
        CollapsedHeaderError: on an unrecognized read name, unless
            ``allow_unknown_as_one``.
    """
    n = 0
    with open_collapsed(path) as f:
        while limit is None or n < limit:
            header = f.readline()
            if not header:
                break
            header = header.strip()
            if not header:
                break
            sequence = f.readline().strip()
            f.readline()                      # '+' line
            quality = f.readline().strip()

            try:
                count, fmt = parse_collapsed_header_count(header)
            except CollapsedHeaderError as exc:
                if not allow_unknown_as_one:
                    raise CollapsedHeaderError(
                        "%s -- in %s, record %d" % (exc, path, n + 1)
                    ) from None
                count, fmt = 1, UNKNOWN_FORMAT

            if format_counter is not None:
                format_counter[fmt] += 1
            n += 1
            yield sequence, count, quality, fmt


def _self_test():
    """Assertions on the contract. Run as ``python -m magestic.utils.collapsed_fastq``."""
    assert parse_collapsed_header_count("@seq_1_count=208") == (208, "count_kv")
    assert parse_collapsed_header_count("@12") == (12, "bare_int")
    assert parse_collapsed_header_count("12") == (12, "bare_int")
    # ordering: a count= name also contains digits, and must not read as bare_int
    assert parse_collapsed_header_count("@seq_7_count=3")[0] == 3
    for bad in ("@0", "@seq_1_count=0", "@read_abc/1", "@", "@seq_1_count=208 extra",
                "@M01234:5:000-ABCDE:1:1101:15589:1332"):
        try:
            parse_collapsed_header_count(bad)
        except CollapsedHeaderError:
            pass
        else:
            raise AssertionError("should have raised on %r" % (bad,))
    print("collapsed_fastq self-test OK")


if __name__ == "__main__":
    _self_test()
