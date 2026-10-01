"""Regression gate: subpool reverse primers must be reverse-complemented at read.

The oligo carries the reverse complement of its subpool's reverse primer at its 3' end,
because the primer anneals to the bottom strand. The original generator did this at load
time; the port to the packaged CLI dropped it, which would have shipped subpools that
cannot be amplified with their intended primer.

This test pins the behavior against **real shipped data** rather than a fixture: the NRD1
SpCas9 pool that was actually ordered. If someone removes the reverse complement again,
this fails with a sequence comparison instead of a silent library failure months later.
"""

from __future__ import annotations

import pathlib

import pytest

from magestic.design._cli import _load_rev_primers

REPO = pathlib.Path(__file__).resolve().parents[2]
PRIMER_FILE = REPO / "NRD1" / "rev_subpool_priming_sequences.txt"
ORDERED_POOL = REPO / "NRD1" / "NRD1_SpCas9_NGG_guide_donor_200bp_oligo_pool.tsv"

#: Line 1 of the primer file, as written 5'->3'.
PRIMER_1_AS_WRITTEN = "TGCTGTCCAGTCGTCTATCG"
#: What every oligo of the real ordered NRD1 pool actually ends with.
ORDERED_POOL_3P_TAIL = "CGATAGACGACTGGACAGCA"


def test_rev_comp_is_applied_at_read():
    primers = _load_rev_primers(PRIMER_FILE)
    assert primers, "no primers loaded"
    assert primers[0] == ORDERED_POOL_3P_TAIL.lower()
    assert primers[0] != PRIMER_1_AS_WRITTEN.lower(), (
        "primer was loaded as written; oligos would not amplify with this subpool primer"
    )


def test_every_ordered_nrd1_oligo_ends_with_the_loaded_primer():
    """The decisive check: 2,308 real ordered oligos, one expected 3' tail."""
    if not ORDERED_POOL.exists():
        pytest.skip(f"{ORDERED_POOL} not present")

    primer_0 = _load_rev_primers(PRIMER_FILE)[0].upper()

    tails = set()
    n = 0
    with open(ORDERED_POOL) as fh:
        for line in fh:
            for field in line.rstrip("\n").split("\t"):
                if len(field) == 200:
                    tails.add(field.upper()[-len(primer_0):])
                    n += 1

    assert n == 2308, f"expected 2308 oligos in the ordered pool, found {n}"
    assert tails == {primer_0}, (
        f"ordered pool 3' tails {tails} do not match the loaded primer {primer_0}"
    )


def test_primer_file_is_not_already_reverse_complemented():
    """Guard the fix against being applied twice.

    If someone pre-reverse-complements the stored file, applying rev_comp at read would
    undo it. The stored file must contain the primer as written.
    """
    first = PRIMER_FILE.read_text().splitlines()[0].strip().upper()
    assert first == PRIMER_1_AS_WRITTEN, (
        "the stored primer file changed orientation; _load_rev_primers assumes it holds "
        "primers as written 5'->3'"
    )
