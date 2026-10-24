"""An amplicon window is checked against the reference, like the barcode flanks.

A window is a raw coordinate pair. If the barcode locus is re-pointed at another reference
and the windows are not, the filter silently looks in the wrong place. So a configured window
must carry the reference sequence it covers, and `validate_tracks` compares that sequence with
the FASTA every run, exactly as it already checks the N-run and its junction flanks.
"""
import pysam
import pytest

from magestic.pipelines.wgs.core.locus_config import LocusConfigError, validate_tracks
from trace_crispr.wgs.bc1_consensus import TrackConfig

LEFT, RIGHT = "ACCCTGTTAT", "GAGCATGACC"
PREFIX = "TCAATGGCTTGCGCTAACTGCGAACAGAGTGCCCTATGAAATAGGGGATACCGTTCGTATAG"[:60]
SUFFIX = "TTTCCAAACTGTCTCACGACGTTTTGAAACGGAAAGTCAATGGCTTGCGCTAACTGCGAACA"[:60]
CONTIG = PREFIX + LEFT + "N" * 20 + RIGHT + SUFFIX   # barcode at 1-based 71-90
LEFT_WINDOW, RIGHT_WINDOW = (5, 20), (101, 110)
LEFT_SEQ = CONTIG[LEFT_WINDOW[0] - 1:LEFT_WINDOW[1]]
RIGHT_SEQ = CONTIG[RIGHT_WINDOW[0] - 1:RIGHT_WINDOW[1]]


@pytest.fixture
def fasta(tmp_path):
    path = tmp_path / "ref.fa"
    path.write_text(">chrII\n" + CONTIG + "\n")
    pysam.faidx(str(path))
    return path


def track(left=LEFT_WINDOW, right=RIGHT_WINDOW):
    return TrackConfig("genome", "chrII", 71, 90, LEFT, RIGHT, True,
                       amplicon_left_window=left, amplicon_right_window=right)


def test_windows_whose_sequences_match_the_reference_pass(fasta):
    validate_tracks([track()], fasta, amplicon_window_sequences=(LEFT_SEQ, RIGHT_SEQ))


def test_no_windows_need_no_sequences(fasta):
    validate_tracks([track(left=None, right=None)], fasta)


def test_a_window_without_its_sequence_is_refused_even_by_the_default_call(fasta):
    with pytest.raises(LocusConfigError, match="amplicon_left_window"):
        validate_tracks([track()], fasta)
    with pytest.raises(LocusConfigError, match="amplicon_right_window"):
        validate_tracks([track()], fasta, amplicon_window_sequences=(LEFT_SEQ, None))


@pytest.mark.parametrize("shift", [1, -1])
def test_a_drifted_window_is_refused(fasta, shift):
    drifted = (LEFT_WINDOW[0] + shift, LEFT_WINDOW[1] + shift)
    with pytest.raises(LocusConfigError, match="drift"):
        validate_tracks([track(left=drifted)], fasta, amplicon_window_sequences=(LEFT_SEQ, RIGHT_SEQ))


def test_a_sequence_without_its_window_is_refused(fasta):
    with pytest.raises(LocusConfigError, match="amplicon_right_window_sequence"):
        validate_tracks([track(right=None)], fasta, amplicon_window_sequences=(LEFT_SEQ, RIGHT_SEQ))


def test_the_comparison_ignores_case(fasta):
    validate_tracks([track()], fasta, amplicon_window_sequences=(LEFT_SEQ.lower(), RIGHT_SEQ))
