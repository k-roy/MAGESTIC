"""
Step 02 core logic: Parse collapsed FASTQs into guide-donor-BC0 counts.

Processes collapsed FASTQ files from Illumina sequencing of plasmid libraries.
Extracts guide, donor, BC0, and SPS sequences using anchor-based parsing and
accumulates read counts per unique sequence combination.

Author: Kevin R. Roy
"""

import timeit
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pysam

from magestic.utils.collapsed_fastq import parse_collapsed_header_count
from magestic.utils.sequence_anchors import (
    extract_donor_bc0_fragment,
    find_anchor_from_end,
    find_re_site,
    DEFAULT_DONOR_BC0_FRAGMENT_LENGTH,
    DEFAULT_RE_SITE_MAX_MISMATCHES,
    RE_SITE_PATTERNS,
)


# ============================================================================
# Sequence Constants
# ============================================================================

@dataclass(frozen=True)
class SequenceConstants:
    """Fixed sequence lengths and anchor patterns for parsing guide-donor-BC0 reads."""
    guide_length: int = 20
    RE_site: str = "GTTTGAAGAGC"        # SpCas9/SpG restriction enzyme site
    RE_site_lbcas12a: str = "TTTCGAAGAGC"  # LbCas12a restriction enzyme site
    RE_site_length: int = 11
    RE_site_search_end: int = 150       # the site is searched in seq[:RE_site_search_end]
    RE_site_max_mismatches: int = DEFAULT_RE_SITE_MAX_MISMATCHES  # substitutions tolerated (1 == old min_score 85)
    msd_max_distance_from_end: int = 30
    bc0_length: int = 10
    SPS_length: int = 20
    msd: str = "AGGAAAACAGAC"           # Constant 3' anchor sequence
    msd_length: int = 12
    donor_bc0_fragment_length: int = 60


SEQ_CONST = SequenceConstants()


# ============================================================================
# Sample Keyfile
# ============================================================================

@dataclass
class SampleInfo:
    """Metadata for a single sequencing sample from the keyfile."""
    sample_number: str
    library_ID: str
    inner_primers: str
    outer_primers: str
    sequencing_dates: List[str]
    polymerase: str


def load_sample_keyfile(keyfile_path: Path) -> Tuple[Dict[str, List[SampleInfo]], Set[str]]:
    """
    Load the step-02 sample keyfile and return a per-library sample map.

    Keyfile format (tab-separated, one header row):
        sample_number  library_ID  inner_primers  outer_primers  sequencing_dates  polymerase

    The ``sequencing_dates`` column may be comma-separated for samples sequenced
    on multiple dates.

    Parameters
    ----------
    keyfile_path : Path
        Path to the TSV sample keyfile.

    Returns
    -------
    library_to_samples : dict[str, list[SampleInfo]]
        Maps each library_ID to the list of samples associated with it.
    all_dates : set[str]
        Union of all sequencing dates found in the keyfile.  Used to
        pre-initialise count columns so every output row has the same columns.
    """
    library_to_samples: Dict[str, List[SampleInfo]] = {}
    all_dates: Set[str] = set()

    with open(keyfile_path) as fh:
        header = fh.readline().strip().split("\t")
        if len(header) < 5:
            raise ValueError(
                f"Expected at least 5 columns in {keyfile_path}, got {len(header)}"
            )

        for line_num, line in enumerate(fh, start=2):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 5:
                continue

            sequencing_dates = [d.strip() for d in fields[4].split(",") if d.strip()]
            all_dates.update(sequencing_dates)

            sample = SampleInfo(
                sample_number=fields[0],
                library_ID=fields[1],
                inner_primers=fields[2],
                outer_primers=fields[3],
                sequencing_dates=sequencing_dates,
                polymerase=fields[5] if len(fields) > 5 else "Q5",
            )
            library_to_samples.setdefault(fields[1], []).append(sample)

    return library_to_samples, all_dates


# ============================================================================
# Sequence Utilities
# ============================================================================

def rev_comp(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(
        str.maketrans("ACGTacgtRYMKrymkVBHDvbhdNn ", "TGCAtgcaYRKMyrkmBVDHbvdhNn ")
    )[::-1]


def calculate_avg_phred(qual_string: str) -> float:
    """Return the mean Phred quality score for an ASCII quality string."""
    if not qual_string:
        return 0.0
    return sum(ord(c) - 33 for c in qual_string) / len(qual_string)


# ============================================================================
# Sequence Parsing
# ============================================================================

PARSE_STAT_KEYS = (
    "reads", "no_msd", "no_re_site", "bad_frame", "parsed",
    "re_site_SpCas9_exact", "re_site_LbCas12a_exact",
    "re_site_SpCas9_mm1", "re_site_LbCas12a_mm1", "re_site_ambiguous_mm1",
)
"""Keys of the per-read parse accounting written to every step-02 provenance sidecar.
``reads = no_msd + no_re_site + bad_frame + parsed`` and ``parsed = sum(re_site_*)``."""


def new_parse_stats() -> Dict[str, int]:
    """A zeroed parse-accounting dict (every key present so sidecars have a fixed schema)."""
    return {k: 0 for k in PARSE_STAT_KEYS}


def parse_sequence(
    seq: str,
    cas_type: str = "SpCas9",
    stats: Optional[Dict[str, int]] = None,
) -> Tuple[str, str, str, str, str, str, str]:
    """
    Parse a read into its guide-donor-BC0 components using anchor-based extraction.

    Anchors (in order of reliability):
    1. ``msd`` found from the 3' end — most stable anchor.
    2. RE (cloning) site found from the 5' end, across ALL nuclease patterns at once
       (:func:`magestic.utils.sequence_anchors.find_re_site`): exact first, then
       substitutions-only within ``SEQ_CONST.RE_site_max_mismatches``.

    Positions relative to anchors::

        5'─ leader ─[guide]─[RE_site]─[donor]─[SPS]─[bc0]─[msd]─3'

    Parameters
    ----------
    seq : str
        Full read sequence (forward orientation).
    cas_type : str
        Retained for API compatibility only. The search is nuclease-agnostic: every
        pattern is tried exactly before any is tried approximately, so the caller's
        expectation can no longer steer the anchor onto the wrong window (GDB-01: the
        SpCas9-first search mis-cut every G-ending Cas12a guide by one base).
    stats : dict, optional
        Accounting dict from :func:`new_parse_stats`; incremented in place per read.

    Returns
    -------
    tuple of (leader, guide, donor, SPS, bc0, msd, donor_bc0_fragment)
        All elements are ``"NA"`` if the anchors cannot be located.
    """
    if stats is not None:
        stats["reads"] += 1
    msd_pos = find_anchor_from_end(seq, SEQ_CONST.msd,
                                   max_distance_from_end=SEQ_CONST.msd_max_distance_from_end)
    if msd_pos is None:
        if stats is not None:
            stats["no_msd"] += 1
        return ("NA",) * 7

    hit = find_re_site(seq, end=min(SEQ_CONST.RE_site_search_end, len(seq)),
                       max_mismatches=SEQ_CONST.RE_site_max_mismatches)
    if hit is None:
        if stats is not None:
            stats["no_re_site"] += 1
        return ("NA",) * 7
    re_pos = hit.pos

    guide_start = re_pos - SEQ_CONST.guide_length
    donor_start = re_pos + SEQ_CONST.RE_site_length
    bc0_end = msd_pos
    bc0_start = bc0_end - SEQ_CONST.bc0_length
    sps_end = bc0_start
    sps_start = sps_end - SEQ_CONST.SPS_length
    donor_end = sps_start

    if guide_start < 0 or donor_end <= donor_start:
        if stats is not None:
            stats["bad_frame"] += 1
        return ("NA",) * 7

    if stats is not None:
        stats["parsed"] += 1
        key = f"re_site_{hit.kind}_" + ("exact" if hit.mismatches == 0 else f"mm{hit.mismatches}")
        stats[key] = stats.get(key, 0) + 1

    leader = seq[:guide_start]
    guide = seq[guide_start:re_pos]
    donor = seq[donor_start:donor_end]
    SPS = seq[sps_start:sps_end]
    bc0 = seq[bc0_start:bc0_end]
    msd = seq[msd_pos: msd_pos + SEQ_CONST.msd_length]

    donor_bc0_fragment = extract_donor_bc0_fragment(
        seq,
        msd_pos=msd_pos,
        msd_pattern=SEQ_CONST.msd,
        fragment_length=SEQ_CONST.donor_bc0_fragment_length,
        bc0_length=SEQ_CONST.bc0_length,
        sps_length=SEQ_CONST.SPS_length,
    )
    if donor_bc0_fragment is None:
        donor_bc0_fragment = "NA"

    return leader, guide, donor, SPS, bc0, msd, donor_bc0_fragment


def parse_parameters() -> Dict[str, object]:
    """Every constant that affects parse_sequence's output, for provenance sidecars."""
    return {
        "guide_length": SEQ_CONST.guide_length,
        "re_site_patterns": {k: v for k, v in RE_SITE_PATTERNS},
        "re_site_length": SEQ_CONST.RE_site_length,
        "re_site_search_end": SEQ_CONST.RE_site_search_end,
        "re_site_max_mismatches": SEQ_CONST.RE_site_max_mismatches,
        "re_site_search": "all patterns exact first (earliest position wins), then substitutions-only "
                          "(Hamming) within re_site_max_mismatches, lowest mismatches then earliest position",
        "msd": SEQ_CONST.msd,
        "msd_max_distance_from_end": SEQ_CONST.msd_max_distance_from_end,
        "msd_search": "exact then rapidfuzz ratio >= 85 (find_anchor_from_end; unchanged)",
        "bc0_length": SEQ_CONST.bc0_length,
        "sps_length": SEQ_CONST.SPS_length,
        "donor_bc0_fragment_length": SEQ_CONST.donor_bc0_fragment_length,
        "reverse_primer_pairs": sorted(_REVERSE_PRIMER_PAIRS),
        "reverse_read_handling": "reverse-complement; KR1882-KR1953 reads get 'GAC' appended",
    }


# ============================================================================
# FASTQ Processing
# ============================================================================

# Reverse primer pairs require rc of the read and optional 3-nt extension.
_REVERSE_PRIMER_PAIRS = {"KR1884-KR1953", "KR1882-KR1953"}


def process_collapsed_fastq(
    fastq_path: Path,
    inner_primers: str,
    sequencing_run_date: str,
    library_ID: str,
    bc0_counts: Dict[tuple, Dict[str, int]],
    all_dates: Set[str],
    min_quality: int = 0,
    quality_stats: Optional[Dict[str, int]] = None,
    parse_stats: Optional[Dict[str, int]] = None,
) -> Tuple[int, int]:
    """
    Process one collapsed FASTQ file and accumulate per-unique-sequence counts.

    Parameters
    ----------
    fastq_path : Path
        Collapsed FASTQ carrying either supported count-header format. Format is
        resolved per record by the canonical parser; unknown/zero counts raise.
    inner_primers : str
        Primer pair identifier (e.g. ``"KR1884-KR1953"``).  Determines
        whether the read should be reverse-complemented before parsing.
    sequencing_run_date : str
        Date string (e.g. ``"20240215"``); last 4 digits become the column suffix.
    library_ID : str
        Library identifier used only for contextual logging.
    bc0_counts : dict
        Accumulator mapping ``bc0_info`` tuples to per-column count dicts.
        **Modified in place.**
    all_dates : set[str]
        All sequencing dates in the keyfile.  Used to pre-initialise count
        columns so every row in the output has the same schema.
    min_quality : int
        Minimum mean Phred quality to keep a read.  0 disables filtering.
    quality_stats : dict, optional
        Accumulator for quality-filtering statistics.  **Modified in place.**
    parse_stats : dict, optional
        Per-read anchor accounting (:func:`new_parse_stats`).  **Modified in place.**

    Returns
    -------
    total_unique_seqs : int
    total_counts : int
    """
    total_unique_seqs = 0
    total_counts = 0
    filtered_low_quality = 0
    filtered_counts = 0

    is_reverse = inner_primers in _REVERSE_PRIMER_PAIRS
    seq_direction = "rev" if is_reverse else "fwd"
    counts_col = f"{seq_direction}_{sequencing_run_date[-4:]}"

    with pysam.FastxFile(str(fastq_path)) as fin:
        for entry in fin:
            total_unique_seqs += 1
            count, _ = parse_collapsed_header_count(entry.name)
            total_counts += count
            seq = entry.sequence
            qual = entry.quality or ""

            if min_quality > 0 and qual:
                if calculate_avg_phred(qual) < min_quality:
                    filtered_low_quality += 1
                    filtered_counts += count
                    continue

            if is_reverse:
                seq = rev_comp(seq)
                qual = qual[::-1]
                if inner_primers == "KR1882-KR1953":
                    seq += "GAC"
                    qual += "EEE"

            leader, guide, donor, SPS, bc0, msd, donor_bc0_frag = parse_sequence(seq, stats=parse_stats)
            bc0_info = (leader, guide, donor, SPS, bc0, msd, donor_bc0_frag, seq)

            if bc0_info not in bc0_counts:
                bc0_counts[bc0_info] = {}
                for run_date in all_dates:
                    for direction in ("fwd", "rev"):
                        bc0_counts[bc0_info][f"{direction}_{run_date[-4:]}"] = 0

            bc0_counts[bc0_info][counts_col] += count

    if quality_stats is not None:
        quality_stats["filtered_low_quality"] = (
            quality_stats.get("filtered_low_quality", 0) + filtered_low_quality
        )
        quality_stats["filtered_counts"] = (
            quality_stats.get("filtered_counts", 0) + filtered_counts
        )

    return total_unique_seqs, total_counts


def calculate_totals(bc0_counts: Dict[tuple, Dict[str, int]]) -> None:
    """
    Add ``total_counts``, ``total_fwd_counts``, ``total_rev_counts`` to each entry.

    Modifies *bc0_counts* in place.
    """
    for counts_dict in bc0_counts.values():
        total = sum(counts_dict.values())
        total_fwd = sum(v for k, v in counts_dict.items() if k.startswith("fwd_"))
        total_rev = sum(v for k, v in counts_dict.items() if k.startswith("rev_"))
        counts_dict["total_counts"] = total
        counts_dict["total_fwd_counts"] = total_fwd
        counts_dict["total_rev_counts"] = total_rev


def write_counts_file(
    output_path: Path,
    library_ID: str,
    bc0_counts: Dict[tuple, Dict[str, int]],
    run_columns: Optional[List[str]] = None,
) -> None:
    """
    Write the accumulated counts dictionary to a tab-separated file.

    Parameters
    ----------
    output_path : Path
        Destination TSV file (created or overwritten).
    library_ID : str
        Written as the first column of every row.
    bc0_counts : dict
        Maps ``(leader, guide, donor, SPS, bc0, msd, donor_bc0_fragment, seq)``
        tuples to per-column count dicts.
    run_columns : list[str], optional
        Per-run columns to preserve when there are no retained sequences.
    """
    sample_counts = next(iter(bc0_counts.values()), {})
    total_cols = ["total_counts", "total_fwd_counts", "total_rev_counts"]
    per_run_cols = sorted(k for k in (sample_counts if bc0_counts else (run_columns or []))
                          if k not in total_cols)
    seq_cols = ["leader", "guide", "donor", "SPS", "bc0", "msd", "donor_bc0_fragment", "seq"]

    header = ["library_ID"] + total_cols + per_run_cols + seq_cols

    with open(output_path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for bc0_info, counts_dict in bc0_counts.items():
            row = [library_ID]
            row.extend(str(counts_dict.get(col, 0)) for col in total_cols)
            row.extend(str(counts_dict.get(col, 0)) for col in per_run_cols)
            row.extend(str(val) for val in bc0_info)
            fh.write("\t".join(row) + "\n")


# ============================================================================
# Per-Library Orchestration
# ============================================================================

def validate_collapsed_fastq_inputs(samples: List[SampleInfo], collapsed_fastq_dir: Path) -> None:
    """Refuse missing expected FASTQs before any result is published.

    An existing zero-record FASTQ is valid. A missing file is not an empty sample.
    Callers can preflight all selected libraries together before processing any.
    """
    missing = []
    for sample in samples:
        for seq_date in sample.sequencing_dates:
            path = collapsed_fastq_dir / f"{seq_date}_{sample.outer_primers}_{sample.sample_number}_collapsed.fastq"
            if not path.is_file():
                missing.append(str(path))
    if missing:
        raise FileNotFoundError("Missing expected collapsed FASTQ files: " + ", ".join(missing))


def process_library(
    library_ID: str,
    samples: List[SampleInfo],
    collapsed_fastq_dir: Path,
    counts_dir: Path,
    all_dates: Set[str],
    min_quality: int = 0,
) -> Dict[str, object]:
    """
    Process all sequencing samples for one library and write a counts TSV.

    Parameters
    ----------
    library_ID : str
        Library identifier used for output file naming.
    samples : list[SampleInfo]
        All samples belonging to this library.
    collapsed_fastq_dir : Path
        Directory containing collapsed FASTQ files from Step 01.
    counts_dir : Path
        Directory where the output ``{library_ID}_guide_donor_bc0_counts.tsv``
        is written.
    all_dates : set[str]
        All sequencing dates across the entire keyfile (column schema).
    min_quality : int
        Minimum mean Phred quality filter (0 = no filtering).

    Returns
    -------
    dict
        Processing statistics: files processed, unique sequences, total counts,
        missing files, elapsed seconds, output file path, filtered counts.
    """
    validate_collapsed_fastq_inputs(samples, collapsed_fastq_dir)
    start_time = timeit.default_timer()
    bc0_counts: Dict[tuple, Dict[str, int]] = {}
    quality_stats: Dict[str, int] = {"filtered_low_quality": 0, "filtered_counts": 0}
    parse_stats: Dict[str, int] = new_parse_stats()
    files_processed: List[Dict[str, object]] = []

    total_files = 0
    total_seqs = 0
    total_counts = 0
    missing_files: List[str] = []

    for sample in samples:
        for seq_date in sample.sequencing_dates:
            sample_name = f"{seq_date}_{sample.outer_primers}_{sample.sample_number}"
            fastq_path = collapsed_fastq_dir / f"{sample_name}_collapsed.fastq"

            total_files += 1
            n_seqs, n_counts = process_collapsed_fastq(
                fastq_path,
                sample.inner_primers,
                seq_date,
                library_ID,
                bc0_counts,
                all_dates,
                min_quality=min_quality,
                quality_stats=quality_stats,
                parse_stats=parse_stats,
            )
            total_seqs += n_seqs
            total_counts += n_counts
            files_processed.append({"path": str(fastq_path), "sample": sample.sample_number,
                                    "sequencing_date": seq_date, "inner_primers": sample.inner_primers,
                                    "unique_seqs": n_seqs, "reads": n_counts})

    calculate_totals(bc0_counts)

    output_path = counts_dir / f"{library_ID}_guide_donor_bc0_counts.tsv"
    run_columns = [f"{direction}_{date[-4:]}" for date in all_dates for direction in ("fwd", "rev")]
    write_counts_file(output_path, library_ID, bc0_counts, run_columns=run_columns)

    return {
        "library_ID": library_ID,
        "files_processed": total_files,
        "unique_sequences": len(bc0_counts),
        "total_seqs_in_files": total_seqs,
        "total_counts": total_counts,
        "missing_files": missing_files,
        "elapsed_seconds": timeit.default_timer() - start_time,
        "output_file": str(output_path),
        "filtered_low_quality": quality_stats["filtered_low_quality"],
        "filtered_counts": quality_stats["filtered_counts"],
        "parse_stats": parse_stats,
        "files": files_processed,
    }
