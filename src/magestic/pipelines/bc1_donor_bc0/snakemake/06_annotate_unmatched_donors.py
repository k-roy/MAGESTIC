#!/usr/bin/env python3
"""
Annotate unmatched donors with variant information via genome alignment.

This script processes bc1 clones whose donors don't match any designed oligo,
attempting to:
1. Use existing gdb columns from bc1 reference table (guide, gdb_oligo_name, etc.)
2. Find closest designed donor by edit distance (using yeast_sublibrary-specific oligo designs)
3. Align donors to reference genome using triple-aligner approach
4. Extract and annotate variants
5. Detect chimeric donors (potential double-edits)
6. Classify rescue status for hit calling

PARALLELIZATION STRATEGY:
- Large reference files (oligo designs) are loaded ONCE in parent
- Data is stored in module-level globals for true copy-on-write memory sharing after fork
- bc1 reference table is split into chunks for parallel processing
- Each worker inherits shared data via fork (no pickling overhead)
- Progress is tracked via per-chunk log files for monitoring
- Results are merged at the end

INCORPORATION PROBABILITY MODEL (from Li et al., Science Advances, 2025):
- Distance from adjacent variant determines incorporation rate:
  - ≤2 bp: ~100%
  - 3-5 bp: ~80-90%
  - 10 bp: ~50%
  - 20 bp: ~10%
  - 30 bp: ~5%

Usage:
    python 06_annotate_unmatched_donors.py \\
        --input <bc1_reference_table.tsv> \\
        --output-dir <output_directory> \\
        --spg-oligo-design <2024_spg_oligo_design.tsv> \\
        --spcas9-oligo-design <2021_spcas9_oligo_design.tsv> \\
        --reference-genome <MAGESTIC_background_strain.fasta> \\
        --annotation-gff <MAGESTIC_background_strain_annotations.gff> \\
        [--harmonized-oligo-key <harmonized_oligo_key.tsv>] \\
        [--n-workers 8] \\
        [--n-chunks 8]

Author: Kevin R. Roy
"""

from pathlib import Path
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple, Set
from dataclasses import dataclass, field
from collections import defaultdict
import subprocess
import tempfile
import os
import sys
import time
import argparse
from datetime import datetime, timedelta
from functools import partial
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import logging

try:
    from magestic.utils.variant_annotation import (
        classify_variant_type,
        compute_indel_length,
        extract_variant_positions,
        format_aa_change,
        CODON_VARIANT_TYPES,
    )
    HAS_VARIANT_ANNOTATION = True
except ImportError:
    HAS_VARIANT_ANNOTATION = False
    print("Warning: variant_annotation utils not found. Using local classification.")

# Optional imports with fallbacks
try:
    import Levenshtein
    HAS_LEVENSHTEIN = True
except ImportError:
    HAS_LEVENSHTEIN = False
    print("Warning: Levenshtein not installed. Using slower edit distance calculation.")

try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False
    print("Warning: pysam not installed. Alignment parsing will fail.")

# =============================================================================
# Configuration
# =============================================================================

# Parameters
MAX_EDIT_DISTANCE = 20  # Max ED to consider for closest donor
RESCUE_DISTANCE_THRESHOLD = 20  # bp from designed edit - beyond this, delta variants unlikely incorporated
MIN_MAPQ = 20  # Minimum mapping quality

# Parallelization parameters
DEFAULT_N_WORKERS = 8  # Default number of parallel workers
DEFAULT_N_CHUNKS = 8   # Default number of file chunks (can be different from workers)

# Oligo structure constants
OLIGO_GUIDE_START = 20
OLIGO_GUIDE_END = 40
OLIGO_DONOR_START = 51
OLIGO_DONOR_END = 180
DONOR_LENGTH = 129
MIDDLE_DONOR_START = 30
MIDDLE_DONOR_END = 90
DONOR_BC0_DONOR_LENGTH = 60  # The donor portion in donor_bc0_fragment is 60 bp

# Yeast_sublibrary to oligo design mapping
# NOTE: SpG_NGG uses SpCas9 with SpG nuclease, so uses 2021 oligo design
YEAST_SUBLIBRARY_TO_DESIGN = {
    "SpG_NGNG": "spg",    # 2024 SpG design
    "SpG_NGNH": "spg",    # 2024 SpG design
    "SpG_NGG": "spcas9",  # 2021 SpCas9 design (NOT SpG!)
    "SpCas9_NGG": "spcas9",
    "LbCas12a": "spcas9",     # 2021 pool — spcas9_lookup contains LbCas12a oligos
    "impLbCas12a": "spcas9",  # 2021 pool — spcas9_lookup contains impLbCas12a oligos
}

# V_library to oligo design — fallback when yeast_sublibrary is NaN
# HANDOFF_20260527 forward-port: used for 2021 QTL screens where step05
# trafo_key was not provided, so yeast_sublibrary cannot be derived.
V_LIBRARY_TO_DESIGN = {
    "V634": "spcas9",  # SpCas9/SpG NGG, 2021 pool
    "V450": "spcas9",  # LbCas12a CRISPEY, 2021 pool
    "V454": "spcas9",  # LbCas12a non-editing controls, 2021 pool
}


def _resolve_design(yeast_sublibrary: str) -> str:
    """Resolve oligo design from yeast_sublibrary or V_library key.

    HANDOFF_20260527 forward-port: NaN guard + V_library fallback. Without this,
    all bc1s routed to default 'spg' k-mer index → 100% tier_5_no_match for
    SpCas9/LbCas12a donors during rescue (V450/V454/V634).
    """
    design = YEAST_SUBLIBRARY_TO_DESIGN.get(yeast_sublibrary)
    if design is None:
        design = V_LIBRARY_TO_DESIGN.get(yeast_sublibrary, "spg")
    return design


# =============================================================================
# GFF Loading and Variant Consequence Annotation
# =============================================================================

@dataclass
class CdsInfo:
    """Information about a CDS region."""
    chrom: str
    start: int
    end: int
    strand: str
    systematic_name: str
    common_name: str


def load_gff_cds_info(gff_path: Path) -> Dict[str, List[CdsInfo]]:
    """
    Load CDS features from GFF file and build chromosome-indexed lookup.

    Args:
        gff_path: Path to GFF annotation file

    Returns:
        Dict mapping chromosome -> list of CdsInfo, sorted by start position
    """
    cds_by_chrom: Dict[str, List[CdsInfo]] = defaultdict(list)

    if not gff_path.exists():
        print(f"Warning: GFF file not found: {gff_path}")
        return {}

    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            if feature_type != "CDS":
                continue

            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]

            # Parse attributes
            attrs = {}
            for attr in parts[8].split(";"):
                if "=" in attr:
                    key, value = attr.split("=", 1)
                    attrs[key.strip()] = value.strip()

            systematic_name = attrs.get("Name", attrs.get("ID", ""))
            common_name = attrs.get("gene", systematic_name)

            cds_by_chrom[chrom].append(CdsInfo(
                chrom=chrom,
                start=start,
                end=end,
                strand=strand,
                systematic_name=systematic_name,
                common_name=common_name,
            ))

    # Sort each chromosome's CDS by start position
    for chrom in cds_by_chrom:
        cds_by_chrom[chrom].sort(key=lambda x: x.start)

    total_cds = sum(len(v) for v in cds_by_chrom.values())
    print(f"  Loaded {total_cds:,} CDS features from GFF")

    return dict(cds_by_chrom)


def find_overlapping_cds(
    chrom: str,
    pos: int,
    cds_lookup: Dict[str, List[CdsInfo]]
) -> Optional[CdsInfo]:
    """Find CDS that overlaps with a given position."""
    if chrom not in cds_lookup:
        return None

    cds_list = cds_lookup[chrom]
    for cds in cds_list:
        if cds.start <= pos <= cds.end:
            return cds
        if cds.start > pos:
            break

    return None


def annotate_variant_consequence(
    variants: List[Tuple[str, int, str, str]],
    cds_lookup: Dict[str, List[CdsInfo]]
) -> Tuple[str, Optional[str], Optional[str]]:
    """
    Annotate a list of variants with their coding consequence.

    Args:
        variants: List of (chrom, pos, ref, alt) tuples
        cds_lookup: CDS lookup from load_gff_cds_info()

    Returns:
        Tuple of (codon_variant_type, systematic_gene_name, common_gene_name)
    """
    if not variants:
        return "non_coding", None, None

    # Compute net indel length
    net_indel = 0
    for chrom, pos, ref, alt in variants:
        net_indel += len(alt) - len(ref)

    # Find overlapping CDS for any variant
    overlapping_cds = None
    for chrom, pos, ref, alt in variants:
        cds = find_overlapping_cds(chrom, pos, cds_lookup)
        if cds:
            overlapping_cds = cds
            break

    if not overlapping_cds:
        return "non_coding", None, None

    # Classify based on indel length
    if net_indel != 0:
        if HAS_VARIANT_ANNOTATION:
            variant_type = classify_variant_type(
                ref_aas=None,
                alt_aas=None,
                indel_length=net_indel,
                cds_start=overlapping_cds.start,
                cds_end=overlapping_cds.end,
            )
        else:
            # Fallback classification
            if net_indel % 3 == 0:
                variant_type = "in_frame_indel"
            else:
                variant_type = "frameshift"
    else:
        variant_type = "coding_SNP"

    return variant_type, overlapping_cds.systematic_name, overlapping_cds.common_name


# =============================================================================
# Module-level globals for copy-on-write memory sharing across workers
# =============================================================================
_WORKER_OLIGO_LOOKUP_SPG: Dict = {}
_WORKER_OLIGO_LOOKUP_SPCAS9: Dict = {}
_WORKER_INDIVIDUAL_VARIANTS_LOOKUP: Dict = {}
_WORKER_REFERENCE_FASTA: Path = None
_WORKER_CDS_LOOKUP: Dict = {}
_WORKER_OUTPUT_DIR: Path = None
_WORKER_PROGRESS_DIR: Path = None


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class AnnotationResult:
    """Results of annotating an unmatched donor."""
    bc1: str
    yeast_sublibrary: str = ""

    # Donor-bc0 lookup results (from existing columns)
    gdb_lookup_found: bool = False
    gdb_lookup_guide: Optional[str] = None
    gdb_lookup_donor: Optional[str] = None
    gdb_lookup_oligo_name: Optional[str] = None
    gdb_lookup_perfect_donor: bool = False
    gdb_lookup_perfect_middle: bool = False

    # Closest designed donor
    closest_oligo_name: Optional[str] = None
    closest_donor_edit_distance: int = -1
    closest_middle_donor_edit_distance: int = -1
    inferred_guide: Optional[str] = None

    # Alignment results
    alignment_chrom: Optional[str] = None
    alignment_start: int = -1
    alignment_end: int = -1
    alignment_mapq: int = 0
    alignment_score: int = -1
    aligner_used: Optional[str] = None

    # Variant annotation
    observed_individual_variants: Optional[str] = None
    designed_individual_variants: Optional[str] = None
    delta_variants: Optional[str] = None
    delta_variant_count: int = 0
    systematic_ORF_name: Optional[str] = None
    gene_name: Optional[str] = None
    aa_change: Optional[str] = None
    codon_variant_type: Optional[str] = None

    # Chimera detection
    is_chimera: bool = False
    chimera_left_parent: Optional[str] = None
    chimera_right_parent: Optional[str] = None
    chimera_junction_pos: int = -1
    chimera_gap_bp: int = -1

    # Extended donor detection
    is_extended: bool = False
    extension_non_bc0_side_bp: int = 0
    extension_bc0_side_bp: int = 0
    observed_donor_length: int = 0

    # Rescue assessment
    rescue_status: Optional[str] = None
    closest_delta_to_designed: int = -1
    incorporation_probability: float = 0.0
    confidence_tier: Optional[str] = None


@dataclass
class OligoInfo:
    """Information about a designed oligo."""
    oligo_name: str
    guide: str
    donor: str
    middle_donor: str
    oligo_seq: str

    systematic_ORF_name: Optional[str] = None
    gene_name: Optional[str] = None
    chrom: Optional[str] = None
    strand: Optional[str] = None
    pam_pos: Optional[int] = None
    variant_start: Optional[int] = None
    variant_end: Optional[int] = None


# =============================================================================
# Edit Distance Functions
# =============================================================================

def edit_distance(s1: str, s2: str) -> int:
    """Calculate Levenshtein edit distance between two strings."""
    if HAS_LEVENSHTEIN:
        return Levenshtein.distance(s1, s2)
    else:
        if len(s1) < len(s2):
            return edit_distance(s2, s1)
        if len(s2) == 0:
            return len(s1)
        previous_row = range(len(s2) + 1)
        for i, c1 in enumerate(s1):
            current_row = [i + 1]
            for j, c2 in enumerate(s2):
                insertions = previous_row[j + 1] + 1
                deletions = current_row[j] + 1
                substitutions = previous_row[j] + (c1 != c2)
                current_row.append(min(insertions, deletions, substitutions))
            previous_row = current_row
        return previous_row[-1]


# =============================================================================
# Oligo Design Loading
# =============================================================================

def load_oligo_designs(filepath: Path, library_name: str = "") -> Tuple[pd.DataFrame, Dict[str, OligoInfo]]:
    """Load oligo design file and build lookup structures."""
    print(f"Loading oligo designs from {filepath.name} ({library_name})...")
    df = pd.read_csv(filepath, sep="\t")

    df["guide"] = df["oligo_seq"].str[OLIGO_GUIDE_START:OLIGO_GUIDE_END]
    df["donor"] = df["oligo_seq"].str[OLIGO_DONOR_START:OLIGO_DONOR_END].str.upper()
    df["middle_donor"] = df["donor"].str[MIDDLE_DONOR_START:MIDDLE_DONOR_END]

    oligo_lookup = {}
    for _, row in df.iterrows():
        info = OligoInfo(
            oligo_name=row["oligo_name"],
            guide=row["guide"],
            donor=row["donor"],
            middle_donor=row["middle_donor"],
            oligo_seq=row["oligo_seq"],
        )
        info = parse_oligo_name(info)
        oligo_lookup[row["oligo_name"]] = info

    print(f"  Loaded {len(oligo_lookup):,} oligo designs")
    return df, oligo_lookup


def parse_oligo_name(info: OligoInfo) -> OligoInfo:
    """Parse structured oligo name to extract annotation fields."""
    name = info.oligo_name
    parts = name.split("_")

    try:
        if len(parts) >= 3:
            info.systematic_ORF_name = parts[1] if parts[1].startswith("Y") else None
            info.gene_name = parts[2] if not parts[2].isdigit() else None

        for i, part in enumerate(parts):
            if part.startswith("chr"):
                info.chrom = part
                if i + 1 < len(parts) and parts[i + 1] in ["+", "-"]:
                    info.strand = parts[i + 1]
                if i + 2 < len(parts) and parts[i + 2].isdigit():
                    info.pam_pos = int(parts[i + 2])
                break

        if len(parts) >= 2:
            try:
                info.variant_end = int(parts[-1])
                info.variant_start = int(parts[-2])
            except ValueError:
                pass

    except Exception:
        pass

    return info


def load_individual_variants_lookup(filepath: Path) -> Dict[str, str]:
    """Load harmonized oligo key and create oligo_name -> individual_variants lookup."""
    if not filepath or not filepath.exists():
        print(f"  WARNING: Harmonized oligo key not found: {filepath}")
        print("    Delta variant computation will be skipped.")
        return {}

    print(f"Loading harmonized oligo key from {filepath.name}...")
    df = pd.read_csv(filepath, sep="\t", usecols=["oligo_name", "individual_variants"])

    lookup = {}
    for _, row in df.iterrows():
        oligo_name = row["oligo_name"]
        indiv_vars = row["individual_variants"]
        if pd.notna(oligo_name) and pd.notna(indiv_vars):
            lookup[oligo_name] = str(indiv_vars)

    print(f"  Loaded {len(lookup):,} oligo -> individual_variants mappings")
    return lookup


def get_oligo_lookup_for_yeast_sublibrary(
    yeast_sublibrary: str,
    spg_lookup: Dict[str, OligoInfo],
    spcas9_lookup: Dict[str, OligoInfo]
) -> Dict[str, OligoInfo]:
    """Get the appropriate oligo lookup dict for a given yeast_sublibrary."""
    design = _resolve_design(yeast_sublibrary)
    if design == "spg":
        return spg_lookup
    else:
        return spcas9_lookup


def get_kmer_index_for_yeast_sublibrary(
    yeast_sublibrary: str,
    spg_kmer_index: Dict[str, set],
    spcas9_kmer_index: Dict[str, set]
) -> Dict[str, set]:
    """Get the appropriate kmer index for a given yeast_sublibrary."""
    design = _resolve_design(yeast_sublibrary)
    if design == "spg":
        return spg_kmer_index
    else:
        return spcas9_kmer_index


# =============================================================================
# K-mer Index for Fast Candidate Finding
# =============================================================================

KMER_SIZE = 12  # 12-mers for indexing
TOP_KMER_CANDIDATES = 30  # Number of top k-mer matches to evaluate

def build_kmer_index(oligo_lookup: Dict[str, OligoInfo]) -> Dict[str, set]:
    """
    Build inverted index: k-mer -> set of oligo names containing it.

    This enables O(1) lookup of candidate oligos that share k-mers with
    an observed donor, reducing candidate search from O(n) to O(k) where
    k is the number of unique k-mers in the observed donor.
    """
    kmer_to_oligos = {}
    for oligo_name, info in oligo_lookup.items():
        donor = info.donor.upper()
        for i in range(len(donor) - KMER_SIZE + 1):
            kmer = donor[i:i + KMER_SIZE]
            if kmer not in kmer_to_oligos:
                kmer_to_oligos[kmer] = set()
            kmer_to_oligos[kmer].add(oligo_name)
    return kmer_to_oligos


def find_closest_designed_donor(
    observed_donor: str,
    oligo_lookup: Dict[str, OligoInfo],
    kmer_index: Dict[str, set] = None,
    max_edit_distance: int = MAX_EDIT_DISTANCE
) -> Tuple[Optional[str], int, int, Optional[str]]:
    """
    Find designed donor with minimum edit distance to observed donor.

    Uses k-mer pre-filtering for ~10,000x speedup:
    1. Extract all k-mers from observed donor
    2. Score oligos by number of shared k-mers
    3. Only compute expensive edit distance for top candidates

    Returns:
        (closest_oligo_name, full_edit_distance, middle_edit_distance, inferred_guide)
    """
    if not observed_donor or len(observed_donor) < MIDDLE_DONOR_END:
        return None, -1, -1, None

    observed_donor = observed_donor.upper()
    observed_middle = observed_donor[MIDDLE_DONOR_START:MIDDLE_DONOR_END]

    best_oligo = None
    best_full_ed = float('inf')
    best_middle_ed = float('inf')

    # Use k-mer index if available (fast path)
    if kmer_index:
        # Extract k-mers from observed donor
        observed_kmers = set()
        for i in range(len(observed_donor) - KMER_SIZE + 1):
            observed_kmers.add(observed_donor[i:i + KMER_SIZE])

        # Count k-mer matches per oligo
        oligo_scores = {}
        for kmer in observed_kmers:
            if kmer in kmer_index:
                for oligo_name in kmer_index[kmer]:
                    oligo_scores[oligo_name] = oligo_scores.get(oligo_name, 0) + 1

        # Get top candidates by k-mer overlap
        if oligo_scores:
            candidates = sorted(oligo_scores.items(), key=lambda x: -x[1])[:TOP_KMER_CANDIDATES]
        else:
            candidates = []

        # Compute edit distance only for top k-mer candidates
        for oligo_name, kmer_count in candidates:
            info = oligo_lookup[oligo_name]
            full_ed = edit_distance(observed_donor, info.donor)
            if full_ed < best_full_ed:
                best_full_ed = full_ed
                best_middle_ed = edit_distance(observed_middle, info.middle_donor)
                best_oligo = oligo_name
    else:
        # Fallback: old method using middle donor screening (slower)
        candidates = []
        for oligo_name, info in oligo_lookup.items():
            middle_ed = edit_distance(observed_middle, info.middle_donor)
            if middle_ed <= max_edit_distance:
                candidates.append((oligo_name, middle_ed))

        candidates.sort(key=lambda x: x[1])

        for oligo_name, middle_ed in candidates[:100]:
            info = oligo_lookup[oligo_name]
            full_ed = edit_distance(observed_donor, info.donor)
            if full_ed < best_full_ed:
                best_full_ed = full_ed
                best_middle_ed = middle_ed
                best_oligo = oligo_name

    if best_full_ed == float('inf') or best_full_ed > max_edit_distance:
        return None, -1, -1, None

    inferred_guide = oligo_lookup[best_oligo].guide if best_oligo else None
    return best_oligo, int(best_full_ed), int(best_middle_ed), inferred_guide


# =============================================================================
# Alignment Functions
# =============================================================================

def write_donors_to_fasta(donors: Dict[str, str], output_path: Path) -> None:
    """Write donor sequences to FASTA file."""
    with open(output_path, "w") as f:
        for bc1, donor_seq in donors.items():
            f.write(f">{bc1}\n{donor_seq}\n")


def check_aligner_available(aligner: str) -> bool:
    """Check if an aligner is available in PATH."""
    try:
        if aligner == "bwa":
            result = subprocess.run(["bwa"], capture_output=True)
            return True
        elif aligner == "bbmap":
            result = subprocess.run(["bbmap.sh", "--version"], capture_output=True)
            return True
        elif aligner == "minimap2":
            result = subprocess.run(["minimap2", "--version"], capture_output=True)
            return result.returncode == 0
    except FileNotFoundError:
        return False
    return False


def run_bwa(donors_fasta: Path, reference_fasta: Path, output_sam: Path, threads: int = 1) -> bool:
    """Run bwa mem alignment."""
    try:
        index_file = Path(str(reference_fasta) + ".bwt")
        if not index_file.exists():
            print("    Creating bwa index...")
            subprocess.run(["bwa", "index", str(reference_fasta)], check=True, capture_output=True)

        with open(output_sam, "w") as outf:
            subprocess.run(
                ["bwa", "mem", "-t", str(threads), str(reference_fasta), str(donors_fasta)],
                stdout=outf, stderr=subprocess.DEVNULL, check=True,
            )
        return True
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"    bwa failed: {e}")
        return False


def run_bbmap(donors_fasta: Path, reference_fasta: Path, output_sam: Path, threads: int = 1) -> bool:
    """Run bbmap alignment."""
    try:
        subprocess.run([
            "bbmap.sh",
            f"in={donors_fasta}",
            f"out={output_sam}",
            f"ref={reference_fasta}",
            "nmtag=t", "mdtag=t", "local=f",
            f"threads={threads}",
        ], check=True, capture_output=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"    bbmap failed: {e}")
        return False


def run_minimap2(donors_fasta: Path, reference_fasta: Path, output_sam: Path, threads: int = 1) -> bool:
    """Run minimap2 alignment."""
    try:
        with open(output_sam, "w") as outf:
            subprocess.run(
                ["minimap2", "-a", "-x", "sr", "-t", str(threads), str(reference_fasta), str(donors_fasta)],
                stdout=outf, stderr=subprocess.DEVNULL, check=True,
            )
        return True
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"    minimap2 failed: {e}")
        return False


def score_alignment(read) -> int:
    """Score an alignment based on total bp in mismatches + indels + soft-clips."""
    if read.is_unmapped:
        return float('inf')

    score = 0
    cigar_ops = read.cigartuples or []

    for op, length in cigar_ops:
        if op == 4:  # S = soft-clip
            score += length
        elif op == 1:  # I = insertion
            score += length
        elif op == 2:  # D = deletion
            score += length
        elif op == 8:  # X = mismatch
            score += length

    if read.has_tag('NM'):
        nm = read.get_tag('NM')
        indel_bp = sum(length for op, length in cigar_ops if op in [1, 2])
        x_bp = sum(length for op, length in cigar_ops if op == 8)
        if x_bp == 0:
            score += nm - indel_bp

    return score


def run_triple_aligners(
    donors_fasta: Path,
    reference_fasta: Path,
    output_dir: Path,
    threads_per_aligner: int = 1,
) -> Dict[str, Tuple]:
    """
    Run all three aligners and select best alignment for each donor.

    Returns:
        Dict[bc1] = (best_alignment, aligner_name, score)
    """
    if not HAS_PYSAM:
        print("Error: pysam not available for alignment parsing")
        return {}

    output_dir.mkdir(parents=True, exist_ok=True)

    aligners = [
        ("bwa", run_bwa, output_dir / "bwa.sam"),
        ("bbmap", run_bbmap, output_dir / "bbmap.sam"),
        ("minimap2", run_minimap2, output_dir / "minimap2.sam"),
    ]

    successful_aligners = []
    for name, func, output_sam in aligners:
        if check_aligner_available(name.split(".")[0]):
            if func(donors_fasta, reference_fasta, output_sam, threads=threads_per_aligner):
                successful_aligners.append((name, output_sam))

    if not successful_aligners:
        print("Error: No aligners succeeded")
        return {}

    best_alignments = {}

    for aligner_name, sam_path in successful_aligners:
        try:
            samfile = pysam.AlignmentFile(str(sam_path), "r")
            for read in samfile:
                bc1 = read.query_name
                score = score_alignment(read)

                if bc1 not in best_alignments:
                    best_alignments[bc1] = (read, aligner_name, score)
                else:
                    current_score = best_alignments[bc1][2]
                    current_aligner = best_alignments[bc1][1]

                    if score < current_score:
                        best_alignments[bc1] = (read, aligner_name, score)
                    elif score == current_score:
                        priority = {"bwa": 0, "bbmap": 1, "minimap2": 2}
                        if priority.get(aligner_name, 3) < priority.get(current_aligner, 3):
                            best_alignments[bc1] = (read, aligner_name, score)

            samfile.close()
        except Exception as e:
            print(f"    Error parsing {sam_path}: {e}")

    return best_alignments


# =============================================================================
# Variant Extraction
# =============================================================================

def extract_variants_from_alignment(
    read,
    reference_fasta_path: Path
) -> List[Tuple[str, int, str, str]]:
    """Extract variants from alignment as (chrom, pos, ref, alt) tuples."""
    if not HAS_PYSAM:
        return []

    if read.is_unmapped:
        return []

    variants = []
    ref_fasta = pysam.FastaFile(str(reference_fasta_path))
    chrom = read.reference_name

    try:
        pairs = read.get_aligned_pairs(with_seq=True)
    except Exception:
        return []

    i = 0
    while i < len(pairs):
        query_pos, ref_pos, ref_base = pairs[i]

        if ref_pos is None:
            # Insertion
            ins_bases = ""
            while i < len(pairs) and pairs[i][1] is None:
                if pairs[i][0] is not None:
                    ins_bases += read.query_sequence[pairs[i][0]]
                i += 1
            if ins_bases:
                if ref_pos is None and i > 0:
                    prev_ref_pos = None
                    for j in range(i - 1, -1, -1):
                        if pairs[j][1] is not None:
                            prev_ref_pos = pairs[j][1]
                            break
                    if prev_ref_pos is not None:
                        ref_base = ref_fasta.fetch(chrom, prev_ref_pos, prev_ref_pos + 1)
                        variants.append((chrom, prev_ref_pos + 1, ref_base, ref_base + ins_bases))
            continue

        if query_pos is None:
            # Deletion
            del_start = ref_pos
            del_bases = ""
            while i < len(pairs) and pairs[i][0] is None:
                if pairs[i][2] is not None:
                    del_bases += pairs[i][2]
                i += 1
            if del_bases:
                if del_start > 0:
                    anchor = ref_fasta.fetch(chrom, del_start - 1, del_start)
                    variants.append((chrom, del_start, anchor + del_bases, anchor))
            continue

        # Mismatch
        if ref_base is not None and query_pos is not None:
            query_base = read.query_sequence[query_pos]
            if ref_base.upper() != query_base.upper():
                variants.append((chrom, ref_pos + 1, ref_base.upper(), query_base.upper()))

        i += 1

    ref_fasta.close()
    return variants


def format_variants_string(variants: List[Tuple[str, int, str, str]]) -> str:
    """Format variants as a comma-separated string."""
    if not variants:
        return ""
    return ",".join(f"{pos}_{ref}_{alt}" for _, pos, ref, alt in variants)


def parse_individual_variants_string(indiv_vars_str: str) -> List[Tuple[str, int, str, str]]:
    """Parse individual_variants string from harmonized oligo key."""
    if not indiv_vars_str or pd.isna(indiv_vars_str):
        return []

    variants = []
    for var in str(indiv_vars_str).split(","):
        var = var.strip()
        if not var:
            continue

        parts = var.split("_")
        if len(parts) < 3:
            continue

        try:
            if parts[0].startswith("chr"):
                chrom = parts[0]
                pos = int(parts[1])
                ref = parts[2]
                alt = parts[3] if len(parts) > 3 else ""
            else:
                chrom = ""
                pos = int(parts[0])
                ref = parts[1]
                alt = parts[2]

            if pos > 0 and ref and alt:
                variants.append((chrom, pos, ref, alt))

        except (ValueError, IndexError):
            continue

    return variants


def compute_delta_variants(
    observed_variants: List[Tuple[str, int, str, str]],
    designed_variants: List[Tuple[str, int, str, str]],
) -> List[Tuple[str, int, str, str]]:
    """Compute delta variants (observed but not designed)."""
    if not observed_variants:
        return []

    designed_positions = {pos for _, pos, _, _ in designed_variants}
    delta = [
        (chrom, pos, ref, alt)
        for chrom, pos, ref, alt in observed_variants
        if pos not in designed_positions
    ]

    return delta


# =============================================================================
# Extended Donor Detection
# =============================================================================

def detect_extended_donor(
    observed_donor: str,
    closest_designed_donor: Optional[str],
    expected_length: int = DONOR_LENGTH
) -> Dict:
    """Detect if a donor is extended beyond the expected 129bp."""
    result = {
        "is_extended": False,
        "extension_non_bc0_side_bp": 0,
        "extension_bc0_side_bp": 0,
        "observed_length": len(observed_donor) if observed_donor else 0,
    }

    if not observed_donor:
        return result

    observed_len = len(observed_donor)

    if observed_len <= expected_length:
        return result

    result["is_extended"] = True
    extra_bp = observed_len - expected_length

    if closest_designed_donor and len(closest_designed_donor) == expected_length:
        designed_middle = closest_designed_donor[MIDDLE_DONOR_START:MIDDLE_DONOR_END]
        middle_start_in_observed = observed_donor.find(designed_middle)

        if middle_start_in_observed >= 0:
            result["extension_non_bc0_side_bp"] = middle_start_in_observed - MIDDLE_DONOR_START
            result["extension_bc0_side_bp"] = extra_bp - result["extension_non_bc0_side_bp"]
        else:
            result["extension_non_bc0_side_bp"] = extra_bp // 2
            result["extension_bc0_side_bp"] = extra_bp - result["extension_non_bc0_side_bp"]
    else:
        result["extension_non_bc0_side_bp"] = extra_bp // 2
        result["extension_bc0_side_bp"] = extra_bp - result["extension_non_bc0_side_bp"]

    result["extension_non_bc0_side_bp"] = max(0, result["extension_non_bc0_side_bp"])
    result["extension_bc0_side_bp"] = max(0, result["extension_bc0_side_bp"])

    return result


# =============================================================================
# Rescue Assessment
# =============================================================================

def calculate_incorporation_probability(distance_bp: int) -> float:
    """Calculate incorporation probability based on distance from adjacent variant."""
    if distance_bp <= 2:
        return 1.0
    elif distance_bp <= 5:
        return 0.85
    else:
        return np.exp(-0.0693 * distance_bp)


def classify_rescue_status(
    observed_variants: List[Tuple[str, int, str, str]],
    designed_variants: List[Tuple[str, int, str, str]],
    oligo_info: OligoInfo
) -> Dict:
    """Classify whether this donor can be rescued for hit calling."""
    result = {
        "rescue_status": None,
        "closest_delta_to_designed": -1,
        "incorporation_probability": 0.0,
        "delta_variant_count": 0,
    }

    if not observed_variants:
        result["rescue_status"] = "NO_VARIANTS"
        return result

    designed_positions = {pos for _, pos, _, _ in designed_variants}
    delta_variants = [
        (chrom, pos, ref, alt) for chrom, pos, ref, alt in observed_variants
        if pos not in designed_positions
    ]
    result["delta_variant_count"] = len(delta_variants)

    if not delta_variants:
        result["rescue_status"] = "NO_DELTA"
        result["incorporation_probability"] = 0.0
        return result

    min_distance = float('inf')
    for _, delta_pos, _, _ in delta_variants:
        for _, designed_pos, _, _ in designed_variants:
            dist = abs(delta_pos - designed_pos)
            if dist < min_distance:
                min_distance = dist

    result["closest_delta_to_designed"] = min_distance if min_distance != float('inf') else -1

    if min_distance != float('inf'):
        result["incorporation_probability"] = calculate_incorporation_probability(int(min_distance))

    if min_distance > RESCUE_DISTANCE_THRESHOLD:
        result["rescue_status"] = "CAN_ATTRIBUTE"
    elif min_distance >= 15:
        result["rescue_status"] = "NEEDS_REVIEW"
    else:
        result["rescue_status"] = "CONFOUNDED"

    return result


def assign_confidence_tier(annotation: AnnotationResult) -> str:
    """Assign confidence tier based on annotation results."""
    ed = annotation.closest_donor_edit_distance
    rescue = annotation.rescue_status

    if ed < 0 or annotation.closest_oligo_name is None:
        return "tier_5_no_match"

    if ed <= 3:
        if rescue in ["CAN_ATTRIBUTE", "NO_DELTA", "NO_VARIANTS"]:
            return "tier_1_low_ed_attributable"
        elif rescue == "NEEDS_REVIEW":
            return "tier_2_low_ed_review"
        else:
            return "tier_2_low_ed_confounded"
    elif ed <= 10:
        if rescue in ["CAN_ATTRIBUTE", "NO_DELTA", "NO_VARIANTS"]:
            return "tier_3_mod_ed_attributable"
        else:
            return "tier_4_mod_ed_confounded"
    else:
        return "tier_5_high_ed"


# =============================================================================
# Chunk Management
# =============================================================================

def get_chunk_paths(n_chunks: int, chunks_dir: Path) -> List[Path]:
    """Get paths for chunk files."""
    return [chunks_dir / f"bc1_unmatched_chunk_{i:03d}.tsv" for i in range(n_chunks)]


def split_bc1_reference_table(
    input_file: Path,
    n_chunks: int,
    output_dir: Path,
    force: bool = False,
) -> Tuple[List[Path], int]:
    """Split bc1 reference table into chunks for parallel processing."""
    output_dir.mkdir(parents=True, exist_ok=True)
    chunk_paths = get_chunk_paths(n_chunks, output_dir)

    print(f"  Loading {input_file.name} for splitting...")
    start = time.time()

    df = pd.read_csv(input_file, sep="\t", low_memory=False)
    load_time = time.time() - start
    print(f"  Loaded {len(df):,} rows in {load_time:.1f}s")

    # Detect the appropriate column for matching status
    if "bdb_oligo_name" in df.columns:
        match_col = "bdb_oligo_name"
    elif "oligo_name" in df.columns:
        match_col = "oligo_name"
    else:
        print("  Warning: No oligo_name column found. Processing all rows.")
        unmatched = df
        match_col = None

    if match_col:
        unmatched = df[df[match_col].isna() | (df[match_col] == "")]

    total_unmatched = len(unmatched)
    print(f"  Unmatched donors: {total_unmatched:,} ({100*total_unmatched/len(df):.1f}%)")

    if total_unmatched == 0:
        print("  No unmatched donors to process")
        return [], 0

    chunk_size = (total_unmatched + n_chunks - 1) // n_chunks
    print(f"  Splitting into {n_chunks} chunks of ~{chunk_size:,} rows each...")

    for i, chunk_path in enumerate(chunk_paths):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, total_unmatched)

        if start_idx >= total_unmatched:
            chunk_df = unmatched.iloc[:0]
        else:
            chunk_df = unmatched.iloc[start_idx:end_idx]

        chunk_df.to_csv(chunk_path, sep="\t", index=False)
        print(f"    Chunk {i}: {len(chunk_df):,} rows -> {chunk_path.name}")

    return chunk_paths, total_unmatched


# =============================================================================
# Worker Functions
# =============================================================================

def init_worker(
    spg_lookup: Dict[str, OligoInfo],
    spcas9_lookup: Dict[str, OligoInfo],
    spg_kmer_index: Dict[str, set],
    spcas9_kmer_index: Dict[str, set],
    individual_variants_lookup: Dict[str, str],
    reference_fasta: Path,
    cds_lookup: Dict[str, List[CdsInfo]],
    output_dir: Path,
    progress_dir: Path,
):
    """Initialize worker process with shared data."""
    global _WORKER_OLIGO_LOOKUP_SPG, _WORKER_OLIGO_LOOKUP_SPCAS9
    global _WORKER_KMER_INDEX_SPG, _WORKER_KMER_INDEX_SPCAS9
    global _WORKER_INDIVIDUAL_VARIANTS_LOOKUP
    global _WORKER_REFERENCE_FASTA, _WORKER_CDS_LOOKUP
    global _WORKER_OUTPUT_DIR, _WORKER_PROGRESS_DIR

    _WORKER_OLIGO_LOOKUP_SPG = spg_lookup
    _WORKER_OLIGO_LOOKUP_SPCAS9 = spcas9_lookup
    _WORKER_KMER_INDEX_SPG = spg_kmer_index
    _WORKER_KMER_INDEX_SPCAS9 = spcas9_kmer_index
    _WORKER_INDIVIDUAL_VARIANTS_LOOKUP = individual_variants_lookup
    _WORKER_REFERENCE_FASTA = reference_fasta
    _WORKER_CDS_LOOKUP = cds_lookup
    _WORKER_OUTPUT_DIR = output_dir
    _WORKER_PROGRESS_DIR = progress_dir


def log_progress(chunk_idx: int, message: str, progress_dir: Path = None):
    """Write progress message to chunk-specific log file."""
    if progress_dir is None:
        progress_dir = _WORKER_PROGRESS_DIR
    if progress_dir is None:
        return

    progress_dir.mkdir(parents=True, exist_ok=True)
    log_file = progress_dir / f"chunk_{chunk_idx:03d}.log"
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(log_file, "a") as f:
        f.write(f"[{timestamp}] {message}\n")


def process_chunk(chunk_path: Path, chunk_idx: int) -> Tuple[int, Path]:
    """Process a single chunk of unmatched donors."""
    spg_lookup = _WORKER_OLIGO_LOOKUP_SPG
    spcas9_lookup = _WORKER_OLIGO_LOOKUP_SPCAS9
    spg_kmer_index = _WORKER_KMER_INDEX_SPG
    spcas9_kmer_index = _WORKER_KMER_INDEX_SPCAS9
    individual_variants_lookup = _WORKER_INDIVIDUAL_VARIANTS_LOOKUP
    reference_fasta = _WORKER_REFERENCE_FASTA
    cds_lookup = _WORKER_CDS_LOOKUP
    output_dir = _WORKER_OUTPUT_DIR
    progress_dir = _WORKER_PROGRESS_DIR

    chunk_output_dir = output_dir / f"chunk_{chunk_idx:03d}"
    chunk_output_dir.mkdir(parents=True, exist_ok=True)

    log_progress(chunk_idx, f"Starting chunk processing: {chunk_path.name}", progress_dir)

    chunk_df = pd.read_csv(chunk_path, sep="\t")
    n_rows = len(chunk_df)
    log_progress(chunk_idx, f"Loaded {n_rows:,} rows", progress_dir)

    if n_rows == 0:
        result_path = chunk_output_dir / "results.tsv"
        pd.DataFrame().to_csv(result_path, sep="\t", index=False)
        log_progress(chunk_idx, "Empty chunk, skipping", progress_dir)
        return chunk_idx, result_path

    results = []
    donors_for_alignment = {}

    # Detect donor column
    donor_col = "bdb_donor" if "bdb_donor" in chunk_df.columns else "donor"
    # Check for yeast_sublibrary column (or legacy "sublibrary" for backward compatibility)
    yeast_sublibrary_col = None
    if "yeast_sublibrary" in chunk_df.columns:
        yeast_sublibrary_col = "yeast_sublibrary"
    elif "sublibrary" in chunk_df.columns:
        yeast_sublibrary_col = "sublibrary"
    # V_library / library_ID column: HANDOFF_20260527 fallback when yeast_sublibrary is NaN
    v_library_col = None
    for _col in ("V_library", "library_ID"):
        if _col in chunk_df.columns:
            v_library_col = _col
            break

    log_progress(chunk_idx, "Step 1: Reading existing gdb columns and finding closest donors...", progress_dir)

    for row_idx, (idx, row) in enumerate(chunk_df.iterrows()):
        bc1 = row["bc1"]
        donor = row.get(donor_col, "")
        yeast_sublibrary = row.get(yeast_sublibrary_col, None) if yeast_sublibrary_col else None
        # Fall back to V_library when yeast_sublibrary is NaN/empty (HANDOFF_20260527)
        if not yeast_sublibrary or pd.isna(yeast_sublibrary):
            v_lib = row.get(v_library_col, None) if v_library_col else None
            if v_lib and not pd.isna(v_lib) and v_lib in V_LIBRARY_TO_DESIGN:
                yeast_sublibrary = v_lib  # resolved via V_LIBRARY_TO_DESIGN in _resolve_design
            else:
                yeast_sublibrary = "SpG_NGNG"

        result = AnnotationResult(bc1=bc1, yeast_sublibrary=yeast_sublibrary)
        oligo_lookup = get_oligo_lookup_for_yeast_sublibrary(yeast_sublibrary, spg_lookup, spcas9_lookup)
        kmer_index = get_kmer_index_for_yeast_sublibrary(yeast_sublibrary, spg_kmer_index, spcas9_kmer_index)

        # Use existing gdb columns
        guide_val = row.get("guide", "")
        result.gdb_lookup_found = pd.notna(guide_val) and guide_val != ""
        result.gdb_lookup_guide = guide_val if pd.notna(guide_val) else None
        gdb_oligo = row.get("gdb_oligo_name", "")
        result.gdb_lookup_oligo_name = gdb_oligo if pd.notna(gdb_oligo) and gdb_oligo != "" else None
        result.gdb_lookup_perfect_donor = bool(row.get("gdb_perfect_donor", False))
        result.gdb_lookup_perfect_middle = bool(row.get("gdb_perfect_middle_donor_region", False))

        # Closest designed donor (using k-mer index for ~10,000x speedup)
        closest_designed_seq = None
        if donor and len(donor) >= MIDDLE_DONOR_END:
            closest_name, full_ed, middle_ed, inferred_guide = find_closest_designed_donor(
                donor, oligo_lookup, kmer_index=kmer_index
            )
            result.closest_oligo_name = closest_name
            result.closest_donor_edit_distance = full_ed
            result.closest_middle_donor_edit_distance = middle_ed
            result.inferred_guide = inferred_guide

            if closest_name and closest_name in oligo_lookup:
                closest_designed_seq = oligo_lookup[closest_name].donor

        # Extended donor detection
        if donor:
            ext_result = detect_extended_donor(donor, closest_designed_seq)
            result.is_extended = ext_result["is_extended"]
            result.extension_non_bc0_side_bp = ext_result["extension_non_bc0_side_bp"]
            result.extension_bc0_side_bp = ext_result["extension_bc0_side_bp"]
            result.observed_donor_length = ext_result["observed_length"]

        results.append(result)

        donor_to_align = result.gdb_lookup_donor if result.gdb_lookup_donor else donor
        if donor_to_align and len(donor_to_align) >= 50:
            donors_for_alignment[bc1] = donor_to_align

        if (row_idx + 1) % 1000 == 0:
            log_progress(chunk_idx, f"  Processed {row_idx + 1:,}/{n_rows:,} rows", progress_dir)

    log_progress(chunk_idx, f"Step 1 complete. {len(donors_for_alignment):,} donors for alignment.", progress_dir)

    # Step 2: Run alignment
    if donors_for_alignment and reference_fasta and reference_fasta.exists():
        log_progress(chunk_idx, "Step 2: Running triple-aligner...", progress_dir)
        donors_fasta = chunk_output_dir / "donors.fa"
        write_donors_to_fasta(donors_for_alignment, donors_fasta)

        best_alignments = run_triple_aligners(
            donors_fasta, reference_fasta, chunk_output_dir, threads_per_aligner=1
        )
        log_progress(chunk_idx, f"  Got {len(best_alignments):,} alignments", progress_dir)

        # Step 3: Extract variants and annotate
        log_progress(chunk_idx, "Step 3: Extracting variants and annotating...", progress_dir)
        bc1_to_result = {r.bc1: r for r in results}

        for bc1, (alignment, aligner, score) in best_alignments.items():
            if bc1 not in bc1_to_result:
                continue

            result = bc1_to_result[bc1]
            oligo_lookup = get_oligo_lookup_for_yeast_sublibrary(
                result.yeast_sublibrary, spg_lookup, spcas9_lookup
            )

            if not alignment.is_unmapped:
                result.alignment_chrom = alignment.reference_name
                result.alignment_start = alignment.reference_start
                result.alignment_end = alignment.reference_end
                result.alignment_mapq = alignment.mapping_quality
                result.alignment_score = score
                result.aligner_used = aligner

                observed_variants = extract_variants_from_alignment(alignment, reference_fasta)
                result.observed_individual_variants = format_variants_string(observed_variants)

                if observed_variants and cds_lookup:
                    var_type, sys_gene, com_gene = annotate_variant_consequence(
                        observed_variants, cds_lookup
                    )
                    result.codon_variant_type = var_type
                    if not result.systematic_ORF_name:
                        result.systematic_ORF_name = sys_gene
                    if not result.gene_name:
                        result.gene_name = com_gene

                if result.closest_oligo_name and result.closest_oligo_name in oligo_lookup:
                    closest_info = oligo_lookup[result.closest_oligo_name]
                    if closest_info.systematic_ORF_name:
                        result.systematic_ORF_name = closest_info.systematic_ORF_name
                    if closest_info.gene_name:
                        result.gene_name = closest_info.gene_name

                    designed_indiv_vars_str = individual_variants_lookup.get(result.closest_oligo_name, "")
                    result.designed_individual_variants = designed_indiv_vars_str if designed_indiv_vars_str else None

                    designed_variants = parse_individual_variants_string(designed_indiv_vars_str)
                    delta_variants = compute_delta_variants(observed_variants, designed_variants)
                    result.delta_variants = format_variants_string(delta_variants) if delta_variants else None

                    rescue_result = classify_rescue_status(
                        observed_variants, designed_variants, closest_info
                    )
                    result.rescue_status = rescue_result["rescue_status"]
                    result.closest_delta_to_designed = rescue_result["closest_delta_to_designed"]
                    result.incorporation_probability = rescue_result["incorporation_probability"]
                    result.delta_variant_count = rescue_result["delta_variant_count"]

            result.confidence_tier = assign_confidence_tier(result)

    log_progress(chunk_idx, "Step 4: Converting results to DataFrame...", progress_dir)

    # Convert results to DataFrame
    result_records = []
    for result in results:
        record = {
            "bc1": result.bc1,
            "yeast_sublibrary": result.yeast_sublibrary,
            "gdb_lookup_found": result.gdb_lookup_found,
            "gdb_lookup_guide": result.gdb_lookup_guide,
            "gdb_lookup_oligo_name": result.gdb_lookup_oligo_name,
            "gdb_lookup_perfect_donor": result.gdb_lookup_perfect_donor,
            "gdb_lookup_perfect_middle": result.gdb_lookup_perfect_middle,
            "closest_oligo_name": result.closest_oligo_name,
            "closest_donor_edit_distance": result.closest_donor_edit_distance,
            "closest_middle_donor_edit_distance": result.closest_middle_donor_edit_distance,
            "inferred_guide": result.inferred_guide,
            "alignment_chrom": result.alignment_chrom,
            "alignment_start": result.alignment_start,
            "alignment_end": result.alignment_end,
            "alignment_mapq": result.alignment_mapq,
            "alignment_score": result.alignment_score,
            "aligner_used": result.aligner_used,
            "observed_individual_variants": result.observed_individual_variants,
            "designed_individual_variants": result.designed_individual_variants,
            "delta_variants": result.delta_variants,
            "delta_variant_count": result.delta_variant_count,
            "systematic_ORF_name": result.systematic_ORF_name,
            "gene_name": result.gene_name,
            "aa_change": result.aa_change,
            "codon_variant_type": result.codon_variant_type,
            "is_chimera": result.is_chimera,
            "chimera_left_parent": result.chimera_left_parent,
            "chimera_right_parent": result.chimera_right_parent,
            "chimera_junction_pos": result.chimera_junction_pos,
            "chimera_gap_bp": result.chimera_gap_bp,
            "is_extended": result.is_extended,
            "extension_non_bc0_side_bp": result.extension_non_bc0_side_bp,
            "extension_bc0_side_bp": result.extension_bc0_side_bp,
            "observed_donor_length": result.observed_donor_length,
            "rescue_status": result.rescue_status,
            "closest_delta_to_designed": result.closest_delta_to_designed,
            "incorporation_probability": result.incorporation_probability,
            "confidence_tier": result.confidence_tier,
        }
        result_records.append(record)

    result_df = pd.DataFrame(result_records)
    result_path = chunk_output_dir / "results.tsv"
    result_df.to_csv(result_path, sep="\t", index=False)

    log_progress(chunk_idx, f"COMPLETE. Saved {len(result_df):,} results to {result_path.name}", progress_dir)

    return chunk_idx, result_path


def process_chunk_wrapper(args):
    """Wrapper for multiprocessing - unpacks arguments."""
    return process_chunk(*args)


# =============================================================================
# Main Pipeline
# =============================================================================

def print_summary(result_df: pd.DataFrame) -> None:
    """Print summary statistics of annotation results."""
    print("\n" + "=" * 70)
    print("ANNOTATION STATISTICS")
    print("=" * 70)

    total = len(result_df)
    print(f"\nTotal unmatched donors processed: {total:,}")

    # Edit distance distribution
    has_ed = result_df["closest_donor_edit_distance"] >= 0
    if has_ed.sum() > 0:
        print(f"\nEdit distance to closest designed donor (full 129bp):")
        ed = result_df.loc[has_ed, "closest_donor_edit_distance"]
        print(f"  Mean: {ed.mean():.1f}")
        print(f"  Median: {ed.median():.1f}")
        print(f"  ED ≤ 5: {(ed <= 5).sum():,} ({100*(ed <= 5).sum()/has_ed.sum():.1f}%)")

    # Middle donor edit distance
    has_middle_ed = result_df["closest_middle_donor_edit_distance"] >= 0
    if has_middle_ed.sum() > 0:
        print(f"\nMiddle donor edit distance (core 60bp, positions 30-90):")
        mid_ed = result_df.loc[has_middle_ed, "closest_middle_donor_edit_distance"]
        n_perfect = (mid_ed == 0).sum()
        print(f"  Middle ED = 0 (ATTRIBUTABLE): {n_perfect:,} ({100*n_perfect/has_middle_ed.sum():.1f}%)")

    # Alignment success
    aligned = result_df["alignment_chrom"].notna().sum()
    print(f"\nAlignment results:")
    print(f"  Successfully aligned: {aligned:,} ({100*aligned/total:.1f}%)")

    # Confidence tiers
    print(f"\nConfidence tiers:")
    for tier in result_df["confidence_tier"].value_counts().index:
        count = (result_df["confidence_tier"] == tier).sum()
        print(f"  {tier}: {count:,} ({100*count/total:.1f}%)")

    print("=" * 70)


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Annotate unmatched donors with variant information via genome alignment"
    )

    # Required arguments
    parser.add_argument("--input", "-i", required=True, type=Path,
                        help="Input bc1 reference table (TSV)")
    parser.add_argument("--output-dir", "-o", required=True, type=Path,
                        help="Output directory for results")

    # Oligo design files
    parser.add_argument("--spg-oligo-design", required=True, type=Path,
                        help="SpG oligo design file (2024 Twist pool)")
    parser.add_argument("--spcas9-oligo-design", required=True, type=Path,
                        help="SpCas9 oligo design file (2021 Twist pool)")

    # Reference files
    parser.add_argument("--reference-genome", required=True, type=Path,
                        help="Reference genome FASTA")
    parser.add_argument("--annotation-gff", required=True, type=Path,
                        help="GFF annotation file")

    # Optional files
    parser.add_argument("--harmonized-oligo-key", type=Path, default=None,
                        help="Harmonized oligo key file for delta variant computation")

    # Processing options
    parser.add_argument("--n-workers", type=int, default=DEFAULT_N_WORKERS,
                        help=f"Number of parallel workers (default: {DEFAULT_N_WORKERS})")
    parser.add_argument("--n-chunks", type=int, default=DEFAULT_N_CHUNKS,
                        help=f"Number of chunks to split input into (default: {DEFAULT_N_CHUNKS})")
    parser.add_argument("--force-split", action="store_true",
                        help="Force re-splitting even if chunk files exist")
    parser.add_argument("--n-samples", type=int, default=None,
                        help="Process only N samples (for testing)")
    parser.add_argument("--sequential", action="store_true",
                        help="Use sequential (single-process) mode")

    args = parser.parse_args()

    # Validate inputs
    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)

    # Create output directories
    args.output_dir.mkdir(parents=True, exist_ok=True)
    chunks_dir = args.output_dir / "chunks"
    progress_dir = args.output_dir / "progress"
    chunks_dir.mkdir(exist_ok=True)
    progress_dir.mkdir(exist_ok=True)

    start_time = time.time()

    print("=" * 80)
    print("UNMATCHED DONOR ANNOTATION PIPELINE")
    print("=" * 80)
    print(f"\nConfiguration:")
    print(f"  Input: {args.input}")
    print(f"  Output: {args.output_dir}")
    print(f"  Workers: {args.n_workers}")
    print(f"  Chunks: {args.n_chunks}")

    # ==========================================================================
    # Step 1: Load shared data
    # ==========================================================================
    print("\n" + "=" * 80)
    print("LOADING SHARED DATA")
    print("=" * 80)

    # Oligo designs
    print("\n1. Loading oligo designs...")
    _, spg_lookup = load_oligo_designs(args.spg_oligo_design, "SpG")
    _, spcas9_lookup = load_oligo_designs(args.spcas9_oligo_design, "SpCas9")
    print(f"  SpG oligos (NGNG/NGNH): {len(spg_lookup):,}")
    print(f"  SpCas9 oligos (NGG): {len(spcas9_lookup):,}")

    # Build k-mer indexes for fast donor matching
    print("\n1b. Building k-mer indexes for fast donor matching...")
    spg_kmer_index = build_kmer_index(spg_lookup)
    spcas9_kmer_index = build_kmer_index(spcas9_lookup)
    print(f"  SpG k-mer index: {len(spg_kmer_index):,} unique {KMER_SIZE}-mers")
    print(f"  SpCas9 k-mer index: {len(spcas9_kmer_index):,} unique {KMER_SIZE}-mers")

    # Individual variants lookup
    print("\n2. Loading harmonized oligo key...")
    individual_variants_lookup = {}
    if args.harmonized_oligo_key:
        individual_variants_lookup = load_individual_variants_lookup(args.harmonized_oligo_key)

    # GFF annotation
    print("\n3. Loading GFF annotation...")
    cds_lookup = load_gff_cds_info(args.annotation_gff)

    load_time = time.time() - start_time
    print(f"\n  Shared data loaded in {load_time:.1f}s")

    # ==========================================================================
    # Step 2: Split input data
    # ==========================================================================
    print("\n" + "=" * 80)
    print("SPLITTING INPUT DATA")
    print("=" * 80)

    chunk_paths, total_unmatched = split_bc1_reference_table(
        args.input, args.n_chunks, chunks_dir, force=args.force_split
    )

    if total_unmatched == 0:
        print("\nNo unmatched donors to process. Exiting.")
        return

    # ==========================================================================
    # Step 3: Process chunks
    # ==========================================================================
    print("\n" + "=" * 80)
    print("PROCESSING CHUNKS IN PARALLEL")
    print("=" * 80)

    # Clear old progress logs
    for log_file in progress_dir.glob("*.log"):
        log_file.unlink()

    chunk_args = [(chunk_path, i) for i, chunk_path in enumerate(chunk_paths)]

    print(f"\n  Processing {len(chunk_paths)} chunks with {args.n_workers} workers...")
    print(f"  Progress logs: {progress_dir}/chunk_*.log")
    process_start = time.time()

    with Pool(
        processes=args.n_workers,
        initializer=init_worker,
        initargs=(spg_lookup, spcas9_lookup, spg_kmer_index, spcas9_kmer_index,
                  individual_variants_lookup, args.reference_genome, cds_lookup,
                  args.output_dir, progress_dir)
    ) as pool:
        results = list(tqdm(
            pool.imap(process_chunk_wrapper, chunk_args),
            total=len(chunk_args),
            desc="  Chunks"
        ))

    process_time = time.time() - process_start
    print(f"\n  Chunk processing completed in {process_time:.1f}s")

    # ==========================================================================
    # Step 4: Merge results
    # ==========================================================================
    print("\n" + "=" * 80)
    print("MERGING RESULTS")
    print("=" * 80)

    result_dfs = []
    for chunk_idx, result_path in sorted(results):
        if result_path.exists():
            df = pd.read_csv(result_path, sep="\t")
            result_dfs.append(df)
            print(f"  Chunk {chunk_idx}: {len(df):,} rows")

    if result_dfs:
        merged_df = pd.concat(result_dfs, ignore_index=True)
        print(f"\n  Total merged: {len(merged_df):,} rows")
    else:
        merged_df = pd.DataFrame()

    # Save merged results
    output_file = args.output_dir / "unmatched_donor_annotations.tsv"
    merged_df.to_csv(output_file, sep="\t", index=False)
    print(f"\n  Saved to: {output_file}")

    # Print summary
    if len(merged_df) > 0:
        print_summary(merged_df)

    total_time = time.time() - start_time
    print(f"\nTotal time: {timedelta(seconds=int(total_time))}")


if __name__ == "__main__":
    main()
