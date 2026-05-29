"""magestic.pipelines.redi.readers — IO wrappers for REDI inputs."""

from .fastq import extract_bc1_from_seq, extract_redi_bc_vectorized
from .keyfiles import combine_keyfiles, load_redi_bc_keys
from .plate_db import load_plate_db, resolve_plate, trace_lineage
from .wgs_outcome import load_wgs_outcome, load_wgs_plate_key

__all__ = [
    "extract_bc1_from_seq",
    "extract_redi_bc_vectorized",
    "combine_keyfiles",
    "load_redi_bc_keys",
    "load_plate_db",
    "resolve_plate",
    "trace_lineage",
    "load_wgs_outcome",
    "load_wgs_plate_key",
]
