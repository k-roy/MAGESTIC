"""magestic.pipelines.wgs.readers — input loaders and the REDI output seam.

`bam_inventory` enumerates BAMs from the manifest (never a glob) and reads the `@PG`
alignment provenance the D1 per-sample reference cross-check consumes.
`oligo_tables` is the Stage-2/Stage-3 designed-oligo seam; both of its public functions
raise today, naming the gate, because the inputs do not exist for this cohort.
`wgs_calls` is the REDI seam: one loader normalizing a TRACE emit table or the legacy
wide table into one schema that carries `has_call` and never pads a barcode.
"""

from .bam_inventory import (
    BamInventoryError,
    index_path,
    observed_references,
    reference_contig_lengths,
    require_indexed,
)
from .oligo_tables import (
    MissingOligoInputError,
    donor_samplesheet_rows,
    load_donor_sequences,
    read_oligo_table,
)
from .wgs_calls import CALL_COLUMNS, LEGACY_BC1_COLUMN, WgsCallsError, load_wgs_bc1_calls

__all__ = [
    "CALL_COLUMNS",
    "LEGACY_BC1_COLUMN",
    "BamInventoryError",
    "MissingOligoInputError",
    "WgsCallsError",
    "donor_samplesheet_rows",
    "index_path",
    "load_donor_sequences",
    "load_wgs_bc1_calls",
    "observed_references",
    "read_oligo_table",
    "reference_contig_lengths",
    "require_indexed",
]
