"""
magestic.pipelines.redi — REDI (Restriction Enzyme Diagnostic Indexing) clone-array engine.

Operates on paired-end reads from REDI mating plates and yields:
- Per-(REDI_bc, MAGESTIC_bc1) read-count tables
- Colony purity scores (top-bc1 / total within top 4)
- Cross-array (SHA435 vs yT182/yT183) agreement metrics
- Annotation against the bc1-donor-bc0 reference table (from pipelines.bc1_donor_bc0)
- Type-(i) bc1 <-> edit correspondence (WGS) and type-(ii) REDI-bc <-> MAGESTIC-bc1
  correspondence (re-REDI) via a unified `core.correspondence` API.

Author: Kevin R. Roy
"""

__all__ = ["config", "core", "readers", "pipeline", "qc"]
