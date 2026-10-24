# Stage 3 functional and restart regressions

Run from the MAGESTIC repository with its WGS dependencies and the compatible local TRACE installed:

```sh
python -m pytest tests/test_pipelines/wgs/test_stage3_safety.py -q
```

The suite contains 23 parametrized cases and requires `pytest`, `pysam`, `PyYAML`, and `samtools` on `PATH`, in addition to the package dependencies. Each case generates its FASTA, BAM, index, manifest and tables in pytest temporary directories. There is no private-data or planning-workspace dependency. Concurrent same-sample calls use two bounded threads.

Tests exercise actual Stage 3 Half A/B entry points, externally specified junction geometry, donor-derived reference boundaries, explicit unmeasured status without eligible contigs, cleaned-BAM membership, BAI/CSI companion recovery, changed assignments, malformed inputs, completion publication, and duplicate computation suppression. Related TRACE primitives and standalone filtering tests live in TRACE's `tests/wgs/test_sv_plasmid_safety.py`.

Synthetic contract passes are not biological accuracy estimates. The suite does not establish network-filesystem locking behavior. Empty persistent `.lock` files prevent lock-inode races and must not be treated as completion markers; cache metadata is written after complete output publication.
