# REDI pipeline reference

The reusable engine is `src/magestic/pipelines/redi/`. The Python/umbrella CLI runs B through H; A processes one named FASTQ pair. This document describes the repaired engine contract. A completed test or run needs its numerical output and provenance, not just exit status. Validation results belong in the audit ledger; documentation is not an accuracy result.

## Configuration and completion

Configuration binds `project_root` to the YAML directory and other data paths to that root. Key filenames are relative to `paths.keyfiles`. Unknown/duplicate settings and invalid values fail before assay processing. Explicit assay selections and dated physical registrations are inputs, not deductions from barcode agreement.

Stages publish through the REDI artifact writer. It records input hashes before reading, checks inputs and source again before publication, atomically replaces individual outputs, and writes one completion sidecar last. A set of files is not a filesystem transaction: consumers must use the primary output's `.provenance.json` as completion evidence. That sidecar includes hashes for every output, including Stage H PNGs. Stage A can reuse only its complete matching five-output attestation when explicitly called with `overwrite=False`; other stage invocations recompute their outputs.

## A — Trim, merge and count

`magestic-redi-trim --config screen.yaml --sample <dated-prefix>` validates paired FASTQ structure, mate identity and pair counts. It uses the installed BBTools Java implementation with explicit heap/thread limits: BBDuk adapter/end trimming and BBMerge overlap merging. The default BBDuk settings include `ftl=20`, `ktrim=r`, `k=23`, `mink=11`, `hdist=1`, `tpe`, `tbo`, and the TruSeq adapter. These are not an additional base-quality threshold.

The native disk-backed collapse counts full merged consensus sequences. Output `R1_seq` is that merged consensus and `R2_seq` its reverse complement, not two independent original reads. Unmerged pairs are counted in provenance but cannot enter assignment. The count table is `<sample>_R1_R2_counts.tsv`; the sidecar distinguishes raw pairs, trimmed pairs, merged reads, retained reads and retained unique sequences. A header-only result can represent an empty input and must be interpreted with those denominators.

## B — Register samples and assays

`magestic-redi-combine-keys` joins the declared sample, gDNA assay and REDI physical-plate keys into `combined_key.tsv`. Literal identifiers are retained. Per-well primer declarations restrict assay expansion where present; the resulting outer-primer/sample identity must be unique. Unregistered or missing physical identifiers fail. Only an explicit `empty` sample is a declared empty control, and a contradictory registration status fails. The engine does not turn malformed plate IDs into a common missing value to make a join succeed.

## C — Extract bc1 and REDI barcode

`magestic-redi-extract-bc` selects the configured bc1 read side from the declared inner-primer rules and searches configured sequence anchors. Both anchor orientation and the configured barcode-length contract matter. The default current cassette is `CAGGTCATGCTC–bc1–CCCTAGGATACG`; older cassettes need their documented per-screen anchors. A positional fallback is disabled by default and, if explicitly enabled, cannot establish cassette identity.

REDI barcoder prefix/length and library come from registered assay metadata. Barcode-to-position keys are validated for unique library/barcode assignment, coordinate shape and compatible prefix lengths. Measured/key columns cannot be overwritten by sample metadata. Explicit sequencing dates must agree with the configuration. Stage C requires all expected count tables unless partial processing is explicitly enabled; missing inputs remain in the recorded denominator.

Outputs are the grouped registered `(bc1, REDI_bc, sample_name)` counts and the separately labelled bc1-only observations. Counts and extraction/assignment categories reconcile to input read counts. Multiple candidate designs do not establish a unique edit.

## D — Observed barcode read fractions

`magestic-redi-purity` produces one observation per sample/barcoder and retains the full observed bc1 count distribution. Its primary fraction is the largest observed barcode count divided by **all observed barcode counts** in that observation. It is a barcode read fraction, not a cell-purity estimate. Close sequence neighbors are not automatically sequencing errors, and truncation to the top four does not define the primary denominator.

Legacy distance-filtered/top-N quantities are separate diagnostics. Tied winners remain ambiguous. An observed off-length winner can remain visible in Stage D; retaining it as evidence does not make it callable in F/G. `all_bc1_counts_json` allows later pooling to reconstruct and validate the full distribution.

## E — Candidate design annotation

`magestic-redi-annotate` requires an explicit bc1 reference table and attaches candidate design information without multiplying measured observation rows. Multiple reference matches are retained as ambiguity with candidate records; they cannot silently duplicate counts or select one design. Design metadata cannot overwrite measured quantities, physical keys or physical-validation flags. A reference match establishes a candidate relationship, not successful editing or independent clone identity.

## F — Compare explicitly selected assays

`magestic-redi-cross-array` processes each configured array pair. Its selectors must distinguish the observations belonging to A and B; the same observation cannot validate itself. Physical joins use the explicit array barcode, 384-well color, row and column. Optional geometry metadata must be uniquely registered.

When pooling observations at a position, the engine combines their complete barcode distributions first and then applies support/callability thresholds. Multiple top-only summaries cannot reconstruct a pooled distribution. Missing, tied, low-support and invalid-length winners remain uncalled. Arrays are outer-joined over physical positions, and agreement is nullable outside the jointly callable denominator. Repeated quadrants or matching barcodes do not establish independent biological replicates.

## G — Registered WGS comparison

`magestic-redi-wgs-check` compares native WGS and REDI observations at independently registered physical positions. WGS is another observation; an expected identity requires independent dated physical/transfer evidence. Barcode-set membership is not a well-level expectation.

The shared WGS call reader accepts the supported native/legacy schemas without padding or truncating observed barcodes. Explicit sample namespace translations may be supplied. A WGS sample must map once to physical coordinates; unregistered samples and ambiguous keys fail. Missing and invalid calls remain visible but unassessed. Every configured array pair is reported, with callable and agreement counts separately recorded. A conflict between valid A/B identities cannot be resolved by automatically favoring B.

Legacy plate-level WGS registration requires an explicit `sample_number` in the sequencing-index geometry. The NNS schema uses its declared `plate_96` key; another supported legacy schema uses `plate_384`. Neither geometry row order nor agreement with an observed bc1 can fill a missing sample-index-to-well map.

## H — Plot registered positions

The `magestic-redi-rearray` compatibility command currently invokes plotting. It does not generate a rearray robot worklist. Stage H handles every configured pair and emits numerical `bc1_registered_position_counts.tsv`, `filtered_registered_positions.tsv`, and PNGs. Multiple pairs receive distinct filename prefixes. The legacy PNG name `bc1_clones_distribution.png` is retained, but its quantity is **registered positions per barcode**, not independent clones.

Defaults are at least20 supporting counts and an observed barcode read fraction of0.95 in a supporting assay. `plotting.min_replicates=0` explicitly means **no independent-replicate gate**. A positive requested gate requires an independently supplied `num_independent_replicates`; absence never means it passed. The engine does not infer independence from quadrant copies or repeated barcodes. Conflicting valid array identities remain unassessed. Candidate reference annotation cannot multiply plotted physical positions.

The completion sidecar belongs to the primary numerical table and records the numerical/PNG output hashes, filters, registered/callable/conflicting/passing position counts and number of distinct passing barcodes. Plot test doubles check numerical data only; visual rendering requires separate inspection.

## I — Selection primitive and remaining robot-output gap

`pipeline/step_i_select.py` exposes `select_top_n_per_group`: it selects up to an explicit positive N per group using declared sort columns, rejects missing group identities and invalid parameters, and preserves input order for exact sorting ties. This is a generic ranking primitive. It does not establish clone eligibility, resolve candidate edits, register physical transfer lineage, assign robot targets or write a robot worklist.

No screen-driver implementation is bundled at the previously advertised `REDI/screens/<screen>/driver.py` path in the audited checkout. Bounded searches were made under both `REDI/` and `src/magestic/` for `driver.py`, `screens/`, robot and rearray implementations. The engine source describes I as “generic clone-selection helpers”; the CLI executes `run_step_h(config)`. The full collection-specific selection/robot-output path therefore remains an implementation gap, not a completed Stage I run.

A concrete implementation needs a declared collection policy (eligibility, per-edit/per-gene caps, tie rules and how ambiguity is handled), explicit source coordinates/lineage, and the robot's target/worklist schema. Estimated work is1–2 engineering days once those artifacts are supplied, including a small independently declared source-to-target fixture and a generated worklist review. This estimate does not count laboratory validation. Generic helper checks cannot substitute for that acceptance test.

## Plate snapshots and workflow status

The PlateDB reader consumes snapshots, never edits the shared workbook, accepts archived primary-ID header variants, rejects duplicate/conflicting IDs and follows only an explicit single-parent column. Cycles, nonintegral IDs, unresolved parent records and depth truncation fail. Branched mating/consolidation lineage and side-by-side spreadsheet tables require explicit role-qualified registration; row adjacency is not a transfer.

The repaired exporter maps both Misc eras explicitly, validates required/unique sheet mappings and primary identities before publication, and attests its TSV/parquet output set. If the new workbook omits an optional stage whose parquet remains in an old output directory, export requires a new destination rather than mixing eras. Legacy snapshot reads do not yet enforce modern completion attestations; that consumer boundary remains under review.

The repaired Snakemake source schedules A through H using the shared typed config and dated registration, declares stage companions, quotes paths, and invokes the same stage implementations as the CLI. Its final target verifies source/input/output content attestations. A missing companion is a scheduling dependency; corrupted/stale attestations fail completion and require an affected-stage rebuild. The declared two-assay fixture verifies numerical CLI/Snakemake parity, unchanged restart behavior, changed-WGS recomputation and a deleted plot companion. Source corrections and prepared assertions alone are not a passing result; consult the audit run ledger for executed outcomes.
