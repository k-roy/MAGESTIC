# Stage 4 functional tests

From the MAGESTIC checkout, with the matching local TRACE checkout available,
Python WGS dependencies and pytest installed, and samtools on PATH:

```sh
python -m pytest -q tests/test_pipelines/wgs/test_reconcile_functional.py
```

These tests generate their own FASTA, BAM, manifest, assignment, bridge, window,
and blacklist inputs. They do not import audit-directory fixtures or require a
cluster or external sequencing files. They exercise the actual Stage 4 entry
point, including refreshed Stage 1 prerequisites and Stage 3 donor cleaning.

The independent raw/cleaned oracle contains four donor-template ALT fragments and
four genomic REF fragments: raw evidence has 4/8 intended ALT, whereas cleaned
evidence has 0/4. Barcode evidence must continue to come from the raw BAM. Further
cases check explicit assignments, conflicting/missing inputs, single-sample and
cohort paths, cache identity, zero coverage after cleaning, and changing the
evidence policy. Tests fail rather than silently skip when samtools is absent.

These synthetic expected-versus-observed cases validate the stated contracts;
they do not supply empirical cohort accuracy. TRACE separately owns the molecular
allele-counting tests in `tests/wgs/test_reconcile_functional.py`.
