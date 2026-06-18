# MAGESTIC TODO

## Preserve local-only work, then it's safe to re-clone-and-delete the M1 copy
**Priority:** Medium (uncommitted work at risk) · **Added:** 2026-06-16 (M1-cleanup pass)
**Context:** `~/repos/magestic` (62 M) has **27 uncommitted/untracked changes** not on
`k-roy/MAGESTIC` — i.e. local-only work that would be lost if the M1 dies or the clone is deleted.
Notable:
- modified: `README.md`, `src/magestic/utils/oligo_annotations.py`, `src/magestic/utils/path_utils.py`
- untracked: `docs/ARCHITECTURE.md`, `.readthedocs.yaml`, `dev/`, and assorted `.DS_Store`
**Action:**
1. Review the 27 changes (`git status`, `git diff`); add `.DS_Store` to `.gitignore`.
2. Commit the real work with explicit `git add <paths>` (not `-A`) and **push to `k-roy/MAGESTIC`**.
3. Once pushed clean, the M1 clone is a re-clonable dead weight and can be removed if desired.
**Note:** MAGESTIC is **SGTC/Stanford** — its canonical home is `k-roy/MAGESTIC` (GitHub) + Oak,
NOT the personal gmail Drive. (Flagged for the CLEANUP_SGTC arm.)
