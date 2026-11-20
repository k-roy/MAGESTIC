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

---

**RESOLVED 2026-11-20.** Verified rather than assumed: `git status --porcelain
--untracked-files=all` is empty, all six files named above are tracked, `.gitignore` covers
`.DS_Store`, and `feat/pool-composition-module` is identical on `origin` (k-roy/MAGESTIC) and
`backup`. Nothing in this clone is local-only, so the "work at risk" this item was opened for is
gone.

⚠️ **The Oak checkout at `Users/kevinroy/software/MAGESTIC` CANNOT be pushed to.** Its git is
1.8.3.1 and its own `receive.denyCurrentBranch = updateInstead` is a value that version does not
understand, so every push dies at receive-pack with `fatal: bad config value`. The working route is
push to GitHub, then `git fetch origin` on Oak. Noted here because the obvious command fails in a
way that looks like a permissions problem and is not.
