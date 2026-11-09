"""Remove a restriction site from a donor by a SYNONYMOUS codon substitution.

Why this exists
---------------
``construct_donor`` builds a donor centred on the variant and, when the donor contains a
BspQI/SapI site, tries to shift the homology arms until the site falls outside them. When no shift
keeps ``minimum_homology`` on both sides it gives up and the design is lost with
``donor_with_restriction_site_not_allowing_shifting``.

For a donor that lies inside a coding sequence that is an avoidable loss. The site can be destroyed
in place by changing a codon to a synonymous alternative: the protein is unchanged, the intended
edit is untouched, and the donor becomes manufacturable. MAGESTIC donors already carry a deliberate
second change of exactly this kind (the PAM-disrupting substitution), so a site-removing synonymous
change is the same class of edit, not a new one.

Measured motivation: on the TRIUMPH yKR888 library 193 designs were rejected for an unshiftable
site. 83 of them lie in a CDS, and every one of those 83 could be rescued this way -- including
*MIP1/YOR330C*, a read-confirmed (two-sided n=1,552) coding SNV that is one of the DHY214 QTL
fixes. A parameter sweep over ``minimum_homology`` x ``max_dist_to_pam`` recovered none of them,
which is what establishes that this is neither a guide nor a homology problem.

What this module guarantees
---------------------------
Every substitution it returns is checked, not assumed:

* the replacement codon encodes the **same amino acid** (compared after translation);
* the replacement **actually removes** at least one site occurrence and **introduces none**;
* the replacement never touches a base of the intended edit;
* a codon that cannot be changed synonymously (Met, Trp) or whose alternatives all fail is
  reported by returning ``None`` -- the caller then rejects the donor exactly as before.

The function is pure: no file or genome access, so it is testable in isolation.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Sequence, Tuple

__all__ = ["recode_donor_to_remove_sites", "CODON_TO_AA", "AA_TO_CODONS"]

_BASES = "TCAG"
_AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

CODON_TO_AA: Dict[str, str] = {
    _BASES[i // 16] + _BASES[(i // 4) % 4] + _BASES[i % 4]: _AMINO[i] for i in range(64)
}

AA_TO_CODONS: Dict[str, List[str]] = {}
for _codon, _aa in CODON_TO_AA.items():
    AA_TO_CODONS.setdefault(_aa, []).append(_codon)

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _count_sites(seq: str, restriction_sites: Iterable[str]) -> int:
    upper = seq.upper()
    total = 0
    for site in restriction_sites:
        site = site.upper()
        total += upper.count(site) + upper.count(_reverse_complement(site))
    return total


def _cds_offset_of_genomic(
    position: int, cds_blocks: Sequence[Tuple[int, int]], strand: str
) -> Optional[int]:
    """0-based offset of a 1-based genomic position within the spliced CDS, or None."""
    cumulative = 0
    blocks = list(cds_blocks) if strand == "+" else list(reversed(cds_blocks))
    for start, end in blocks:
        if start <= position <= end:
            return cumulative + ((position - start) if strand == "+" else (end - position))
        cumulative += end - start + 1
    return None


def _genomic_of_cds_offset(
    offset: int, cds_blocks: Sequence[Tuple[int, int]], strand: str
) -> Optional[int]:
    """Inverse of :func:`_cds_offset_of_genomic`; returns a 1-based genomic position."""
    cumulative = 0
    blocks = list(cds_blocks) if strand == "+" else list(reversed(cds_blocks))
    for start, end in blocks:
        length = end - start + 1
        if cumulative <= offset < cumulative + length:
            within = offset - cumulative
            return (start + within) if strand == "+" else (end - within)
        cumulative += length
    return None


def recode_donor_to_remove_sites(
    donor: str,
    donor_start: int,
    cds_blocks: Sequence[Tuple[int, int]],
    strand: str,
    cds_sequence: str,
    restriction_sites: Iterable[str],
    protected: Iterable[int] = (),
) -> Tuple[Optional[str], List[Tuple[int, str, str]]]:
    """Destroy restriction sites in ``donor`` with synonymous codon substitutions.

    Parameters
    ----------
    donor:
        The donor sequence as it sits on the PLUS strand of the genome.
    donor_start:
        1-based plus-strand genomic coordinate of ``donor[0]``.
    cds_blocks:
        The gene's CDS blocks as 1-based inclusive ``(start, end)`` pairs, ascending.
    strand:
        ``"+"`` or ``"-"`` -- the gene's strand.
    cds_sequence:
        The spliced CDS in translation orientation (already reverse-complemented for a
        minus-strand gene). Used to read the codon being replaced.
    restriction_sites:
        Recognition sequences to remove. Both orientations are searched.
    protected:
        1-based genomic positions that must not be altered -- the bases of the intended edit.

    Returns
    -------
    ``(recoded_donor, changes)`` where ``changes`` is a list of
    ``(genomic_position, from_base, to_base)``, or ``(None, [])`` when no synonymous
    substitution removes a site without introducing another.
    """
    if _count_sites(donor, restriction_sites) == 0:
        return donor, []

    protected_set = {int(p) for p in protected}
    sites = [s.upper() for s in restriction_sites]
    best_donor = donor
    changes: List[Tuple[int, str, str]] = []

    # one pass per remaining site occurrence; each pass must strictly reduce the count
    for _ in range(8):
        remaining = _count_sites(best_donor, sites)
        if remaining == 0:
            break
        upper = best_donor.upper()
        occurrences: List[int] = []
        for site in sites:
            for pattern in (site, _reverse_complement(site)):
                start = 0
                while True:
                    idx = upper.find(pattern, start)
                    if idx < 0:
                        break
                    occurrences.append(idx)
                    start = idx + 1
        improved = False
        for offset_in_donor in sorted(occurrences):
            site_len = max(len(s) for s in sites)
            for within in range(site_len):
                genomic = donor_start + offset_in_donor + within
                if genomic in protected_set:
                    continue
                cds_offset = _cds_offset_of_genomic(genomic, cds_blocks, strand)
                if cds_offset is None:
                    continue
                codon_index, frame = divmod(cds_offset, 3)
                codon = cds_sequence[codon_index * 3: codon_index * 3 + 3].upper()
                if len(codon) != 3:
                    continue
                amino = CODON_TO_AA.get(codon)
                if amino is None or amino == "*":
                    continue
                alternatives = [c for c in AA_TO_CODONS.get(amino, []) if c != codon]
                if not alternatives:
                    continue  # Met and Trp have no synonym
                codon_genomic = [
                    _genomic_of_cds_offset(codon_index * 3 + k, cds_blocks, strand)
                    for k in range(3)
                ]
                if any(g is None for g in codon_genomic):
                    continue
                if any(g in protected_set for g in codon_genomic):
                    continue
                for candidate in alternatives:
                    if CODON_TO_AA.get(candidate) != amino:
                        continue  # checked by translation, not trusted from the table
                    trial = list(best_donor)
                    ok = True
                    trial_changes: List[Tuple[int, str, str]] = []
                    for k, genomic_pos in enumerate(codon_genomic):
                        index = genomic_pos - donor_start
                        if not 0 <= index < len(trial):
                            ok = False
                            break
                        base = candidate[k] if strand == "+" else _reverse_complement(candidate[k])
                        if trial[index].upper() != base.upper():
                            trial_changes.append((genomic_pos, trial[index].upper(), base.upper()))
                        trial[index] = base
                    if not ok or not trial_changes:
                        continue
                    trial_seq = "".join(trial)
                    if _count_sites(trial_seq, sites) < remaining:
                        best_donor = trial_seq
                        changes.extend(trial_changes)
                        improved = True
                        break
                if improved:
                    break
            if improved:
                break
        if not improved:
            return None, []

    if _count_sites(best_donor, sites) != 0:
        return None, []
    return best_donor, changes


# ---------------------------------------------------------------------------------------------
# Self-test, runnable without pytest: `python -m magestic.design.synonymous_recode` or
# `python src/magestic/design/synonymous_recode.py`. This mirrors the `--selftest` convention the
# adversarial panel scripts already use, and matters because the module must be checkable in
# environments that have no test runner installed (the cluster envs here do not).
# Every success is paired with a PLANTED case that must be REFUSED.
# ---------------------------------------------------------------------------------------------
def _selftest() -> int:  # pragma: no cover - exercised by tests/test_design and by hand
    site = "GCTCTTC"

    def _rc(s: str) -> str:
        return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

    def _tr(s: str) -> str:
        return "".join(CODON_TO_AA.get(s[i: i + 3], "?") for i in range(0, len(s) - 2, 3))

    failures: List[str] = []

    def check(name: str, cond: bool, detail: str = "") -> None:
        if not cond:
            failures.append(f"{name}: {detail}")

    # 1. a donor with no site is returned untouched
    cds = "ATG" + "GAAGAAGAAGAAGAA" + "TAA"
    out, ch = recode_donor_to_remove_sites(cds, 1, [(1, len(cds))], "+", cds, [site])
    check("no-site passthrough", out == cds and ch == [], f"{out!r} {ch!r}")

    # 2. an in-frame site is removed, the protein is preserved
    cds = "ATG" + "GAA" + site + "GAAGAAGAA" + "TAA"
    cds += "A" * ((3 - len(cds) % 3) % 3)
    before = _tr(cds)
    out, ch = recode_donor_to_remove_sites(cds, 1, [(1, len(cds))], "+", cds, [site])
    check("site removed", out is not None and site not in (out or "").upper(), f"{out!r}")
    if out:
        check("protein preserved", _tr(out) == before, f"{_tr(out)} != {before}")
        check("changes reported", bool(ch), "no changes reported for a successful recode")
        check("no new site", _count_sites(out, [site]) == 0, "a new site was introduced")

    # 3. PLANTED: every base protected -> must refuse
    out, ch = recode_donor_to_remove_sites(
        cds, 1, [(1, len(cds))], "+", cds, [site], protected=range(1, len(cds) + 1)
    )
    check("refuses when all bases protected", out is None and ch == [], f"{out!r}")

    # 4. PLANTED: a pattern made only of Met/Trp codons has no synonym -> must refuse
    only_mw = "ATG" + "ATGTGG" + "TAA"
    out, ch = recode_donor_to_remove_sites(
        only_mw, 1, [(1, len(only_mw))], "+", only_mw, ["ATGTGG"]
    )
    check("refuses with no synonym available", out is None and ch == [], f"{out!r}")

    # 5. minus-strand gene: the protein is preserved after re-deriving from the plus strand
    coding = "ATG" + "GAA" + site + "GAAGAAGAA" + "TAA"
    coding += "A" * ((3 - len(coding) % 3) % 3)
    genomic = _rc(coding)
    before = _tr(coding)
    out, ch = recode_donor_to_remove_sites(genomic, 1, [(1, len(genomic))], "-", coding, [site])
    check("minus strand removed", out is not None and site not in (out or "").upper(), f"{out!r}")
    if out:
        check("minus strand protein preserved", _tr(_rc(out)) == before, f"{_tr(_rc(out))}")

    # 6. reported positions are GENOMIC, not donor offsets
    out, ch = recode_donor_to_remove_sites(coding, 101, [(101, 100 + len(coding))], "+",
                                           coding, [site])
    if out:
        check("positions are genomic",
              all(101 <= p <= 100 + len(coding) for p, _f, _t in ch), f"{ch!r}")

    # 7. the table and its inverse agree
    check("table consistent",
          all(c in AA_TO_CODONS[a] for c, a in CODON_TO_AA.items()) and len(CODON_TO_AA) == 64)
    check("only M and W lack synonyms",
          {a for a, c in AA_TO_CODONS.items() if len(c) == 1} == {"M", "W"})

    if failures:
        print(f"SELFTEST FAILED ({len(failures)}):")
        for f in failures:
            print(f"  {f}")
        return 1
    print("synonymous_recode selftest: 7 checks passed, including 2 planted refusals")
    return 0


if __name__ == "__main__":  # pragma: no cover
    import sys as _sys

    _sys.exit(_selftest())
