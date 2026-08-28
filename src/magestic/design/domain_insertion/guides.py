"""
NGG guide enumeration around a domain-insertion point.

SpCas9 with NGG PAMs, per project decision (WT enzyme for dependable, robust efficiency).

All coordinates are indices into a TRANSCRIPT-ORIENTED genomic window `W`. A "junction index" j
means the position between W[j-1] and W[j]; both the Cas9 blunt cut and the domain-insertion point
are junctions, which keeps the distance arithmetic free of off-by-ones.

The blunt cut lies 6 bp upstream of the END of the PAM, equivalently 3 bp into the spacer from
its PAM-proximal end. In window (sense) coordinates that is `pam_index - 3` for a sense guide and
`pam_index + 6` for an antisense guide -- the asymmetry is real, because the spacer extends in
opposite directions in the two cases.

    sense guide      [--------- 20 nt protospacer ---------][NGG]
                                          ^ cut, 3 nt into the spacer

    antisense guide  [CCN][--------- 20 nt protospacer (rc) ---------]
                                ^ cut, 3 nt into the spacer
"""

from dataclasses import dataclass

from .genome import rc
from .spec import DISALLOWED_IN_DONOR, GUIDE_LENGTH

# A run of >= this many T's terminates Pol III transcription of the sgRNA.
DEFAULT_MAX_T_RUN = 6


@dataclass(frozen=True)
class Guide:
    protospacer: str          # 20 nt, 5'->3' as the sgRNA sees it
    pam: str                  # 3 nt NGG
    strand: str               # '+' = same orientation as the transcript window
    pam_index: int            # index in W where the PAM begins (sense) / the CCN begins (antisense)
    cut_junction: int         # junction index in W of the blunt cut

    @property
    def with_pam(self) -> str:
        return self.protospacer + self.pam


def has_t_run(seq: str, n: int = DEFAULT_MAX_T_RUN) -> bool:
    return "T" * n in seq.upper()


def guide_is_clonable(protospacer: str) -> bool:
    """The guide is synthesised inside the oligo, so it must not carry a cloning site."""
    p = protospacer.upper()
    for site in DISALLOWED_IN_DONOR:
        if site in p or rc(site) in p:
            return False
    return True


def find_guides(window: str, guide_length: int = GUIDE_LENGTH,
                max_t_run: int = DEFAULT_MAX_T_RUN) -> list:
    """Enumerate every NGG guide fully contained in `window`, on both strands."""
    W = window.upper()
    out = []

    # sense: protospacer then NGG
    for p in range(guide_length, len(W) - 2):
        if W[p + 1] == "G" and W[p + 2] == "G":
            proto = W[p - guide_length:p]
            if "N" in proto:
                continue
            out.append(Guide(protospacer=proto, pam=W[p:p + 3], strand="+",
                             pam_index=p, cut_junction=p - 3))

    # antisense: CCN then (rc of the next 20 nt) is the protospacer
    for p in range(0, len(W) - 2 - guide_length):
        if W[p] == "C" and W[p + 1] == "C":
            proto = rc(W[p + 3:p + 3 + guide_length])
            if "N" in proto:
                continue
            # Blunt cut sits 6 bp upstream of the END of the PAM, i.e. 3 bp into the spacer
            # (Kevin, 2026-07-25). For an antisense guide the spacer runs to INCREASING sense
            # coordinates from the PAM, so the cut is 3 nt past the spacer's PAM-proximal end:
            # junction = p + 3 + 3. An earlier version used p + 3, which placed the cut at the
            # PAM/spacer boundary and under-reported every antisense cut distance by 3 nt.
            out.append(Guide(protospacer=proto, pam=rc(W[p:p + 3]), strand="-",
                             pam_index=p, cut_junction=p + 6))

    return [g for g in out
            if not has_t_run(g.protospacer, max_t_run) and guide_is_clonable(g.protospacer)]


def guides_near(window: str, insertion_junction: int, max_distance: int, **kw) -> list:
    """Guides whose blunt cut is within `max_distance` nt of the insertion junction."""
    return sorted(
        (g for g in find_guides(window, **kw)
         if abs(g.cut_junction - insertion_junction) <= max_distance),
        key=lambda g: (abs(g.cut_junction - insertion_junction), g.strand, g.pam_index),
    )


def count_offtargets(protospacer: str, genome: dict, require_pam: bool = True) -> int:
    """
    Exact genome-wide occurrences of the protospacer followed by an NGG PAM.

    A perfect-match count of 1 means the guide is unique. This is the strict screen; the legacy
    MAGESTIC pipeline additionally scores 1- and 2-mismatch neighbours from a precomputed
    Azimuth/BLAST table, which is the better filter once wired in.
    """
    p = protospacer.upper()
    p_rc = rc(p)
    n = 0
    for seq in genome.values():
        for probe, is_rc in ((p, False), (p_rc, True)):
            i = seq.find(probe)
            while i != -1:
                if not require_pam:
                    n += 1
                else:
                    if not is_rc:
                        pam = seq[i + len(probe):i + len(probe) + 3]
                        if len(pam) == 3 and pam[1:] == "GG":
                            n += 1
                    else:
                        pam = seq[i - 3:i]
                        if len(pam) == 3 and rc(pam)[1:] == "GG":
                            n += 1
                i = seq.find(probe, i + 1)
    return n


# ---------------------------------------------------------------------------
# Precomputed off-target screen against the MAGESTIC background strain
# ---------------------------------------------------------------------------
#
# IMPORTANT: guide off-target status for this project is NOT computed against S288C. The
# MAGESTIC background carries designed fixes plus SNVs/indels called by WGS (yKR831), and the
# allowed-guide set was recomputed against THAT genome. Screening against S288C would both admit
# guides that have off-targets in the real background and reject guides that are fine in it.
#
# Canonical file (Oak):
#   .../guide_donor_design/common_scripts_and_reference_files/
#       MAGESTIC_base_strain_allowed_Cas9_20mer_guides_NGG.txt
# Columns: chrom, PAM_coord, strand, full_guide, full_PAM  (876,281 allowed NGG guides)

def load_allowed_guides(path: str, by_site: bool = False):
    """
    Load the precomputed allowed-guide set for the MAGESTIC background strain.

    Returns a set of 20-mer protospacers, or -- with `by_site` -- a dict keyed by
    (chrom, PAM_coord, strand) for position-exact screening.
    """
    import csv as _csv

    opener = open
    if str(path).endswith(".gz"):
        import gzip as _gzip
        opener = _gzip.open

    with opener(path, "rt", newline="") as fh:
        reader = _csv.DictReader(fh, delimiter="\t")
        if by_site:
            return {(r["chrom"], int(r["PAM_coord"]), r["strand"]): r["full_guide"].upper()
                    for r in reader}
        return {r["full_guide"].upper() for r in reader}


def guide_passes_offtarget(protospacer: str, allowed: set) -> bool:
    """Membership in the precomputed allowed set IS the off-target screen."""
    return protospacer.upper() in allowed


# ---------------------------------------------------------------------------
# SCORE: genome-wide ML predictions for every unique NGG sgRNA (Li et al. 2025)
# ---------------------------------------------------------------------------
#
# Efficacy is derived as **1 - SCORE.NE.minmax**, the inverse of the predicted UNEDITED
# likelihood (Kevin, 2026-07-24). Azimuth is deliberately NOT used: it was trained on NHEJ
# indel propensity in human cells, which is the wrong readout for HDR-based donor editing in
# yeast. SCORE was trained genome-wide on this editing system.
#
# Source table (Oak): common/published_data/li_et_al_2025/s288c_NGG_PAM_sites_editing_SCORE.tsv
#   (636 MB, 873,420 sites, 66 columns). A 10-column subset is kept locally.
# Columns used:
#   SCORE.NE.minmax     predicted likelihood the site is left UNEDITED  -> efficacy = 1 - NE
#   SCORE.SV.minmax     predicted structural-variant formation potential (lower is better)
#   SCORE.LD.minmax     predicted large-deletion potential
#   SCORE.NRecT.minmax  predicted non-reciprocal translocation potential
#
# NOTE the table is built on S288C NGG sites, while the library is designed against the MAGESTIC
# background. Joining on the GUIDE SEQUENCE (not coordinates) sidesteps the offset entirely, and
# is the same convention `magestic.utils.harmonize_oligos.assign_sv_scores` already uses.

def load_score_predictions(path: str) -> dict:
    """guide 20-mer -> dict(efficacy, unedited, sv, ld, nrect)."""
    import csv as _csv
    import gzip as _gzip

    opener = _gzip.open if str(path).endswith(".gz") else open
    out = {}
    with opener(path, "rt", newline="") as fh:
        for r in _csv.DictReader(fh, delimiter="\t"):
            try:
                ne = float(r["SCORE.NE.minmax"])
            except (KeyError, ValueError, TypeError):
                continue

            def _f(k):
                try:
                    return float(r[k])
                except (KeyError, ValueError, TypeError):
                    return None

            out[r["guide"].upper()] = {
                "unedited": ne,
                "efficacy": 1.0 - ne,
                "sv": _f("SCORE.SV.minmax"),
                "ld": _f("SCORE.LD.minmax"),
                "nrect": _f("SCORE.NRecT.minmax"),
            }
    return out


def guide_quality(protospacer: str, scores: dict, sv_weight: float = 1.0):
    """
    Composite ranking value: higher is better.

        quality = efficacy - sv_weight * SV
                = (1 - unedited_likelihood) - sv_weight * SV_formation_potential

    Returns None when the guide has no prediction, so callers can decide whether to keep it.
    """
    s = scores.get(protospacer.upper())
    if s is None:
        return None
    sv = s["sv"] or 0.0
    return s["efficacy"] - sv_weight * sv
