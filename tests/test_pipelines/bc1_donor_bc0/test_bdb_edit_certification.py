"""Edit certification for the bc1_donor_bc0 matcher (planning 218; production ledger 216 section 3, item 10).

THE DEFECT (BDB-H05, DONOR-WINDOW-R1). The sliding-window branch of
``magestic.pipelines.bc1_donor_bc0.core.matching`` sets ``perfect_middle_region=True`` and can report a unique
installed ``variant`` as soon as ANY 60-mer of the observed donor equals a 60-mer of a designed donor
(``matching.py:655-670``). Nothing asks whether the matched windows contain the edit the design installs. For a
129-nt donor a window that avoids any single offset always exists, so the UNEDITED donor of a design still
matches that design; the audit measured this on 98 real registry rows (70 oligo names), every one of which kept
its oligo, ``perfect_middle_region`` and closest-middle ED0 after its annotated SNVs were reverted to reference.

THE CONTRACT UNDER TEST. A read certifies a design only when the complete designed donor is present
(``full_donor``) or every intended allele of that design lies inside an exactly matched window
(``edit_positions_matched``). Everything else stays a CANDIDATE: names, candidate lists, counts and
``match_method`` are untouched, and only ``perfect_middle_region``, ``is_unambiguous`` and ``variant`` change,
downwards. Design and rationale: ``reviews/bdb_edit_certification/REPAIR_DESIGN.md``.

HOW TO READ THIS FILE.
  ``test_defect_*``    fail on the current matcher BECAUSE of the defect and pass after the repair. They
                       assert only fields that exist today, so the failure is the defect itself.
  ``test_control_*``   must pass BEFORE and AFTER. These are the behaviours the repair may not break.
  ``test_contract_*``  describe the new evidence fields, so before the repair they fail with
                       AttributeError/KeyError. That is expected and is not evidence of the defect.
  ``test_guideonly_*`` encode Decision D-2 (a designed guide is construct evidence, never edit evidence).
                       Drop THIS section, not the others, if Kevin rules the other way.

New names (``oligo_to_edit_positions``, ``edit_evidence``, ``num_edit_certified_variants``) are reached through
attributes rather than imports, and statuses are compared as string literals, so a missing symbol reds only the
test that needs it and never a control at collection time.

No Oak path, no network, no fixture outside tmp_path.
"""

import random

import pandas as pd

from magestic.pipelines.bc1_donor_bc0.core.matching import (
    DEFAULT_WINDOW_SIZE,
    DONOR_VARIABLE_REGION_ANCHORS,
    build_lookup_tables,
    load_and_build_lookup,
    match_dataframe,
    match_sequence,
)

WINDOW = 60
assert WINDOW == DEFAULT_WINDOW_SIZE, "the window arithmetic in these tests assumes a 60-nt window"

COMPLEMENT = str.maketrans("ACGT", "TGCA")
ANCHOR = DONOR_VARIABLE_REGION_ANCHORS[0]
DONOR_START = 1000           # zero-based, half-open: the donor occupies [DONOR_START, DONOR_START + len)
rng = random.Random(21810)


def _background(n):
    """An A/C-only stretch.

    Both pre-trim anchors (GTTTGAAGAGC, TTTCGAAGAGC) carry at least six G/T positions, so no A/C window comes
    within the tolerant anchor search's Hamming budget of 1 and ``trim_to_variable_region`` is a no-op on these
    donors. The tests therefore exercise the matcher, not the trimmer, while running with production defaults.
    """
    return "".join(rng.choice("AC") for _ in range(n))


def _put(seq, offset, base):
    return seq[:offset] + base + seq[offset + 1:]


def _design(name, donor, edits, strand="+", start=DONOR_START, guide=None):
    """One design row in the harmonized-registry schema.

    ``edits`` is ``[(donor offset, the base the REFERENCE allele carries there), ...]``; the designed donor
    already carries the alternate allele at each offset. Coordinates follow the registry convention verified
    in REPAIR_DESIGN.md section 2 on all 98 audited rows: the donor occupies the zero-based half-open interval
    ``[donor_start_coord, donor_end_coord)``, variant positions are one-based, and a '-' donor is the reverse
    complement of that interval.
    """
    end = start + len(donor)
    tokens = []
    for offset, reference_base in edits:
        alternate_base = donor[offset]
        assert alternate_base != reference_base, "a design must differ from the reference allele"
        if strand == "+":
            position, ref, alt = start + offset + 1, reference_base, alternate_base
        else:
            position = end - offset
            ref, alt = reference_base.translate(COMPLEMENT), alternate_base.translate(COMPLEMENT)
        tokens.append(f"{position}_{ref}_{alt}_SNP")
    return {
        "oligo_name": name,
        "guide": guide if guide else _background(20),
        "donor": donor,
        "individual_variants": ",".join(tokens),
        "donor_start_coord": float(start),
        "donor_end_coord": float(end),
        "donor_strand": strand,
    }


def _tables(*rows):
    return build_lookup_tables(pd.DataFrame(list(rows)))


# --- a locus whose design installs an EDGE edit: offset 92, the geometry of the audited outlier rows -------
EDGE = 92
REF_EDGE = _background(129)                        # the reference-allele sequence of this locus
DESIGN_EDGE = _put(REF_EDGE, EDGE, "G")            # the designed donor: G installed at offset 92
EDGE_ROW = _design("design_edge", DESIGN_EDGE, [(EDGE, REF_EDGE[EDGE])])

# --- a locus whose design installs a CENTRED edit, inside the fixed 30:90 middle ---------------------------
CENTRE = 64
REF_CENTRE = _background(129)
DESIGN_CENTRE = _put(REF_CENTRE, CENTRE, "T")
CENTRE_ROW = _design("design_centre", DESIGN_CENTRE, [(CENTRE, REF_CENTRE[CENTRE])])

# --- a locus with TWO intended alleles, the shape of the 13 two-token outlier rows -------------------------
REF_TWO = _background(129)
DESIGN_TWO = _put(_put(REF_TWO, 92, "G"), 94, "T")
TWO_ROW = _design("design_two_alleles", DESIGN_TWO, [(92, REF_TWO[92]), (94, REF_TWO[94])])

# --- two designs at ONE locus, installing different edits --------------------------------------------------
REF_AMB = _background(129)
AMB_A_DONOR = _put(REF_AMB, 30, "T")
AMB_B_DONOR = _put(REF_AMB, 92, "G")
AMB_A_ROW = _design("design_a_edit30", AMB_A_DONOR, [(30, REF_AMB[30])])
AMB_B_ROW = _design("design_b_edit92", AMB_B_DONOR, [(92, REF_AMB[92])])

# --- two designs that share their first 70 nt (one 60-mer window in common) --------------------------------
REF_SHARED = _background(129)
P_DONOR = _put(REF_SHARED, 100, "G")
Q_TAIL = _background(59)
Q_DONOR = _put(P_DONOR[:70] + Q_TAIL, 110, "T")
P_ROW = _design("p_design", P_DONOR, [(100, REF_SHARED[100])])
Q_ROW = _design("q_design", Q_DONOR, [(110, Q_TAIL[110 - 70])])

# --- a donor that matches no design, for the guide-only and no-match paths ---------------------------------
UNRELATED = _background(129)
GUIDE_ONLY_GUIDE = "ACACACACACACACACACAC"
GUIDE_ONLY_ROW = _design("design_guide_only", DESIGN_CENTRE, [(CENTRE, REF_CENTRE[CENTRE])],
                         guide=GUIDE_ONLY_GUIDE)

# --- a design with NO geometry, as legacy design files and today's loader provide ---------------------------
LEGACY_ROW = {"oligo_name": "legacy_design", "guide": "ACACACACACACACACACAA", "donor": DESIGN_CENTRE,
              "individual_variants": "1065_A_T_SNP"}

# --- one REAL registry row pair, so the audited geometry itself is tested ----------------------------------
# Rows 4213 and 4214 of the shipped registry (md5 2dce503d): two guides for one edit, donor on the '-' strand
# of chrIV:[344980, 345109), variant 345017 C>T, so donor offset 345109 - 345017 = 92. Sequences copied from
# 216/preflight/donor_geometry/all_window_outlier_full_rows.json; no Oak path is read here.
USO1_DONOR = ("CCACAGTTTGGCAAGGGCAGATATGATAGGAAGCAGTAACGGCAAGGAAGGATAAGAACATCATTGAGGGAGTCTGTGG"
              "CAGTTTAGCACATACTTTGGACCATTAAAGGGTTACGTAGAGGAGAAGAA")
assert len(USO1_DONOR) == 129 and USO1_DONOR[92] == "A", "the copied registry donor is corrupt"
USO1_ROWS = [
    {"oligo_name": "SpG_Cas9_YDL058W_USO1_chrIV_345017_C_T_AGGG_PAM_-_PAM_strand_345002_PAM_coord"
                   "_guide_disruption_position_15",
     "guide": "CACATGCTTTGGACCATTAA", "donor": USO1_DONOR, "individual_variants": "345017_C_T_SNP",
     "donor_start_coord": 344980.0, "donor_end_coord": 345109.0, "donor_strand": "-"},
    {"oligo_name": "SpG_Cas9_YDL058W_USO1_chrIV_345017_C_T_TGTG_PAM_+_PAM_strand_345018_PAM_coord"
                   "_guide_disruption_position_2",
     "guide": "CCCTTTAATGGTCCAAAGCA", "donor": USO1_DONOR, "individual_variants": "345017_C_T_SNP",
     "donor_start_coord": 344980.0, "donor_end_coord": 345109.0, "donor_strand": "-"},
]
# the audit's counterfactual: the alternate allele T (donor base A) reverted to the reference C (donor base G)
USO1_REVERTED = _put(USO1_DONOR, 92, "G")


# ==========================================================================================================
# test_defect_*  --  fail on the current matcher, pass after the repair
# ==========================================================================================================

def test_defect_unedited_donor_of_an_edge_design_is_not_certified():
    """The read is the design's donor with its one intended allele reverted to reference (offset 92, the
    geometry of the 98 audited rows). 33 windows of the design (starts 0..32) avoid offset 92 and still match,
    so the read is a candidate for the design and no window vouches for the edit.

    Expected: the candidate is kept unchanged -- oligo_name 'design_edge', num_candidates 1,
    num_candidate_variants 1, match_method 'sliding_window', perfect_donor False -- and the certification is
    withdrawn: perfect_middle_region False, is_unambiguous False, variant None.
    Current matcher: perfect_middle_region True, is_unambiguous True and variant '1093_<ref>_G_SNP', i.e. the
    plasmid is certified as carrying an edit the read shows is absent.
    """
    r = match_sequence(REF_EDGE, None, _tables(EDGE_ROW))
    assert r.oligo_name == "design_edge" and r.match_method == "sliding_window"
    assert (r.num_candidate_oligos, r.num_candidate_variants) == (1, 1)
    assert r.perfect_donor is False
    assert r.perfect_middle_region is False
    assert r.is_unambiguous is False
    assert r.variant is None


def test_defect_unedited_donor_of_a_centred_design_is_not_certified():
    """The same with the edit at offset 64, INSIDE the fixed 30:90 middle, which shows the defect is not
    confined to the audited shifted designs: the 10 windows starting 0..4 and 65..69 avoid offset 64, so the
    unedited donor still matches its own design.

    Expected: candidate kept, not certified (perfect_middle_region False, is_unambiguous False, variant None).
    Current matcher: certified.
    """
    r = match_sequence(REF_CENTRE, None, _tables(CENTRE_ROW))
    assert r.oligo_name == "design_centre" and r.match_method == "sliding_window"
    assert r.perfect_middle_region is False
    assert r.is_unambiguous is False
    assert r.variant is None


def test_defect_real_registry_row_uso1_counterfactual_is_not_certified():
    """The audit's own counterfactual on a REAL design pair: USO1 rows 4213/4214, one edit designed with two
    guides, donor offset 92 on the '-' strand, with the alternate allele reverted to reference. 33 design
    windows still match.

    Expected: both names remain candidates (num_candidates 2, num_candidate_variants 1 -- one edit, two
    guides), oligo_name is still the sorted-first candidate, and nothing is certified.
    Current matcher: the audit measured the same oligo, perfect_middle_region and middle ED0 on 98/98 such
    rows (216 register DONOR-WINDOW-R1).
    """
    r = match_sequence(USO1_REVERTED, None, _tables(*USO1_ROWS))
    assert (r.num_candidate_oligos, r.num_candidate_variants) == (2, 1)
    assert r.oligo_name == min(row["oligo_name"] for row in USO1_ROWS)
    assert r.perfect_middle_region is False
    assert r.is_unambiguous is False
    assert r.variant is None


def test_defect_read_truncated_before_the_edit_is_not_certified():
    """The read is the first 85 nt of the designed donor, so the molecule never reaches the intended allele at
    offset 92; 26 windows (starts 0..25) match.

    Expected: candidate kept, not certified. Current matcher: certified, from a molecule that does not even
    contain the position in question.
    """
    r = match_sequence(DESIGN_EDGE[:85], None, _tables(EDGE_ROW))
    assert r.oligo_name == "design_edge" and r.match_method == "sliding_window"
    assert r.perfect_middle_region is False
    assert r.is_unambiguous is False
    assert r.variant is None


def test_defect_design_with_only_one_of_two_alleles_present_is_not_certified():
    """Two intended alleles at offsets 92 and 94 (the shape of the 13 two-token audited rows). The read
    carries the first and the reference base at the second. Windows 33 and 34 cover offset 92 exactly, so a
    rule asking "is ANY intended allele seen" would certify this read.

    Expected: not certified -- EVERY intended allele of a design must be seen. Current matcher: certified.
    """
    r = match_sequence(_put(DESIGN_TWO, 94, REF_TWO[94]), None, _tables(TWO_ROW))
    assert r.oligo_name == "design_two_alleles"
    assert r.perfect_middle_region is False
    assert r.is_unambiguous is False
    assert r.variant is None


def test_defect_match_dataframe_does_not_certify_the_unedited_donor():
    """The same counterfactual through match_dataframe, the entry point step 04 calls
    (04_map_to_oligos.py:171) and whose frame step 04 writes verbatim.

    Expected on that row: oligo_name 'design_edge' kept, perfect_donor False, perfect_middle_region False,
    is_unambiguous False, variant null. Current matcher: True for both flags and a variant reported.
    """
    out = match_dataframe(pd.DataFrame({"donor": [REF_EDGE]}), _tables(EDGE_ROW))
    row = out.iloc[0]
    assert row["oligo_name"] == "design_edge"
    assert not bool(row["perfect_donor"])
    assert not bool(row["perfect_middle_region"])
    assert not bool(row["is_unambiguous"])
    assert pd.isna(row["variant"])


# ==========================================================================================================
# test_guideonly_*  --  Decision D-2: a designed guide is construct evidence, never edit evidence
# ==========================================================================================================

def test_guideonly_match_does_not_certify_a_variant():
    """The donor on this row matched nothing at all (no perfect donor, no window) and only the guide is a
    designed guide, so the matcher falls through to the guide-only branch (matching.py:674-683).

    Expected: the guide answer is kept -- oligo_name 'design_guide_only', match_method 'guide_only' -- but no
    edit is certified: is_unambiguous False, variant None. Current matcher: is_unambiguous True with a variant
    reported, on a row whose donor evidence is nil. The prior contract already states the rule this encodes
    ("Shared 60mer, guide-only ... cannot [support a plasmid-design claim]",
    reviews/bdb_attribution_contract/PLAN_AND_REVIEW.md).
    """
    r = match_sequence(UNRELATED, GUIDE_ONLY_GUIDE, _tables(GUIDE_ONLY_ROW))
    assert r.match_method == "guide_only" and r.oligo_name == "design_guide_only"
    assert r.guide_match is True                       # the guide evidence itself is untouched
    assert r.is_unambiguous is False
    assert r.variant is None


# ==========================================================================================================
# test_control_*  --  must pass BEFORE and AFTER the repair
# ==========================================================================================================

def test_control_perfect_donor_is_certified():
    """An exact designed donor is the certification the repair must leave alone: match_method 'perfect',
    perfect_donor and perfect_middle_region True, is_unambiguous True, the installed variant reported."""
    r = match_sequence(DESIGN_EDGE, None, _tables(EDGE_ROW))
    assert r.match_method == "perfect"
    assert r.perfect_donor is True and r.perfect_middle_region is True
    assert r.is_unambiguous is True
    assert r.variant == EDGE_ROW["individual_variants"]


def test_control_vectorized_fast_path_perfect_donor_is_certified():
    """match_dataframe's anchor-split fast path (matching.py:865-894) carries the perfect-donor majority of
    production rows. A designed donor behind the SpCas9/SpG anchor must still come out certified, by the fast
    path and not the fallback."""
    out = match_dataframe(pd.DataFrame({"donor": ["ACGTACGT" + ANCHOR + DESIGN_EDGE]}), _tables(EDGE_ROW))
    row = out.iloc[0]
    assert row["oligo_name"] == "design_edge" and row["match_method"] == "perfect"
    assert bool(row["perfect_donor"]) and bool(row["perfect_middle_region"])
    assert bool(row["is_unambiguous"])


def test_control_edit_inside_a_matched_window_stays_certified():
    """A read error at offset 5, far from the edit at offset 64: windows 6..64 still match exactly and cover
    the edit, so the design's allele IS seen in exact context. This is the case the repair must NOT withdraw
    -- tolerating an error elsewhere in the donor is why the sliding window exists."""
    r = match_sequence(_put(DESIGN_CENTRE, 5, "G"), None, _tables(CENTRE_ROW))
    assert r.match_method == "sliding_window" and r.perfect_donor is False
    assert r.perfect_middle_region is True
    assert r.is_unambiguous is True
    assert r.variant == CENTRE_ROW["individual_variants"]


def test_control_five_prime_truncation_that_keeps_the_edit_stays_certified():
    """A 20-nt 5' truncation. The read still carries the edit inside exactly matched windows (design starts
    20..64), so it stays certified."""
    r = match_sequence(DESIGN_CENTRE[20:], None, _tables(CENTRE_ROW))
    assert r.perfect_middle_region is True
    assert r.is_unambiguous is True


def test_control_extended_donor_is_certified_through_the_full_designed_donor():
    """A 3'-extended or untrimmed read contains the WHOLE designed donor, which certifies it without needing
    any coordinates. The repo's own test_matching.py::test_extended_donor asserts exactly this
    (perfect_middle_region True for design + extension), so this control protects that test too."""
    r = match_sequence(DESIGN_CENTRE + "ACGT" * 12, None, _tables(CENTRE_ROW))
    assert r.perfect_donor is False
    assert r.perfect_middle_region is True
    assert r.is_unambiguous is True


def test_control_design_without_geometry_is_certified_when_the_full_donor_is_present():
    """Legacy design files (load_oligo_design) and any design whose geometry cannot be projected -- 35.3% of
    the shipped registry -- can never be certified from a window. A read containing the full designed donor
    still certifies them, which is the fallback that keeps those libraries usable."""
    r = match_sequence(DESIGN_CENTRE + "ACGT" * 12, None, build_lookup_tables(pd.DataFrame([LEGACY_ROW])))
    assert r.perfect_middle_region is True
    assert r.is_unambiguous is True


def test_control_variant_ambiguous_window_match_keeps_both_candidates():
    """Two designs share their first 70 nt and install different edits, and the read is a 65-nt fragment of
    that shared region. It is a candidate for both. The repair must not narrow, drop or reorder candidates:
    only the certification claim changes."""
    r = match_sequence(P_DONOR[:65], None, _tables(P_ROW, Q_ROW))
    assert r.match_method == "sliding_window_ambiguous"
    assert (r.num_candidate_oligos, r.num_candidate_variants) == (2, 2)
    assert r.candidate_oligos == ["p_design", "q_design"]
    assert r.oligo_name == "p_design"
    assert r.is_unambiguous is False


def test_control_same_edit_two_guides_still_collapses_to_one_variant():
    """BDB-H06 boundary. Two designs share a donor and install the same edit, so a perfect match is
    variant-unique with num_candidates 2 and no serialized candidate list (format_candidates emits names only
    when there is more than one variant). This patch deliberately does not change that: construct-level
    serialization belongs to the BDB-to-screen projection (ledger item 11). The control exists so that any
    later change to it is visible rather than incidental."""
    r = match_sequence(USO1_DONOR, None, _tables(*USO1_ROWS))
    assert r.match_method == "perfect" and r.is_unambiguous is True
    assert (r.num_candidate_oligos, r.num_candidate_variants) == (2, 1)
    out = match_dataframe(pd.DataFrame({"donor": [USO1_DONOR]}), _tables(*USO1_ROWS))
    assert pd.isna(out.iloc[0]["candidate_oligos"])


def test_control_unmatched_donor_still_reports_nothing():
    """A donor that matches nothing, with no guide: oligo_name None, match_method 'none', not unambiguous."""
    r = match_sequence(UNRELATED, None, _tables(EDGE_ROW))
    assert r.oligo_name is None and r.match_method == "none"
    assert r.is_unambiguous is False


# ==========================================================================================================
# test_contract_*  --  the new evidence fields (these fail before the repair with AttributeError/KeyError)
# ==========================================================================================================

def test_contract_edit_evidence_reports_what_the_read_showed():
    """One status per situation, read off the matcher's own result object:

      full_donor                  the complete designed donor is present (exact, or inside a longer read)
      edit_positions_matched      every intended allele lies inside an exactly matched window
      edit_positions_mismatched   an intended allele is spanned and every matched diagonal shows another base
      edit_positions_not_covered  the read is truncated before the intended allele
      edit_positions_unconfirmed  the read spans it with the design's base but without exact 60-nt context
      not_assessable              this design's intended positions cannot be placed in donor coordinates
      no_donor_evidence           guide-only: the donor showed nothing
    """
    edge, centre = _tables(EDGE_ROW), _tables(CENTRE_ROW)
    legacy = build_lookup_tables(pd.DataFrame([LEGACY_ROW]))
    assert match_sequence(DESIGN_EDGE, None, edge).edit_evidence == "full_donor"
    assert match_sequence(DESIGN_EDGE + "ACGT" * 12, None, edge).edit_evidence == "full_donor"
    assert match_sequence(_put(DESIGN_CENTRE, 5, "G"), None, centre).edit_evidence == "edit_positions_matched"
    assert match_sequence(REF_EDGE, None, edge).edit_evidence == "edit_positions_mismatched"
    assert match_sequence(USO1_REVERTED, None, _tables(*USO1_ROWS)).edit_evidence == "edit_positions_mismatched"
    assert match_sequence(DESIGN_EDGE[:85], None, edge).edit_evidence == "edit_positions_not_covered"
    assert match_sequence(DESIGN_CENTRE[:85], None, legacy).edit_evidence == "not_assessable"
    assert match_sequence(UNRELATED, GUIDE_ONLY_GUIDE, _tables(GUIDE_ONLY_ROW)).edit_evidence == \
        "no_donor_evidence"


def test_contract_unconfirmed_is_the_accepted_conservative_false_negative():
    """The edit IS present at offset 92, but a read error at offset 80 falls inside every window that could
    cover offset 92 (windows 33..69 all span offset 80), so no exact 60-nt context vouches for the allele; the
    only matched windows are starts 0..20.

    The row is reported edit_positions_unconfirmed and is NOT certified. Certifying it would mean trusting an
    extrapolated diagonal, which an undetected indel turns into a one-in-four guess. REPAIR_DESIGN.md section
    4 names this cost explicitly; this test is here so the choice cannot be reversed by accident.
    """
    r = match_sequence(_put(DESIGN_EDGE, 80, "T"), None, _tables(EDGE_ROW))
    assert r.edit_evidence == "edit_positions_unconfirmed"
    assert r.perfect_middle_region is False
    assert r.is_unambiguous is False


def test_contract_build_lookup_tables_projects_intended_alleles_to_donor_offsets():
    """build_lookup_tables must place each design's intended alleles in DONOR coordinates, and refuse to guess
    when it cannot. '+' offset = position - 1 - donor_start_coord; '-' offset = donor_end_coord - position,
    with the designed donor carrying the COMPLEMENT of the alternate allele. Each row is validated against its
    own donor, so a stale or mis-projected annotation yields None instead of pointing at the wrong base.

    The five None cases are the ones that actually occur: 41,452 shipped registry rows have no donor_strand,
    10,193 carry a multi-base or indel allele, and a gapped design fails the interval-length check.
    """
    assert _tables(EDGE_ROW).oligo_to_edit_positions["design_edge"] == (92,)
    assert _tables(TWO_ROW).oligo_to_edit_positions["design_two_alleles"] == (92, 94)
    assert set(_tables(*USO1_ROWS).oligo_to_edit_positions.values()) == {(92,)}   # 345109 - 345017, '-' strand

    def _positions(**changes):
        row = dict(EDGE_ROW)
        row.update(changes)
        return build_lookup_tables(pd.DataFrame([row])).oligo_to_edit_positions["design_edge"]

    assert _positions(donor_strand=float("nan")) is None                    # no donor_strand
    assert _positions(individual_variants="345017_CA_GT_SNP") is None       # multi-base allele
    assert _positions(individual_variants="1093_C_-_INDEL") is None         # indel allele
    assert _positions(donor_end_coord=float(DONOR_START + 130)) is None     # interval disagrees with donor
    assert _positions(individual_variants="1083_A_G_SNP") is None           # design lacks that alt allele


def test_contract_match_dataframe_emits_the_evidence_columns():
    """Step 04 writes whatever match_dataframe returns (04_map_to_oligos.py:171,186), so the evidence has to
    be columns, on both paths: the vectorized fast path (row 0, an anchor-borne perfect donor) and the per-row
    fallback (rows 1 and 2). num_edit_certified_variants counts the distinct installed variants the read
    actually vouches for: every candidate of a perfect donor, and none of an uncertified window match."""
    out = match_dataframe(
        pd.DataFrame({"donor": ["ACGTACGT" + ANCHOR + DESIGN_EDGE, REF_EDGE, DESIGN_EDGE[:85]]}),
        _tables(EDGE_ROW))
    assert out["edit_evidence"].tolist() == ["full_donor", "edit_positions_mismatched",
                                             "edit_positions_not_covered"]
    assert out["num_edit_certified_variants"].tolist() == [1, 0, 0]
    assert out["perfect_middle_region"].tolist() == [True, False, False]
    assert out["is_unambiguous"].tolist() == [True, False, False]


def test_contract_an_ambiguous_row_says_how_many_variants_it_certified():
    """Two designs at one locus, edits at offsets 30 and 92. The read carries the second design's allele, the
    reference base where the first design's allele would be, and one read error at offset 120. Both designs
    stay candidates (2 variants), so the row stays ambiguous -- but exactly one variant is certified, which is
    the fact a later one-row-per-barcode projection (ledger item 11) needs.

    The reported oligo_name is still the sorted-first candidate, here the design the read CONTRADICTS: this
    patch does not narrow by edit evidence (Decision D-3). perfect_middle_region is True because SOME
    candidate's edit-bearing region matched; it certifies oligo_name only when is_unambiguous is True.
    """
    r = match_sequence(_put(AMB_B_DONOR, 120, "T"), None, _tables(AMB_A_ROW, AMB_B_ROW))
    assert (r.num_candidate_oligos, r.num_candidate_variants) == (2, 2)
    assert r.num_edit_certified_variants == 1
    assert r.edit_evidence == "edit_positions_matched"
    assert r.perfect_middle_region is True
    assert r.is_unambiguous is False and r.variant is None
    assert r.oligo_name == "design_a_edit30"


def test_contract_loader_delivers_donor_geometry_to_the_lookup(tmp_path, monkeypatch):
    """The projection is only possible if the BDB loader REQUESTS the geometry. COLUMNS_BC1_DONOR_BC0
    (utils/oligo_annotations.py:103) is handed to pd.read_csv(usecols=...) at oligo_annotations.py:168, so a
    column absent from that list never reaches build_lookup_tables.

    This is the only test that can tell the two halves of the repair apart, and it needs a POSITIVE case to do
    it: a "not certified" assertion would pass without the loader change too, because missing geometry also
    fails closed. Here the read carries the edit inside an exactly matched window, so it must come out
    certified -- which is impossible unless the three geometry columns arrive.
    """
    from magestic.utils import oligo_annotations

    table = {column: "" for column in oligo_annotations.COLUMNS_BC1_DONOR_BC0}
    table.update({key: CENTRE_ROW[key] for key in
                  ("oligo_name", "guide", "donor", "individual_variants",
                   "donor_start_coord", "donor_end_coord", "donor_strand")})
    path = tmp_path / "harmonized_oligos.tsv"
    pd.DataFrame([table]).to_csv(path, sep="\t", index=False)
    monkeypatch.setattr(oligo_annotations, "HARMONIZED_OLIGO_TABLE", path)

    tables = load_and_build_lookup(use_harmonized_table=True)
    assert tables.oligo_to_edit_positions["design_centre"] == (64,), \
        "the loader did not deliver donor_start_coord / donor_end_coord / donor_strand to the lookup"
    r = match_sequence(_put(DESIGN_CENTRE, 5, "G"), None, tables)
    assert r.edit_evidence == "edit_positions_matched"
    assert r.is_unambiguous is True
