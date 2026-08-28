"""
Regression + negative-control tests for the AsLOV2 domain-insertion oligo design.

The positive fixture is the hand-designed reference oligo for POL3 (`pol3-example-oligo-250-nt`,
Wilson lab / K. Roy 2026-07), which inserts AsLOV2 between POL3 K794 and H795.

Every negative test corrupts exactly one property of that known-good oligo and asserts the
validator rejects it FOR THE RIGHT REASON (checked by error code). A validator that only ever
accepts is worthless, so the negative controls are the substance of this file.
"""

import pytest

from magestic.design.domain_insertion.spec import ASLOV2_AA, ASLOV2_CDS, OligoSpec
from magestic.design.domain_insertion.validate import (
    locate_paqci,
    rc,
    translate,
    validate_oligo,
)

# --- the reference oligo, transcribed from the GenBank record ---------------
# positions 41..290 of POL3_example_oligo_250_nt (the "designed oligo" feature)
POL3_OLIGO = (
    "TGATCGGAGCTGCGATTGGCAGGCGCGCC"[-20:]          # 41..60  5' constant, ends in AscI
    + "TTTTCAAATTCTAAGTTAAT"                        # 61..80  guide
    + "GTTTGAAGAGC"                                 # 81..91  scaffold + BspQI
    + "ACAGATTTAAAGGAAGCTATGGATCTTGGTACCGAAGCTGCCAAATATGTCTCCACTCTATTCAAA"  # 92..157 left arm
    + "GGTTCT"                                      # 158..163 G-S linker
    + "GAGTGCAGGTGGTAGTGCTAACAGGATCACCTGCGCTA"      # 164..201 PaqCI stuffer
    + "TCCGGA"                                      # 202..207 S-G linker
    + "CATCCGATCAATTTGGAGTTCGAGAAAGCATACTTCCCTTACCTTTTGATAAATAAAAAGCGT"     # 208..270 right arm
    + "CGATAGACGACTGGACAGCA"                        # 271..290 subpool rev primer
)

POL3_INSERTION_AFTER = 794          # insert between K794 and H795
POL3_LEFT_ARM_AA = "TDLKEAMDLGTEAAKYVSTLFK"
POL3_RIGHT_ARM_AA = "HPINLEFEKAYFPYLLINKKR"
POL3_NATIVE_WINDOW = POL3_LEFT_ARM_AA + POL3_RIGHT_ARM_AA


@pytest.fixture
def spec():
    return OligoSpec(
        oligo_length=250,
        fwd_const=POL3_OLIGO[:20],
        rev_primer=POL3_OLIGO[-20:],
    )


@pytest.fixture
def native_protein():
    """A synthetic protein containing the POL3 window at the right offset.

    Avoids depending on a genome download in unit tests; the genome-backed check lives in
    the integration test below.
    """
    pad = "A" * (POL3_INSERTION_AFTER - len(POL3_LEFT_ARM_AA))
    return pad + POL3_NATIVE_WINDOW


# --------------------------------------------------------------------------
# structural facts
# --------------------------------------------------------------------------

def test_reference_oligo_is_250_nt():
    assert len(POL3_OLIGO) == 250


def test_spec_budget(spec):
    assert spec.fixed_overhead == 121
    assert spec.arm_budget == 129
    # left arm must be a whole number of codons (consideration 4); matches the hand design
    assert spec.split_arms() == (66, 63)


def test_aslov2_payload_is_consistent():
    assert len(ASLOV2_CDS) == 408
    assert len(ASLOV2_AA) == 136
    assert translate(ASLOV2_CDS) == ASLOV2_AA


def test_aslov2_is_free_of_cloning_sites():
    for site in ("GCTCTTC", "CACCTGC", "GGCGCGCC", "GCGGCCGC"):
        assert site not in ASLOV2_CDS, f"{site} present in AsLOV2"
        assert rc(site) not in ASLOV2_CDS, f"rc({site}) present in AsLOV2"


def test_paqci_overhangs_are_constant(spec):
    """The whole pooled step-3 strategy rests on these two 4-mers."""
    from magestic.design.domain_insertion.validate import extract_donor

    sites = locate_paqci(extract_donor(POL3_OLIGO, spec))
    assert sites.left_overhang == "TTCT"
    assert sites.right_overhang == "TCCG"


# --------------------------------------------------------------------------
# positive control
# --------------------------------------------------------------------------

def test_reference_oligo_validates(spec, native_protein):
    r = validate_oligo(POL3_OLIGO, spec, native_protein=native_protein,
                       insertion_after_residue=POL3_INSERTION_AFTER)
    assert r.ok, r.errors
    assert r.info["residues_added"] == 140
    assert r.info["left_arm_nt"] == 66
    assert r.info["right_arm_nt"] == 63
    assert r.info["edited_protein"] == (
        POL3_LEFT_ARM_AA + "GS" + ASLOV2_AA + "SG" + POL3_RIGHT_ARM_AA)


def test_edited_protein_length(spec, native_protein):
    r = validate_oligo(POL3_OLIGO, spec, native_protein=native_protein,
                       insertion_after_residue=POL3_INSERTION_AFTER)
    assert len(r.info["edited_protein"]) == 22 + 2 + 136 + 2 + 21


# --------------------------------------------------------------------------
# negative controls -- each corrupts one property and must fail with its code
# --------------------------------------------------------------------------

def _codes(result):
    return {e.split(":")[0] for e in result.errors}


def test_C1_second_bspqi_site_rejected(spec, native_protein):
    """A BspQI site inside the donor would be cut during step-1 cloning."""
    bad = POL3_OLIGO[:120] + "GCTCTTC" + POL3_OLIGO[127:]
    r = validate_oligo(bad, spec, native_protein=native_protein,
                       insertion_after_residue=POL3_INSERTION_AFTER)
    assert not r.ok
    assert "C1" in _codes(r)


def test_C2_broken_left_overhang_rejected(spec, native_protein):
    """If the overhang is not TTCT the pooled LOV insert cannot ligate."""
    i = POL3_OLIGO.find("GGTTCT")
    bad = POL3_OLIGO[:i] + "GGTAGT" + POL3_OLIGO[i + 6:]      # G-S -> G-S, but overhang TAGT
    r = validate_oligo(bad, spec, native_protein=native_protein,
                       insertion_after_residue=POL3_INSERTION_AFTER)
    assert not r.ok
    assert "C2" in _codes(r)


def test_C2_extra_paqci_site_rejected(spec, native_protein):
    """A PaqCI site in an arm is cut during step 3 -- the defect the legacy screen misses."""
    bad = POL3_OLIGO[:100] + "CACCTGC" + POL3_OLIGO[107:]
    r = validate_oligo(bad, spec, native_protein=native_protein,
                       insertion_after_residue=POL3_INSERTION_AFTER)
    assert not r.ok
    assert "C2" in _codes(r)


def test_C3_asci_site_in_donor_rejected(spec, native_protein):
    bad = POL3_OLIGO[:100] + "GGCGCGCC" + POL3_OLIGO[108:]
    r = validate_oligo(bad, spec, native_protein=native_protein,
                       insertion_after_residue=POL3_INSERTION_AFTER)
    assert not r.ok
    assert "C3" in _codes(r)


def test_C4_frameshifted_donor_rejected(spec, native_protein):
    """Dropping one nt from the left arm puts the donor out of frame."""
    i = POL3_OLIGO.find("ACAGATTTAAAGG")
    bad = POL3_OLIGO[:i] + POL3_OLIGO[i + 1:]
    r = validate_oligo(bad, spec, native_protein=native_protein,
                       insertion_after_residue=POL3_INSERTION_AFTER)
    assert not r.ok
    assert _codes(r) & {"C4", "C5", "C7"}


def test_C5_nonsynonymous_arm_change_rejected(spec, native_protein):
    """A non-synonymous edit in an arm changes the protein and must be caught."""
    i = POL3_OLIGO.find("ACAGATTTAAAGG")
    bad = POL3_OLIGO[:i] + "CCC" + POL3_OLIGO[i + 3:]          # T -> P at the first codon
    r = validate_oligo(bad, spec, native_protein=native_protein,
                       insertion_after_residue=POL3_INSERTION_AFTER)
    assert not r.ok
    assert "C5" in _codes(r)


def test_C5_synonymous_recoding_is_accepted(spec, native_protein):
    """Synonymous recoding is the whole point of consideration 2 -- it must NOT be rejected."""
    i = POL3_OLIGO.find("ACAGATTTAAAGG")
    assert POL3_OLIGO[i:i + 3] == "ACA"                         # Thr
    good = POL3_OLIGO[:i] + "ACG" + POL3_OLIGO[i + 3:]          # Thr -> Thr
    r = validate_oligo(good, spec, native_protein=native_protein,
                       insertion_after_residue=POL3_INSERTION_AFTER)
    assert r.ok, r.errors


def test_C5_wrong_insertion_residue_rejected(spec, native_protein):
    """Correct oligo, but claimed to insert at the wrong residue."""
    r = validate_oligo(POL3_OLIGO, spec, native_protein=native_protein,
                       insertion_after_residue=POL3_INSERTION_AFTER + 5)
    assert not r.ok
    assert "C5" in _codes(r)


def test_C7_frameshifting_payload_rejected(spec, native_protein):
    """A payload that is not a codon multiple must not silently pass."""
    r = validate_oligo(POL3_OLIGO, spec, native_protein=native_protein,
                       insertion_after_residue=POL3_INSERTION_AFTER,
                       lov_cds=ASLOV2_CDS[:-1], lov_aa=ASLOV2_AA)
    assert not r.ok


def test_C8_guide_still_cutting_donor_rejected(spec, native_protein):
    """
    Place the guide protospacer + NGG inside the right arm. The donor would be re-cut,
    destroying the repaired allele -- this is what the disruption score exists to prevent.
    """
    guide = "TTTTCAAATTCTAAGTTAAT"
    i = POL3_OLIGO.find("CATCCGATCAATT")
    bad = POL3_OLIGO[:i] + guide + "TGG" + POL3_OLIGO[i + len(guide) + 3:]
    r = validate_oligo(bad, spec, native_protein=native_protein,
                       insertion_after_residue=POL3_INSERTION_AFTER)
    assert not r.ok
    assert "C8" in _codes(r)


# --------------------------------------------------------------------------
# genome-backed integration test
# --------------------------------------------------------------------------

@pytest.mark.integration
def test_reference_oligo_against_real_pol3():
    """Full check against SGD R64 POL3 rather than a synthetic protein."""
    import gzip
    import os

    fa = os.environ.get(
        "SGD_ORF_TRANS",
        "/Users/kevinroy/work/UCSB/wilson_lab/planning/data/orf_trans_all.fasta.gz")
    if not os.path.exists(fa):
        pytest.skip("SGD orf_trans_all.fasta.gz not available")

    seqs, name, buf = {}, None, []
    with gzip.open(fa, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).rstrip("*")
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
        seqs[name] = "".join(buf).rstrip("*")

    pol3 = seqs["YDL102W"]
    assert pol3[POL3_INSERTION_AFTER - 22:POL3_INSERTION_AFTER] == POL3_LEFT_ARM_AA
    assert pol3[POL3_INSERTION_AFTER:POL3_INSERTION_AFTER + 21] == POL3_RIGHT_ARM_AA

    s = OligoSpec(oligo_length=250, fwd_const=POL3_OLIGO[:20], rev_primer=POL3_OLIGO[-20:])
    r = validate_oligo(POL3_OLIGO, s, native_protein=pol3,
                       insertion_after_residue=POL3_INSERTION_AFTER)
    assert r.ok, r.errors
    assert r.info["native_context"] == "KYVSTLFK^HPINLEFE"


# --------------------------------------------------------------------------
# Cas9 cut geometry -- 6 bp upstream of the END of the PAM = 3 bp into the spacer
# --------------------------------------------------------------------------

def test_sense_guide_cut_is_3nt_into_the_spacer():
    from magestic.design.domain_insertion.guides import find_guides

    # 20-nt spacer then TGG; PAM starts at index 20, so the cut junction is 17
    W = "A" * 5 + "GTCAGTCAGTCAGTCAGTCA" + "TGG" + "T" * 5
    g = [x for x in find_guides(W) if x.strand == "+" and x.pam == "TGG"][0]
    assert g.pam_index == 25
    assert g.cut_junction == g.pam_index - 3
    # 3 bp into the spacer from its PAM-proximal end
    assert g.protospacer[-3:] == W[g.cut_junction:g.pam_index]


def test_antisense_guide_cut_is_3nt_into_the_spacer():
    """Regression: an earlier version used pam_index + 3, placing the cut at the PAM/spacer
    boundary and under-reporting every antisense cut distance by 3 nt."""
    from magestic.design.domain_insertion.guides import find_guides, rc

    # CCA on the sense strand = TGG PAM on the antisense strand, then a 20-nt spacer
    W = "T" * 5 + "CCA" + "GTCAGTCAGTCAGTCAGTCA" + "A" * 5
    g = [x for x in find_guides(W) if x.strand == "-"][0]
    assert g.pam_index == 5
    assert g.pam == rc("CCA")
    assert g.cut_junction == g.pam_index + 6
    # the 3 spacer bases between the PAM and the cut
    assert W[g.pam_index + 3:g.cut_junction] == "GTC"


def test_pol3_reference_guide_cut_distance():
    """The POL3 reference guide is antisense; its cut sits 9 nt from the insertion point."""
    import gzip
    import os

    import pytest as _pytest

    from magestic.design.domain_insertion.genome import (
        load_chromosomes, load_gene_models)
    from magestic.design.domain_insertion.guides import find_guides

    base = os.environ.get(
        "MAGESTIC_REF_DIR",
        "/Users/kevinroy/work/UCSB/wilson_lab/planning/data/magestic_ref")
    gff = os.path.join(base, "MAGESTIC_background_strain_annotations.gff")
    fsa = os.path.join(base, "MAGESTIC_background_strain.fasta")
    if not (os.path.exists(gff) and os.path.exists(fsa)):
        _pytest.skip("MAGESTIC reference not available")

    models = load_gene_models(gff)
    chroms = load_chromosomes(fsa)
    m = models["YDL102W"]
    W, c = m.genomic_window(chroms[m.chrom], m.cds_to_genomic(794 * 3), 120)
    ins = c + 1
    g = [x for x in find_guides(W) if x.protospacer == "TTTTCAAATTCTAAGTTAAT"][0]
    assert g.strand == "-"
    assert g.cut_junction - ins == 9


def test_linker_can_reconstitute_a_pam_at_the_arm_junction():
    """
    The G-S linker starts with GG. A sense guide whose protospacer ends exactly at the
    arm/linker boundary therefore gets a functional NGG PAM handed to it by the linker, and the
    finished donor would be re-cut. C8 must catch this.
    """
    from magestic.design.domain_insertion.spec import LEFT_LINKER

    assert LEFT_LINKER.startswith("GG"), "this hazard depends on the linker's leading GG"

    guide = "CAAGAACAAAAAATTCTGAA"
    left_arm = ("ACAGATTTAAAGGAAGCTATGGATCTTGGTACCGAAGCTGCC" + guide + "C")[-66:]
    oligo = (POL3_OLIGO[:20] + guide + "GTTTGAAGAGC" + left_arm + "GGTTCT"
             + "GAGTGCAGGTGGTAGTGCTAACAGGATCACCTGCGCTA" + "TCCGGA"
             + "CATCCGATCAATTTGGAGTTCGAGAAAGCATACTTCCCTTACCTTTTGATAAATAAAAAGCGT"
             + POL3_OLIGO[-20:])
    s = OligoSpec(oligo_length=250, fwd_const=POL3_OLIGO[:20], rev_primer=POL3_OLIGO[-20:])
    r = validate_oligo(oligo, s)
    assert not r.ok
    assert any(e.startswith("C8") for e in r.errors), r.errors


def test_spans_insertion_requires_splitting_the_protospacer_not_the_pam():
    from magestic.design.domain_insertion.designer import spans_insertion
    from magestic.design.domain_insertion.guides import Guide

    # sense guide: protospacer [80,100), PAM [100,103), cut at 97
    g = Guide(protospacer="A" * 20, pam="TGG", strand="+", pam_index=100, cut_junction=97)
    assert spans_insertion(g, 20, 90) is True        # inside the protospacer
    assert spans_insertion(g, 20, 101) is False      # only splits the PAM -> not disruptive
