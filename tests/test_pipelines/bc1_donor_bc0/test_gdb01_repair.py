"""GDB-01 repair in bdb 04: Cas12a guides that gdb Step 02 extracted 1 nt early are restored (planning 30_)."""
import random

import pandas as pd

from magestic.pipelines.bc1_donor_bc0.core.matching import (
    build_lookup_tables, match_dataframe, repair_gdb01_shifted_guides)

rng = random.Random(9)


def _seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


G_LB = _seq(19) + "G"          # LbCas12a design ending in G
G_SP = _seq(19) + "G"          # SpCas9 design ending in G (never shifted by the defect)
DONOR = _seq(129)
OLIGOS = pd.DataFrame({
    "oligo_name": ["LbCas12a_x", "SpCas9_y", "impLbCas12a_z"],
    "guide": [G_LB, G_SP, _seq(20)],
    # x and y share a donor but install different edits, so only the guide can resolve the match
    "donor": [DONOR, DONOR, _seq(129)],
    "individual_variants": ["10_A_C_SNP", "20_G_T_SNP", "30_C_G_SNP"],
})
TABLES = build_lookup_tables(OLIGOS)


def test_repair_rule():
    raw = pd.Series(["T" + G_LB[:19], "T" + G_SP[:19], G_LB, "A" + G_LB[:19], None])
    out, flag = repair_gdb01_shifted_guides(raw, TABLES)
    assert out.tolist()[:4] == [G_LB, "T" + G_SP[:19], G_LB, "A" + G_LB[:19]]
    assert flag.tolist() == [True, False, False, False, False]


def test_designed_shifted_string_is_never_touched():
    shifted = "T" + G_LB[:19]
    olig = pd.concat([OLIGOS, pd.DataFrame({"oligo_name": ["SpCas9_w"], "guide": [shifted], "donor": [_seq(129)],
                                            "individual_variants": ["40_T_A_SNP"]})], ignore_index=True)
    out, flag = repair_gdb01_shifted_guides(pd.Series([shifted]), build_lookup_tables(olig))
    assert out.iloc[0] == shifted and not flag.iloc[0]


def test_repaired_guide_resolves_the_match():
    shifted = "T" + G_LB[:19]
    df = pd.DataFrame({"donor": [DONOR, DONOR], "guide": [shifted, shifted]})
    before = match_dataframe(df.iloc[:1], TABLES)
    assert before.iloc[0].match_method == "perfect_ambiguous"           # shifted guide narrows nothing
    fixed, _ = repair_gdb01_shifted_guides(df["guide"], TABLES)
    after = match_dataframe(df.assign(guide=fixed).iloc[:1], TABLES)
    r = after.iloc[0]
    assert (r.oligo_name, r.match_method, bool(r.guide_match)) == ("LbCas12a_x", "perfect", True)
