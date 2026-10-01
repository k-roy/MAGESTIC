#!/usr/bin/env python3
"""Driver for the NN-240 known-answer round trip. Each experiment appends its result to
the results JSONL immediately, so a killed session keeps everything already measured.

Usage:  roundtrip_experiments.py <workdir> [experiment ...]
"""
from __future__ import annotations

import json
import sys
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd  # noqa: E402

import synth_bc1_library as S  # noqa: E402
import roundtrip_harness as H  # noqa: E402

RESULTS: Path = None  # set in main


def emit(name: str, payload: dict) -> None:
    rec = {"experiment": name, **payload}
    with open(RESULTS, "a") as f:
        f.write(json.dumps(rec, default=str) + "\n")
    print(f"\n#### {name}\n" + json.dumps(payload, indent=1, default=str)[:4000], flush=True)


def base_design(n_loci=8, vpl=3, gpv=2):
    d = S.make_design(n_loci=n_loci, variants_per_locus=vpl, guides_per_variant=gpv)
    S.plant_bc1s(d)
    return d


def mols_for(d, depth=5):
    return [S.Molecule(bc1=b, oligo=d.by_name(n[0]), depth=depth) for b, n in d.links.items()]


def summarize(r, d, stages=("04", "05")):
    out = {"survival": {k: {"planted": v[0], "spurious": v[1], "rows": v[2]}
                        for k, v in r.survival.items()}}
    for st in stages:
        if st in r.stage_frames:
            sc = H.score_table(r.stage_frames[st], d, st)
            out[f"score_{st}"] = {"denominator": sc.denominator, "recovery": sc.recovery,
                                  "precision_oligo": sc.precision_oligo,
                                  "precision_variant": sc.precision_variant,
                                  "counts": sc.counts, "spurious_bc1": sc.spurious_bc1,
                                  "detail": dict(list(sc.detail.items())[:8])}
    return out


# =================================================================================
# A. clean baseline + non-vacuity control
# =================================================================================
def exp_clean(W):
    d = base_design()
    r = H.run_route_a(d, mols_for(d), W / "clean", n_samples=2)
    res = summarize(r, d)
    # per-stage reported counts, straight from the stage logs
    rep = {}
    t = (W / "clean" / "log_02.txt").read_text()
    for line in t.splitlines():
        if "links," in line and "bc1" in line:
            rep.setdefault("02_lines", []).append(line.strip())
    res["reported_02"] = rep
    res["match_method"] = r.stage_frames["04"]["match_method"].value_counts().to_dict()
    res["guide_match"] = {str(k): int(v) for k, v in
                          r.stage_frames["04"]["guide_match"].value_counts(dropna=False).items()}
    emit("A1_clean_baseline", res)

    # NON-VACUITY: feed 04 the collapsed table (NO guide bridge). The donor alone is 2-way
    # ambiguous at oligo_name because each variant has 2 guides, so if the bridge were
    # inert, precision_oligo would still be 100% -- and the whole test would be vacuous.
    mapped2 = W / "clean" / "mapped_noguide.tsv"
    H.run([H.PY, H.SNK / "04_map_to_oligos.py", "--input", r.collapsed, "--output", mapped2],
          W / "clean" / "common", W / "clean" / "log_04_noguide.txt")
    f2 = pd.read_csv(mapped2, sep="\t", low_memory=False)
    sc = H.score_table(f2, d, "04_no_bridge")
    emit("A2_no_bridge_control", {"line": sc.line(), "counts": sc.counts,
                                  "recovery": sc.recovery,
                                  "precision_oligo": sc.precision_oligo,
                                  "precision_variant": sc.precision_variant,
                                  "num_candidates": {str(k): int(v) for k, v in
                                                     f2["num_candidates"].value_counts().items()}})


# =================================================================================
# B. read-depth curve
# =================================================================================
def exp_depth(W):
    for depth in (1, 2, 5, 20, 100):
        d = base_design(n_loci=5, vpl=2, gpv=2)
        try:
            r = H.run_route_a(d, mols_for(d, depth=depth), W / f"depth{depth}", n_samples=2)
            res = summarize(r, d)
            res["depth_per_bc1_per_sample"] = depth
            res["total_depth_per_bc1"] = depth * 2
            if "04" in r.stage_frames:
                f = r.stage_frames["04"]
                res["support_tier"] = {str(k): int(v) for k, v in
                                       f["support_tier"].value_counts(dropna=False).items()}
                res["num_PCR_replicates"] = {str(k): int(v) for k, v in
                                             f["num_PCR_replicates"].value_counts().items()}
            emit(f"B_depth_{depth}", res)
        except Exception as e:
            emit(f"B_depth_{depth}", {"ERROR": str(e)[:2000]})


# =================================================================================
# C. bc0 collision: two oligos share a bc0 AND the same donor[-30:]
#    -> their 60-nt fragments are identical, so 03's 40-bp linking key collides.
# =================================================================================
def exp_bc0_collision(W, permute_names=False):
    d = base_design(n_loci=4, vpl=2, gpv=1)
    # pick two oligos at the SAME locus (identical flanks -> identical donor[-30:]) and
    # force the same bc0: the fragment then collides exactly.
    a, b = d.oligos[0], d.oligos[1]
    assert a.locus == b.locus and a.donor[-30:] == b.donor[-30:], "need same-locus pair"
    b.bc0 = a.bc0
    if permute_names:
        a.oligo_name, b.oligo_name = b.oligo_name, a.oligo_name
    d.notes["bc0_collision"] = f"{a.oligo_name} / {b.oligo_name} share bc0 {a.bc0}"
    assert a.fragment == b.fragment
    r = H.run_route_a(d, mols_for(d), W / ("bc0col_perm" if permute_names else "bc0col"),
                      n_samples=2)
    res = summarize(r, d)
    f = r.stage_frames["04"]
    coll = f[f["bc1"].isin([k for k, v in d.links.items() if v[0] in (a.oligo_name, b.oligo_name)])]
    res["collided_rows"] = coll[["bc1", "guide", "oligo_name", "guide_link_n_guides",
                                 "guide_status", "is_unambiguous",
                                 "num_candidates"]].to_dict("records")
    res["flag_guide_link_n_guides_survives_to_04"] = "guide_link_n_guides" in f.columns
    res["flag_survives_to_05"] = ("guide_link_n_guides" in r.stage_frames["05"].columns
                                  if "05" in r.stage_frames else None)
    res["collision_note"] = d.notes["bc0_collision"]
    emit("C_bc0_collision" + ("_name_permuted" if permute_names else ""), res)


# =================================================================================
# D. bc1 collision: one bc1 carries two different oligos (a real cloning outcome).
#    HARD form: two variants at the SAME locus, donors differ by 1 nt -- well inside
#    02's --donor-edit-distance 5, so a silent collapse is the prime suspect.
#    CONTROL: two oligos at DIFFERENT loci (donors differ everywhere).
# =================================================================================
def exp_bc1_collision(W):
    for form, same_locus in (("hard_same_locus_1nt", True), ("control_diff_locus", False)):
        d = base_design(n_loci=4, vpl=2, gpv=1)
        if same_locus:
            o1, o2 = d.oligos[0], d.oligos[1]          # same locus, donors differ by 1 nt
        else:
            o1, o2 = d.oligos[0], d.oligos[2]          # different loci
        shared = list(d.links)[0]
        # re-point the shared bc1 at BOTH oligos, and drop o2's own bc1
        own2 = [k for k, v in d.links.items() if v[0] == o2.oligo_name][0]
        del d.links[own2]
        d.links[shared] = [o1.oligo_name, o2.oligo_name]
        mols = [S.Molecule(bc1=b, oligo=d.by_name(n[0]), depth=5) for b, n in d.links.items()]
        mols.append(S.Molecule(bc1=shared, oligo=o2, depth=5))   # the second clone
        hd = sum(x != y for x, y in zip(o1.donor, o2.donor))
        r = H.run_route_a(d, mols, W / f"bc1col_{form}", n_samples=2)
        res = summarize(r, d)
        f = r.stage_frames["04"]
        rows = f[f["bc1"] == shared]
        res["donor_hamming_between_the_two_oligos"] = hd
        res["shared_bc1"] = shared
        res["n_rows_for_shared_bc1"] = len(rows)
        res["rows"] = rows[["donor", "bc0", "oligo_name", "counts", "num_distinct_donors",
                            "bc1_purity", "is_unambiguous", "num_candidates",
                            "match_method"]].assign(donor=lambda x: x["donor"].str[-20:]
                                                    ).to_dict("records")
        res["both_oligos_reported"] = bool({o1.oligo_name, o2.oligo_name} <= set(
            rows["oligo_name"].dropna().astype(str)))
        emit(f"D_bc1_collision_{form}", res)


# =================================================================================
# E. synthesis error in the donor: 1 mm, 2, 3, truncation
# =================================================================================
def exp_synth_error(W):
    for label, kw in (("mm1", {"synth_mismatches": 1}), ("mm2", {"synth_mismatches": 2}),
                      ("mm3", {"synth_mismatches": 3}), ("mm8", {"synth_mismatches": 8}),
                      ("trunc10", {"synth_deletion": 10}), ("trunc40", {"synth_deletion": 40})):
        d = base_design(n_loci=4, vpl=2, gpv=2)
        mols = [S.Molecule(bc1=b, oligo=d.by_name(n[0]), depth=5, **kw)
                for b, n in d.links.items()]
        try:
            r = H.run_route_a(d, mols, W / f"synth_{label}", n_samples=2)
            res = summarize(r, d)
            f = r.stage_frames["04"]
            res["match_method"] = {str(k): int(v) for k, v in
                                   f["match_method"].value_counts(dropna=False).items()}
            res["perfect_donor"] = {str(k): int(v) for k, v in
                                    f["perfect_donor"].value_counts(dropna=False).items()}
            res["guide_status"] = {str(k): int(v) for k, v in
                                   f["guide_status"].value_counts(dropna=False).items()}
            res["donor_len"] = {str(k): int(v) for k, v in
                                f["donor"].astype(str).str.len().value_counts().items()}
            emit(f"E_synth_{label}", res)
        except Exception as e:
            emit(f"E_synth_{label}", {"ERROR": str(e)[:1500],
                                      "TRACE": traceback.format_exc()[-1500:]})


# =================================================================================
# F. mixed well: two bc1s at 50/50, 80/20, 95/5 in the SAME sample
# =================================================================================
def exp_mixed_well(W):
    for major, minor in ((50, 50), (80, 20), (95, 5)):
        d = base_design(n_loci=4, vpl=2, gpv=1)
        keys = list(d.links)
        bA, bB = keys[0], keys[1]
        oA, oB = d.by_name(d.links[bA][0]), d.by_name(d.links[bB][0])
        mols = [S.Molecule(bc1=b, oligo=d.by_name(n[0]), depth=5)
                for b, n in d.links.items() if b not in (bA, bB)]
        mols += [S.Molecule(bc1=bA, oligo=oA, depth=major),
                 S.Molecule(bc1=bB, oligo=oB, depth=minor)]
        r = H.run_route_a(d, mols, W / f"mixed_{major}_{minor}", n_samples=2)
        res = summarize(r, d)
        f = r.stage_frames["04"]
        res["ratio"] = f"{major}/{minor}"
        res["minor_clone_reported"] = bool((f["bc1"] == bB).any())
        sub = f[f["bc1"].isin([bA, bB])]
        res["rows"] = sub[["bc1", "oligo_name", "counts", "bc1_purity", "support_tier",
                           "num_PCR_replicates"]].to_dict("records")
        emit(f"F_mixed_well_{major}_{minor}", res)


# =================================================================================
# G. sequencing error inside the bc1 -- does a 1-error read become a new lineage?
# =================================================================================
def exp_bc1_seq_error(W):
    for nerr in (1, 2):
        d = base_design(n_loci=4, vpl=2, gpv=1)
        mols = []
        for b, n in d.links.items():
            o = d.by_name(n[0])
            mols.append(S.Molecule(bc1=b, oligo=o, depth=20))
            mols.append(S.Molecule(bc1=b, oligo=o, depth=2, bc1_seq_errors=nerr))
        r = H.run_route_a(d, mols, W / f"bc1err{nerr}", n_samples=2)
        res = summarize(r, d)
        res["errors_per_read"] = nerr
        res["parsed_unique_bc1"] = int(r.stage_frames["01"]["bc1"].nunique())
        res["collapsed_unique_bc1"] = int(r.stage_frames["02"]["bc1"].nunique())
        res["planted"] = len(d.links)
        res["richness_inflation_at_01"] = res["parsed_unique_bc1"] / len(d.links)
        res["richness_inflation_at_02"] = res["collapsed_unique_bc1"] / len(d.links)
        emit(f"G_bc1_seq_error_{nerr}", res)


# =================================================================================
# H. edge cases -- assert BEHAVIOR, not merely that it runs
# =================================================================================
def exp_edges(W):
    # H1 empty well: zero reads
    d = base_design(n_loci=2, vpl=1, gpv=1)
    wd = W / "edge_empty"
    wd.mkdir(parents=True, exist_ok=True)
    common = wd / "common"
    S.write_common_dir(d, common)
    for n in ("e_R1.fastq", "e_R2.fastq"):
        (wd / n).write_text("")
    out = wd / "empty_bc1_donor_bc0.tsv"
    p = H.run([H.PY, H.SNK / "01_parse_bc1_donor_bc0.py", wd / "e_R1.fastq",
               wd / "e_R2.fastq", out], common, wd / "log.txt", check=False)
    hdr = out.read_text().splitlines() if out.exists() else []
    meta = json.loads((wd / "empty_bc1_donor_bc0.tsv.parse_meta.json").read_text()) \
        if (wd / "empty_bc1_donor_bc0.tsv.parse_meta.json").exists() else {}
    emit("H1_empty_well", {"rc": p.returncode, "output_exists": out.exists(),
                           "n_lines": len(hdr), "header": hdr[0] if hdr else None,
                           "schema_matches_real": (hdr[0].split("\t") if hdr else None),
                           "parse_meta_note": meta.get("note"),
                           "reads_processed": meta.get("reads_processed")})

    # H2 reads that match no oligo at all
    d2 = base_design(n_loci=4, vpl=1, gpv=1)
    alien = S.make_design(seed=999999, n_loci=2, variants_per_locus=1, guides_per_variant=1)
    S.plant_bc1s(alien, seed=4242)
    mols = mols_for(d2) + [S.Molecule(bc1=b, oligo=alien.by_name(n[0]), depth=5)
                           for b, n in alien.links.items()]
    try:
        r = H.run_route_a(d2, mols, W / "edge_nomatch", n_samples=2)
        f = r.stage_frames["04"]
        unnamed = f[f["oligo_name"].isna()]
        emit("H2_reads_match_no_oligo", {
            "planted_in_design": len(d2.links), "alien_bc1s": len(alien.links),
            **summarize(r, d2),
            "rows_with_no_oligo_name": len(unnamed),
            "alien_bc1s_present_unnamed": int(unnamed["bc1"].isin(alien.links).sum()),
            "alien_bc1s_given_a_name": int(f[f["oligo_name"].notna()]["bc1"]
                                           .isin(alien.links).sum()),
            "match_method_of_alien": {str(k): int(v) for k, v in
                                      f[f["bc1"].isin(alien.links)]["match_method"]
                                      .value_counts(dropna=False).items()}})
    except Exception as e:
        emit("H2_reads_match_no_oligo", {"ERROR": str(e)[:1500]})

    # H3 bc1 that is not exactly N20 (19 and 21)
    for L in (19, 21):
        d3 = base_design(n_loci=3, vpl=1, gpv=1)
        mols = []
        for i, (b, n) in enumerate(d3.links.items()):
            o = d3.by_name(n[0])
            mols.append(S.Molecule(bc1=b, oligo=o, depth=10,
                                   bc1_len_override=(L if i == 0 else None)))
        try:
            r = H.run_route_a(d3, mols, W / f"edge_len{L}", n_samples=2)
            f01 = r.stage_frames["01"]
            lens = f01["bc1"].astype(str).str.len().value_counts().to_dict()
            emit(f"H3_bc1_length_{L}", {"bc1_lengths_at_01": {str(k): int(v) for k, v in lens.items()},
                                        "planted": len(d3.links),
                                        "unique_bc1_at_01": int(f01["bc1"].nunique()),
                                        "unique_bc1_at_02": int(r.stage_frames["02"]["bc1"].nunique()),
                                        **summarize(r, d3)})
        except Exception as e:
            emit(f"H3_bc1_length_{L}", {"ERROR": str(e)[:1500]})

    # H4 both collapsed-fastq producer header formats
    for style in ("plain", "int", "count"):
        d4 = base_design(n_loci=3, vpl=1, gpv=1)
        try:
            r = H.run_route_a(d4, mols_for(d4, depth=7), W / f"edge_hdr_{style}",
                              n_samples=2, header_style=style)
            f01 = r.stage_frames["01"]
            emit(f"H4_fastq_header_{style}", {
                "total_counts_at_01": int(f01["counts"].astype(int).sum()),
                "expected_if_weight_1": 7 * len(d4.links) * 2,
                **summarize(r, d4)})
        except Exception as e:
            emit(f"H4_fastq_header_{style}", {"ERROR": str(e)[:1500]})

    # H5 duplicate oligo names in the design file
    d5 = base_design(n_loci=3, vpl=2, gpv=1)
    d5.oligos[1].oligo_name = d5.oligos[0].oligo_name
    d5.links = {b: [n[0]] for b, n in d5.links.items()}
    try:
        r = H.run_route_a(d5, mols_for(d5), W / "edge_dupnames", n_samples=2)
        emit("H5_duplicate_oligo_names", {"ran": True, **summarize(r, d5),
                                          "note": "two design rows share one oligo_name"})
    except Exception as e:
        emit("H5_duplicate_oligo_names", {"CRASH": str(e)[:2500]})

    # H6 library_ID / contig names containing underscores (03 does f.name.split("_")[0])
    d6 = base_design(n_loci=3, vpl=1, gpv=1)
    for o in d6.oligos:
        o.library = "V545_sub_pool_A"
    try:
        r = H.run_route_a(d6, mols_for(d6), W / "edge_underscore", n_samples=2)
        f03 = r.stage_frames["03"]
        emit("H6_underscore_library_id", {
            "library_ID_in_linked": sorted(set(f03.get("library_ID", pd.Series(dtype=str))
                                               .dropna().astype(str)))[:4],
            **summarize(r, d6)})
    except Exception as e:
        emit("H6_underscore_library_id", {"CRASH": str(e)[:2500]})

    # H7 established_offset mismatch BETWEEN samples (layout C vs layout A)
    d7 = base_design(n_loci=4, vpl=1, gpv=1)
    try:
        r = H.run_route_a(d7, mols_for(d7, depth=10), W / "edge_offset",
                          n_samples=2, tails=[S.TAIL_C, S.TAIL_A])
        metas = [json.loads(Path(str(p) + ".parse_meta.json").read_text()) for p in r.parsed]
        emit("H7_offset_mismatch_between_samples", {
            "offsets": [m.get("established_offset") for m in metas],
            "layouts": [m.get("layout") for m in metas],
            "purities": [m.get("offset_purity") for m in metas],
            "emitting": [m.get("reads_emitting_bc1") for m in metas],
            "processed": [m.get("reads_processed") for m in metas],
            **summarize(r, d7)})
    except Exception as e:
        emit("H7_offset_mismatch_between_samples", {"ERROR": str(e)[:1500]})

    # H8 MIXED layout inside ONE sample (purity < 0.9)
    d8 = base_design(n_loci=4, vpl=1, gpv=1)
    wd = W / "edge_mixedlayout"
    wd.mkdir(parents=True, exist_ok=True)
    common = wd / "common"
    S.write_common_dir(d8, common)
    half = len(d8.links) // 2
    mols = []
    for i, (b, n) in enumerate(d8.links.items()):
        mols.append(S.Molecule(bc1=b, oligo=d8.by_name(n[0]), depth=10,
                               tail=(S.TAIL_C if i < half else S.TAIL_A)))
    reads = S.emit_reads(mols, seed=5)
    S.write_fastq_pair(reads, wd / "m_R1.fastq", wd / "m_R2.fastq")
    out = wd / "mixed_bc1_donor_bc0.tsv"
    p = H.run([H.PY, H.SNK / "01_parse_bc1_donor_bc0.py", wd / "m_R1.fastq",
               wd / "m_R2.fastq", out], common, wd / "log.txt", check=False)
    meta = json.loads(Path(str(out) + ".parse_meta.json").read_text())
    f = pd.read_csv(out, sep="\t", dtype=str)
    emit("H8_mixed_layout_one_sample", {
        "rc": p.returncode, "established_offset": meta.get("established_offset"),
        "offset_purity": meta.get("offset_purity"), "layout": meta.get("layout"),
        "reads_processed": meta.get("reads_processed"),
        "reads_emitting_bc1": meta.get("reads_emitting_bc1"),
        "planted": len(d8.links),
        "planted_recovered_at_01": len(set(f["bc1"]) & set(d8.links)),
        "spurious_at_01": len(set(f["bc1"]) - set(d8.links)),
        "warned_mixed": "MIXED sample" in (wd / "log.txt").read_text()})


# =================================================================================
# I. stage 08: split the clean reads across two runs; integrated == single-run
# =================================================================================
def exp_integrate_08(W):
    d = base_design(n_loci=4, vpl=2, gpv=2)
    one = H.run_route_a(d, mols_for(d), W / "i08_single", n_samples=2, run_05=False)
    wd = W / "i08_two"
    wd.mkdir(parents=True, exist_ok=True)
    common = wd / "common"
    S.write_common_dir(d, common)
    halves = []
    keys = list(d.links)
    for ri, chunk in enumerate((keys[::2], keys[1::2])):
        sub = S.Design(oligos=d.oligos, links={k: d.links[k] for k in chunk})
        rr = H.run_route_a(sub, mols_for(sub), wd / f"run{ri}", n_samples=2, run_05=False)
        halves.append(rr.collapsed)
    out = wd / "integrated.tsv"
    p = H.run([H.PY, H.SNK / "08_integrate_multiple_runs.py", "--input-files", *halves,
               "--output", out, "--n-jobs", "1"], common, wd / "log_08.txt", check=False)
    res = {"rc": p.returncode}
    if p.returncode == 0 and out.exists():
        fi = pd.read_csv(out, sep="\t", low_memory=False)
        f1 = one.stage_frames["02"]
        res.update({
            "integrated_unique_bc1": int(fi["bc1"].nunique()),
            "single_run_unique_bc1": int(f1["bc1"].nunique()),
            "planted": len(d.links),
            "integrated_planted": len(set(fi["bc1"].astype(str)) & set(d.links)),
            "integrated_spurious": len(set(fi["bc1"].astype(str)) - set(d.links)),
            "bc1_sets_equal": set(fi["bc1"].astype(str)) == set(f1["bc1"].astype(str)),
        })
        # does the integrated table still carry the donor so 03/04 can run?
        res["has_donor_col"] = "donor" in fi.columns
    else:
        res["stderr_tail"] = (wd / "log_08.txt").read_text()[-1500:]
    emit("I_stage08_two_runs", res)


# =================================================================================
# J. Route B: collapse_gdb_bridge ON vs OFF, scored against truth
# =================================================================================
def exp_collapse_onoff(W):
    from magestic.pipelines.bc1_donor_bc0.readers.data_loader import (
        load_gdb_data, collapse_gdb_bridge)
    d = base_design(n_loci=6, vpl=2, gpv=2)
    wd = W / "collapse"
    wd.mkdir(parents=True, exist_ok=True)
    # an ambiguous fragment whose WRONG row is FIRST in file order with LOWER support
    victim = d.oligos[0]
    wrong = d.oligos[5].oligo_name
    gdb = wd / "gdb"
    S.write_gdb_dir(d, gdb, reads_per_fragment=3,
                    ambiguous={victim.fragment: wrong})
    res = {}
    for label, flag in (("on", True), ("off", False)):
        frame, info = load_gdb_data(gdb, collapse_bridge=flag)
        res[label] = {"rows": len(frame),
                      "unique_fragments": int(frame["donor_bc0_fragment"].nunique())}
        sel = frame[frame["donor_bc0_fragment"] == victim.fragment]
        res[label]["victim_oligo_names"] = sorted(set(sel["oligo_name"].astype(str)))
        res[label]["victim_rows"] = len(sel)
        res[label]["victim_kept_oligo"] = (str(sel.iloc[0]["oligo_name"])
                                           if len(sel) else None)
    # the support column the collapse actually used
    raw, _ = load_gdb_data(gdb, collapse_bridge=False)
    from magestic.pipelines.bc1_donor_bc0.readers.data_loader import _GDB_SUPPORT_COLUMNS
    res["real_gdb_columns"] = list(raw.columns)
    res["_GDB_SUPPORT_COLUMNS"] = list(_GDB_SUPPORT_COLUMNS)
    res["support_column_found"] = next((c for c in _GDB_SUPPORT_COLUMNS
                                        if c in raw.columns), None)
    _, st = collapse_gdb_bridge(raw)
    res["collapse_stats"] = st
    res["truth_for_victim"] = victim.oligo_name
    res["wrong_oligo_planted_first"] = wrong
    res["collapse_kept_the_truth"] = res["on"]["victim_kept_oligo"] == victim.oligo_name
    emit("J1_collapse_on_off_real_schema", res)

    # same thing with the column the unit tests use, to show the contrast
    gdb2 = wd / "gdb_count_col"
    S.write_gdb_dir(d, gdb2, reads_per_fragment=3, count_col="count",
                    ambiguous={victim.fragment: wrong})
    raw2, _ = load_gdb_data(gdb2, collapse_bridge=False)
    coll2, st2 = collapse_gdb_bridge(raw2)
    kept = coll2[coll2["donor_bc0_fragment"] == victim.fragment]
    emit("J2_collapse_with_count_column", {
        "support_column_found": next((c for c in _GDB_SUPPORT_COLUMNS
                                      if c in raw2.columns), None),
        "victim_kept_oligo": str(kept.iloc[0]["oligo_name"]) if len(kept) else None,
        "truth_for_victim": victim.oligo_name,
        "collapse_kept_the_truth": (str(kept.iloc[0]["oligo_name"]) == victim.oligo_name
                                    if len(kept) else None),
        "stats": st2})

    # J3: gdb table with NO donor_bc0_fragment column -> the 50-nt fallback key
    gdb3 = wd / "gdb_nofrag"
    S.write_gdb_dir(d, gdb3, reads_per_fragment=2, omit_fragment_column=True)
    raw3, _ = load_gdb_data(gdb3, collapse_bridge=False)
    fr = raw3["donor_bc0_fragment"].dropna().astype(str)
    emit("J3_gdb_without_fragment_column", {
        "fallback_built": "donor_bc0_fragment" in raw3.columns,
        "fallback_lengths": {str(k): int(v) for k, v in fr.str.len().value_counts().items()},
        "stage01_fragment_length": 60,
        "matches_stage01_key": bool((fr.str.len() == 60).all()) if len(fr) else None,
        "truth_fragment_example": d.oligos[0].fragment,
        "fallback_example": fr.iloc[0] if len(fr) else None})


# =================================================================================
# K. Route A vs Route B on identical clean input
# =================================================================================
def exp_route_b(W):
    d = base_design(n_loci=5, vpl=2, gpv=2)
    wd = W / "routeb"
    proj = wd / "project"
    (proj / "processed_data" / "bc1_donor_bc0").mkdir(parents=True, exist_ok=True)
    common = wd / "common"
    S.write_common_dir(d, common)
    # Route B reads parsed 2x300 / 2x150 tables from the project tree; build the parsed
    # tables with Route A's 01 so BOTH routes see byte-identical input.
    ra = H.run_route_a(d, mols_for(d), wd / "prep", n_samples=2, run_05=False)
    src = proj / "processed_data" / "bc1_donor_bc0" / "2x300"
    src.mkdir(parents=True, exist_ok=True)
    for p in ra.parsed:
        (src / p.name).write_text(p.read_text())
    gdbdir = proj / "processed_data" / "guide_donor_bc0"
    S.write_gdb_dir(d, gdbdir, reads_per_fragment=3)
    for label, extra in (("default", []), ):
        outd = wd / f"out_{label}"
        log = wd / f"log_routeb_{label}.txt"
        p = H.run([H.PY, "-m", "magestic.pipelines.bc1_donor_bc0.pipeline.run_pipeline",
                   "--project-dir", proj, "--output-dir", outd, "--common-dir", common,
                   "--gdb-dir", gdbdir, "--n-jobs", "1"], common, log, check=False)
        res = {"rc": p.returncode, "log_tail": log.read_text()[-2500:]}
        cands = sorted(outd.rglob("*.tsv")) if outd.exists() else []
        res["outputs"] = [str(c.relative_to(outd)) for c in cands]
        for c in cands:
            try:
                fb = pd.read_csv(c, sep="\t", low_memory=False)
            except Exception:
                continue
            if "bc1" in fb.columns and "oligo_name" in fb.columns:
                sc = H.score_table(fb, d, f"route_b_{c.name}")
                res[f"score_{c.name}"] = {"recovery": sc.recovery,
                                          "precision_oligo": sc.precision_oligo,
                                          "precision_variant": sc.precision_variant,
                                          "counts": sc.counts,
                                          "spurious": sc.spurious_bc1}
        emit(f"K_route_b_{label}", res)


EXPERIMENTS = {
    "clean": exp_clean, "depth": exp_depth, "bc0col": exp_bc0_collision,
    "bc0col_perm": lambda W: exp_bc0_collision(W, permute_names=True),
    "bc1col": exp_bc1_collision, "synth": exp_synth_error, "mixed": exp_mixed_well,
    "bc1err": exp_bc1_seq_error, "edges": exp_edges, "i08": exp_integrate_08,
    "collapse": exp_collapse_onoff, "routeb": exp_route_b,
}


def main():
    global RESULTS
    W = Path(sys.argv[1])
    W.mkdir(parents=True, exist_ok=True)
    RESULTS = W / "results.jsonl"
    names = sys.argv[2:] or list(EXPERIMENTS)
    for n in names:
        print(f"\n{'=' * 70}\n=== {n}\n{'=' * 70}", flush=True)
        try:
            EXPERIMENTS[n](W)
        except Exception as e:
            emit(f"{n}_DRIVER_ERROR", {"ERROR": str(e)[:2000],
                                       "TRACE": traceback.format_exc()[-2500:]})


if __name__ == "__main__":
    main()
