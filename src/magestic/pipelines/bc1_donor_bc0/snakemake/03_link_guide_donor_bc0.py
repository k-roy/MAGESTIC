#!/usr/bin/env python3
"""
Link bc1-donor-bc0 data with guide-donor-bc0 data.

This script matches bc1 entries to their corresponding guides by:
1. Loading guide-donor-bc0 data (with donor_bc0_fragment)
2. Creating a 40bp linking key from both datasets:
   - bc1_donor_bc0: donor[-30:] + bc0[:10]
   - guide_donor_bc0: donor_bc0_fragment[:30] + donor_bc0_fragment[-10:]
3. Adding guide and SPS information to bc1 reference

Usage:
    python 03_link_guide_donor_bc0.py \\
        --input <bc1_data.tsv> \\
        --output <output.tsv> \\
        --guide-dir <guide_donor_bc0_dir>          # or: --no-guide-data

2026-09-10 (G3, planning 27_/28_): the guide-donor-bc0 source is now EXPLICIT and
loud. Previously a missing/nonexistent --guide-dir, or a glob matching zero files,
silently copied the input to the output and exited 0 -- so a rebuild could "succeed"
with no guide link at all. Now: --guide-dir is required unless --no-guide-data is
given; a missing dir, zero matched files, or one library_ID appearing in two files
(e.g. both V454_guide_donor_bc0_counts.tsv and V454_matched_guide_donor_bc0_counts.tsv,
as in the superseded in-place 20210822 tree) is an ERROR. Every file loaded is logged
with bytes + mtime and written to <output>.guide_sources.tsv.

Author: Kevin R. Roy
"""

import argparse
import os
import sys
import time
import pandas as pd
from pathlib import Path
from tqdm import tqdm

from magestic.utils.provenance import StepProvenance, read_sidecar

GUIDE_USECOLS = ["library_ID", "guide", "donor_bc0_fragment", "counts", "total_counts", "SPS"]


def _die(msg):
    print(f"ERROR: {msg}", file=sys.stderr, flush=True)
    sys.exit(2)


def _guide_key_series(frag: pd.Series) -> pd.Series:
    """Vectorized create_linking_key_guide: fragment[:30] + fragment[-10:] for 60+ nt fragments."""
    f = frag.fillna("").astype(str).str.upper()
    ok = (f.str.len() >= 60) & (f != "NA")
    key = f.str[:30] + f.str[-10:]
    return key.where(ok)


def load_guide_table(guide_dir: Path, pattern: str, fallback_pattern: str, output: str):
    """Select guides by summed sequence-row support within each source library.

    Keep a real highest-count representative row for SPS/source metadata. Named
    source libraries remain separate support candidates; selection across them
    retains the existing strongest-source policy rather than silently pooling.
    """
    primary, fallback = sorted(guide_dir.glob(pattern)), sorted(guide_dir.glob(fallback_pattern))
    files = primary or fallback
    if primary and fallback:
        lib = lambda f: f.name.split("_")[0]
        orphans = sorted({lib(f) for f in fallback} - {lib(f) for f in primary})
        if orphans:
            _die(f"libraries {orphans} exist only as {fallback_pattern!r} in {guide_dir} while others "
                 f"match {pattern!r}: mixing would silently skip them (review F2). Point --guide-dir "
                 f"at a directory holding ONE generation, or pass --guide-pattern explicitly.")
    if not files:
        _die(f"no guide-donor-bc0 files match {pattern!r} (or {fallback_pattern!r}) in {guide_dir}")
    print(f"\nLoading guide-donor-bc0 data from: {guide_dir} ({len(files)} files)")
    srcs, per_file, lib_to_file = [], [], {}
    for f in tqdm(files, desc="Loading guide data"):
        st = f.stat()
        head = pd.read_csv(f, sep="\t", nrows=0).columns
        use = [c for c in GUIDE_USECOLS if c in head]
        df = pd.read_csv(f, sep="\t", usecols=use, dtype=str, keep_default_na=False)
        cnt_col = "total_counts" if "total_counts" in df.columns else "counts"
        df["_n"] = pd.to_numeric(df[cnt_col], errors="coerce").fillna(0).astype("int64")
        libs = sorted(set(df["library_ID"])) if "library_ID" in df.columns else []
        for lib in libs:
            if lib in lib_to_file:
                _die(f"library_ID {lib} appears in BOTH {lib_to_file[lib]} and {f.name}: "
                     f"two generations of the same table would be double-loaded")
            lib_to_file[lib] = f.name
        df["linking_key"] = _guide_key_series(df["donor_bc0_fragment"])
        df = df[df["linking_key"].notna() & (df["guide"] != "") & (df["guide"] != "NA")]
        # Repeated full-sequence rows can carry the same guide/linking key.
        # Sum their support within the source library before choosing a guide.
        support_keys = ["linking_key", "guide"] + (["library_ID"] if "library_ID" in df.columns else [])
        support = df.groupby(support_keys, sort=False, dropna=False)["_n"].agg(
            guide_support="sum", representative_index="idxmax").reset_index()
        df = df.loc[support["representative_index"]].copy()
        df["_n"] = support["guide_support"].to_numpy()
        # Keep every source-guide candidate until global ambiguity is counted.
        # A discarded per-file runner-up may be absent from every other file.
        keep = ["linking_key", "guide", "_n"] + [c for c in ("library_ID", "SPS") if c in df.columns]
        per_file.append(df[keep])
        srcs.append({"file": str(f), "bytes": st.st_size,
                     "mtime": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(st.st_mtime)),
                     "library_IDs": ";".join(libs), "keys_kept": int(df["linking_key"].nunique()),
                     "guide_candidates": len(df)})
        print(f"  {f.name}: {st.st_size:,} B, mtime {srcs[-1]['mtime']}, libs {libs}, "
              f"{srcs[-1]['keys_kept']:,} keys, {len(df):,} source-guide candidates")
    pd.DataFrame(srcs).to_csv(f"{output}.guide_sources.tsv", sep="\t", index=False)
    g = pd.concat(per_file, ignore_index=True)
    # across files: a key seen in several libraries -> keep best-supported, record ambiguity
    n_guides_x = g.groupby("linking_key")["guide"].nunique()
    g = g.sort_values("_n", ascending=False, kind="mergesort").drop_duplicates("linking_key")
    g["guide_link_n_guides"] = g["linking_key"].map(n_guides_x).astype("int64")
    g = g.rename(columns={"_n": "guide_link_counts"})
    return g, srcs


def _upstream_md5(path: str, nbytes: int):
    """The md5 a guide-donor-bc0 table's OWN sidecar recorded for it, if that sidecar exists and its byte
    count still matches -- so the 11 multi-GB inputs are not re-hashed on every link, yet a table that was
    rewritten since its sidecar (size changed) falls back to a fresh md5. Returns None to force a recompute."""
    sc = Path(f"{path}.provenance.json")
    if not sc.is_file():
        return None
    try:
        for rec in read_sidecar(sc).get("outputs") or []:
            if str(rec.get("path")) == str(path) and int(rec.get("bytes", -1)) == int(nbytes) and rec.get("md5"):
                return rec["md5"]
    except Exception:
        return None
    return None


def create_linking_key_bc1(row):
    """
    Create 40bp linking key from bc1-donor-bc0 data.
    Key = last 30bp of donor + first 10bp of bc0
    """
    donor = row.get("donor") or row.get("top_donor")
    bc0 = row.get("bc0") or row.get("top_bc0")

    if pd.isna(donor) or pd.isna(bc0):
        return None

    donor = str(donor).upper()
    bc0 = str(bc0).upper()

    if len(donor) < 30 or len(bc0) < 10:
        return None

    return donor[-30:] + bc0[:10]


def create_linking_key_guide(row):
    """
    Create 40bp linking key from guide-donor-bc0 data.
    Key = first 30bp of donor_bc0_fragment + last 10bp of donor_bc0_fragment
    (This skips the 20bp SPS in the middle)
    """
    fragment = row.get("donor_bc0_fragment")

    if pd.isna(fragment) or fragment == "NA":
        return None

    fragment = str(fragment).upper()

    if len(fragment) < 60:
        return None

    # donor_bc0_fragment structure: donor[-30:] + SPS[20] + bc0[:10] = 60bp
    # Extract: first 30bp (donor part) + last 10bp (bc0 part)
    return fragment[:30] + fragment[-10:]


def main():
    parser = argparse.ArgumentParser(description="Link bc1 to guide data")
    parser.add_argument("--input", required=True, help="Input bc1 collapsed TSV")
    parser.add_argument("--output", required=True, help="Output linked TSV")
    parser.add_argument("--guide-dir", help="Directory with guide-donor-bc0 files (REQUIRED unless --no-guide-data)")
    parser.add_argument("--no-guide-data", action="store_true",
                        help="Explicitly run WITHOUT guide linkage (recorded in the output log)")
    parser.add_argument("--guide-pattern", default="*_guide_donor_bc0_counts.tsv")
    parser.add_argument("--guide-fallback-pattern", default="*_guide_donor_bc0_purity_filtered.tsv")

    args = parser.parse_args()
    if args.no_guide_data and args.guide_dir:
        _die("give --guide-dir or --no-guide-data, not both")
    if not args.no_guide_data:
        if not args.guide_dir:
            _die("--guide-dir is required (or pass --no-guide-data to skip guide linkage explicitly)")
        if not Path(args.guide_dir).is_dir():
            _die(f"--guide-dir does not exist: {args.guide_dir}")

    # Load bc1 data
    print("Loading bc1-donor-bc0 data...")
    bc1_df = pd.read_csv(args.input, sep="\t", low_memory=False)
    print(f"Loaded {len(bc1_df)} rows with {bc1_df['bc1'].nunique()} unique bc1s")

    if args.no_guide_data:
        print("--no-guide-data: guide linkage SKIPPED by request. Copying input to output.")
        bc1_df.to_csv(args.output, sep="\t", index=False)
        return

    guide_dedup, guide_srcs = load_guide_table(Path(args.guide_dir), args.guide_pattern,
                                               args.guide_fallback_pattern, args.output)
    print(f"Unique guide linking keys: {len(guide_dedup):,} "
          f"(ambiguous, >1 guide: {(guide_dedup['guide_link_n_guides'] > 1).sum():,})")

    print("\nCreating bc1 linking keys...")
    donor = bc1_df["donor"] if "donor" in bc1_df.columns else bc1_df.get("top_donor")
    bc0 = bc1_df["bc0"] if "bc0" in bc1_df.columns else bc1_df.get("top_bc0")
    d = donor.fillna("").astype(str).str.upper()
    b = bc0.fillna("").astype(str).str.upper()
    bc1_df["linking_key"] = (d.str[-30:] + b.str[:10]).where((d.str.len() >= 30) & (b.str.len() >= 10))
    print(f"  BC1 entries with valid linking keys: {bc1_df['linking_key'].notna().sum():,}/{len(bc1_df):,}")

    transfer = ["linking_key", "guide", "guide_link_counts", "guide_link_n_guides"] + \
               [c for c in ("library_ID", "SPS") if c in guide_dedup.columns]
    merged = bc1_df.merge(guide_dedup[transfer], on="linking_key", how="left",
                          suffixes=("", "_from_guide"))
    if len(merged) != len(bc1_df):
        _die(f"row count changed in the guide merge ({len(bc1_df)} -> {len(merged)})")
    matched = merged["guide"].notna().sum()
    print(f"Matched {matched:,}/{len(merged):,} ({100*matched/max(len(merged),1):.1f}%) links to guide data")
    if "donor" in merged.columns:
        dl = merged["donor"].fillna("").astype(str).str.len()
        for cls, m in (("donor<150 (2x150)", dl < 150), ("donor>=150 (2x300)", dl >= 150)):
            n = int(m.sum())
            if n:
                print(f"  {cls}: {merged.loc[m, 'guide'].notna().sum():,}/{n:,} linked")
    if matched == 0:
        _die("ZERO links matched guide data -- wrong --guide-dir, or the linking key is broken")
    if "library_ID" in merged.columns:
        print("\nMatches by library_ID:")
        print(merged["library_ID"].value_counts(dropna=False).head(15))
    merged = merged.drop(columns=["linking_key"])
    merged.to_csv(args.output, sep="\t", index=False)
    print(f"\nSaved linked data to: {args.output}")

    # --- step sidecar (PIPELINE_CORRECTNESS_POLICY.md Sec. 2). The guide tables' md5s come from their own
    # sidecars when present (size-checked), so the provenance chain is joined by checksum, not by path.
    prov = StepProvenance(pipeline="bc1_donor_bc0", step="03_link_guide_donor_bc0", output=args.output,
                          description="bc1-donor-bc0 links joined to the guide-donor-bc0 tables on the 40-nt key "
                                      "donor[-30:]+bc0[:10] == donor_bc0_fragment[:30]+donor_bc0_fragment[-10:]; "
                                      "best-supported guide per key, ambiguity counted")
    prov.add_input(args.input, role="bc1_donor_bc0_collapsed")
    for s in guide_srcs:
        up = _upstream_md5(s["file"], s["bytes"])
        prov.add_input(s["file"], role="guide_donor_bc0_table", md5=up,
                       md5_source="upstream sidecar" if up else "computed here",
                       library_IDs=s["library_IDs"], keys_kept=int(s["keys_kept"]))
    prov.set_parameters(guide_dir=str(args.guide_dir), guide_pattern=args.guide_pattern,
                        guide_fallback_pattern=args.guide_fallback_pattern, guide_usecols=GUIDE_USECOLS,
                        guide_support_scope="sum sequence rows per linking_key/guide/library_ID within each file; source libraries not pooled",
                        guide_metadata_row="highest-count row within selected source-guide group; first appearance on ties",
                        guide_ambiguity_scope="distinct guides across all source-guide candidates, before winner selection",
                        linking_key="donor[-30:] + bc0[:10]  <->  donor_bc0_fragment[:30] + donor_bc0_fragment[-10:]")
    prov.set_counts(links_in=int(len(bc1_df)), links_out=int(len(merged)), unique_bc1=int(bc1_df["bc1"].nunique()),
                    links_with_valid_key=int(bc1_df["linking_key"].notna().sum()), links_guide_linked=int(matched),
                    guide_keys=int(len(guide_dedup)),
                    guide_keys_ambiguous=int((guide_dedup["guide_link_n_guides"] > 1).sum()),
                    guide_files=len(guide_srcs))
    prov.add_output(args.output, rows=int(len(merged)))
    prov.add_output(f"{args.output}.guide_sources.tsv")
    print(f"sidecar: {prov.write()}")


if __name__ == "__main__":
    main()
