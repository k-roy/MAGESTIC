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
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm

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
    """Load, reduce per file (best-supported guide per linking key), fail loudly."""
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
        # best-supported guide per key within this file; remember how many distinct guides
        n_guides = df.groupby("linking_key")["guide"].nunique()
        df = df.sort_values("_n", ascending=False, kind="mergesort").drop_duplicates("linking_key")
        df["_n_guides"] = df["linking_key"].map(n_guides).astype("int64")
        keep = ["linking_key", "guide", "_n", "_n_guides"] + [c for c in ("library_ID", "SPS") if c in df.columns]
        per_file.append(df[keep])
        srcs.append({"file": str(f), "bytes": st.st_size,
                     "mtime": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(st.st_mtime)),
                     "library_IDs": ";".join(libs), "keys_kept": len(df)})
        print(f"  {f.name}: {st.st_size:,} B, mtime {srcs[-1]['mtime']}, libs {libs}, {len(df):,} keys")
    pd.DataFrame(srcs).to_csv(f"{output}.guide_sources.tsv", sep="\t", index=False)
    g = pd.concat(per_file, ignore_index=True)
    # across files: a key seen in several libraries -> keep best-supported, record ambiguity
    n_guides_x = g.groupby("linking_key")["guide"].nunique()
    g = g.sort_values("_n", ascending=False, kind="mergesort").drop_duplicates("linking_key")
    g["guide_link_n_guides"] = np.maximum(g["linking_key"].map(n_guides_x).to_numpy(),
                                          g["_n_guides"].to_numpy())
    g = g.rename(columns={"_n": "guide_link_counts"}).drop(columns=["_n_guides"])
    return g


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

    guide_dedup = load_guide_table(Path(args.guide_dir), args.guide_pattern,
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


if __name__ == "__main__":
    main()
