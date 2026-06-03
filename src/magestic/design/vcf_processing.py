"""
magestic.design.vcf_processing — VCF variant loading and pre-processing.

This module handles the middle steps of the MAGESTIC guide-donor design pipeline:

1. :func:`load_haplotype_vcf` — read the per-sample haplotype TSV produced by
   ``get_haplotype_blocks_from_vcf_with_coord_offset.py``.
2. :func:`combine_linked_variants` — cluster nearby SNPs/INDELs on the same
   haplotype into a single combined variant entry that can be edited by one
   guide-donor pair.
3. :func:`filter_variants` — remove background-strain variants and optionally
   restrict to a gene list (for suppressing non-QTL subpools).

The upstream step (extracting haplotype blocks from the raw multi-sample VCF
with coordinate-offset adjustment for the MAGESTIC background strain) is
project-specific and is documented in the reproducibility wrappers under
``reproduce/``.

Ported from:
  - ``projects/QTL/20240404_SpG_QTL_guide_donor_libraries/final_scripts/
    combine_linked_variants.py``
  - ``projects/QTL/20240404_SpG_QTL_guide_donor_libraries/final_scripts/
    remove_shared_variants_and_filter.py``

Author: Kevin R. Roy <kevinrjroy@gmail.com>
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------------

# Individual variant: (POS, REF, ALT, variant_description, ORF, NAME)
_VarTuple = Tuple[int, str, str, str, str, str]

# variants_by_sample[SAMPLE][HAPLOTYPE][CHROM][POS] = [_VarTuple, ...]
_VariantsBySample = Dict[str, Dict[str, Dict[str, Dict[int, List[_VarTuple]]]]]

# Combined variant dict key: (CHROM, POS, REF, ALT, variant_description, ORF, NAME)
_CombinedKey = Tuple[str, int, str, str, str, str, str]

# Combined variant dict value: list of SAMPLE strings that carry this variant
_CombinedVariants = Dict[_CombinedKey, List[str]]

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Chromosome names for S. cerevisiae nuclear genome (Roman numerals).
YEAST_CHROMS: List[str] = [
    f"chr{n}"
    for n in "I,II,III,IV,V,VI,VII,VIII,IX,X,XI,XII,XIII,XIV,XV,XVI".split(",")
]

#: Map from ``chromosomeN`` (integer) to ``chrROMAN`` format.
_CHROM_NUM_TO_NUMERAL: Dict[str, str] = {
    f"chromosome{i + 1}": f"chr{num}"
    for i, num in enumerate("I,II,III,IV,V,VI,VII,VIII,IX,X,XI,XII,XIII,XIV,XV,XVI".split(","))
}


# ---------------------------------------------------------------------------
# 1. Load haplotype VCF
# ---------------------------------------------------------------------------

def load_haplotype_vcf(filename: str | Path) -> _VariantsBySample:
    """
    Load the per-sample haplotype TSV produced by the haplotype-block extraction
    step (``get_haplotype_blocks_from_vcf_with_coord_offset.py``).

    Expected column order::

        CHROM  POS  REF  ALT  SAMPLE  HAPLOTYPE  PID  genotypes
        all_ALT_alleles  ORF  NAME

    Coordinate trimming (removing shared prefix/suffix bases from REF/ALT) is
    applied here so that variant positions are always the first differing base.

    Parameters
    ----------
    filename : str or Path
        Path to the haplotype block TSV.

    Returns
    -------
    dict
        Nested structure: ``variants_by_sample[SAMPLE][HAPLOTYPE][CHROM][POS]``
        contains a list of ``(POS, REF, ALT, variant_description, ORF, NAME)``
        tuples.
    """
    filename = Path(filename)
    variants_by_sample: _VariantsBySample = {}

    with open(filename) as infile:
        _header = infile.readline()
        for line in infile:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 11:
                logger.debug("Skipping short line: %s", line.rstrip())
                continue
            chrom_raw, pos_str, ref, alt, sample, haplotype = parts[:6]
            orf, name = parts[9], parts[10]

            # Normalise chromosome name format
            chrom = _CHROM_NUM_TO_NUMERAL.get(chrom_raw, chrom_raw)

            # Skip mitochondrial and non-autosomal entries
            if "M" in chrom or chrom not in YEAST_CHROMS:
                continue

            # Skip structural variants
            if alt in ("*", "DUP", "DEL"):
                continue

            pos = int(pos_str)

            # Trim shared suffix bases (left-align the variant)
            while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
                ref = ref[:-1]
                alt = alt[:-1]
            # Trim shared prefix bases
            while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
                pos += 1
                ref = ref[1:]
                alt = alt[1:]

            # Classify variant type
            if len(ref) == 1 and len(alt) == 1:
                vtype = "SNP"
            elif len(ref) == len(alt):
                vtype = "MNP"
            else:
                vtype = "INDEL"

            description = f"{pos}_{ref}_{alt}_{vtype}"
            var_tuple: _VarTuple = (pos, ref, alt, description, orf, name)

            # Insert into nested dict
            variants_by_sample.setdefault(sample, {})
            variants_by_sample[sample].setdefault(haplotype, {})
            variants_by_sample[sample][haplotype].setdefault(chrom, {})

            existing = variants_by_sample[sample][haplotype][chrom]
            if pos in existing:
                logger.warning(
                    "Duplicate variant at %s %s:%d — appending",
                    sample, chrom, pos,
                )
            existing.setdefault(pos, []).append(var_tuple)

    n_samples = len(variants_by_sample)
    n_variants = sum(
        len(pos_lst)
        for s in variants_by_sample.values()
        for h in s.values()
        for c in h.values()
        for pos_lst in c.values()
    )
    logger.info(
        "Loaded %d variants across %d samples from %s",
        n_variants, n_samples, filename,
    )
    return variants_by_sample


# ---------------------------------------------------------------------------
# 2. Combine linked variants
# ---------------------------------------------------------------------------

def combine_linked_variants(
    variants_by_sample: _VariantsBySample,
    genome_seq: dict,
    max_bp_between_snvs: int = 3,
    max_bp_between_indels: int = 10,
) -> _CombinedVariants:
    """
    Cluster variants that are close enough to be edited by a single guide-donor
    pair into combined multi-variant entries.

    SNPs/MNPs within ``max_bp_between_snvs`` bp of each other (or of an INDEL)
    are merged; INDELs within ``max_bp_between_indels`` bp of each other are
    merged.  The combined REF/ALT alleles are assembled from the reference
    genome sequence spanning the cluster.

    Parameters
    ----------
    variants_by_sample : dict
        Output of :func:`load_haplotype_vcf`.
    genome_seq : dict
        Genome sequence dict ``{chrom: sequence_string}`` as returned by
        :func:`magestic.design.guide_donor_functions.load_genome`.
    max_bp_between_snvs : int
        Maximum inter-variant gap (bp) for combining SNPs/MNPs. Default 3.
    max_bp_between_indels : int
        Maximum inter-variant gap (bp) for combining when either variant is an
        INDEL. Default 10.

    Returns
    -------
    dict
        ``combined_variants[(CHROM, POS, REF, ALT, description, ORF, NAME)]``
        = ``[SAMPLE, SAMPLE, ...]`` listing every sample carrying this variant
        (with duplicates if multiple haplotypes in the same sample carry it).

    Notes
    -----
    Single-variant entries are also included in the output so that the dict
    is a complete replacement for the per-sample structure for downstream use.
    """
    combined: _CombinedVariants = {}

    for sample, hap_dict in variants_by_sample.items():
        for haplotype, chrom_dict in hap_dict.items():
            for chrom, pos_dict in chrom_dict.items():
                if chrom not in genome_seq:
                    logger.warning("Chromosome %s not in genome_seq, skipping", chrom)
                    continue

                # Always register single variants first
                for pos, var_list in pos_dict.items():
                    for var_tuple in var_list:
                        key = (chrom, var_tuple[0], var_tuple[1], var_tuple[2],
                               var_tuple[3], var_tuple[4], var_tuple[5])
                        combined.setdefault(key, []).append(sample)

                # Now find and register clusters of ≥2 nearby variants
                sorted_positions = sorted(pos_dict.keys())
                n_pos = len(sorted_positions)

                for left_idx, left_pos in enumerate(sorted_positions):
                    for left_var in pos_dict[left_pos]:
                        left_ref, left_alt, left_desc, orf, name = (
                            left_var[1], left_var[2], left_var[3],
                            left_var[4], left_var[5],
                        )
                        # Use a queue-based BFS to build all valid clusters
                        # starting at this variant
                        initial_cluster = [(left_pos, left_ref, left_alt, left_desc)]
                        initial_end = left_pos + len(left_ref)
                        queue = [(initial_cluster, left_pos, initial_end, left_idx)]

                        while queue:
                            cluster, cluster_start, cluster_end, cur_idx = queue.pop()

                            if len(cluster) >= 2:
                                _register_cluster(
                                    cluster, cluster_start, cluster_end,
                                    chrom, orf, name, genome_seq, combined, sample,
                                )

                            next_idx = cur_idx + 1
                            if next_idx >= n_pos:
                                continue

                            next_pos = sorted_positions[next_idx]
                            gap = next_pos - cluster_end
                            # Choose linking distance based on whether either variant is INDEL
                            is_indel = any("INDEL" in v[3] for v in cluster)
                            max_gap = (
                                max_bp_between_indels if is_indel else max_bp_between_snvs
                            )
                            if gap < max_gap:
                                for next_var in pos_dict[next_pos]:
                                    n_ref, n_alt, n_desc = (
                                        next_var[1], next_var[2], next_var[3]
                                    )
                                    new_cluster = cluster + [
                                        (next_pos, n_ref, n_alt, n_desc)
                                    ]
                                    new_end = next_pos + len(n_ref)
                                    queue.append(
                                        (new_cluster, cluster_start, new_end, next_idx)
                                    )

    n_combined = sum(1 for k in combined if "," in k[4])
    logger.info(
        "combine_linked_variants: %d total entries (%d multi-variant clusters)",
        len(combined), n_combined,
    )
    return combined


def _register_cluster(
    cluster: list,
    cluster_start: int,
    cluster_end: int,
    chrom: str,
    orf: str,
    name: str,
    genome_seq: dict,
    combined: _CombinedVariants,
    sample: str,
) -> None:
    """
    Assemble the combined REF/ALT for a cluster of variants and register it.

    Parameters
    ----------
    cluster : list
        List of ``(POS, REF, ALT, description)`` tuples in genomic order.
    cluster_start, cluster_end : int
        Genomic coordinates spanning the cluster (1-based).
    chrom, orf, name, sample : str
        Metadata forwarded to the combined-variants dict.
    genome_seq : dict
        Genome sequence for REF allele reconstruction.
    combined : dict
        Mutable combined-variants dict (updated in-place).
    """
    # REF is the reference-genome sequence spanning the cluster
    ref_seq = genome_seq[chrom][cluster_start - 1: cluster_end - 1]
    alt_chars = list(ref_seq)

    # Apply each individual variant as a substitution within the REF window
    offset = cluster_start
    descriptions: list[str] = []
    for pos, ind_ref, ind_alt, ind_desc in cluster:
        rel_start = pos - offset
        alt_chars[rel_start: rel_start + len(ind_ref)] = list(ind_alt)
        offset += len(ind_ref) - len(ind_alt)
        descriptions.append(ind_desc)

    alt_seq = "".join(alt_chars)
    combined_desc = ",".join(descriptions)

    key: _CombinedKey = (
        chrom, cluster_start, ref_seq, alt_seq, combined_desc, orf, name,
    )
    combined.setdefault(key, []).append(sample)


# ---------------------------------------------------------------------------
# 3. Filter variants
# ---------------------------------------------------------------------------

def filter_variants(
    combined_variants: _CombinedVariants,
    gene_list: Optional[Set[str]] = None,
    exclude_background: bool = True,
    exclude_gene_names: Optional[Set[str]] = None,
) -> _CombinedVariants:
    """
    Filter the combined-variant dict before oligo design.

    Parameters
    ----------
    combined_variants : dict
        Output of :func:`combine_linked_variants`.
    gene_list : set of str, optional
        If provided, keep only variants whose ORF or NAME is in this set.
        Pass ``None`` to keep all genes. Used to suppress non-QTL subpools.
    exclude_background : bool
        Remove entries whose ORF column is ``"background_variant"`` (i.e.
        variants that distinguish the MAGESTIC background strain from the
        S288c reference). Default ``True``.
    exclude_gene_names : set of str, optional
        Names to exclude (e.g. ``{"HO"}``).

    Returns
    -------
    dict
        Filtered copy of *combined_variants*.
    """
    exclude_gene_names = exclude_gene_names or set()
    filtered: _CombinedVariants = {}
    n_excluded = 0

    for key, samples in combined_variants.items():
        chrom, pos, ref, alt, desc, orf, name = key

        if exclude_background and orf == "background_variant":
            n_excluded += 1
            continue
        if name in exclude_gene_names or orf in exclude_gene_names:
            n_excluded += 1
            continue
        if gene_list is not None and orf not in gene_list and name not in gene_list:
            n_excluded += 1
            continue

        filtered[key] = samples

    logger.info(
        "filter_variants: kept %d / %d variants (excluded %d)",
        len(filtered), len(combined_variants), n_excluded,
    )
    return filtered


def load_gene_list(
    filename: str | Path,
    max_fdr: Optional[float] = None,
) -> Set[str]:
    """
    Load a gene list from a tab-separated or newline-separated file.

    Handles two formats:

    - **Simple two-column TSV** (``systematic_name\\tcommon_name``): both
      columns are collected from every row including the header.
    - **Multi-column TSV with named headers**: if the header contains an
      ``ORF`` column and/or a ``NAME`` column, those columns are extracted
      from every data row.  This covers ``Bloom_QTL_gene_level.tsv`` and
      similar per-gene tables.  If a ``localFDR`` column is present and
      *max_fdr* is given, rows with ``localFDR >= max_fdr`` are skipped.
    - **Single-column list** (no header or only one column): the first column
      of every row is treated as a gene identifier.

    Parameters
    ----------
    filename : str or Path
    max_fdr : float, optional
        If provided and the file has a ``localFDR`` column, only include genes
        whose localFDR is strictly less than this threshold. Typical value:
        ``0.10`` (matches legacy ``remove_shared_variants_and_filter.py``).

    Returns
    -------
    set of str
    """
    genes: Set[str] = set()
    with open(filename) as fh:
        header = fh.readline().strip().split("\t")
        # Check for named ORF / NAME columns (e.g. Bloom_QTL_gene_level.tsv)
        orf_idx = header.index("ORF") if "ORF" in header else None
        name_idx = header.index("NAME") if "NAME" in header else None
        fdr_idx = header.index("localFDR") if "localFDR" in header else None
        if orf_idx is not None or name_idx is not None:
            # Named multi-column table — extract ORF and NAME columns only
            for line in fh:
                parts = line.strip().split("\t")
                if fdr_idx is not None and max_fdr is not None:
                    try:
                        if float(parts[fdr_idx]) >= max_fdr:
                            continue
                    except (ValueError, IndexError):
                        pass
                if orf_idx is not None and orf_idx < len(parts):
                    genes.add(parts[orf_idx])
                if name_idx is not None and name_idx < len(parts):
                    genes.add(parts[name_idx])
        else:
            # Simple two-column or single-column file
            two_col = len(header) == 2
            if two_col:
                genes.update(h for h in header if h)
            else:
                if header[0]:
                    genes.add(header[0])
            for line in fh:
                parts = line.strip().split("\t")
                if parts:
                    genes.add(parts[0])
                if two_col and len(parts) >= 2:
                    genes.add(parts[1])
    genes.discard("")
    logger.info("Loaded %d gene names from %s", len(genes), filename)
    return genes


# ---------------------------------------------------------------------------
# 4. Write processed VCF
# ---------------------------------------------------------------------------

def write_processed_vcf(
    combined_variants: _CombinedVariants,
    output_filename: str | Path,
) -> None:
    """
    Write the combined-variant dict to a tab-separated file for downstream use
    by :func:`magestic.design.oligo_generation.design_oligos_from_vcf`.

    Output columns::

        CHROM  POS  REF  ALT  variant_description  SAMPLES  ORF  NAME

    Parameters
    ----------
    combined_variants : dict
        Output of :func:`combine_linked_variants` or :func:`filter_variants`.
    output_filename : str or Path
        Destination file path.
    """
    rows = []
    for key, samples in combined_variants.items():
        chrom, pos, ref, alt, desc, orf, name = key
        rows.append(
            {
                "CHROM": chrom,
                "POS": pos,
                "REF": ref,
                "ALT": alt,
                "variant_description": desc,
                "SAMPLES": ",".join(sorted(set(samples))),
                "ORF": orf,
                "NAME": name,
            }
        )
    df = pd.DataFrame(rows).sort_values(["CHROM", "POS"])
    df.to_csv(output_filename, sep="\t", index=False)
    logger.info("Wrote %d variants to %s", len(df), output_filename)
