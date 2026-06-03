"""
magestic.design._cli — CLI entry points for MAGESTIC oligo design.

Entry points registered in pyproject.toml:
    magestic-design-vcf        -> design_vcf_main()
    magestic-design-saturation -> design_saturation_main()
    magestic-pam-search        -> pam_search_main()

Author: Kevin R. Roy <kevinrjroy@gmail.com>
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s — %(message)s",
    datefmt="%H:%M:%S",
)


def pam_search_main() -> None:
    """CLI entry point: magestic-pam-search

    Enumerate all PAM sites genome-wide and write a 0/1-mismatch off-target
    analysis file.  Run once per genome/nuclease combination; the output is a
    permanent reference file consumed by magestic-design-vcf.
    """
    parser = argparse.ArgumentParser(
        prog="magestic-pam-search",
        description=(
            "Enumerate all PAM sites in a genome and write a 0/1-mismatch "
            "off-target analysis file."
        ),
    )
    parser.add_argument("--genome", required=True, metavar="FASTA",
                        help="Reference genome FASTA.")
    parser.add_argument("--nuclease", required=True,
                        choices=["SpCas9", "SpG", "LbCas12a", "impLbCas12a"],
                        help="Nuclease system.")
    parser.add_argument("--output", required=True, metavar="TSV",
                        help="Destination for the off-target analysis TSV.")
    parser.add_argument("--seed-length", type=int, default=19, metavar="N",
                        help="Seed region length for 1-mismatch search. (default: 19)")
    args = parser.parse_args()

    from magestic.design.guide_donor_functions import load_genome
    from magestic.design.nuclease import get_nuclease
    from magestic.design.pam_search import search_pam_sites_genomewide

    nuclease = get_nuclease(args.nuclease)
    genome_seq = load_genome(args.genome)
    search_pam_sites_genomewide(
        genome_seq, nuclease, args.output, seed_length=args.seed_length,
    )


def guide_selection_main() -> None:
    """CLI entry point: magestic-guide-selection

    Classify guides from the off-target analysis file into allowed (unique) and
    disallowed (multiple targets) sets.
    """
    parser = argparse.ArgumentParser(
        prog="magestic-guide-selection",
        description="Classify guides from the off-target analysis into allowed/disallowed sets.",
    )
    parser.add_argument("--off-target-file", required=True, metavar="TSV",
                        help="Off-target analysis TSV from magestic-pam-search.")
    parser.add_argument("--allowed-out", required=True, metavar="TSV",
                        help="Output file for allowed guides.")
    parser.add_argument("--disallowed-out", required=True, metavar="TSV",
                        help="Output file for disallowed guides.")
    args = parser.parse_args()

    from magestic.design.guide_selection import assign_allowable_guides

    n_allowed, n_disallowed = assign_allowable_guides(
        args.off_target_file, args.allowed_out, args.disallowed_out,
    )
    print(f"Allowed guides: {n_allowed:,}")
    print(f"Disallowed guides: {n_disallowed:,}")


def design_vcf_main() -> None:
    """CLI entry point: magestic-design-vcf

    Design guide-donor oligo pairs for natural variants from a processed VCF
    file.  Requires pre-computed guide databases from magestic-pam-search and
    magestic-guide-selection.
    """
    parser = argparse.ArgumentParser(
        prog="magestic-design-vcf",
        description=(
            "Design guide-donor 200-bp oligos for natural variants. "
            "Input VCF must be pre-processed by magestic-vcf-process."
        ),
    )
    parser.add_argument("--vcf", required=True, metavar="TSV",
                        help="Processed variant TSV from magestic-vcf-process.")
    parser.add_argument("--genome", required=True, metavar="FASTA",
                        help="Reference genome FASTA.")
    parser.add_argument("--allowed-guides", required=True, metavar="TSV",
                        help="Allowed guides file from magestic-guide-selection.")
    parser.add_argument("--disallowed-guides", required=True, metavar="TSV",
                        help="Disallowed guides file from magestic-guide-selection.")
    parser.add_argument("--plus-bedgraph", required=True, metavar="BG",
                        help="TIF-seq plus-strand bedgraph (for donor orientation).")
    parser.add_argument("--minus-bedgraph", required=True, metavar="BG",
                        help="TIF-seq minus-strand bedgraph.")
    parser.add_argument("--rev-primers", required=True, metavar="TXT",
                        help="File of reverse subpool priming sequences (one per line).")
    parser.add_argument("--nuclease", required=True,
                        choices=["SpCas9", "SpG", "LbCas12a", "impLbCas12a"],
                        help="Nuclease system.")
    parser.add_argument("--output-dir", required=True, metavar="DIR",
                        help="Output directory for oligo pool, info, and untargetable TSVs.")
    parser.add_argument("--output-prefix", default="oligo_design", metavar="STR",
                        help="Filename prefix for output files. (default: oligo_design)")
    parser.add_argument("--subpool-idx", type=int, default=0, metavar="N",
                        help="Starting subpool index into the rev primer list. (default: 0)")
    parser.add_argument("--oligo-length", type=int, default=200, metavar="N",
                        help="Target total oligo length in bases. (default: 200)")
    parser.add_argument("--minimum-homology", type=int, default=30, metavar="N",
                        help="Minimum donor homology arm length. (default: 30)")
    parser.add_argument("--max-guides-per-donor", type=int, default=4, metavar="N",
                        help="Maximum guide-donor pairs per variant. (default: 4)")
    parser.add_argument("--max-dist-to-pam", type=int, default=20, metavar="N",
                        help="Maximum variant-to-PAM distance in bp. (default: 20)")
    parser.add_argument("--fasta-exclusion", nargs="*", metavar="FASTA",
                        help="FASTA files whose guides must be excluded from the design.")
    args = parser.parse_args()

    from magestic.design.guide_donor_functions import load_genome, load_transcript_coverage
    from magestic.design.nuclease import get_nuclease, expand_pam_motif
    from magestic.design.oligo_generation import design_oligos_from_vcf

    nuclease = get_nuclease(args.nuclease)
    genome_seq = load_genome(args.genome)

    plus_bg = load_transcript_coverage(args.plus_bedgraph)
    minus_bg = load_transcript_coverage(args.minus_bedgraph)

    rev_primers = []
    with open(args.rev_primers) as fh:
        for line in fh:
            rev_primers.append(line.strip().lower())

    # Default PAM preference: NGNG (higher activity) before NGNH for SpG
    if args.nuclease == "SpG":
        pam_preference = [expand_pam_motif("NGNG"), expand_pam_motif("NGNH")]
    else:
        pam_preference = [nuclease.all_pams]

    out_dir = Path(args.output_dir)
    prefix = args.output_prefix

    design_oligos_from_vcf(
        processed_vcf_file=args.vcf,
        genome_seq=genome_seq,
        nuclease=nuclease,
        allowed_guides_file=args.allowed_guides,
        disallowed_guides_file=args.disallowed_guides,
        plus_strand_bedgraph=plus_bg,
        minus_strand_bedgraph=minus_bg,
        rev_primer_sequences=rev_primers,
        output_oligo_pool_file=out_dir / f"{prefix}_oligo_pool.tsv",
        output_info_file=out_dir / f"{prefix}_info.tsv",
        output_untargetable_file=out_dir / f"{prefix}_untargetable.tsv",
        fasta_exclusion_files=args.fasta_exclusion or [],
        subpool_idx=args.subpool_idx,
        pam_preference_order=pam_preference,
        oligo_length=args.oligo_length,
        minimum_homology=args.minimum_homology,
        max_guides_per_donor=args.max_guides_per_donor,
        max_dist_to_pam=args.max_dist_to_pam,
        nuclease_label=args.nuclease,
    )


def design_vcf_preprocess_main() -> None:
    """CLI entry point: magestic-vcf-process

    Load a haplotype VCF, combine linked variants, filter, and write a
    processed TSV for use by magestic-design-vcf.
    """
    parser = argparse.ArgumentParser(
        prog="magestic-vcf-process",
        description="Load, combine, and filter variants from a haplotype VCF for oligo design.",
    )
    parser.add_argument("--haplotype-vcf", required=True, metavar="TSV",
                        help="Per-sample haplotype TSV from get_haplotype_blocks_from_vcf_with_coord_offset.py.")
    parser.add_argument("--genome", required=True, metavar="FASTA",
                        help="Reference genome FASTA (needed for combined REF/ALT assembly).")
    parser.add_argument("--output", required=True, metavar="TSV",
                        help="Output processed variant TSV.")
    parser.add_argument("--gene-list", metavar="TSV",
                        help="Optional gene list TSV to restrict to QTL genes.")
    parser.add_argument("--max-fdr", type=float, default=None, metavar="FDR",
                        help="If --gene-list has a localFDR column, only include genes with "
                             "localFDR < MAX_FDR. Default: no FDR filter. Typical value: 0.10.")
    parser.add_argument("--exclude-background", action="store_true", default=True,
                        help="Exclude background_variant ORF entries. (default: True)")
    parser.add_argument("--exclude-genes", nargs="*", metavar="GENE",
                        help="Gene names to exclude (e.g. HO).")
    parser.add_argument("--max-bp-between-snvs", type=int, default=3, metavar="N",
                        help="Max gap between SNPs to link them. (default: 3)")
    parser.add_argument("--max-bp-between-indels", type=int, default=10, metavar="N",
                        help="Max gap between INDELs to link them. (default: 10)")
    args = parser.parse_args()

    from magestic.design.guide_donor_functions import load_genome
    from magestic.design.vcf_processing import (
        load_haplotype_vcf,
        combine_linked_variants,
        filter_variants,
        load_gene_list,
        write_processed_vcf,
    )

    genome_seq = load_genome(args.genome)
    variants_by_sample = load_haplotype_vcf(args.haplotype_vcf)
    combined = combine_linked_variants(
        variants_by_sample, genome_seq,
        max_bp_between_snvs=args.max_bp_between_snvs,
        max_bp_between_indels=args.max_bp_between_indels,
    )

    gene_list = None
    if args.gene_list:
        gene_list = load_gene_list(args.gene_list, max_fdr=args.max_fdr)

    exclude_genes = set(args.exclude_genes) if args.exclude_genes else None
    filtered = filter_variants(
        combined,
        gene_list=gene_list,
        exclude_background=args.exclude_background,
        exclude_gene_names=exclude_genes,
    )

    write_processed_vcf(filtered, args.output)


def design_saturation_main() -> None:
    """CLI entry point: magestic-design-saturation

    Design guide-donor oligo pairs for all amino acid substitutions in a
    target ORF (saturation mutagenesis).
    """
    parser = argparse.ArgumentParser(
        prog="magestic-design-saturation",
        description=(
            "Design saturation mutagenesis guide-donor oligos covering every "
            "amino acid substitution in a target ORF."
        ),
    )
    parser.add_argument("--orf", required=True, nargs="+", metavar="ORF",
                        help="Systematic ORF name(s) to design (e.g. YOR332W PDR5).")
    parser.add_argument("--genome", required=True, metavar="FASTA",
                        help="Reference genome FASTA.")
    parser.add_argument("--gff", required=True, metavar="GFF",
                        help="GFF annotation file for ORF coordinates.")
    parser.add_argument("--nuclease", required=True,
                        choices=["SpCas9", "SpG", "LbCas12a", "impLbCas12a"],
                        help="Nuclease system.")
    parser.add_argument("--allowed-guides", required=True, metavar="TXT",
                        help="Allowed guide sequences from magestic-guide-selection.")
    parser.add_argument("--rev-primers", required=True, metavar="TXT",
                        help="File of reverse subpool priming sequences (one per line).")
    parser.add_argument("--output-dir", required=True, metavar="DIR",
                        help="Output directory for oligo pool TSV(s).")
    parser.add_argument("--output-prefix", default="saturation", metavar="STR",
                        help="Filename prefix for output files. (default: saturation)")
    parser.add_argument("--codon-table", metavar="TXT",
                        help="Codon table file (codon AA, one per line). "
                             "Defaults to bundled S. cerevisiae table.")
    parser.add_argument("--suboptimal-codon-table", metavar="TXT",
                        help="Codon table with suboptimal codons removed. "
                             "Defaults to bundled S. cerevisiae table (frequency > 10%%).")
    parser.add_argument("--subpool-idx", type=int, default=0, metavar="N",
                        help="Starting subpool index into the rev primer list. (default: 0)")
    parser.add_argument("--oligo-length", type=int, default=230, metavar="N",
                        help="Target total oligo length in bases. (default: 230)")
    parser.add_argument("--minimum-homology", type=int, default=20, metavar="N",
                        help="Minimum donor homology arm length. (default: 20)")
    parser.add_argument("--min-disruption-score", type=int, default=4, metavar="N",
                        help="Minimum guide disruption score. (default: 4)")
    parser.add_argument("--subpool-window-size", type=int, default=250, metavar="N",
                        help="Chromosomal window size (bp) for subpool assignment. (default: 250)")
    parser.add_argument("--fold-pseudo-wt-controls", type=int, default=5, metavar="N",
                        help="Copies of synonymous controls to write. (default: 5)")
    parser.add_argument("--no-stop-codons", dest="stop_codons_included",
                        action="store_false", default=True,
                        help="Exclude stop codon as a target mutation.")
    parser.add_argument("--mutate-initiating-methionine", action="store_true",
                        help="Include M1 in the saturation scan.")
    args = parser.parse_args()

    import importlib.resources

    from magestic.design.guide_donor_functions import (
        load_genome,
        load_ORF_coordinates,
        load_codon_table,
    )
    from magestic.design.guide_selection import load_allowed_guides
    from magestic.design.nuclease import get_nuclease
    from magestic.design.saturation_design import design_saturation_oligos

    nuclease = get_nuclease(args.nuclease)
    genome_seq = load_genome(args.genome)
    orf_dict = load_ORF_coordinates(args.gff)

    # --- Codon tables ---
    _data = Path(__file__).parent / "data"
    codon_table_path = args.codon_table or str(_data / "DNA_codon_to_AA.txt")
    suboptimal_table_path = (
        args.suboptimal_codon_table
        or str(_data / "DNA_codon_to_AA_suboptimal_codons_removed.txt")
    )

    codon_to_aa, _ = load_codon_table(codon_table_path)
    _, aa_to_codon_subopt = load_codon_table(suboptimal_table_path)

    # Build aa_to_codon_fraction from suboptimal table (equal weights — used
    # only for highest-frequency-codon subpool assignment)
    aa_to_codon_fraction: dict = {}
    for aa, codons in aa_to_codon_subopt.items():
        n = len(codons)
        aa_to_codon_fraction[aa] = [(1.0 / n, c) for c in codons]

    # suboptimal_removed_aa_to_codon is the aa→[codon] map from suboptimal table
    suboptimal_removed_aa_to_codon = aa_to_codon_subopt

    # --- Allowed guides ---
    allowed_guides = load_allowed_guides(args.allowed_guides)
    logging.getLogger(__name__).info(
        "Loaded %d allowed guides.", len(allowed_guides)
    )

    # --- Reverse primers ---
    rev_primers: list = []
    with open(args.rev_primers) as fh:
        for line in fh:
            rev_primers.append(line.strip().lower())

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for orf_name in args.orf:
        out_file = out_dir / f"{args.output_prefix}_{orf_name}_oligo_pool.tsv"
        logging.getLogger(__name__).info("Designing saturation oligos for %s → %s", orf_name, out_file)
        design_saturation_oligos(
            orf_name=orf_name,
            genome_seq=genome_seq,
            orf_dict=orf_dict,
            codon_to_aa=codon_to_aa,
            aa_to_codon_fraction=aa_to_codon_fraction,
            suboptimal_removed_aa_to_codon=suboptimal_removed_aa_to_codon,
            nuclease=nuclease,
            rev_primer_sequences=rev_primers,
            output_oligo_pool_file=out_file,
            allowed_guides=allowed_guides,
            starting_subpool_idx=args.subpool_idx,
            oligo_length=args.oligo_length,
            minimum_homology=args.minimum_homology,
            min_disruption_score=args.min_disruption_score,
            subpool_window_size=args.subpool_window_size,
            fold_pseudo_wt_controls=args.fold_pseudo_wt_controls,
            stop_codons_included=args.stop_codons_included,
            mutate_initiating_methionine=args.mutate_initiating_methionine,
        )
        print(f"Wrote {out_file}")
