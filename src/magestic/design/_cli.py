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


def _load_rev_primers(path: str | Path) -> list[str]:
    """Read subpool reverse-priming sequences and return them **reverse-complemented**.

    🔴 The file stores each subpool's reverse primer in its own 5'->3' orientation, but an
    oligo carries that primer's **reverse complement** at its 3' end, because the primer
    anneals to the bottom strand. ``design_oligos_from_vcf`` and
    ``design_saturation_oligos`` both simply concatenate what they are given, so the
    orientation has to be fixed here, at the read.

    Verified against the real ordered NRD1 pool
    (``NRD1/NRD1_SpCas9_NGG_guide_donor_200bp_oligo_pool.tsv``): all **2,308** oligos end
    in ``CGATAGACGACTGGACAGCA``, which is ``rev_comp`` of line 1 of
    ``rev_subpool_priming_sequences.txt`` (``TGCTGTCCAGTCGTCTATCG``) and matches no line
    as written. The original generator did this at load time
    (``NRD1/generate_guide_donors_for_aa_substitutions_Cas9.py:161``); the port to this CLI
    dropped it, so every emitted subpool would have been unamplifiable with its intended
    primer. Caught by the EMBL/TRIUMPH design review, 2026-09-30, where a collaborator had
    independently reverse-complemented the primers by hand.
    """
    from magestic.design.guide_donor_functions import rev_comp

    primers: list[str] = []
    with open(path) as fh:
        for line in fh:
            seq = line.strip()
            if seq:
                primers.append(rev_comp(seq.lower()))
    return primers


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

    rev_primers = _load_rev_primers(args.rev_primers)

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
    parser.add_argument("--exclude-background", action=argparse.BooleanOptionalAction,
                        default=True,
                        help="Exclude background_variant ORF entries, i.e. variants that "
                             "distinguish the background strain from S288c and are therefore not "
                             "edits. Pass --no-exclude-background to keep them; the 2024-04-11 "
                             "Twist order deliberately included 479 such oligos. "
                             "(default: exclude)")
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

    # --- Reverse primers (reverse-complemented at read; see _load_rev_primers) ---
    rev_primers: list = _load_rev_primers(args.rev_primers)

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


def strain_background_main() -> None:
    """``magestic-strain-background`` — build a strain's design reference from a VCF.

    Packages the step that was previously project-specific
    (``get_haplotype_blocks_from_vcf_with_coord_offset.py``): substitute one strain's
    variants into a reference, emit the coordinate-offset table, remap the annotation,
    and optionally project other strains' alleles onto the new coordinates so they can
    be designed as edits.
    """
    from .strain_background import (
        annotate_variants_with_features,
        build_strain_reference,
        load_strain_variants,
        project_variants,
        read_fasta,
        remap_gff,
        write_fasta,
    )

    parser = argparse.ArgumentParser(
        description="Build a de novo background-strain reference and coordinate map.")
    parser.add_argument("--vcf", required=True, metavar="VCF",
                        help="Multi-sample VCF of strain variants (may be .gz).")
    parser.add_argument("--recipient", required=True, metavar="SAMPLE",
                        help="Sample whose variants define the new reference (the strain "
                             "that will be transformed).")
    parser.add_argument("--genome", required=True, metavar="FASTA",
                        help="Starting reference FASTA, e.g. S288C R64 .chrnames.fsa. "
                             "Every contig is retained, including chrM and the 2-micron.")
    parser.add_argument("--gff", metavar="GFF",
                        help="Annotation to remap onto the new coordinates.")
    parser.add_argument("--donors", nargs="*", default=[], metavar="SAMPLE",
                        help="Other samples whose alleles should be projected onto the "
                             "recipient reference as candidate edits. For an all-by-all "
                             "swap, pass the other N-1 strains.")
    parser.add_argument("--output-dir", required=True, metavar="DIR")
    parser.add_argument("--min-qual", type=float, default=None, metavar="Q",
                        help="Drop records with QUAL below Q. Default: no QUAL filter.")
    parser.add_argument("--keep-filtered", action="store_true",
                        help="Keep records whose FILTER is not PASS. Default: drop them.")
    parser.add_argument("--allow-heterozygous", action="store_true",
                        help="Apply the non-reference allele of a heterozygous call. Default: "
                             "skip het calls, because writing one allele of a het site into a "
                             "design reference asserts more than the data supports.")
    parser.add_argument("--flank", type=int, default=2000, metavar="N",
                        help="Distance from a gene within which a non-coding variant is "
                             "labeled upstream:/downstream: rather than intergenic. "
                             "Default 2000.")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    reference = read_fasta(args.genome)
    print(f"Reference: {len(reference)} contigs, "
          f"{sum(len(s) for s in reference.values()):,} bp")

    recipient_vars, rejected = load_strain_variants(
        args.vcf, args.recipient,
        require_pass=not args.keep_filtered,
        min_qual=args.min_qual,
        require_homozygous=not args.allow_heterozygous,
        chroms=set(reference),
    )
    print(f"{args.recipient}: {len(recipient_vars):,} variants loaded, "
          f"{len(rejected):,} rejected at load")

    strain_seqs, cmap, build_rejected = build_strain_reference(reference, recipient_vars)
    all_rejected = rejected + build_rejected
    print(f"{args.recipient}: {len(build_rejected):,} rejected during application")

    ref_path = out_dir / f"{args.recipient}_reference.fasta"
    write_fasta(strain_seqs, ref_path)
    print(f"Wrote {ref_path}")

    offsets_path = out_dir / f"{args.recipient}_coordinate_offsets.tsv"
    with open(offsets_path, "w") as fh:
        fh.write("CHROM\tREF_LENGTH\tSTRAIN_LENGTH\tNET_LENGTH_CHANGE\n")
        for chrom in reference:
            fh.write(f"{chrom}\t{len(reference[chrom])}\t{len(strain_seqs[chrom])}"
                     f"\t{cmap.total_length_change(chrom)}\n")
    print(f"Wrote {offsets_path}")

    # Rejections are written, never just counted: a silent drop is the failure mode
    # this whole module exists to prevent.
    rej_path = out_dir / f"{args.recipient}_rejected_variants.tsv"
    with open(rej_path, "w") as fh:
        fh.write("CHROM\tPOS\tREF\tALT\tGENOTYPE\tREASON\n")
        for r in all_rejected:
            v = r.variant
            fh.write(f"{v.chrom}\t{v.pos}\t{v.ref}\t{v.alt}\t{v.genotype}\t{r.reason}\n")
    print(f"Wrote {rej_path} ({len(all_rejected):,} rows)")

    if args.gff:
        gff_out = out_dir / f"{args.recipient}_annotation.gff"
        counts = remap_gff(args.gff, cmap, gff_out)
        print(f"Wrote {gff_out} — {counts['written']:,} features, "
              f"{counts['spans_indel']:,} span an indel, {counts['dropped']:,} dropped")
        if counts["spans_indel"]:
            print(f"  NOTE: {counts['spans_indel']:,} features change length in this strain. "
                  f"CDS-relative indexing for those genes must use the coordinate map, "
                  f"not the remapped endpoints.")

    for donor in args.donors:
        donor_vars, donor_rej = load_strain_variants(
            args.vcf, donor,
            require_pass=not args.keep_filtered,
            min_qual=args.min_qual,
            require_homozygous=not args.allow_heterozygous,
            chroms=set(reference),
        )
        projected, proj_rej = project_variants(donor_vars, cmap, strain_seqs)
        if args.gff:
            projected = annotate_variants_with_features(
                projected, out_dir / f"{args.recipient}_annotation.gff", flank=args.flank)

        edits_path = out_dir / f"edits_{donor}_into_{args.recipient}.tsv"
        with open(edits_path, "w") as fh:
            fh.write("CHROM\tPOS\tREF\tALT\tSAMPLE\tORF\tNAME\n")
            for v in projected:
                fh.write(f"{v.chrom}\t{v.pos}\t{v.ref}\t{v.alt}\t{v.sample}"
                         f"\t{v.orf}\t{v.name}\n")
        n_noop = sum(1 for r in proj_rej if "already carries" in r.reason)
        print(f"Wrote {edits_path} — {len(projected):,} edits "
              f"({n_noop:,} dropped as already present in {args.recipient}, "
              f"{len(donor_rej):,} rejected at load)")


def verify_designs_main() -> None:
    """``magestic-verify-designs`` — check a design table against a de novo assembly.

    The recommended guard whenever a design reference was built by substituting variants
    into a distant reference: a guide or homology arm built on locally wrong sequence
    fails at the bench in a way that looks exactly like a biological negative.
    """
    import csv

    from .strain_background import read_fasta, verify_designs_against_assembly

    parser = argparse.ArgumentParser(
        description="Verify designed guides and homology arms occur verbatim in an assembly.")
    parser.add_argument("--designs", required=True, metavar="TSV")
    parser.add_argument("--assembly", required=True, metavar="FASTA",
                        help="De novo assembly for the strain the designs target.")
    parser.add_argument("--passed-out", required=True, metavar="TSV")
    parser.add_argument("--failed-out", required=True, metavar="TSV")
    parser.add_argument("--guide-column", default="guide", metavar="COL")
    parser.add_argument("--left-arm-column", default="left_homology_arm", metavar="COL")
    parser.add_argument("--right-arm-column", default="right_homology_arm", metavar="COL")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    with open(args.designs) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        print(f"No rows in {args.designs}", file=sys.stderr)
        raise SystemExit(1)

    assembly = read_fasta(args.assembly)
    passed, failed = verify_designs_against_assembly(
        rows, assembly,
        guide_key=args.guide_column,
        left_arm_key=args.left_arm_column,
        right_arm_key=args.right_arm_column,
    )

    for path, recs in ((args.passed_out, passed), (args.failed_out, failed)):
        if not recs:
            Path(path).write_text("")
            continue
        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(recs[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(recs)

    total = len(passed) + len(failed)
    pct = 100.0 * len(passed) / total if total else 0.0
    print(f"{len(passed):,} / {total:,} designs verified against the assembly ({pct:.1f}%)")
    print(f"Wrote {args.passed_out} and {args.failed_out}")
