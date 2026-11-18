from __future__ import print_function

import logging
from typing import Any, Iterable, NamedTuple, Optional, Sequence

from genmod.vcf_tools import parse_annotations

INTERESTING_SO_TERMS = {
    "transcript_ablation",
    "splice_donor_variant",
    "splice_acceptor_variant",
    "stop_gained",
    "start_lost",
    "frameshift_variant",
    "stop_lost",
    "initiator_codon_variant",
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
    "missense_variant",
    "splice_region_variant",
    "incomplete_terminal_codon_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_exon_variant",
    "non_coding_transcript_exon_variant",
    "nc_transcript_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    # Terms used by SnpEff but not by VEP
    "conservative_inframe_deletion",
    "conservative_inframe_insertion",
    "disruptive_inframe_deletion",
    "disruptive_inframe_insertion",
    "exon_loss_variant",
    "rare_amino_acid_variant",
    "5_prime_UTR_premature_start_codon_gain_variant",
    "start_retained_variant",
    "splice_branch_variant",
}


class ConsequenceSource(NamedTuple):
    """Where a consequence annotation tool stores its per-transcript annotations.

    Attributes:
        info_key: The INFO field with the annotations.
        gene_column: The annotation field with the gene id.
        consequence_column: The annotation field with the ``&``-joined SO terms.
    """

    info_key: str
    gene_column: str
    consequence_column: str


VEP = ConsequenceSource(info_key="CSQ", gene_column="Gene", consequence_column="Consequence")
SNPEFF = ConsequenceSource(info_key="ANN", gene_column="Gene_ID", consequence_column="Annotation")


def get_annotated_genes(
    annotations: Iterable[dict[str, str]], source: ConsequenceSource
) -> set[str]:
    """
    Return a set with the genes that a variant is annotated with.

    VEP and SnpEff annotate all variants but we are only interested in the ones
    that belong to a gene. The terms are specified in INTERESTING_SO_TERMS.

    Arguments:
        annotations (list[dict]): The parsed annotations, see parse_annotations
        source (ConsequenceSource): Which fields hold the gene and the consequence

    Returns:
        annotation (set): A set with genes
    """
    genes = set()
    for annotation in annotations:
        gene = annotation.get(source.gene_column, "")
        if not gene:
            continue
        consequences = annotation.get(source.consequence_column, "").split("&")
        # These are the SO terms that indicate that the variant belongs to a gene
        if INTERESTING_SO_TERMS.intersection(consequences):
            genes.add(gene)
    return genes


def get_annotation(
    variant: dict[str, Any],
    annotation_key: str = "Annotation",
    source: Optional[ConsequenceSource] = None,
    columns: Sequence[str] = (),
) -> set[str]:
    """
    Return the features that a variant belongs to.

    Arguments:
        variant (dict): A variant dictionary
        annotation_key (str): The name of the info field to search
        source (ConsequenceSource): If given, take the genes from the VEP or SnpEff
            annotation instead of from annotation_key
        columns (list[str]): The field names of the source annotation

    Returns:
        annotations (set): A set with annotated features
    """
    logger = logging.getLogger(__name__)

    annotation = set()
    variant_id = variant.get("variant_id", "")
    logger.debug("Checking variant annotation for %s", variant_id)
    info_dict = variant.get("info_dict", {})
    if source:
        logger.debug("Using %s annotation.", source.info_key)
        annotation = get_annotated_genes(
            annotations=parse_annotations(info_dict.get(source.info_key), columns),
            source=source,
        )

    elif info_dict.get(annotation_key, None):
        annotation = set(info_dict[annotation_key].split(","))

    logger.debug("Annotations found for %s: %s", variant_id, ",".join(annotation))
    return annotation
