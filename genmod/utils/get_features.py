from __future__ import print_function

import logging

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
}


def check_vep_annotation(variant):
    """
    Return a set with the genes that vep has annotated this variant with.

    Vep annotates all variants but we are only interested in the exonic ones.
    The terms are specified in INTERESTING_SO_TERMS

    Arguments:
        variant (dict): A variant dictionary

    Returns:
        annotation (set): A set with genes
    """

    annotation = set()
    # vep_info is a dictionary with genes as key and annotation as values
    ##TODO use extract_vcf to get the annotation here

    vep_info = variant.get("vep_info", {})

    for allele in vep_info:
        for vep_annotation in variant["vep_info"][allele]:
            for consequence in vep_annotation.get("Consequence", {}).split("&"):
                # These are the SO terms that indicate that the variant
                # belongs to a gene
                if consequence in INTERESTING_SO_TERMS:
                    annotation.add(vep_annotation.get("Gene", ""))
    return annotation


def get_annotation(variant, annotation_key="Annotation", vep=False):
    """
    Return the features that a variant belongs to.

    Arguments:
        variant (dict): A variant dictionary
        annotation_key (str): The name of the info field to search
        vep (bool): If variants are annotated with vep

    Returns:
        annotations (set): A set with annotated features
    """
    logger = logging.getLogger(__name__)
    ##TODO use extract_vcf to get the annotation here

    annotation = set()
    variant_id = variant.get("variant_id", "")
    logger.debug("Checking variant annotation for {0}".format(variant_id))
    # If the variant has already been annotated by genmod we do not need to
    # check again
    if vep:
        logger.debug("Using vep annotation.")
        annotation = check_vep_annotation(variant)

    else:
        info_dict = variant.get("info_dict", {})
        if info_dict.get(annotation_key, None):
            annotation = set(info_dict[annotation_key].split(","))

    logger.debug("Annotations found for {0}: {1}".format(variant_id, ",".join(annotation)))
    return annotation
