#!/usr/bin/env python
# encoding: utf-8
"""
compound_model.py

Checks if the Autosomal Recessive Compound model is followed.


Created by Måns Magnusson on 2013-02-12.
Copyright (c) 2013 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function
from genmod.vcf_tools.genotype import Genotype
from intervaltree import IntervalTree
import logging
from typing import Union

# TODO: What is family?
# TODO: Remove intervals
def check_compounds(variant_1: dict, variant_2: dict, family, intervals: IntervalTree, phased: bool) -> bool:
    """
    Check if two variants of a pair follow the compound heterozygous model.

    At this stage we know:
        - None of the individuals are homozygote alternative for the variants
        - All affected individuals are heterozygote for both variants.

    We do not allow healthy individuals to be heterozygote for both variants
    in the pair
    (ref.
    http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0070151)

    If the individuals are phased we will only consider pairs that are on
    different alleles in affected individuals.

    Args:
        variant_1, variant_2: Variants in a potential compound pair
        family: A family object with the individuals
        intervals: A interval tree that describes the phased intervals
        phased: A bool that tells if the individuals are phased

    Returns:
        bool: depending on if the pair follow the rules stated above

    """
    # Check in all individuals what genotypes that are in the trio based of the individual picked.
    logger = logging.getLogger(__name__)
    for individual_id, individual in family.individuals.items():
        logger.debug("Check compounds for individual {0}".format(individual_id))
        individual = family.individuals[individual_id]

        # If the individual has parents we can check if the parents are healthy and have both variants, if they do they can not be a compound pair, since the variants could be on different alleles.
        # The parents would then be carriers of the disease, which is not possible if the parents are healthy.
        if individual.has_parents:
            for parent_id in (individual.mother, individual.father):
                if parent_id != "0":
                    parent = family.individuals[parent_id]
                    parent_genotypes = [get_genotype(variant, parent) for variant in (variant_1, variant_2)]
                    if parent.healthy and all(genotype.has_variant for genotype in parent_genotypes):
                        return False

        genotype_1 = get_genotype(variant_1, individual_id)
        genotype_2 = get_genotype(variant_2, individual_id)

        # If a healthy individual is not phased and has both variants it can not be a compound pair, since the variants could be on different alleles.
        # The individual would then be a carrier of the disease, which is not possible if the individual is healthy.
        if individual.healthy and (genotype_1.heterozygote and genotype_2.heterozygote):
            return False

        # If the individual is affected and phased we can say that the variants need to be on different alleles, otherwise it can not be a compound pair.
        if individual.affected and phased and variants_on_same_allele(individual.individual_id, variant_1, variant_2):
            return False

    # If the individual is affected and not phased we can not say anything about the phase, so we say it is a compound pair, since it could be a compound pair.
    return True

# TODO: Write a test for this
# TODO: Fix KeyError
def get_genotype(variant: dict, individual_id: str) -> Genotype:
    """
    Return the Genotype object for a variants for a given individual.

    Args:
        variant (dict): A dictionary representing a variant (with per-sample genotypes)
        individual_id (str): Sample/individual ID

    Returns:
        Genotype: The genotype object for the individual at the given variant
    """
    return variant["genotypes"][individual_id]

# TODO: Write a test for this
def get_phase_set(variant_dict: dict, sample_id: str) -> Union[str, None]:
    """
    Extracts the PS (phase set) field for a given sample from a VCF-like variant dictionary.

    Parameters:
        variant_dict (dict): A dictionary representing a variant (with FORMAT and per-sample fields)
        sample_id (str): The sample/individual ID to extract the PS for

    Returns:
        str or None: The PS value for the sample, or None if missing
    """
    format_fields = variant_dict["FORMAT"].split(":")
    sample_values = variant_dict.get(sample_id, "").split(":")
    field_map = dict(zip(format_fields, sample_values))
    return field_map.get("PS")

# TODO: Write a test for this
def variants_on_same_allele(individual_id: str, variant_1: dict, variant_2: dict) -> bool:
    """
    Determine if two phased variants are on the same haplotype for a given individual.

    This function assumes that both variants are phased and that the PS (phase set)
    field is present in the genotype information. Two variants are considered to be
    on the same allele if they belong to the same phase set and have at least one
    allele in the same haplotype position (allele_1 with allele_1 or allele_2 with allele_2).

    Parameters:
        individual_id (str): ID of the individual/sample.
        variant_1 (dict): First variant dictionary containing per-sample genotypes.
        variant_2 (dict): Second variant dictionary containing per-sample genotypes.

    Returns:
        bool: True if the variants are on the same allele in the same phase set,
              False otherwise.
    """
    genotype_1 = variant_1["genotypes"][individual_id]
    genotype_2 = variant_2["genotypes"][individual_id]

    same_phase = get_phase_set(variant_1, individual_id) == get_phase_set(variant_2, individual_id)
    overlapping_allele = (
        genotype_1.allele_1 == genotype_2.allele_1 or
        genotype_1.allele_2 == genotype_2.allele_2
    )

    # Variants are on the same allele if they overlap in the same haplotype and are in the same phase set
    return same_phase and overlapping_allele

def main():
    pass


if __name__ == "__main__":
    main()
