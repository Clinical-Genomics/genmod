#!/usr/bin/env python
# encoding: utf-8
"""
compound_model.py

Checks if the Autosomal Recessive Compound model is followed.


Created by Måns Magnusson on 2013-02-12.
Copyright (c) 2013 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function

import os
import sys
import logging

def check_compounds(variant_1, variant_2, family, intervals, phased):
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
    # Check in all individuals what genotypes that are in the trio based of the
    # individual picked.
    logger = logging.getLogger(__name__)
    
    for individual_id in family.individuals:
        logger.debug("Check compounds for individual {0}".format(individual_id))
        individual = family.individuals[individual_id]
        
        genotype_1 = variant_1['genotypes'][individual_id]
        genotype_2 = variant_2['genotypes'][individual_id]
        
        if individual.has_parents:
            mother_id = individual.mother
            father_id = individual.father
            
            if mother_id != '0':
                # mother_genotype_1 = variant_1['genotypes'][mother_id]
                # mother_genotype_2 = variant_2['genotypes'][mother_id]
                mother = family.individuals[mother_id]
    
            if father_id != '0':
                # father_genotype_1 = variant_1['genotypes'][father_id]
                # father_genotype_2 = variant_2['genotypes'][father_id]
                father = family.individuals[father_id]
            
            if mother_id != '0' and mother.healthy:
                if (variant_1['genotypes'][mother_id].has_variant and 
                    variant_2['genotypes'][mother_id].has_variant):
                    return False
            
            if father_id != '0' and father.healthy:
                if (variant_1['genotypes'][father_id].has_variant and 
                    variant_2['genotypes'][father_id].has_variant):
                    return False
        
        #check if variants are in the same phased interval:
        
        if phased:
            variant_1_interval = intervals[individual_id].find_range(
                                                    [
                                                        int(variant_1['POS']),
                                                        int(variant_1['POS'])
                                                    ]
                                                    )
            variant_2_interval = intervals[individual_id].find_range(
                                                    [
                                                        int(variant_2['POS']),
                                                        int(variant_2['POS'])
                                                    ]
                                                    )
        
        # If phased a healthy individual can have both variants if they are on
        # the same haploblock
        # if not phased:
        
        if individual.healthy:
            if genotype_1.heterozygote and genotype_2.heterozygote:
                return False
        # The case where the individual is affected
        # We know since ealier that all affected are heterozygotes 
        #for these variants
        # So we only need to know if the variants are on the same phase
        elif individual.affected:
            #If the individual is sick and phased it has to have one variant on
            # each allele
            if phased:
                # Variants need to be in the same phased interval, othervise we
                # do not have any extra info
                if variant_1_interval == variant_2_interval:
                # If they are in the same interval they can not be on same 
                #allele
                    if ((genotype_1.allele_1 == genotype_2.allele_1) or 
                        (genotype_1.allele_2 == genotype_2.allele_2)):
                        return False
    
    return True

def main():
    pass

if __name__ == '__main__':
    main()

