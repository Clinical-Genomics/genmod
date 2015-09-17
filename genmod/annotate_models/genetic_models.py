#!/usr/bin/env python
# encoding: utf-8
"""
genetic_models.py

Genetic models take a family object with individuals and variants and annotates
for each variant which models they follow in this family.

The following models are checked:

- Autosomal Dominant(AD)
- Autosomal Dominant De Novo(AD_DN)
- Autosomal Recessive(AR_hom)
- Autosomal Recessive De Novo(AR_DN)
- Autosomal Recesive Compound(AR_comp).

In this model a variant must imply affected status, otherwise it can not be 
dominant. All sick has to be ay least heterozygote for the variant and all 
healthy can not have it.

We will assume that each individual that we have information about is present
among the individual in self.family.individuals.


is the individual sick?

    - If the individual is homozygote alternative then AD/AD-denovo 
        and AR/AR-denovo are ok

    - If the individual is is heterozygote then AD/AD-denovo are ok 
        but AR/AR-denovo are not ok

    - If the individual is homozygote reference no model is ok

    - If the individual has no call all models are ok



is the individual healthy?

    - If the individual is homozygote alternative then no model is ok

    - If the individual is heterozygote then AR/AR-denove are ok but 
        AD/AD-denovo are not ok

    - If the individual is homozygote referense all models are ok

    - If there is no call all models are ok



is there no known phenotype?

    - All models are ok for all variants



Created by Måns Magnusson on 2013-02-12.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

from __future__ import print_function, absolute_import

import os
import sys
import logging
from datetime import datetime

from .models import (check_dominant, check_recessive, check_compounds, 
check_X_recessive, check_X_dominant)

from genmod.utils import generate_pairs

def check_genetic_models(variant_batch, families, phased = False, 
                         strict = False):
    """
    Check and annotate which genetic models that are followed for the variants 
    in a batch
    
    See more description in header of file or documentation.
    
    Arguments:
        variant_batch (dict): A dictionary with variant ids as keys and variant
                              dictionaries as values
        families (dict): A dictionary with family ids as keys and Family 
                        objects as values
        phased (bool): If the variants are phased
        strict (bool): If the strict mode should be used when checking the 
                       genetic models
        
    """
    # A variant batch is a dictionary on the form 
    # {variant_id:variant_dict, variant_2_id:variant_dict_2, ...}
    logger = logging.getLogger(__name__)
    intervals = variant_batch.pop('haploblocks', {})
    
    # We check the genetic models for one family at a time
    for family_id in families:
        logger.debug("Checking genetic models for family {0}".format(
            family_id
        ))
        family = families[family_id]
        individuals = family.individuals
        
        compound_candidates = []
        compound_pairs = []
        
        for variant_id in variant_batch:
            inheritance_models = {
                'XR' : False,
                'XR_dn' : False,
                'XD' : False,
                'XD_dn' : False,
                'AD' : False,
                'AD_dn' : False,
                'AR_hom' : False,
                'AR_hom_dn' : False,
                'AR_comp' : False,
                'AR_comp_dn' : False
            }
            
            variant = variant_batch[variant_id]
            # save the compound pairs for a variant in a set
            if 'compounds' in variant:
                variant['compounds'][family_id] = set()
            else:
                variant['compounds'] = {family_id : set()}
            
            # Add information of models followed:
            if 'inheritance_models' in variant:
                variant['inheritance_models'][family_id] = inheritance_models
            else:
                variant['inheritance_models'] = {family_id: inheritance_models}
                
                
            # If the variant is in a genetic region we check for compound 
            # candidates
            if variant.get('compound_candidate',True):
                
                if check_compound_candidate(variant, family, strict):
                    compound_candidates.append(variant_id)
            
            # Only check X-linked for the variants in the X-chromosome:
            # For X-linked we do not need to check the other models
            if variant['CHROM'] == 'X':
                if check_X_recessive(variant, family, strict):
                    variant['inheritance_models'][family_id]['XR'] = True
                    for individual_id in individuals:
                        individual = individuals[individual_id]
                        if individual.has_parents:
                            check_parents(
                                'X_recessive', 
                                individual_id, 
                                family, 
                                variant
                            )
            
                if check_X_dominant(variant, family, strict):
                    variant['inheritance_models'][family_id]['XD'] = True
                    for individual_id in family.individuals:
                        individual = individuals[individual_id]
                        if individual.has_parents:
                            check_parents(
                                'X_dominant', 
                                individual_id, 
                                family, 
                                variant
                            )
            # If variant is not on X:
            else:
                # Check the dominant model:
                if check_dominant(variant, family, strict):
                    variant['inheritance_models'][family_id]['AD'] = True
                    for individual_id in individuals:
                        individual = individuals[individual_id]
                        if individual.has_parents:
                            check_parents(
                                'dominant', 
                                individual_id, 
                                family, 
                                variant
                            )
                
                # Check the recessive model:
                if check_recessive(variant, family, strict):
                    variant['inheritance_models'][family_id]['AR_hom'] = True
                    for individual_id in individuals:
                        individual = individuals[individual_id]
                        if individual.has_parents:
                            check_parents(
                                'recessive', 
                                individual_id, 
                                family, 
                                variant
                            )
            
        # Now check the compound models:
            
        if len(compound_candidates) > 1:
            for pair in generate_pairs(compound_candidates):
            # If the variants in the pair belong to the same gene we check for compounds:
                variant_1 = variant_batch[pair[0]]
                variant_2 = variant_batch[pair[1]]
                # Check that the pair is in the same feature:
                if variant_1['annotation'].intersection(variant_2['annotation']):
                    if len(individuals) == 1:
                        variant_1['compounds'][family_id].add(pair[1])
                        variant_2['compounds'][family_id].add(pair[0])
                        variant_1['inheritance_models'][family_id]['AR_comp'] = True
                        variant_2['inheritance_models'][family_id]['AR_comp'] = True
                # We know from check_compound_candidates that all variants are present in all affected
                    elif check_compounds(variant_1, variant_2, family, intervals, phased):
                        parents_found = False
                        for individual_id in individuals:
                            individual = individuals[individual_id]
                            if individual.has_parents:
                                check_parents('compound', individual_id, family, variant_1, variant_2)
                                parents_found = True
                        if not parents_found:
                            variant_1['inheritance_models'][family_id]['AR_comp'] = True
                            variant_2['inheritance_models'][family_id]['AR_comp'] = True
                                
                        variant_1['compounds'][family_id].add(pair[1])
                        variant_2['compounds'][family_id].add(pair[0])
    return

def check_compound_candidate(variant, family, strict):
    """
    Sort out the variants that are potential compound candidates. 
    This function is used to reduce the number of potential candidates
    for the future analysis. It will go through all variants in a 
    batch(gene or other feature) and filter out those variants that not 
    fit the model. Returns a bool depending on if the variant is a 
    potential compound candidate.
    
    Cases:
        Affected:
            - If individual is affected it needs to be heterozygpte 
            otherwise it can not be a compound candidate
            
        Healthy:
            - Can not be hom. alt for any variant in a potential 
            compound pair.
    
    If strict:
        Affected must be heterozygote
    
    Args:
        variant : A variant dictionary.
        family  : A family object with information about the family 
                    members for this analysis
        
    Returns:
        bool: depending on if the variant is a potential compound 
                candidate according to therules stated above
    
    """
    # This is the case when the variant is located in an uninteresting region(non gene region):
    
    for individual_id in family.individuals:
        individual = family.individuals[individual_id]
        individual_genotype = variant['genotypes'][individual_id]
        
        # No individuals can be homo_alt
        if individual_genotype.homo_alt:
            return False
        
        if individual.affected:
            # Affected have to be heterozygote for compounds
            if not individual_genotype.heterozygote:
                return False
            # If both parents are healthy none of them can have both variants
            
            mother_id = individual.mother
            father_id = individual.father
            
            if mother_id != '0':
                mother_genotype = variant['genotypes'][mother_id]
                mother = family.individuals[mother_id]
                
            if father_id != '0':
                father_genotype = variant['genotypes'][father_id]
                father = family.individuals[father_id]
            
            # If both parents exist and both are healthy, 
            # both can not have the variant
            if mother_id != '0' and father_id != '0':
                if mother.healthy and father.healthy:
                    if (mother_genotype.has_variant and 
                            father_genotype.has_variant):
                        return False
            # We have now significantly reduced the number
            # of compound candidates.
            # In the next step we check if pairs of compounds
            # follow the compound inheritance pattern.
    
    return True


def check_parents(model, individual_id, family, variant, variant_2={}, 
                    strict = False):
    """
    Check if information in the parents can tell us if model is 
    de novo or not. 
    Model IN ['recessive', 'compound', 'dominant', 'X_recessive', 'X_dominant'].
    If the expected pattern of a variant is followed in the family, 
    de novo will be False. Otherwise de novo will be True.
    
    
    If only one parent is present then we can never exclude denovo for 
    heterozygous inheritance patterns.
    If strict and one parent we will never say it is denovo
    
    Args:
        model  : String, one of 'recessive', 'compound', 'dominant', 
                    'X_recessive', 'X_dominant'
        individual_id   : String that represents the individual id
        family  : A family object
        variant : A dictionary that represents the variant
        variant_2   : If compound pair this is the second variant
        strict  : Bool
    
    """
    sex = family.individuals[individual_id].sex
    family_id = family.family_id
    
    mother_genotype = False
    father_genotype = False
    
    parent_genotypes = []
    mother_id = family.individuals[individual_id].mother
    father_id = family.individuals[individual_id].father
    
    if mother_id != '0':
        mother_genotype = variant['genotypes'][mother_id]
        mother_phenotype = family.get_phenotype(mother_id)
        parent_genotypes.append(mother_genotype)
    
    if father_id != '0':
        father_genotype = variant['genotypes'][father_id]
        father_phenotype = family.get_phenotype(father_id)
        parent_genotypes.append(father_genotype)

    if model == 'recessive':
        # If a parent is homozygote or if both parents are heterozygote 
        # the variant is not denovo.
        # If strict we know from before that both parents are genotyped
        if len(parent_genotypes) == 2:
            if not (mother_genotype.has_variant and 
                    father_genotype.has_variant):
                variant['inheritance_models'][family_id]['AR_hom_dn'] = True
            # If both parents are called but none of the above is 
            # fullfilled it is pure denovo
                if (mother_genotype.genotyped and father_genotype.genotyped):
                    variant['inheritance_models'][family_id]['AR_hom'] = False
        elif not strict:
            variant['inheritance_models'][family_id]['AR_hom_dn'] = True
                    
    elif model == 'dominant':
        # If none of the parents carry variant it is de novo
        if len(parent_genotypes) == 2:
            if not (mother_genotype.has_variant or father_genotype.has_variant):
                variant['inheritance_models'][family_id]['AD_dn'] = True
        # If both parents are called but none of them carry the variant it is denovo
                if mother_genotype.genotyped and father_genotype.genotyped:
                    variant['inheritance_models'][family_id]['AD'] = False
        else:
            for parent in parent_genotypes:
                if not parent.has_variant:
                    variant['inheritance_models'][family_id]['AD_dn'] = True
                    variant['inheritance_models'][family_id]['AD'] = False
            
    elif model == 'X_recessive':
        #If the individual is a male we only need to check if the mother carry the variant:
        if sex == 1:
            if mother_genotype:
                if not mother_genotype.has_variant:
                    variant['inheritance_models'][family_id]['XR_dn'] = True
                    if mother_genotype.genotyped:
                        variant['inheritance_models'][family_id]['XR'] = False
            elif not strict:
                variant['inheritance_models'][family_id]['XR_dn'] = True
                
        #If female, both parents must have the variant otherwise denovo is true
        elif sex == 2:
            if len(parent_genotypes) == 2:
                if not (mother_genotype.has_variant and father_genotype.has_variant):
                    variant['inheritance_models'][family_id]['XR_dn'] = True
        #If both parents are genotyped but they both are not carriers XR is not true
                    if (mother_genotype.genotyped and father_genotype.genotyped):
                        variant['inheritance_models'][family_id]['XR'] = False
            elif not strict:
                variant['inheritance_models'][family_id]['XR_dn'] = True
                
    elif model == 'X_dominant':
        #If the individual is a male we only need to look at the mother:
        if sex == 1:
            if mother_genotype:
                if not mother_genotype.has_variant:
                    variant['inheritance_models'][family_id]['XD_dn'] = True
                    if mother_genotype.genotyped:
                        variant['inheritance_models'][family_id]['XD'] = False
        #If female, one of the parents must have the variant otherwise denovo is true
        elif sex == 2:
            if len(parent_genotypes) == 2:
                if not (mother_genotype.has_variant or father_genotype.has_variant):
                    variant['inheritance_models'][family_id]['XD_dn'] = True
                    if mother_genotype.genotyped and father_genotype.genotyped:
                        variant['inheritance_models'][family_id]['XD'] = False
            elif not strict:
                variant['inheritance_models'][family_id]['XD_dn'] = True
    
    elif model == 'compound':
        
        mother_genotype_2 = None
        father_genotype_2 = None
        parent_genotypes_2 = []
        
        if mother_id != '0':
            mother_genotype_2 = variant_2['genotypes'][mother_id]
            parent_genotypes_2.append(mother_genotype_2)
        if father_id != '0':
            father_genotype_2 = variant_2['genotypes'][father_id]
            parent_genotypes_2.append(father_genotype_2)
        # One of the variants must come from father and one from mother
        if (len(parent_genotypes) == 2 and len(parent_genotypes_2) == 2):
            # If both parents are genotyped and one of them are homozygote reference for both variants
            # the pair will be considered AR compound de novo
            if ((mother_genotype.genotyped and mother_genotype_2.genotyped) and
                father_genotype.genotyped and father_genotype_2.genotyped):
                
                # if not both parents have one of the variants it is de novo
                if not ((mother_genotype.has_variant or mother_genotype_2.has_variant) and 
                        (father_genotype.has_variant or father_genotype_2.has_variant)):
                    variant['inheritance_models'][family_id]['AR_comp_dn'] = True
                    variant_2['inheritance_models'][family_id]['AR_comp_dn'] = True
                
                else:
                    
                    variant['inheritance_models'][family_id]['AR_comp'] = True
                    variant_2['inheritance_models'][family_id]['AR_comp'] = True
            
        elif not strict:
            variant['inheritance_models'][family_id]['AR_comp_dn'] = True
            variant_2['inheritance_models'][family_id]['AR_comp_dn'] = True
            variant['inheritance_models'][family_id]['AR_comp'] = True
            variant_2['inheritance_models'][family_id]['AR_comp'] = True
            
    return
            
        

def main():
    from ped_parser import family, individual
    from interval_tree import interval_tree
    from vcf_parser import genotype
    
    duo_family = family.Family(family_id = '1')
    sick_son = individual.Individual(ind='1', family='1',mother='3', father='0', sex=1, phenotype=2)
    healthy_mother = individual.Individual(ind='3', family='1',mother='0', father='0', sex=2, phenotype=1)
    duo_family.add_individual(sick_son)
    duo_family.add_individual(healthy_mother)
    
    pp(duo_family.individuals)
    intervals = {ind_id:interval_tree.IntervalTree([[1,100, '1']], 1, 100) for ind_id in duo_family.individuals}

    pp(intervals)
    
    #Setup two variants with autosomal recessive compound pattern
    recessive_comp_simple_1 = {'CHROM':'1', 'POS':'5', 'ALT':'A', 'REF':'C', 'ID':'rs2230749',
                             '1':'0|1', '3':'0|0'}
    genotypes = {'1':genotype.Genotype('0|1'), '3':genotype.Genotype('0|0')}
    recessive_comp_simple_1['genotypes'] = genotypes
    
    recessive_comp_simple_2 = {'CHROM':'1', 'POS':'10', 'ALT':'C', 'REF':'T', 'ID':'.', 
                            '1':'1|0', '3':'0|1'}
    genotypes = {'1':genotype.Genotype('1|0'), '3':genotype.Genotype('0|1')}
    recessive_comp_simple_2['genotypes'] = genotypes
    
    
    batch = {'ABC':{'1_5_A_C':recessive_comp_simple_1, '1_10_C_T':recessive_comp_simple_2}}
    batch['haploblocks'] = intervals
    
    check_genetic_models(batch, duo_family, phased=True)
    for gene in batch:
        # pp(batch[gene])
        for variant in batch[gene]:
            pp(batch[gene][variant])
            for ind in batch[gene][variant]['genotypes']:
                print(ind, batch[gene][variant]['genotypes'][ind].__dict__)
    

if __name__ == '__main__':
    main()

