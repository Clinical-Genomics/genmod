#!/usr/bin/env python
# encoding: utf-8
"""
genetic_models.py

Genetic models take a family object with individuals and variants and annotates for each variant which models they follow in this family.

The following models are checked:

- Autosomal Dominant(AD)
- Autosomal Dominant De Novo(AD_DN)
- Autosomal Recessive(AR_hom)
- Autosomal Recessive De Novo(AR_DN)
- Autosomal Recesive Compound(AR_comp).

In this model a variant must imply affected status, otherwise it can not be dominant. All sick has to be ay least heterozygote for the variant and all healthy can not have it.

We will assume that each individual that we have information about is present among the individual in self.family.individuals.


is the individual sick?

    - If the individual is homozygote alternative then AD/AD-denovo and AR/AR-denovo are ok

    - If the individual is is heterozygote then AD/AD-denovo are ok but AR/AR-denovo are not ok

    - If the individual is homozygote reference no model is ok

    - If the individual has no call all models are ok



is the individual healthy?

    - If the individual is homozygote alternative then no model is ok

    - If the individual is heterozygote then AR/AR-denove are ok but AD/AD-denovo are not ok

    - If the individual is homozygote referense all models are ok

    - If there is no call all models are ok



is there no known phenotype?

    - All models are ok for all variants



Created by Måns Magnusson on 2013-02-12.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import os
import sys
from datetime import datetime
from pprint import pprint as pp

from genmod.variants import genotype
from genmod.utils import pair_generator

def check_genetic_models(variant_batch, family, verbose = False, phased = False, strict = False, proc_name = None):
    #A variant batch is a dictionary on the form {gene_id: {variant_id:variant_dict}}
    # Start by getting the genotypes for each variant:
    individuals = family.individuals.values()
    intervals = variant_batch.pop('haploblocks', {})
    for gene in variant_batch:
        
        compound_candidates = []
        compound_pairs = []
        for variant_id in variant_batch[gene]:
            genotypes = {}
            variant = variant_batch[gene][variant_id]
            for individual in family.individuals:
                try:
                    # print(variant_batch[gene][variant_id])
                    gt_info = variant[individual].split(':')[0]
                except KeyError:# If individual not in variant file
                    if verbose:
                        print('Warning! Individual %s is not in variant file!' % individual)
                    gt_info = './.'
                
                individual_genotype = genotype.Genotype(GT=gt_info)
                genotypes[individual] = individual_genotype
            variant['Genotypes'] = genotypes
            variant['Compounds'] = {}
            # Add information of models followed:
            variant['Inheritance_model'] = {'XR' : False, 'XR_dn' : False, 'XD' : False, 
                                            'XD_dn' : False, 'AD' : False, 'AD_dn' : False, 
                                            'AR_hom' : False, 'AR_hom_dn' : False, 'AR_comp' : False, 
                                            'AR_comp_dn' : False
                                            }
            if gene != '-':
                if check_compound_candidates(variant, family):
                    compound_candidates.append(variant_id)
            # Only check X-linked for the variants in the X-chromosome:
            # For X-linked we do not need to check the other models
            if variant['CHROM'] == 'X':
                
                if check_X_recessive(variant, family, strict):
                    variant['Inheritance_model']['XR'] = True
                    for individual in family.individuals:
                        if family.individuals[individual].has_parents:
                            check_parents('X_recessive', individual, family, variant)
                                
                if check_X_dominant(variant, family, strict):
                    variant['Inheritance_model']['XD'] = True
                    for individual in family.individuals:
                        if family.individuals[individual].has_parents:
                            check_parents('X_dominant', individual, family, variant)
            else:
                # Check the dominant model:
                if check_dominant(variant, family, strict):
                    variant['Inheritance_model']['AD'] = True
                    for individual in family.individuals:
                        if family.individuals[individual].has_parents:
                            print('Dominant', variant['variant_id'], individual)
                            check_parents('dominant', individual, family, variant)
                    
                # Check the recessive model:
                if check_recessive(variant, family, strict):
                    variant['Inheritance_model']['AR_hom'] = True
                    for individual in family.individuals:
                        if family.individuals[individual].has_parents:
                            check_parents('recessive', individual, family, variant)
                    
                
            
    
        # Now check the compound models:
        if len(compound_candidates) > 1:
            
            compound_pairs = pair_generator.Pair_Generator(compound_candidates)
            
            for pair in compound_pairs.generate_pairs():
                variant_1 = variant_batch[gene][pair[0]]
                variant_2 = variant_batch[gene][pair[1]]
                # We know from check_compound_candidates that all variants are present in all affected
                
                if check_compounds(variant_1, variant_2, family, intervals, phased):
                    for individual in family.individuals:
                        if family.individuals[individual].has_parents:
                            check_parents('compound', individual, family, variant_1, variant_2)
                    variant_1['Compounds'][pair[1]] = 0
                    variant_2['Compounds'][pair[0]] = 0
                    variant_1['Inheritance_model']['AR_comp'] = True
                    variant_2['Inheritance_model']['AR_comp'] = True
                    
    return

def check_compound_candidates(variant, family):
    """Sort out the variants that are potential compound candidates. 
        This function is used to reduce the number of potential candidates for the future analysis.
        It will go through all variants in a batch(gene or other feature) and filter out those variants that not fit the model.
        Returns a bool depending on if the variant is a potential compound candidate.
        
        Cases:
            Affected:
                - If individual is affected it needs to be heterozygpte otherwise it can not be a compound candidate
                
            Healthy:
                - Can not be hom. alt for any variant in a potential compound pair.
                
        Args:
            variant: A variant dictionary.
            family: A family object with information about the family members for this analysis
            
        Returns:
            bool: depending on if the variant is a potential compound candidate according to the
            rules stated above
        
    """
    # This is the case when the variant is located in an uninteresting region:
    if not variant.get('comp_candidate',True):
        return False
    
    for individual in family.individuals:
        individual_genotype = variant['Genotypes'].get(individual, genotype.Genotype())
        # No individuals can be homo_alt
        if individual_genotype.homo_alt:
            return False
        if family.individuals[individual].affected:
        # Affected can not be hom. ref. for compounds
            if individual_genotype.homo_ref:
                return False
    
    return True

def check_compounds(variant_1, variant_2, family, intervals, phased):
    """Check if two variants of a pair follow the compound heterozygous model. 
        At this stage we know: 
            - None of the individuals are homozygote alternative for the variants
            - All affected individuals are heterozygote for both variants.
            
        We do not allow healthy individuals to be heterozygote for both variants in the pair(ref. 
        http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0070151).
        
        If the individuals are phased we will only consider pairs that are on different alleles in affected individuals.
        
        Args:
            variant_1, variant_2: Variants in a potential compound pair
            family: A family object with the individuals
            intervals: A interval tree that describes the phased intervals
            phased: A bool that tells if the individuals are phased
        
        Returns:
            bool: depending on if the pair follow the rules stated above
            
    """
    # Check in all individuals what genotypes that are in the trio based of the individual picked.
    
    for individual in family.individuals:

        genotype_1 = variant_1['Genotypes'].get(individual, genotype.Genotype())
        genotype_2 = variant_2['Genotypes'].get(individual, genotype.Genotype())
        #check if variants are in the same phased interval:
        if phased:
            variant_1_interval = intervals[individual].find_range([int(variant_1['POS']),int(variant_1['POS'])])
            variant_2_interval = intervals[individual].find_range([int(variant_1['POS']),int(variant_1['POS'])])
        
        if family.individuals[individual].healthy:
            if genotype_1.heterozygote and genotype_2.heterozygote:
                return False
        # If the family is phased we need to check if a healthy individual have both variants on same allele
        # This individual intervals can not overlap
        # If the variants are not in the same phased interval we can not say that the model is not followed
            # if phased:
            #     if variant_1_interval == variant_2_interval:
            #         # If the variants are on different alleles it can not be a compound pair:
            #         if genotype_1.has_variant or genotype_2.has_variant:
            #             if genotype_1.allele_1 != '0':
            #                 if genotype_2.allele_2 != '0':
            #                     return False
            #                 if genotype_1.allele_2 != '0':
            #                     if genotype_2.allele_1 != '0':
            #                         return False        
            # The case where the individual is affected
        elif family.individuals[individual].affected:
            #If the individual is sick and phased it has to have one variant on each allele
            if phased:
                if variant_1_interval == variant_2_interval:
                    if (genotype_1.allele_1 == genotype_2.allele_1) or (genotype_1.allele_2 == genotype_2.allele_2):
                        return False
    
    return True


def check_dominant(variant, family, strict):
    """Check if the variant follows the autosomal dominant (AD) pattern in this family.
        A variant is following the dominant patttern if:
        Healthy:
            - Can not have the variant in any form.
            - If no call we can not exclude dominant.
            if strict:
                - Have to be homozygote reference
                - No call will return false
        
        Affected:
            - Has to be heterozygote for this position.
            - If no call we can not exclude dominant.
            if strict:
                - Have to be heterozygote
                - No call will return false
        
        No affection status:
            We can not tell if variant follows the model or not.
        
        Args:
            variant: variant dictionary.
            family: A family object with the individuals
            strict: A boolean that tells if strict analyzis should be performed.
        
        Return:
            bool: depending on if the model is followed in these indivduals

            """
    
    for individual in family.individuals: 
        # Check in all individuals what genotypes that are in the trio based of the individual picked.
        individual_genotype = variant['Genotypes'].get(individual, genotype.Genotype())
        if family.individuals[individual].healthy:# The case where the individual is healthy
            if individual_genotype.genotyped:
                if not individual_genotype.homo_ref:
                    return False
            elif strict:
                return False
        
        elif family.individuals[individual].affected:
            # The case when the individual is sick
            if individual_genotype.genotyped:
                if not individual_genotype.heterozygote:
                    return False
            elif strict:
                return False
    
    return True

def check_recessive(variant, family, strict):
    """Check if the variant follows the autosomal recessive homozygote (AR_hom) pattern in this family.
        A variant is following the AR_hom pattern if:
        Healthy:
            - Can not be homozygote alternative.
            - If no call we can not exclude AR
            if strict:
                - Have to be homozygote reference or heterozygote.
                - No call will return False
        
        Affected:
            - Have to be homozygote alternative.
            - If no call we can not exclude AR
            if strict:
                - Have to be homozygote alternative
                - No call will return false
        
        No affection status:
            We can not tell if variant follows the model or not.
        
        Args:
            variant: variant dictionary.
            family: A family object with the individuals
            strict: A boolean that tells if strict analyzis should be performed.
        
        Return:
            bool: depending on if the model is followed in these indivduals
        
            """
    for individual in family.individuals:
        individual_genotype = variant['Genotypes'].get(individual, genotype.Genotype())
        # The case where the individual is healthy:
        if family.individuals[individual].healthy:
        # If the individual is healthy and homozygote alt the model is broken.
            if individual_genotype.genotyped:
                if individual_genotype.homo_alt:
                    return False
            elif strict:
                return False
                
        # The case when the individual is sick:
        elif family.individuals[individual].affected:
        # In the case of a sick individual it must be homozygote alternative for Autosomal recessive to be true.
        # Also, we can not exclude the model if no call.
            if individual_genotype.genotyped:
                if not individual_genotype.homo_alt:
                    return False
            elif strict:
                return False
    
    return True

def check_X_recessive(variant, family, strict):
    """Check if the variant follows the x linked heterozygous (XR) pattern of inheritance in this family.
        A variant is following the XR pattern if:
        
        Healthy:
            - Can not be homozygote alternative
            - If no call we can not exclude XR
            - Males can not have variant at all. This is added since sometimes males 
                get called as heterozygotes but this should not be possible since 
                they only have one copy of the X chromosome.
            if strict:
                - Have to be homozygote reference or heterozygote.
                - No call will return False
        
        Affected:
            - Have to be homozygote alternative.
            - If no call we can not exclude AR
            if strict:
                - Have to be homozygote alternative
                - No call will return false
        
        No affection status:
                We can not tell if variant follows the model or not.
        
        Args:
            variant: variant dictionary.
            family: A family object with the individuals
            strict: A boolean that tells if strict analyzis should be performed.
        
        Return:
            bool: depending on if the model is followed in these indivduals
        
            """
    for individual in family.individuals:
        # Get the genotype for this variant for this individual
        individual_genotype = variant['Genotypes'].get(individual, genotype.Genotype())

        # The case where the individual is healthy
        if family.individuals[individual].healthy:
            # If individual is healthy and homozygote alternative the variant can not be deleterious:
            if individual_genotype.genotyped:
                if individual_genotype.homo_alt:
                    return False
                if family.individuals[individual].sex == 1:
                    if individual_genotype.has_variant:
                        return False
            elif strict:
                return False
        
        # The case when the individual is sick
        elif family.individuals[individual].affected:
        #If the individual is sick and homozygote ref it can not be x-recessive
            if individual_genotype.genotyped:
                if individual_genotype.homo_ref:
                    return False
        # Women have to be hom alt to be sick (almost allways carriers)
                elif family.individuals[individual].sex == 2:
                    if not individual_genotype.homo_alt:
                        return False
            elif strict:
                return False
    return True

def check_X_dominant(variant, family, strict):
    """Check if the variant follows the x linked dominant (XD) pattern of inheritance in this family.
        A variant is following the XD pattern if:
        
        Healthy:
            - Can not be homozygote alternative
            - Healthy females can be heterozygotes. This is possible since there are several
             documented diseases where only one allele is expressed at a time during development.
            - If no call we can not exclude XR
            if strict:
                - Have to be homozygote reference (or heterozygote womens).
                - No call will return False
        
        Affected:
            - Have to be heterozygote.
            - If no call we can not exclude AR
            if strict:
                - Have to be heterozygote
                - No call will return false
        
        No affection status:
                We can not tell if variant follows the model or not.
        
        Args:
            variant: variant dictionary.
            family: A family object with the individuals
            strict: A boolean that tells if strict analyzis should be performed.
        
        Return:
            bool: depending on if the model is followed in these indivduals
        
        """
    for individual in family.individuals:
        # Get the genotype for this variant for this individual
        individual_genotype = variant['Genotypes'].get(individual, genotype.Genotype())
        # The case where the individual is healthy
        if family.individuals[individual].healthy:
        # Healthy womans can be carriers but not homozygote:
            if individual_genotype.genotyped:
                if family.individuals[individual].sex == 2:
                    if individual_genotype.homo_alt:
                        return False
                # Males can not carry the variant:
                elif family.individuals[individual].sex == 1:
                    if individual_genotype.has_variant:
                        return False
            elif strict:
                return False
        # The case when the individual is sick
        elif family.individuals[individual].affected:
        #If the individual is sick and homozygote ref it can not be x-linked-dominant
            if individual_genotype.genotyped:
                if individual_genotype.homo_ref:
                    return False
            elif strict:
                return False
    return True

def check_parents(model, individual, family, variant, variant_2={}, strict = False):
    """Check if information in the parents can tell us if model is de novo or not. 
        Model in ['recessive', 'compound', 'dominant'].
        This function is to check if variant should be annotated with de novo or not. If the expected pattern of a variant is
        followed in the family, de novo will be False. Otherwise de novo will be True.
         
         """
    sex = family.individuals[individual].sex
    individual_genotype = variant['Genotypes'].get(individual, genotype.Genotype())

    mother_id = family.individuals[individual].mother
    # If mother does not exist we pretend that she allways have the variant
    if mother_id == '0':
        mother_genotype = variant['Genotypes'].get(mother_id, genotype.Genotype('0/1'))
    
    mother_genotype = variant['Genotypes'].get(mother_id, genotype.Genotype())
    mother_phenotype = family.get_phenotype(mother_id)

    father_id = family.individuals[individual].father
    if father_id == '0':
        father_genotype = variant['Genotypes'].get(father_id, genotype.Genotype('0/1'))

    father_genotype = variant['Genotypes'].get(father_id, genotype.Genotype())
    father_phenotype = family.get_phenotype(father_id)


    if model == 'recessive':
        # If a parnent is homozygote or if both parents are heterozygote the variant is not denovo
        if not (mother_genotype.has_variant and father_genotype.has_variant):
            variant['Inheritance_model']['AR_hom_dn'] = True
            # If both parents are called but none of the above is fullfilled it is pure denovo
            if mother_genotype.genotyped and father_genotype.genotyped:
                variant['Inheritance_model']['AR_hom'] = False
                    
    elif model == 'dominant':
        # If one or both parents are affected it is de novo if none of them have a variant
        if not (mother_genotype.has_variant or father_genotype.has_variant):
            variant['Inheritance_model']['AD_dn'] = True
        # If both parents are called but none of them carry the variant it is denovo
            if mother_genotype.genotyped and father_genotype.genotyped:
                variant['Inheritance_model']['AD'] = False
            
    elif model == 'X_recessive':
        #If the individual is a male we only need to check if the mother carry the variant:
        if sex == 1:
            if not mother_genotype.has_variant:
                variant['Inheritance_model']['XR_dn'] = True
                if mother_genotype.genotyped:
                    variant['Inheritance_model']['XR'] = False
        #If female, both parents must have the variant otherwise denovo is true
        elif sex == 2:
            if not (mother_genotype.has_variant and father_genotype.has_variant):
                variant['Inheritance_model']['XR_dn'] = True
        #If both parents are genotyped but they both are not carriers XR is not true
                if mother_genotype.genotyped and father_genotype.genotyped:
                    variant['Inheritance_model']['XR'] = False
    
    elif model == 'X_dominant':
        #If the individual is a male we only need to look at the mother:
        if sex == 1:
            if not mother_genotype.has_variant:
                variant['Inheritance_model']['XD_dn'] = True
                if mother_genotype.genotyped:
                    variant['Inheritance_model']['XD'] = False
        #If female, one of the parents must have the variant otherwise denovo is true
        elif sex == 2:
            if not (mother_genotype.has_variant or father_genotype.has_variant):
                variant['Inheritance_model']['XD_dn'] = True
                if mother_genotype.genotyped and father_genotype.genotyped:
                    variant['Inheritance_model']['XD'] = False
    
    elif model == 'compound':
        individual_genotype_2 = variant_2['Genotypes'].get(individual, genotype.Genotype())
        mother_genotype_2 = variant_2['Genotypes'].get(mother_id, genotype.Genotype())
        father_genotype_2 = variant_2['Genotypes'].get(father_id, genotype.Genotype())
        # One of the variants must come from father and one from mother
        if not (mother_genotype.has_variant or mother_genotype_2.has_variant):
            variant['Inheritance_model']['AR_comp_dn'] = True
            variant_2['Inheritance_model']['AR_comp_dn'] = True
            if mother_genotype.genotyped and mother_genotype_2.genotyped:
                variant['Inheritance_model']['AR_comp'] = False
                variant_2['Inheritance_model']['AR_comp'] = False
        if not (father_genotype.has_variant or father_genotype_2.has_variant):
            variant['Inheritance_model']['AR_comp_dn'] = True
            variant_2['Inheritance_model']['AR_comp_dn'] = True
            if father_genotype.genotyped and father_genotype_2.genotyped:
                variant['Inheritance_model']['AR_comp'] = False
                variant_2['Inheritance_model']['AR_comp'] = False
    return
            
        

def main():
    from ped_parser import family, individual
    from interval_tree import interval_tree
    
    duo_family = family.Family(family_id = '1')
    sick_son = individual.Individual(ind='1', family='1',mother='3', father='2', sex=1, phenotype=2)
    healthy_father = individual.Individual(ind='2', family='1',mother='0', father='0', sex=1, phenotype=1)
    healthy_mother = individual.Individual(ind='3', family='1',mother='0', father='0', sex=2, phenotype=1)
    duo_family.add_individual(sick_son)
    duo_family.add_individual(healthy_mother)
    duo_family.add_individual(healthy_father)
    
    pp(duo_family.individuals)
    intervals = {ind_id:interval_tree.IntervalTree([[1,100, '1']], 1, 100) for ind_id in duo_family.individuals}

    pp(intervals)
    
    #Setup two variants with autosomal recessive compound pattern
    recessive_comp_simple_1 = {'CHROM':'1', 'POS':'5', 'ALT':'A', 'REF':'C', 'ID':'rs2230749',
                                             '1':'0|1', '2':'1|0', '3':'0|0'}
    
    recessive_comp_simple_2 = {'CHROM':'1', 'POS':'10', 'ALT':'C', 'REF':'T', 'ID':'.', 
                                            '1':'1|0', '2':'1|0', '3':'0|0'}
    
    
    batch = {'ABC':{'1_5_A_C':recessive_comp_simple_1, '1_10_C_T':recessive_comp_simple_2}}
    batch['haploblocks'] = intervals
    
    check_genetic_models(batch, duo_family, phased=True)
    for gene in batch:
        pp(batch[gene])
    

if __name__ == '__main__':
    main()

