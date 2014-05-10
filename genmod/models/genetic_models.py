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

def check_genetic_models(variant_batch, family, verbose = False, phased = False, proc_name = None):
    #A variant batch is a dictionary on the form {gene_id: {variant_id:variant_dict}}
    # Start by getting the genotypes for each variant:
    individuals = family.individuals.values()
    intervals = variant_batch.pop('haploblocks', {})
    for gene in variant_batch:
        for variant_id in variant_batch[gene]:
            genotypes = {}
            for individual in family.individuals:
                try:
                    # print(variant_batch[gene][variant_id])
                    gt_info = variant_batch[gene][variant_id][individual].split(':')[0]
                except KeyError:# If individual not in variant file
                    if verbose:
                        print('Warning! Individual %s is not in variant file!' % individual)
                    gt_info = './.'
                
                individual_genotype = genotype.Genotype(GT=gt_info)
                genotypes[individual] = individual_genotype
            variant_batch[gene][variant_id]['Genotypes'] = genotypes
            variant_batch[gene][variant_id]['Compounds'] = {}
            # Add information of models followed:
            variant_batch[gene][variant_id]['Inheritance_model'] = {'XR' : True, 'XR_dn' : True, 'XD' : True, 
            'XD_dn' : True, 'AD' : True, 'AD_dn' : True, 'AR_hom' : True, 'AR_hom_dn' : True, 'AR_comp' : False, 
            'AR_comp_dn' : False}
    # Now check the genetic models:
    for gene in variant_batch:
        compound_candidates = []
        compound_pairs = []
        # We look at compounds only when variants are in genes:
        if gene != '-':
            # First remove all variants that can't be compounds to reduce the number of lookup's:
            compound_candidates = check_compound_candidates(variant_batch[gene], family)
        
        for variant_id in variant_batch[gene]:
            
            # Only check X-linked for the variants in the X-chromosome:
            # For X-linked we do not need to check the other models
            if variant_batch[gene][variant_id]['CHROM'] == 'X':
                check_X_recessive(variant_batch[gene][variant_id], family)
                check_X_dominant(variant_batch[gene][variant_id], family)
                variant_batch[gene][variant_id]['Inheritance_model']['AD'] = False
                variant_batch[gene][variant_id]['Inheritance_model']['AD_dn'] = False
                variant_batch[gene][variant_id]['Inheritance_model']['AR_hom'] = False
                variant_batch[gene][variant_id]['Inheritance_model']['AR_hom_dn'] = False
            else:
                variant_batch[gene][variant_id]['Inheritance_model']['XR'] = False
                variant_batch[gene][variant_id]['Inheritance_model']['XR_dn'] = False
                variant_batch[gene][variant_id]['Inheritance_model']['XD'] = False
                variant_batch[gene][variant_id]['Inheritance_model']['XD_dn'] = False
            # Check the dominant model:
                check_dominant(variant_batch[gene][variant_id], family)
            # Check the recessive model:
                check_recessive(variant_batch[gene][variant_id], family)
            
        if len(compound_candidates) > 1:
            
            compound_pairs = pair_generator.Pair_Generator(compound_candidates)
            for pair in compound_pairs.generate_pairs():
                # Add the compound pair id to each variant
                if check_compounds(variant_batch[gene][pair[0]], variant_batch[gene][pair[1]], family, phased, intervals):
                    variant_batch[gene][pair[0]]['Compounds'][pair[1]] = 0
                    variant_batch[gene][pair[1]]['Compounds'][pair[0]] = 0
                    variant_batch[gene][pair[0]]['Inheritance_model']['AR_comp'] = True
                    variant_batch[gene][pair[1]]['Inheritance_model']['AR_comp'] = True
    return

def check_compound_candidates(variants, family):
    """Sort out the compound candidates, this function is used to reduce the number of potential candidates."""
    #Make a copy of the dictionary to not change the original one. {variant_id:variant_dict}
    comp_candidates = {variant_id: variants[variant_id] for variant_id in variants if 
                        variants[variant_id].get('comp_candidate',True)}
    for individual in family.individuals:
        individual_variants = {}
        for variant_id in dict((variant_id, comp_candidates[variant_id]) for variant_id in comp_candidates):
            individual_genotype = variants[variant_id]['Genotypes'].get(individual, genotype.Genotype())
            if individual_genotype.homo_alt:
                comp_candidates.pop(variant_id,0)
            # If an individual is affected:
            else:
                if family.individuals[individual].affected():
                    # It has to be heterozygote for the variant to be a candidate
                    if not individual_genotype.heterozygote:
                        comp_candidates.pop(variant_id, 0)
                    else: # Now we have a potential candidate:
                        individual_variants[variant_id] = ''
        #If the individual is sick then all potential compound candidates of a gene must exist in that individual.
        #So we remove all variants that the sick individual don't have
        if family.individuals[individual].affected():
            if len(individual_variants) > 1:
                for variant_id in comp_candidates:
                    if variant_id not in individual_variants:
                        comp_candidates.pop(variant_id,0)
            else:
                # If a sick individual dont have any compounds pairs there are no compound candidates.
                comp_candidates = {}
    return list(comp_candidates.keys())

def check_compounds(variant_1, variant_2, family, phased, intervals):
    """Check which variants in the list that follow the compound heterozygous model. At this stage we\
        know that none of the individuals are homozygote alternative for the variants."""
    
    # Check in all individuals what genotypes that are in the trio based of the individual picked.
    for individual in family.individuals:
        genotype_1 = variant_1['Genotypes'].get(individual, genotype.Genotype())
        genotype_2 = variant_2['Genotypes'].get(individual, genotype.Genotype())
        #check if variants are in the same phased interval:
        if phased:
            variant_1_interval = intervals[individual].find_range([int(variant_1['POS']),int(variant_1['POS'])])
            variant_2_interval = intervals[individual].find_range([int(variant_1['POS']),int(variant_1['POS'])])
        if family.individuals[individual].phenotype != 2:
        # If the individual is not sick and have both variants it can not be compound
            if genotype_1.has_variant and genotype_2.has_variant:
                if phased:
                # If the family is phased we need to check if a healthy individual have both variants on same allele
                    if variant_1_interval == variant_2_interval:
                        # If the variants are on different alleles it can not be a compound pair:
                        if genotype_1.allele_1 == '0' and genotype_2.allele_1 != '0':
                            return False        
                #In this case we can not tell if the variants are on the same haplotype so we assume that compound is not ok
                    else:
                        return False
                else:
                    return False
        # The case where the individual is affected
        else:
            if phased:
                #If the individual is sick and phased it has to have one variant on each allele
                if variant_1_interval == variant_2_interval:
                    if genotype_1.allele_1 == genotype_2.allele_1 or genotype_1.allele_2 == genotype_2.allele_2:
                        return False
                    
            elif family.individuals[individual].has_parents:
                mother_id = family.individuals[individual].mother
                mother_genotype_1 = variant_1['Genotypes'].get(mother_id, genotype.Genotype())
                mother_genotype_2 = variant_2['Genotypes'].get(mother_id, genotype.Genotype())
                mother_phenotype = family.get_phenotype(mother_id)
                
                father_id = family.individuals[individual].father
                father_genotype_1 = variant_1['Genotypes'].get(father_id, genotype.Genotype())
                father_genotype_2 = variant_2['Genotypes'].get(father_id, genotype.Genotype())
                father_phenotype = family.get_phenotype(father_id)
                # If a parent has both variants and is unaffected it can not be a compound.
                # This will change when we get the phasing information.
                if ((mother_genotype_1.heterozygote and mother_genotype_2.heterozygote and mother_phenotype == 1) or 
                        (father_genotype_1.heterozygote and father_genotype_2.heterozygote and father_phenotype == 1)):
                    return False
    return True


def check_dominant(variant, family):
    """Check if the variant follows the dominant pattern in this family."""
    for individual in family.individuals: 
        # Check in all individuals what genotypes that are in the trio based of the individual picked.
        individual_genotype = variant['Genotypes'].get(individual, genotype.Genotype())
        if family.individuals[individual].phenotype == 1:# The case where the individual is healthy
            if individual_genotype.has_variant:
                # If the individual is healthy and have a variation on one or both alleles it can not be dominant.
                variant['Inheritance_model']['AD'] = False
                variant['Inheritance_model']['AD_dn'] = False
                return
        elif family.individuals[individual].phenotype == 2:
            # The case when the individual is sick
            if individual_genotype.genotyped:
                if not individual_genotype.heterozygote:
                # Individual has to be heterozygote i AD can be true
                    variant['Inheritance_model']['AD'] = False
                    variant['Inheritance_model']['AD_dn'] = False
                    return
            # Now the ind is sick and have a variant ≠ ref, check parents for de novo
            if family.individuals[individual].has_parents:
                check_parents('dominant', individual, variant, family)
    return

def check_recessive(variant, family):
    """Check if the variant follows the autosomal recessive pattern in this family."""
    for individual in family.individuals:
        individual_genotype = variant['Genotypes'].get(individual, genotype.Genotype())
        # The case where the individual is healthy:
        if family.individuals[individual].phenotype == 1:
        # If the individual is healthy and homozygote alt the model is broken.
            if individual_genotype.homo_alt:
                variant['Inheritance_model']['AR_hom'] = False
                variant['Inheritance_model']['AR_hom_dn'] = False
                return
        # The case when the individual is sick:
        elif family.individuals[individual].phenotype == 2:
        # In the case of a sick individual it must be homozygote alternative for Autosomal recessive to be true.
        # Also, we can not exclude the model if no call.
            if not individual_genotype.homo_alt:
                variant['Inheritance_model']['AR_hom'] = False
                variant['Inheritance_model']['AR_hom_dn'] = False
                return
            #Models are followed but we need to check the parents to see if de novo is followed or not.
            elif family.individuals[individual].has_parents:
                check_parents('recessive', individual, variant, family)
    return

def check_X_recessive(variant, family):
    """Check if the variant follows the x linked heterozygous pattern of inheritance in this family."""
    for individual in family.individuals:
        # Get the genotype for this variant for this individual
        individual_genotype = variant['Genotypes'].get(individual, genotype.Genotype())
        
        # The case where the individual is healthy
        if not family.individuals[individual].affected():
            # If individual is healthy and homozygote alternative the variant can not be deleterious:
            if individual_genotype.homo_alt:
                variant['Inheritance_model']['XR'] = False
                variant['Inheritance_model']['XR_dn'] = False
                return
        #The case where the individual is a male
            if family.individuals[individual].sex == 1:
                if individual_genotype.has_variant:
        # If the individual is healthy, male and have a variation it can not be x-linked-recessive.
                    variant['Inheritance_model']['XR'] = False
                    variant['Inheritance_model']['XR_dn'] = False
                    return
        
        # The case when the individual is sick
        elif family.individuals[individual].affected():
        #If the individual is sick and homozygote ref it can not be x-recessive
            if individual_genotype.homo_ref:
                variant['Inheritance_model']['XR'] = False
                variant['Inheritance_model']['XR_dn'] = False
                return
        # Women have to be hom alt to be sick (almost allways carriers)
            elif family.individuals[individual].sex == 2:
                if individual_genotype.genotyped:
                    if not individual_genotype.homo_alt:
                        variant['Inheritance_model']['XR'] = False
                        variant['Inheritance_model']['XR_dn'] = False
                        return
            if family.individuals[individual].has_parents:
                check_parents('X_recessive', individual, variant, family)
    return

def check_X_dominant(variant, family):
    """Check if the variant follows the x linked dominant pattern of inheritance in this family."""
    for individual in family.individuals:
        # Get the genotype for this variant for this individual
        individual_genotype = variant['Genotypes'].get(individual, genotype.Genotype())
        # The case where the individual is healthy
        if not family.individuals[individual].affected():
        # Healthy womans can be carriers but not homozygote:
            if family.individuals[individual].sex == 2:
                if individual_genotype.homo_alt:
                    variant['Inheritance_model']['XD'] = False
                    variant['Inheritance_model']['XD_dn'] = False
                    return
        # Males can not carry the variant:
            elif family.individuals[individual].sex == 1:
                if individual_genotype.has_variant:
                    variant['Inheritance_model']['XD'] = False
                    variant['Inheritance_model']['XD_dn'] = False
                    return
        # The case when the individual is sick
        elif family.individuals[individual].affected():
        #If the individual is sick and homozygote ref it can not be x-linked-dominant
            if individual_genotype.homo_ref:
                variant['Inheritance_model']['XD'] = False
                variant['Inheritance_model']['XD_dn'] = False
                return
            elif individual_genotype.has_variant:
                if family.individuals[individual].has_parents:
                    check_parents('X_dominant', individual, variant, family)
    return

def check_parents(model, individual, variant, family):
    """Check if information in the parents can tell us if model is de novo or not. Model in ['recessive', 'compound', 'dominant']."""
    sex = family.individuals[individual].sex
    individual_genotype = variant['Genotypes'].get(individual, genotype.Genotype())

    mother_id = family.individuals[individual].mother
    mother_genotype = variant['Genotypes'].get(mother_id, genotype.Genotype())
    mother_phenotype = family.get_phenotype(mother_id)

    father_id = family.individuals[individual].father
    father_genotype = variant['Genotypes'].get(father_id, genotype.Genotype())
    father_phenotype = family.get_phenotype(father_id)


    if model == 'recessive':
        # If a parnent is homozygote or if both parents are heterozygote the variant is not denovo
        if ((mother_genotype.homo_alt or father_genotype.homo_alt) or
                (mother_genotype.has_variant and father_genotype.has_variant)):
            variant['Inheritance_model']['AR_hom_dn'] = False
        # If both parents are called but none of the above is fullfilled it is denovo
        elif mother_genotype.genotyped and father_genotype.genotyped:
                variant['Inheritance_model']['AR_hom'] = False
                    
    elif model == 'dominant':
        # If one or both parents are affected it is de novo if none of them have a variant
        if mother_genotype.has_variant or father_genotype.has_variant:
            variant['Inheritance_model']['AD_dn'] = False
        # If both parents are called but none of them carry the variant it is denovo
        elif mother_genotype.genotyped and father_genotype.genotyped:
            variant['Inheritance_model']['AD'] = False
            
    elif model == 'X_recessive':
        #If the individual is a male we only need if the mother carry the variant:
        if sex == 1:
            if mother_genotype.has_variant:
                variant['Inheritance_model']['XR_dn'] = False
            elif mother_genotype.genotyped:
                variant['Inheritance_model']['XR'] = False
        #If female, both parents must have the variant otherwise denovo is true
        elif sex == 2:
            if (mother_genotype.has_variant and father_genotype.has_variant):
                variant['Inheritance_model']['XR_dn'] = False
        #If both parents are genotyped but they both are not carriers XR is not true
            elif mother_genotype.genotyped and father_genotype.genotyped:
                variant['Inheritance_model']['XR'] = False
    
    elif model == 'X_dominant':
        #If the individual is a male we only need to look at the mother:
        if sex == 1:
            if mother_genotype.has_variant:
                variant['Inheritance_model']['XD_dn'] = False
            elif mother_genotype.genotyped:
                variant['Inheritance_model']['XD'] = False
        #If female, one of the parents must have the variant otherwise denovo is true
        elif sex == 2:
            if (mother_genotype.has_variant or father_genotype.has_variant):
                variant['Inheritance_model']['XD_dn'] = False
            elif mother_genotype.genotyped and father_genotype.genotyped:
                variant['Inheritance_model']['XD'] = False
        

def main():
    from ped_parser import family, individual
    from genmod.utils import interval_tree
    
    duo_family = family.Family(family_id = '1')
    sick_son = individual.Individual(ind='1', family='1',mother='3', father='0', sex=1, phenotype=2)
    healthy_mother = individual.Individual(ind='3', family='1',mother='0', father='0', sex=2, phenotype=1)
    duo_family.add_individual(sick_son)
    duo_family.add_individual(healthy_mother)
    
    pp(duo_family.individuals)
    intervals = {ind_id:interval_tree.IntervalTree([[1,100, '1']], 0, 1, 1, 100) for ind_id in duo_family.individuals}

    pp(intervals)
    
    #Setup two variants with autosomal recessive compound pattern
    recessive_comp_simple_1 = {'CHROM':'1', 'POS':'5', 'ALT':'A', 'REF':'C', 'ID':'rs2230749',
                                             '1':'0|1', '3':'0|0'}
    
    recessive_comp_simple_2 = {'CHROM':'1', 'POS':'10', 'ALT':'C', 'REF':'T', 'ID':'.', 
                                            '1':'1|0', '3':'0|1'}
    
    
    batch = {'ABC':{'1_5_A_C':recessive_comp_simple_1, '1_10_C_T':recessive_comp_simple_2}}
    batch['intervals'] = intervals
    
    check_genetic_models(batch, duo_family, phased=True)
    for gene in batch:
        pp(batch[gene])
    

if __name__ == '__main__':
    main()

