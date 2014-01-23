#!/usr/bin/env python
# encoding: utf-8
"""
genetic_models.py

Genetic models take a family object with individuals and variants and annotates for each variant which models they follow in this family.

The following models are checked:

- Autosomal Dominant(AD)
- Autosomal Dominant De Novo(AD_DN)
- Autosomal Recessive(AR)
- Autosomal Recessive De Novo(AR_DN)
- Autosomal Recesive Compound.

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
from collections import Counter
from datetime import datetime
from genmod.variants import genotype
from genmod.utils import pair_generator
from time import sleep

def check_genetic_models(variant_batch, family, verbose = False):
    #A variant batch is a dictionary on the form {gene_id: {variant_id:variant_obj}}
    for gene, variants in variant_batch.items():
        compound_candidates = []
        compound_pairs = []
        if gene != '-':
            # First remove all variants that can't be compounds to reduce the number of lookup's:
            compound_candidates = [variant for variant_id, variant in check_compound_candidates(variants, family).items()]
            if len(compound_candidates) > 1:
            # Now check the compound candidates:
                compound_pairs = check_compound(compound_candidates, family)
            else:
                compound_pairs = []
        
        for variant_id, variant in variants.items():
                                
            # Only check X-linked for the variants in the X-chromosome:
            # For X-linked we do not need to check the other models
            if variant.chr == 'X':
                check_x_linked(variant, family)
            else:
            # Check the dominant model:
                check_dominant(variant, family)
            # Check the recessive model:
                check_recessive(variant, family)
            for pair in compound_pairs:
                variant_1 = variants[pair[0]]
                variant_2 = variants[pair[1]]
                # Add the compound pair id to each variant
                variant_1.ar_comp_variants[variant_2.variant_id] = 0
                variant_2.ar_comp_variants[variant_1.variant_id] = 0
                variant_1.ar_comp = True
                variant_2.ar_comp = True
            variant.check_models()
    return

def check_compound_candidates(variants, family):
    """Sort out the compound candidates, this function is used to reduce the number of potential candidates."""
    
    comp_candidates = dict(variants)#Make a copy of the dictionary to not change the original one.
    for individual in family.individuals:
        individual_variants = {}
        for variant_id in variants:
            variant = variants[variant_id]
            genotype = variant.get_genotype(individual.individual_id)
            if individual.affected():
                if not genotype.heterozygote:
                    if variant_id in comp_candidates:
                        del comp_candidates[variant_id]
                else:
                    individual_variants[variant_id] = variant
            else:#If individual is healthy or not known
                if genotype.homo_alt:
                    if variant_id in comp_candidates:
                        del comp_candidates[variant_id]
                elif genotype.heterozygote:
                    individual_variants[variant_id] = variant
        if individual.affected():
            #If the individual is sick then all potential compound candidates of a gene must exist in that individual.
            if len(individual_variants) > 1:
                for variant_id in comp_candidates:
                    if variant_id not in individual_variants:
                        del comp_candidates[variant_id]
            else:
                # If a sick individual dont have any compounds pairs there are no compound candidates.
                comp_candidates = {}
        else:
            #If an individual is healthy and have compound pairs they can not be deleterious:
            if len(individual_variants) > 1:
                for variant_id in individual_variants:
                    if variant_id in comp_candidates:
                        del comp_candidates[variant_id]
    return comp_candidates

def check_compound(list_of_variants, family):
    """Check which variants in the list that follow the compound heterozygous model. 
    We need to go through all variants and sort them into their corresponding genes 
    to see which that are candidates for compound heterozygotes first. 
    The cheapest way to store them are in a hash table. After this we need to go
     through all pairs, if both variants of a pair is found in a healthy individual
      the pair is not a deleterious compound heterozygote."""
    true_variants = []
    false_variants = []
            
    def add_variant_pair(variant_pair, variant_list):
        """Add the pairs that where found to be not true."""
        in_list = False
        for pair in variant_list:
            if Counter(variant_pair) == Counter(pair):
                in_list = True
        if not in_list:
            variant_list.append(variant_pair)
        
    # Returns a generator with all possible (unordered) pairs for this individual:
    my_pairs = pair_generator.Pair_Generator(list_of_variants)
    for pair in my_pairs.generate_pairs():
        variant_1 = pair[0]
        variant_2 = pair[1]
    # Check in all individuals what genotypes that are in the trio based of the individual picked.
        for individual in family.individuals:
            genotype_1 = variant_1.get_genotype(individual.individual_id)
            genotype_2 = variant_2.get_genotype(individual.individual_id)
            if individual.phenotype == 2:# The case where the individual is affected
                mother_id = individual.mother
                mother_genotype_1 = variant_1.get_genotype(mother_id)
                mother_genotype_2 = variant_2.get_genotype(mother_id)
                mother_phenotype = family.get_phenotype(mother_id)
                
                father_id = individual.father
                father_genotype_1 = variant_1.get_genotype(father_id)
                father_genotype_2 = variant_2.get_genotype(father_id)
                father_phenotype = family.get_phenotype(father_id)
                # If a parent has both variants and is unaffected it can not be a compound.
                # This will change when we get the phasing information.
                if ((mother_genotype_1.heterozygote and mother_genotype_2.heterozygote and mother_phenotype == 1) or (father_genotype_1.heterozygote and father_genotype_2.heterozygote and father_phenotype == 1)):
                    add_variant_pair([variant_1.variant_id, variant_2.variant_id], false_variants)
                else:
                    add_variant_pair([variant_1.variant_id, variant_2.variant_id], true_variants)
    compound_pairs = []
    pair_to_add = True
    for variant_pair in true_variants:
        for false_pair in false_variants:
            # If the variant pair is among the false ones do not add it
            if Counter(variant_pair) == Counter(false_pair):
                pair_to_add = False
        if pair_to_add:
            compound_pairs.append(variant_pair)
    return compound_pairs
                 
    

def check_x_linked(variant, family):
    """Check if the variant follows the x linked patter of inheritance in this family."""
    for individual in family.individuals:
        # Get the genotype for this variant for this individual
        genotype = variant.get_genotype(individual.individual_id)
    
        # The case where the individual is healthy
        if individual.phenotype == 1:
        #The case where the individual is a male
            if individual.sex == 1:
                if genotype.has_variant:
        # If the individual is healthy, male and have a variation it can not be x-linked.
                    variant.x_linked = False
                    variant.x_linked_dn = False
        
            #The case where the individual is a female
            elif individual.sex == 2:
                # If the individual is HEALTHY, female and is homozygote alternative it can not be x - linked.
                if genotype.homo_alt:
                    variant.x_linked = False
                    variant.x_linked_dn = False
    
        # The case when the individual is sick
        elif individual.phenotype == 2:
        #If the individual is sick and homozygote ref it can not be x-linked
            if genotype.homo_ref:
                variant.x_linked = False
                variant.x_linked_dn = False
            elif genotype.has_variant:
                check_parents('x_linked', individual, variant, family)
        # Else if phenotype is unknown we can not say anything about the model

def check_dominant(variant, family):
    """Check if the variant follows the dominant pattern in this family."""
    for individual in family.individuals: 
        # Check in all individuals what genotypes that are in the trio based of the individual picked.
    
        genotype = variant.get_genotype(individual.individual_id)
        if individual.phenotype == 1:# The case where the individual is healthy
            if genotype.has_variant:
                # If the individual is healthy and have a variation on one or both alleles it can not be dominant.
                variant.ad = False
                variant.ad_dn = False
        elif individual.phenotype == 2:
            # The case when the individual is sick
            if genotype.homo_ref:
                variant.ad = False
                variant.ad_dn = False
            else: 
            # Now the ind is sick and have a variant ≠ ref, check parents for de novo
                check_parents('dominant', individual, variant, family)
            # Else if phenotype is unknown we can not say anything about the model

def check_recessive(variant, family):
    """Check if the variant follows the autosomal recessive pattern in this family."""
    for individual in family.individuals:
        genotype = variant.get_genotype(individual.individual_id)
        # The case where the individual is healthy:
        if individual.phenotype == 1:
        # If the individual is healthy and homozygote alt the model is broken.
            if genotype.homo_alt:
                variant.ar = False
                variant.ar_dn = False
        # The case when the individual is sick:
        elif individual.phenotype == 2:
        # In the case of a sick individual it must be homozygote alternative for compound heterozygote to be true.
        # Also, we can not exclude the model if no call.
            if genotype.homo_ref or genotype.heterozygote:
                variant.ar = False
                variant.ar_dn = False
            else:
            #Models are followed but we need to check the parents to see if de novo is followed or not.
                check_parents('recessive', individual, variant, family)
                
def check_parents(model, individual, variant, family):
    """Check if information in the parents can tell us if model is de novo or not. Model in ['recessive', 'compound', 'dominant']."""
    sex = individual.sex
    individual_genotype = variant.get_genotype(individual.individual_id)

    mother_id = individual.mother
    mother_genotype = variant.get_genotype(mother_id)
    mother_phenotype = family.get_phenotype(mother_id)

    father_id = individual.father
    father_genotype = variant.get_genotype(father_id)
    father_phenotype = family.get_phenotype(father_id)


    if model == 'recessive':
        # If any of the parents doesent exist de novo will be true as the model is specified
        if mother_id != '0' and father_id != '0':
        # If both parents have the variant or if one of the parents are homozygote alternative, the de novo model is NOT followed, otherwise de novo is true.
            if ((mother_genotype.homo_alt or father_genotype.homo_alt) or 
                (mother_genotype.has_variant and father_genotype.has_variant)):
                variant.ar_dn = False
        if variant.ar_dn:# If de novo is true then the it is only de novo
            variant.ar = False
                
    elif model == 'dominant':
    # If one of the parents have the variant on any form the de novo model is NOT followed.
        if mother_genotype.has_variant or father_genotype.has_variant:
            variant.ad_dn = False
        if variant.ad_dn:# If variant is ad de novo then it is not ad
            variant.ad = False
            
    elif model == 'x_linked':
        #If the individual is a male:
        if sex == 1:
        #If any of the parents have the variant it is not dn
            if mother_genotype.has_variant or father_genotype.has_variant:
                variant.x_linked_dn = False
        #If female 
        elif sex == 2:
            #If 
            if individual_genotype.homo_alt:
                if ((mother_genotype.homo_alt or father_genotype.homo_alt) 
                    or (mother_genotype.has_variant and father_genotype.has_variant)):
                    variant.x_linked_dn = False
            elif individual_genotype.heterozygote:
                if mother_genotype.has_variant or father_genotype.has_variant:
                    variant.x_linked_dn = False
        if variant.x_linked_dn:
            variant.x_linked = False          
    

def main():
    pass


if __name__ == '__main__':
    main()

