#!/usr/bin/env python
# encoding: utf-8
"""
x_models.py

Checks is the X linked models are followed.

Created by Måns Magnusson on 2013-02-12.
Copyright (c) 2013 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function

import os
import sys

def check_X_recessive(variant, family, strict):
    """
    Check if the variant follows the x linked heterozygous (XR) pattern of 
    inheritance in this family.
    
    A variant is following the XR pattern if:
    
    Healthy:
        - Can not be homozygote alternative
        - If no call we can not exclude XR
        - Males can not have variant at all. This is added since sometimes males
            get called as heterozygotes but this should not be possible since 
            they only have one copy of the X chromosome.
        if strict:
            - Have to be homozygote reference(if male) or heterozygote(if female).
            - No call will return False
        
    Affected:
        - Have to be homozygote alternative(or heterozygote if male).
        - If no call we can not exclude AR
        if strict:
            - Have to be homozygote alternative(or heterozygote if male)
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
        individual_genotype = variant['genotypes'][individual]
        
        if strict:
            if not individual_genotype.genotyped:
                return False
        # The case where the individual is healthy
        if family.individuals[individual].healthy:
            # If individual is healthy and homozygote alternative 
            # the variant can not be deleterious:
            if individual_genotype.genotyped:
                if individual_genotype.homo_alt:
                    return False
                # If individual is male it can not have the variant at all
                if family.individuals[individual].sex == 1:
                    if individual_genotype.has_variant:
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
    return True

def check_X_dominant(variant, family, strict):
    """
    Check if the variant follows the x linked dominant (XD) pattern of 
    inheritance in this family.
    A variant is following the XD pattern if:
    
    Healthy:
        - Can not be homozygote alternative
        - Healthy females can be heterozygotes. This is possible since there 
            are several documented diseases where only one allele at a time is
            expressed during development.
        - If no call we can not exclude XR
        if strict:
            - Have to be homozygote reference (or heterozygote womens).
            - No call will return False
    
    Affected:
        - Have to be heterozygote.
        - If no call we can not exclude AR
        if strict:
            - Have to be heterozygote or homozygote(for males)
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
        individual_genotype = variant['genotypes'][individual]
        
        if strict:
            if not individual_genotype.genotyped:
                return False
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
        
        # The case when the individual is sick
        elif family.individuals[individual].affected:
        # If the individual is sick and homozygote ref it 
        # can not be x-linked-dominant
            if individual_genotype.genotyped:
                if individual_genotype.homo_ref:
                    return False
    return True

def main():
    pass

if __name__ == '__main__':
    main()

