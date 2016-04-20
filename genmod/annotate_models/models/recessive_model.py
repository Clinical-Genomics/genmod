#!/usr/bin/env python
# encoding: utf-8
"""
recessive_model.py

Checks is the Autosomal Homozygote recessive model is followed.

Created by Måns Magnusson on 2013-02-12.
Copyright (c) 2013 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function

import os
import sys

def check_recessive(variant, family, strict):
    """
    Check if the variant follows the autosomal recessive homozygote (AR_hom) 
    pattern in this family.
    
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
        individual_genotype = variant['genotypes'][individual]
        if strict:
            if not individual_genotype.genotyped:
                return False
        # The case where the individual is healthy:
        if family.individuals[individual].healthy:
        # If the individual is healthy and homozygote alt the model is broken.
            if individual_genotype.genotyped:
                if individual_genotype.homo_alt:
                    return False
                
        # The case when the individual is sick:
        elif family.individuals[individual].affected:
        # In the case of a sick individual it must be homozygote alternative 
        # for Autosomal recessive to be true.
        # Also, we can not exclude the model if no call.
            if individual_genotype.genotyped:
                if not individual_genotype.homo_alt:
                    return False
    
    return True


