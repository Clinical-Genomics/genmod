#!/usr/bin/env python
# encoding: utf-8
"""
dominant_model.py

Checks is the Autosomal Dominant model is followed.


Created by Måns Magnusson on 2013-02-12.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

from __future__ import print_function

import logging

logger = logging.getLogger(__name__)

def check_dominant(variant, family, strict=False):
    """
    Check if the variant follows the autosomal dominant (AD) pattern in 
    this family.
    
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
    
    SPECIAL CASE:
        If the variants is annotated with incomplete penetrance we allow healthy
        individuals to be carriers (i.e. healthy individuals can be het.)
    
    Args:
        variant: variant dictionary.
        family: A family object with the individuals
        strict: A boolean that tells if strict analyzis should be performed.
    
    Return:
        bool: depending on if the model is followed in these indivduals
    
    """
    
    for individual in family.individuals: 
        # Check in all individuals what genotypes that are in the trio based 
        # of the individual picked.
        logger.debug("Checking autosomal dominant pattern for variant {0},"\
        " individual: {1}".format(
            variant.get('variant_id', None),
            individual)
        )
        individual_genotype = variant['genotypes'][individual]
        if strict:
            if not individual_genotype.genotyped:
                return False
        # The case where the individual is healthy
        if family.individuals[individual].healthy:
            logger.debug("Individual {0} is healthy".format(individual))
            if individual_genotype.has_variant:
                if variant.get('reduced_penetrance', False):
                    if individual_genotype.homo_alt:
                        return False
                else:
                    return False
        
        elif family.individuals[individual].affected:
            logger.debug("Individual {0} is affected".format(individual))
            # The case when the individual is sick
            if individual_genotype.genotyped:
                if not individual_genotype.heterozygote:
                    return False
    
    return True

