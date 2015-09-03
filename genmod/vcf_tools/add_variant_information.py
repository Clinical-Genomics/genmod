#!/usr/bin/env python
# encoding: utf-8
"""
add_variant_information.py

Add information to a variant line in the vcf format

Created by Måns Magnusson on 2015-08-25.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function

import logging

def add_vcf_info(variant_line, keyword, annotation=None):
    """
    Add information to the info field of a vcf variant line.
    
    Arguments:
        variant_line (str): A vcf formatted variant line
        keyword (str): The info field key
        annotation (str): If the annotation is a key, value pair
                          this is the string that represents the value
    
    Returns:
        variant_line (str): A annotated variant line
    """
    logger = logging.getLogger(__name__)
    splitted_variant = variant_line.rstrip('\n').split('\t')
    
    if annotation:
        new_info = '='.join([keyword, annotation])
    else:
        new_info = keyword

    logger.debug("New info: {0}".format(new_info))
    logger.debug("Adding new info to variant line")
        
    splitted_variant[7] = "{0};{1}".format(splitted_variant[7], new_info)
    
    
    return '\t'.join(splitted_variant)
