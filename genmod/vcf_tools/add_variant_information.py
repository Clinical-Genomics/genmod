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

def replace_vcf_info(keyword, annotation, variant_line=None, variant_dict=None):
    """Replace the information of a info field of a vcf variant line.
    
    
    Arguments:
        variant_line (str): A vcf formatted variant line
        variant_dict (dict): A variant dictionary
        keyword (str): The info field key
        annotation (str): If the annotation is a key, value pair
                          this is the string that represents the value
    
    Returns:
        variant_line (str): A annotated variant line
    """
    logger = logging.getLogger(__name__)
    
    new_info = '{0}={1}'.format(keyword, annotation)
    
    logger.debug("Replacing the variant information {0}".format(new_info))
    
    fixed_variant = None
    new_info_list = []
    
    if variant_line:
        logger.debug("Adding information to a variant line")
        splitted_variant = variant_line.rstrip('\n').split('\t')
        logger.debug("Adding information to splitted variant line")
        old_info = splitted_variant[7]
        if old_info == '.':
            new_info_string = new_info
        else:
            splitted_info_string = old_info.split(';')
            for info in splitted_info_string:
                splitted_info_entry = info.split('=')
                if splitted_info_entry[0] == keyword:
                    new_info_list.append(new_info)
                else:
                    new_info_list.append(info)
            new_info_string = ';'.join(new_info_list)
        
        splitted_variant[7] = new_info_string
        
        fixed_variant = '\t'.join(splitted_variant)
    
    elif variant_dict:
        logger.debug("Adding information to a variant dict")
        old_info = variant_dict['INFO']
        
        if old_info == '.':
            variant_dict['INFO'] = new_info
        else:
            for info in old_info.split(';'):
                splitted_info_entry = info.split('=')
                if splitted_info_entry[0] == keyword:
                    new_info_list.append(new_info)
                else:
                    new_info_list.append(info)
            new_info_string = ';'.join(new_info_list)
        
        variant_dict['INFO'] = new_info_string
        fixed_variant = variant_dict
    
    return fixed_variant

def add_vcf_info(keyword, variant_line=None, variant_dict=None, annotation=None):
    """
    Add information to the info field of a vcf variant line.
    
    Arguments:
        variant_line (str): A vcf formatted variant line
        keyword (str): The info field key
        annotation (str): If the annotation is a key, value pair
                          this is the string that represents the value
    
    Returns:
        fixed_variant : str if variant line, or dict if variant_dict
    """
    logger = logging.getLogger(__name__)
    
    if annotation:
        annotation = str(annotation)
        new_info = '{0}={1}'.format(keyword, annotation)
    else:
        new_info = keyword
    
    logger.debug("Adding new variant information {0}".format(new_info))
    
    fixed_variant = None
    
    if variant_line:
        logger.debug("Adding information to a variant line")
        splitted_variant = variant_line.rstrip('\n').split('\t')
        logger.debug("Adding information to splitted variant line")
        old_info = splitted_variant[7]
        if old_info == '.':
            splitted_variant[7] = new_info
        else:
            splitted_variant[7] = "{0};{1}".format(splitted_variant[7], new_info)
        
        fixed_variant = '\t'.join(splitted_variant)
    
    elif variant_dict:
        logger.debug("Adding information to a variant dict")
        old_info = variant_dict['INFO']
        if old_info == '.':
            variant_dict['INFO'] = new_info
        else:
            variant_dict['INFO'] = "{0};{1}".format(old_info, new_info)
        fixed_variant = variant_dict
    
    return fixed_variant
