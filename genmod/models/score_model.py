#!/usr/bin/env python
# encoding: utf-8
"""
score_model.py

Script that takes a batch of variants as input and generates a
rank score calculated using a weigthed sum matrix.

Created by Henrik Stranneheim on 2015-01-12.
Copyright (c) 2013 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function, unicode_literals

import sys
import os
import numbers
import re
from collections import defaultdict
from pprint import pprint as pp

# Import third party library
# https://github.com/mitsuhiko/logbook
from logbook import Logger, StderrHandler
log = Logger('Logbook')
log_handler = StderrHandler()


def convert_to_number(record_element):
    """
    Check if record_element is of type unicode and tries to make record_element
    into float.
    
    Args:
        record_element  (element) : Variant record
    
    Return:
        float:  If able
        None:   If unable to turn into float
    """
    try:
        float(record_element)
        return float(record_element)
    except ValueError:  # Not a float
        return None

def number_to_list(record_list, separator=','):
    """
    Checks if record list contains numbers and
    appends these to list and returns it.
    
    Args:
        record_list  (list) : Variant record
        separator    (string, optional): The separator to use in split
    
    Return:
        list:  Empty or with elements of type floats
    """
    score_list = []  # Create score list
    
    for score in record_list.split(separator):  # Split record list
        if isinstance(score, numbers.Number):  # Score is a number
            score_list.append(score)
        else:  # Try to convert into number
            score_element = convert_to_number(score)
            if score_element is not None:  # Score is a number
                score_list.append(score_element)
    
    return score_list

def string_to_dict(term_list):
    """
    Adds term_list as strings to a dictionnary. Saves
    term string as lower case since config file keys are
    automatically converted to lower case.
    
    Args:
        term_list  (list) : Variant record terms list
    
    Return:
        dict:  Empty or keys as strings
    """
    term_dict = {}  # Create term dict
    
    for term in term_list:
        if term is not None:
            term_dict[str(term).lower()] = str(term.lower())
    
    return term_dict


def split_record(record_list, separator_list):
    """
    Splits a list on multiple separators
    
    Args:
        record_list     (list) : Variant record
        separator_list  (lsit) : Separators to use in split
    
    Return:
        dict:  Empty or keys as strings
    """
    
    term = "["  # Create string for resplit
    ## Collect separator(s)
    for separator in separator_list:
        term += separator
    term += "]"  # Closing the string
    term_list = filter(None, re.split(term, record_list))
    term_dict = string_to_dict(term_list)
    
    return term_dict


def number_equal(element, value, score):
    """
    Check if element equals config value
    
    Args:
        element  (float) : Usually vcf record
        value    (float) : Config value
        score    (integer) : config score
    
    Return:
        Float:  Score
    """
    if element == value:
        return score
    return None


def number_less_than(element, value, score):
    """
    Check if element is lower than config value.
    
    Args:
        element  (float) : Usually vcf record
        value    (float) : Config value
        score    (integer) : config score
    
    Return:
        Float:  Score
    """
    if element < value:
        return score
    return None

def number_less_equal(element, value, score):
    """
    Check if element is lower than or equal config value.
    
    Args:
        element  (float) : Usually vcf record
        value    (float) : Config value
        score    (integer) : config score
    
    Return:
        Float:  Score
    """
    if element <= value:
        return score
    return None


def number_greater_than(element, value, score):
    """
    Check if element is greater than config value.
    
    Args:
        element  (float) : Usually vcf record
        value    (float) : Config value
        score    (integer) : config score

    Return:
        Float:  Score
    """
    if element > value:
        return score
    return None

def number_greater_equal(element, value, score):
    """
    Check if element is greater than or equals config value.
    
    Args:
        element  (float) : Usually vcf record
        value    (float) : Config value
        score    (integer) : config score
    
    Return:
        Float:  Score
    """
    if element >= value:
        return score
    return None

def score_float(alt, record_list, value_dict, operation_dict,
                score_dict, record_aggregate=max):
    """
    Calculates the performance score for each element.
    
    Args:
        alt               (string) : Usually vcf key
        record_list       (list) : list of record elements
        value_dict        (dict) : Dictionnary of alt config value
        operation_dict    (dict) : Dictionnary of alt comparisons
        score_dict        (dict) : Dictionnary of alt config score
        record_aggregate  (string, optional) : Method of record aggregation
    
    Return:
        Float:  Alternative final score for variant
    """
    score_list = []  # Collects all scores for vcf record
    final_score = 0
    
    ## Can be multiple fields within vcf record
    for element in record_list:
    
        for key in value_dict[alt]:
            score = None
            if operation_dict[alt][key] == "lt":
                score = number_less_than(element,
                                  float(value_dict[alt][key]),
                                  float(score_dict[alt][key])
                                  )
            elif operation_dict[alt][key] == "le":
                score = number_less_equal(element,
                                  float(value_dict[alt][key]),
                                  float(score_dict[alt][key])
                                  )
            elif operation_dict[alt][key] == "gt":
                score = number_greater_than(element,
                                  float(value_dict[alt][key]),
                                  float(score_dict[alt][key])
                                  )
            elif operation_dict[alt][key] == "ge":
                score = number_greater_equal(element,
                                  float(value_dict[alt][key]),
                                  float(score_dict[alt][key])
                                  )
            elif operation_dict[alt][key] == "e":
                score = number_equal(element,
                                  float(value_dict[alt][key]),
                                  float(score_dict[alt][key])
                                  )
            if score:
                score_list.append(score)
    
    if score_list:
        
        if record_aggregate == "max":  # Keep highest score
            final_score = max(score_list)
        
        if record_aggregate == "min":  # Keep lowest score
            final_score = min(score_list)
    
    return final_score

def score_string(alt, term_dict, value_dict, operation_dict,
                 score_dict, record_aggregate=max):
    """
    Calculates the performance score for each element.
    
    Args:
        alt               (string) : Usually vcf key
        term_dict         (dict) : Dictionnary of alt config term
        value_dict        (dict) : Dictionnary of alt config value
        operation_dict    (dict) : Dictionnary of alt comparisons
        score_dict        (dict) : Dictionnary of alt config score
        record_aggregate  (string, optional) : Method of record aggregation

    Return:
        Float:  Alternative final score for variant
    """
    score_list = []  # Collects all scores for vcf record
    final_score = 0  # Final score to return
    
    ## Can be multiple fields within vcf record
    for term in term_dict:
        
        if term in value_dict[alt]:  # Term present in record
            score_list.append(float(score_dict[alt][term]))
        
        elif term == ".":  # No information on element
            score_list.append(float(score_dict[alt]['notreported']))
        
        if score_list:
            
            ## Remove None type
            score_list = [i for i in score_list if i is not None]
            if record_aggregate == "max":  # Keep highest score
                final_score = max(score_list)
            
            if record_aggregate == "min":  # Keep lowest score
                final_score = min(score_list)
    
    return final_score


def evaluate_float(alt, category, value_dict, operation_dict,
                   score_dict, perf_score_dict, separator_list=[],
                   record=None, record_aggregate=max, verbose=False):
    """
    Evalutes vcf data type floats. Fileds can be separated on supplied
    delimiter.
    
    Args:
        alt               (string) : Usually vcf key
        category          (string) : alt config category
        value_dict        (dict) : Dictionnary of alt config value
        operation_dict    (dict) : Dictionnary of alt comparisons
        score_dict        (dict) : Dictionnary of alt config score
        perf_score_dict   (dict) : Dictionnary of category scores
        separator_list    (list) : List of vcf record field delimiters
        record            (unicode, optional) : vcf record
        record_aggregate  (string, optional) : Method of record aggregation
        verbose           (boolean, optional) : Enable informative print
    
    Return:
        Float:  Alternative final score for variant
    """
    if len(separator_list) > 0:  # Multiple elements
        record_list = number_to_list(record, separator_list[0])
    
    else:
        record_list = number_to_list(record)
    
    if len(record_list) > 0:
        if verbose == 2:
            log.info("Record elements:" + str(record_list))
        perf_score_dict[category].append(score_float(alt,
                                         record_list,
                                         value_dict,
                                         operation_dict,
                                         score_dict,
                                         record_aggregate))
    
    else:  # For instance "."
        if (score_dict[alt]['notreported']):
            perf_score_dict[category].append(float(
                                             score_dict[alt]
                                             ['notreported']))


def score_variants(batch, predicted_models=[], alt_dict=None, score_dict=None,
                   value_dict=None, operation_dict=None, verbose=False):
    """
    Score a variant object according to a weigthed sum model and inserts
    score in variant object info field.
    
    Args:
        batch             (object) : Variant object
        predicted_models  (list, optional): List of predicted genetic models
        alt_dict          (dict, optional) : Dictionnary of config alternatives
        score_dict        (dict, optional) : Dictionnary of alt config score
        value_dict        (dict, optional) : Dictionnary of alt config value
        operation_dict    (dict, optional) : Dictionnary of alt comparisons
        verbose           (boolean, optional) : Enable informative print

    Return:
        None:
    """
    for variant_id in batch:
        
        variant = batch[variant_id]  # Get variant
        if verbose == 2:
            log.info("Variant Line: " + str(variant))
        variant_score = 0  # Variant score
        perf_score_dict = defaultdict(list)  # Performance score dict
        info_dict = variant.get('info_dict', {})  # Create dict of info field
        category_dict = {}  # Collects all alternative categories
        
        for alt in alt_dict:  # Config alternatives
            
            category = ""
            record_aggregate = "max"  # Set default
            separator_list = []
            
            if verbose == 2:
                log.info("Alternative: " + alt)
            
            category = alt_dict[alt]['category']  # Alias
            
            if 'category_aggregate' in alt_dict[alt]:
                category_dict[category] = alt_dict[alt]['category_aggregate']
            
            if 'field_separators' in alt_dict[alt]:
                separator_list = alt_dict[alt]['field_separators'].split('_')
            
            if 'record_aggregate' in alt_dict[alt]:
                record_aggregate = alt_dict[alt]['record_aggregate']
            
            if alt in variant:  # Not a vcf INFO key
                record = variant.get(alt, None)
            
            else:  # Vcf INFO key
                record_list = info_dict.get(alt, None)
                record = None
                
                ## Create string from list
                if record_list:
                    record = ','.join(record_list)
                
            if record is not None:  # Record exists in vcf
                if verbose == 2:
                    log.info("Record:" + record)
                
                ## Number comparisons
                if alt_dict[alt]['data_type'] == 'float':
                    evaluate_float(alt, category,
                                   value_dict, operation_dict,
                                   score_dict, perf_score_dict,
                                   separator_list, record,
                                   record_aggregate, verbose)
                ## String comparisons
                if alt_dict[alt]['data_type'] == 'string':
                    ## Split vcf info field
                    if (separator_list) and (len(separator_list) > 0):
                        term_dict = split_record(record, separator_list)
                    
                    else:
                        term_dict = {}  # Create term dict
                        term_dict[str(record).lower()] = str(record.lower())
                    
                    perf_score_dict[category].append(score_string(alt,
                                                     term_dict,
                                                     value_dict,
                                                     operation_dict,
                                                     score_dict,
                                                     record_aggregate))
            else:  # No record found in vcf for this alternative
                if (score_dict[alt]['notreported']):
                    perf_score_dict[category].append(float(score_dict[alt]
                                                           ['notreported']))
            
            if verbose == 2:
                log.info(category + ": " + str(perf_score_dict[category]))
            
        for category in category_dict:  # All should be numbers from here on
            category_score = 0
            
            if category_dict[category] == 'min':
                category_score = min(perf_score_dict[category])
            
            elif category_dict[category] == 'max':
                category_score = max(perf_score_dict[category])
            
            elif category_dict[category] == 'sum':
                category_score = sum(perf_score_dict[category])
            
            variant_score += category_score
            
        if verbose == 2:
            log.info("Variant score " + str(variant_score))
        
        # Added the individual rank score to the variant.
        variant['Individual_rank_score'] = int(variant_score)
        
    return


def main():
    pass


if __name__ == '__main__':
    main()
