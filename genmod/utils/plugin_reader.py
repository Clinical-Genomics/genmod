#!/usr/bin/env python
# coding: utf-8
"""
plugin_reader.py

Read a config file with plugins that define a scoring model.


Created by Henrik Stranneheim on 2015-01-12.
Copyright (c) 2015 __MyCompanyName__. All rights reserved.
"""

from __future__ import print_function, unicode_literals

import sys
import configparser

from collections import defaultdict
from pprint import pprint as pp

# Import third party library
# https://github.com/mitsuhiko/logbook
from logbook import Logger, StderrHandler
log = Logger('Logbook')
log_handler = StderrHandler()


def read_config(config_file):
    """
    Read the config file.
    
    Args:
        config_file (file) : Plugin config file
    
    Returns:
        object: config object
    """
    config = configparser.ConfigParser()
    
    config.read(config_file)
    
    return config


def check_vcf_config(key, not_info_list, my_vcf_parser, config):
    """
    Check that the alternatives in config is present within vcf INFO keys.
    
    Args:
        key            (string) : Plugin config alternative key
        not_info_list  (list) : List of none vcf INFO keys
        my_vcf_parser  (object) : vcf_parser object
        config         (object) : Config object
    
    Returns:
        None:
    """
    if key in not_info_list:
        pass
    elif key not in my_vcf_parser.metadata.info_dict:
        
        log.notice('Plug-in.' +
                   config['Plug-in']['name'] +
                   ': Could not find key="' + key +
                   '" in INFO fields of variant file'
                   + '\n'
                   )
        log.notice('Aborting ranking' + '\n')
        sys.exit()
    return

def collectKeys(config_file, my_vcf_parser, verbose):
    """
    Collects database keys determined by plug-in and saves them in dict.
    
    Args:
        config_file (file) : plugin config file
        my_vcf_parser (object) : vcf_parser object
    
    Returns:
        dict: dict[collection][list_of_keys]
    """
    ## Create config object
    config = read_config(config_file)
    
    ## Create dictionnaries to store config info
    alt_dict = dict()
    score_dict = dict()
    value_dict = dict()
    operation_dict = dict()
    ## Alternatives not found in vcf INFO
    not_info_list = ["FILTER"]
    
    ## Collect keys for plug-in
    for section in config.items():
        
        section_name = section[0]  # Alias
        if verbose:
            log.info("Config alternative: " + str(config[section_name]))
            if 'version' in config[section_name]:
                log.info("Plugin version: " + str(config[section_name]['version']))
        
        ## Only vcf records
        if 'category' in config[section_name]:
            
            ## Test for presence of vcf INFO keys
            check_vcf_config(section_name, not_info_list, my_vcf_parser, config)
            
            alt_dict[section_name] = dict()
            score_dict[section_name] = dict()
            value_dict[section_name] = dict()
            operation_dict[section_name] = dict()
            
            ## Save vcf_info_key and general parameters
            alt_dict[section_name]['category'] = config[section_name]['category']
            alt_dict[section_name]['data_type'] = config[section_name]['data_type']
            
            if 'category_aggregate' in config[section_name]:
            
                alt_dict[section_name]['category_aggregate'] = config[section_name]['category_aggregate']
            if 'field_separators' in config[section_name]:
            
                alt_dict[section_name]['field_separators'] = config[section_name]['field_separators']
            if 'record_aggregate' in config[section_name]:
            
                alt_dict[section_name]['record_aggregate'] = config[section_name]['record_aggregate']
            
            for alternative in config[section_name]:
                
                order = alternative.split('-')
                
                for term in order:
                    
                    if term == 'score':
                        
                        score_key = order[1]  # Alias
                        score_dict[section_name][score_key] = config[section_name][alternative]
                        
                    if term == 'value':
                        
                        value_key = order[1]  # Alias
                        ## Split to seperate operation from number
                        operation_list = config[section_name][alternative].split(':')
                        
                        ## Check that a operation was given
                        if len(operation_list) > 1:
                            
                            operation_dict[section_name][value_key] = operation_list[0]
                            value_dict[section_name][value_key] = operation_list[1]
                        else:
                            
                            value_dict[section_name][value_key] = operation_list[0]
    # Return dicts with vcf_info_key as 1st key and alternative as 2nd key
    return alt_dict, score_dict, value_dict, operation_dict
