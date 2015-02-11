#!/usr/bin/env python
# coding: utf-8
"""
plugin_reader.py

Read a config file with plugins that define a scoring model.
Perform config format validation tests.


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


def LoggErrorandExist(text, key, config):
    """
    Log errors in config file and exist

    Args:
        text   (string) : Error message to log
        key    (unicode) : Section key
        config (object) : Config object

    Returns:
        None:
    """
    log.critical('Plugin.' +
                 config['Plugin']['name'] +
                 ': ' + text + '"' + key +
                 '" in config file'
                 + '\n'
                 )
    log.critical('Aborting ranking' + '\n')
    sys.exit()
    return


def check_vcf_config(section_name, not_info_list, my_vcf_parser, config):
    """
    Check that the alternatives in config is present within vcf INFO keys.

    Args:
        section_name  (unicode) : Section name
        not_info_list (list) : List of none vcf INFO keys
        my_vcf_parser (object) : vcf_parser object
        config        (object) : Config object

    Returns:
        None:
    """
    if section_name in not_info_list:
        pass
    elif section_name not in my_vcf_parser.metadata.info_dict:

        LoggErrorandExist("Could not find config key " +
                          " in VCF INFO header for key=",
                          section_name, config
                          )
    return


def checkMandatoryConfigAlternatives(section_name, config):
    """
    Check mandatory alternatives(s) in the config file

    Args:
        section_name (unicode) : Section name
        config       (object) : Config object

    Returns:
        None:
    """
    ## Create mandatory element list
    config_elements = ["data_type", "category",
                       "category_aggregate", "record_aggregate"]

    ## Check alternatives
    for alternative in config_elements:

        if alternative not in config[section_name]:

            LoggErrorandExist("Could not find mandatory key alternative=",
                              alternative, config)
    return


def CheckScoreValueExist(section_name, score_switch, value_switch, config):
    """"
    Check config section for score and value alterantives

    Args:
        section_name (unicode) : Section name
        score_switch (boolean) : Seen score alternative switch
        value_switch (boolean) : Seen value alternative switch
        config       (object) : Config object

    Returns:
        None:
    """
    if value_switch:

        LoggErrorandExist("Could not find value alternative in section",
                          section_name, config)
    if score_switch:

        LoggErrorandExist("Could not find score alternative in section",
                          section_name, config)
    return


def CheckScoreValuePairs(section_name, score_dict, value_dict, config):
    """
    Check that we have identical keys for score and value pairs

    Args:
        section_name (unicode) : Section name
        score_dict   (dict) : Score dictionnary
        value_dict   (dict) : Value dictionnary
        config       (object) : Config object

    Returns:
        None:
    """

    for key in score_dict[section_name].keys():

        if not key in value_dict[section_name]:

            LoggErrorandExist("Could not find corresponding score-key: '" +
                              key + "' with value-key in section",
                              section_name, config
                              )


def collectKeys(config_file, my_vcf_parser, verbose):
    """
    Collects database keys determined by Plugin and saves them in dict.

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
    ## Section that should not go through vcf data tests
    section_metadata_list = ["DEFAULT", "Meta_information", "Plugin"]
    ## Alternatives not found in vcf INFO
    not_info_list = ["FILTER"]

    ## Check that we have a plugin section
    if "Plugin" not in config:
        log.critical("Could not find mandatory Plugin section in config file")
        log.critical('Aborting ranking' + '\n')
        sys.exit()

    ## Collect keys for Plugin
    for section_name in config.sections():

        if verbose:
            log.info("Config alternative: " + str(config[section_name]))
            if 'version' in config[section_name]:
                log.info("Plugin version: " +
                         str(config[section_name]['version']))

        ## Only vcf data fields
        if section_name not in section_metadata_list:

            ## Enable test for score and value alternative existance downstream
            score_switch = True
            value_switch = True

            ## Test for mandatory config parameters
            checkMandatoryConfigAlternatives(section_name, config)

            ## Test for presence of vcf INFO keys
            check_vcf_config(section_name, not_info_list,
                             my_vcf_parser, config)

            alt_dict[section_name] = dict()
            score_dict[section_name] = dict()
            value_dict[section_name] = dict()
            operation_dict[section_name] = dict()

            ## Save vcf_info_key and mandatory general alternative
            alt_dict[section_name]['data_type'] = config[section_name]['data_type']
            alt_dict[section_name]['category'] = config[section_name]['category']
            alt_dict[section_name]['category_aggregate'] = config[section_name]['category_aggregate']
            alt_dict[section_name]['record_aggregate'] = config[section_name]['record_aggregate']

            ## Optional alternative
            if 'field_separators' in config[section_name]:

                alt_dict[section_name]['field_separators'] = config[section_name]['field_separators']

            for alternative in config[section_name]:

                order = alternative.split('-')

                for term in order:

                    if term == 'score':

                        score_switch = False
                        score_key = order[1]  # Alias
                        score_dict[section_name][score_key] = config[section_name][alternative]

                    if term == 'value':

                        value_switch = False
                        value_key = order[1]  # Alias
                        ## Split to seperate operation from number
                        operation_list = config[section_name][alternative].split(':')

                        ## Check that a operation was given
                        if len(operation_list) > 1:

                            operation_dict[section_name][value_key] = operation_list[0]
                            value_dict[section_name][value_key] = operation_list[1]
                        else:

                            value_dict[section_name][value_key] = operation_list[0]
            CheckScoreValueExist(section_name, score_switch,
                                 value_switch, config)
            CheckScoreValuePairs(section_name, score_dict, value_dict, config)
    # Return dicts with vcf_info_key as 1st key and alternative as 2nd key
    return alt_dict, score_dict, value_dict, operation_dict
