#!/usr/bin/env python
# encoding: utf-8
"""
config_parser.py

Read a config file with rules for how to score variants

[Version]
  name = config_name # This is a string
  version = config_version # This is a float

Each plugin section will look like:

[Plugin_name]
  section = section_name # str in ['CHROM','POS','ID','REF','ALT', 'FILTER',
                          'QUAL', 'FILTER','INFO','FORMAT','sample_id']
  data_type = data_type # str in ['integer','float','flag','string']
  record_rule = record_rule # str in ['min', 'max']


Created by Måns Magnusson on 2015-04-16.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function

import logging
import click
import configobj

from six import string_types
from validate import ValidateError

from extract_vcf import Plugin

from genmod.score_variants import ScoreFunction

class ConfigParser(configobj.ConfigObj):
    """
    Class for holding information from config file.
    
    """
    def __init__(self, config_file, indent_type='  ', encoding='utf-8'):
        super(ConfigParser, self).__init__(
                                        infile=config_file, 
                                        indent_type=indent_type, 
                                        encoding=encoding,
                                        )
        self.logger = logging.getLogger(__name__)
        self.logger = logging.getLogger("genmod.score_varaints.config_parser")
        
        self.categories = {}
        self.plugins = {}
        self.score_functions = {}
        
        self.vcf_columns = ['CHROM','POS','ID','REF','ALT', 'FILTER','QUAL',
                            'FILTER','INFO','FORMAT','sample_id']
        
        self.data_types = ['integer','float','flag','character','string']
        
        self.logger.info("Checking version and name")
        self.version_check()
        self.version = float(self['Version']['version'])
        self.logger.debug("Set version to {0}".format(self.version))
        
        self.name = self['Version']['name']
        self.logger.debug("Set name to {0}".format(self.name))

        self.logger.info("Config name: {0}".format(self['Version']['name']))
        self.logger.info("Config version: {0}".format(self['Version']['version']))
        
        for plugin in self.keys():
            if plugin != 'Version' and plugin != 'Categories':
                self.plugins[plugin] = None
                self.score_functions[plugin] = None
        
        self.logger.info("Found plugins: {0}".format(
            ', '.join(list(self.plugins.keys()))))
        
        self.logger.info("Checking categories")
        for category in self['Categories']:
            self.logger.debug("Found category {0}".format(category))
            self.categories[category] = {}
            aggregation = self['Categories'][category].get('category_aggregation', 'min')
            self.logger.debug("Setting aggregation to {0}".format(aggregation))
            self.categories[category]['category_aggregation'] = aggregation
            self.categories[category]['plugins'] = []
            
        self.logger.info("Checking plugins")
        
        for plugin in self.keys():
            if plugin != 'Version' and plugin != 'Categories':
                self.logger.debug("Checking plugin: {0}".format(plugin))
                self.check_plugin(plugin)
                self.logger.debug("Plugin {0} is ok".format(plugin))
                plugin_info = self[plugin]
                
                string_rules = {}
                if plugin_info['data_type'] == 'string':
                    self.logger.info("Checking string rules for plugin {0}".format(
                        plugin
                    ))
                    string_rules = self.get_string_dict(plugin_info)
                
                self.logger.debug("Adding plugin {0} to ConfigParser".format(plugin))
                
                category = plugin_info.get('category', None)
                
                self.plugins[plugin] = Plugin(
                    name=plugin, 
                    field=plugin_info['field'], 
                    data_type=plugin_info['data_type'], 
                    separators=plugin_info.get('separators',[]), 
                    info_key=plugin_info.get('info_key',None),
                    category=category,
                    csq_key=plugin_info.get('csq_key', None),
                    record_rule=plugin_info.get('record_rule', 'max'),
                    string_rules=string_rules
                )
                
                if category:
                    if category in self.categories:
                        self.categories[category]['plugins'].append(plugin)
                        self.logger.debug("Adding {0} to category {1}".format(
                            plugin, category))
                    else:
                        self.logger.error("{0} have an undefined category {1}".format(
                            plugin, category))
                        self.logger.info("Please add specifications for category {0}".format(
                            category
                        ))
                        raise ValidateError("Plugins must have a category"\
                                        " defined in 'Categories' section")
                else:
                    self.logger.error("{0} is missing category".format(plugin))
                    raise ValidateError("Plugins must have a category")
                
                self.logger.info("Check score function for plugin {0}".format(plugin))
                self.score_functions[plugin] = self.get_score_function(plugin_info)
                self.logger.debug("Added score function for plugin {0}".format(plugin))

                    
    #
    #     logger.info("Checking plugin scores")
    #     for plugin in self.plugins:
    #         logger.debug("Checking plugin score: {0}".format(plugin))
    #         self[plugin] = self.vcf_score_check(self[plugin], plugin)
    #
    def get_score_function(self, plugin_info):
        """Convert the scoring information
        
        If data_type = String we use the string dict to get the score
        If data_type = Flag we only need to check if we have an annotation
        If data_type = Float or Int we can build a interval tree with score as data
        
        Arguments:
            plugin_info (dict): A dictionary with plugin information
        
        """
        
        score_dict = {}
        
        for key in plugin_info:
            try:
                score_dict[key] = dict(plugin_info[key])
            except ValueError:
                pass

        match_type = plugin_info['data_type']
        
        score_function = ScoreFunction(match_type=match_type)
        
        for score_name in score_dict:
            raw_info = score_dict[score_name]
            
            if score_name == 'not_reported':
                not_reported_score = float(raw_info['score'])
                score_function.set_not_reported(not_reported_score)
            
            elif match_type == 'flag':
                reported_score = float(raw_info['score'])
                score_function.set_reported(reported_score)
            
            elif match_type == 'string':
                string = raw_info['string']
                score = float(raw_info['score'])
                score_function.add_string_rule(key=string, score=score)
            
            else:
                if raw_info['score'] == 'eq':
                    score_function.set_equal()
                elif raw_info.get('value'):
                    score_function.add_value(
                        value = raw_info['value'],
                        score = raw_info['score']
                    )
                else:
                    lower_bound = float(raw_info['lower'])
                    upper_bound = float(raw_info['upper'])
                    score = float(raw_info['score'])
                    score_function.add_interval(
                        lower=lower_bound,
                        upper=upper_bound,
                        score=score
                    )
        
        return score_function
                
        
    def get_string_dict(self, plugin_info):
        """Convert a section with information of priorities to a string dict.
        
        Arguments:
            plugin_info (dict): A dictionary with plugin information
        
        Return:
            string_dict (dict): A dictionary with strings as keys and integer
                                that specifies their priorities as values
        """
        string_info = []
        string_dict = {}
        
        for key in plugin_info:
            try:
                if key != 'not_reported':
                    string_info.append(dict(plugin_info[key]))
            except ValueError:
                pass
        
        string_rules = {}
        
        for raw_info in string_info:
            try:
                string = raw_info['string']
            except KeyError:
                
                raise ValidateError("String information has to have a 'string'")
            try:
                priority = raw_info['priority']
            except KeyError:
                raise ValidateError("String information has to have a 'priority'")
            try:
                priority = int(priority)
            except ValueError:
                raise ValidateError("'priority' has to be an integer")
            
            string_dict[string] = priority

        if len(string_dict) == 0:
            raise ValidateError("'string' entrys must have string rules defined")
            
        return string_dict
            
    
    
    def version_check(self):
        """
        Check if the version entry is in the proper format
        """
        try:
            version_info = self['Version']
        except KeyError:
            raise ValidateError('Config file has to have a Version section')
        try:
            float(version_info['version'])
        except KeyError:
            raise ValidateError('Config file has to have a version section')
        except ValueError:
            raise ValidateError('Version has to be a float.')
        
        try:
            version_info['name']
        except KeyError:
            raise ValidateError("Config file has to have a name")
        return
   
    def check_plugin(self, plugin):
        """
        Check if the section is in the proper format vcf format.

        Args:
            vcf_section (dict): The information from a vcf section

        Returns:
            True is it is in the proper format

        """
        
        vcf_section = self[plugin]
        
        try:
            vcf_field = vcf_section['field']
            if not  vcf_field in self.vcf_columns:
                raise ValidateError(
                        "field has to be in {0}\n"
                        "Wrong field name in plugin: {1}".format(
                        self.vcf_columns, plugin
                    ))
            if vcf_field == 'INFO':
                try:
                    info_key = vcf_section['info_key']

                    if info_key == 'CSQ':
                        try:
                            csq_key = vcf_section['csq_key']
                        except KeyError:
                            try:
                                csq_key = vcf_section['vep_key']
                            except KeyError:
                                raise ValidateError(
                                "CSQ entrys has to refer to an csq field.\n"
                                "Refer with keyword 'csq_key'\n"
                                "csq_key is missing in section: {0}".format(
                                plugin
                                )
                            )


                except KeyError:
                    raise ValidateError(
                        "INFO entrys has to refer to an INFO field.\n"
                        "Refer with keyword 'info_key'\n"
                        "info_key is missing in section: {0}".format(
                            plugin
                            )
                        )
        except KeyError:
            raise ValidateError(
                "Vcf entrys have to refer to a field in the VCF with keyword"
                " 'field'.\nMissing keyword 'field' in plugin: {0}".format(
                  plugin
                ))

        try:
            data_type = vcf_section['data_type']
            if not data_type in self.data_types:
                raise ValidateError(
                    "data_type has to be in {0}\n"
                    "Wrong data_type in plugin: {1}".format(
                        self.data_types, plugin)
                    )
        except KeyError:
            raise ValidateError(
                "Vcf entrys have to refer to a data type in the VCF with "
                "keyword 'data_type'.\n"
                "Missing data_type in plugin: {0}".format(plugin)
                )

        
        separators = vcf_section.get('separators', None)
        
        if separators:
            if len(separators) == 1:
                self[plugin]['separators'] = list(separators)
        else:
            if data_type != 'flag':
                raise ValidateError(
                    "If data_type != flag the separators have to be defined.\n"
                    "Missing separators in plugin: {0}".format(plugin)
                    )
                
        
        record_rule = vcf_section.get('record_rule', None)
        
        if record_rule:
            if not record_rule in ['min', 'max']:
                raise ValidateError(
                    "Record rules have to be in {0}\n"
                    "Wrong record_rule in plugin: {1}".format(
                        ['min', 'max'], plugin)
                )
        else:
            self.logger.info("Setting record rule to default: 'max'")
                
        return True

@click.command()
@click.argument('config_file',
                nargs=1,
                type=click.Path(exists=True)
)
@click.option('-out', '--outfile',
                nargs=1,
                type=click.File('w')
)
@click.option('-l', '--loglevel',
                type=click.Choice(['DEBUG', 'INFO', 'WARNING']),
                default = 'INFO'
)
def read_config(config_file, outfile, loglevel):
    """Parse the config file and print it to the output."""

    from genmod import logger
    from genmod.log import init_log
    init_log(logger, loglevel='DEBUG')

    logger.info("Reading Config File: {0}".format(config_file))

    config_reader = ConfigParser(config_file)

    for plugin in config_reader.plugins:
        logger.info("Found plugin:{0}".format(plugin))
        logger.info("{0}: {1}".format(
            plugin,config_reader.plugins[plugin])
            )

    for category in config_reader.categories:
        logger.info("Category {0}: {1}".format(
            category, config_reader.categories[category]
        ))


if __name__ == '__main__':
    read_config()
