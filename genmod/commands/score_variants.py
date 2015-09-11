#!/usr/bin/env python
# encoding: utf-8
"""
score_variants.py

Script for scoring genetic variants in VCF files.

Created by Måns Magnusson on 2015-09-03.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function

import sys
import os
import click
import logging


from codecs import open
from datetime import datetime
from validate import ValidateError

from ped_parser import FamilyParser
from vcf_parser import VCFParser

from genmod.vcf_tools import (add_metadata, print_variant_dict, add_vcf_info)
from genmod.score_variants import (ConfigParser, score_variant)

from genmod import __version__

@click.command()
@click.argument('variant_file',
                nargs=1,
                type=click.Path(),
                metavar='<vcf_file> or -'
)
@click.option('-f', '--family_file',
                nargs=1, 
                type=click.File('r'),
                metavar='<ped_file>'
)
@click.option('-t' ,'--family_type', 
                type=click.Choice(['ped', 'alt', 'cmms', 'mip']), 
                default='ped',
                help='If the analysis use one of the known setups, please specify which one.'
)
@click.option('-s', '--silent',
                is_flag=True,
                help='Do not print the variants.'
)
@click.option('-o', '--outfile',
                type=click.Path(exists=False),
                help='Specify the path to a file where results should be stored.'
)
@click.option('-c', '--score_config',
              type=click.Path(exists=True),
              help="The plug-in config file(.ini)"
)
@click.option('-v', '--verbose',
              count=True,
              help='Increase output verbosity. If -vv all scores will be printed'
)

def score(family_file, variant_file, family_type,
          score_config, silent, outfile, verbose):
    """
    Score variants in a vcf file using Weighted Sum Model.
    The specific scores should be defined in a config file, see examples on 
    github.
    """
    from genmod import logger as root_logger
    from genmod.log import init_log, LEVELS
    loglevel = LEVELS.get(min(verbose,2), "WARNING")
    init_log(root_logger, loglevel=loglevel)
    
    logger = logging.getLogger(__name__)
    logger = logging.getLogger("genmod.commands.score")
    
    logger.info('Running GENMOD score, version: {0}'.format(__version__))
    
    ## Check the score config:
    if not score_config:
        logger.warning("Please provide a score config file.")
        sys.exit(1)
    
    logger.debug("Parsing config file")
    try:
        config_parser = ConfigParser(score_config)
    except ValidateError as e:
        logger.error("Something wrong in plugin file, please see log")
        logger.info("Exiting")
        sys.exit(1)

    logger.debug("Config parsed succesfully")
    
    if variant_file == '-':
        variant_parser = VCFParser(
            fsock = sys.stdin, 
            )
    else:
        variant_parser = VCFParser(
            infile = variant_file, 
            )
    
    head = variant_parser.metadata
    header_line = head.header
    
    add_metadata(
        head,
        'info',
        'RankScore',
        annotation_number='1', 
        entry_type='Integer', 
        description="Combined rank score for the variant in this family."'GeneticModels'
    )
    
    for variant in variant_parser:
        rank_score = score_variant(variant, config_parser)
        variant = add_vcf_info(
            keyword = 'RankScore', 
            variant_dict=variant,
            annotation=rank_score)
        
        print_variant_dict(variant, header_line)
    

if __name__ == '__main__':
    score()
