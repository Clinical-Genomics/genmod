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
import itertools

from codecs import open
from datetime import datetime
from validate import ValidateError

from ped_parser import FamilyParser

from genmod.vcf_tools import (add_metadata, print_variant, add_vcf_info,
print_headers, HeaderParser, get_variant_dict, get_info_dict)
from genmod.score_variants import (ConfigParser, score_variant)

from genmod import __version__

@click.command()
@click.argument('variant_file',
                nargs=1,
                type=click.File('r'),
                metavar='<vcf_file> or -'
)
@click.option('-i   ', '--family_id',
                default='1', 
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
                type=click.File('w'),
                help='Specify the path to a file where results should be stored.'
)
@click.option('-c', '--score_config',
              type=click.Path(exists=True),
              help="The plug-in config file(.ini)"
)
def score(variant_file, family_id, family_file, family_type, score_config, 
silent, outfile):
    """
    Score variants in a vcf file using Weighted Sum Model.
    
    The specific scores should be defined in a config file, see examples on 
    github.
    """
    # from genmod import logger as root_logger
    # from genmod.log import init_log, LEVELS
    # loglevel = LEVELS.get(min(verbose,2), "WARNING")
    # init_log(root_logger, loglevel=loglevel)
    
    logger = logging.getLogger(__name__)
    # logger = logging.getLogger("genmod.commands.score")
    
    logger.info('Running GENMOD score, version: {0}'.format(__version__))
    
    logger.info("Checking family id")
    
    if family_file:
        logger.info("Setting up a family parser")
        family_parser = FamilyParser(family_file, family_type)
        logger.debug("Family parser done")
        family_id = list(family_parser.families.keys())[0]
    
    logger.info("Family used in analysis: {0}".format(family_id))
    
    
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

    logger.info("Initializing a Header Parser")
    head = HeaderParser()

    for line in variant_file:
        line = line.rstrip()
        if line.startswith('#'):
            if line.startswith('##'):
                head.parse_meta_data(line)
            else:
                head.parse_header_line(line)
        else:
            break
    
    #Add the first variant to the iterator
    variant_file = itertools.chain([line], variant_file)
    header_line = head.header
    
    if "RankScore" in head.info_dict:
        logger.warning("Variants already scored according to VCF header")
        logger.info("Please check VCF file")
        logger.info("Exiting...")
        sys.exit(1)
    
    add_metadata(
        head,
        'info',
        'RankScore',
        annotation_number='.', 
        entry_type='String', 
        description="The rank score for this variant in this family. family_id:rank_score."
    )
    print_headers(
        head=head,
        outfile=outfile,
        silent=silent
    )
    start_scoring = datetime.now()
    last_twenty = datetime.now()
    nr_of_variants = 1

    for line in variant_file:
        if not line.startswith('#'):
            variant = get_variant_dict(line, header_line)
            variant['info_dict'] = get_info_dict(variant['INFO'])

            rank_score = score_variant(variant, config_parser)
            
            variant = add_vcf_info(
                keyword = 'RankScore',
                variant_dict=variant,
                annotation="{0}:{1}".format(family_id, rank_score)
            )

            print_variant(
                variant_dict=variant,
                header_line=header_line,
                outfile=outfile,
                silent=silent
            )

            nr_of_variants += 1

            if nr_of_variants % 20000 == 0:
                logger.info("{0} variants scored.".format(nr_of_variants))
                logger.info("Last 20000 took {0} to score.".format(datetime.now()-last_twenty))
                last_twenty = datetime.now()

    logger.info("Variants scored. Number of variants: {0}".format(nr_of_variants))
    logger.info("Time to score variants: {0}".format(datetime.now()-start_scoring))
        

if __name__ == '__main__':
    score()
