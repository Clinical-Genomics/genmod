#!/usr/bin/env python
# encoding: utf-8
"""
filter_variants.py

Command line tool for filtering vcf variants.
Filter variants based on some annotation.

Created by Måns Magnusson on 2015-08-25.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import (print_function)

import sys
import os
import logging
import itertools

import click

from codecs import open
from datetime import datetime

from genmod import __version__

from genmod.vcf_tools import (HeaderParser, get_variant_dict, get_info_dict,
print_variant, print_headers)
from extract_vcf import Plugin

@click.command()
@click.argument('variant_file',
                    nargs=1,
                    type=click.File('r'),
                    metavar="<vcf_file> or -"
)
@click.option('-a', '--annotation', 
                default='1000GAF',
                help="Specify the info annotation to search for."\
                " Default 1000GAF"
)
@click.option('-t', '--threshold',
                    default=0.05, 
                    help="""Threshold for filter variants. Default 0.05"""
)
@click.option('-d', '--discard',
                    is_flag=True, 
                    help="If variants without the annotation should be"\
                    " discarded"
)
@click.option('-g', '--greater',
                    is_flag=True, 
                    help="If greater than threshold should be used instead of"\
                    " less thatn threshold."
)
@click.option('-s', '--silent',
                is_flag=True,
                help='Do not print the variants.'
)
@click.option('-o', '--outfile', 
                    type=click.File('w'),
                    help="Specify the path to a file where results should be stored."
)
def filter(variant_file, annotation, threshold, discard, greater, silent, outfile):
    """
    Filter vcf variants.
    
    Filter variants based on their annotation
    """

    logger = logging.getLogger(__name__)
    #For testing
    logger = logging.getLogger("genmod.commands.filter")
    
    logger.info("Running genmod filter version {0}".format(__version__))
    
    start_time_analysis = datetime.now()
    
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
    
    if not annotation in head.info_dict:
        logger.warning("Annotation {0} not specified in header".format(annotation))
        logger.info("Please check VCF file")
        logger.info("Exiting...")
        sys.exit(1)
    
    logger.info("Building a plugin from extract_vcf for {0}".format(annotation))
    annotation_plugin = Plugin(
        name=annotation, 
        field='INFO',
        info_key=annotation,
        separators = [','],
        record_rule = 'min',
        data_type = 'float'
    )
    logger.debug("Plugin=(field={0},info_key={1},separators={2},record_rule={3}"\
    ",data_type={4})".format('INFO', annotation, "','", 'min', 'float'))
    
    print_headers(head=head, outfile=outfile, silent=silent)
    
    nr_of_variants = 0
    nr_of_passed_variants = 0
    for variant in variant_file:
        nr_of_variants += 1
        keep_variant = False
        value = annotation_plugin.get_value(variant_line=variant)
        logger.debug("Found value {0}".format(value))
        if value:
            if greater:
                if value > threshold:
                    keep_variant = True
            else:
                if value < threshold:
                    keep_variant = True
        else:
            if not discard:
                keep_variant = True
        
        if keep_variant:
            logger.debug("Keeping variant")
            nr_of_passed_variants += 1
            print_variant(
                variant_line=variant, 
                outfile=outfile, 
                mode='vcf', 
                silent=silent
            )
        else:
            logger.debug("Discarding variant")
            

    logger.info("Number of variants in file {0}".format(nr_of_variants))
    logger.info("Number of variants passing filter {0}".format(nr_of_passed_variants))
    logger.info("Number of variants filtered {0}".format(
        nr_of_variants - nr_of_passed_variants))

if __name__ == '__main__':
    from genmod.log import init_log
    from genmod import logger
    init_log(logger, loglevel="DEBUG")
    filter()
