#!/usr/bin/env python
# encoding: utf-8
"""
analyze.py

Command line tool for analyzing variants that are annotated with genmod.

Created by Måns Magnusson on 2014-09-03.
Copyright (c) 2014 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function

import logging
import itertools

from datetime import datetime
from pprint import pprint as pp

import click

from genmod.vcf_tools import (HeaderParser, get_variant_dict, 
                              get_info_dict, get_genetic_models)

from genmod.commands.utils import (variant_file, get_file_handle)

logger = logging.getLogger(__name__)

GENETIC_MODELS = ['AD','AD_dn','AR','AR_dn','AR_comp','AR_comp_dn','XD',
                  'XD_dn','XR','XR_dn']

@click.command()
@variant_file
@click.option('--gq_treshold', '-gq',
                    default=20, 
                    nargs=1,
                    help="""Specify genotype quality treshold for variants 
                            to be considered. Default 20."""
)
@click.option('-v', '--verbose', 
                is_flag=True,
                help='Increase output verbosity.'
)
def analyse(variant_file, gq_treshold, verbose):
    """Analyse variants annotated by genmod in a VCF file. 

    """

    start_time_analysis = datetime.now()
    variants = get_file_handle(variant_file)
    
    # Dictionary where we collect information about the analysis
    # per family
    families = {}

    logger.info("Initializing a Header Parser")
    head = HeaderParser()

    line = None
    for line in variants:
        line = line.rstrip()

        if line.startswith('#'):
            if line.startswith('##'):
                head.parse_meta_data(line)
            else:
                head.parse_header_line(line)
        else:
            break
    
    #Add the first variant back to the iterator
    # If the vcf has no variants the last line will be a header
    if not line.startswith('#'):
        variants = itertools.chain([line], variants)
    else:
        print_headers(head, outfile, silent)
        sys.exit(0)
    
    header_line = head.header
    
    for variant in variants:
        variant_dict = get_variant_dict(variant, header_line)
        info = get_info_dict(variant_dict['INFO'])
        genetic_models = get_genetic_models(info)
        for family_id in genetic_models:
            for genetic_model in genetic_models[family_id]:
                if family_id in families:
                    if genetic_model in families[family_id]:
                        families[family_id][genetic_model] += 1
                    else:
                        families[family_id][genetic_model] = 1
                else:
                    families[family_id] = {}
                    families[family_id][genetic_model] = 1
    
    click.echo('\t'.join(['family_id']+GENETIC_MODELS))
    for family_id in families:
        out_line = [family_id]
        for model in GENETIC_MODELS:
            out_line.append(str(families[family_id].get(model, 0)))
        click.echo('\t'.join(out_line))

if __name__ == '__main__':
    analyse()