#!/usr/bin/env python
# encoding: utf-8
"""
genmod_sort.py

Command line tool for sorting variants with genmod.

Created by Måns Magnusson on 2015-08-21.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function

import sys
import os
import click
import logging

from codecs import open
from tempfile import NamedTemporaryFile
from datetime import datetime


from genmod.vcf_tools import (print_variant_for_sorting, sort_variants, 
get_info_dict, print_variant, HeaderParser, print_headers)

from genmod.utils import (get_chromosome_priority, get_rank_score)

from genmod import __version__

from .utils import (variant_file, outfile, silent, temp_dir)

@click.command()
@variant_file
@outfile
@click.option('-f', '--family_id',
                    type=str,
                    help='Specify the family id for sorting.'
)
@silent
@temp_dir
@click.option('-p', '--position',
                is_flag=True,
                help='If variants should be sorted by position.'
)
def sort(variant_file, outfile, family_id, silent, position, temp_dir):
    """
    Sort a VCF file based on rank score.
    """    
    logger = logging.getLogger(__name__)
    head = HeaderParser()

    logger.info("Running GENMOD sort version {0}".format(__version__))
    start = datetime.now()
    # Create a temporary variant file for sorting
    logger.debug("Creating temporary file for sorting")
    if temp_dir:
        temp_file = NamedTemporaryFile(delete=False, dir=temp_dir)
    else:
        temp_file = NamedTemporaryFile(delete=False)
    temp_file.close()
    # Open the temp file with codecs
    temp_file_handle = open(
                                temp_file.name,
                                mode='w',
                                encoding='utf-8',
                                errors='replace'
                                )
    logger.debug("Temp file created")
    logger.info("Printing variants to temp file")
    nr_variants = 0
    # Print the variants with rank score in first column
    for line in variant_file:
        line = line.rstrip()
        if line.startswith('#'):
            if line.startswith('##'):
                head.parse_meta_data(line)
            else:
                head.parse_header_line(line)
        else:
            nr_variants += 1
            priority = '0'
            
            if position:
                chrom = line.split()[0]
                priority = get_chromosome_priority(chrom)
            else:
                priority = get_rank_score(line)
            
            print_variant(
                variant_line=line, 
                priority=priority, 
                outfile=temp_file_handle
            )
    
    temp_file_handle.close()
    
    logger.info("Variants printed to temp file")
    logger.info("Nr or variants in VCF file: {0}".format(nr_variants))
    
    sort_mode = 'rank'
    
    if position:
        sort_mode = 'chromosome'
    
    logger.info("Sorting variants")
    sort_variants(
        infile = temp_file.name, 
        mode=sort_mode
    )
    logger.info("Variants sorted")

    logger.debug("Printing headers")
    print_headers(
        head = head, 
        outfile = outfile, 
        silent=silent
    )
    logger.debug("Headers printed")
    
    logger.info("Printing variants")
    with open(temp_file.name, mode='r', encoding='utf-8', errors='replace') as f:
        for variant_line in f:
            print_variant(
                variant_line = variant_line, 
                outfile = outfile, 
                mode = 'modified',
                silent=False
                )
    logger.debug("Variants printed")
    
    logger.info("Removing temp file")
    os.remove(temp_file.name)
    logger.debug("Temp file removed")
    
    logger.info("Sorting done, time for sorting: {0}".format(datetime.now()-start))


if __name__ == '__main__':
    sort()