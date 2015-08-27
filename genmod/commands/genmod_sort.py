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
import logging

import click


from codecs import open
from tempfile import NamedTemporaryFile

from genmod.utils import (print_variants, print_variant_for_sorting, 
                          sort_variants)

@click.command()
@click.argument('variant_file', 
                    nargs=1, 
                    type=click.File('rb'),
                    metavar='<vcf_file> or -'
)
@click.option('-o', '--outfile', 
                    type=click.Path(exists=False),
                    help='Specify the path to a file where results should be stored.'
)
@click.option('-f', '--family_id', 
                    help="Specify the family id for sorting. If no family id "\
                    "the first family found in annotation will be used."
)
@click.option('-v', '--verbose', 
                is_flag=True,
                help='Increase output verbosity.'
)
def sort(variant_file, outfile, family_id, verbose):
    """
    Sort a VCF file based on rank score.\n
    """    
    logger = logging.getLogger(__name__)
    #Save the variant lines for printing
    header_lines = []
    
    # Create a temporary variant file for sorting
    temp_file = NamedTemporaryFile(delete=False)
    temp_file.close()
    # Open the temp file with codecs
    temp_file_handle = open(
                                temp_file.name,
                                mode='w',
                                encoding='utf-8',
                                errors='replace'
                                )
    
    # Print the variants with rank score in first column
    for line in variant_file:
        line = line.rstrip()
        if line.startswith('#'):
            header_lines.append(line)
        else:
            print_variant_for_sorting(
                variant_line = line, 
                outfile = temp_file_handle,
                family_id = family_id
            )
    # The tempfile includes the unsorted variants
    temp_file_handle.close()
    
    # Sort the variants based on rank score
    sort_variants(
        infile = temp_file.name, 
        mode='rank'
    )
    
    # Print the headers
    if outfile:
        g = open(outfile, 'w', encoding='utf-8')
    
    for header in header_lines:
        if outfile:
            g.write(header + '\n')
        else:
            print(header)
    if outfile:
        g.close()
    
    # Print the variants
    print_variants(
        variant_file = temp_file.name,
        outfile = outfile,
        mode = 'modified'
    )
    
    logger.info("Removing temp file")
    os.remove(temp_file.name)
    logger.debug("Temp file removed")
    


if __name__ == '__main__':
    sort()