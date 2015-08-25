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


from codecs import open
from datetime import datetime
from tempfile import mkdtemp, TemporaryFile, NamedTemporaryFile
from pprint import pprint as pp

import shutil
import pkg_resources

# from genmod import (sort_variants, print_headers)
from genmod.utils import (print_variant_for_sorting, sort_variants, 
                          print_variants)



@click.command()
@click.argument('variant_file', 
                    nargs=1, 
                    type=click.File('rb'),
                    metavar='<vcf_file> or -'
)
@click.option('-o', '--outfile', 
                    type=click.File('w'),
                    help='Specify the path to a file where results should be stored.'
)
@click.option('-f', '--family_id', 
                    help='Specify the family id for sorting.'
)
@click.option('-v', '--verbose', 
                is_flag=True,
                help='Increase output verbosity.'
)
def sort(variant_file, outfile, family_id, verbose):
    """
    Sort a VCF file.\n
    """    
    
    header = []
    
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
        if line.startswith('#'):
            line = line.rstrip()
        #     print(line)
            header.append(line)
        else:
            print_variant_for_sorting(
                variant_line = line, 
                outfile = temp_file_handle,
                family_id = family_id
            )
    temp_file_handle.close()
    
    # Sort the variants based on rank score
    sort_variants(
        infile = temp_file.name, 
        mode='rank'
    )
    
    # Print the headers
    for line in header:
        if outfile:
            outfile.write(line+'\n')
        else:
            print(line)
    
    # Print the variants
    print_variants(
        variant_file = temp_file.name,
        outfile = outfile
    )



if __name__ == '__main__':
    sort()