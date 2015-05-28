#!/usr/bin/env python
# encoding: utf-8
"""
print_variants.py

Print the variants of a file to vcf file.

There are two modes, 'vcf' or 'modified'.
If 'vcf' we expect plain vcf variants and print them as they came in.
If 'modified' the first column has been used for sorting so we skip that one.

If a outfile is provided the variants will be printed to this one.

Created by Måns Magnusson on 2015-01-22.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function, unicode_literals

from codecs import open


def print_variants(variant_file, outfile=None, mode='modified', silent=False):
    """
    Print the variants.
    
    If a result file is provided the variante will be appended to the file, 
    otherwise they are printed to stdout.
    
    There are two modes, 'vcf' or 'modified'.
    If 'vcf' we expect plain vcf variants and print them as they came in.
    If 'modified' the first column has been used for sorting so we skip 
    that one.
    
    Args:
        variants_file (str): A string with the path to a file
        outfile (str): Path to outfile or None
        mode (str): 'vcf' or 'modified'
        silent (bool): Bool. If nothing should be printed.
    
    """
    
    with open(variant_file, mode='r', encoding='utf-8') as f:
        if outfile:
            g = open(outfile, 'a', encoding='utf-8')
        
        for line in f:
            if not line.startswith('#'):
                line = line.rstrip().split('\t')
                if mode == 'modified':
                    line = line[1:]
                
                if outfile:
                    g.write('\t'.join(line)+'\n')
                
                else:
                    if not silent:
                        print('\t'.join(line))
    return

