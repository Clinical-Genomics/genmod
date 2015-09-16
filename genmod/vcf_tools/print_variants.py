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

from __future__ import print_function

from codecs import open


def print_variant_dict(variant, header_line, outfile=None, silent=False):
    """Print a variant dictionary
    
        Prints a variant dictionary to a file or stdout
        
        Args:
            variant (dict): A variant dictionary
            header_line (list): A list with the header columns
            outfile (FileHandle): A file handle
        
    """
    print_line = [variant.get(entry, '.') for entry in header_line]
    
    if outfile:
        outfile.write('\t'.join(print_line) + '\n')
    else:
        if not silent:
            print('\t'.join(print_line))
    
def print_variant(variant_line, outfile=None, mode='vcf', silent=False):
    """
    Print a variant line.
    
    If a result file is provided the variante will be appended to the file, 
    otherwise they are printed to stdout.
    
    There are two modes, 'vcf' or 'modified'.
    If 'vcf' we expect plain vcf variants and print them as they came in.
    If 'modified' the first column has been used for sorting so we skip 
    that one.
    
    Args:
        variants_file (str): A string with the path to a file
        outfile (FileHandle): An opened file_handle
        mode (str): 'vcf' or 'modified'
        silent (bool): Bool. If nothing should be printed.
    
    """
    
    if not variant_line.startswith('#'):
        splitted_line = variant_line.rstrip().split('\t')
        if mode == 'modified':
            splitted_line = splitted_line[1:]
            
        if outfile:
            outfile.write('\t'.join(splitted_line)+'\n')
        
        else:
            if not silent:
                print('\t'.join(splitted_line))
    return

def print_variant_for_sorting(variant_line, priority, outfile):
    """
    Print the variants for sorting
    
    Arguments:
        variant_line (str): A vcf variant line
        prority (str): The priotiy for this variant
        outfile (file_handle): A filehandle to the temporary variant file
        family_id (str): The family Id for sorting on rank score
    """
    variant_line = variant_line.rstrip().split("\t")
    
    outfile.write("{0}\t{1}\n".format(priority, '\t'.join(variant_line)))
                    
