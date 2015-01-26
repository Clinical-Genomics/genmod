#!/usr/bin/env python
# encoding: utf-8
"""
print_headers.py

Print the lines in a header object to a file or stdout.

Created by Måns Magnusson on 2015-01-22.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function, unicode_literals

import sys
import os
import click

from codecs import open

def print_headers(head, outfile=None, silent=False):
    """
    Print the vcf headers.
    
    If a result file is provided headers will be printed here, otherwise
    they are printed to stdout.
    
    Args:
        head    : A vcf header object
        outfile : Path to outfile or None
        silent  : Bool. If nothing should be printed.
        
    """
    if outfile:
        with open(outfile, 'w', encoding='utf-8') as f:
            for head_count in head.print_header():
                f.write(head_count+'\n')
    else:
        if not silent:
            for line in head.print_header():
                print(line)
    return

@click.command()
@click.argument('vcf_file', 
                nargs=1, 
                type=click.Path(),
                metavar='<vcf_file> or -'
)
@click.option('-s' ,'--split_variants', 
                    is_flag=True,
                    help='If the variants should be splitted.'
)
@click.option('-o' ,'--outfile', 
                type=click.Path(exists=False),
                help='Print to output.'
)
@click.option('--silent',
              is_flag=True,
              help='Only error messages printed.'
)
@click.option('-v', '--verbose',
              is_flag=True,
              help='Increase output verbosity.'
)
def cli(vcf_file, split_variants, outfile, silent, verbose):
    from vcf_parser import VCFParser
    
    if vcf_file == '-':
        variant_parser = VCFParser(
                            fsock = sys.stdin, 
                            split_variants=split_variants
                        )
    else:
        variant_parser = VCFParser(
                            infile = vcf_file, 
                            split_variants=split_variants
                        )
    
    head = variant_parser.metadata
    
    print_headers(head, outfile, silent)

if __name__ == '__main__':
    cli()

