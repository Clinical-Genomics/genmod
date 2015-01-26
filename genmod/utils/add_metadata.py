#!/usr/bin/env python
# encoding: utf-8
"""
add_metadata.py

Add a new metadata string to a vcf header object.

Created by Måns Magnusson on 2015-01-22.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function, unicode_literals

import sys
import os
import click

from datetime import datetime

from genmod.utils import print_headers

def add_metadata(head, metadata_type, annotation_id, annotation_number='.',
                entry_type=None, description=None, version=None,
                command_line_string=''):
    """
    Add genmod metadata to the vcf header.
    
    Args:
        head            : A vcf header object
        annotate_models : Bool. If genetic models should be annotated
        vep             : Bool. If variants already annotated with vep the they
                         will not get any region annotation by genmod.
        cadd_annotation : Bool. If cadd scores are annoated.
        cadd_raw        : Bool. If raw cadd scores should be annotated
        thousand_g      : Bool. If thousand genome frequencies are annotated
        exac            : Bool. If exac frequencies are added
        command_line_string : A string with the command that envoked genmod
    
    """
    # Update INFO headers
    if metadata_type == 'info':
        head.add_info(
            annotation_id,
            annotation_number,
            entry_type,
            description
        )
    elif metadata_type == 'version':
        head.add_version_tracking(
                        annotation_id, 
                        version, 
                        datetime.now().strftime("%Y-%m-%d %H:%M"), 
                        command_line_string
                    )
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
    
    add_metadata(
        head,
        'info',
        'GeneticModels', 
        annotation_number='.', 
        entry_type='String', 
        description="':'-separated list of genetic models for this variant."
        
        )
    # Test if metadata was added properly.
    print_headers(head)

if __name__ == '__main__':
    cli()

