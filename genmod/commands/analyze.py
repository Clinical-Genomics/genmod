#!/usr/bin/env python
# encoding: utf-8
"""
analyze.py

Command line tool for analyzing variants that are annotated with genmod.

Created by Måns Magnusson on 2014-09-03.
Copyright (c) 2014 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function
from __future__ import unicode_literals

import sys
import os
import click
import inspect

try:
    import cPickle as pickle
except:
    import pickle

from codecs import open
from datetime import datetime
from tempfile import NamedTemporaryFile
from pprint import pprint as pp
from configobj import ConfigObj

import pkg_resources

from vcf_parser import parser as vcf_parser

import genmod
from genmod import variant_sorter, variant_printer, warning


def print_headers(head, outfile=None, silent=False):
    """Print the headers to a results file."""
    if outfile:
        with open(outfile, 'w', encoding='utf-8') as f:
            for head_count in head.print_header():
                f.write(head_count+'\n')
    else:
        if not silent:
            for line in head.print_header():
                print(line)
    return


def print_variants(sorted_variants, outfile, silent=False):
    """Print the variants to a results file or stdout."""
    
    with open(sorted_variants, mode='r', encoding='utf-8') as f:
        if outfile:
            with open(outfile, 'a', encoding='utf-8') as g:
                for variant in f:
                    g.write(variant)
        else:
            if not silent:
                for line in f:
                    print(line.rstrip())
    return



###           This is for analyzing the variants       ###


def make_models(list_of_models):
    """Make a dictionary of the prefered models.
        If no models are specified all are considered interesting."""
    model_set = set()
    # If no models are specified we allow all models
    if len(list_of_models) == 0:
        list_of_models = ['AR', 'AD', 'X']
    
    for model in list_of_models:
        if 'AR' in model:
            model_set.add('AR_hom')
            model_set.add('AR_hom_dn')
            model_set.add('AR_comp')
            model_set.add('AR_comp_dn')
        if 'AD' in model:
            model_set.add('AD')
            model_set.add('AD_dn')
        if 'X' in model:
            model_set.add('XR')
            model_set.add('XR_dn')
            model_set.add('XD')
            model_set.add('XD_dn')
    return model_set

@click.command()
@click.argument('variant_file',
                    nargs=1,
                    type=click.Path(exists=True),
                    metavar='<vcf_file> or "-"'
)
@click.option('-c', '--config_file',
                    type=click.Path(exists=True),
                    help="""Specify the path to a config file."""
)
@click.option('--frequency', '-freq',
                    type=float, 
                    nargs=1,
                    help='Specify the treshold for variants to be considered. Default 0.05'
)
@click.option('-p', '--patterns',
                    type=click.Choice(['AR', 'AD', 'X']),
                    multiple=True,
                    help='Specify the inheritance patterns. Default is all patterns'
)
@click.option('-o', '--outfile', 
                    type=click.Path(exists=False),
                    help='Specify the path to a file where results should be stored.'
)
@click.option('-s', '--silent', 
                is_flag=True,
                help='Do not output variants.'
)
@click.option('-v', '--verbose', 
                is_flag=True,
                help='Increase output verbosity.'
)
def analyze(variant_file, frequency, patterns, config_file, outfile, silent,verbose):
    """Analyze the annotated variants in a VCF file."""
    configs = ConfigObj(config_file)
    
    freq_treshold = 0.05
    freq_keyword = '1000GMAF'
    prefered_models = make_models([])
    inheritance_keyword = 'GeneticModels'
    if config_file:
        freq_treshold = float(configs.get('frequency', {}).get('rare', freq_treshold))
        freq_keyword = configs.get('frequency', {}).get('keyword', freq_keyword)    
        inheritance_patterns = [pattern for pattern in configs.get('inheritance', {}).get('patterns',[])]
        inheritance_keyword = configs.get('inheritance', {}).get('keyword',inheritance_keyword)
        prefered_models = make_models(inheritance_patterns)
    
    # Command line overrides config file
    if frequency:
        freq_treshold = frequency
    if patterns:
        prefered_models = make_models(patterns)
    
    if variant_file == '-':
        variant_parser = vcf_parser.VCFParser(fsock = sys.stdin)
    else:
        variant_parser = vcf_parser.VCFParser(infile = variant_file)
    
    dominant_file = NamedTemporaryFile(delete=False)
    dominant_file.close()
    de_novo_file = NamedTemporaryFile(delete=False)
    de_novo_file.close()
    
    with open(dominant_file.name, mode='w', encoding = 'utf-8') as f:
        
        for variant in variant_parser:
            models_found = set(variant['info_dict'].get(inheritance_keyword, '').split(','))
            # print(models_found)
            
            maf = min([float(frequency) for frequency in variant['info_dict'].get(freq_keyword, '0').split(',')])
            cadd_score = max([float(cscore) for cscore in variant['info_dict'].get('CADD', '0').split(',')])
            
            # Check if variant models overlap with prefered:
            if models_found.intersection(set(['AD'])):
                # Check if cadd score is available:
                # print('hej')
                # pp(variant)
                if cadd_score > 0:
                    # Check if MAF is below treshold:
                    if maf < freq_treshold:
                        print_line = [variant.get(entry, '-') for entry in variant_parser.header]
                        f.write('\t'.join(print_line)+'\n')
    
    # print_headers(variant_parser.metadata, outfile=outfile)
    
    dominant_results = NamedTemporaryFile(delete=False)
    dominant_results.close()
    
    var_sorter = variant_sorter.FileSort(dominant_file.name, mode='cadd', outfile=dominant_results.name)
    var_sorter.sort()
    
    print_variants(dominant_results.name, outfile, silent)
    
    

if __name__ == '__main__':
    analyze()