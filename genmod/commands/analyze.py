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

def print_variants(variant_dict, variant_file, header):
    """print the variants to a file"""
    with open(variant_file, mode='w', encoding='utf-8') as f:
        for variant_id in variant_dict:
            print_line = [variant_dict[variant_id].get(entry, '-') for entry in header]
            f.write('\t'.join(print_line)+'\n')
    return

def print_results(variant_dict, mode = 'homozygote', outfile=None, silent=False):
    """Print the variants to a results file or stdout."""
    
    score_key = 'CADD'
    score_dict = {} # A dictionary with {variant_id: score}. Score is usually cadd score or rank score
    # for variant_id, variant in sorted(variant_dict.items(), key = lambda sort_key: float(sort_key[1]['info_dict'].get('CADD', '0')), reverse=True):
    column_width = 12
    length_of_output = 40
    for variant_id in variant_dict:
        score = float(variant_dict[variant_id]['info_dict'].get(score_key, '0'))
        if mode == 'compound':
            for variant_2_id in variant_dict[variant_id]['info_dict'].get('Compounds', '').split(','):
                if variant_2_id in variant_dict:
                    score_2 = float(variant_dict[variant_2_id]['info_dict'].get(score_key, '0'))
                    if score_2 > 10:
                        # print(variant_dict[variant_2_id])
                        variant_pair = (variant_id, variant_2_id)
                        score = (score + score_2)/2
                        already_scored = [set(var_pair) for var_pair in list(score_dict.keys())]
                        if set(variant_pair) not in already_scored:
                            score_dict[variant_pair] = score
        else:
            score_dict[variant_id] = score
            print('score_dict:' % score_dict)
            print('variant_dict:' % variant_dict)
    if mode == 'compound':
        print('\nCompound analysis:\n')
        variant_header = ['Variant 1 in pair', 'Variant 2 in pair']
        header = ['Chrom', 
                    'Position', 
                    'Reference', 
                    'Alternative', 
                    'Cadd score', 
                    '1000GMAF', 
                    'Position', 
                    'Reference', 
                    'Alternative', 
                    'Cadd score', 
                    '1000GMAF', 
                    'Annotation']
        print(''.join(word.ljust(column_width*6) for word in variant_header))
        print(''.join(word.ljust(column_width) for word in header))
        
    else:
        if mode == 'dominant':
            print('\nDominant analysis:\n')            
        if mode == 'homozygote':
            print('\nHomozygote analysis:\n')
        if mode == 'denovo':
            print('\nDe novo analysis:\n')
        if mode == 'xlinked':
            print('\nX-linked analysis:\n')
        header = ['Chrom', 'Position', 'Reference', 'Alternative', 'Cadd score', '1000GMAF', 'Annotation']
    
        print(''.join(word.ljust(column_width) for word in header))
    
    i = 0
    for variant_id in sorted(score_dict, key=score_dict.get, reverse=True):
        if mode == 'compound':
            print_line = [variant_dict[variant_id[0]]['CHROM'],
                            variant_dict[variant_id[0]]['POS'],
                            variant_dict[variant_id[0]]['REF'],
                            variant_dict[variant_id[0]]['ALT'],
                            variant_dict[variant_id[0]]['info_dict'].get('CADD', '-'),
                            variant_dict[variant_id[0]]['info_dict'].get('1000GMAF', '-'),
                            variant_dict[variant_id[1]]['POS'],
                            variant_dict[variant_id[1]]['REF'],
                            variant_dict[variant_id[1]]['ALT'],
                            variant_dict[variant_id[1]]['info_dict'].get('CADD', '-'),
                            variant_dict[variant_id[1]]['info_dict'].get('1000GMAF', '-'),
                            variant_dict[variant_id[0]]['info_dict'].get('Annotation', '-')
                            ]
            
        else:
            print_line = [variant_dict[variant_id]['CHROM'],
                            variant_dict[variant_id]['POS'],
                            variant_dict[variant_id]['REF'],
                            variant_dict[variant_id]['ALT'],
                            variant_dict[variant_id]['info_dict'].get('CADD', '-'),
                            variant_dict[variant_id]['info_dict'].get('1000GMAF', '-'),
                            variant_dict[variant_id]['info_dict'].get('Annotation', '-')
                            ]
        if i < length_of_output:
            print(''.join(word.ljust(column_width) for word in print_line))
        i += 1
            
    
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



def remove_inacurate_compounds(compound_dict):
    """If the second variant in a compound pair does not meet the requirements they should not be considered."""
    
    for variant_id in list(compound_dict.keys()):
        # Get the compounds for the variant
        compounds = compound_dict[variant_id]['info_dict'].get('Compounds', '').split(',')
        compound_set = set(compounds) 
        high_maf_compounds = 0
        for compound in compounds:
            # If requrements are not met it has never been placed in compound dict
            if compound not in compound_dict:
                compound_set.discard(compound)
        # If no compounds in the pair upfills the requirements we remove the pair
        if len(compound_set) == 0:
            compound_dict.pop(variant_id)
    return


def covered_in_all(variant, coverage_treshold = 7):
    """Check if the variant is covered in all individuals."""
    for individual in variant['genotypes']:
        if variant['genotypes'][individual].quality_depth < coverage_treshold:
            return False
    return True
        

def get_interesting_variants(variant_parser, dominant_dict, homozygote_dict, compound_dict, x_linked_dict, 
    dominant_dn_dict, freq_treshold, freq_keyword, cadd_treshold, cadd_keyword):
    """Collect the interesting variants in their dictionarys. add RankScore."""
    
    inheritance_keyword = 'GeneticModels'
    
    gq_treshold = 100
    
    de_novo_set = set(['AD_dn', 'AR_hom_dn', 'AR_comp_dn', 'XD_dn', 'XR_dn'])
    dominant_set = set(['AD'])
    homozygote_set = set(['AR_hom'])
    compound_set = set(['AR_comp'])
    x_linked_set = set(['XD', 'XR'])
    dominant_dn_set = set(['XD_dn'])
    
    for variant in variant_parser:
        models_found = set(variant['info_dict'].get(inheritance_keyword, '').split(','))
        
        maf = min([float(frequency) for frequency in variant['info_dict'].get(freq_keyword, '0').split(',')])
        cadd_score = max([float(cscore) for cscore in variant['info_dict'].get(cadd_keyword, '0').split(',')])
        
        variant_id = variant['variant_id']
        
        if covered_in_all(variant):
        
            if variant['FILTER'] == 'PASS' and float(variant['QUAL']) > gq_treshold:
                # Check if cadd score is available:
                if cadd_score > cadd_treshold:
                    # Check if MAF is below treshold:
                    if maf < freq_treshold:
                        # First we look at the variants that are not dn:
                        if models_found.intersection(dominant_set):
                            dominant_dict[variant_id] = variant
                        if models_found.intersection(homozygote_set):
                            homozygote_dict[variant_id] = variant
                        if models_found.intersection(compound_set):
                            compound_dict[variant_id] = variant
                        if models_found.intersection(x_linked_set):
                            x_linked_dict[variant_id] = variant
                        elif models_found.intersection(dominant_dn_set):
                                dominant_dn_dict[variant_id] = variant
    return


@click.command()
@click.argument('variant_file',
                    nargs=1,
                    type=click.Path(exists=True),
                    metavar='<vcf_file> or "-"'
)
# @click.option('-c', '--config_file',
#                     type=click.Path(exists=True),
#                     help="""Specify the path to a config file."""
# )
@click.option('--frequency_treshold', '-freq',
                    default=0.03, 
                    nargs=1,
                    help='Specify maf treshold for variants to be considered. Default 0.03'
)
@click.option('--frequency_keyword', '-freqkey',
                    default='1000GMAF', 
                    nargs=1,
                    help='Specify keyword for frequency in vcf. Default 1000GMAF'
)
@click.option('--cadd_treshold', '-cadd',
                    default=10.0, 
                    nargs=1,
                    help='Specify the cadd treshold for variants to be considered. Default 10.0'
)
@click.option('--cadd_keyword', '-caddkey',
                    default='CADD', 
                    nargs=1,
                    help='Specify keyword for CADD scores in vcf. Default CADD'
)
# @click.option('-p', '--patterns',
#                     type=click.Choice(['AR', 'AD', 'X']),
#                     multiple=True,
#                     help='Specify the inheritance patterns. Default is all patterns'
# )
# @click.option('-o', '--outfile',
#                     type=click.Path(exists=False),
#                     help='Specify the path to a file where results should be stored.'
# )
@click.option('-s', '--silent', 
                is_flag=True,
                help='Do not output variants.'
)
@click.option('-v', '--verbose', 
                is_flag=True,
                help='Increase output verbosity.'
)
def analyze(variant_file, frequency_treshold, frequency_keyword, cadd_treshold, cadd_keyword, silent, verbose):
    """Analyze the annotated variants in a VCF file."""
    
    start_time_analysis = datetime.now()
    
    # configs = ConfigObj(config_file)        
    # prefered_models = make_models([])
    
    inheritance_keyword = 'GeneticModels'
    
    gq_treshold = 100
    
    # if config_file:
    #     frequency_treshold = float(configs.get('frequency', {}).get('rare', frequency_treshold))
    #     freq_keyword = configs.get('frequency', {}).get('keyword', freq_keyword)
    #     inheritance_patterns = [pattern for pattern in configs.get('inheritance', {}).get('patterns',[])]
    #     inheritance_keyword = configs.get('inheritance', {}).get('keyword',inheritance_keyword)
    #     prefered_models = make_models(inheritance_patterns)
    
    
    if variant_file == '-':
        variant_parser = vcf_parser.VCFParser(fsock = sys.stdin)
    else:
        variant_parser = vcf_parser.VCFParser(infile = variant_file)
    
    head = variant_parser.metadata
    
    
    dominant_dict = {}
    homozygote_dict = {}
    compound_dict = {}
    x_linked_dict = {}
    dominant_dn_dict = {}
    
    
    get_interesting_variants(variant_parser, dominant_dict, homozygote_dict, compound_dict, x_linked_dict, dominant_dn_dict,
                                frequency_treshold, frequency_keyword, cadd_treshold, cadd_keyword)
    
    if len(dominant_dict) > 0:
        print_results(dominant_dict, mode='dominant')
                
    if len(homozygote_dict) > 0:
        
        print_results(homozygote_dict, mode='homozygote')
        
    if len(compound_dict) > 0:
        print_results(compound_dict, mode='compound')
    #     compound_file = NamedTemporaryFile(delete=False)
    #     print_variants(compound_dict, compound_file.name)
    #     compound_file.close()
    #     homozygote_sorted = NamedTemporaryFile(delete=False, suffix='.vcf')
    #     homozygote_sorted.close()
    #     print_headers(head, homozygote_sorted.name)
    #     var_sorter = variant_sorter.FileSort(homozygote_file.name, mode='cadd', outfile=homozygote_sorted.name)
    #     var_sorter.sort()
    if len(x_linked_dict) > 0:
        print_results(x_linked_dict, mode='xlinked')

    if len(dominant_dn_dict) > 0:
        print_results(dominant_dn_dict, mode='denovo')
    #     de_novo_file = NamedTemporaryFile(delete=False)
    #     print_variants(dominant_dn_dict, de_novo_file.name)
    #     de_novo_file.close()
    #     homozygote_sorted = NamedTemporaryFile(delete=False, suffix='.vcf')
    #     homozygote_sorted.close()
    #     print_headers(head, homozygote_sorted.name)
    #     var_sorter = variant_sorter.FileSort(homozygote_file.name, mode='cadd', outfile=homozygote_sorted.name)
    #     var_sorter.sort()
        
    
    
        
    
    # with open(homozygote_sorted.name, mode='r', encoding='utf-8') as f:
    #     for line in f:
    #         print(line.rstrip())
    
    # print_variants(homozygote_sorted.name)
    
    remove_inacurate_compounds(compound_dict)
    
    print('')
    
    print('Dominant variants: %s' % len(dominant_dict))
    print('Homozygote variants: %s' %len(homozygote_dict))
    print('Compound variants: %s' %len(compound_dict))
    print('X-linked variants: %s' %len(x_linked_dict))
    
    # pp(compound_dict)
    
    print('Time for analysis: %s' % str(datetime.now() - start_time_analysis))
    # print_headers(variant_parser.metadata, outfile=outfile)
    
    # dominant_results = NamedTemporaryFile(delete=False)
    # dominant_results.close()
    #
    # var_sorter = variant_sorter.FileSort(dominant_file.name, mode='cadd', outfile=dominant_results.name)
    # var_sorter.sort()
    #
    # print(dominant_results)
    # print(outfile)
    # print_variants(dominant_results.name, outfile, silent)
    
    

if __name__ == '__main__':
    analyze()