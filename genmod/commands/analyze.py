#!/usr/bin/env python
# encoding: utf-8
"""
analyze.py

Command line tool for analyzing variants that are annotated with genmod.

Created by Måns Magnusson on 2014-09-03.
Copyright (c) 2014 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function

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

from vcf_parser import VCFParser

import genmod

from genmod.utils import print_headers
from genmod.errors import warning

# This is an ad hoc solution to remove huge mostly uninteresting genes. 
# Please modify this set for your own needs
PROBLEMATIC_GENES = set(['MIR6077-1',
                            'MIR6077-2',
                            'MIR4315-1',
                            'MIR4315-2',
                            'LINC00623',
                            'LINC00869',
                            'NBPF8',
                            'NBPF9',
                            'NBPF20',
                            'PPIAL4A',
                            'PPIAL4B',
                            'PPIAL4C',
                            'PDE4DIP',
                            'LOC100132057',
                            'LOC100288162',
                            'SRGAP2D',
                            'FAM272C',
                            'SNAR-A3',
                            'SNAR-A4',
                            'SNAR-A5',
                            'SNAR-A6',
                            'SNAR-A7',
                            'SNAR-A8',
                            'SNAR-A9',
                            'SNAR-A10',
                            'SNAR-A11',
                            'SNAR-A14',
                            'GLUD1P7',
                        ])

def check_families(variant_file):
    """Loop through the vcf file and check which families that are found."""
    families = set([])
    if variant_file == '-':
        variant_parser = VCFParser(fsock = sys.stdin)
    else:
        variant_parser = VCFParser(infile = variant_file)
    for variant in variant_parser:
        genetic_models = variant['info_dict'].get('GeneticModels', None)
        if genetic_models:
            for family_models in genetic_models:
                family = family_models.split(':')[0]
                families.add(family)
    return families


def print_results(variant_dict, 
                    outfile, 
                    vcf_header,
                    family_id, 
                    score_key='CADD', 
                    freq_key='1000G_freq', 
                    mode = 'homozygote', 
                    silent=False):
    """Print the variants to a results file or stdout."""
    
    score_dict = {} # A dictionary with {variant_id: score}. Score is usually cadd score or rank score
    # for variant_id, variant in sorted(variant_dict.items(), key = lambda sort_key: float(sort_key[1]['info_dict'].get('CADD', '0')), reverse=True):
    column_width = 12
    length_of_output = 20
    for variant_id in variant_dict:
        # Get the score for each variant:
        max_score = max(
                            [
                                float(score) for score in 
                                variant_dict[variant_id]['info_dict'].get(
                                                                    score_key, 
                                                                    '0')
                            ]
                        )
        if mode == 'compound':
            # If we look at compounds we want to consider the combined score
            family_compounds = compound_dict[variant_id]['info_dict'].get('Compounds', None)
            if compounds:
                for family in family_compounds:
                    splitted_compounds = family.split(':')
                    if splitted_compounds[0] == family_id:
                        compounds = splitted_compounds[1].split('|')
            
                for variant_2_id in compounds:
                    if variant_2_id in variant_dict:
                        max_score_2 = max(
                                        [
                                            float(score) for score in
                                            variant_dict[variant_2_id]['info_dict'].get(
                                                                                score_key, 
                                                                                '0')
                                        ]
                                    )
                    if max_score_2 > 10:
                        # print(variant_dict[variant_2_id])
                        variant_pair = (variant_id, variant_2_id)
                        max_score = (max_score + max_score_2)/2
                        already_scored = [set(var_pair) for var_pair in list(score_dict.keys())]
                        if set(variant_pair) not in already_scored:
                            score_dict[variant_pair] = max_score
        else:
            score_dict[variant_id] = max_score
    
    if mode == 'compound':
        print('\nCompound analysis:\n')
    if mode == 'dominant':
        print('\nDominant analysis:\n')
    if mode == 'homozygote':
        print('\nHomozygote analysis:\n')
    if mode == 'denovo':
        print('\nDe novo analysis:\n')
    if mode == 'xlinked':
        print('\nX-linked analysis:\n')
    header = ['Chrom', 
                'Position', 
                'Reference', 
                'Alternative', 
                'Cadd score', 
                '1000GMAF', 
                'Annotation'
            ]
    
    print(''.join(word.ljust(column_width) for word in header))
    
    i = 0
    
    with open(outfile , mode='a', encoding='utf-8') as f:
        for variant_id in sorted(score_dict, key=score_dict.get, reverse=True):
            if mode == 'compound':
                if i < length_of_output:
                    print('Pair %s' % (i+1))
                for compound_id in variant_id:
                    print_line = [variant_dict[compound_id]['CHROM'],
                                    variant_dict[compound_id]['POS'],
                                    variant_dict[compound_id]['REF'],
                                    variant_dict[compound_id]['ALT'],
                                    variant_dict[compound_id]['info_dict'].get(score_key, '-'),
                                    variant_dict[compound_id]['info_dict'].get(freq_key, '-'),
                                    variant_dict[compound_id]['info_dict'].get('Annotation', '-')
                                ]
                    if i < length_of_output:
                        print(''.join(word.ljust(column_width) for word in print_line))
                    print_line = [variant_dict[compound_id].get(entry, '-') for entry in vcf_header]
                    f.write('\t'.join(print_line)+'\n')

            else:
                print_line = [variant_dict[variant_id]['CHROM'],
                                variant_dict[variant_id]['POS'],
                                variant_dict[variant_id]['REF'],
                                variant_dict[variant_id]['ALT'],
                                variant_dict[variant_id]['info_dict'].get(score_key, ['-'])[0],
                                variant_dict[variant_id]['info_dict'].get(freq_key, ['-'])[0],
                                variant_dict[variant_id]['info_dict'].get('Annotation', ['-'])[0]
                            ]
                # Print the highest ranked variants to screen:
                if i < length_of_output:
                    print(''.join(word.ljust(column_width) for word in print_line))
                print_line = [variant_dict[variant_id].get(entry, '-') for entry in vcf_header]
                f.write('\t'.join(print_line)+'\n')
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



def remove_inacurate_compounds(compound_dict, family_id):
    """If the second variant in a compound pair does not meet the requirements they should not be considered."""
    
    for variant_id in list(compound_dict.keys()):
        # Get the compounds for the variant
        family_compounds = compound_dict[variant_id]['info_dict'].get('Compounds', None)
        if compounds:
            for family in family_compounds.split(','):
                splitted_compounds = family.split(':')
                if splitted_compounds[0] == family_id:
                    compounds = splitted_compounds[1].split('|')
                    compound_set = set(compounds) 
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
        

def get_interesting_variants(variant_parser, family_id, dominant_dict, 
        homozygote_dict, compound_dict, x_linked_dict, dominant_dn_dict, 
        freq_treshold, freq_keyword, cadd_treshold, cadd_keyword, gq_treshold, 
        coverage, exclude_problematic):
    """Collect the interesting variants in their dictionarys. add RankScore."""
    
    inheritance_keyword = 'GeneticModels'
    
    
    de_novo_set = set(['AD_dn', 'AR_hom_dn', 'AR_comp_dn', 'XD_dn', 'XR_dn'])
    dominant_set = set(['AD'])
    homozygote_set = set(['AR_hom'])
    compound_set = set(['AR_comp'])
    x_linked_set = set(['XD', 'XR'])
    dominant_dn_set = set(['AD_dn'])
    
    
    for variant in variant_parser:
        annotation = set(variant['info_dict'].get('Annotation', ''))
        models_found = set([])
        
        family_models = variant['info_dict'].get(inheritance_keyword, None)
        if family_models:
            #This is a string on the form 'fam_1:AR_hom,fam_2:AR_hom|AR_hom_dn
            for family_info in family_models:
                splitted_family = family_info.split(':')
                if splitted_family[0] == family_id:
                    models_found = set(splitted_family[1].split('|'))
        
        maf = min(
                    [
                        float(frequency) for frequency in 
                        variant['info_dict'].get(
                                            freq_keyword, 
                                            '0'
                                            )
                    ]
                )
        cadd_score = max(
                            [
                                float(cscore) for cscore in 
                                variant['info_dict'].get(
                                                cadd_keyword, 
                                                '0'
                                                )
                            ]
                        )
        
        variant_id = variant['variant_id']
        
        # There is a list of huge genes that becomes problematic when analysing single individuals
        
        interesting = True
        
        if not models_found:
            interesting = False
        
        if exclude_problematic:
            if annotation.intersection(PROBLEMATIC_GENES):
                interesting = False
            
        # if not covered_in_all(variant, coverage):
        #         interesting = False
            
        if not variant['FILTER'] == 'PASS':
            interesting = False
        
        if not float(variant['QUAL']) > gq_treshold:
            interesting = False
        
        
        if interesting:
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
                    if models_found.intersection(dominant_dn_set):
                        dominant_dn_dict[variant_id] = variant
                            
    return


@click.command()
@click.argument('variant_file',
                    nargs=1,
                    type=click.Path(exists=True),
                    metavar='<vcf_file> or "-"'
)
@click.option('-t' ,'--family_type', 
                type=click.Choice(['ped', 'alt', 'cmms', 'mip']), 
                default='ped',
                help="""If the analysis use one of the known setups, 
                      please specify which one."""
)
# @click.option('-c', '--config_file',
#                     type=click.Path(exists=True),
#                     help="""Specify the path to a config file."""
# )
@click.option('--frequency_treshold', '-freq',
                    default=0.02, 
                    nargs=1,
                    help="""Specify maf treshold for variants to be considered.
                            Default 0.02"""
)
@click.option('--frequency_keyword', '-freqkey',
                    default='1000G_freq', 
                    nargs=1,
                    help="""Specify keyword for frequency in vcf. 
                            Default 1000G_freq"""
)
@click.option('--cadd_treshold', '-cadd',
                    default=12.0, 
                    nargs=1,
                    help="""Specify the cadd treshold for variants to be 
                            considered. Default 12.0"""
)
@click.option('--cadd_keyword', '-caddkey',
                    default='CADD', 
                    nargs=1,
                    help="""Specify keyword for CADD scores in vcf. 
                            Default CADD"""
)
@click.option('--coverage', '-cov',
                    default=7, 
                    nargs=1,
                    help="""Specify minimum read depth in all individuals for
                           variant to be considered. Default 7"""
)
@click.option('--gq_treshold', '-gq',
                    default=20, 
                    nargs=1,
                    help="""Specify genotype quality treshold for variants 
                            to be considered. Default 20."""
)
# @click.option('-p', '--patterns',
#                     type=click.Choice(['AR', 'AD', 'X']),
#                     multiple=True,
#                     help='Specify the inheritance patterns. Default is all patterns'
# )
@click.option('-o', '--outdir',
                    type=click.Path(exists=True),
                    default=os.getcwd(),
                    help="""Specify the path to a directory where results 
                            should be stored. Default is ./"""
)
@click.option('-s', '--silent', 
                is_flag=True,
                help='Do not output variants.'
)
@click.option('-exclude', '--exclude_problematic', 
                is_flag=True,
                help="""Exclude problematic genes. This flag is preferable 
                        if analysis of only one individual."""
)
@click.option('-v', '--verbose', 
                is_flag=True,
                help='Increase output verbosity.'
)
def analyze(variant_file, family_type, frequency_treshold, frequency_keyword, 
            cadd_treshold, cadd_keyword, coverage, gq_treshold, outdir, silent,
            exclude_problematic, verbose):
    """Analyze the annotated variants in a VCF file. 
        
        If there are multiple families in the ped one analysis per family will
        be done. The variants are analyzed in five different categories based 
        on what inheritance patterns that are followed.
        The differen analysies are: 
        
                AR compound\n
                AR homozygote\n
                Dominant\n
                X linked\n
                Dominant dn\n
        
        Which variants to be considered are specified in the command line. 
        Defaults are (based on a rare disease assumption):
        
            MAF < 0.02\n
            CADD score > 12\n
            Coverage in all individuals > 7\n
            Call quality > 20\n
        
        The highest scoring variants of each category is printed to screen.
        The full list of each category is printed to new vcf files in a 
        directory specified by the user. Default current dir.
        File names are the same like the input vcf with the name of the 
        analysis appended.
    
    """
    
    start_time_analysis = datetime.now()
    
    # configs = ConfigObj(config_file)        
    # prefered_models = make_models([])
    
    inheritance_keyword = 'GeneticModels'
    families = check_families(variant_file)
    file_name = os.path.splitext(os.path.split(variant_file)[-1])[0]
    
    
    # if config_file:
    #     frequency_treshold = float(configs.get('frequency', {}).get('rare', frequency_treshold))
    #     freq_keyword = configs.get('frequency', {}).get('keyword', freq_keyword)
    #     inheritance_patterns = [pattern for pattern in configs.get('inheritance', {}).get('patterns',[])]
    #     inheritance_keyword = configs.get('inheritance', {}).get('keyword',inheritance_keyword)
    #     prefered_models = make_models(inheritance_patterns)
    
    if variant_file == '-':
        variant_parser = VCFParser(fsock = sys.stdin)
    else:
        variant_parser = VCFParser(infile = variant_file)
    
    for family_id in families:
        print('Analysis for family: %s' % family_id)
    
        head = variant_parser.metadata
            
        dominant_dict = {}
        homozygote_dict = {}
        compound_dict = {}
        x_linked_dict = {}
        dominant_dn_dict = {}
        
        
        get_interesting_variants(variant_parser,
                                    family_id,
                                    dominant_dict, 
                                    homozygote_dict, 
                                    compound_dict, 
                                    x_linked_dict, 
                                    dominant_dn_dict, 
                                    frequency_treshold, 
                                    frequency_keyword, 
                                    cadd_treshold, 
                                    cadd_keyword,
                                    gq_treshold,
                                    coverage, 
                                    exclude_problematic)
        
        remove_inacurate_compounds(compound_dict, family_id)
        
        if len(dominant_dict) > 0:
            dominant_file = os.path.join(
                                    outdir, 
                                    file_name+'_dominant_analysis.vcf'
                                    )
            
            print_headers(head, dominant_file)
            
            print_results(
                      dominant_dict, 
                      dominant_file, 
                      family_id, 
                      variant_parser.header, 
                      cadd_keyword, 
                      frequency_keyword, 
                      mode='dominant'
                      )
                    
        if len(homozygote_dict) > 0:
            homozygote_file = os.path.join(
                                      outdir, 
                                      file_name+'_homozygote_analysis.vcf'
                                      )
            print_headers(head, homozygote_file)
            
            print_results(
                      homozygote_dict, 
                      homozygote_file, 
                      family_id, 
                      variant_parser.header, 
                      cadd_keyword, 
                      frequency_keyword, 
                      mode='homozygote'
                      )
            
        if len(compound_dict) > 0:
            compound_file = os.path.join(
                      outdir, 
                      file_name+'_compound_analysis.vcf'
                      )
            print_headers(head, compound_file)
            
            print_results(
                      compound_dict, 
                      compound_file, 
                      family_id, 
                      variant_parser.header, 
                      cadd_keyword, 
                      frequency_keyword, 
                      mode='compound'
                      )
        
        if len(x_linked_dict) > 0:
            xlinked_file = os.path.join(
                                  outdir, 
                                  file_name+'_x_linked_analysis.vcf'
                                  )
            print_headers(head, xlinked_file)
            
            print_results(
                      x_linked_dict, 
                      xlinked_file, 
                      family_id, 
                      variant_parser.header, 
                      cadd_keyword, 
                      frequency_keyword, 
                      mode='xlinked'
                      )
        
        if len(dominant_dn_dict) > 0:
            dominant_dn_file = os.path.join(
                      outdir, 
                      file_name+'_ad_denovo_analysis.vcf'
                      )
            print_headers(head, dominant_dn_file)
        
            print_results(
                      dominant_dn_dict, 
                      dominant_dn_file, 
                      family_id, 
                      variant_parser.header, 
                      cadd_keyword, 
                      frequency_keyword, 
                      mode='denovo'
                      )
        
        print('')
        
        print('Number of interesting Dominant variants: %s' % 
                len(dominant_dict))
        print('Number of interesting Homozygote variants: %s' %
                len(homozygote_dict))
        print('Number of interesting Compound variants: %s' %
                len(compound_dict))
        print('Number of interesting X-linked variants: %s' %
                len(x_linked_dict))
        print('Number of interesting Autosomal Dominant de novo variants: %s' %
                len(dominant_dn_dict))
        
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