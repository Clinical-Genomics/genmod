#!/usr/bin/env python
# encoding: utf-8
"""
analyze_variants.py

Analyze the the variants in a vcf, the following will be printed:
    
    - How many variants found
    - How many mendelian violations
    - How many variants where not covered in all individuals. (Default depth 7)
    - How many variants did not satisfy the base call 
        quality treshold. (Default 10)
    - How many variants followed each model:
        - AR_hom
        - AR_comp
        - AR_hom_dn
        - AR_comp_dn
        - AD
        - AD_dn
        - XD
        - XD_dn
        - XR
        - XR_dn
    - How many variants in genetic regions
    - How many rare variants (Default maf < 0.02)
    - How many high scored cadd. (Default cadd = 0)
    - How many rare + high score cadd
    - How many follow a genetic model + rare + high cadd

Created by Måns Magnusson on 2014-09-08.
Copyright (c) 2014 __MoonsoInc__. All rights reserved.
"""

from __future__ import (print_function, division)

import sys
import os
import logging
import json
import click

from codecs import open
from datetime import datetime

import pkg_resources

from vcf_parser import VCFParser
from ped_parser import FamilyParser

from genmod.vcf_tools import HeaderParser


def get_inheritance_models(variant, family_id, inheritance_keyword):
    """Return the genetic models found for this family in this variant"""
    models_found = set([])
    family_models = variant['info_dict'].get(inheritance_keyword, None)
    if family_models:
        #This is a string on the form 'fam_1:AR_hom,fam_2:AR_hom|AR_hom_dn
        for family_info in family_models:
            splitted_family = family_info.split(':')
            if splitted_family[0] == family_id:
                models_found = set(splitted_family[1].split('|'))
    return models_found
    

###           This is for analyzing the variants       ###

@click.command()
@click.argument('variant_file',
                    nargs=1,
                    type=click.Path(exists=True),
                    metavar='<vcf_file> or "-"'
)
@click.option('-f', '--family_file',
                    nargs=1, 
                    type=click.File('r'),
                    metavar='<ped_file>'
)

# @click.option('-c', '--config_file',
#                     type=click.Path(exists=True),
#                     help="""Specify the path to a config file."""
# )
@click.option('--frequency_keyword', '-freqkey',
                    default='1000G_freq', 
                    nargs=1,
                    help='Specify keyword for frequency in vcf. Default 1000G_freq'
)
@click.option('--frequency_treshold', '-freq',
                    default=0.05,
                    nargs=1,
                    help='Specify the ferquency treshold for variants to be considered. Default=0.05'
)
@click.option('--cadd_keyword', '-caddkey',
                    default='CADD', 
                    nargs=1,
                    help='Specify keyword for CADD scores in vcf. Default CADD'
)
@click.option('--gq_treshold', '-gq',
                    default=20.0,
                    nargs=1,
                    help='Specify the genotype quality treshold for variants to be considered. Default=50'
)
@click.option('--read_depth_treshold', '-depth',
                    default=10.0,
                    nargs=1,
                    help="""Specify the genotype quality treshold for variants to be considered. Default=10.
                            The read deth is taken from AD, so it is the sum of the quality reads from reference and alternative alleles."""
)
@click.option('--cadd_treshold', '-cadd',
                    default=12.0, 
                    nargs=1,
                    help='Specify the cadd treshold for variants to be considered. Default 12.0'
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
# @click.option('-v', '--verbose',
#                 is_flag=True,
#                 help='Increase output verbosity.'
# )
def summarize(variant_file, family_file, frequency_treshold, frequency_keyword,
                cadd_treshold, cadd_keyword, gq_treshold, read_depth_treshold):
    """
    Summarize the the variants in a vcf.
    
    There will be one result line per individual.
    
    - How many variants found\n
    - How many variants did not satisfy the base call 
        quality treshold. (Default 20)\n
    - How many variants where not covered in all individuals. 
        (Default depth 10)\n
    - How many variants followed each model in each family:\n
            - AR_hom\n
            - AR_comp\n
            - AR_hom_dn\n
            - AR_comp_dn\n
            - AD\n
            - AD_dn\n
            - XD\n
            - XD_dn\n
            - XR\n
            - XR_dn\n
        - How many rare variants (Default maf < 0.02)\n
        - How many high scored cadd. (Default cadd = 0)\n
        - How many rare + high score cadd\n
        - How many no cadd score\n
        - How many indels\n
        - How many indels without cadd score\n
    
    """
    logger = logging.getLogger(__name__)
    logger = logging.getLogger("genmod.commands.summarize_variants")
    
    head = HeaderParser()
    
    nr_of_variants = 0
    
    header = ['sample_id', 'nr_of_variants']
    
    samples = {}
    
    logger.debug("Setting up a variant parser")
    if variant_file == '-':
        variant_parser = VCFParser(
            fsock = sys.stdin,
            check_info=False
            )
    else:
        variant_parser = VCFParser(
            infile = variant_file,
            check_info=False
            )
    logger.debug("Variant parser setup")
    
    head = variant_parser.metadata
    
    for sample_id in head.individuals:
        samples[sample_id] = {}
        samples[sample_id]["nr_of_variants"] = 0
    
    
    for variant in variant_parser:
        for sample_id in samples:
            samples[sample_id]["nr_of_variants"] += 1
            print(variant['genotypes'][sample_id].depth_of_coverage)
    
    print(json.dumps(samples))
    
    
    # inheritance_models = [
    #                         'AR_hom',
    #                         'AR_hom_dn',
    #                         'AR_comp',
    #                         'AR_comp_dn',
    #                         'AD',
    #                         'AD_dn',
    #                         'XD',
    #                         'XD_dn',
    #                         'XR',
    #                         'XR_dn'
    #                     ]
    #
    # family_dict = {}
    # for family_id in families:
    #     family_dict[family_id] = {}
    #     for inheritance_model in inheritance_models:
    #         family_dict[family_id][inheritance_model] = 0
    #
    #
    # number_of_variants = 0
    # interesting_variants = 0
    # rare_variants = 0
    # high_cadd_scores = 0
    # no_cadd_score = 0
    # high_cadd_and_rare = 0
    # high_gq = 0
    # covered_in_all = 0
    # indels = 0
    # indel_no_cadd = 0
    # true_de_novos = 0
    # low_genotype = 0
    # low_coverage = 0
    # low_genotype_and_low_coverage = 0
    #
    # analysis_start = datetime.now()
    #
    # if variant_file == '-':
    #     variant_parser = VCFParser(fsock = sys.stdin)
    # else:
    #     variant_parser = VCFParser(infile = variant_file)
    #
    # for variant in variant_parser:
    #
    #     maf = min(
    #             [
    #                 float(frequency) for frequency in
    #                     variant['info_dict'].get(
    #                                         frequency_keyword,
    #                                         '0'
    #                                         )
    #             ]
    #             )
    #     cadd_score = max(
    #             [
    #                 float(cscore) for cscore in
    #                 variant['info_dict'].get(
    #                                     cadd_keyword,
    #                                     '0'
    #                                     )
    #             ]
    #             )
    #     reference = variant['REF']
    #     alternative = variant['ALT']
    #
    #     number_of_variants += 1
    #     genotypes = variant.get('genotypes', {})
    #
    #     correct_genotype = True
    #     adequate_depth = True
    #     high_cadd = True
    #     rare = True
    #
    #     for individual in genotypes:
    #         if genotypes[individual].genotype_quality < gq_treshold:
    #             correct_genotype = False
    #
    #         #If any individual has depth below "depth" we do not consider the variant
    #         if genotypes[individual].quality_depth < read_depth_treshold:
    #             adequate_depth = False
    #
    #     if not correct_genotype:
    #         low_genotype += 1
    #         if not adequate_depth:
    #             low_genotype_and_low_coverage += 1
    #     if not adequate_depth:
    #         low_coverage += 1
    #     # We are most interested in the variants that meet the criterias of read depth and genotype quality:
    #     if correct_genotype and adequate_depth:
    #         interesting_variants += 1
    #         # Check the cadd score:
    #
    #         if cadd_score >= cadd_treshold:
    #             high_cadd_scores += 1
    #         else:
    #             high_cadd = False
    #         if cadd_score == 0:
    #             no_cadd_score += 1
    #
    #         # Check the frequency of the variants:
    #         if maf <= frequency_treshold:
    #             rare_variants += 1
    #             if high_cadd:
    #                 high_cadd_and_rare += 1
    #         else:
    #             rare = False
    #
    #         # Check if indel:
    #
    #         if len(reference) > 1 or len(alternative) > 1:
    #             indels += 1
    #             if cadd_score == 0:
    #                 indel_no_cadd += 1
    #
    #     for family_id in families:
    #         models_found = get_inheritance_models(variant, family_id, inheritance_keyword)
    #         if models_found:
    #             for model in models_found:
    #                 family_dict[family_id][model] += 1
    #
    # # pp(inheritance_dict)
    # print("\n\nSUMMARY OF VARIANTS FOUND IN %s.vcf" % vcf_file_name)
    # print("===========================================================================================\n")
    # print('Number of variants in file: %s' % number_of_variants)
    # print('Number of variants with low genotype quality (gq<%s): %s' % (gq_treshold, low_genotype))
    # print('Number of variants with low coverage (cov<%s): %s' % (read_depth_treshold, low_coverage))
    # print('Number of variants with low coverage AND low genotype quality: %s \n\n' % low_genotype_and_low_coverage)
    #
    # print("The following statistics are for the variants that meet the criterias for genotype quality and read depth.\n"
    #         "This means that the variants are covered in all individuals.\n"
    #         "-----------------------------------------------------------------------------------------\n")
    #
    # print("Number of variants to be considered in the analysis(according to the statement above): %s \n" %
    #          (interesting_variants))
    #
    # for family_id in families:
    #     print('Models followed for family %s \n' % family_id)
    #     for model in inheritance_models:
    #         print("%s = %s" % (model, family_dict[family_id][model]))
    #
    #
    # print('\nNumber of rare (maf<%s): %s. Frequency of all: %.2f'
    #         % (frequency_treshold, rare_variants, rare_variants/interesting_variants))
    # print('Number of high cadd scores (cadd >= %s): %s. Frequency of all: %.2f'
    #         % (cadd_treshold, high_cadd_scores, high_cadd_scores/interesting_variants))
    # print('Number of high cadd scores and rare: %s. Frequency of all: %.2f'
    #         % (high_cadd_and_rare, high_cadd_and_rare/interesting_variants))
    # print('Number of no cadd scores: %s. Frequency of all: %.2f \n'
    #         % (no_cadd_score, no_cadd_score/interesting_variants))
    # print('Number of indels: %s. Frequency of all: %.2f'
    #         % (indels, indels/number_of_variants))
    # print('Number of indels and no cadd score: %s. Frequency of all: %.2f \n'
    #             % (indel_no_cadd, indel_no_cadd/number_of_variants))
    # print('Time for analysis: %s' % str(datetime.now()-analysis_start))

if __name__ == '__main__':
    summarize()