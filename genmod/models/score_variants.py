#!/usr/bin/env python
# encoding: utf-8
"""
score_variant.py

Script that takes a variant as input and modify it with a score depending on its different values.

Possible names for the list of genetic models are:

AD, AD_denovo, AR, AR_denovo, AR_compound, X, X_denovo


Created by Måns Magnusson on 2013-08-14.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from genmod.utils import is_number


def score_variant(variant_object, genetic_models = ['AR_hom', 'AD']):
    """Score a variant object according to Henriks score model. Input: A variant object and a list of genetic models."""
    
    score = 0
    # Models of inheritance
    variant_models = variant_object.models
    
    # Predictors
    mutation_taster = variant_object.all_info.get('Mutation_taster', None)
    avsift = variant_object.all_info.get('SIFT', None)
    poly_phen = variant_object.all_info.get('Poly_phen', None)
    
    # Annotations:
    functional_annotation = variant_object.all_info.get('Functional_annotation', None)
    gene_annotation = variant_object.all_info.get('Gene_annotation', None)
    
    # Frequency in databases:
    thousand_genomes_frequency = variant_object.all_info.get('1000G', None)
    dbsnp_frequency = variant_object.all_info.get('Dbsnp129', None)
    dbsnp_id = variant_object.all_info.get('Dbsnp_nonflagged' ,None)
    
    # Filter
    
    filt = variant_object.all_info.get('GT_call_filter', None)
    
    # Conservation scores:
        # Base
    gerp_base = variant_object.all_info.get('GERP', None)
        # Region
    mce64way = variant_object.all_info.get('Phast_cons_lements', None)
    gerp_region = variant_object.all_info.get('GERP_elements', None)
    
    
    phylop = variant_object.all_info.get('Phylo_p', None)
    
    segdup = variant_object.all_info.get('Genomic_super_dups', None)
    
    hgmd = variant_object.all_info.get('HGMD', None)
    
    score += check_inheritance(variant_models, genetic_models)
    score += check_predictions(mutation_taster, avsift, poly_phen)
    score += check_functional_annotation(functional_annotation)
    score += check_gene_annotation(gene_annotation)
    score += check_frequency_score(thousand_genomes_frequency, dbsnp_frequency, dbsnp_id)
    score += check_filter(filt)
    score += check_region_conservation(mce64way, gerp_region)
    score += check_base_conservation(gerp_base)
    score += check_phylop_score(phylop)
    score += check_segmental_duplication(segdup)
    score += check_hgmd(hgmd)
    variant_object.rank_score = score
    
    # return variant_object
    
def check_inheritance(variant_models, genetic_models):
    """Check if the models of inheritance are followed for the variant."""
    model_score = 0
    #If any of the prefered models are followed:
    for model_followed in variant_models:
        if model_followed in genetic_models:
            model_score = 3
    #Else if any model is followed
    if model_score != 3:
        if variant_models != ['Na']:
            model_score = 1
        else:
            model_score = -12
    return model_score
    
def check_predictions(mutation_taster = None, avsift = None, poly_phen = None):
    """Score the variant based on the scores from prediction databases."""
    prediction_score = 0
    if is_number.is_number(avsift):
        if float(avsift) <= 0.05:
            prediction_score += 1
    if is_number.is_number(mutation_taster):
        if float(mutation_taster) >= 0.05:
            prediction_score += 1
    if is_number.is_number(poly_phen):
        if float(poly_phen) >= 0.85:
            prediction_score += 1
    return prediction_score
    
def check_functional_annotation(functional_annotation = None):
    """Score the variant based on its functional annotation"""
    functional_annotation_score = 0
    if functional_annotation in ['frameshift deletion', 'frameshift insertion', 'nonframeshift deletion', 'nonframeshift insertion', 'stopgain SNV', 'stoploss SNV', 'nonsynonymous SNV']:
        functional_annotation_score += 3
    elif functional_annotation in ['unknown', 'synonymous SNV']:
        functional_annotation_score += 1
    return functional_annotation_score
    
def check_gene_annotation(gene_annotation = None):
    """Score the variant based onits gene annotation."""
    gene_annotation_score = 0
    if gene_annotation in ['exonic', 'exonic;splicing',  'splicing']:
        gene_annotation_score += 3
    elif gene_annotation in ['intronic', 'UTR3', 'UTR5', 'UTR5;UTR3', 'upstream', 'downstream', 'upstream;downstream']:
        gene_annotation_score += 1
    return gene_annotation_score
    
def check_frequency_score(thousand_genomes_frequency = None, dbsnp_frequency = None, dbsnp_id = None):
    """Score the variant based on the frequency in population."""

    frequency_score = 0
    freq_scores = []
    def get_freq_score(frequency):
        """Returns a score depending on the frequency"""
        if is_number.is_number(frequency):
            if float(frequency) <= 0.005:
                return 2
            elif float(frequency) <= 0.02:
                return 1
            #If common variant:
            else:
                    return -12
        else:# If not existing in database
            return 3
    
    freq_scores.append(get_freq_score(thousand_genomes_frequency))
    freq_scores.append(get_freq_score(dbsnp_frequency))
    common = False
    # If the variant if common in any database(if a score is negative) we give a low score:
    for freq_score in freq_scores:
        if freq_score < 0:
            common = True
    if common:
        frequency_score = -12
    else:
        frequency_score += sum(freq_scores) / 2
    # If variant has no ID in dbSNP it get an extra score
        if dbsnp_id == '-':
            frequency_score += 1
    return frequency_score
    
def check_filter(filt):
    """Check if variant has passed the filter process."""
    filter_score = 0
    if filt == 'PASS':
        filter_score = 3
    elif filt == 'PRES':
        filter_score = 1
    return filter_score
    
def check_region_conservation(mce64way = None, gerp_region = None):
    """Score the variant based on what annotations it has for the region conservations"""
    region_conservation_score = 0
    if mce64way != '-' and gerp_region != '-':
        region_conservation_score += 2
    elif mce64way != '-' or gerp_region != '-':
        region_conservation_score += 1
    return region_conservation_score
    
def check_base_conservation(gerp_base_score = None):
    """Score the variant based on the base level conservation."""
    base_conservation_score = 0
    if is_number.is_number(gerp_base_score):
        if float(gerp_base_score) >= 4:
            base_conservation_score += 2
        elif float(gerp_base_score) >= 2:
            base_conservation_score += 1
    return base_conservation_score
    
def check_phylop_score(phylop = None):
    """Score the variant based on the Phylop score."""
    phylop_score = 0
    if is_number.is_number(phylop):
        if float(phylop) >= 0.9984188612:
            phylop_score += 2
        elif float(phylop) >= 0.95:
            phylop_score += 1
    return phylop_score
    
def check_segmental_duplication(segdup):
    """Check if there are any annotations for segmental duplication"""
    segdup_score = 0
    if segdup != '-':
        segdup_score -= 2
    return segdup_score
    
def check_hgmd(hgmd):
    """Check if the variant have any annotation from hgmd"""
    hgmd_score = 0
    if hgmd != '-':
        hgmd_score += 1
    return hgmd_score
    
    



def main():
    pass


if __name__ == '__main__':
    main()

