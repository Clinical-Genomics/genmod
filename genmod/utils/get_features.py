from __future__ import (print_function, unicode_literals, absolute_import)

import logging

from . import INTERESTING_SO_TERMS

def check_vep_annotation(variant):
    """
    Return a set with the genes that vep has annotated this variant with.
    
    Arguments:
        variant (dict): A variant dictionary
    
    Returns:
        annotation (set): A set with genes
    """
    
    annotation = set()
    # vep_info is a dictionary with genes as key and annotation as values
    
    for allele in variant.get('vep_info',{}):
        if allele != 'gene_ids':
            for vep_annotation in variant['vep_info'][allele]:
                for consequence in vep_annotation.get('Consequence', {}).split('&'):
                    # These are the SO terms that indicate that the variant 
                    # belongs to a gene
                    if consequence in INTERESTING_SO_TERMS:
                        annotation.add(vep_annotation.get('SYMBOL', ''))
    return annotation

def get_annotation(variant, gene_trees):
    """
    Return the features that a variant belongs to.
    
    Arguments: 
        variant (dict): A variant dictionary
        gene_trees (dict): A dictionary with interval trees
        
    Returns: 
        annotations (set): A set with annotated features
    """
    logger = logging.getLogger(__name__)
    
    annotation = set()
    variant_id = variant.get('variant_id', '')
    # If the variant has already been annotated by genmod we do not need to 
    # check again
    if variant.get('info_dict', {}).get('Annotation', None):
        
        logger.debug("Variant already annotated.")
        return set(variant['info_dict']['Annotation'])
    #If annotated with vep we do not need to check interval trees
    if variant.get("CSQ", None):
        
        logger.debug("Using vep annotation.")        
        return check_vep_annotation(variant)

    chrom = variant['CHROM']
    # Internally we never use 'chr' in the chromosome names:
    if chrom.startswith('chr'):
        chrom = chrom[3:]
    position = int(variant['POS'])
    # When checking what features that are overlapped we use the longest 
    # alternative
    longest_alt = max([len(alternative) for 
                        alternative in variant['ALT'].split(',')])

    variant_interval = [position, (position + longest_alt-1)]
    logger.debug("Checking annotations for variant on chrom:{0} and"\
                " position: {1}".format(chrom, position))


    try:
        annotation = set(gene_trees[chrom].find_range(variant_interval))
        logger.debug("Annotations found for variant {0}: {1}".format(
            variant_id, annotation))
        
    except KeyError:
        logger.warning("Chromosome {0} is not in annotation file".format(chrom))
    
    return annotation

