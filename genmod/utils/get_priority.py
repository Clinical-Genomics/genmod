import logging
from genmod.score_variants import RANK_SCORE_TYPE_NAMES

def get_chromosome_priority(chrom, chrom_dict={}):
    """
    Return the chromosome priority
    
    Arguments:
        chrom (str): The cromosome name from the vcf
        chrom_dict (dict): A map of chromosome names and theis priority
    
    Return:
        priority (str): The priority for this chromosom
    """
    priority = '0'
    
    chrom = chrom.lstrip('chr')
    
    if chrom_dict:
        priority = chrom_dict.get(chrom, '0')
    
    else:
        try:
            if int(chrom) < 23:
                priority = chrom
        except ValueError:
            if chrom == 'X':
                priority = '23'
            elif chrom == 'Y':
                priority = '24'
            elif chrom == 'MT':
                priority = '25'
            else:
                priority = '26'
    
    return priority

def get_rank_score(variant_line=None, variant_dict=None, family_id=0, rank_score_type: str = 'RankScore'):
    """
    Return the rank score priority for a certain family.
    
    If no family is given the first family found is used
    
    Arguments:
        variant_line (str): A vcf variant line
        variant_dict (dict): A variant dictionary
        family_id (str): A family id
        rank_score_type(str): Return rank score based on raw or normalized format
        See the genmod.score_variants.rank_score_variant_definitions for more info.
    
    Return:
        rank_score (str): The rank score for this variant
    """
    if rank_score_type not in RANK_SCORE_TYPE_NAMES:
        raise ValueError('Unknown rank_score_type', rank_score_type)
    
    rank_score = -100
    raw_entry = None
    
    if variant_line:
        variant_line = variant_line.split("\t")
    
        for info_annotation in variant_line[7].split(';'):
            info_annotation = info_annotation.split('=')
            key = None
            if len(info_annotation) == 2:
                key = info_annotation[0]
                value = info_annotation[1]
            if key == rank_score_type:
                raw_entry = value
                break
    
    elif variant_dict:
        raw_entry: str = variant_dict['info_dict'].get(rank_score_type)
    
    if raw_entry:
        for family_annotation in raw_entry.split(','):
            family_annotation = family_annotation.split(':')
            if family_id:
                # If we should sort on a certain family we look for the
                # correct id
                if family_id == family_annotation[0]:
                    rank_score = float(family_annotation[1])
            else:
            # If no family id is given we choose the first family found
                rank_score = float(family_annotation[1])

    return str(rank_score)
    