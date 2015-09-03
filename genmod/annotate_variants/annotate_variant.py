from __future__ import absolute_import, print_function

import logging

from genmod.vcf_tools import add_vcf_info

from . import get_frequency

def annotate_thousand_g(variant_line, thousand_g_handle, keyword="1000GAF"):
    """
    Annotate a variant line with the given annotations
    """
    splitted_line = variant_line.rstrip().split('\t')
    chrom = splitted_line[0]
    pos = splitted_line[1]
    alt = splitted_line[4]
    
    thousand_g_freq = get_frequency(
        tabix_reader=thousand_g, 
        chrom=chrom, 
        start=pos, 
        alt=alt
    )
    
    if thousand_g_freq:
        variant_line = add_vcf_info(
            variant_line=variant_line, 
            keyword=keyword, 
            annotation=str(thousand_g_freq))
    
    return variant_line
    