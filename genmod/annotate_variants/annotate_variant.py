import logging

from genmod.annotate_variants import get_frequency

def annotate_frequency(chrom, position, alternatives, 
tabix_handle=None):
    """
    Annotate a variant line with the given annotations
    """
    frequencies = []
    for alternative in alternatives.split(','):
        frequency = get_frequency(
            tabix_reader=tabix_handle, 
            chrom=chrom, 
            start=position, 
            alt=alternative
            )
        if frequency:
            frequencies.append(frequency)
        else:
            frequencies.append('0')
    
    return frequencies
    