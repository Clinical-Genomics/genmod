import logging

from genmod.annotate_variants import get_frequency

def annotate_frequency(chrom, position, alternative, 
tabix_handle=None):
    """Return the frequency
    
        Args:
            chrom (str): The chromosome
            position (int): Start position for variant
            alternative (str): The variant
            tabix_handle (Tabix.Handle)
        
        Returns:
            frequency (str): The frequency found
    """
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
    