from genmod.annotate_variants import (get_frequencies)

def test_get_frequencies(thousand_g_handle):
    chrom = '1'
    start = '879537'
    alt = 'C'
    frequencies = get_frequencies(thousand_g_handle, chrom, start, alt)
    assert frequencies
    assert frequencies['AF'] == '0.000199681'
    assert frequencies['MAX_AF'] == '0.001'

def test_get_frequencies_wrong_alt(thousand_g_handle):
    chrom = '1'
    start = '879537'
    alt = 'T'
    frequencies = get_frequencies(thousand_g_handle, chrom, start, alt)
    assert frequencies == {}
    

def test_get_non_existing_tabix_record(thousand_g_handle):
    chrom = '1'
    start = '10'
    alt = 'C'
    
    frequencies = get_frequencies(thousand_g_handle, chrom, start, alt)
    
    assert frequencies == {}