from genmod.variant_annotation import get_haploblocks

def get_variant_batch():
    """Return a small variant batch"""
    variant_1 = {
        'CHROM': '1',
        'POS': '12',
        'ALT': 'A',
        'REF': 'T',
        '1': '0|1'
    }
    variant_2 = {
        'CHROM': '1',
        'POS': '13',
        'ALT': 'T',
        'REF': 'G',
        '1': '0|1'
    }
    variant_3 = {
        'CHROM': '1',
        'POS': '14',
        'ALT': 'A',
        'ALT': 'C',
        '1': '0|1'
    }
    variant_batch = {
        '1_12_T_A': variant_1,
        '1_13_G_T': variant_2,
        '1_14_C_A': variant_3,
    }
    return variant_batch

def test_simple():
    """Test if get_haploblocks behave as suspected"""
    variant_batch = get_variant_batch()
    
    haploblocks = get_haploblocks(variant_batch, ['1'])
    
    assert list(haploblocks.keys()) == ['1']
    assert haploblocks['1'].find_range([13,13]) == ['1']
