from genmod.annotate_models import get_haploblocks

def get_variant(**kwargs):
    """
    Construct a variant and return it
    """
    variant = {
        'CHROM': kwargs.get('CHROM', '1'),
        'POS': kwargs.get('POS', '12'),
        'REF': kwargs.get('REF', 'A'),
        'ALT': kwargs.get('REF', 'C'),
        'FILTER': 'PASS'
    }
    for indivivdual_id in kwargs.get('genotype_calls',{'1':'0|1'}):
        variant[indivivdual_id] = kwargs.get(
            'genotype_calls',{'1':'0|1'})[indivivdual_id]
    return variant

def get_variant_batch():
    """Return a small variant batch"""
    from collections import OrderedDict
    variant_batch = OrderedDict()
    variant_1 = get_variant(**{
        'POS': '12',
        'genotype_calls':{
            '1': '0|1',
            '2': '0/1',
        }
    })
    variant_2 = get_variant(**{
        'POS': '13',
        'genotype_calls':{
            '1': '0|1',
            '2': '0/1',
        }
    })
    variant_3 = get_variant(**{
        'POS': '14',
        'genotype_calls':{
            '1': '0|1',
            '2': '0|1',
        }
    })
    variant_4 = get_variant(**{
        'POS': '15',
        'genotype_calls':{
            '1': '0/1',
            '2': '0|1',
        }
    })
    variant_5 = get_variant(**{
        'POS': '16',
        'genotype_calls':{
            '1': '0|1',
            '2': '0|1',
        }
    })
    variant_6 = get_variant(**{
        'POS': '17',
        'genotype_calls':{
            '1': '0|1',
            '2': '0/1',
        }
    })
    variant_batch['1_12_A_C'] = variant_1
    variant_batch['1_13_A_C'] = variant_2
    variant_batch['1_14_A_C'] = variant_3
    variant_batch['1_15_A_C'] = variant_4
    variant_batch['1_16_A_C'] = variant_5
    variant_batch['1_17_A_C'] = variant_6
    
    return variant_batch

def test_simple():
    """Test if get_haploblocks behave as suspected"""
    variant_batch = get_variant_batch()
    
    haploblocks = get_haploblocks(variant_batch, ['1', '2'])
    
    assert set(haploblocks.keys()) == set(['1','2'])
    
    assert haploblocks['1'].find_range([12,12]) == ['1']
    assert haploblocks['2'].find_range([12,12]) == []
    
    assert haploblocks['1'].find_range([15,15]) == []
    assert haploblocks['2'].find_range([15,15]) == ['2']

    assert haploblocks['1'].find_range([17,17]) == ['3']
    assert haploblocks['2'].find_range([17,17]) == []
    

