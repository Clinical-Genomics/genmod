from genmod.utils import get_chromosome_priority


def test_get_chromosome_priority():
    """docstring for test_get_chromosome_priority"""
    assert get_chromosome_priority(chrom='1', chrom_dict={}) == '1'

def test_get_X_priority():
    """docstring for test_get_X_priority"""
    assert get_chromosome_priority(chrom='X', chrom_dict={}) == '23'

def test_get_Y_priority():
    """docstring for test_get_Y_priority"""
    assert get_chromosome_priority(chrom='Y', chrom_dict={}) == '24'

def test_get_MT_priority():
    """docstring for test_get_MT_priority"""
    assert get_chromosome_priority(chrom='MT', chrom_dict={}) == '25'

def test_get_OTHER_priority():
    """docstring for test_get_MT_priority"""
    assert get_chromosome_priority(chrom='GL37', chrom_dict={}) == '26'

def test_get_chr_prority():
    """docstring for test_get_chr_prority"""
    assert get_chromosome_priority(chrom='chr1', chrom_dict={}) == '1'

def test_get_custom_prority():
    """docstring for test_get_chr_prority"""
    assert get_chromosome_priority(chrom='AHA_1', chrom_dict={'AHA_1':2, 'AHA_2':3}) == 2
    
