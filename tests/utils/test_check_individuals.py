from genmod.utils import check_individuals

import pytest

def test_correct():
    """
    Test if return True when the individuals exist
    """
    ped_individuals = ['1','2','3']
    vcf_individuals = ['1','2','3']
    
    assert check_individuals(ped_individuals, vcf_individuals)


def test_wrong():
    """
    Test if raise error when the individuals not exist
    """
    ped_individuals = ['1','2','3', '4']
    vcf_individuals = ['1','2','3']
    
    with pytest.raises(IOError):
        check_individuals(ped_individuals, vcf_individuals)