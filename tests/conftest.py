import pytest
import tabix

@pytest.fixture(scope='function')
def thousand_g_path(request):
    """Return the path to a bgzipped 1000G file"""
    thousand_g = "tests/fixtures/annotate_variant/small_1000G_maxAF.vcf.gz"
    return thousand_g

@pytest.fixture(scope='function')
def thousand_g_handle(request, thousand_g_path):
    """Return a tabix handle with a 1000G file"""
    thousand_g = tabix.open(thousand_g_path)
    return thousand_g
