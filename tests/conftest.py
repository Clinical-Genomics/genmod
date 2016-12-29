import os
import pytest
import tabix
from tempfile import NamedTemporaryFile

from genmod.annotations import ensembl_path

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


@pytest.fixture(scope='function')
def bed_lines(request):
    """Return a iterable with bed lines"""
    lines = [
        "#Chromosome/scaffold name\tStart (bp)\tGene End (bp)\tHGNC symbol\tGene ID",
        "13\t23708313\t23708703\tHMGA1P6\tENSG00000233440\n",
        "13\t23726725\t23726825\tRNY3P4\tENSG00000207157\n",
        "13\t23743974\t23744736\tLINC00362\tENSG00000229483\n"
    ]
    return lines

@pytest.fixture(scope='function')
def ensembl_file(request):
    """Return the path to ensembl file with region defenitions"""
    return ensembl_path