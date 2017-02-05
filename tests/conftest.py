import os
import pytest
import tabix
from tempfile import NamedTemporaryFile

from genmod.annotations import ensembl_path_37 as ensembl_path

from genmod.vcf_tools.header_parser import HeaderParser

thousandg_path = "tests/fixtures/annotate_variant/small_1000G_maxAF.vcf.gz"
thousandg_chr_path = "tests/fixtures/annotate_variant/small_1000G_chr.vcf.gz"
vcf = "tests/fixtures/test_vcf.vcf"

@pytest.fixture(scope='function')
def thousand_g_path(request):
    """Return the path to a bgzipped 1000G file"""
    return thousandg_path

@pytest.fixture(scope='function')
def thousand_g_chr_path(request):
    """Return the path to a bgzipped 1000G file"""
    return thousandg_chr_path

@pytest.fixture(scope='function')
def thousand_g_handle(request, thousand_g_path):
    """Return a tabix handle with a 1000G file"""
    thousand_g = tabix.open(thousand_g_path)
    return thousand_g

@pytest.fixture(scope='function')
def thousand_g_chr_handle(request, thousand_g_chr_path):
    """Return a tabix handle with a 1000G file"""
    thousand_g = tabix.open(thousand_g_chr_path)
    return thousand_g

@pytest.fixture(scope='function')
def vcf_path(request):
    """Return the path to a vcf file"""
    return vcf

@pytest.fixture(scope='function')
def header(request, vcf_path):
    """Return the a header object"""
    head = HeaderParser()

    with open(vcf_path, 'r') as variant_file:
        for line in variant_file:
            line = line.rstrip()
        
            if line.startswith('#'):
                if line.startswith('##'):
                    head.parse_meta_data(line)
                else:
                    head.parse_header_line(line)
            else:
                break
    
    return head


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