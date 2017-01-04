from genmod.commands.base import cli
from click.testing import CliRunner

VCF_FILE = "tests/fixtures/test_vcf_regions.vcf"
THOUSAND_G_FILE = "tests/fixtures/annotate_variant/small_1000G.vcf.gz"
CADD_FILE = "tests/fixtures/annotate_variant/small_CADD.tsv.gz"
CADD_1000G_FILE = "tests/fixtures/annotate_variant/small_1000G_CADD.tsv.gz"

from genmod import logger
from genmod.log import init_log
init_log(logger, loglevel="INFO")


def test_genmod_annotate_features():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(
        cli, [
            'annotate',
            VCF_FILE,
            '-r'
        ])

    assert result.exit_code == 0

def test_genmod_annotate_thousand_g():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(
        cli, [
            'annotate',
            VCF_FILE,
            '--thousand-g',
            THOUSAND_G_FILE
        ])
    
    assert result.exit_code == 0

def test_genmod_annotate_cadd():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(
        cli, [
            'annotate',
            VCF_FILE,
            '--cadd-file',
            CADD_FILE
        ])
    
    assert result.exit_code == 0

def test_genmod_annotate_multiple_cadd():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(
        cli, [
           'annotate',
            VCF_FILE,
            '--cadd-file',
            CADD_FILE,
            '--cadd-file',
            CADD_1000G_FILE
            
        ])
    
    assert result.exit_code == 0

