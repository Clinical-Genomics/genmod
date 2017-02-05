from genmod.commands import models_command
from click.testing import CliRunner

ANNOTATED_VCF_FILE = "tests/fixtures/test_vcf_annotated.vcf"
VCF_FILE = "tests/fixtures/test_vcf_regions.vcf"
FAMILY_FILE = "tests/fixtures/recessive_trio.ped"
BAD_FAMILY_FILE = "tests/fixtures/annotate_models/one_ind.ped"
EMPTY_VCF_FILE = "tests/fixtures/empty.vcf"

from genmod import logger
from genmod.log import init_log
init_log(logger, loglevel="INFO")


def test_genmod_annotate_models_no_family():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(models_command, [VCF_FILE])
    
    # This should fail since there is no family file
    assert result.exit_code == 1

def test_genmod_annotate_models():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(models_command, [
        VCF_FILE, 
        '-f',
        FAMILY_FILE 
        ]
    )
    
    print(result.output)
    assert result.exit_code == 0

def test_genmod_annotate_models_empty_vcf():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(models_command, [
        EMPTY_VCF_FILE, 
        '-f',
        FAMILY_FILE 
        ]
    )
    
    print(result.output)
    assert result.exit_code == 0

def test_annotate_models_already_annotated():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(models_command, [
        ANNOTATED_VCF_FILE, 
        '-f',
        FAMILY_FILE 
        ]
    )
    
    assert result.exit_code == 1

def test_annotate_models_lacking_ind():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(models_command, [
        VCF_FILE, 
        '-f',
        BAD_FAMILY_FILE 
        ]
    )
    
    assert result.exit_code == 1