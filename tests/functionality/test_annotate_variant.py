from genmod.commands import annotate_variant_command
from click.testing import CliRunner

VCF_FILE = "tests/fixtures/test_vcf_regions.vcf"
THOUSAND_G_FILE = "tests/fixtures/annotate_variant/small_1000G.vcf.gz"
CADD_FILE = "tests/fixtures/annotate_variant/small_CADD.tsv.gz"
CADD_1000G_FILE = "tests/fixtures/annotate_variant/small_1000G_CADD.tsv.gz"

from genmod import logger
from genmod.log import init_log
init_log(logger, loglevel="INFO")


# def test_genmod_annotate_features():
#     """docstring for test_genmod_annotate_models"""
#     runner = CliRunner()
#     result = runner.invoke(
#         annotate_variant_command, [
#             VCF_FILE,
#             '-r'
#         ])
#
#     # This should fail since there is no family file
#     assert result.exit_code == 0

def test_genmod_annotate_thousand_g():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(
        annotate_variant_command, [
            VCF_FILE,
            '--thousand_g',
            THOUSAND_G_FILE
        ])
    
    # This should fail since there is no family file
    assert result.exit_code == 0

def test_genmod_annotate_cadd():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(
        annotate_variant_command, [
            VCF_FILE,
            '--cadd_file',
            CADD_FILE
        ])
    
    # This should fail since there is no family file
    assert result.exit_code == 0

def test_genmod_annotate_multiple_cadd():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(
        annotate_variant_command, [
            VCF_FILE,
            '--cadd_file',
            CADD_FILE,
            '--cadd_1000g',
            CADD_1000G_FILE
            
        ])
    
    # This should fail since there is no family file
    assert result.exit_code == 0

# def test_genmod_annotate_models():
#     """docstring for test_genmod_annotate_models"""
#     runner = CliRunner()
#     result = runner.invoke(models_command, [
#         VCF_FILE,
#         '-f',
#         FAMILY_FILE
#         ]
#     )
#
#     print(result.output)
#     assert result.exit_code == 0
#
# def test_annotate_models_already_annotated():
#     """docstring for test_genmod_annotate_models"""
#     runner = CliRunner()
#     result = runner.invoke(models_command, [
#         ANNOTATED_VCF_FILE,
#         '-f',
#         FAMILY_FILE
#         ]
#     )
#
#     assert result.exit_code == 1
#
# def test_annotate_models_lacking_ind():
#     """docstring for test_genmod_annotate_models"""
#     runner = CliRunner()
#     result = runner.invoke(models_command, [
#         VCF_FILE,
#         '-f',
#         BAD_FAMILY_FILE
#         ]
#     )
#
#     assert result.exit_code == 1