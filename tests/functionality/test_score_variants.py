from genmod.commands import score_command
from click.testing import CliRunner

ANNOTATED_VCF_FILE = "tests/fixtures/test_vcf_annotated.vcf"
VCF_FILE = "tests/fixtures/test_vcf_regions.vcf"
EMPTY_VCF_FILE = "tests/fixtures/test_vcf_annotated_empty.vcf"
SCORED_VCF = "tests/fixtures/test_vcf_annotated_scored.vcf"
SCORE_CONFIG = "tests/fixtures/score_variants/genmod_example.ini"
BAD_FAMILY_FILE = "tests/fixtures/annotate_models/one_ind.ped"

from genmod import logger
from genmod.log import init_log
init_log(logger, loglevel="INFO")


def test_genmod_score_no_config():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(score_command, [ANNOTATED_VCF_FILE])
    
    # This should fail since there is no family file
    assert result.exit_code == 1

def test_genmod_score():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(score_command, [
        ANNOTATED_VCF_FILE,
        '-c',
        SCORE_CONFIG
        ]
    )

    print(result.output)
    assert result.exit_code == 0

def test_genmod_score_empty_vcf():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(score_command, [
        EMPTY_VCF_FILE,
        '-c',
        SCORE_CONFIG
        ]
    )

    print(result.output)
    assert result.exit_code == 0

def test_annotate_models_already_scored():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(score_command, [
        SCORED_VCF,
        '-c',
        SCORE_CONFIG
        ]
    )

    assert result.exit_code == 1
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