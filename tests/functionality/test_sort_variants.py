from click.testing import CliRunner
from genmod import logger
from genmod.commands import sort_command
from genmod.log import init_log

ANNOTATED_VCF_FILE = "tests/fixtures/test_vcf_annotated.vcf"
VCF_FILE = "tests/fixtures/test_vcf_regions.vcf"
EMPTY_VCF_FILE = "tests/fixtures/empty.vcf"
SCORED_VCF = "tests/fixtures/test_vcf_annotated_scored.vcf"
SCORE_CONFIG = "tests/fixtures/score_variants/genmod_example.ini"
BAD_FAMILY_FILE = "tests/fixtures/annotate_models/one_ind.ped"

init_log(logger, loglevel="INFO")


# def test_genmod_score_no_config():
#     """docstring for test_genmod_annotate_models"""
#     runner = CliRunner()
#     result = runner.invoke(score_command, [ANNOTATED_VCF_FILE])
#
#     # This should fail since there is no family file
#     assert result.exit_code == 1


def test_genmod_sort():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(
        sort_command,
        [
            SCORED_VCF,
        ],
    )

    print(result.output)
    assert result.exit_code == 0


def test_genmod_sort_empty():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(
        sort_command,
        [
            EMPTY_VCF_FILE,
        ],
    )

    print(result.output)
    assert result.exit_code == 0


def test_genmod_sort_not_scored():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(
        sort_command,
        [
            ANNOTATED_VCF_FILE,
        ],
    )

    print(result.output)
    assert result.exit_code == 0


def test_genmod_sort_position():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(sort_command, [SCORED_VCF, "-p"])

    print(result.output)
    assert result.exit_code == 0


# def test_annotate_models_already_scored():
#     """docstring for test_genmod_annotate_models"""
#     runner = CliRunner()
#     result = runner.invoke(score_command, [
#         SCORED_VCF,
#         '-c',
#         SCORE_CONFIG
#         ]
#     )

# assert result.exit_code == 1
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
