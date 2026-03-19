from tempfile import NamedTemporaryFile

from click.testing import CliRunner
from genmod import logger
from genmod.commands import score_command
from genmod.log import init_log
from test_utils import generate_variants_from_file

ANNOTATED_VCF_FILE = "tests/fixtures/test_vcf_annotated.vcf"
VCF_FILE = "tests/fixtures/test_vcf_regions.vcf"
EMPTY_VCF_FILE = "tests/fixtures/test_vcf_annotated_empty.vcf"
SCORED_VCF = "tests/fixtures/test_vcf_annotated_scored.vcf"
SCORE_CONFIG = "tests/fixtures/score_variants/genmod_example.ini"
BAD_FAMILY_FILE = "tests/fixtures/annotate_models/one_ind.ped"
SV_SAME_POS_VCF_FILE = "tests/fixtures/test_vcf_sv_same_pos_annotated.vcf"

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
    result = runner.invoke(score_command, [ANNOTATED_VCF_FILE, "-c", SCORE_CONFIG])

    print(result.output)
    assert result.exit_code == 0


def test_genmod_score_empty_vcf():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(score_command, [EMPTY_VCF_FILE, "-c", SCORE_CONFIG])

    print(result.output)
    assert result.exit_code == 0


def test_annotate_models_already_scored():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(score_command, [SCORED_VCF, "-c", SCORE_CONFIG])

    assert result.exit_code == 1


def test_genmod_score_same_pos_sv_keeps_distinct_end_variants():
    """Test scoring keeps symbolic SVs at same locus as distinct variants."""
    runner = CliRunner()
    result = runner.invoke(score_command, [SV_SAME_POS_VCF_FILE, "-c", SCORE_CONFIG])

    assert result.exit_code == 0

    with NamedTemporaryFile(delete=False) as temp_file:
        temp_file.write(result.stdout_bytes)
        temp_file.seek(0)
        output_variants = list(generate_variants_from_file(temp_file.name))

    assert len(output_variants) == 2
    assert {variant["info_dict"].get("END") for variant in output_variants} == {"1001", "101"}
    assert all("RankScore" in variant["info_dict"] for variant in output_variants)


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
