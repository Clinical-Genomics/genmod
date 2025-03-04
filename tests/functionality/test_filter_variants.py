from click.testing import CliRunner
from genmod import logger
from genmod.commands import filter_command
from genmod.log import init_log

ANNOTATED_VCF_FILE = "tests/fixtures/test_vcf_annotated.vcf"

init_log(logger, loglevel="INFO")


def test_genmod_filter():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(filter_command, [ANNOTATED_VCF_FILE])

    assert result.exit_code == 0
