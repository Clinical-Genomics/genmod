from click.testing import CliRunner
from tempfile import NamedTemporaryFile

from genmod import logger
from genmod.commands import score_command, score_compounds_command
from genmod.log import init_log

SCORED_VCF = "tests/fixtures/test_vcf_annotated_scored.vcf"
SCORE_CONFIG = "tests/fixtures/score_variants/genmod_example.ini"

init_log(logger, loglevel="INFO")

def test_rescore_with_annotation_suffix():
    """
    Test for asserting secondary Genmod annotation is present in VCF.
    This is intended for secondary analysis using a rank config,
    and where score annotations are added with a suffix.

    Expected:
        - All previously existing GENMOD annotations are present and unaltered
        - New scoring annotations are added that's based on an additional scoring config
    """

    # Variant Scoring
    runner = CliRunner()
    result = runner.invoke(score_command,
                           [SCORED_VCF,
                            "-c", SCORE_CONFIG,
                            "--skip_is_previously_scored_check"])
    assert result.exit_code == 0
    temporary_file = NamedTemporaryFile()
    with open(temporary_file.name, "w") as file:
        file.write(result.stdout_bytes.decode("utf-8"))  # Save processed VCF to file

    # Compound Scoring
    runner = CliRunner()
    # WHEN computing compound score
    result = runner.invoke(
        score_compounds_command,
        [temporary_file.name],
    )
    assert result.exit_code == 0
    temporary_file = NamedTemporaryFile()
    with open(temporary_file.name, "w") as file:
        file.write(result.stdout_bytes.decode("utf-8"))  # Save processed VCF to file