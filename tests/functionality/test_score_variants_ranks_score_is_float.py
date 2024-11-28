from tempfile import NamedTemporaryFile

import pytest
from click.testing import CliRunner
from genmod.commands import score_command, score_compounds_command
from test_utils import generate_variants_from_file

ANNOTATED_VCF_FILE = "tests/fixtures/test_vcf_annotated.vcf"
SCORE_CONFIG = "tests/fixtures/score_variants/genmod_example.ini"


def _generate_rank_score_strings_from_file(file_path: str) -> str:
    """
    Yield rank score strings from VCF.
    :param file_path: VCF to be read
    """
    for variant in generate_variants_from_file(file_path=file_path):
        rank_score_entry: str = variant["info_dict"].get("RankScore", "")
        for family_rank_score in rank_score_entry.split(","):
            family_rank_score = family_rank_score.split(":")
            rank_score: str = family_rank_score[1]
            yield rank_score


def _check_rankscore_string_is_float(rank_score_string: str):
    """
    Check to see if a rank score string contains a floating point value.
    :param rank_score_string: A string to be checked
    :raises AssertionError: In case rank_score_string is not a float-type value
    """
    # Check decimal point presence
    assert "." in rank_score_string
    # Check all numerical digits, might contain period and minus sign.
    assert rank_score_string.replace(".", "").replace("-", "").isdigit()
    # Check successfully parsed to float
    assert isinstance(float(rank_score_string), float)


def test_check_rankscore_string_is_float():
    """
    Test for checking integrity of rank score float check method.
    """
    # GIVEN some rank score strings
    # WHEN running the method under test
    # THEN expect it to behave as expected in positive case
    _check_rankscore_string_is_float("0.0")
    _check_rankscore_string_is_float("132.1")
    _check_rankscore_string_is_float("132.10")
    _check_rankscore_string_is_float("-10.0")

    # THEN expect it to behave as expected in negative case
    with pytest.raises(AssertionError):
        _check_rankscore_string_is_float("")
    with pytest.raises(AssertionError):
        _check_rankscore_string_is_float("132")
    with pytest.raises(AssertionError):
        _check_rankscore_string_is_float("b")
    with pytest.raises(AssertionError):
        _check_rankscore_string_is_float("0")
    with pytest.raises(AssertionError):
        _check_rankscore_string_is_float(".")
    with pytest.raises(AssertionError):
        _check_rankscore_string_is_float("-")


def test_rankscore_is_float_type():
    """
    Test to make sure RankScore value is float
    provided by both the CLI command 'score' and 'compound'.
    """
    # GIVEN some VCF file to be ranked
    runner = CliRunner()
    # WHEN computing rank score
    result = runner.invoke(score_command, [ANNOTATED_VCF_FILE, "-c", SCORE_CONFIG])
    assert result.exit_code == 0
    temporary_file = NamedTemporaryFile()
    with open(temporary_file.name, "w") as file:
        file.write(result.stdout_bytes.decode("utf-8"))  # Save processed VCF to file
    # THEN expect all rank scores to be float type
    for rank_score_string in _generate_rank_score_strings_from_file(file_path=temporary_file.name):
        _check_rankscore_string_is_float(rank_score_string)

    # GIVEN some ranked VCF file, run compound scoring (which modify the RankScore)
    runner = CliRunner()
    # WHEN computing compound score
    result = runner.invoke(
        score_compounds_command,
        [
            temporary_file.name,
        ],
    )
    assert result.exit_code == 0
    temporary_file = NamedTemporaryFile()
    with open(temporary_file.name, "w") as file:
        file.write(result.stdout_bytes.decode("utf-8"))  # Save processed VCF to file
    # THEN expect all rank scores (including modified compound scores) to be float type
    for rank_score_string in _generate_rank_score_strings_from_file(file_path=temporary_file.name):
        _check_rankscore_string_is_float(rank_score_string)
