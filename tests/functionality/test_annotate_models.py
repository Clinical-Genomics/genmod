from tempfile import NamedTemporaryFile

from click.testing import CliRunner
from genmod import logger
from genmod.commands import models_command
from genmod.log import init_log
from test_utils import generate_variants_from_file

ANNOTATED_VCF_FILE = "tests/fixtures/test_vcf_annotated.vcf"
VCF_FILE = "tests/fixtures/test_vcf_regions.vcf"
VCF_FILE_WITH_CHR = "tests/fixtures/test_vcf_regions_with_chr.vcf"
FAMILY_FILE = "tests/fixtures/recessive_trio.ped"
BAD_FAMILY_FILE = "tests/fixtures/annotate_models/one_ind.ped"
EMPTY_VCF_FILE = "tests/fixtures/empty.vcf"

init_log(logger, loglevel="INFO")


def _generate_genetic_models_string_from_file(file_path: str) -> str:
    """
    Yield genetic model string from VCF.
    :param file_path: VCF to be read
    """
    for variant in generate_variants_from_file(file_path=file_path):
        genetic_models_entry: str = variant["info_dict"].get("GeneticModels", "")
        for family_genetic_models in genetic_models_entry.split(","):
            family_genetic_models = family_genetic_models.split(":")
            if len(family_genetic_models) > 1:  # Not all variants will have a model
                genetic_models: str = family_genetic_models[1]
                yield genetic_models


def _run_model_command(vcf_file):
    """Helper function to run models_command and return output as a list."""
    runner = CliRunner()
    result = runner.invoke(models_command, [vcf_file, "-f", FAMILY_FILE])
    assert result.exit_code == 0, f"Command failed with exit code: {result.exit_code}"

    with NamedTemporaryFile(delete=False) as temp_file:
        temp_file.write(result.stdout_bytes)
        temp_file.seek(0)  # Move back to the start of the file
        return list(_generate_genetic_models_string_from_file(temp_file.name))


def test_genmod_annotate_models_no_family():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(models_command, [VCF_FILE])

    # This should fail since there is no family file
    assert result.exit_code == 1


def test_genmod_annotate_models():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(models_command, [VCF_FILE, "-f", FAMILY_FILE])

    print(result.output)
    assert result.exit_code == 0


def test_genmod_annotate_models_empty_vcf():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(models_command, [EMPTY_VCF_FILE, "-f", FAMILY_FILE])

    print(result.output)
    assert result.exit_code == 0


def test_annotate_models_already_annotated():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(models_command, [ANNOTATED_VCF_FILE, "-f", FAMILY_FILE])

    assert result.exit_code == 1


def test_annotate_models_lacking_ind():
    """docstring for test_genmod_annotate_models"""
    runner = CliRunner()
    result = runner.invoke(models_command, [VCF_FILE, "-f", BAD_FAMILY_FILE])

    assert result.exit_code == 1


def test_annotate_models_chr_prefix():
    """Test that genetic models are identical for VCF with and without 'chr' prefix."""

    # Get models from both VCF files
    models_list = _run_model_command(VCF_FILE)
    models_list_with_chr = _run_model_command(VCF_FILE_WITH_CHR)

    # Assert that the lists of models are identical
    assert len(models_list) > 0 and len(models_list_with_chr) > 0, "No models in VCFs"
    assert models_list == models_list_with_chr, "Models differ between VCF files."
