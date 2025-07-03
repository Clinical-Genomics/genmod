import pytest as pt

from tempfile import NamedTemporaryFile
from typing import Dict, List

from click.testing import CliRunner
from genmod.commands import score_command, score_compounds_command
from genmod.vcf_tools import HeaderParser, get_info_dict, get_variant_dict, get_variant_id


SCORED_VCF = "tests/fixtures/test_vcf_annotated_scored.vcf"
SCORE_CONFIG = "tests/fixtures/score_variants/rank_model_-v1.39-mivmir.ini"
SUFFIX = "MIVMIR"

def test_mivmir_minimal_score_config():
    runner = CliRunner()
    result = runner.invoke(
        score_command,
        [
            SCORED_VCF,
            "-c",
            SCORE_CONFIG,
            "--skip_is_previously_scored_check",
            "--annotation_suffix",
            SUFFIX,
        ],
    )
    assert result.exit_code == 0
    print(result.stdout_bytes.decode("utf-8"))
    assert False