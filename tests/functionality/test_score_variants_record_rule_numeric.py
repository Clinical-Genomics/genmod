"""Numeric record_rule can choose by score instead of by value (genmod issue #86).

The variant has ``CLNSIG=4|5``. The rank model maps value 4 to score 100 and value 5 to
score 1. By default ``record_rule`` compares the annotation values, so ``max`` takes value 5
(score 1). With ``record_rule_on = score`` it compares the scores, so ``max`` gives 100 and
``min`` gives 1.
"""

from pathlib import Path

from click.testing import CliRunner
from genmod.commands import score_command
from test_utils import generate_variants_from_file

MINIMAL_VCF = "tests/fixtures/score_variants/record_rule_max_numeric.vcf"
MINIMAL_CONFIG = "tests/fixtures/score_variants/record_rule_max_numeric.ini"


def write_config(tmp_path: Path, record_rule: str, record_rule_on: str | None) -> str:
    """Write the rank model with the given record_rule and optional record_rule_on."""
    config = (
        Path(MINIMAL_CONFIG)
        .read_text()
        .replace("record_rule = max", f"record_rule = {record_rule}")
    )
    if record_rule_on:
        config = config.replace(
            f"record_rule = {record_rule}",
            f"record_rule = {record_rule}\n  record_rule_on = {record_rule_on}",
        )
    path = tmp_path / "rankmodel.ini"
    path.write_text(config)
    return str(path)


def rank_score(config: str, tmp_path: Path) -> float:
    """Run `genmod score` for family 'test' and return the rank score of the single variant."""
    result = CliRunner().invoke(score_command, ["-i", "test", MINIMAL_VCF, "-c", config])
    assert result.exit_code == 0, f"{result.output}\n{result.exception!r}"

    scored_vcf = tmp_path / "scored.vcf"
    scored_vcf.write_bytes(result.stdout_bytes)
    variants = list(generate_variants_from_file(str(scored_vcf)))
    assert len(variants) == 1
    return float(variants[0]["info_dict"]["RankScore"].split(":")[1])


def test_record_rule_defaults_to_comparing_values(tmp_path: Path) -> None:
    # GIVEN a numeric plugin with record_rule = max and no record_rule_on
    config = write_config(tmp_path, "max", None)

    # WHEN scoring the variant
    score = rank_score(config, tmp_path)

    # THEN the highest value (5, score 1) is used, as before record_rule_on existed
    assert score == 1


def test_record_rule_max_on_score_uses_the_highest_score(tmp_path: Path) -> None:
    # GIVEN a numeric plugin with record_rule = max and record_rule_on = score
    config = write_config(tmp_path, "max", "score")

    # WHEN scoring the variant
    score = rank_score(config, tmp_path)

    # THEN the highest score (100, from value 4) is used
    assert score == 100


def test_record_rule_min_on_score_uses_the_lowest_score(tmp_path: Path) -> None:
    # GIVEN a numeric plugin with record_rule = min and record_rule_on = score
    config = write_config(tmp_path, "min", "score")

    # WHEN scoring the variant
    score = rank_score(config, tmp_path)

    # THEN the lowest score (1, from value 5) is used
    assert score == 1


def test_record_rule_on_with_a_bad_value_is_rejected(tmp_path: Path) -> None:
    # GIVEN a rank model with an invalid record_rule_on
    config = write_config(tmp_path, "max", "nonsense")

    # WHEN scoring
    result = CliRunner().invoke(score_command, ["-i", "test", MINIMAL_VCF, "-c", config])

    # THEN genmod refuses the config
    assert result.exit_code != 0
