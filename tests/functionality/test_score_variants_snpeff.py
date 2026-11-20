"""Scoring SnpEff ANN fields with a rank model (genmod issue #96).

A plugin with ``info_key = ANN`` scores the ANN field named by ``csq_key``, like a CSQ plugin
scores a VEP column.
"""

from pathlib import Path
from typing import Dict

from click.testing import CliRunner
from genmod.commands import score_command
from test_utils import generate_variants_from_file

SNPEFF_VCF = "tests/fixtures/snpeff/snpeff_trio.vcf"
RANK_MODEL = "tests/fixtures/snpeff/snpeff_rank_model.ini"
MISSING_FIELD_RANK_MODEL = "tests/fixtures/snpeff/snpeff_missing_ann_field.ini"

# Rank score per position: Annotation_Impact score + Annotation score
EXPECTED_SCORES = {
    "879537": 6.0,  # MODERATE (3) + missense_variant (3), the MODIFIER annotation is not the max
    "879595": 7.0,  # MODERATE (3) + splice_region_variant (4) from missense&splice_region
    "880086": 2.0,  # LOW (1) + synonymous_variant (1)
    "900000": 0.0,  # MODIFIER (0) + intergenic_region, not in the rank model (0)
    "950000": 10.0,  # HIGH (5) + stop_gained (5), from the second ALT allele
    "960000": 5.0,  # MODERATE (3) + conservative_inframe_deletion (2)
}


def rank_scores(args: list[str], tmp_path: Path) -> Dict[str, float]:
    """Run `genmod score` and return the rank score of each variant keyed by position."""
    result = CliRunner().invoke(score_command, args)
    assert result.exit_code == 0, f"{result.output}\n{result.exception!r}"

    scored_vcf = tmp_path / "scored.vcf"
    scored_vcf.write_bytes(result.stdout_bytes)
    return {
        variant["POS"]: float(variant["info_dict"]["RankScore"].split(":")[1])
        for variant in generate_variants_from_file(str(scored_vcf))
    }


def test_score_snpeff_ann_fields(tmp_path: Path) -> None:
    # GIVEN a rank model scoring the ANN fields Annotation_Impact and Annotation
    # WHEN scoring the SnpEff annotated VCF
    scores = rank_scores([SNPEFF_VCF, "-c", RANK_MODEL], tmp_path)

    # THEN each variant gets the score of its most severe impact and consequence
    assert scores == EXPECTED_SCORES


def test_score_without_skip_plugin_check_aborts_on_missing_ann_field() -> None:
    # GIVEN a rank model scoring an ANN field that is not in the ANN header
    # WHEN scoring without --skip_plugin_check
    result = CliRunner().invoke(score_command, [SNPEFF_VCF, "-c", MISSING_FIELD_RANK_MODEL])

    # THEN genmod aborts before scoring
    assert result.exit_code != 0


def test_score_with_skip_plugin_check_scores_missing_ann_field_as_not_reported(
    tmp_path: Path,
) -> None:
    # GIVEN a rank model scoring an ANN field that is not in the ANN header
    # WHEN scoring with --skip_plugin_check
    scores = rank_scores(
        [SNPEFF_VCF, "-c", MISSING_FIELD_RANK_MODEL, "--skip_plugin_check"], tmp_path
    )

    # THEN the missing field adds nothing to the score of any variant
    assert scores == EXPECTED_SCORES
