"""Scoring with a CSQ column that is missing from the VEP annotation (genmod issue #85).

The rank model scores a CSQ column, ``NOT_IN_VEP``, that is not in the ``Format:`` of the CSQ
INFO header of the VCF. Without ``--skip_plugin_check`` genmod refuses to score. With
``--skip_plugin_check`` the user accepts missing fields, so the column should score as
not reported, instead of crashing with ``UnboundLocalError: csq_column``.
"""

from pathlib import Path

from click.testing import CliRunner
from genmod.commands import score_command
from genmod.plugin import Plugin
from test_utils import generate_variants_from_file

CLINVAR_VCF = "tests/fixtures/score_variants/clinvar_csq.vcf"
EXACT_CONFIG = "tests/fixtures/score_variants/clinvar_string_exact.ini"
MISSING_FIELD_CONFIG = "tests/fixtures/score_variants/clinvar_missing_csq_field.ini"

CSQ_FORMAT = [
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Feature",
    "SIFT",
    "PolyPhen",
    "CLNVID",
    "CLNSIG",
    "CLNVID_2",
    "CLNREVSTAT",
]
CSQ_ENTRY = "G|missense_variant|MODERATE|VCAN|ENST00000265077|tolerated|benign|1|Pathogenic|1|x"


def rank_scores(args: list[str], tmp_path: Path) -> dict[str, float]:
    """Run `genmod score` and return the rank score of each variant keyed by CHROM:POS."""
    result = CliRunner().invoke(score_command, args)
    assert result.exit_code == 0, f"{result.output}\n{result.exception!r}"

    scored_vcf = tmp_path / "scored.vcf"
    scored_vcf.write_bytes(result.stdout_bytes)
    return {
        f"{variant['CHROM']}:{variant['POS']}": float(
            variant["info_dict"]["RankScore"].split(":")[1]
        )
        for variant in generate_variants_from_file(str(scored_vcf))
    }


def missing_csq_plugin() -> Plugin:
    """A string plugin for a CSQ column that is not in the CSQ format of the test VCF."""
    return Plugin(
        name="NOT_IN_VEP",
        field="INFO",
        data_type="string",
        info_key="CSQ",
        csq_key="NOT_IN_VEP",
        separators=["&"],
        record_rule="max",
        string_rules={"present": 1},
    )


def test_plugin_value_is_none_when_csq_key_is_not_in_the_csq_format() -> None:
    # GIVEN a plugin for a CSQ column that is not in the CSQ format, and a variant with CSQ
    plugin = missing_csq_plugin()
    variant = {"info_dict": {"CSQ": CSQ_ENTRY}}

    # WHEN getting the value
    value = plugin.get_value(variant_dict=variant, csq_format=CSQ_FORMAT)

    # THEN no value is found, instead of an UnboundLocalError on `csq_column`
    assert value is None


def test_plugin_entry_is_empty_when_csq_key_is_not_in_the_csq_format() -> None:
    # GIVEN a plugin for a CSQ column that is not in the CSQ format
    plugin = missing_csq_plugin()

    # WHEN getting the entry of a variant that has a CSQ annotation
    entry = plugin.get_entry(raw_entry=CSQ_ENTRY, csq_format=CSQ_FORMAT)

    # THEN the entry is empty, as for any other entry that was not found
    assert entry == []


def test_score_without_skip_plugin_check_aborts_on_missing_csq_field() -> None:
    # GIVEN a rank model scoring a CSQ column that is not in the VCF header
    # WHEN scoring without --skip_plugin_check
    result = CliRunner().invoke(score_command, [CLINVAR_VCF, "-c", MISSING_FIELD_CONFIG])

    # THEN genmod aborts before scoring, and does not crash
    assert result.exit_code != 0
    assert not isinstance(result.exception, UnboundLocalError)


def test_score_with_skip_plugin_check_scores_missing_csq_field_as_not_reported(
    tmp_path: Path,
) -> None:
    # GIVEN the ClinVar scores without the missing CSQ column in the rank model
    expected = rank_scores([CLINVAR_VCF, "-c", EXACT_CONFIG], tmp_path)
    assert expected

    # WHEN scoring with the missing CSQ column in the rank model and --skip_plugin_check
    scores = rank_scores([CLINVAR_VCF, "-c", MISSING_FIELD_CONFIG, "--skip_plugin_check"], tmp_path)

    # THEN the missing column adds nothing to the score of any variant
    assert scores == expected
