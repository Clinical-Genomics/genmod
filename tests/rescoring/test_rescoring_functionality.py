from tempfile import NamedTemporaryFile
from typing import Dict, List

from click.testing import CliRunner
from genmod.commands import score_command, score_compounds_command
from genmod.vcf_tools import HeaderParser, get_info_dict, get_variant_dict, get_variant_id

SCORED_VCF = "tests/fixtures/test_vcf_annotated_scored.vcf"
SCORE_CONFIG = "tests/fixtures/score_variants/genmod_example.ini"


def _flatten_dict(d: Dict) -> Dict:
    """
    Flatten dictionary
    """
    output = d.copy()
    for key, value in d.items():
        if isinstance(value, dict):
            output = d | value
    return output


def _yield_variants_from_vcf(vcf_file_path: str):
    head = HeaderParser()
    # Parse file header
    with open(vcf_file_path, "r") as fp:
        for line in fp:
            line = line.rstrip()
            if line.startswith("#"):
                if line.startswith("##"):
                    head.parse_meta_data(line)
                else:
                    head.parse_header_line(line)
            else:
                break
    # Yield variants
    with open(vcf_file_path, "r") as fp:
        try:
            for line in fp:
                if line.startswith("#"):
                    continue
                variant = get_variant_dict(line, head.header)
                variant_id = get_variant_id(variant)
                variant["variant_id"] = variant_id
                variant["info_dict"] = get_info_dict(variant["INFO"])
                variant = _flatten_dict(variant)
                yield variant
        except Exception as e:
            print(f"Failed parsing VCF line; '{line}'")
            raise e


def _sort_variants_on_id(variants: List[Dict]) -> List[Dict]:
    return sorted(variants, key=lambda variant: variant["variant_id"])


def test_rescore_with_annotation_suffix():
    """
    Test for asserting secondary Genmod annotation is present in VCF.
    This is intended for secondary scoring using a different rank config,
    and where score annotations are added with a suffix appended.

    Applies only to:
        - Variant scoring
        - Compound scoring
    remaining annotations, models etc are out of scope for this test.

    Expected:
        - All previously existing GENMOD annotations are present and unaltered
        - New scoring annotations are added that's based on an additional scoring config
    """
    SUFFIX = "Suffix"

    # Variant Scoring with annotation suffix
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
    rescored_variants = NamedTemporaryFile(suffix=".vcf")
    with open(rescored_variants.name, "w") as file:
        file.write(result.stdout_bytes.decode("utf-8"))  # Save processed VCF to file

    # Compound Scoring with annotation suffix
    runner = CliRunner()
    # WHEN computing compound score
    result_compounds = runner.invoke(
        score_compounds_command,
        [rescored_variants.name, "--annotation_suffix", SUFFIX],
    )
    assert result_compounds.exit_code == 0
    rescored_variants_with_compounds = NamedTemporaryFile(suffix=".vcf")
    with open(rescored_variants_with_compounds.name, "w") as file:
        file.write(result_compounds.stdout_bytes.decode("utf-8"))  # Save processed VCF to file

    # WHEN comparing the variants pre and post suffix annotation
    original_variants = []
    for variant in _yield_variants_from_vcf(SCORED_VCF):
        original_variants.append(variant)
    original_variants = _sort_variants_on_id(original_variants)

    rescored_variants = []
    for variant in _yield_variants_from_vcf(rescored_variants_with_compounds.name):
        rescored_variants.append(variant)
    rescored_variants = _sort_variants_on_id(rescored_variants)

    assert len(original_variants) == len(rescored_variants), (
        len(original_variants),
        len(rescored_variants),
    )

    unaltered_annotations = ["RankScore", "RankScoreNormalized", "RankScoreMinMax", "Compounds"]
    added_annotations = [f"{name}{SUFFIX}" for name in unaltered_annotations]

    # THEN expect that the previously annotated information is unchanged ...
    for unalted_annotation in unaltered_annotations:
        for original_variant, rescored_variant in zip(original_variants, rescored_variants):
            # Compounds annotation not present for all variants (continue if missing)
            if (
                unalted_annotation == "Compounds"
                and unalted_annotation not in original_variant.keys()
            ):
                continue
            assert original_variant[unalted_annotation] == rescored_variant[unalted_annotation], (
                original_variant["variant_id"],
                rescored_variant["variant_id"],
                unalted_annotation,
                original_variant[unalted_annotation],
                rescored_variant[unalted_annotation],
            )

    # THEN expect that new annotations are added
    for added_annotation in added_annotations:
        for original_variant, rescored_variant in zip(original_variants, rescored_variants):
            # Compounds annotation not present for all variants (continue if missing)
            if (
                added_annotation.replace(SUFFIX, "") == "Compounds"
                and added_annotation not in original_variant.keys()
            ):
                continue
            assert added_annotation not in original_variant.keys(), (
                added_annotation,
                original_variant,
            )
            assert added_annotation in rescored_variant.keys(), (added_annotation, rescored_variant)
