"""Genetic models for SnpEff annotated VCFs (genmod issue #96).

`genmod models --snpeff` takes the genes of each variant from the SnpEff ``ANN`` field, as
`--vep` does from the VEP ``CSQ`` field. ``snpeff_trio.vcf`` and ``vep_trio.vcf`` hold the
same synthetic variants, annotated as SnpEff ANN and as VEP CSQ.
"""

from pathlib import Path
from typing import Dict

from click.testing import CliRunner
from genmod.commands import models_command
from test_utils import generate_variants_from_file

SNPEFF_VCF = "tests/fixtures/snpeff/snpeff_trio.vcf"
CSQ_VCF = "tests/fixtures/snpeff/vep_trio.vcf"
NO_ANNOTATION_VCF = "tests/fixtures/test_vcf_regions.vcf"
FAMILY_FILE = "tests/fixtures/recessive_trio.ped"


def sort_pairs(value: str) -> str:
    """Sort the '|'-separated items of a 'family:item|item' value, their order is not stable."""
    family, _, items = value.partition(":")
    return f"{family}:{'|'.join(sorted(items.split('|')))}"


def model_annotations(args: list[str], tmp_path: Path) -> Dict[str, Dict[str, str]]:
    """Run `genmod models` and return GeneticModels and Compounds of each variant by position."""
    result = CliRunner().invoke(models_command, args)
    assert result.exit_code == 0, f"{result.output}\n{result.exception!r}"

    annotated_vcf = tmp_path / "annotated.vcf"
    annotated_vcf.write_bytes(result.stdout_bytes)
    models = {
        variant["POS"]: {
            key: sort_pairs(variant["info_dict"][key])
            for key in ("GeneticModels", "Compounds")
            if key in variant["info_dict"]
        }
        for variant in generate_variants_from_file(str(annotated_vcf))
    }
    # Guard against comparing empty outputs
    assert len(models) == 6, models
    return models


def test_snpeff_finds_compound_in_shared_gene(tmp_path: Path) -> None:
    # GIVEN a trio where the proband has one paternal and one maternal variant in SAMD11
    # WHEN annotating models with the genes from the SnpEff ANN field
    models = model_annotations([SNPEFF_VCF, "-f", FAMILY_FILE, "--snpeff"], tmp_path)

    # THEN the two variants are a compound pair
    assert "AR_comp" in models["879537"]["GeneticModels"]
    assert "1_879595_C_T" in models["879537"]["Compounds"]
    assert models["879595"] == {"GeneticModels": "1:AR_comp", "Compounds": "1:1_879537_T_C"}


def test_snpeff_intergenic_variant_has_no_compounds(tmp_path: Path) -> None:
    # GIVEN an intergenic variant, where Gene_ID holds the two closest genes
    # WHEN annotating models with the SnpEff ANN field
    models = model_annotations([SNPEFF_VCF, "-f", FAMILY_FILE, "--snpeff"], tmp_path)

    # THEN it does not belong to any gene, so it is in no compound pair
    assert models["900000"] == {"GeneticModels": "1:AR_hom"}


def test_snpeff_gives_same_models_as_vep(tmp_path: Path) -> None:
    # GIVEN the same variants annotated by SnpEff (ANN) and VEP (CSQ)
    # WHEN annotating models with each annotation
    snpeff = model_annotations([SNPEFF_VCF, "-f", FAMILY_FILE, "--snpeff"], tmp_path)
    vep = model_annotations([CSQ_VCF, "-f", FAMILY_FILE, "--vep"], tmp_path)

    # THEN the genetic models and compounds are identical
    assert snpeff == vep


def test_snpeff_without_ann_header_aborts() -> None:
    # GIVEN a VCF without an ANN INFO header
    # WHEN annotating models with --snpeff
    result = CliRunner().invoke(models_command, [NO_ANNOTATION_VCF, "-f", FAMILY_FILE, "--snpeff"])

    # THEN genmod aborts
    assert result.exit_code != 0


def test_vep_and_snpeff_together_is_a_usage_error() -> None:
    # GIVEN both --vep and --snpeff
    # WHEN annotating models
    result = CliRunner().invoke(
        models_command, [SNPEFF_VCF, "-f", FAMILY_FILE, "--vep", "--snpeff"]
    )

    # THEN it is a usage error
    assert result.exit_code == 2
    assert "--vep and --snpeff can not be used together" in result.output
