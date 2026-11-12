"""Compound pairs where each parent is only genotyped for one of the variants (issue #99)"""

from dataclasses import dataclass

import pytest
from genmod.annotate_models.genetic_models import check_genetic_models
from genmod.vcf_tools import Genotype
from ped_parser import FamilyParser

TRIO_PED = [
    "#FamilyID\tSampleID\tFather\tMother\tSex\tPhenotype\n",
    "1\tchild\tfather\tmother\t1\t2\n",
    "1\tfather\t0\t0\t1\t1\n",
    "1\tmother\t0\t0\t2\t1\n",
]


@dataclass(frozen=True)
class TrioCalls:
    """The genotype calls of the trio for one variant, e.g. ``0/1`` or ``./.``"""

    child: str
    father: str
    mother: str

    def to_variant(self) -> dict:
        """Return a variant in gene GENE1 with these calls, as used by the genetic models"""
        return {
            "CHROM": "1",
            "annotation": {"GENE1"},
            "genotypes": {
                "child": Genotype(GT=self.child),
                "father": Genotype(GT=self.father),
                "mother": Genotype(GT=self.mother),
            },
        }


def run_models(variant_1: dict, variant_2: dict, strict: bool = False) -> dict:
    """Check the genetic models of two variants in the trio and return the batch"""
    families = FamilyParser(TRIO_PED).families
    batch = {"variant_1": variant_1, "variant_2": variant_2}
    check_genetic_models(batch, families, strict=strict)
    return batch


def test_compound_when_parents_are_called_for_both_variants():
    # GIVEN a healthy trio where the mother carries variant 1 and the father variant 2
    variant_1 = TrioCalls(child="0/1", father="0/0", mother="0/1").to_variant()
    variant_2 = TrioCalls(child="0/1", father="0/1", mother="0/0").to_variant()

    # WHEN checking the genetic models
    batch = run_models(variant_1, variant_2)

    # THEN the variants are an AR compound pair
    assert batch["variant_1"]["inheritance_models"]["1"]["AR_comp"] is True
    assert batch["variant_2"]["inheritance_models"]["1"]["AR_comp"] is True
    assert batch["variant_1"]["compounds"]["1"] == {"variant_2"}
    assert batch["variant_2"]["compounds"]["1"] == {"variant_1"}


@pytest.mark.parametrize("strict", [False, True])
def test_compound_when_each_parent_has_a_no_call_for_one_variant(strict):
    # GIVEN a healthy trio where the mother carries variant 1 but has no call for variant 2,
    # and the father carries variant 2 but has no call for variant 1 (issue #99)
    variant_1 = TrioCalls(child="0/1", father="./.", mother="0/1").to_variant()
    variant_2 = TrioCalls(child="0/1", father="0/1", mother="./.").to_variant()

    # WHEN checking the genetic models
    batch = run_models(variant_1, variant_2, strict=strict)

    # THEN the variants are an AR compound pair, same as when the missing calls are 0/0
    assert batch["variant_1"]["inheritance_models"]["1"]["AR_comp"] is True
    assert batch["variant_2"]["inheritance_models"]["1"]["AR_comp"] is True
    assert batch["variant_1"]["compounds"]["1"] == {"variant_2"}
    assert batch["variant_2"]["compounds"]["1"] == {"variant_1"}


@pytest.mark.parametrize("strict", [False, True])
def test_de_novo_compound_when_one_parent_is_reference_and_the_other_has_no_calls(strict):
    # GIVEN a healthy trio where the father is 0/0 for both variants and the mother has no calls
    variant_1 = TrioCalls(child="0/1", father="0/0", mother="./.").to_variant()
    variant_2 = TrioCalls(child="0/1", father="0/0", mother="./.").to_variant()

    # WHEN checking the genetic models
    batch = run_models(variant_1, variant_2, strict=strict)

    # THEN the father cannot have passed on a variant, so the pair is a de novo compound
    for variant in batch.values():
        assert variant["inheritance_models"]["1"]["AR_comp_dn"] is True
        assert variant["inheritance_models"]["1"]["AR_comp"] is False


def test_compound_when_a_parent_has_no_calls_and_the_other_carries_one_variant():
    # GIVEN a healthy trio where the father carries variant 2 and the mother has no calls
    variant_1 = TrioCalls(child="0/1", father="./.", mother="./.").to_variant()
    variant_2 = TrioCalls(child="0/1", father="0/1", mother="./.").to_variant()

    # WHEN checking the genetic models without strict mode
    batch = run_models(variant_1, variant_2, strict=False)

    # THEN the pair may be a compound, since the mother could carry variant 1
    for variant in batch.values():
        assert variant["inheritance_models"]["1"]["AR_comp"] is True
        assert variant["inheritance_models"]["1"]["AR_comp_dn"] is False


def test_no_compound_in_strict_mode_when_a_parent_carrier_status_is_unknown():
    # GIVEN the same trio, with the mother's carrier status unknown
    variant_1 = TrioCalls(child="0/1", father="./.", mother="./.").to_variant()
    variant_2 = TrioCalls(child="0/1", father="0/1", mother="./.").to_variant()

    # WHEN checking the genetic models in strict mode
    batch = run_models(variant_1, variant_2, strict=True)

    # THEN the pair cannot be confirmed, so no compound model is set
    for variant in batch.values():
        assert variant["inheritance_models"]["1"]["AR_comp"] is False
        assert variant["inheritance_models"]["1"]["AR_comp_dn"] is False
