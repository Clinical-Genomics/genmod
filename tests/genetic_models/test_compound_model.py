from genmod.annotate_models.models.compound_model import (
    check_compounds,
    get_genotype,
    get_phase_set,
    variants_on_same_allele,
)
from genmod.vcf_tools import Genotype
from ped_parser import FamilyParser


def get_family(family_lines):
    return FamilyParser(family_lines)


def make_variant(*, genotypes, format_str="GT:PS", sample_fields=None):
    variant = {"genotypes": dict(genotypes), "FORMAT": format_str}
    if sample_fields:
        variant.update(sample_fields)
    return variant


def test_get_genotype_returns_expected_object():
    genotype = Genotype(**{"GT": "0/1"})
    variant = {"genotypes": {"sample": genotype}}
    assert get_genotype(variant, "sample") is genotype


def test_get_phase_set_returns_ps_when_present():
    variant = {"FORMAT": "GT:AD:PS", "sample": "0/1:10,5:12345"}
    assert get_phase_set(variant, "sample") == "12345"


def test_get_phase_set_returns_none_when_ps_missing_in_format():
    variant = {"FORMAT": "GT:AD", "sample": "0/1:10,5"}
    assert get_phase_set(variant, "sample") is None


def test_get_phase_set_returns_none_when_sample_missing_or_unset():
    variant_missing_sample = {"FORMAT": "GT:PS"}
    assert get_phase_set(variant_missing_sample, "sample") is None

    variant_empty_sample = {"FORMAT": "GT:PS", "sample": ""}
    assert get_phase_set(variant_empty_sample, "sample") is None

    variant_short_sample = {"FORMAT": "GT:PS", "sample": "0|1"}
    assert get_phase_set(variant_short_sample, "sample") is None


def test_variants_on_same_allele_true_same_ps_and_overlap():
    variant_1 = make_variant(
        genotypes={"sample": Genotype(**{"GT": "0|1"})},
        sample_fields={"sample": "0|1:PS1"},
    )
    variant_2 = make_variant(
        genotypes={"sample": Genotype(**{"GT": "0|1"})},
        sample_fields={"sample": "0|1:PS1"},
    )
    assert variants_on_same_allele("sample", variant_1, variant_2) is True


def test_variants_on_same_allele_false_different_phase_set():
    variant_1 = make_variant(
        genotypes={"sample": Genotype(**{"GT": "0|1"})},
        sample_fields={"sample": "0|1:1"},
    )
    variant_2 = make_variant(
        genotypes={"sample": Genotype(**{"GT": "0|1"})},
        sample_fields={"sample": "0|1:2"},
    )
    assert variants_on_same_allele("sample", variant_1, variant_2) is False


def test_variants_on_same_allele_false_same_phase_set_no_overlap():
    variant_1 = make_variant(
        genotypes={"sample": Genotype(**{"GT": "0|1"})},
        sample_fields={"sample": "0|1:1"},
    )
    variant_2 = make_variant(
        genotypes={"sample": Genotype(**{"GT": "1|0"})},
        sample_fields={"sample": "1|0:1"},
    )
    assert variants_on_same_allele("sample", variant_1, variant_2) is False


def test_check_compounds_false_when_affected_phased_same_allele():
    family_lines = [
        "#FamilyID\tSampleID\tFather\tMother\tSex\tPhenotype\n",
        "1\tchild\tfather\tmother\t1\t2\n",
        "1\tfather\t0\t0\t1\t1\n",
        "1\tmother\t0\t0\t2\t1\n",
    ]
    family = get_family(family_lines)

    variant_1 = make_variant(
        genotypes={
            "child": Genotype(**{"GT": "0|1"}),
            "father": Genotype(**{"GT": "0|0"}),
            "mother": Genotype(**{"GT": "0|0"}),
        },
        sample_fields={"child": "0|1:1"},
    )
    variant_2 = make_variant(
        genotypes={
            "child": Genotype(**{"GT": "0|1"}),
            "father": Genotype(**{"GT": "0|0"}),
            "mother": Genotype(**{"GT": "0|0"}),
        },
        sample_fields={"child": "0|1:1"},
    )

    assert check_compounds(variant_1, variant_2, family=family, phased=True) is False


def test_check_compounds_false_when_healthy_parent_has_both_variants():
    family_lines = [
        "#FamilyID\tSampleID\tFather\tMother\tSex\tPhenotype\n",
        "1\tchild\tfather\tmother\t1\t2\n",
        "1\tfather\t0\t0\t1\t1\n",
        "1\tmother\t0\t0\t2\t1\n",
    ]
    family = get_family(family_lines)

    variant_1 = make_variant(
        genotypes={
            "child": Genotype(**{"GT": "0/1"}),
            "father": Genotype(**{"GT": "0/1"}),
            "mother": Genotype(**{"GT": "0/0"}),
        }
    )
    variant_2 = make_variant(
        genotypes={
            "child": Genotype(**{"GT": "0/1"}),
            "father": Genotype(**{"GT": "0/1"}),
            "mother": Genotype(**{"GT": "0/0"}),
        }
    )

    assert check_compounds(variant_1, variant_2, family=family, phased=False) is False
