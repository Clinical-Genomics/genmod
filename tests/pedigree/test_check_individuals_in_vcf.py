"""The individuals of a ped file must all be samples in the VCF."""

import pytest
from genmod.pedigree import FamilyParser, check_individuals_in_vcf


def test_all_ped_individuals_in_vcf_passes() -> None:
    # GIVEN ped individuals that are all in the VCF
    ped_individuals = ["1", "2", "3"]
    vcf_individuals = ["1", "2", "3"]

    # WHEN checking
    # THEN nothing is raised
    check_individuals_in_vcf(ped_individuals, vcf_individuals)


def test_vcf_individuals_missing_from_ped_are_allowed() -> None:
    # GIVEN a VCF with more samples than the ped file
    # WHEN checking
    # THEN it passes, only the ped individuals must be in the VCF
    check_individuals_in_vcf(["1"], ["1", "2", "3"])


def test_ped_individual_missing_from_vcf_raises() -> None:
    # GIVEN a ped individual that is not in the VCF
    # WHEN checking
    # THEN an IOError is raised
    with pytest.raises(IOError, match="must exist in VCF"):
        check_individuals_in_vcf(["1", "2", "3", "4"], ["1", "2", "3"])


def test_error_names_the_missing_individuals() -> None:
    # GIVEN two ped individuals that are not in the VCF
    # WHEN checking
    with pytest.raises(IOError) as error:
        check_individuals_in_vcf(["1", "x", "y"], ["1"])

    # THEN the message tells which ones are missing
    assert "x" in str(error.value)
    assert "y" in str(error.value)
    assert "'1'" not in str(error.value)


def test_no_ped_individuals_passes() -> None:
    # GIVEN no ped individuals
    # THEN there is nothing to be missing
    check_individuals_in_vcf([], ["1"])


def test_accepts_the_individuals_of_a_family_parser() -> None:
    # GIVEN the individuals dict of a parsed trio, keyed on individual id
    parser = FamilyParser(
        ["fam\tkid\tdad\tmom\t1\t2\n", "fam\tdad\t0\t0\t1\t1\n", "fam\tmom\t0\t0\t2\t1\n"]
    )

    # WHEN checking against a VCF that has them all
    # THEN it passes
    check_individuals_in_vcf(parser.individuals, ["kid", "dad", "mom"])
    with pytest.raises(IOError):
        check_individuals_in_vcf(parser.individuals, ["kid", "dad"])
