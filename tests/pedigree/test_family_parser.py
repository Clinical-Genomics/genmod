"""Tests for the pedigree parser (ped, fam and the alternative ped formats)."""

from pathlib import Path

import pytest
from genmod.pedigree import (
    FAMILY_HEADER,
    Family,
    FamilyParser,
    Individual,
    PedigreeError,
    WrongLineFormat,
)

FAMILY = "healthyParentsAffectedSon"

TRIO_LINES = [
    "#Standard trio\n",
    "#FamilyID\tSampleID\tFather\tMother\tSex\tPhenotype\n",
    f"{FAMILY}\tproband\tfather\tmother\t1\t2\n",
    f"{FAMILY}\tmother\t0\t0\t2\t1\n",
    f"{FAMILY}\tfather\t0\t0\t1\t1\n",
]


def test_parse_standard_trio() -> None:
    # GIVEN the lines of a standard trio in ped format
    # WHEN parsing them
    parser = FamilyParser(TRIO_LINES)

    # THEN the family and its three individuals are found
    assert set(parser.families) == {FAMILY}
    family = parser.families[FAMILY]
    assert isinstance(family, Family)
    assert family.family_id == FAMILY
    assert set(family.individuals) == {"proband", "mother", "father"}

    # AND the parent ids, sex and phenotype are read from the right columns
    proband = family.individuals["proband"]
    assert isinstance(proband, Individual)
    assert (proband.father, proband.mother) == ("father", "mother")
    assert (proband.sex, proband.phenotype) == (1, 2)
    assert proband.family == FAMILY
    assert proband.individual_id == "proband"


def test_header_columns_are_ped_order() -> None:
    # THEN the columns of a ped file are family, sample, father, mother, sex, phenotype
    assert FAMILY_HEADER == ["family_id", "sample_id", "father_id", "mother_id", "sex", "phenotype"]


def test_affected_and_healthy_follow_phenotype() -> None:
    # GIVEN a parsed trio
    individuals = FamilyParser(TRIO_LINES).families[FAMILY].individuals

    # THEN phenotype 2 is affected and 1 is healthy
    assert individuals["proband"].affected and not individuals["proband"].healthy
    assert individuals["mother"].healthy and not individuals["mother"].affected


def test_parents_of_individuals() -> None:
    # GIVEN a parsed trio
    individuals = FamilyParser(TRIO_LINES).families[FAMILY].individuals

    # THEN the child has both parents and the founders have none
    assert individuals["proband"].has_parents
    assert individuals["proband"].has_both_parents
    assert not individuals["mother"].has_parents
    assert not individuals["mother"].has_both_parents


def test_duo_has_parent_but_not_both_parents() -> None:
    # GIVEN a child with a mother but no father
    lines = [
        f"{FAMILY}\tchild\t0\tmother\t1\t2\n",
        f"{FAMILY}\tmother\t0\t0\t2\t1\n",
    ]

    # WHEN parsing
    child = FamilyParser(lines).families[FAMILY].individuals["child"]

    # THEN the child has parents but not both
    assert child.has_parents
    assert not child.has_both_parents
    assert child.father == "0"


def test_individuals_are_available_on_the_parser() -> None:
    # GIVEN two families
    lines = TRIO_LINES + ["fam2\tsolo\t0\t0\t2\t2\n"]

    # WHEN parsing
    parser = FamilyParser(lines)

    # THEN individuals of all families can be looked up on the parser
    assert set(parser.families) == {FAMILY, "fam2"}
    assert set(parser.individuals) == {"proband", "mother", "father", "solo"}
    assert parser.individuals["solo"].family == "fam2"


def test_extra_daughter() -> None:
    # GIVEN a trio with an extra daughter
    lines = TRIO_LINES + [f"{FAMILY}\tdaughter\tfather\tmother\t2\t1\n"]

    # WHEN parsing
    family = FamilyParser(lines).families[FAMILY]

    # THEN all four are in the family
    assert set(family.individuals) == {"proband", "mother", "father", "daughter"}


def test_comments_and_empty_lines_are_skipped() -> None:
    # GIVEN a ped file with comments and empty or blank lines
    lines = ["# comment\n", "\n", "   \n"] + TRIO_LINES

    # WHEN parsing
    parser = FamilyParser(lines)

    # THEN only the individuals are found
    assert len(parser.individuals) == 3


def test_whitespace_separated_lines_are_accepted() -> None:
    # GIVEN a ped line separated by spaces instead of tabs
    lines = [f"{FAMILY} proband father mother 1 2\n", f"{FAMILY} mother 0 0 2 1\n"]
    lines.append(f"{FAMILY} father 0 0 1 1\n")

    # WHEN parsing
    individuals = FamilyParser(lines).families[FAMILY].individuals

    # THEN the line is split on whitespace
    assert individuals["proband"].mother == "mother"


@pytest.mark.parametrize("sex, expected", [("1", 1), ("2", 2), ("0", 0), ("-9", 0), ("other", 0)])
def test_sex_outside_1_and_2_is_unknown(sex: str, expected: int) -> None:
    # GIVEN a line with a sex code
    # WHEN parsing
    individual = FamilyParser([f"fam\tind\t0\t0\t{sex}\t1\n"]).individuals["ind"]

    # THEN only 1 and 2 are kept, anything else is unknown (0)
    assert individual.sex == expected


@pytest.mark.parametrize("phenotype, expected", [("1", 1), ("2", 2), ("0", 0), ("-9", 0), (".", 0)])
def test_phenotype_outside_1_and_2_is_missing(phenotype: str, expected: int) -> None:
    # GIVEN a line with a phenotype
    # WHEN parsing
    individual = FamilyParser([f"fam\tind\t0\t0\t1\t{phenotype}\n"]).individuals["ind"]

    # THEN only 1 and 2 are kept, anything else is missing (0)
    assert individual.phenotype == expected
    assert individual.affected == (expected == 2)
    assert individual.healthy == (expected == 1)


def test_dot_parents_are_unknown_parents() -> None:
    # GIVEN an individual with '.' for both parents
    # WHEN parsing
    individual = FamilyParser(["fam\tind\t.\t.\t1\t1\n"]).individuals["ind"]

    # THEN they are stored as '0', the missing-parent value
    assert (individual.father, individual.mother) == ("0", "0")
    assert not individual.has_parents


def test_parse_from_file_handle(tmp_path: Path) -> None:
    # GIVEN a ped file on disk
    ped_file = tmp_path / "trio.ped"
    ped_file.write_text("".join(TRIO_LINES))

    # WHEN parsing the open file handle
    with open(ped_file) as handle:
        parser = FamilyParser(handle)

    # THEN the individuals are found
    assert set(parser.individuals) == {"proband", "mother", "father"}


def test_fam_type_is_parsed_as_ped() -> None:
    # WHEN parsing with family type 'fam'
    parser = FamilyParser(TRIO_LINES, "fam")

    # THEN it is read like a ped file
    assert set(parser.individuals) == {"proband", "mother", "father"}


def test_unknown_family_type_raises() -> None:
    # WHEN parsing with a family type that does not exist
    # THEN an error is raised instead of silently returning nothing
    with pytest.raises(ValueError, match="family type"):
        FamilyParser(TRIO_LINES, "madeline")


def test_empty_input_gives_no_families() -> None:
    # WHEN parsing nothing
    parser = FamilyParser([])

    # THEN no families are found
    assert parser.families == {}
    assert parser.individuals == {}


@pytest.mark.parametrize(
    "line",
    [f"{FAMILY}\tproband\tfather\tmother\t1\n", f"{FAMILY}\tproband\tfather\tmother\t0\t1\t2\n"],
)
def test_wrong_number_of_columns_raises(line: str) -> None:
    # GIVEN a ped line with 5 or 7 columns
    # WHEN parsing
    # THEN an error about the line format is raised
    with pytest.raises(WrongLineFormat):
        FamilyParser(TRIO_LINES[:2] + [line])


def test_missing_father_raises() -> None:
    # GIVEN a child whose father is not in the file
    lines = [f"{FAMILY}\tproband\tfather\tmother\t1\t2\n", f"{FAMILY}\tmother\t0\t0\t2\t1\n"]

    # WHEN parsing
    # THEN a pedigree error names the family and the missing parent
    with pytest.raises(PedigreeError, match="not in family") as error:
        FamilyParser(lines)
    assert error.value.family_id == FAMILY
    assert error.value.individual_id == "father"


def test_missing_mother_raises() -> None:
    # GIVEN a child whose mother is not in the file
    lines = [f"{FAMILY}\tproband\tfather\tmother\t1\t2\n", f"{FAMILY}\tfather\t0\t0\t1\t1\n"]

    # WHEN parsing
    # THEN a pedigree error is raised
    with pytest.raises(PedigreeError, match="not in family"):
        FamilyParser(lines)


def test_parents_with_wrong_sex_raise() -> None:
    # GIVEN a trio where the mother is male and the father female
    lines = [
        f"{FAMILY}\tproband\tfather\tmother\t1\t2\n",
        f"{FAMILY}\tmother\t0\t0\t1\t1\n",
        f"{FAMILY}\tfather\t0\t0\t2\t1\n",
    ]

    # WHEN parsing
    # THEN a pedigree error is raised
    with pytest.raises(PedigreeError):
        FamilyParser(lines)


def test_parent_in_another_family_is_not_found() -> None:
    # GIVEN a child whose parents are listed under a different family id
    lines = [
        f"{FAMILY}\tproband\tfather\tmother\t1\t2\n",
        "other\tmother\t0\t0\t2\t1\n",
        "other\tfather\t0\t0\t1\t1\n",
    ]

    # WHEN parsing
    # THEN the parents are not found in the child's family
    with pytest.raises(PedigreeError, match="not in family"):
        FamilyParser(lines)


# Duplicate individual ids


def test_duplicate_individual_id_in_a_family_raises() -> None:
    # GIVEN a family where the id "ind" is on two lines with different sex and phenotype
    lines = [f"{FAMILY}\tind\t0\t0\t1\t2\n", f"{FAMILY}\tind\t0\t0\t2\t1\n"]

    # WHEN parsing
    # THEN an error names the family and individual, instead of the last line replacing the first
    with pytest.raises(PedigreeError, match="more than once") as error:
        FamilyParser(lines)
    assert error.value.family_id == FAMILY
    assert error.value.individual_id == "ind"


def test_identical_duplicate_line_also_raises() -> None:
    # GIVEN the exact same line twice
    line = f"{FAMILY}\tind\t0\t0\t1\t2\n"

    # THEN it is still reported, since it is most likely a copy-paste error
    with pytest.raises(PedigreeError, match="more than once"):
        FamilyParser([line, line])


def test_same_individual_id_in_two_families_raises() -> None:
    # GIVEN two families that both have an individual "ind"
    lines = ["fam1\tind\t0\t0\t1\t2\n", "fam2\tind\t0\t0\t2\t1\n"]

    # WHEN parsing
    # THEN an error explains that genmod needs ids that are unique across the file,
    # since they are matched against the sample names of the VCF
    with pytest.raises(PedigreeError, match="unique across") as error:
        FamilyParser(lines)
    assert error.value.family_id == "fam2"
    assert error.value.individual_id == "ind"


def test_same_individual_id_in_two_families_raises_in_alternative_format() -> None:
    # GIVEN an alternative format file with the same individual in two families
    lines = [ALT_HEADER, "fam1\tind\t0\t0\t1\t2\tYes\tx\n", "fam2\tind\t0\t0\t2\t1\tNo\ty\n"]

    # THEN it is rejected as for ped files
    with pytest.raises(PedigreeError, match="unique across"):
        FamilyParser(lines, "alt")


# Alternative format: a header line and extra columns

ALT_HEADER = "#FamilyID\tSampleID\tFather\tMother\tSex\tPhenotype\tProband\tNotes\n"
ALT_LINES = [
    ALT_HEADER,
    f"{FAMILY}\tproband\tfather\tmother\t1\t2\tYes\tfirst\n",
    f"{FAMILY}\tmother\t0\t0\t2\t1\tNo\tsecond\n",
    f"{FAMILY}\tfather\t0\t0\t1\t1\tNo\tthird\n",
]


@pytest.mark.parametrize("family_type", ["alt", "cmms", "mip"])
def test_alternative_formats_are_parsed(family_type: str) -> None:
    # GIVEN a ped file with a header and extra columns
    # WHEN parsing it as an alternative type
    family = FamilyParser(ALT_LINES, family_type).families[FAMILY]

    # THEN the first six columns are read as in a ped file
    proband = family.individuals["proband"]
    assert (proband.father, proband.mother, proband.sex, proband.phenotype) == (
        "father",
        "mother",
        1,
        2,
    )
    # AND the extra columns are kept under their header names
    assert proband.extra_info == {"Proband": "Yes", "Notes": "first"}
    assert family.individuals["father"].extra_info["Notes"] == "third"


def test_alternative_format_requires_a_header() -> None:
    # GIVEN data lines without a header line
    lines = TRIO_LINES[2:]

    # WHEN parsing as the alternative format
    # THEN an error is raised asking for a header
    with pytest.raises(WrongLineFormat, match="header"):
        FamilyParser(lines, "alt")


@pytest.mark.parametrize(
    "line",
    [
        f"{FAMILY}\tproband\tfather\tmother\t1\t2\tYes\n",
        f"{FAMILY}\tproband\tfather\tmother\t1\t2\tYes\tfirst\textra\n",
    ],
)
def test_alternative_line_must_match_header_length(line: str) -> None:
    # GIVEN a data line with fewer or more columns than the header
    # WHEN parsing
    # THEN an error is raised
    with pytest.raises(WrongLineFormat):
        FamilyParser([ALT_HEADER, line], "alt")


def test_alternative_format_checks_parents() -> None:
    # GIVEN an alternative file where the father is missing
    lines = [ALT_HEADER, ALT_LINES[1], ALT_LINES[2]]

    # WHEN parsing
    # THEN the pedigree is checked like for ped files
    with pytest.raises(PedigreeError):
        FamilyParser(lines, "alt")


def test_alternative_header_and_lines_can_be_whitespace_separated() -> None:
    # GIVEN an alternative format file where the header and the data lines are separated by spaces
    lines = [
        "#FamilyID SampleID Father Mother Sex Phenotype Proband\n",
        f"{FAMILY} proband father mother 1 2 Yes\n",
        f"{FAMILY} mother 0 0 2 1 No\n",
        f"{FAMILY} father 0 0 1 1 No\n",
    ]

    # WHEN parsing as the alternative format
    family = FamilyParser(lines, "alt").families[FAMILY]

    # THEN the extra column is read under its header name
    assert family.individuals["proband"].extra_info == {"Proband": "Yes"}


def test_alternative_line_keeps_empty_trailing_column() -> None:
    # GIVEN an alternative format file where the last column is empty for every individual
    lines = [
        ALT_HEADER,
        f"{FAMILY}\tproband\tfather\tmother\t1\t2\tYes\t\n",
        f"{FAMILY}\tmother\t0\t0\t2\t1\tNo\t\n",
        f"{FAMILY}\tfather\t0\t0\t1\t1\tNo\t\n",
    ]

    # WHEN parsing as the alternative format
    family = FamilyParser(lines, "alt").families[FAMILY]

    # THEN the empty column is kept as an empty value
    assert family.individuals["proband"].extra_info == {"Proband": "Yes", "Notes": ""}


def test_ped_line_with_trailing_tab_is_still_accepted() -> None:
    # GIVEN a ped line with a trailing tab after its six columns
    lines = [f"{FAMILY}\tproband\t0\t0\t1\t2\t\n"]

    # WHEN parsing as a ped file
    # THEN the trailing tab is ignored and the individual is read
    family = FamilyParser(lines, "ped").families[FAMILY]
    assert family.individuals["proband"].extra_info == {}
