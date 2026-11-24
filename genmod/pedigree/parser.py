"""Parse pedigree files into families and individuals.

Two formats are read. Both are tab separated (whitespace separated is also accepted):

``ped`` and ``fam``: exactly six columns on each line::

    Family_ID  Individual_ID  Paternal_ID  Maternal_ID  Sex  Phenotype

``alt``, ``cmms`` and ``mip``: a header line starting with ``#`` that names every column,
and data lines with as many columns as the header. The first six columns are the same as
in a ped file, the other columns are kept in ``Individual.extra_info``.

Individual ids must be unique across the whole file, as they are matched against the VCF
sample names (a ped file allows the same id in different families, genmod does not).

In both formats ``0`` or ``.`` mean that a parent is unknown, sex is 1 (male) or 2 (female)
and phenotype is 1 (unaffected) or 2 (affected). Other sex and phenotype codes (such as
``0`` and ``-9``) are read as unknown and missing (0). Lines starting with ``#`` and empty
lines are ignored in ped files.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable

from genmod.pedigree.exceptions import PedigreeError, WrongLineFormat
from genmod.pedigree.family import Family
from genmod.pedigree.individual import NO_PARENT, Individual

logger = logging.getLogger(__name__)

# The first six columns of every supported format
FAMILY_HEADER = ["family_id", "sample_id", "father_id", "mother_id", "sex", "phenotype"]
PED_TYPES = ["ped", "fam"]
ALTERNATIVE_TYPES = ["alt", "cmms", "mip"]


class FamilyParser:
    """Parse a ped file into families.

    Attributes:
        families: The families, keyed on family id
        individuals: All individuals of all families, keyed on individual id
    """

    def __init__(self, family_info: Iterable[str], family_type: str = "ped") -> None:
        """Parse the lines of a pedigree file.

        Args:
            family_info: The lines of the file, e.g. an open file handle or a list of strings
            family_type: ``ped`` or ``fam`` for six column files, ``alt``, ``cmms`` or ``mip``
                for files with a header and extra columns

        Raises:
            ValueError: If the family type is not known
            WrongLineFormat: If a line has the wrong number of columns or the header is missing
            PedigreeError: If a parent is missing from the family or has the wrong sex, or an
                individual id is used more than once
        """
        self.families: dict[str, Family] = {}
        self.individuals: dict[str, Individual] = {}
        logger.debug(f"Family type: {family_type}")

        if family_type in PED_TYPES:
            self._parse_ped(family_info)
        elif family_type in ALTERNATIVE_TYPES:
            self._parse_alternative(family_info)
        else:
            raise ValueError(
                f"Unknown family type '{family_type}', use one of {PED_TYPES + ALTERNATIVE_TYPES}"
            )

        for family in self.families.values():
            family.check()

    def _parse_ped(self, lines: Iterable[str]) -> None:
        for line in lines:
            if line.startswith("#") or not line.strip():
                continue
            columns = self._split_line(line)
            if len(columns) != len(FAMILY_HEADER):
                raise WrongLineFormat(message="WRONG FORMATED PED LINE!", ped_line=line.rstrip())
            self._add_individual(columns)

    def _parse_alternative(self, lines: Iterable[str]) -> None:
        header: list[str] | None = None
        for line in lines:
            if line.startswith("#"):
                header = self._split_alternative_line(line[1:])
                logger.debug(f"Alternative header found: {line}")
            elif line.strip():
                if not header:
                    raise WrongLineFormat(
                        message="Alternative ped files must have a header! "
                        "Please add a header line.",
                        ped_line=line.rstrip(),
                    )
                columns = self._split_alternative_line(line)
                if len(columns) != len(header):
                    raise WrongLineFormat(
                        message=(
                            f"Number of columns ({len(columns)}) differs from the header "
                            f"({len(header)})"
                        ),
                        ped_line=line.rstrip(),
                    )
                extra_info = dict(zip(header[len(FAMILY_HEADER) :], columns[len(FAMILY_HEADER) :]))
                self._add_individual(columns[: len(FAMILY_HEADER)], extra_info)

    @staticmethod
    def _split_line(line: str) -> list[str]:
        """Split on tabs, or on any whitespace if the line is not tab separated."""
        columns = line.rstrip().split("\t")
        if len(columns) < len(FAMILY_HEADER):
            columns = line.rstrip().split()
        return columns

    @staticmethod
    def _split_alternative_line(line: str) -> list[str]:
        """Split a header or data line of an alternative file, keeping empty trailing columns.

        Only the line break is removed before splitting on tabs, so a row ending in an empty
        extra column keeps that column. Lines that are not tab separated are split on whitespace.
        """
        stripped = line.rstrip("\r\n")
        columns = stripped.split("\t")
        if len(columns) < len(FAMILY_HEADER):
            columns = stripped.split()
        return [column.strip() for column in columns]

    def _add_individual(self, columns: list[str], extra_info: dict[str, str] | None = None) -> None:
        family_id, sample_id, father_id, mother_id, sex, phenotype = columns
        individual = Individual(
            individual_id=sample_id,
            family=family_id,
            mother=NO_PARENT if mother_id == "." else mother_id,
            father=NO_PARENT if father_id == "." else father_id,
            sex=int(sex) if sex in ("1", "2") else 0,
            phenotype=int(phenotype) if phenotype in ("1", "2") else 0,
            extra_info=extra_info or {},
        )
        family = self.families.setdefault(family_id, Family(family_id))
        family.add_individual(individual)
        if sample_id in self.individuals:
            # Legal in a ped file, but genmod matches the ids against the sample names in the
            # VCF, which are unique, so the individual would be ambiguous.
            raise PedigreeError(
                family_id,
                sample_id,
                f"Individual {sample_id} is in both family {self.individuals[sample_id].family} "
                f"and family {family_id}. Individual ids must be unique across the file.",
            )
        self.individuals[sample_id] = individual
