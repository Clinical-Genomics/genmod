"""Pedigree handling: parse ped files into families and individuals.

Based on the ``ped_parser`` package, trimmed to what genmod uses.
"""

from genmod.pedigree.exceptions import PedigreeError, WrongLineFormat
from genmod.pedigree.family import Family
from genmod.pedigree.individual import Individual
from genmod.pedigree.parser import FAMILY_HEADER, FamilyParser

__all__ = [
    "FAMILY_HEADER",
    "Family",
    "FamilyParser",
    "Individual",
    "PedigreeError",
    "WrongLineFormat",
]
