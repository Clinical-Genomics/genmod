"""An individual in a pedigree."""

from __future__ import annotations

from dataclasses import dataclass, field

# Value of mother and father when the parent is not known
NO_PARENT = "0"


@dataclass
class Individual:
    """A person in a family, as described by one line in a ped file.

    Attributes:
        individual_id: Id of the individual, unique within the family
        family: Id of the family the individual belongs to
        mother: Id of the mother, ``"0"`` if unknown
        father: Id of the father, ``"0"`` if unknown
        sex: 1 for male, 2 for female, 0 for unknown
        phenotype: 1 for unaffected, 2 for affected, 0 for missing
        extra_info: Values of columns after the sixth in alternative ped files,
            keyed on the column name in the header
    """

    individual_id: str
    family: str
    mother: str = NO_PARENT
    father: str = NO_PARENT
    sex: int = 0
    phenotype: int = 0
    extra_info: dict[str, str] = field(default_factory=dict)

    @property
    def has_parents(self) -> bool:
        """True if at least one parent is known."""
        return self.mother != NO_PARENT or self.father != NO_PARENT

    @property
    def has_both_parents(self) -> bool:
        """True if both parents are known."""
        return self.mother != NO_PARENT and self.father != NO_PARENT

    @property
    def affected(self) -> bool:
        """True if the phenotype is affected (2)."""
        return self.phenotype == 2

    @property
    def healthy(self) -> bool:
        """True if the phenotype is unaffected (1)."""
        return self.phenotype == 1
