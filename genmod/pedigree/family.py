"""A family: a set of related individuals."""

from __future__ import annotations

from dataclasses import dataclass, field

from genmod.pedigree.exceptions import PedigreeError
from genmod.pedigree.individual import NO_PARENT, Individual

MALE = 1
FEMALE = 2


@dataclass
class Family:
    """The individuals of one family, keyed on individual id.

    Attributes:
        family_id: Id of the family
        individuals: The members of the family
    """

    family_id: str
    individuals: dict[str, Individual] = field(default_factory=dict)

    def add_individual(self, individual: Individual) -> None:
        """Add an individual to the family.

        Args:
            individual: The individual to add

        Raises:
            PedigreeError: If the individual belongs to another family or is already in
                this family
        """
        if individual.family != self.family_id:
            raise PedigreeError(
                self.family_id,
                individual.individual_id,
                "Family id of individual is not the same as family id for Family object",
            )
        if individual.individual_id in self.individuals:
            raise PedigreeError(
                self.family_id,
                individual.individual_id,
                f"Individual {individual.individual_id} is in family {self.family_id} "
                "more than once.",
            )
        self.individuals[individual.individual_id] = individual

    def check(self) -> None:
        """Check that the family members do not contradict each other.

        Every known parent must be in the family, the father must be male and the
        mother female.

        Raises:
            PedigreeError: If a parent is missing from the family or has the wrong sex
        """
        for individual in self.individuals.values():
            self._check_parent(individual.father, expected_sex=MALE, role="Father")
            self._check_parent(individual.mother, expected_sex=FEMALE, role="Mother")

    def _check_parent(self, parent_id: str, expected_sex: int, role: str) -> None:
        if parent_id == NO_PARENT:
            return
        if parent_id not in self.individuals:
            raise PedigreeError(self.family_id, parent_id, "Parent is not in family.")
        if self.individuals[parent_id].sex != expected_sex:
            label = "male" if expected_sex == MALE else "female"
            raise PedigreeError(self.family_id, parent_id, f"{role} is not specified as {label}.")
