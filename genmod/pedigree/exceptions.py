"""Errors raised when a pedigree (ped) file is malformed or inconsistent."""

from __future__ import annotations


class PedigreeError(Exception):
    """The individuals of a family contradict each other, e.g. a parent is missing."""

    def __init__(self, family_id: str, individual_id: str, message: str = "") -> None:
        """Create the error.

        Args:
            family_id: The family where the problem was found
            individual_id: The individual the problem is about
            message: A description of the problem
        """
        super().__init__(message)
        self.family_id = family_id
        self.individual_id = individual_id
        self.message = message

    def __str__(self) -> str:
        return self.message


class WrongLineFormat(Exception):
    """A line in a ped file has the wrong number of columns, or the header is missing."""

    def __init__(self, message: str = "", ped_line: str = "") -> None:
        """Create the error.

        Args:
            message: A description of the problem
            ped_line: The line that has the wrong format
        """
        super().__init__(message)
        self.message = message
        self.ped_line = ped_line

    def __str__(self) -> str:
        return self.message
