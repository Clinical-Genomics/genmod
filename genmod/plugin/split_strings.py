"""Vendored from extract_vcf (https://github.com/moonso/extract_vcf, MIT License, Måns Magnusson)."""

from __future__ import annotations

import logging
from collections.abc import Sequence

logger = logging.getLogger(__name__)


def split_strings(string: str, separators: Sequence[str]) -> list[str]:
    """Split a string on several separators, applied in order, into a flat list.

    The string is first split on the first separator, each piece is then split on the
    second separator and so on. The pieces from the last split are returned in order.
    A plain string is accepted as ``separators`` and is treated as a sequence of
    single-character separators. With no separators the string is returned unsplit.

    Args:
        string: The string to split, e.g. ``"a:1|2,b:2"``.
        separators: The separators to split on, in order, e.g. ``[",", ":", "|"]``.

    Returns:
        The flat list of split values, e.g. ``["a", "1", "2", "b", "2"]``.
    """
    logger = logging.getLogger("genmod.plugin.split_strings")
    logger.debug("splitting string '{0}' with separators {1}".format(string, separators))
    results = []

    def recursion(recursive_string: str, separators: Sequence[str], i: int = 1) -> None:
        """Split a string on separator ``i`` and recurse on the pieces with the next one.

        The values from the last separator are appended to the enclosing ``results`` list.

        Args:
            recursive_string: The string (or piece of the original string) to split.
            separators: The separators to split on, in order.
            i: One-based index of the separator to split on at this level.
        """
        if i == len(separators):
            for value in recursive_string.split(separators[i - 1]):
                logger.debug("Adding {0} to results".format(value))
                results.append(value)
        else:
            for value in recursive_string.split(separators[i - 1]):
                recursion(value, separators, i + 1)

    if len(separators) > 0:
        recursion(string, separators)
    else:
        results = [string]

    return results
