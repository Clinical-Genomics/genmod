#!/usr/bin/env python
# encoding: utf-8
"""
print_variants.py

Print the variants of a file to vcf file.

There are two modes, 'vcf' or 'modified'.
If 'vcf' we expect plain vcf variants and print them as they came in.
If 'modified' the first column has been used for sorting so we skip that one.

If a outfile is provided the variants will be printed to this one.

Created by Måns Magnusson on 2015-01-22.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Literal, TextIO


def print_variant_dict(
    variant: Mapping[str, str],
    header_line: Sequence[str],
    outfile: TextIO | None = None,
    silent: bool = False,
) -> None:
    """Print a variant dictionary.

    Prints a variant dictionary to a file or stdout, with the columns in the order of the
    header. Columns that the variant lacks are printed as ``.``.

    Args:
        variant: A variant dictionary
        header_line: The header columns
        outfile: A file handle to write to
        silent: If nothing should be printed to stdout
    """
    print_line = [variant.get(entry, ".") for entry in header_line]

    if outfile is not None:
        outfile.write("\t".join(print_line) + "\n")
    elif not silent:
        print("\t".join(print_line))


def print_variant(
    variant_line: str | None = None,
    variant_dict: Mapping[str, str] | None = None,
    header_line: Sequence[str] | None = None,
    priority: str | None = None,
    outfile: TextIO | None = None,
    mode: Literal["vcf", "modified"] = "vcf",
    silent: bool = False,
) -> None:
    """Print a variant line.

    If a result file is provided the variant will be appended to the file,
    otherwise it is printed to stdout.

    There are two modes, 'vcf' or 'modified'.
    If 'vcf' we expect plain vcf variants and print them as they came in.
    If 'modified' the first column has been used for sorting so we skip
    that one.

    Args:
        variant_line: A vcf formatted variant line
        variant_dict: A variant dictionary, used instead of ``variant_line`` if given
        header_line: The header columns, needed to print a ``variant_dict``
        priority: The priority for this variant, added as first column in ``vcf`` mode
        outfile: A file handle to write to
        mode: ``vcf`` or ``modified``
        silent: If nothing should be printed to stdout

    Raises:
        ValueError: If neither a variant line nor a variant dictionary is given
        IOError: If a variant dictionary is given without a header line
    """
    if variant_dict is not None:
        if not header_line:
            raise IOError("Print line needs a header_line when printing variant dict.")
        print_line = [variant_dict.get(entry, ".") for entry in header_line]
    elif variant_line is not None:
        print_line = variant_line.rstrip().split("\t")
    else:
        raise ValueError("print_variant needs a variant_line or variant_dict")

    if mode == "modified":
        print_line = print_line[1:]
    elif priority:
        print_line = [priority] + print_line

    print_string = "\t".join(print_line)

    if outfile is not None:
        outfile.write(print_string + "\n")
    elif not silent:
        print(print_string)


def print_variant_for_sorting(variant_line: str, priority: str, outfile: TextIO) -> None:
    """Print a variant with its priority as first column, for sorting.

    Args:
        variant_line: A vcf variant line
        priority: The priority for this variant
        outfile: A file handle to the temporary variant file
    """
    columns = variant_line.rstrip().split("\t")

    outfile.write("{0}\t{1}\n".format(priority, "\t".join(columns)))
