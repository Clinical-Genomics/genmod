#!/usr/bin/env python
# encoding: utf-8
"""
print_headers.py

Print the lines in a header object to a file or stdout.

Created by Måns Magnusson on 2015-01-22.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import annotations

from typing import TextIO

from genmod.vcf_tools.header_parser import HeaderParser


def print_headers(head: HeaderParser, outfile: TextIO | None = None, silent: bool = False) -> None:
    """Print the vcf headers.

    If a result file is provided headers will be printed here, otherwise
    they are printed to stdout.

    Args:
        head: A vcf header object
        outfile: A file handle to write to
        silent: If nothing should be printed to stdout. Has no effect when writing to a file.
    """
    for header_line in head.print_header():
        if outfile is not None:
            outfile.write(header_line + "\n")
        elif not silent:
            print(header_line)
