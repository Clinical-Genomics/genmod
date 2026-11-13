from __future__ import annotations

import gzip
import sys
from multiprocessing import cpu_count
from typing import TextIO

import click

variant_file = click.argument("variant_file", type=click.Path(), metavar="<vcf_file> or -")

outfile = click.option(
    "-o",
    "--outfile",
    type=click.File("w"),
    help="Specify the path to a file where results should be stored.",
)

silent = click.option("-s", "--silent", is_flag=True, help="Do not print the variants.")

processes = click.option(
    "-p",
    "--processes",
    default=min(4, cpu_count()),
    help="Define how many processes that should be use for annotation.",
)

temp_dir = click.option("--temp_dir", type=click.Path(exists=True), help="Path to tempdir")

family_file = click.option("-f", "--family_file", type=click.File("r"), metavar="<ped_file>")

family_type = click.option(
    "-t",
    "--family_type",
    type=click.Choice(["ped", "alt", "cmms", "mip"]),
    default="ped",
    help="If the analysis use one of the known setups, please specify which one.",
)


def get_file_handle(path: str) -> TextIO:
    """Open a file for reading as text.

    Args:
        path: Path to the file. ``-`` means standard input and a path ending in ``.gz`` is
            read as a gzipped file.

    Returns:
        A text file handle. Bytes that are not valid utf-8 in gzipped files are replaced.

    Raises:
        FileNotFoundError: If the file does not exist
    """
    if path == "-":
        return sys.stdin

    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")

    return open(path, "r")
