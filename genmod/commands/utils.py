import gzip
import sys
from codecs import getreader, open
from multiprocessing import cpu_count

import click

from typing import Optional

from genmod.utils import SNPEFF, VEP, ConsequenceSource

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


vep = click.option(
    "--vep", is_flag=True, help="If variants are annotated with the Variant Effect Predictor."
)

snpeff = click.option(
    "--snpeff", is_flag=True, help="If variants are annotated with SnpEff (the ANN field)."
)


def get_consequence_source(vep: bool, snpeff: bool) -> Optional[ConsequenceSource]:
    """Return where to take the genes of a variant from, given the --vep and --snpeff flags

    Args:
        vep: If variants are annotated with VEP (CSQ)
        snpeff: If variants are annotated with SnpEff (ANN)

    Returns:
        VEP, SNPEFF or None if neither flag is set

    Raises:
        click.UsageError: If both flags are set
    """
    if vep and snpeff:
        raise click.UsageError("--vep and --snpeff can not be used together")
    if vep:
        return VEP
    if snpeff:
        return SNPEFF
    return None


def get_file_handle(path):
    """Get a file handle"""
    if path == "-":
        if sys.version_info < (3, 0):
            sys.stdin = getreader("utf-8")(sys.stdin)

        file_handle = sys.stdin

    elif path.endswith(".gz"):
        file_handle = getreader("utf-8")(gzip.open(path, "r"), errors="replace")

    else:
        file_handle = open(path, "r")

    return file_handle
