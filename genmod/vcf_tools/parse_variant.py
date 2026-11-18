"""Parse a variant line in different ways"""

import logging
from collections.abc import Sequence
from typing import Union

logging = logging.getLogger(__name__)


def get_variant_dict(variant_line, header_line):
    """Parse a variant line

    Split a variant line and map the fields on the header columns

    Args:
        variant_line (str): A vcf variant line
        header_line (list): A list with the header columns
    Returns:
        variant_dict (dict): A variant dictionary
    """
    # logger.debug("Building variant dict from {0}".format(variant_line))

    return dict(zip(header_line, variant_line.rstrip().split("\t")))


def get_info_dict(info_line):
    """Parse a info field of a variant

    Make a dictionary from the info field of a vcf variant.
    Keys are the info keys and values are the raw strings from the vcf

    Args:
        info_line (str): The info field of a vcf variant
    Returns:
        info_dict (dict): A INFO dictionary
    """

    variant_info = {}
    for raw_info in info_line.split(";"):
        splitted_info = raw_info.split("=")
        if len(splitted_info) == 2:
            variant_info[splitted_info[0]] = splitted_info[1]
        else:
            variant_info[splitted_info[0]] = []

    return variant_info


def get_variant_id(variant_dict):
    """Build a variant id

    Build CHROM_POS_REF_ALT and sanitize ALT by stripping '<>[]:'.
    For structural variants with symbolic ALT values (for example <DEL>),
    append _END{END} when END is available. If END is missing but SVLEN is
    present, append _SVLEN{SVLEN} instead.

    Args:
        variant_dict (dict): A variant dictionary

    Returns:
        variant_id (str)
    """
    chrom = variant_dict["CHROM"]
    pos = variant_dict["POS"]
    ref = variant_dict["REF"]
    alt_raw = variant_dict["ALT"]
    # There are several symbols in structural variant calls that make
    # things hard. We will strip those symbols
    bad_chars = "<>[]:"
    alt = "".join(c for c in alt_raw if c not in bad_chars)
    base_id = "_".join([chrom, pos, ref, alt])

    # For SVs with symbolic ALT values, we append END/SVLEN if present to reduce id collisions
    if "<" in alt_raw or ">" in alt_raw:
        info = variant_dict.get("info_dict", {})
        end = info.get("END")
        svlen = info.get("SVLEN")
        if end:
            return "{}_END{}".format(base_id, end)
        if svlen:
            return "{}_SVLEN{}".format(base_id, svlen)

    return base_id


def parse_annotations(
    annotation_string: Union[str, list, None], columns: Sequence[str]
) -> list[dict[str, str]]:
    """Parse a consequence annotation INFO value, VEP ``CSQ`` or SnpEff ``ANN``

    Both formats separate the annotations (one per transcript or feature) with ``,``
    and the fields of one annotation with ``|``, in the order of ``columns``.

    Args:
        annotation_string: The raw INFO value. Empty, ``None`` (key not in the variant) and
            ``[]`` (``get_info_dict`` value of a flag without value) give no annotations
        columns: The field names from the INFO header line

    Return:
        One dictionary per annotation, field name to value
    """
    if not annotation_string:
        return []
    return [
        dict(zip(columns, annotation.split("|"))) for annotation in annotation_string.split(",")
    ]
