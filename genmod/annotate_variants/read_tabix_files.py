from __future__ import annotations

import logging
import os
from typing import Optional

import pysam

logger = logging.getLogger(__name__)


class TabixError(Exception):
    """Raised when a tabix file is unusable (not bgzipped, missing index, unreadable)"""


def get_tabixhandle(path: str) -> pysam.TabixFile:
    """Open a bgzipped, tabix indexed file.

    Args:
        path: Path to a ``.gz`` file with a ``.tbi`` index next to it.

    Returns:
        A handle to the file.

    Raises:
        TabixError: If the file is not ``.gz``, has no index or cannot be opened.
    """
    if not path.endswith(".gz"):
        raise TabixError(f"File {path} does not end with '.gz'")
    index_file = path + ".tbi"
    if not os.path.isfile(index_file):
        raise TabixError(f"No index could be found for {path}")

    try:
        return pysam.TabixFile(path)
    except (OSError, ValueError) as err:
        raise TabixError(f"Could not open {path}: {err}")


def get_contig_name(tabix_reader: pysam.TabixFile, chrom: str) -> Optional[str]:
    """Find the name a chromosome has in a tabix file, with or without a ``chr`` prefix.

    Args:
        tabix_reader: A tabix handle.
        chrom: The chromosome, e.g. ``1`` or ``chr1``.

    Returns:
        The matching contig name in the file, or None if the file has no such contig.
    """
    if chrom.startswith("chr"):
        alternative = chrom.removeprefix("chr")
    else:
        alternative = "chr" + chrom
    for contig in (chrom, alternative):
        if contig in tabix_reader.contigs:
            return contig
    return None


def get_tabix_records(
    tabix_reader: pysam.TabixFile, chrom: str, start: str | int
) -> list[tuple[str, ...]]:
    """Get the records overlapping a position, whether or not the file uses ``chr`` prefixes.

    Args:
        tabix_reader: A tabix handle.
        chrom: The chromosome of the position, with or without ``chr`` prefix.
        start: The 1-based position.

    Returns:
        The overlapping records as tuples of columns. Empty if there are none or the
        chromosome is not in the file.
    """
    logger.debug(f"Looking for record in {tabix_reader}")
    logger.debug(f"Looking for records with chr:{chrom}, pos:{start}")
    tabix_key = int(start)
    contig = get_contig_name(tabix_reader, str(chrom))
    if contig is None:
        logger.info(f"Chromosome {chrom} does not seem to exist in {tabix_reader}")
        return []

    return list(tabix_reader.fetch(contig, tabix_key - 1, tabix_key, parser=pysam.asTuple()))


def get_frequencies(
    tabix_reader: pysam.TabixFile, chrom: str, start: str | int, alt: str
) -> dict[str, str]:
    """Return the ``AF`` and ``MAX_AF`` frequencies of an alternative from an indexed VCF.

    Args:
        tabix_reader: A tabix handle.
        chrom: The chromosome of the position.
        start: The 1-based position.
        alt: The alternative sequence.

    Returns:
        The frequencies found, keyed by INFO key.
    """
    records = get_tabix_records(tabix_reader, chrom, start)

    frequencies = {}
    for record in records:
        row = "\t".join(record)
        logger.debug(f"Found record: {row}")
        # We can get multiple rows so need to check each one
        # We also need to check each one of the alternatives per row
        for i, alternative in enumerate(record[4].split(",")):
            if alternative == alt:
                for info in record[7].split(";"):
                    key_value = info.split("=")
                    if key_value[0] in ["AF", "MAX_AF"]:
                        freqs = key_value[-1].split(",")
                        frequencies[key_value[0]] = freqs[i]

    return frequencies


def get_spidex_score(
    tabix_reader: pysam.TabixFile, chrom: str, start: str | int, alt: str
) -> Optional[float]:
    """Return the SPIDEX z score of an alternative from a SPIDEX file.

    Args:
        tabix_reader: A tabix handle.
        chrom: The chromosome of the position.
        start: The 1-based position.
        alt: The alternative sequence.

    Returns:
        The z score, or None if there is no match.
    """
    records = get_tabix_records(tabix_reader, chrom, start)
    spidex_score = None

    for record in records:
        if record[3] == alt:
            spidex_score = float(record[5])

    logger.debug(f"Found spidex score: {spidex_score}")

    return spidex_score


def get_cosmic(tabix_reader: pysam.TabixFile, chrom: str, start: str | int, alt: str) -> bool:
    """Return whether an alternative at a position is in COSMIC.

    Args:
        tabix_reader: A tabix handle.
        chrom: The chromosome of the position.
        start: The 1-based position.
        alt: The alternative sequence.

    Returns:
        True if a COSMIC record matches the alternative.
    """
    records = get_tabix_records(tabix_reader, chrom, start)

    return any(record[4] == alt for record in records)


def get_cadd_scores(
    tabix_reader: pysam.TabixFile, chrom: str, start: str | int, alt: str
) -> dict[str, str]:
    """Return the CADD scores of an alternative from a CADD file.

    Args:
        tabix_reader: A tabix handle.
        chrom: The chromosome of the position.
        start: The 1-based position.
        alt: The alternative sequence.

    Returns:
        ``cadd_raw`` and ``cadd_phred`` as strings. Empty if there is no match.
    """
    cadd_scores = {}
    records = get_tabix_records(tabix_reader, chrom, start)
    # CADD values are only for snps:
    for record in records:
        if record[3] == alt:
            # We need to send both cadd values
            cadd_scores["cadd_raw"] = record[-2]
            cadd_scores["cadd_phred"] = record[-1]

    return cadd_scores
