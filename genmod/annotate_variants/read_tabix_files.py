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
        raise TabixError("File {0} does not end with '.gz'".format(path))
    index_file = path + ".tbi"
    if not os.path.isfile(index_file):
        raise TabixError("No index could be found for {0}".format(path))

    try:
        return pysam.TabixFile(path)
    except (OSError, ValueError) as err:
        raise TabixError("Could not open {0}: {1}".format(path, err))


def get_tabix_records(
    tabix_reader: pysam.TabixFile, chrom: str, start: str | int
) -> list[tuple[str, ...]]:
    """Get the records overlapping a position, trying the contig with and without ``chr``.

    Args:
        tabix_reader: A tabix handle.
        chrom: The chromosome of the position.
        start: The 1-based position.

    Returns:
        The overlapping records as tuples of columns. Empty if there are none or the
        chromosome is not in the file.
    """
    logger.debug("Looking for record in %s" % (tabix_reader))
    logger.debug("Looking for records with chr:%s, pos:%s" % (chrom, start))
    tabix_key = int(start)
    chrom = str(chrom)
    # pysam raises ValueError for contigs that are missing from the index.
    # Records are materialised so that errors surface here and not in the caller.
    for contig in (chrom, "chr" + chrom):
        try:
            return list(
                tabix_reader.fetch(contig, tabix_key - 1, tabix_key, parser=pysam.asTuple())
            )
        except ValueError:
            continue
    logger.info("Chromosome {0} does not seem to exist in {1}".format(chrom, tabix_reader))

    return []


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
        logger.debug("Found record: %s" % "\t".join(record))
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

    logger.debug("Found spidex score: %s" % str(spidex_score))

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
    in_cosmic = False

    for record in records:
        if record[4] == alt:
            in_cosmic = True
    if in_cosmic:
        logger.debug("Variant was found in COSMIC")
    else:
        logger.debug("Variant was not found in COSMIC")

    return in_cosmic


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
