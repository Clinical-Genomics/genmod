from __future__ import print_function

import logging
from collections import OrderedDict
from datetime import datetime

from genmod.vcf_tools import get_info_dict, get_variant_dict, get_variant_id, get_vep_dict

from .get_features import get_annotation

logger = logging.getLogger(__name__)


def get_batches(
    variants, batch_queue, header, vep=False, results_queue=None, annotation_keyword="Annotation"
):
    """
    Create variant batches based on their annotation and put them into the
    batch queue.

    Variant batches are are dictionaries with variant_id as key and
    variant_dict as value.

    get_batches will then use the annotation to search for sequences of variants
    with overlapping annotations. These are collected into one batch and gets
    put into a queue.

    Variants that are in between features will be in their own batch.

    Arguments:
         variants (Iterator):An iterator that returns variant dictionaries
         batch_queue (Queue): A queue where the batches will be putted
         header (HeaderParser): A HeaderParser object
         vep (bool): If variant is annotated with vep
         compound_mode (bool): If only compounds should be used
         results_queue (Queue): A queue where variants can be put for printing

    Returns:
         Does not return but put the results in a queue
    """

    logger.debug("Set beginning to True")
    beginning = True
    logger.debug("Create first empty batch")
    # A batch is a ordered dictionary with variants
    batch = OrderedDict()
    new_chrom = None
    current_chrom = None
    current_features = set()
    chromosomes = []

    start_parsing_time = datetime.now()
    start_chrom_time = start_parsing_time
    start_twenty_time = start_parsing_time

    nr_of_variants = 0
    nr_of_batches = 0

    header_line = header.header
    vep_header = header.vep_columns
    logger.info("Start parsing the variants")

    for line in variants:
        if not line.startswith("#"):
            variant = get_variant_dict(line, header_line)
            variant_id = get_variant_id(variant)
            variant["variant_id"] = variant_id
            variant["info_dict"] = get_info_dict(variant["INFO"])

            if vep:
                variant["vep_info"] = get_vep_dict(
                    vep_string=variant["info_dict"]["CSQ"],
                    vep_header=vep_header,
                    allele=variant["ALT"].split(",")[0],
                )

            logger.debug("Checking variant {0}".format(variant_id))

            nr_of_variants += 1
            new_chrom = variant["CHROM"]
            if new_chrom.startswith("chr"):
                new_chrom = new_chrom[3:]

            logger.debug("Update new chrom to {0}".format(new_chrom))

            new_features = get_annotation(
                variant=variant, vep=vep, annotation_key=annotation_keyword
            )
            logger.debug("Adding {0} to variant {1}".format(", ".join(new_features), variant_id))

            variant["annotation"] = new_features

            if nr_of_variants % 20000 == 0:
                logger.info("{0} variants parsed".format(nr_of_variants))
                logger.info(
                    "Last 20.000 took {0} to parse.".format(str(datetime.now() - start_twenty_time))
                )
                start_twenty_time = datetime.now()

            if beginning:
                logger.debug("First variant.")
                current_features = new_features

                logger.debug("Adding %s to variant batch" % variant_id)
                batch[variant_id] = variant

                logger.debug("Updating current chrom to {0}".format(new_chrom))
                current_chrom = new_chrom

                chromosomes.append(current_chrom)
                logger.debug("Adding chr {0} to chromosomes".format(new_chrom))

                beginning = False
                logger.debug("Updating beginning to False")

            else:
                # If we should put the batch in the queue:
                logger.debug("Updating send to True")
                send = True

                # Check if the variant ovelapps any features
                if len(new_features) != 0:
                    # Check if the features overlap the previous variants features
                    if new_features.intersection(current_features):
                        logger.debug("Set send to False since variant features overlap")
                        send = False
                # If we are at a new chromosome we finish the current batch:
                if new_chrom != current_chrom:
                    if current_chrom not in chromosomes:
                        chromosomes.append(current_chrom)
                    logger.debug("Adding chr {0} to chromosomes".format(new_chrom))
                    # New chromosome means new batch
                    send = True
                    logger.info(
                        "Chromosome {0} parsed. Time to parse chromosome: {1}".format(
                            current_chrom, datetime.now() - start_chrom_time
                        )
                    )
                    start_chrom_time = datetime.now()
                    current_chrom = new_chrom

                if send:
                    # Put the job in the queue
                    if len(batch) > 0:
                        logger.debug("Adding batch in queue")
                        batch_queue.put(batch)
                        nr_of_batches += 1
                    # Reset the variables
                    current_features = new_features
                    logger.debug("Initializing empty batch")
                    batch = {}
                else:
                    current_features = current_features.union(new_features)

                # Add variant to batch
                batch[variant_id] = variant

    if current_chrom not in chromosomes:
        logger.debug("Adding chr {0} to chromosomes".format(current_chrom))
        chromosomes.append(current_chrom)

    logger.info(
        "Chromosome {0} parsed. Time to parse chromosome: {0}".format(
            current_chrom,
        )
    )

    if len(batch) > 0:
        nr_of_batches += 1
        batch_queue.put(batch)
        logger.debug("Adding batch to queue")

    logger.info(
        "Variants parsed. Time to parse variants: {0}".format(
            str(datetime.now() - start_parsing_time)
        )
    )

    logger.info("Number of variants in variant file: {0}".format(nr_of_variants))
    logger.info("Number of batches created: {0}".format(nr_of_batches))

    return chromosomes
