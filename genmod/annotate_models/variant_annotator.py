#!/usr/bin/env python
# encoding: utf-8
"""
variant_consumer.py

Consumes batches of variants and annotates them. Each batch is a dictionary
with variant_id:s as keys and dictionaries with variant information.
The variants will get different annotations depending on input

Created by Måns Magnusson on 2013-03-01.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

from __future__ import absolute_import, division, print_function

import logging
import queue
from multiprocessing import Process, queues
from typing import Any, Optional, Union

from ped_parser import Family

from genmod.vcf_tools import get_genotypes

from .fix_variant import make_print_version
from .genetic_models import check_genetic_models

# A variant batch maps variant ids to variant dictionaries. ``None`` tells a worker to stop.
VariantBatch = dict[str, dict[str, Any]]


class VariantAnnotator(Process):
    """
    Annotates variant in batches from the task queue and puts the result in
    the results queue.
    """

    def __init__(
        self,
        task_queue: queues.JoinableQueue,
        results_queue: Union[queues.Queue, queue.Queue],
        families: dict[str, Family],
        individuals: list[str],
        phased: bool = False,
        strict: bool = False,
        reduced_penetrance_genes: Optional[set[str]] = None,
    ) -> None:
        """
        Initialize the VariantAnnotator

        Consume variant batches from the task queue, annotate them with the
        genetic inheritance patterns that they follow and put them in the
        results queue.

        Arguments:
            task_queue: Variant batches to annotate, and one ``None`` per worker to stop
            results_queue: Where the annotated variants are put for printing, a
                multiprocessing queue or a ``Manager().Queue()`` proxy
            families: The families to annotate, by family id
            individuals: The ids of the individuals in the analysis
            phased: If the genotype calls are phased (PS tags)
            strict: If a model is only annotated when all individuals are called
            reduced_penetrance_genes: Genes with reduced penetrance, as the gene ids
                of the variant annotation
        """
        Process.__init__(self)
        self.logger = logging.getLogger(__name__)

        self.proc_name = self.name

        self.logger.info("Setting up variant_annotator: {0}".format(self.proc_name))

        self.logger.debug("Setting up task queue")
        self.task_queue = task_queue

        self.logger.debug("Setting up results queue")
        self.results_queue = results_queue

        # The families that should be annotated
        self.families = families
        self.logger.debug("Families found: {0}".format(self.families))
        self.individuals = individuals
        self.logger.debug("Individuals found: {0}".format(self.individuals))

        # Settings for the annotation
        self.phased = phased
        self.logger.debug("Setting phased to {0}".format(self.phased))
        self.strict = strict
        self.logger.debug("Setting strict to {0}".format(self.strict))
        self.reduced_penetrance: set[str] = reduced_penetrance_genes or set()

    def run(self) -> None:
        """Run the consuming"""
        self.logger.info("%s: Starting!" % self.proc_name)
        # Check if there are any batches in the queue
        while True:
            # A batch is a dictionary with varints on the form {variant_id:variant_dict}
            self.logger.debug("Getting task from task_queue")
            variant_batch: Optional[VariantBatch] = self.task_queue.get()

            if variant_batch is None:
                self.logger.info("No more batches")
                self.task_queue.task_done()
                self.logger.info("{0}: Exiting".format(self.proc_name))
                break

            # We are now going to check the genetic models for the variants in
            # the batch

            for variant_id in variant_batch:
                variant = variant_batch[variant_id]
                variant["genotypes"] = get_genotypes(variant, self.individuals)

                # Check if the variant is in a gene with reduced penetrance
                if variant.get("annotation", set()).intersection(self.reduced_penetrance):
                    self.logger.debug(
                        "Setting reduced_penetrance to True for variant: {0}".format(variant_id)
                    )

                    variant["reduced_penetrance"] = True

            if len(variant_batch) > 1:
                # We only need to check compound candidates if there is
                # more than one variant in the batch
                for variant_id in variant_batch:
                    self.logger.debug("Check compound candidates")
                    variant = variant_batch[variant_id]

                    variant["compound_candidate"] = False

                    if variant.get("annotation"):
                        variant["compound_candidate"] = True
                        self.logger.debug("Set compound_candidate to True")

            # Check the genetic models for all variants in the batch
            check_genetic_models(
                variant_batch=variant_batch,
                families=self.families,
                phased=self.phased,
                strict=self.strict,
            )

            # # Now we want to make versions of the variants that are ready for printing.
            for variant_id in variant_batch:
                variant = make_print_version(
                    variant=variant_batch[variant_id], families=self.families
                )
                self.logger.debug("Putting variant in results_queue")
                self.results_queue.put(variant)

            self.task_queue.task_done()

        return
