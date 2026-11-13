#!/usr/bin/env python
# encoding: utf-8
"""
add_metadata.py

Add a new metadata string to a vcf header object.

Created by Måns Magnusson on 2015-01-22.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import annotations

import logging
from datetime import datetime
from typing import Literal

from genmod import __version__
from genmod.vcf_tools.header_parser import HeaderParser

logger = logging.getLogger(__name__)


def add_metadata(
    head: HeaderParser,
    metadata_type: Literal["info"],
    annotation_id: str,
    annotation_number: str = ".",
    entry_type: str | None = None,
    description: str | None = None,
) -> None:
    """Add genmod metadata to the vcf header.

    Args:
        head: The vcf header to add to
        metadata_type: The kind of header line to add, only ``info`` is supported
        annotation_id: The ID of the INFO key, e.g. ``GeneticModels``
        annotation_number: The Number of the INFO key: an integer, ``A``, ``R``, ``G`` or ``.``
        entry_type: The Type of the INFO key: Integer, Float, Flag, Character or String
        description: The Description of the INFO key

    Raises:
        ValueError: If the metadata type is not supported
    """
    if metadata_type != "info":
        raise ValueError(f"Unknown metadata type '{metadata_type}', only 'info' is supported")
    logger.debug(f"Updating INFO header with {annotation_id}")
    head.add_info(annotation_id, annotation_number, entry_type, description)


def add_version_header(head: HeaderParser, command_line_string: str = "") -> None:
    """Add Version information to the header.

    Args:
        head: The vcf header to add to
        command_line_string: The command that invoked genmod
    """
    head.add_version_tracking(
        "genmod", __version__, datetime.now().strftime("%Y-%m-%d %H:%M"), command_line_string
    )


def add_annotation_header(head: HeaderParser) -> None:
    """Add the Annotation information to a vcf header."""
    add_metadata(
        head,
        "info",
        "Annotation",
        annotation_number=".",
        entry_type="String",
        description="Annotates what feature(s) this variant belongs to.",
    )


def add_exonic_header(head: HeaderParser) -> None:
    """Add the Exonic information to a vcf header."""
    add_metadata(
        head,
        "info",
        "Exonic",
        annotation_number="0",
        entry_type="Flag",
        description="Indicates if the variant is exonic.",
    )


def add_genetic_models_header(head: HeaderParser) -> None:
    """Add Genetic Models to a vcf header."""
    add_metadata(
        head,
        "info",
        "GeneticModels",
        annotation_number=".",
        entry_type="String",
        description="':'-separated list of genetic models for this variant.",
    )


def add_model_score_header(head: HeaderParser) -> None:
    """Add Model Score to a vcf header."""
    add_metadata(
        head,
        "info",
        "ModelScore",
        annotation_number="1",
        entry_type="Integer",
        description="PHRED score for genotype models.",
    )


def add_compounds_header(head: HeaderParser) -> None:
    """Add compounds to a vcf header."""
    add_metadata(
        head,
        "info",
        "Compounds",
        annotation_number=".",
        entry_type="String",
        description=(
            "List of compound pairs for this variant."
            "The list is splitted on ',' family id is separated with compounds"
            "with ':'. Compounds are separated with '|'."
        ),
    )
