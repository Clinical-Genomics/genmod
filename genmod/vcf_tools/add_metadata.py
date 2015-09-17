#!/usr/bin/env python
# encoding: utf-8
"""
add_metadata.py

Add a new metadata string to a vcf header object.

Created by Måns Magnusson on 2015-01-22.
Copyright (c) 2015 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function

import logging

from datetime import datetime
from genmod import __version__

def add_metadata(head, metadata_type, annotation_id, annotation_number='.',
                entry_type=None, description=None, version=None,
                command_line_string=''):
    """
    Add genmod metadata to the vcf header.
    
    Args:
        head            : A vcf header object
        metadata_type(str): 'info' or 'version'
        annotate_models : Bool. If genetic models should be annotated
        vep             : Bool. If variants already annotated with vep the they
                         will not get any region annotation by genmod.
        cadd_annotation : Bool. If cadd scores are annoated.
        cadd_raw        : Bool. If raw cadd scores should be annotated
        thousand_g      : Bool. If thousand genome frequencies are annotated
        exac            : Bool. If exac frequencies are added
        command_line_string : A string with the command that envoked genmod
    
    """
    logger = logging.getLogger(__name__)
    
    if metadata_type == 'info':
        logger.debug("Updating INFO header with {0}".format(
            annotation_id
        ))
        head.add_info(
            annotation_id,
            annotation_number,
            entry_type,
            description
        )
    elif metadata_type == 'version':
        logger.debug("Updating version header with {0}".format(
            annotation_id
        ))
    return

def add_version_header(head, command_line_string = ""):
    """Add Version information to the header"""
    head.add_version_tracking(
                    'genmod',
                    __version__,
                    datetime.now().strftime("%Y-%m-%d %H:%M"),
                    command_line_string
                )
    return
    

def add_annotation_header(head):
    """
    Add the Annotation information to a vcf header
    """
    add_metadata(
        head,
        'info',
        'Annotation',
        annotation_number='.',
        entry_type='String',
        description='Annotates what feature(s) this variant belongs to.'
    )
    return

def add_exonic_header(head):
    """
    Add the Exonic information to a vcf header
    """
    add_metadata(
        head,
        'info',
        'Exonic',
        annotation_number='0',
        entry_type='Flag',
        description='Indicates if the variant is exonic.'
    )
    return
    
def add_genetic_models_header(head):
    """
    Add Genetic Models to vcf header
    """
    add_metadata(
        head,
        'info',
        'GeneticModels',
        annotation_number='.',
        entry_type='String',
        description="':'-separated list of genetic models for this variant."
    )
    return

def add_model_score_header(head):
    """
    Add Model Score to vcf header
    """
    add_metadata(
        head,
        'info',
        'ModelScore',
        annotation_number='1',
        entry_type='Integer',
        description="PHRED score for genotype models."
    )
    return

def add_compounds_header(head):
    """
    Add compounds to vcf header
    """
    add_metadata(
        head,
        'info',
        'Compounds',
        annotation_number='.',
        entry_type='String',
        description=("List of compound pairs for this variant."
        "The list is splitted on ',' family id is separated with compounds"
        "with ':'. Compounds are separated with '|'.")
    )
    return
