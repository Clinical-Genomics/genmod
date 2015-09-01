from genmod.utils import (get_annotation, check_exonic)

import pkg_resources
import pytest

def get_variant(chrom='1', pos='1', ref='A', alt='G', annotation=["ADK"], exonic=True):
    """
    Return a variant dictionary
    """
    variant_id = '_'.join([chrom, pos, ref, alt])
    variant = {
        "CHROM":chrom,
        "POS":pos,
        "INFO":"Annotation={0}".format(annotation),
        'info_dict':{
            "Annotation":[region for region in annotation],
        },
        "variant_id": variant_id
    }
    
    if exonic:
        variant['info_dict']["Exonic"] = []

    return variant

def get_vep_variant(chrom='1', pos='1', ref='A', alt='G', annotation="ADK"):
    """
    Return a variant dictionary
    """
    variant_id = '_'.join([chrom, pos, ref, alt])
    variant = {
        "CHROM":chrom,
        "POS":pos,
        "INFO":"Annotation={0}".format(annotation),
        'vep_info':{
            'A': [{
                "SYMBOL": annotation,
                "Consequence": 'transcript_ablation' 
            }]
        },
        "variant_id": variant_id
    }
    return variant


def test_empty():
    """Test if get_features behave as suspected"""
    variant = {
        'CHROM': '1',
        'POS': '12',
        'ALT': 'A'
    }
    
    assert get_annotation(variant, {}) == set()

def test_simple():
    """Test if get_annotation behave as suspected"""
    variant = get_variant()
    assert get_annotation(variant) == set(['ADK'])

def test_double_region():
    """Test if get_annotation behave as suspected"""
    variant = get_variant(annotation=["ADK", "DDD"])
    assert get_annotation(variant) == set(["ADK", "DDD"])

def test_get_vep_region():
    """docstring for test_get_vep_region"""
    variant = get_vep_variant()
    assert get_annotation(variant, vep=True) == set(["ADK"])

def test_get_exonic():
    """docstring for test_get_exonic"""
    variant = get_variant()
    assert check_exonic(variant) == True

def test_get_exonic_vep():
    """docstring for test_get_exonic"""
    variant = get_vep_variant()
    assert check_exonic(variant, vep=True) == True
    
