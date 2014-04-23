#!/usr/bin/env python
# encoding: utf-8
"""
test_models_compound.py

Test the so that the genetic models behave as suspected.


Created by Måns Magnusson on 2013-03-07.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from ped_parser import family, individual
from genmod.models import genetic_models
from genmod.variants import genotype


class TestModelsCompound(object):
    """Test class for testing how the genetic models behave with combinations of compound variants."""

    def setup_class(self):
        """Setup a simple family with family id 1, sick son id 1,
         healthy father id 2, healthy mother id 3"""
        # Setup family with sick kid, sick father and healthy mother:
        self.recessive_family = family.Family(family_id = '1')
        sick_son = individual.Individual(ind='1', family='1',mother='3', father='2', sex=1, phenotype=2)
        healthy_father = individual.Individual(ind='2', family='1',mother='0', father='0', sex=1, phenotype=1)
        healthy_mother = individual.Individual(ind='3', family='1',mother='0', father='0', sex=2, phenotype=1)
        self.recessive_family.add_individual(healthy_father)
        self.recessive_family.add_individual(sick_son)
        self.recessive_family.add_individual(healthy_mother)
        
        #Setup two variants with only autosomal recessive pattern
        self.recessive_comp_variant_1 = {'CHROM':'1', 'POS':'5', 'ALT':'A', 'REF':'C', 'ID':'rs2230749', '1':'0/1', '2':'0/1', '3':'0/0'}
        self.recessive_comp_variant_2 = {'CHROM':'1', 'POS':'10', 'ALT':'C', 'REF':'T', 'ID':'.', '1':'0/1', '2':'0/0', '3':'0/1'}
        #This one should not be recessive since both parents will have a combination:
        self.not_recessive_comp =  {'CHROM':'1', 'POS':'7', 'ALT':'G', 'REF':'T', 'ID':'.', '1':'0/1', '2':'0/1', '3':'0/1'}

        self.recessive_comp_missing = {'CHROM':'1', 'POS':'15', 'ALT':'C', 'REF':'T', 'ID':'.', '1':'0/1', '2':'0/1', '3':'./.'}
        
        batch = {'ABC':{'1_5_A_C':self.recessive_comp_variant_1, '1_10_C_T':self.recessive_comp_variant_2, '1_7_G_T':self.not_recessive_comp},
                'BBC':{'1_10_C_T':self.recessive_comp_variant_2, '1_15_C_T':self.recessive_comp_missing}}
        
        genetic_models.check_genetic_models(batch, self.recessive_family)
    
    def test_recessive_comp_1(self):
        """The compound inheritnace pattern should be followed for this variant."""
        assert not self.recessive_comp_variant_1['Inheritance_model']['AR_hom']
        assert not self.recessive_comp_variant_1['Inheritance_model']['AR_hom_dn']
        assert not self.recessive_comp_variant_1['Inheritance_model']['AD']
        assert not self.recessive_comp_variant_1['Inheritance_model']['AD_dn']
        assert self.recessive_comp_variant_1['Inheritance_model']['AR_comp']
    
    def test_recessive_comp_2(self):
        """The compound inheritnace pattern should be followed for this variant."""
        assert not self.recessive_comp_variant_2['Inheritance_model']['AR_hom']
        assert not self.recessive_comp_variant_2['Inheritance_model']['AR_hom_dn']
        assert not self.recessive_comp_variant_2['Inheritance_model']['AD']
        assert not self.recessive_comp_variant_2['Inheritance_model']['AD_dn']
        assert self.recessive_comp_variant_2['Inheritance_model']['AR_comp']
    
    def test_not_recessive_comp(self):
        """No inheritance pattern should be followed for this model."""
        assert not self.not_recessive_comp['Inheritance_model']['AR_hom']
        assert not self.not_recessive_comp['Inheritance_model']['AR_hom_dn']
        assert not self.not_recessive_comp['Inheritance_model']['AD']
        assert not self.not_recessive_comp['Inheritance_model']['AD_dn']
        assert not self.not_recessive_comp['Inheritance_model']['AR_comp']
    
    def test_recessive_comp_missing(self):
        """Compound pattern should be followed here."""
        assert not self.recessive_comp_missing['Inheritance_model']['AR_hom']
        assert not self.recessive_comp_missing['Inheritance_model']['AR_hom_dn']
        assert not self.recessive_comp_missing['Inheritance_model']['AD']
        assert not self.recessive_comp_missing['Inheritance_model']['AD_dn']
        assert self.recessive_comp_missing['Inheritance_model']['AR_comp']
    


def main():
    pass


if __name__ == '__main__':
    main()

