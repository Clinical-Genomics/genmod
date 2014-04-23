#!/usr/bin/env python
# encoding: utf-8
"""
test_x_dominant.py

Test the so that the genetic models behave as suspected.


Created by Måns Magnusson on 2013-03-07.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from ped_parser import family, individual
from genmod.models import genetic_models
from genmod.variants import genotype


class TestXDominantModel(object):
    """Test class for testing how the genetic models behave with a dominant variant"""

    def setup_class(self):
        """Setup a simple family with family id 1, sick son id 1,
         healthy father id 2, healthy mother id 3"""
        # Setup family with sick kid, healthy father and sick mother:
        self.dominant_family = family.Family(family_id = '1')
        #The son will only inherit x variants from the mother
        sick_son = individual.Individual(ind='1', family='1',mother='3', father='2', sex=1, phenotype=2)
        healthy_father = individual.Individual(ind='2', family='1',mother='0', father='0', sex=1, phenotype=1)
        sick_mother = individual.Individual(ind='3', family='1',mother='0', father='0', sex=2, phenotype=2)
        
        self.dominant_family.add_individual(healthy_father)
        self.dominant_family.add_individual(sick_son)
        self.dominant_family.add_individual(sick_mother)
        
        self.x_dominant_het_variant = {'CHROM':'X', 'POS':'5', 'ALT':'A', 'REF':'C', 'ID':'rs2230749', '1':'0/1', '2':'0/0', '3':'0/1'}
        
        self.x_dominant_dn_missing = {'CHROM':'X', 'POS':'10', 'ALT':'C', 'REF':'T', 'ID':'.', '1':'1/1', '2':'0/0', '3':'./.'}
        
        self.x_dominant_variant = {'CHROM':'X', 'POS':'13', 'ALT':'A', 'REF':'C', 'ID':'rs2230749', '1':'1/1', '2':'0/0', '3':'0/1'}

        self.not_x_dominant = {'CHROM':'X', 'POS':'15', 'ALT':'C', 'REF':'T', 'ID':'.', '1':'0/1', '2':'0/1', '3':'0/1'}

        
        #This batch simulates one genewith all above variants:
        batch = {'ABC':{'X_5_A_C':self.x_dominant_het_variant, 'X_10_C_T':self.x_dominant_dn_missing,
             'X_13_A_C':self.x_dominant_variant, 'X_15_C_T':self.not_x_dominant}}
        
        genetic_models.check_genetic_models(batch, self.dominant_family, verbose=True)

            
    def test_x_dominant_het_variant(self):
        """This variant should only follow the x dominant pattern."""
        assert self.x_dominant_het_variant['Inheritance_model']['XD']
        assert not self.x_dominant_het_variant['Inheritance_model']['XD_dn']
        # Since mother is sick she has to be 1/1 for the variant to be recessive:
        assert not self.x_dominant_het_variant['Inheritance_model']['XR']
        assert not self.x_dominant_het_variant['Inheritance_model']['XR_dn']
    
    def test_x_dominant_dn_missing(self):
        """This variant should follow all x-patterns."""
        assert self.x_dominant_dn_missing['Inheritance_model']['XD']
        assert self.x_dominant_dn_missing['Inheritance_model']['XD_dn']
        assert self.x_dominant_dn_missing['Inheritance_model']['XR']
        assert self.x_dominant_dn_missing['Inheritance_model']['XR_dn']
    
    def test_x_dominant_variant(self):
        """This variant should only follow the x dominant pattern."""
        assert self.x_dominant_variant['Inheritance_model']['XD']
        assert not self.x_dominant_variant['Inheritance_model']['XD_dn']
        assert not self.x_dominant_variant['Inheritance_model']['XR']
        assert not self.x_dominant_variant['Inheritance_model']['XR_dn']
        
    def test_not_x_dominant(self):
        """This variant should not follow any pattern."""
        assert not self.not_x_dominant['Inheritance_model']['XD']
        assert not self.not_x_dominant['Inheritance_model']['XD_dn']
        assert not self.not_x_dominant['Inheritance_model']['XR']
        assert not self.not_x_dominant['Inheritance_model']['XR_dn']
    


def main():
    pass


if __name__ == '__main__':
    main()

