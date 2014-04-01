#!/usr/bin/env python
# encoding: utf-8
"""
test_x_recessive.py

Test the so that the genetic models behave as suspected.


Created by Måns Magnusson on 2013-03-07.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from ped_parser import family, individual
from genmod.models import genetic_models
from genmod.variants import genotype


class TestXRecessiveModel(object):
    """Test class for testing how the genetic models behave with a x_recessive variant"""

    def setup_class(self):
        """Setup a simple family with family id 1, sick son id 1,
         healthy father id 2, healthy mother id 3"""
        # Setup family with sick kid, sick father and healthy mother:
        self.recessive_family = family.Family(family_id = '1')
        sick_daughter = individual.Individual(ind='1', family='1',mother='3', father='2', sex=2, phenotype=2)
        healthy_father = individual.Individual(ind='2', family='1',mother='0', father='0', sex=1, phenotype=1)
        healthy_mother = individual.Individual(ind='3', family='1',mother='0', father='0', sex=2, phenotype=1)
        
        self.recessive_family.add_individual(healthy_father)
        self.recessive_family.add_individual(sick_daughter)
        self.recessive_family.add_individual(healthy_mother)
        
        self.x_recessive_dn_variant = {'CHROM':'X', 'POS':'5', 'ALT':'A', 'REF':'C', 'ID':'rs2230749', '1':'0/1', '2':'0/0', '3':'0/0'}
        
        self.x_recessive_missing = {'CHROM':'X', 'POS':'10', 'ALT':'C', 'REF':'T', 'ID':'.', '1':'1/1', '2':'0/0', '3':'./.'}
        
        self.x_recessive_variant = {'CHROM':'X', 'POS':'13', 'ALT':'A', 'REF':'C', 'ID':'rs2230749', '1':'1/1', '2':'0/0', '3':'0/1'}

        self.not_x_recessive = {'CHROM':'X', 'POS':'15', 'ALT':'C', 'REF':'T', 'ID':'.', '1':'0/1', '2':'./.', '3':'1/1'}

        
        #This batch simulates one genewith all above variants:
        batch = {'ABC':{'X_5_A_C':self.x_recessive_dn_variant, 'X_10_C_T':self.x_recessive_missing,
             'X_13_A_C':self.x_recessive_variant, 'X_15_C_T':self.not_x_recessive}}
        
        genetic_models.check_genetic_models(batch, self.recessive_family, verbose=True)

            
    def test_x_recessive_dn_variant(self):
        """This variant should follow the x-dominant dn and x-recessive dn patterns."""
        assert not self.x_recessive_dn_variant['Inheritance_model']['XD']
        assert self.x_recessive_dn_variant['Inheritance_model']['XD_denovo']
        assert not self.x_recessive_dn_variant['Inheritance_model']['XR']
        assert self.x_recessive_dn_variant['Inheritance_model']['XR_denovo']
    
    def test_x_recessive_missing(self):
        """This variant should follow the x-dominant and x-recessive patterns."""
        assert self.x_recessive_missing['Inheritance_model']['XD']
        assert self.x_recessive_missing['Inheritance_model']['XD_denovo']
        assert self.x_recessive_missing['Inheritance_model']['XR']
        assert self.x_recessive_missing['Inheritance_model']['XR_denovo']
    
    def test_x_recessive_variant(self):
        """This variant should only follow the x-recessive pattern."""
        assert not self.x_recessive_variant['Inheritance_model']['XD']
        assert not self.x_recessive_variant['Inheritance_model']['XD_denovo']
        assert self.x_recessive_variant['Inheritance_model']['XR']
        assert not self.x_recessive_variant['Inheritance_model']['XR_denovo']
        
    def test_not_x_recessive(self):
        """This variant should not follow any pattern."""
        assert not self.not_x_recessive['Inheritance_model']['XD']
        assert not self.not_x_recessive['Inheritance_model']['XD_denovo']
        assert not self.not_x_recessive['Inheritance_model']['XR']
        assert not self.not_x_recessive['Inheritance_model']['XR_denovo']
    


def main():
    pass


if __name__ == '__main__':
    main()

