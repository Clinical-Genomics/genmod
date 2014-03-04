#!/usr/bin/env python
# encoding: utf-8
"""
test_models_phased_duo.py

Test the so that the genetic models behave as suspected.


Created by Måns Magnusson on 2014-03-04.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from ped_parser import family, individual
from genmod.utils import interval_tree
from genmod.models import genetic_models
from genmod.variants import genotype


class TestModelsCompound(object):
    """Test class for testing how the genetic models behave with a recessive variant"""

    def setup_class(self):
        """Setup a simple family with family id 1, sick son id 1,
         healthy father id 2, healthy mother id 3"""
        # Setup family with sick kid, sick father and healthy mother:
        self.duo_family = family.Family(family_id = '1')
        sick_son = individual.Individual(ind='1', family='1',mother='3', father='0', sex=1, phenotype=2)
        healthy_mother = individual.Individual(ind='3', family='1',mother='0', father='0', sex=2, phenotype=1)
        self.duo_family.add_individual(sick_son)
        self.duo_family.add_individual(healthy_mother)
        
        intervals = {ind_id:interval_tree.intervalTree([[1,100, '1']], 0, 1, 1, 100) for ind_id in self.duo_family.individuals}
        
        #Setup two variants with only autosomal recessive pattern
        self.recessive_comp_simple_1 = {'CHROM':'1', 'POS':'5', 'ALT':'A', 'REF':'C', 'ID':'rs2230749',
                                                 '1':'0|1', '3':'0|0'}
        
        self.recessive_comp_simple_2 = {'CHROM':'1', 'POS':'10', 'ALT':'C', 'REF':'T', 'ID':'.', 
                                                '1':'1|0', '3':'0|1'}

        #Setup two variants that parents have but on same allele
        self.recessive_comp_not_simple_1 = {'CHROM':'1', 'POS':'15', 'ALT':'A', 'REF':'C', 'ID':'.',
                                                    '1':'0|1', '3':'1|0'}
        
        self.recessive_comp_not_simple_2 = {'CHROM':'1', 'POS':'17', 'ALT':'C', 'REF':'T', 'ID':'.', 
                                                    '1':'1|0', '3':'1|0'}
        
        #This one should not be recessive since the sick have both variants on same allele:
        self.not_recessive_comp_1 =  {'CHROM':'1', 'POS':'27', 'ALT':'G', 'REF':'T', 'ID':'.', 
                                        '1':'0|1', '3':'0|0'}

        self.not_recessive_comp_2 = {'CHROM':'1', 'POS':'35', 'ALT':'C', 'REF':'T', 'ID':'.', 
                                        '1':'0|1', '3':'0|1'}
        
        #Setup two variants that parents have but on same allele
        self.recessive_comp_missing_1 = {'CHROM':'1', 'POS':'45', 'ALT':'A', 'REF':'C', 'ID':'.',
                                                    '1':'0|1', '3':'0|0'}
        
        self.recessive_comp_missing_2 = {'CHROM':'1', 'POS':'57', 'ALT':'C', 'REF':'T', 'ID':'.', 
                                                    '1':'1|0', '3':'./.'}
        
        batch = {'ABC':{'1_5_A_C':self.recessive_comp_simple_1, '1_10_C_T':self.recessive_comp_simple_2},
                'BBC':{'1_15_A_C':self.recessive_comp_not_simple_1, '1_17_C_T':self.recessive_comp_not_simple_2},
                'DBC':{'1_45_A_C':self.recessive_comp_missing_1, '1_57_C_T':self.recessive_comp_missing_2},
                'CBC':{'1_27_G_T':self.not_recessive_comp_1, '1_35_C_T':self.not_recessive_comp_2}
                }
        
        batch['intervals'] = intervals
        genetic_models.check_genetic_models(batch, self.duo_family, phased=True)
    
    def test_recessive_comp_simple(self):
        """Check if the genetic models are followed for the heterozygote variant"""
        assert not self.recessive_comp_simple_1['Inheritance_model']['AR_hom']
        assert not self.recessive_comp_simple_1['Inheritance_model']['AR_hom_denovo']
        assert self.recessive_comp_simple_1['Inheritance_model']['AD']
        assert self.recessive_comp_simple_1['Inheritance_model']['AD_denovo']
        assert not self.recessive_comp_simple_1['Inheritance_model']['X']
        assert not self.recessive_comp_simple_1['Inheritance_model']['X_dn']
        assert self.recessive_comp_simple_1['Inheritance_model']['AR_compound']
    
    def test_recessive_comp_not_simple(self):
        """docstring for test_recessive_comp_2"""
        assert not self.recessive_comp_not_simple_2['Inheritance_model']['AR_hom']
        assert not self.recessive_comp_not_simple_2['Inheritance_model']['AR_hom_denovo']
        assert not self.recessive_comp_not_simple_2['Inheritance_model']['AD']
        assert not self.recessive_comp_not_simple_2['Inheritance_model']['AD_denovo']
        assert not self.recessive_comp_not_simple_2['Inheritance_model']['X']
        assert not self.recessive_comp_not_simple_2['Inheritance_model']['X_dn']
        assert self.recessive_comp_not_simple_2['Inheritance_model']['AR_compound']
    
    def test_not_recessive_comp(self):
        """docstring for test_not_recessive_comp"""
        assert not self.not_recessive_comp_1['Inheritance_model']['AR_hom']
        assert not self.not_recessive_comp_1['Inheritance_model']['AR_hom_denovo']
        assert self.not_recessive_comp_1['Inheritance_model']['AD']
        assert self.not_recessive_comp_1['Inheritance_model']['AD_denovo']
        assert not self.not_recessive_comp_1['Inheritance_model']['X']
        assert not self.not_recessive_comp_1['Inheritance_model']['X_dn']
        assert not self.not_recessive_comp_1['Inheritance_model']['AR_compound']
    
    def test_recessive_comp_missing(self):
        """docstring for test_recessive_comp_2"""
        assert not self.recessive_comp_missing_2['Inheritance_model']['AR_hom']
        assert not self.recessive_comp_missing_2['Inheritance_model']['AR_hom_denovo']
        assert self.recessive_comp_missing_2['Inheritance_model']['AD']
        assert self.recessive_comp_missing_2['Inheritance_model']['AD_denovo']
        assert not self.recessive_comp_missing_2['Inheritance_model']['X']
        assert not self.recessive_comp_missing_2['Inheritance_model']['X_dn']
        assert self.recessive_comp_missing_2['Inheritance_model']['AR_compound']


def main():
    pass


if __name__ == '__main__':
    main()

