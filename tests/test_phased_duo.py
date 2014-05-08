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
import pkg_resources

from pprint import pprint as pp

from ped_parser import family, individual
from interval_tree import interval_tree
import genmod
from genmod.models import genetic_models


class TestPhasedDuo(object):
    """Test how phased variants behave in a duo."""

    def setup_class(self):
        """Setup a simple family with family id 1, sick son id 1 and a healthy mother with id 3"""
        # Setup family with sick kid, sick father and healthy mother:
        self.duo_family = family.Family(family_id = '1')
        sick_son = individual.Individual(ind='1', family='1',mother='3', father='0', sex=1, phenotype=2)
        healthy_mother = individual.Individual(ind='3', family='1',mother='0', father='0', sex=2, phenotype=1)
        self.duo_family.add_individual(sick_son)
        self.duo_family.add_individual(healthy_mother)
        print(os.path.dirname(genmod.models.__file__))
        print(pkg_resources.get_distribution('genmod').version)
        print(sys.executable)
        
        # This is a dummy phased interval:
        interval = [1,100, '1']
        
        intervals = {ind_id:interval_tree.IntervalTree([interval], 1, 100) for ind_id in self.duo_family.individuals}
        
        #Setup two variants with autosomal recessive compound pattern
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
        
        batch['haploblocks'] = intervals
        # pp(self.duo_family)
        # pp(batch)
        # pp(sys.path)
        genetic_models.check_genetic_models(batch, self.duo_family, phased=True)
    
    def test_recessive_comp_duo_simple(self):
        """Check if the variant follows dominant pattern and compound pattern of inheritance."""
        assert not self.recessive_comp_simple_1['Inheritance_model']['AR_hom']
        assert not self.recessive_comp_simple_1['Inheritance_model']['AR_hom_dn']
        assert self.recessive_comp_simple_1['Inheritance_model']['AD']
        assert self.recessive_comp_simple_1['Inheritance_model']['AD_dn']
        pp(self.recessive_comp_simple_1['Inheritance_model'])
        assert self.recessive_comp_simple_1['Inheritance_model']['AR_comp']
    
    def test_recessive_comp_duo_not_simple(self):
        """Check in the variant follows the compound inheritance pattern."""
        assert not self.recessive_comp_not_simple_2['Inheritance_model']['AR_hom']
        assert not self.recessive_comp_not_simple_2['Inheritance_model']['AR_hom_dn']
        assert not self.recessive_comp_not_simple_2['Inheritance_model']['AD']
        assert not self.recessive_comp_not_simple_2['Inheritance_model']['AD_dn']
        assert self.recessive_comp_not_simple_2['Inheritance_model']['AR_comp']
    
    def test_not_recessive_comp_duo(self):
        """This variant should not follow the compound inheritance pattern."""
        assert not self.not_recessive_comp_1['Inheritance_model']['AR_hom']
        assert not self.not_recessive_comp_1['Inheritance_model']['AR_hom_dn']
        assert self.not_recessive_comp_1['Inheritance_model']['AD']
        assert self.not_recessive_comp_1['Inheritance_model']['AD_dn']
        assert not self.not_recessive_comp_1['Inheritance_model']['AR_comp']
    
    def test_recessive_comp_duo_missing(self):
        """Check if the variant follows the compound inheritance pattern."""
        assert not self.recessive_comp_missing_2['Inheritance_model']['AR_hom']
        assert not self.recessive_comp_missing_2['Inheritance_model']['AR_hom_dn']
        assert self.recessive_comp_missing_2['Inheritance_model']['AD']
        assert self.recessive_comp_missing_2['Inheritance_model']['AD_dn']
        assert self.recessive_comp_missing_2['Inheritance_model']['AR_comp']


def main():
    pass


if __name__ == '__main__':
    main()

