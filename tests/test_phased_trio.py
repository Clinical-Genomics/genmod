#!/usr/bin/env python
# encoding: utf-8
"""
test_models_phased_compound.py

Test the so that the genetic models behave as suspected.


Created by Måns Magnusson on 2013-03-07.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from ped_parser import family, individual
from interval_tree import interval_tree
from genmod import genetic_models
from vcf_parser import genotype


class TestPhasedTrio(object):
    """Check how the inheritnace patterns are followed in a phased trio."""

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
        
        interval = [1,100, '1']
        intervals = {ind_id:interval_tree.IntervalTree([interval], 1, 100) for ind_id in self.recessive_family.individuals}
        
        #Setup two variants with only autosomal recessive pattern
        self.recessive_comp_simple_1 = {'CHROM':'1', 'POS':'5', 'ALT':'A', 'REF':'C', 'ID':'rs2230749',
                                                 '1':'0|1', '2':'0|1', '3':'0|0'}
        genotypes = {'1':genotype.Genotype('0|1'), '2':genotype.Genotype('0|1'), '3':genotype.Genotype('0|0')}
        self.recessive_comp_simple_1['genotypes'] = genotypes
                
        self.recessive_comp_simple_2 = {'CHROM':'1', 'POS':'10', 'ALT':'C', 'REF':'T', 'ID':'.', 
                                                '1':'1|0', '2':'0|0', '3':'0|1'}
        genotypes = {'1':genotype.Genotype('1|0'), '2':genotype.Genotype('0|0'), '3':genotype.Genotype('0|1')}
        self.recessive_comp_simple_2['genotypes'] = genotypes
        
        #Setup two variants that parents have but on same allele
        self.recessive_comp_not_simple_1 = {'CHROM':'1', 'POS':'15', 'ALT':'A', 'REF':'C', 'ID':'.',
                                                    '1':'0|1', '2':'0|1', '3':'1|0'}
        genotypes = {'1':genotype.Genotype('0|1'), '2':genotype.Genotype('0|1'), '3':genotype.Genotype('1|0')}
        self.recessive_comp_not_simple_1['genotypes'] = genotypes
        
        self.recessive_comp_not_simple_2 = {'CHROM':'1', 'POS':'17', 'ALT':'C', 'REF':'T', 'ID':'.', 
                                                    '1':'1|0', '2':'0|1', '3':'1|0'}
        genotypes = {'1':genotype.Genotype('1|0'), '2':genotype.Genotype('0|1'), '3':genotype.Genotype('1|0')}
        self.recessive_comp_not_simple_2['genotypes'] = genotypes
        
        #This one should not be recessive since the sick have both variants on same allele:
        self.not_recessive_comp_1 =  {'CHROM':'1', 'POS':'27', 'ALT':'G', 'REF':'T', 'ID':'.', 
                                        '1':'0|1', '2':'0|1', '3':'0|0'}
        genotypes = {'1':genotype.Genotype('0|1'), '2':genotype.Genotype('0|1'), '3':genotype.Genotype('0|0')}
        self.not_recessive_comp_1['genotypes'] = genotypes
        
        self.not_recessive_comp_2 = {'CHROM':'1', 'POS':'35', 'ALT':'C', 'REF':'T', 'ID':'.', 
                                        '1':'0|1', '2':'0|0', '3':'0|1'}
        genotypes = {'1':genotype.Genotype('0|1'), '2':genotype.Genotype('0|0'), '3':genotype.Genotype('0|1')}
        self.not_recessive_comp_2['genotypes'] = genotypes
        
        #Setup two variants that parents have but on same allele
        self.recessive_comp_missing_1 = {'CHROM':'1', 'POS':'45', 'ALT':'A', 'REF':'C', 'ID':'.',
                                                    '1':'0|1', '2':'./.', '3':'0|0'}
        genotypes = {'1':genotype.Genotype('0|1'), '2':genotype.Genotype('./.'), '3':genotype.Genotype('0|0')}
        self.recessive_comp_missing_1['genotypes'] = genotypes
        
        self.recessive_comp_missing_2 = {'CHROM':'1', 'POS':'47', 'ALT':'C', 'REF':'T', 'ID':'.', 
                                                    '1':'1|0', '2':'0|1', '3':'1|0'}
        genotypes = {'1':genotype.Genotype('1|0'), '2':genotype.Genotype('0|1'), '3':genotype.Genotype('1|0')}
        self.recessive_comp_missing_2['genotypes'] = genotypes
        
        batch = {'ABC':{'1_5_A_C':self.recessive_comp_simple_1, '1_10_C_T':self.recessive_comp_simple_2},
                'BBC':{'1_15_A_C':self.recessive_comp_not_simple_1, '1_17_C_T':self.recessive_comp_not_simple_2},
                'CBC':{'1_27_G_T':self.not_recessive_comp_1, '1_35_C_T':self.not_recessive_comp_2},
                'DBC':{'1_45_A_C':self.recessive_comp_missing_1, '1_47_C_T':self.recessive_comp_missing_2}
                }
        
        batch['haploblocks'] = intervals
        genetic_models.check_genetic_models(batch, self.recessive_family, phased=True)
    
    def test_recessive_comp_simple(self):
        """Check if the compound inheritance pattern are followed."""
        assert not self.recessive_comp_simple_1['Inheritance_model']['AR_hom']
        assert not self.recessive_comp_simple_1['Inheritance_model']['AR_hom_dn']
        assert not self.recessive_comp_simple_1['Inheritance_model']['AD']
        assert not self.recessive_comp_simple_1['Inheritance_model']['AD_dn']
        assert self.recessive_comp_simple_1['Inheritance_model']['AR_comp']
    
    def test_recessive_comp_not_simple(self):
        """Check if the compound inheritance pattern are followed."""
        assert not self.recessive_comp_not_simple_2['Inheritance_model']['AR_hom']
        assert not self.recessive_comp_not_simple_2['Inheritance_model']['AR_hom_dn']
        assert not self.recessive_comp_not_simple_2['Inheritance_model']['AD']
        assert not self.recessive_comp_not_simple_2['Inheritance_model']['AD_dn']
        assert not self.recessive_comp_not_simple_2['Inheritance_model']['AR_comp']
        
    
    def test_not_recessive_comp(self):
        """This variants should not follow any inheritance patterns."""
        assert not self.not_recessive_comp_1['Inheritance_model']['AR_hom']
        assert not self.not_recessive_comp_1['Inheritance_model']['AR_hom_dn']
        assert not self.not_recessive_comp_1['Inheritance_model']['AD']
        assert not self.not_recessive_comp_1['Inheritance_model']['AD_dn']
        assert not self.not_recessive_comp_1['Inheritance_model']['AR_comp']
    
    def test_recessive_comp_missing(self):
        """This variant should follow the compound inheritnce pattern."""
        assert not self.recessive_comp_missing_2['Inheritance_model']['AR_hom']
        assert not self.recessive_comp_missing_2['Inheritance_model']['AR_hom_dn']
        assert not self.recessive_comp_missing_2['Inheritance_model']['AD']
        assert not self.recessive_comp_missing_2['Inheritance_model']['AD_dn']
        assert self.recessive_comp_missing_2['Inheritance_model']['AR_comp']


def main():
    pass


if __name__ == '__main__':
    main()

