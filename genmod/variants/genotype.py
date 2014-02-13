#!/usr/bin/env python
# encoding: utf-8
"""
genotype.py

This is a class with information about genotypecalls that follows the (GATK) .vcf standard.

The indata, that is the genotype call, is allways on the form x/x, so they look like 0/0, 1/2, 1/1 and so on.
The first sign inidcates what we find on the first allele, the second is a separator on the form '/' or '|' and the third indicates what is seen on the second allele.
The alleles are unordered.

Attributes:

    - genotype STRING
    - allele_1 STRING
    - allele_2 STRING
    - nocall BOOL
    - heterozygote BOOL
    - homo_alt BOOL
    - homo_ref BOOL
    - has_variant BOOL
    - filter STRING
    - ref_depth INT
    - alt_depth INT
    - phred_likelihoods LIST with FLOAT
    - depth_of_coverage INT
    - genotype_quality FLOAT

If a variant is present, that is if homo_alt or heterozygote is true, then has_variant is True
    
Maybe we can include other genotype info here, not sure yet...

When dealing with phased data we will see the '|'-delimiter


#TODO:
Should we allow '1/2', '2/2' and so on? This type of call looses it's point when moving from vcf -> bed since bed files only have one kind of variant on each line.
For now we will only allow './.', '0/0', '0/1', '1/1'   

Created by Måns Magnusson on 2013-02-18.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os


class Genotype(object):
    """Holds information about a genotype"""
    def __init__(self, GT='./.', AD='.,.', DP='0', GQ='0', PL=[],
                 allele_1 = '.', allele_2 = '.', FILTER = '.', original_info = ''):
        super(Genotype, self).__init__()        
        # Genotype call, ./., 1/1, 0/1, 0|1 ...
        self.genotype = GT
        self.allele_1 = '.'
        self.allele_2 = '.'
        self.allele_1_base = allele_1
        self.allele_2_base = allele_2
        self.ref_depth = '.'
        self.alt_depth = '.'
        self.original_info = original_info
        if len(AD) > 2:
            if AD[0].isdigit():
                self.ref_depth = int(AD.split(',')[0])
            if AD[2].isdigit():
                self.alt_depth = int(AD.split(',')[1])
        self.depth_of_coverage = int(DP)
        self.genotype_quality = float(GQ)
        self.phred_likelihoods = []
        if PL :
            for score in PL.split(','):
                self.phred_likelihoods.append(int(score))
        
        # These are the different genotypes:
        self.nocall = True
        self.heterozygote = False
        self.homo_alt = False
        self.homo_ref = False
        self.has_variant = False
        
        if len(self.genotype) < 3: #This is the case when only one allele is present(eg. X-chromosome) and presented like '0' or '1'.
            self.allele_1 = self.genotype
            self.allele_2 = '.'
        else:
            self.allele_1 = self.genotype[0]
            self.allele_2 = self.genotype[-1]
        self.genotype = self.allele_1 +'/'+ self.allele_2 # The genotype should allways be represented on the same form
        if self.genotype != './.':
            self.nocall = False
            self.check_alleles(self.allele_1, self.allele_2)
            self.check_alleles(self.allele_2, self.allele_1)
        if self.heterozygote or self.homo_alt:
            self.has_variant = True
    
    def check_alleles(self, variant1, variant2):
        """Check if the genotype is heterozygote, homozygote etc..."""
        if variant1 == '.':# First is the case with './x' or 'x/.'
            if variant2 == '0':
                self.homo_ref = True
            else:
                self.homo_alt = True
        elif variant1 == '0':
            if variant2 == variant1:
                self.homo_ref = True
            elif variant2 != '.':
                self.heterozygote = True
        else:
            if variant1 == variant2:
                self.homo_alt = True
            elif variant2 != '.':
                self.heterozygote = True
    
    def get_vcf_genotype(self):
        """Returns the genotype in the original vcf-format"""
        
        return self.original_info
                    
            
    
    def __str__(self):
        """Specifies what will be printed when printing the object."""
        return self.allele_1+'/'+self.allele_2

def main():
    pass


if __name__ == '__main__':
    main()

