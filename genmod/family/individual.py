#!/usr/bin/env python
# encoding: utf-8
"""
Individual.py

Holds the information of an individual

Attributes:

ind STRING Can be any id unique within the family
family STRING Can be any unique id within the cohort
mother STRING The ind_id of the mother or [0,-9] if info is missing
father STRING ---------||------ father --------------||---------------
sex INT 1=male 2=female 0=unknown
phenotype INT 1=unaffected, 2=affected, missing = [0,-9]
genotypes DICT Container with genotype information on the form {<variant_id>: <Genotype>}
phasing BOOL If the genotype information includes phasing for this individual

Created by Måns Magnusson on 2012-10-31.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from genmod.variants import genotype



class Individual(object):
    """docstring for Individual"""
    def __init__(self, ind, family='0', mother='0', father='0',sex=0,phenotype=0):
        
        #TODO write test to throw exceptions if malformed input.
        
        self.individual_id = ind #Individual Id STRING
        self.family = family #Family Id STRING
        self.mother = mother #Mother Id STRING
        self.father = father # Father Id STRING
        if type(sex) == type(0):
            self.sex = sex # Sex Integer
        else:
            self.sex = int(sex) # Sex Integer            
        if type(phenotype) == type(0):
            self.phenotype = phenotype # Phenotype INTEGER 
        else:
            self.phenotype = int(phenotype)
        self.genotypes = {} # DICTIONARY Container with genotypes for all variants on the form {<variant_id>:<Genotype object>}
        self.phasing = False # If we have phasing info for this individual BOOL
    
    def add_genotype(self, variant_id, genotype):
        """Add the information of a variant to this individual."""
        self.genotypes[variant_id] = genotype
    
    def get_genotype(self, variant_id):
        """Returns a genotype object for the variant with the given ID and a default(./.) genotype object if not present."""
        if variant_id in self.genotypes:
            return self.genotypes[variant_id]
        else:
            return genotype.Genotype()
    
    def get_info(self):
        """Return the information about this individual as a dictionary with the attributes."""
        return {'ind_id':self.ind_id, 'family':self.family, 'mother': self.mother, 'father':self.father, 'sex':self.sex, 'phenotype':self.phenotype, 'phasing':self.phasing}
    
    def affected(self):
        """Returns true is affected and false if healthy or unknown(?)"""
        if self.phenotype == 2:
            return True
        else:
            return False
    
    def __str__(self):
        """Returns what should be printed if object is printed."""
        ind_info = ['ind:', self.ind_id, 'family:', self.family, 'mother:', self.mother, 'father:', self.father,
         'sex:', str(self.sex), 'phenotype:', str(self.phenotype)]
        return ' '.join(ind_info)

def main():
    pass


if __name__ == '__main__':
    main()

