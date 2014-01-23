#!/usr/bin/env python
# encoding: utf-8
"""
get_family.py


Parse a file with family info, this can be a .ped file, a .fam, a .txt(CMMS style) 
file or a .txt(Broad style) file.
.ped and .fam always have 6 columns, these are

Family_ID Individual_ID Paternal_ID Maternal_ID 
Sex(1=male; 2=female; other=unknown) 
Phenotype(-9 missing, 0 missing, 1 unaffected, 2 affected)

The .txt allways have two columns on the form 

Individual_ID key=value

Where keys can be fid(=Family Id), mom, dad, sex, phenotype


If a pedigree file includes information about several families this must be taken care
 of by the parser by creating several family objects and then add information about the
  familymembers to theese familys. 

Create a family object and its family members from different types of input file
Created by Måns Magnusson on 2013-01-17.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import argparse
from genmod.family import family, individual

class FamilyParser(object):
    """Parses a file with family info and creates a family object with individuals."""
    def __init__(self, infile, family_type):
        super(FamilyParser, self).__init__()
        self.family_type = family_type
        self.families = {}
        self.header = []
        with open(infile, 'r') as f:
            line_count = 0
            for line in f:
                line_count += 1
                individual_line = line.rstrip()
                if line[0] != '#' and len(line) > 1:
                    if family_type == 'CMMS':
                        self.cmms_parser(individual_line, line_count, self.header)
                    elif family_type in ['FAM', 'PED']:
                        self.ped_parser(individual_line, line_count)
                    # elif family_type == 'broad':
                    #     self.broad_parser(individual_line, line_count)
                    my_counter = 1
                else:
                    self.header = line[1:].split()
    
    def ped_parser(self, individual_line, line_count):
        """Parse a .ped ped file."""
        if len(individual_line) < 6:
            print 'Malformed ped file, to few entrys on line ' + str(line_count)
            print individual_line
            sys.exit()
        line = individual_line.split()
        fam_id = line[0]
        if fam_id not in self.families:
            new_family = family.Family(fam_id)
            self.families[fam_id] = new_family
        ind = line[1]
        father = line[2]
        mother = line[3]
        sex = line[4]
        phenotype = line[5]
        my_individual = individual.Individual(ind, fam_id, mother, father, sex, phenotype)
        self.families[my_individual.family].add_individual(my_individual)


    def cmms_parser(self, individual_line, line_count, header):
        """Parse a .ped ped file."""
        if len(individual_line) < 6:
            print 'Malformed cmms file, to few entrys on line ' + str(line_count)
            print individual_line
            sys.exit()
        line = individual_line.split('\t')
        info = {}
        for i in range(len(line)):
            if header[i] == 'Inheritance_model':
                #If inheritance model is specified it is a ';'-separated list of models
                if line[i] != 'NA':
                    info[header[i]] = line[i].split(';')
            else:
                info[header[i]] = line[i]
        fam_id = info.get('FamilyID', '0')
        if fam_id not in self.families:
            self.families[fam_id] = family.Family(fam_id)
        ind = info.get('SampleID', '0')
        father = info.get('Father', '0')
        mother = info.get('Mother', '0')
        sex = info.get('Sex', '0')
        phenotype = info.get('Phenotype', '0')
        my_individual = individual.Individual(ind, fam_id, mother, father, sex, phenotype)
        models_of_inheritance = info.get('Inheritance_model', 'NA')
        if models_of_inheritance != 'NA':
            # this is something we will need to take care about while we change names
            # 'X', 'X_denovo', 'AD', 'AD_denovo', 'AR_hom', 'AR_hom_denovo', 'AR_compound'
            correct_model_names = []
            for model in models_of_inheritance:
                if model == 'AR':
                    model = 'AR_hom'
                elif model == 'AR_denovo':
                    moel = 'AR_hom_denovo'
                correct_model_names.append(model)
            self.families[my_individual.family].models_of_inheritance = correct_model_names
        self.families[my_individual.family].add_individual(my_individual)
    

def main():
    parser = argparse.ArgumentParser(description="Parse different kind of pedigree files.")
    parser.add_argument('pedigree_file', type=str, nargs=1 , help='A file with pedigree information.')
    parser.add_argument('-ped', '--ped', action="store_true", help='Pedigree file is in ped format.')
    parser.add_argument('-fam', '--fam', action="store_true", help='Pedigree file is in fam format.')
    parser.add_argument('-cmms', '--cmms', action="store_true", help='Pedigree file is in fam format.')
    parser.add_argument('-broad', '--broad', action="store_true", help='Pedigree file is in broad format.')
    args = parser.parse_args()
    infile = args.pedigree_file[0]
    file_type = 'cmms'
    if args.cmms:
        file_type = 'cmms'
    if args.ped:
        file_type = 'ped'
    if args.fam:
        file_type = 'ped'
    if args.broad:
        file_type = 'broad'
    my_parser = FamilyParser(infile, file_type)
    print my_parser.families
    for family in my_parser.families:
        print my_parser.families[family]
        print my_parser.families[family].models_of_inheritance
        


if __name__ == '__main__':
    main()
