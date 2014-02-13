#!/usr/bin/env python
# encoding: utf-8
"""
Variants.py

- Hold the information of some kind af variant.

- The type of variant should be checked in the Variant Parser

- A variant must allways have a coordinate in the form of chrom, start(pos), stop(pos)

- This object will have some mandatory attributes and a dictionary with extra info since different input formats will give very different information.

- Right now the formats that I can see that we will use are .vcf and the cmms standard which will probably change a bit from time to time. .bed will be harder to handle since they often do not include genotype info, often it is only a position.

There will be a base class called Variant with the required information and sub classes, Cmms and Vcf, that includes the special information from the different formats.


Attributes: 

Mandatory:

chr STRING 
start INT
stop INT 
ref STRING Reference nucleotide(s)
alt STRING Alternative sequence

Alternative:

identity STRING dbSNP-id
var_info DICT A dictionary with all the original info from the variant file
filter STRING The filter status
genotypes LIST A list with the genotypes found for this variants
ad BOOL If following Autosomal Dominant pattern
ad_dn BOOL If following Autosomal Dominant De novo pattern
ar BOOL If following Autosomal Recessive pattern
ar_dn BOOL If following Autosomal Recessive De nove pattern
ar_comp BOOL If following Autosomal Recessive compound pattern
ar_comp_dn BOOL If following Autosomal Recessive Compound De Novo pattern

Methods:

get_variant(self):
    Returns a dictionary with basic info to stdout

print_model_info(self):
    Print for each variant which patterns of inheritance they follow.    
    
print_vcf_variant(self):
    Print the variant in vcf-format to stdout

print_original_version(self, header_columns):
    Prints the variant in its original format.

check_noncomplete_call(self):
    Check if GATK have missed to report some info.

get_genotype(self):
    Returns the list with genotypes for this variant.

Created by Måns Magnusson on 2012-10-17.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from genmod.utils.is_number import is_number
from genmod.variants import genotype
from collections import OrderedDict

class Variant(object):
    """This class holds the necessary information about a genetic variant."""
    def __init__(self, chrom, start, stop, reference, alternative, identity='.', quality='.',
                filt='.', info=[], format_info=[], genes = [], all_info = OrderedDict()):
        super(Variant, self).__init__()
        # Properties:
        # If chromosome name is like Chr3 or chr5 we change to 3 or 5
        if 'hr' in chrom:
            self.chr = chrom[3:]
        else:
            self.chr = chrom # STRING
        # Make sure that start and stop are integers:
        # TODO raise proper exceptions here
        
        if type(start) == type(0):
            self.start = start # INT 
        elif is_number(start):
            self.start = int(start)
        else:
            print 'Start position', start, 'is not an integer.'
            sys.exit()
        self.ref = reference # Reference nucleotide(s) STRING
        self.alt = alternative # Alternative sequence LIST
        if identity == '-':
            self.identity = '.'
        else:
            self.identity = identity #dbSNP-id STRING
        self.quality = quality
        self.filter = filt
        self.genes = genes
        self.info = info # INFO field for vcf:s
        self.format_info = format_info# FORMAT field for vcf:s
        self.genotypes = OrderedDict() # Dict with {'ind_id':<Genotype>, 'ind_id':<Genotype>, ...}
        
        
        id_values = [self.chr, str(self.start), self.ref, self.alt[0]]
        self.variant_id = '_'.join(id_values)
        
        self.all_info = all_info
        
        # Models:
        
        if self.chr in ['X', 'x']:
            self.x_linked = True
            self.x_linked_dn = True
            self.ad = False
            self.ad_dn = False
            self.ar = False
            self.ar_dn = False
        else:
            self.ad = True #If following Autosomal Dominant pattern BOOL
            self.ad_dn = True #If following Autosomal Dominant De novo pattern BOOL
            self.ar = True #If following Autosomal Recessive pattern BOOL
            self.ar_dn = True #If following Autosomal Recessive De nove pattern BOOL
            self.x_linked = False
            self.x_linked_dn = False
        #If following Autosomal Recessive compound pattern for certain genes. DICT with {<Gene_ID>: [var_id_1, var_id_2, ...]]}
        self.ar_comp = False
        self.ar_comp_dn = False
        self.ar_comp_variants = {}#Dict with {variant_id:pair_score}
        
        # self.ar_comp_dn = True #If following Autosomal Recessive Compound De Novo pattern BOOL
        self.models = ['Na']
        self.rank_score = '-'
    
    def get_vcf_variant(self):
        """Return a list with information in vcf format"""
        vcf_info = [self.chr, str(self.start), self.identity, self.ref, ','.join(self.alt),
                    self.quality, self.filter, ';'.join(self.info), ':'.join(self.format_info)]
        
        for ind_id in self.genotypes:
            vcf_info.append(self.genotypes[ind_id].get_vcf_genotype())
        
        return vcf_info
        
    def get_cmms_variant(self):
        """Return a list with the information in the cmms format."""
        cmms_info = []
        for key in self.all_info:
            cmms_info.append(self.all_info[key])
        cmms_info.append(':'.join(self.models))
        if len(self.ar_comp_variants) > 0:
            compound_variants = []
            for compound_variant_id in self.ar_comp_variants:
                compound_variants.append(compound_variant_id+'='+str(self.ar_comp_variants[compound_variant_id]))
        else:
            compound_variants = ['-']
        cmms_info.append(':'.join(compound_variants))
        cmms_info.append(str(self.rank_score))
        return cmms_info
            
    def get_genotype(self, ind_id):
        """Return the genotype for a certain individual. Return a nocall if not existing."""
        return self.genotypes.get(ind_id, genotype.Genotype())
    
    def check_models(self):
        """Add the models to a list based on if they follow the pattern."""
        self.models = [] # If this function is run several times we empty the list each time
        if self.chr in ['x', 'X']:
            if self.x_linked:
                self.models.append('X')
            if self.x_linked_dn:
                self.models.append('X_denovo')
        else:
            if self.ad:
                self.models.append('AD')
            if self.ad_dn:
                self.models.append('AD_denovo')
            if self.ar:
                self.models.append('AR_hom')
            if self.ar_dn:
                self.models.append('AR_hom_denovo')
            if self.ar_comp:
                self.models.append('AR_compound')
        # if self.ar_comp_dn:  0:
        #     models.append('Autosomal Recessive Compound De Novo')
        if len(self.models) == 0:
            self.models.append('Na')
        
    
    def __str__(self):
        """Print information about the variant"""
        variant_info = [self.chr, str(self.start) , str(self.stop), self.alt, self.ref]
        for ind, geno in self.genotypes.items():
            geno_info = ind+':'+geno.genotype
            variant_info.append(geno_info)
        models = ';'.join(self.models)
        # for entry in self.all_info:
        #     if entry not in ['Chr', 'CHROM', 'Start', 'Stop', 'POS']:
        #         variant_info.append(self.all_info[entry])
        variant_info.append(models)
        variant_info.append(str(self.rank_score))
        variant_info.append(', '.join(self.genes))
        variant_info.append(', '.join(self.ar_comp_variants.keys()))
        return '\t'.join(variant_info)
    



def main():
    ranked_header = 'Chr    Start    Stop    Ref_allele    Alt_allele    IDN:Filter:GT=Genotype:AD=Allelic_depths_for_the_ref_and_alt_alleles:DP=Approximate_read_depth:GQ=Genotype Quality:PL=Normalized_Phred-scaled_likelihoods_for_genotypes    IDN:Filter:GT=Genotype:AD=Allelic_depths_for_the_ref_and_alt_alleles:DP=Approximate_read_depth:GQ=Genotype Quality:PL=Normalized_Phred-scaled_likelihoods_for_genotypes    IDN:Filter:GT=Genotype:AD=Allelic_depths_for_the_ref_and_alt_alleles:DP=Approximate_read_depth:GQ=Genotype Quality:PL=Normalized_Phred-scaled_likelihoods_for_genotypes    IDN:Filter:GT=Genotype:AD=Allelic_depths_for_the_ref_and_alt_alleles:DP=Approximate_read_depth:GQ=Genotype Quality:PL=Normalized_Phred-scaled_likelihoods_for_genotypes    HGNC_symbol    HGNC_approved_name    HGNC_Synonyms    Ensemble_GeneID    OMIM_gene_description    OMIM_morbid_description    Disease_group(dbIEMv.1.2)    HGMD    Gene_annotation    Functional_annotation    HGNC:Transcript_ID:Exon:cDNA_change:AA_change    phastConsElements    gerp++elem    genomicSuperDups    Pseudogene    1000G    dbsnp129    dbsnp132    dbsnp135    dbsnp135_Id    Esp6500    SIFT_Whole-exome    PolyPhen_version_2_HumDiv_Whole-exome    MutationTaster_Whole-exome    GERP++_Whole-exome    LRT_Whole-exome    PhyloP_Whole-exome    Main_Location_HPASubCL(120821)    Other_location_HPASubCL(120821)    Expression_type_HPASubCL(120821)    Reliability_HPASubCL(120821)    Gwas_catalog    Tfbs_transfac_Matrix_Database_(v7.0)    snoRNA_&_miRNA_annotations    GeneModel    Rank_Score    IDN:PASS;TotalFeaturesPASS||IDN:Chr:Start:Stop:Fraction_ME_Ten-Coverage_Bases:;TotalFeaturesPASS    IDN:PASS;TotalFeaturesPASS||IDN:Chr:Start:Stop:Fraction_ME_Ten-Coverage_Bases:;TotalFeaturesPASS    IDN:PASS;TotalFeaturesPASS||IDN:Chr:Start:Stop:Fraction_ME_Ten-Coverage_Bases:;TotalFeaturesPASS    IDN:PASS;TotalFeaturesPASS||IDN:Chr:Start:Stop:Fraction_ME_Ten-Coverage_Bases:;TotalFeaturesPASS'
    header_info = ranked_header.split('\t')
    example_variant = 'chr13    26153950    26153950    G    T    30-1-2A:PASS:0/1:3,0:42:42,0,93    30-1-3A:PASS:0/1:2,1:55:59,0,55    30-2-1U:PASS:0/1:3,0:24:24,0,68    30-2-2U:PASS:0/0:1,0:3:0,3,37    ATP8A2    ATPase, aminophospholipid transporter, class I, type 8A, member 2;    ATPIB, ML-1;    ENSG00000132932;    -;     ATPase, CLASS I, TYPE 8A, MEMBER 2; ATP8A2;    -;    -    exonic    nonsynonymous SNV    ATP8A2:NM_016529:exon22:c.1872G>T:p.L624F,    Score=446;Name=lod=87    Name=461.3    -    processed_pseudogene;    -    -    -    -    -    -    0    0.98    0.99321    5.49    1.0    0.999729    Nucleoli;    Nucleus but not nucleoli;    Staining;    Uncertain;    -    -    -    AR_compound;    25    30-2-2U:chr13:25946208:25946473:0.000:;48/49,     30-1-2A:chr13:25946208:25946473:0.611:;48/49,     30-1-3A:chr13:25946208:25946473:0.570:;48/49,     30-2-1U:chr13:25946208:25946473:0.000:30-2-1U:chr13:26153969:26154172:0.985:;47/49, '
    variant_info = example_variant.split('\t')
    variant_dict = {}
    for i in range(len(header_info)):
        print i, len(header_info), len(variant_info), header_info[i], variant_info[i]
        variant_dict[header_info[i]] = variant_info[i]
    my_variant = Cmms(**variant_dict)
    my_variant.print_original_info(header_info)
        


if __name__ == '__main__':
    main()

