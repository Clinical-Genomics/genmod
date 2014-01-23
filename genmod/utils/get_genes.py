#!/usr/bin/env python
# encoding: utf-8
"""
get_genes.py

Function that parses a string of genes and returns a list of genes.
 
Created by Måns Magnusson on 2014-01-14.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os

def get_genes(gene_string, gene_annotation = 'Ensembl'):
    """Return a list with the genes that this variant is present in.
    gene_annotation in ['Ensembl', 'HGNC']"""
    genes = []
    
    if gene_annotation == 'Ensembl':
        delimiter = ';'
        if len(gene_string) > 1:
            if delimiter == gene_string[-1]:
                genes = gene_string.split(delimiter)[:-1]
        else:
            genes = gene_string.split()
                    
    elif gene_annotation == 'HGNC':
        
        def parse_gene(gene_entry):
            """Parse out the relevant genes and return a list with them."""
            
            def correct_entry(entry):
                """Removes parenthesis and the info in them from an entry."""
                return entry.split('(')[0]
                
            def insert_gene(gene):
                """Insert a gene to the gene list if it is not already there."""
                if gene not in genes:
                    genes.append(gene)
                    
            semi_splitted = gene_entry.split(';')
            if len(semi_splitted) == 1:
                for entry in semi_splitted:
                    entry = correct_entry(entry)
                    comma_splitted = entry.split(',')
                    for entry in comma_splitted:
                        entry = correct_entry(entry)
                        minus_splitted = entry.split('-')
                        for entry in minus_splitted:
                            entry = correct_entry(entry)
                            insert_gene(entry)
            else:
                for entry in semi_splitted:
                    parse_gene(entry)
        
        if 'dist' in gene_string or 'NONE' in gene_string:
            genes = ['-']
        else:
            parse_gene(gene_string)
        
    if len(genes) == 0:
        genes = ['-']
     
    return genes    

def main():
    pass


if __name__ == '__main__':
    main()

