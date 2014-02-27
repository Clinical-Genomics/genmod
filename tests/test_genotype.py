#!/usr/bin/env python
# encoding: utf-8
"""
test_genotype_pytest.py

Attributes of genotypes to be tested:

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
	- phred_likelihoods TUPLE with INT
	- depth_of_coverage INT
	- genotype_quality FLOAT

Created by Måns Magnusson on 2013-02-26.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from genmod.variants import genotype

def test_nocall():
	"""A nocall is when no informations is found on this position for the individual. It should be False on all questions except nocall. Also in the case of haploidity the result should be the same."""
	my_nocall = genotype.Genotype('./.')
	my_short_nocall = genotype.Genotype('.')
	assert my_nocall.genotype == './.' #We never need to look at the alleles since genotype is defined by 'allele_1/allele_2'
	assert my_short_nocall.genotype == './.'
	assert not my_nocall.heterozygote
	assert not my_short_nocall.heterozygote
	assert not my_nocall.homo_ref
	assert not my_short_nocall.homo_ref
	assert not my_nocall.homo_alt
	assert not my_short_nocall.homo_alt
	assert not my_nocall.has_variant
	assert not my_short_nocall.has_variant
	assert my_nocall.nocall
	assert my_short_nocall.nocall

def test_genotype_1():
	"""The vcf format allows calls that look like '1'. This is a haploid call and should be treated as Homozygote Alternative."""
	my_genotype = genotype.Genotype('1')
	assert my_genotype.genotype == '1/.'
	assert not my_genotype.nocall
	assert not my_genotype.heterozygote
	assert not my_genotype.homo_ref
	assert my_genotype.homo_alt
	assert my_genotype.has_variant

def test_genotype_0():
	"""The vcf format allows calls that look like '0'. This is a haploid call and should be treated as Homozygote Reference."""
	my_genotype = genotype.Genotype('0')
	assert my_genotype.genotype == '0/.'
	assert not my_genotype.nocall
	assert not my_genotype.heterozygote
	assert my_genotype.homo_ref
	assert not my_genotype.homo_alt
	assert not my_genotype.has_variant

def test_genotype_0_1():
	"""docstring for test_genotype_0_1. A normal heterozygote call, has_variant and heterozygote is true."""
	my_genotype = genotype.Genotype('0/1')
	assert my_genotype.genotype == '0/1'
	assert not my_genotype.nocall
	assert my_genotype.heterozygote
	assert not my_genotype.homo_ref
	assert not my_genotype.homo_alt
	assert my_genotype.has_variant

def test_genotype_1_2():
	"""docstring for test_genotype_1_2. A normal heterozygote call, has_variant and heterozygote is true."""
	my_genotype = genotype.Genotype('1/2')
	assert my_genotype.genotype == '1/2'
	assert not my_genotype.nocall
	assert my_genotype.heterozygote
	assert not my_genotype.homo_ref
	assert not my_genotype.homo_alt
	assert my_genotype.has_variant

def test_homo_ref():
	"""A homozygote reference call. has_variant and nocall is False and homo_ref is true."""
	my_homo_ref_genotype = genotype.Genotype('0/0')
	assert my_homo_ref_genotype.genotype == '0/0'
	assert not my_homo_ref_genotype.nocall
	assert not my_homo_ref_genotype.heterozygote
	assert my_homo_ref_genotype.homo_ref
	assert not my_homo_ref_genotype.homo_alt
	assert not my_homo_ref_genotype.has_variant

def test_homo_alt():
	"""A homozygote alternative call. has_variant and homo_alt is true."""
	my_genotype = genotype.Genotype('1/1')
	assert my_genotype.genotype == '1/1'
	assert not my_genotype.nocall
	assert not my_genotype.heterozygote
	assert not my_genotype.homo_ref
	assert my_genotype.homo_alt
	assert my_genotype.has_variant

def test_homo_alt_2():
	"""A homozygote alternative call. has_variant and homo_alt is true."""
	my_genotype = genotype.Genotype('3/3')
	assert my_genotype.genotype == '3/3'
	assert not my_genotype.nocall
	assert not my_genotype.heterozygote
	assert not my_genotype.homo_ref
	assert my_genotype.homo_alt
	assert my_genotype.has_variant

def test_phased_data():
	"""Try if the class van handle phased data. In this case a heterozygote."""
	my_genotype = genotype.Genotype('1|0')
	assert my_genotype.genotype == '1/0'# If asked about the genotype, it should still be on the same form.
	assert not my_genotype.nocall
	assert my_genotype.heterozygote
	assert not my_genotype.homo_ref
	assert not my_genotype.homo_alt
	assert my_genotype.has_variant

	
def main():
	pass


if __name__ == '__main__':
	main()

