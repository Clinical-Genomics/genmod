#!/usr/bin/env python
# encoding: utf-8
"""
genotype.py

This is a class with information about genotypecalls that follows the (GATK) .vcf standard.

The indata, that is the genotype call, is typically on the form x/x (unphased) or x|x (phased),
so they look like 0/0, 1/2, 1/1, 0|1 and so on.
The first sign inidcates what we find on the first allele, the second is a separator on the form '/' or '|' and the third indicates what is seen on the second allele.
The alleles are unordered.

Attributes:

    - genotype STRING (Same as in VCF-standard)
    - allele_1 STRING (Base on allele 1)
    - allele_2 STRING (Base on allele 2)
    - nocall BOOL
    - heterozygote BOOL
    - homo_alt BOOL (If individual is homozygote alternative)
    - homo_ref BOOL (If individual is homozygote reference)
    - has_variant BOOL (If individual is called and not homozygote reference)
    - ref_depth INT
    - alt_depth INT
    - phred_likelihoods LIST with FLOAT
    - depth_of_coverage INT
    - genotype_quality FLOAT
    - phased BOOL

If a variant is present, that is if homo_alt or heterozygote is true, then has_variant is True

When dealing with phased data we will see the '|'-delimiter


#TODO:
Should we allow '1/2', '2/2' and so on? This type of call looses it's point when moving from vcf -> bed since bed files only have one kind of variant on each line.
For now we will only allow './.', '0/0', '0/1', '1/1'

Created by Måns Magnusson on 2014-06-30.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""


class Genotype(object):
    """Holds information about a genotype"""

    def __init__(self, GT="./.", AD=".,.", DP="0", GQ="0", PL=None, **kwargs):
        GT = GT or "./."
        AD = AD or ".,."
        DP = DP or "0"
        GQ = GQ or "0"
        PL = PL or None

        self.phased = "|" in GT
        self.separator = "|" if self.phased else "/"

        # Split alleles safely
        if self.separator in GT:
            self.alleles = GT.split(self.separator)
        else:
            # Haploid calls (e.g. "0", "1") or missing (".")
            self.alleles = [GT] if GT else ["."]
        while len(self.alleles) < 2:
            self.alleles.append(".")

        self.allele_1 = self.alleles[0]
        self.allele_2 = self.alleles[1]

        # Flags
        self.genotyped = GT not in ("./.", ".|.", ".")
        self.heterozygote = False
        self.homo_alt = False
        self.homo_ref = False
        self.has_variant = False
        self.depth_of_coverage = 0
        self.quality_depth = 0
        self.genotype_quality = 0.0
        self.genotype = self.separator.join(self.alleles)

        # Check if a variant is homozygote reference, homozygote alternative or heterozygote
        if self.genotyped:
            if all(allele in ("0", ".") for allele in (self.allele_1, self.allele_2)):
                self.homo_ref = True
            # Treat e.g. ./1 or 1/. as heterozygous
            elif self.allele_1 == self.allele_2:
                self.homo_alt = True
            else:
                self.heterozygote = True

        self.has_variant = self.homo_alt or self.heterozygote

        # Parse allele depth
        self.ref_depth = 0
        self.alt_depth = 0
        try:
            allele_depths = [int(depth) if depth.isdigit() else 0 for depth in AD.split(",")]
            self.ref_depth = allele_depths[0] if len(allele_depths) > 0 else 0
            self.alt_depth = sum(allele_depths[1:]) if len(allele_depths) > 1 else 0
        except Exception:
            self.ref_depth = 0
            self.alt_depth = 0

        self.quality_depth = self.ref_depth + self.alt_depth

        # Check the depth of coverage:
        try:
            self.depth_of_coverage = int(DP)
        except ValueError:
            pass

        # Check the genotype quality
        try:
            self.genotype_quality = float(GQ)
        except ValueError:
            pass

        # Check the genotype likelihoods
        self.phred_likelihoods = []
        if PL:
            try:
                self.phred_likelihoods = [int(score) for score in PL.split(",")]
            except ValueError:
                pass

    def __str__(self):
        """Specifies what will be printed when printing the object."""
        return self.genotype
