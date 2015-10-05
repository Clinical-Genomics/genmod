## Conditions for Genetic Models ##

### Short explanation of genotype calls in VCF format:###

Since we only look at humans, that are diploid, the genotypes represent what we see on both alleles in a single position.
``0`` represents the reference sequence, ``1`` is the first of the alternative alleles, ``2`` second alternative and so on.
If no phasing has been done the genotype is an unordered pair on the form ``x/x``, so ``0/1`` means that the individual is heterozygote in this given position with the reference base on one of the alleles and the first of the alternatives on the other.
``2/2`` means that we see the second of the alternatives on both alleles.
Some chromosomes are only present in one copy in humans, here it is allowed to only use a single digit to show the genotype. A ``0`` would mean reference and ``1`` first of alternatives.

If phasing has been done the pairs are not unordered anymore and the delimiter is then changed to '|', so one can be heterozygote in two ways; ``0|1`` or ``1|0``.


### Autosomal Recessive ###

For this model individuals can be carriers so healthy individuals can be heterozygous. Both alleles need to have the variant for an individual to be sick so a healthy individual can not be homozygous alternative and a sick individual *has* to be homozygous alternative.

* Affected individuals have to be homozygous alternative (hom. alt.)
* Healthy individuals cannot be hom. alt.
* Variant is considered _de novo_ if both parents are genotyped and do not carry the variant
* The variants following this modell will be annotated with ``AR_hom``, or ``AR_hom_dn`` for de novo


### Autosomal Dominant ###

* Affected individuals have to be heterozygous (het.)
* Healthy individuals cannot have the alternative variant
* Variant is considered _de novo_ if both parents are genotyped and do not carry the variant
* The variants following this modell will be annotated with ``AD``or ``AD_dn``

**Special Case:** If the variant belongs to a gene that have been annotated 
with reduced penetrance we will allow healthy carriers. (i.e. healthy individuals can be heterozygotes)

### Autosomal Compound Heterozygote ###

This model includes pairs of exonic variants that are present within the same gene.
**The default behaviour of GENMOD is to look for compounds only in exonic/canonical splice sites**.
The reason for this is that since some genes have huge intronic regions the data will be drowned in compound pairs.
If the user wants all variants in genes checked use the flag -gene/--whole_gene.

1. Non-phased data:
	* Affected individuals have to be het. for both variants
	* Healthy individuals can be het. for one of the variants but cannot have both variants
	* Variant is considered _de novo_ if only one or no variant is found in the parents


2. Phased data:
	* All affected individuals have to be het. for both variants **and** the variants has to be on two different alleles
	* Healthy individuals can be heterozygous for one but cannot have both variants
	* If only one or no variant is found in parents it is considered _de novo_


Variant following this patters will be annotated with ``AR_comp`` or ``AR_comp_dn``

### X-Linked Dominant###

These traits are inherited on the x-chromosome, of which men have one allele and women have two.

* Variant has to be on chromosome X
* Affected individuals have to be het. or hom. alt.
* Healthy males cannot carry the variant
* Healthy females can carry the variant (because of X inactivation)
* If sex is male the variant is considered _de novo_ if mother is genotyped and does not carry the variant
* If sex is female variant is considered _de novo_ if none of the parents carry the variant
* The variants following this modell will be annotated with ``XD``or ``XD_dn``


### X Linked Recessive ###

* Variant has to be on chromosome X
* Affected males have to be het. or hom. alt. (het is theoretically not possible in males, but can occur due to Pseudo Autosomal Regions).
* Healthy males cannot carry the variant
* Affected females have to be hom. alt.
* Healthy females cannot be hom. alt.
* If sex is male the variant is considered _de novo_ if mother is genotyped and does not carry the variant
* If sex is female variant is considered _de novo_ if not both parents carry the variant
* The variants following this modell will be annotated with ``XR``or ``XR_dn``

