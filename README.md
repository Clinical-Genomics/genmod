# GENMOD #


Tool for annotating patterns of inheritance Variant Call Format (VCF) files.

Each variant in the VCF-file will be annotated with which genetic models that are followed in the family.
The genetic models that are checked are the following:

* Autsomal Recessive, denoted 'AR'
* Autsomal Recessive denovo, denoted 'AR_dn'
* Autsomal Dominant, 'AD'
* Autsomal Dominant denovo, 'AD_dn'
* Autosomal Compound Heterozygote, 'AR_comp'
* X-linked dominant, 'XD'
* X-linked dominant de novo, 'XD_dn'
* X-linked Recessive, 'XR'
* X-linked Recessive de novo, 'XR_dn'

**GENMOD** will add entrys to the INFO column for the given VCF file. The new entrys are: 
    
- GM: A colon separated list with genetic models followed
- ANN: Colon separated list with features overlapped in the annotation file
- Comp: Colon separated list with compound pairs(if any). These are described like CHR_POS_REF_ALT.
- MS: Model Score, a phred-score based on the genotype qualities to describe the uncertainty of the genetic model.


##Installation:##

genmod works with Python 2.7 and Python 3.

    pip install genmod

or


    git clone git@github.com:moonso/genmod.git
    cd genmod
    python setup.py install

###USAGE:###

    genmod ped_file variant_file annotation_file


## Conditions for Genetic Models ##

### Short explanation of genotype calls in VCF format:###

Since we only look at humans, that are diploid, the genotypes represent what we see on both alleles in a single position.
0 represents the reference sequence, 1 is the first of the alternative alleles, 2 second alternative and so on.
If no phasing has been done the genotype is an unordered pair on the form x/x, so 0/1 means that the individual is heterozygote in this given position with the reference base on one of the alleles and the first of the alternatives on the other.
2/2 means that we see the second of the alternatives on both alleles.
Some chromosomes are only present in one copy in humans, here it is allowed to only use a single digit to show the genotype. A 0 would mean reference and 1 first of alternatives.

If phasing has been done the pairs are not unordered anymore and the delimiter is then changed to '|', so one can be heterozygote in two ways; 0|1 or 1|0.


### Autosomal Recessive ###

For this model individuals can be carriers so healthy individuals can be heterozygous. Both alleles need to have the variant for an individual to be sick so a healthy individual can not be homozygous alternative and a sick individual *has* to be homozygous alternative.

* Affected individuals have to be homozygous alternative (hom. alt.)
* Healthy individuals cannot be hom. alt.
* Variant is considered _de novo_ if both parents are genotyped and do not carry the variant


### Autosomal Dominant ###

* Affected individuals have to be heterozygous (het.)
* Healthy individuals cannot have the alternative variant
* Variant is considered _de novo_ if both parents are genotyped and do not carry the variant


### Autosomal Compound Heterozygote ###

This model includes pairs of variants that are present within the same feature.

1. Non-phased data:
	* Affected individuals have to be het. for both variants
	* Healthy individuals can be het. for one of the variants but cannot have both variants
	* Variant is considered _de novo_ if only one or no variant is found in the parents
  	
  
2. Phased data:
	* All affected individuals have to be het. for both variants and the variants has to be on two different alleles
	* Healthy individuals can be heterozygous for one or both of the variants if they are on the \emph{same} allele
	* If only one or no variant is found in parents it is considered _de novo_


### X-Linked Dominant###

These traits are inherited on the x-chromosome, of which men have one allele and women have two. 

* Variant has to be on chromosome X
* Affected individuals have to be het. or hom. alt.
* Healthy males cannot carry the variant
* Healthy females can carry the variant
* If sex is male the variant is considered _de novo_ if mother is genotyped and does not carry the variant
* If sex is female variant is considered _de novo_ if none of the parents carry the variant
    

### X Linked Recessive ###

* Variant has to be on chromosome X
* Affected males have to be het. or hom. alt.
* Affected females have to be hom. alt.
* Healthy females cannot be hom. alt.
* Healthy males cannot carry the variant
* If sex is male the variant is considered _de novo_ if mother is genotyped and does not carry the variant
* If sex is female variant is considered _de novo_ if not both parents carry the variant



<!-- ## Detailed Structure ##

Here all attributes and methods of the classes will be showed:

### Genotype ###

Store the genotype information of a variant that is specific for an individual

**Attributes:**

* genotype STRING
* allele_1 STRING
* allele_2 STRING
* nocall BOOL
* heterozygote BOOL
* homo_alt BOOL
* homo_ref BOOL
* has_variant BOOL
* filter STRING
* ref_depth INT
* alt_depth INT
* phred_likelihoods TUPLE with INT
* depth_of_coverage INT
* genotype_quality FLOAT


### Variant ###

Holds the info of a variant and it's specific behaviour in a family.

**Attributes**

* chr STRING Have to be string since X, Y
* start INT
* stop INT 
* ref STRING Reference nucleotide(s)
* alt STRING Alternative sequence
* identity STRING dbSNP-id
* var_info DICT A dictionary with all the info from the variant file
* qual STRING A value for the score of the base call
* filter STRING The filter status
* genotypes LIST A list with the genotypes found for this variants
* gene STRING Semicolon separated string with ensemble gene names
* ad BOOL If following Autosomal Dominant pattern
* ad_dn BOOL If following Autosomal Dominant De novo pattern
* ar BOOL If following Autosomal Recessive pattern
* ar_dn BOOL If following Autosomal Recessive De nove pattern
* ar_comp BOOL If following Autosomal Recessive compound pattern
* ar_comp_dn BOOL If following Autosomal Recessive Compound De Novo pattern

**Methods**

* get_variant(self):
    Returns a dictionary with basic info to stdout
* print_model_info(self):
    Print for each variant which patterns of inheritance they follow.    
* print_vcf_variant(self):
    Print the variant in vcf-format to stdout
* print_original_version(self, header_columns):
    Prints the variant in its original format.
* check_noncomplete_call(self):
    Check if GATK have missed to report some info.
* get_genotype(self):
    Returns the list with genotypes for this variant.


### Individual ###

Holds the information about an individual and the individual specific genotypes.

**Attributes**

* ind STRING Can be any id unique within the family
* family STRING Can be any unique id within the cohort
* mother STRING The ind_id of the mother or [0,-9] if info is missing
* father STRING ---------||------ father --------------||---------------
* sex INT 1=male 2=female 0=unknown
* phenotype INT 1=unaffected, 2=affected, missing = [0,-9]
* genotypes DICT Container with genotype information on the form {<variant_id>: <Genotype>}
* phasing BOOL If the genotype information includes phasing for this individual

### Family ###

**Attributes**

* individuals DICT dictionary with family members on the form {<ind_id>:<Individual>}
* variants DICT dictionary with all the variants that exists in the family on the form {<var_id>:<Variant>} -->


[![Bitdeli Badge](https://d2weczhvl823v0.cloudfront.net/moonso/genmod/trend.png)](https://bitdeli.com/free "Bitdeli Badge")

