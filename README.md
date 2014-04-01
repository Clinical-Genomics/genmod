# GENMOD #


Tool for annotating patterns of inheritance Variant Call Format (VCF) files.

Each variant in the VCF-file will be annotated with which genetic models that are followed in the family.
The genetic models that are checked are the following:

    * Autsomal Recessive, denoted 'AR'
    * Autsomal Recessive denovo, denoted 'AR_denovo'
    * Autsomal Dominant, 'AD'
    * Autsomal Dominant denovo, 'AD_denovo'
    * Autosomal Compound Heterozygote, 'AR_compound'
    * X-linked dominant, 'XD'
    * X-linked dominant de novo, 'XD_denovo'
    * X-linked Recessive, 'XR'
    * X-linked Recessive de novo, 'XR_denovo'

**GENMOD** will add entrys to the INFO column for the given VCF file. The new entrys are: 
    
    - GM: A colon separated list with genetic models followed
    - ANN: Colon separated list with features overlapped in the annotation file
    - Comp: Colon separated list with compound pairs(if any). These are described like CHR_POS_REF_ALT.
    - MS: Model Score, a phred-score based on the genotype qualities to describe the uncertainty of the genetic model.


###Installation:###

genmod works with Python 2.7 and Python 3.


```
pip install genmod
´´´

or

```
 git clone git@github.com:moonso/genmod.git
 cd genmod
 python setup.py install

genmod ped_file variant_file annotation_file

```


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
    * Variant is considered \emph{de novo} if both parents are genotyped and do not carry the variant

### Autosomal Recessive De Novo ###

Same as above but with the difference that one or both of the parents are missing the genotypes.

### Autosomal Dominant ###

Here it is enough that one of the alleles have the variant for an individual to be sick.

    


### Autosomal Dominant De Novo ###

Same as above but with the difference that one or both of the parents are missing the genotypes.


### Autosomal Compound Heterozygote ###

Here we are looking at pairs of *heterozygote* variants that are located within the same gene. For two variants to qualify for being compound pair they have to be present together in a affected individual but can **not** be present together in any of the healthy individuals.
So *one and only one*, variant should be found in each parent, otherwise the variant follows the pattern of compound heterozygote de novo.

One of the problems here is that a compound heterozygote can exist in different ways such that they can be deleterious or not.  
If we do not have the phasing information from the individuals, that is if we have not resolved on which alleles the variants that we look at are located on, we can not say if a pair is deleterious or not.

Right now this is how we do this:

1. If individual is healthy:
	* Can be a carrier so 0/1 is ok
	* Can have ref call or no call so 0/0 and ./. is ok
	* Can not be homozygote alternative, so 1/1 is NOT ok.
  	
  
2. If individual is sick:
	* Can be a carrier so 0/1 is not ok.
	* Can not have ref call, 0/0 is not ok
	* Must make sure that each of the parents carry one (and only one) of the variants in the potential pair each.
	* No call, ./., is ok since we can not exclude the model by this.
    
### Autosomal Compound Heterozygote De Novo ###

Same as above but one of the variants are not seen in the parents.

### X Linked ###

These traits are inherited on the x-chromosome, of which men have one allele and women have two. This means that woman's can be carries but men get affected if they get the variant. It is rare that women get these type of disease since they have to inherit a variant from both parents which implies that the father has to be sick. 

1. If individual is healthy:
	* If woman: can be a carrier so 0/1 is ok
	* If man: can not be a carrier so 0/1 is not ok
	* Can have ref call or no call so 0/0 and ./. is ok
	* Can not be homozygote alternative, so 1/1 is NOT ok.
  	
  
2. If individual is sick:
	* If woman: can be a carrier so 0/1 is not ok.
	* Can not have ref call, 0/0 is not ok
	* No call, ./., is ok since we can not exclude the model by this.



### X Linked De Novo ###

Same as above but variant is not seen in any of the parents.

## Detailed Structure ##

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
* variants DICT dictionary with all the variants that exists in the family on the form {<var_id>:<Variant>}
